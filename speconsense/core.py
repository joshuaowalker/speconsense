#!/usr/bin/env python3

from collections import defaultdict
import argparse
import json
import logging
import os
import random
import statistics
import subprocess
import sys
import tempfile
import math
import itertools
from typing import List, Set, Tuple, Optional, Dict, Any
from datetime import datetime

import edlib
from adjusted_identity import score_alignment, AdjustmentParams
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

# Import analysis functions for variance calculation and variant detection
from speconsense.analyze import (
    extract_alignments_from_msa,
    analyze_positional_variation,
    is_variant_position_with_composition
)

try:
    from speconsense import __version__
except ImportError:
    # Fallback for when running as a script directly (e.g., in tests)
    __version__ = "dev"

# IUPAC nucleotide ambiguity codes mapping
# Maps sets of nucleotides to their corresponding IUPAC code
IUPAC_CODES = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N',
}



class SpecimenClusterer:
    def __init__(self, min_identity: float = 0.8,
                 inflation: float = 4.0,
                 min_size: int = 5,
                 min_cluster_ratio: float = 0.2,
                 max_sample_size: int = 500,
                 presample_size: int = 1000,
                 k_nearest_neighbors: int = 20,
                 sample_name: str = "sample",
                 disable_homopolymer_equivalence: bool = False,
                 output_dir: str = "clusters",
                 mean_error_rate: Optional[float] = None):
        self.min_identity = min_identity
        self.inflation = inflation
        self.min_size = min_size
        self.min_cluster_ratio = min_cluster_ratio
        self.max_sample_size = max_sample_size
        self.presample_size = presample_size
        self.k_nearest_neighbors = k_nearest_neighbors
        self.sample_name = sample_name
        self.disable_homopolymer_equivalence = disable_homopolymer_equivalence
        self.output_dir = output_dir
        self.mean_error_rate = mean_error_rate
        self.sequences = {}  # id -> sequence string
        self.records = {}  # id -> SeqRecord object
        self.id_map = {}  # short_id -> original_id
        self.rev_id_map = {}  # original_id -> short_id

        # Calculate outlier threshold from mean error rate
        self.outlier_threshold = None
        if mean_error_rate is not None:
            OUTLIER_THRESHOLD_MULTIPLIER = 3.0
            self.outlier_threshold = mean_error_rate * OUTLIER_THRESHOLD_MULTIPLIER

        # Create output directory and debug subdirectory
        os.makedirs(self.output_dir, exist_ok=True)
        self.debug_dir = os.path.join(self.output_dir, "cluster_debug")
        os.makedirs(self.debug_dir, exist_ok=True)

        # Initialize attributes that may be set later
        self.input_file = None
        self.augment_input = None
        self.algorithm = None
        self.orient_mode = None
        self.primers_file = None

    def write_metadata(self) -> None:
        """Write run metadata to JSON file for use by post-processing tools."""
        metadata = {
            "version": __version__,
            "timestamp": datetime.now().isoformat(),
            "sample_name": self.sample_name,
            "parameters": {
                "algorithm": self.algorithm,
                "min_identity": self.min_identity,
                "inflation": self.inflation,
                "min_size": self.min_size,
                "min_cluster_ratio": self.min_cluster_ratio,
                "max_sample_size": self.max_sample_size,
                "presample_size": self.presample_size,
                "k_nearest_neighbors": self.k_nearest_neighbors,
                "disable_homopolymer_equivalence": self.disable_homopolymer_equivalence,
                "mean_error_rate": self.mean_error_rate,
                "outlier_threshold": self.outlier_threshold,
                "orient_mode": self.orient_mode,
            },
            "input_file": self.input_file,
            "augment_input": self.augment_input,
        }

        # Add primer information if loaded
        if hasattr(self, 'primers') and self.primers:
            metadata["primers_file"] = self.primers_file
            metadata["primers"] = {}

            # Store primer sequences (avoid duplicates from RC versions)
            seen_primers = set()
            for primer_name, primer_seq in self.primers:
                # Skip RC versions (they end with _RC)
                if not primer_name.endswith('_RC') and primer_name not in seen_primers:
                    metadata["primers"][primer_name] = primer_seq
                    seen_primers.add(primer_name)

        # Write metadata file
        metadata_file = os.path.join(self.debug_dir, f"{self.sample_name}-metadata.json")
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logging.info(f"Wrote run metadata to {metadata_file}")

    def write_phasing_stats(self, initial_clusters_count: int, subclusters_count: int,
                           merged_count: int, final_count: int) -> None:
        """Write phasing statistics to JSON file after clustering completes.

        Args:
            initial_clusters_count: Number of clusters from initial clustering
            subclusters_count: Number of sub-clusters after phasing
            merged_count: Number of clusters after merging
            final_count: Number of final clusters after filtering
        """
        phasing_stats = {
            "phasing_enabled": self.outlier_threshold is not None,
            "initial_clusters": initial_clusters_count,
            "phased_subclusters": subclusters_count,
            "after_merging": merged_count,
            "after_filtering": final_count,
            "clusters_split": subclusters_count > initial_clusters_count,
            "clusters_merged": merged_count < subclusters_count,
            "net_change": final_count - initial_clusters_count
        }

        # Write phasing stats to separate JSON file
        stats_file = os.path.join(self.debug_dir, f"{self.sample_name}-phasing_stats.json")
        with open(stats_file, 'w') as f:
            json.dump(phasing_stats, f, indent=2)

        logging.info(f"Wrote phasing statistics to {stats_file}")

    def add_sequences(self, records: List[SeqIO.SeqRecord],
                      augment_records: Optional[List[SeqIO.SeqRecord]] = None) -> None:
        """Add sequences to be clustered, with optional presampling."""
        all_records = records.copy()  # Start with primary records

        # Track the source of each record for potential logging/debugging
        primary_count = len(records)
        augment_count = 0

        # Add augmented records if provided
        if augment_records:
            augment_count = len(augment_records)
            all_records.extend(augment_records)

        if self.presample_size and len(all_records) > self.presample_size:
            logging.info(f"Presampling {self.presample_size} sequences from {len(all_records)} total "
                         f"({primary_count} primary, {augment_count} augmented)")

            # First, sort primary sequences by quality and take as many as possible
            primary_sorted = sorted(
                records,
                key=lambda r: -statistics.mean(r.letter_annotations["phred_quality"])
            )

            # Determine how many primary sequences to include (all if possible)
            primary_to_include = min(len(primary_sorted), self.presample_size)
            presampled = primary_sorted[:primary_to_include]

            # If we still have room, add augmented sequences sorted by quality
            remaining_slots = self.presample_size - primary_to_include
            if remaining_slots > 0 and augment_records:
                augment_sorted = sorted(
                    augment_records,
                    key=lambda r: -statistics.mean(r.letter_annotations["phred_quality"])
                )
                presampled.extend(augment_sorted[:remaining_slots])

            logging.info(f"Presampled {len(presampled)} sequences "
                         f"({primary_to_include} primary, {len(presampled) - primary_to_include} augmented)")
            all_records = presampled

        # Add all selected records to internal storage
        for record in all_records:
            self.sequences[record.id] = str(record.seq)
            self.records[record.id] = record

    def write_mcl_input(self, output_file: str) -> None:
        """Write similarity matrix in MCL input format using k-nearest neighbors approach."""
        self._create_id_mapping()

        # Calculate all pairwise similarities
        similarities = {}
        n = len(self.sequences)
        total_comparisons = (n * (n - 1)) // 2

        with tqdm(total=total_comparisons, desc="Calculating pairwise sequence similarities") as pbar:
            for id1 in self.sequences:
                similarities[id1] = {}
                for id2 in self.sequences:
                    if id1 >= id2:  # Only calculate upper triangle
                        continue
                    sim = self.calculate_similarity(self.sequences[id1], self.sequences[id2])
                    similarities[id1][id2] = sim
                    similarities.setdefault(id2, {})[id1] = sim  # Mirror for easy lookup
                    pbar.update(1)

        # Use k-nearest neighbors for each sequence to create sparse graph
        k = min(self.k_nearest_neighbors, n - 1)  # Connect to at most k neighbors
        min_edges_per_node = 3  # Ensure at least some connections per node

        with open(output_file, 'w') as f:
            for id1 in self.sequences:
                short_id1 = self.rev_id_map[id1]

                # Sort neighbors by similarity
                neighbors = sorted(
                    [(id2, sim) for id2, sim in similarities[id1].items()],
                    key=lambda x: x[1], reverse=True
                )

                # Take top k neighbors with sufficient similarity
                top_neighbors = [
                    (id2, sim) for id2, sim in neighbors[:k]
                    if sim >= self.min_identity
                ]

                # Ensure at least min_edges connections if possible
                if len(top_neighbors) < min_edges_per_node and len(neighbors) >= min_edges_per_node:
                    additional_needed = min_edges_per_node - len(top_neighbors)
                    for id2, sim in neighbors[k:k + additional_needed]:
                        if sim < self.min_identity * 0.9:  # Allow slightly lower threshold
                            break
                        top_neighbors.append((id2, sim))

                # Write edges
                for id2, sim in top_neighbors:
                    short_id2 = self.rev_id_map[id2]
                    # Transform similarity to emphasize differences
                    transformed_sim = sim ** 2  # Square the similarity
                    f.write(f"{short_id1}\t{short_id2}\t{transformed_sim:.6f}\n")

    def run_mcl(self, input_file: str, output_file: str) -> None:
        """Run MCL clustering algorithm with optimized parameters."""
        cmd = [
            "mcl",
            input_file,
            "--abc",  # Input is in ABC format (node1 node2 weight)
            "-I", str(self.inflation),  # Inflation parameter
            "-scheme", "7",  # More advanced flow simulation
            "-pct", "50",  # Prune weakest 50% of connections during iterations
            "-o", output_file  # Output file
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            logging.debug(f"MCL stdout: {result.stdout}")
            logging.debug(f"MCL stderr: {result.stderr}")

        except subprocess.CalledProcessError as e:
            logging.error(f"MCL failed with return code {e.returncode}")
            logging.error(f"Command: {' '.join(cmd)}")
            logging.error(f"Stderr: {e.stderr}")
            raise

    def merge_similar_clusters(self, clusters: List[Dict]) -> List[Dict]:
        """
        Merge clusters whose consensus sequences are identical or homopolymer-equivalent.
        Preserves provenance metadata through the merging process.

        Note: Primer trimming is performed before comparison to ensure clusters that differ
        only in primer regions are properly merged. Trimmed consensuses are used only for
        comparison and are discarded after merging.

        Args:
            clusters: List of cluster dictionaries with 'read_ids' and provenance fields

        Returns:
            List of merged cluster dictionaries with combined provenance
        """
        if not clusters:
            return []

        # Sort clusters by size, largest first
        clusters = sorted(clusters, key=lambda c: len(c['read_ids']), reverse=True)

        # Generate a consensus sequence for each cluster
        logging.info("Generating consensus sequences for cluster merging...")
        consensuses = []
        cluster_to_consensus = {}  # Map from cluster index to its consensus

        with tqdm(total=len(clusters), desc="Generating cluster consensuses") as pbar:
            for i, cluster_dict in enumerate(clusters):
                cluster_reads = cluster_dict['read_ids']

                # Sample from larger clusters to speed up consensus generation
                if len(cluster_reads) > self.max_sample_size:
                    # Sample by quality
                    qualities = []
                    for seq_id in cluster_reads:
                        record = self.records[seq_id]
                        mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                        qualities.append((mean_quality, seq_id))

                    sampled_ids = [seq_id for _, seq_id in
                                   sorted(qualities, reverse=True)[:self.max_sample_size]]
                    sampled_seqs = [self.sequences[seq_id] for seq_id in sampled_ids]
                else:
                    sampled_seqs = [self.sequences[seq_id] for seq_id in cluster_reads]

                consensus, _ = self.run_spoa(sampled_seqs)  # Ignore MSA for merging

                # Skip empty clusters (can occur when all sequences are filtered out)
                if consensus is None:
                    logging.warning(f"Cluster {i} produced no consensus (empty cluster), skipping")
                    pbar.update(1)
                    continue

                # Trim primers before comparison to merge clusters that differ only in primer regions
                # The trimmed consensus is used only for comparison and discarded after merging
                if hasattr(self, 'primers'):
                    consensus, _ = self.trim_primers(consensus)  # Discard found_primers

                consensuses.append(consensus)
                cluster_to_consensus[i] = consensus
                pbar.update(1)

        consensus_to_clusters = defaultdict(list)
        
        if self.disable_homopolymer_equivalence:
            # Only merge exactly identical sequences
            for i, consensus in enumerate(consensuses):
                consensus_to_clusters[consensus].append(i)
        else:
            # Group by homopolymer-equivalent sequences
            for i, consensus in enumerate(consensuses):
                # Find if this consensus is homopolymer-equivalent to any existing group
                found_group = False
                for existing_consensus in consensus_to_clusters.keys():
                    if self.are_homopolymer_equivalent(consensus, existing_consensus):
                        consensus_to_clusters[existing_consensus].append(i)
                        found_group = True
                        break
                
                if not found_group:
                    consensus_to_clusters[consensus].append(i)

        merged = []
        merged_indices = set()

        # Handle clusters with equivalent consensus sequences
        for equivalent_clusters in consensus_to_clusters.values():
            if len(equivalent_clusters) > 1:
                # Merge clusters with equivalent consensus
                merged_read_ids = set()
                merged_from_list = []

                for idx in equivalent_clusters:
                    merged_read_ids.update(clusters[idx]['read_ids'])
                    merged_indices.add(idx)

                    # Track what we're merging from
                    cluster_info = {
                        'initial_cluster_num': clusters[idx]['initial_cluster_num'],
                        'allele_combo': clusters[idx].get('allele_combo'),
                        'size': len(clusters[idx]['read_ids'])
                    }
                    merged_from_list.append(cluster_info)

                # Create merged cluster with provenance
                merged_cluster = {
                    'read_ids': merged_read_ids,
                    'initial_cluster_num': None,  # Multiple sources
                    'allele_combo': None,  # Multiple alleles merged
                    'variant_positions': None,
                    'num_variants': 0,
                    'merged_from': merged_from_list  # Track merge provenance
                }
                merged.append(merged_cluster)

        # Add remaining unmerged clusters
        for i, cluster_dict in enumerate(clusters):
            if i not in merged_indices:
                merged.append(cluster_dict)

        merge_type = "identical" if self.disable_homopolymer_equivalence else "homopolymer-equivalent"
        if len(merged) < len(clusters):
            logging.info(f"Merged {len(clusters) - len(merged)} clusters with {merge_type} consensus sequences")
        else:
            logging.info(f"No clusters were merged (no {merge_type} consensus sequences found)")

        return merged

    def _find_root(self, merged_to: List[int], i: int) -> int:
        """Find the root index of a merged cluster using path compression."""
        if merged_to[i] != i:
            merged_to[i] = self._find_root(merged_to, merged_to[i])
        return merged_to[i]

    def write_cluster_files(self, cluster_num: int, cluster: Set[str],
                            consensus: str, trimmed_consensus: Optional[str] = None,
                            found_primers: Optional[List[str]] = None,
                            mean_var: Optional[float] = None,
                            median_var: Optional[float] = None,
                            p95_var: Optional[float] = None,
                            variant_positions: Optional[List[Dict]] = None,
                            actual_size: Optional[int] = None,
                            consensus_fasta_handle = None,
                            sampled_ids: Optional[Set[str]] = None,
                            msa: Optional[str] = None) -> None:
        """Write cluster files: reads FASTQ, MSA, consensus FASTA, and variant debug files.

        Variance metrics measure per-read edit distance to consensus:
        - Low variance: good clustering and accurate consensus
        - High variance: may indicate read errors OR imperfect clustering/consensus

        Variant positions indicate potential unphased biological variation that may
        require cluster subdivision in future processing.
        """
        cluster_size = len(cluster)
        ric_size = min(actual_size or cluster_size, self.max_sample_size)

        # Create info string with size first
        info_parts = [f"size={cluster_size}", f"ric={ric_size}"]

        # Add variance metrics (as percentages for readability)
        if mean_var is not None:
            info_parts.append(f"var_mean={mean_var*100:.1f}")
        if median_var is not None:
            info_parts.append(f"var_med={median_var*100:.1f}")
        if p95_var is not None:
            info_parts.append(f"var_p95={p95_var*100:.1f}")

        # Add variant flags
        if variant_positions is not None:
            num_variants = len(variant_positions)
            if num_variants > 0:
                info_parts.append(f"has_variants=true")
                info_parts.append(f"num_variants={num_variants}")
            else:
                info_parts.append(f"has_variants=false")

        if found_primers:
            info_parts.append(f"primers={','.join(found_primers)}")
        info_str = " ".join(info_parts)

        # Write reads FASTQ to debug directory with new naming convention
        reads_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-reads.fastq")
        with open(reads_file, 'w') as f:
            for seq_id in cluster:
                SeqIO.write(self.records[seq_id], f, "fastq")

        # Write sampled reads FASTQ (only sequences used for consensus generation)
        if sampled_ids is not None:
            sampled_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-sampled.fastq")
            with open(sampled_file, 'w') as f:
                for seq_id in sampled_ids:
                    SeqIO.write(self.records[seq_id], f, "fastq")

        # Write MSA (multiple sequence alignment) to debug directory
        if msa is not None:
            msa_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-msa.fasta")
            with open(msa_file, 'w') as f:
                f.write(msa)

        # Write variant positions to debug directory
        if variant_positions and len(variant_positions) > 0:
            variant_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-variants.txt")
            with open(variant_file, 'w') as f:
                f.write(f"Variant positions for {self.sample_name}-c{cluster_num}\n")
                f.write(f"Total variants: {len(variant_positions)}\n\n")

                for var in variant_positions:
                    # Format position info
                    if var['consensus_position'] is not None:
                        pos_str = f"MSA:{var['msa_position']}/Consensus:{var['consensus_position']}"
                    else:
                        pos_str = f"MSA:{var['msa_position']}/Consensus:insertion"

                    # Format base composition
                    sorted_bases = sorted(var['base_composition'].items(), key=lambda x: x[1], reverse=True)
                    comp_str = ', '.join([f"{base}:{count}" for base, count in sorted_bases[:4]])

                    # Write variant info
                    f.write(f"Position {pos_str}\n")
                    f.write(f"  Coverage: {var['coverage']}\n")
                    f.write(f"  Error rate: {var['error_rate']*100:.1f}%\n")
                    f.write(f"  Base composition: {comp_str}\n")
                    f.write(f"  Variant bases: {var['variant_bases']}\n")
                    f.write(f"  Reason: {var['reason']}\n\n")

        # Write untrimmed consensus to debug directory
        with open(os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-untrimmed.fasta"),
                  'w') as f:
            f.write(f">{self.sample_name}-c{cluster_num} {info_str}\n")
            f.write(consensus + "\n")

        # Write consensus to main output file if handle is provided
        if consensus_fasta_handle:
            final_consensus = trimmed_consensus if trimmed_consensus else consensus
            consensus_fasta_handle.write(f">{self.sample_name}-c{cluster_num} {info_str}\n")
            consensus_fasta_handle.write(final_consensus + "\n")

    def run_mcl_clustering(self, temp_dir: str) -> List[Set[str]]:
        """Run MCL clustering algorithm and return the clusters.

        Args:
            temp_dir: Path to temporary directory for intermediate files

        Returns:
            List of clusters, where each cluster is a set of sequence IDs
        """
        mcl_input = os.path.join(temp_dir, "input.abc")
        mcl_output = os.path.join(temp_dir, "output.cls")

        self.write_mcl_input(mcl_input)

        logging.info(f"Running MCL algorithm with inflation {self.inflation}...")
        self.run_mcl(mcl_input, mcl_output)
        return self.parse_mcl_output(mcl_output)

    def run_greedy_clustering(self, temp_dir: str) -> List[Set[str]]:
        """Run greedy clustering algorithm and return the clusters.

        This algorithm iteratively finds the sequence with the most connections above
        the similarity threshold and forms a cluster around it.

        Args:
            temp_dir: Path to temporary directory for intermediate files

        Returns:
            List of clusters, where each cluster is a set of sequence IDs
        """
        logging.info("Running greedy clustering algorithm...")

        # Build similarity matrix if not already built
        if not hasattr(self, 'alignments'):
            self.alignments = defaultdict(dict)
            self.build_similarity_matrix()

        # Initial clustering
        clusters = []
        available_ids = set(self.sequences.keys())

        cluster_count = 0

        while available_ids:
            center, members = self.find_cluster_center(available_ids)
            available_ids -= members

            clusters.append(members)
            cluster_count += 1

        return clusters

    def build_similarity_matrix(self) -> None:
        """Calculate all pairwise similarities between sequences."""
        logging.info("Calculating pairwise sequence similarities...")

        seq_ids = list(self.sequences.keys())
        total = len(seq_ids) * (len(seq_ids) - 1) // 2

        with tqdm(total=total, desc="Building similarity matrix") as pbar:
            for i, id1 in enumerate(seq_ids):
                for id2 in seq_ids[i + 1:]:
                    sim = self.calculate_similarity(
                        self.sequences[id1],
                        self.sequences[id2]
                    )

                    if sim >= self.min_identity:
                        self.alignments[id1][id2] = sim
                        self.alignments[id2][id1] = sim

                    pbar.update(1)

    def find_cluster_center(self, available_ids: Set[str]) -> Tuple[str, Set[str]]:
        """
        Find the sequence with most similar sequences above threshold,
        and return its ID and the IDs of its cluster members.
        """
        best_center = None
        best_members = set()
        best_count = -1

        for seq_id in available_ids:
            # Get all sequences that align with current sequence
            members = {other_id for other_id in self.alignments.get(seq_id, {})
                       if other_id in available_ids}

            if len(members) > best_count:
                best_count = len(members)
                best_center = seq_id
                best_members = members

        if best_center is None:
            # No alignments found, create singleton cluster with first available sequence
            singleton_id = next(iter(available_ids))
            return singleton_id, {singleton_id}

        best_members.add(best_center)  # Include center in cluster
        return best_center, best_members

    def cluster(self, algorithm: str = "graph") -> None:
        """Perform complete clustering process with variant phasing and write output files.

        Pipeline: Splitting (clustering + outlier removal + phasing) → Merging → Filtering → Output

        Args:
            algorithm: Clustering algorithm to use ('graph' for MCL or 'greedy')
        """
        # Create temporary directory outside the clustering method to allow
        # flexibility in choosing different clustering algorithms
        with tempfile.TemporaryDirectory() as temp_dir:
            # ============================================================
            # PHASE 1: INITIAL CLUSTERING
            # ============================================================
            if algorithm == "graph":
                try:
                    initial_clusters = self.run_mcl_clustering(temp_dir)
                except (subprocess.SubprocessError, FileNotFoundError) as e:
                    logging.error(f"MCL clustering failed: {str(e)}")
                    logging.error("You may need to install MCL: https://micans.org/mcl/")
                    logging.error("Falling back to greedy clustering algorithm...")
                    initial_clusters = self.run_greedy_clustering(temp_dir)
            elif algorithm == "greedy":
                initial_clusters = self.run_greedy_clustering(temp_dir)
            else:
                raise ValueError(f"Unknown clustering algorithm: {algorithm}")

            # Sort initial clusters by size (largest first)
            initial_clusters.sort(key=lambda c: len(c), reverse=True)

            logging.info(f"Initial clustering produced {len(initial_clusters)} clusters")

            # ============================================================
            # PHASE 2: DETECTION + PHASING (Splitting Phase)
            # ============================================================
            # Process each initial cluster to detect variants and phase reads
            all_subclusters = []  # Will hold all phased sub-clusters with provenance

            logging.info("Processing clusters for variant detection and phasing...")

            for initial_idx, cluster in enumerate(initial_clusters, 1):
                cluster_size = len(cluster)

                # Sample sequences for consensus generation if needed
                if cluster_size > self.max_sample_size:
                    logging.debug(f"Initial cluster {initial_idx}: Sampling {self.max_sample_size} from {cluster_size} reads")
                    qualities = []
                    for seq_id in cluster:
                        record = self.records[seq_id]
                        mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                        qualities.append((mean_quality, seq_id))
                    sampled_ids = {seq_id for _, seq_id in
                                   sorted(qualities, reverse=True)[:self.max_sample_size]}
                else:
                    sampled_ids = cluster

                # Generate consensus and MSA
                sampled_seqs = [self.sequences[seq_id] for seq_id in sampled_ids]
                consensus, msa = self.run_spoa(sampled_seqs)

                if not consensus or not msa:
                    logging.warning(f"Initial cluster {initial_idx}: Failed to generate consensus, skipping")
                    continue

                # Calculate per-read variance from MSA
                mean_var, median_var, p95_var = self.calculate_read_variance(msa, consensus)

                # Optional: Remove outlier reads and regenerate consensus
                if self.outlier_threshold is not None:
                    keep_ids, outlier_ids = self.identify_outlier_reads(
                        msa, consensus, sampled_ids, self.outlier_threshold
                    )

                    if outlier_ids:
                        logging.info(f"Initial cluster {initial_idx}: Removing {len(outlier_ids)}/{len(sampled_ids)} outlier reads, "
                                   f"regenerating consensus")

                        # Update sampled_ids to exclude outliers
                        sampled_ids = keep_ids

                        # CRITICAL: Also remove outliers from the full cluster
                        # This ensures outliers don't reappear in phased haplotypes
                        cluster = cluster - outlier_ids

                        # Regenerate consensus with filtered reads
                        sampled_seqs = [self.sequences[seq_id] for seq_id in sampled_ids]
                        consensus, msa = self.run_spoa(sampled_seqs)

                        # Recalculate variance
                        if consensus and msa:
                            mean_var, median_var, p95_var = self.calculate_read_variance(msa, consensus)

                # Detect variant positions
                variant_positions = []
                if consensus and msa and self.outlier_threshold is not None:
                    variant_positions = self.detect_variant_positions(
                        msa, consensus, self.outlier_threshold, min_alt_freq=0.20
                    )

                    if variant_positions:
                        logging.info(f"Initial cluster {initial_idx}: Detected {len(variant_positions)} variant positions")

                # Phase reads into haplotypes
                # Note: Uses the FULL cluster (not just sampled_ids) for phasing
                phased_haplotypes = self.phase_reads_by_variants(
                    msa, consensus, cluster, variant_positions
                )

                # Store each haplotype as a sub-cluster with provenance
                for haplotype_idx, (allele_combo, haplotype_reads) in enumerate(phased_haplotypes):
                    subcluster = {
                        'read_ids': haplotype_reads,
                        'initial_cluster_num': initial_idx,
                        'allele_combo': allele_combo,
                        'variant_positions': variant_positions if allele_combo else None,
                        'num_variants': len(variant_positions) if variant_positions else 0
                    }
                    all_subclusters.append(subcluster)

            logging.info(f"After phasing, created {len(all_subclusters)} sub-clusters from {len(initial_clusters)} initial clusters")

            # ============================================================
            # PHASE 3: MERGING
            # ============================================================
            # Merge sub-clusters with identical/homopolymer-equivalent consensus sequences
            logging.info("Merging homopolymer-equivalent sub-clusters...")
            merged_subclusters = self.merge_similar_clusters(all_subclusters)

            logging.info(f"After merging, have {len(merged_subclusters)} clusters")

            # ============================================================
            # PHASE 4: FILTERING
            # ============================================================
            # Filter by absolute size
            large_clusters = [c for c in merged_subclusters if len(c['read_ids']) >= self.min_size]

            if len(large_clusters) < len(merged_subclusters):
                filtered_count = len(merged_subclusters) - len(large_clusters)
                logging.info(f"Filtered {filtered_count} clusters below minimum size ({self.min_size})")

            # Filter by relative size ratio
            if large_clusters and self.min_cluster_ratio > 0:
                largest_size = max(len(c['read_ids']) for c in large_clusters)
                before_ratio_filter = len(large_clusters)
                large_clusters = [c for c in large_clusters
                                 if len(c['read_ids']) / largest_size >= self.min_cluster_ratio]

                if len(large_clusters) < before_ratio_filter:
                    filtered_count = before_ratio_filter - len(large_clusters)
                    logging.info(f"Filtered {filtered_count} clusters below minimum ratio ({self.min_cluster_ratio})")

            # Sort by size and renumber as c1, c2, c3...
            large_clusters.sort(key=lambda c: len(c['read_ids']), reverse=True)

            total_sequences = len(self.sequences)
            sequences_covered = sum(len(c['read_ids']) for c in large_clusters)

            logging.info(f"Final: {len(large_clusters)} clusters covering {sequences_covered} sequences "
                        f"({sequences_covered / total_sequences:.1%} of total)")

            # ============================================================
            # PHASE 5: OUTPUT
            # ============================================================
            # Generate final consensus and write output files
            consensus_output_file = os.path.join(self.output_dir, f"{self.sample_name}-all.fasta")
            with open(consensus_output_file, 'w') as consensus_fasta_handle:
                for final_idx, cluster_dict in enumerate(large_clusters, 1):
                    cluster = cluster_dict['read_ids']
                    actual_size = len(cluster)

                    # Sample sequences for final consensus generation if needed
                    if len(cluster) > self.max_sample_size:
                        logging.info(f"Cluster {final_idx}: Sampling {self.max_sample_size} from {len(cluster)} reads for final consensus")
                        qualities = []
                        for seq_id in cluster:
                            record = self.records[seq_id]
                            mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                            qualities.append((mean_quality, seq_id))
                        sampled_ids = {seq_id for _, seq_id in
                                       sorted(qualities, reverse=True)[:self.max_sample_size]}
                    else:
                        sampled_ids = cluster

                    # Generate final consensus and MSA
                    sampled_seqs = [self.sequences[seq_id] for seq_id in sampled_ids]
                    consensus, msa = self.run_spoa(sampled_seqs)

                    # Calculate final variance metrics
                    mean_var, median_var, p95_var = None, None, None
                    if consensus and msa:
                        mean_var, median_var, p95_var = self.calculate_read_variance(msa, consensus)

                    # Note: No outlier removal or variant detection in output phase
                    # (already done in detection phase)

                    if consensus:
                        # Perform primer trimming
                        trimmed_consensus = None
                        found_primers = None
                        if hasattr(self, 'primers'):
                            trimmed_consensus, found_primers = self.trim_primers(consensus)

                        # Write output files (no variant positions in output phase)
                        self.write_cluster_files(
                            cluster_num=final_idx,
                            cluster=cluster,
                            consensus=consensus,
                            trimmed_consensus=trimmed_consensus,
                            found_primers=found_primers,
                            mean_var=mean_var,
                            median_var=median_var,
                            p95_var=p95_var,
                            variant_positions=None,  # Don't report variants in final output
                            actual_size=actual_size,
                            consensus_fasta_handle=consensus_fasta_handle,
                            sampled_ids=sampled_ids,
                            msa=msa
                        )

            # ============================================================
            # WRITE PHASING STATISTICS
            # ============================================================
            self.write_phasing_stats(
                initial_clusters_count=len(initial_clusters),
                subclusters_count=len(all_subclusters),
                merged_count=len(merged_subclusters),
                final_count=len(large_clusters)
            )

    def _create_id_mapping(self) -> None:
        """Create short numeric IDs for all sequences."""
        for i, seq_id in enumerate(self.sequences.keys()):
            short_id = str(i)
            self.id_map[short_id] = seq_id
            self.rev_id_map[seq_id] = short_id

    def calculate_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity using edlib alignment."""
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0

        max_dist = int((1 - self.min_identity) * max(len(seq1), len(seq2)))
        result = edlib.align(seq1, seq2, task="distance", k=max_dist)

        if result["editDistance"] == -1:
            return 0.0

        return 1.0 - (result["editDistance"] / max(len(seq1), len(seq2)))

    def run_spoa(self, sequences: List[str]) -> Tuple[Optional[str], Optional[str]]:
        """Run SPOA to generate consensus sequence and MSA.

        Returns:
            Tuple of (consensus_sequence, msa_fasta) where msa_fasta includes all
            aligned sequences and the consensus
        """
        if not sequences:
            return None, None

        try:
            # Create temporary input file
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
                for i, seq in enumerate(sequences):
                    f.write(f">seq{i}\n{seq}\n")
                temp_input = f.name

            # Construct SPOA command with parameters
            cmd = [
                "spoa",
                temp_input,
                "-r", "2",  # Result mode 2: MSA + consensus (consensus is last)
                "-l", "1",  # Global alignment mode
                "-m", "5",  # Match score
                "-n", "-4",  # Mismatch penalty
                "-g", "-8",  # Gap opening penalty
                "-e", "-6",  # Gap extension penalty
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # Clean up input file
            os.unlink(temp_input)

            # Parse SPOA output - capture both MSA and consensus
            msa_output = result.stdout
            consensus = None
            current_seq = []
            is_consensus = False
            for line in result.stdout.split('\n'):
                if line.startswith('>'):
                    # Check if this is the consensus sequence
                    is_consensus = 'Consensus' in line
                    current_seq = []
                elif line.strip() and is_consensus:
                    current_seq.append(line.strip())
            # Extract the consensus sequence and remove gaps
            # (consensus is part of MSA so will have gaps, but we want ungapped sequence)
            if current_seq:
                consensus = ''.join(current_seq).replace('-', '')

            if not consensus:
                logging.warning("SPOA did not generate consensus sequence")
                return sequences[0], msa_output  # Fall back to first sequence

            return consensus, msa_output

        except subprocess.CalledProcessError as e:
            logging.error(f"SPOA failed with return code {e.returncode}")
            logging.error(f"Command: {' '.join(cmd)}")
            logging.error(f"Stderr: {e.stderr}")
            return sequences[0], None  # Fall back to first sequence

        except Exception as e:
            logging.error(f"Error running SPOA: {str(e)}")
            return sequences[0], None  # Fall back to first sequence

    def identify_outlier_reads(self, msa_string: str, consensus_seq: str, sampled_ids: Set[str],
                              threshold: float) -> Tuple[Set[str], Set[str]]:
        """Identify outlier reads that exceed error rate threshold.

        Args:
            msa_string: MSA in FASTA format from SPOA
            consensus_seq: Ungapped consensus sequence
            sampled_ids: Set of read IDs that were sampled (in order written to SPOA)
            threshold: Error rate threshold (e.g., mean × 3.0)

        Returns:
            Tuple of (keep_ids, outlier_ids) where both are sets of actual read IDs
        """
        if not msa_string or not consensus_seq:
            return sampled_ids, set()

        try:
            # Extract alignments from MSA string directly (no temp file needed)
            alignments, _, _ = extract_alignments_from_msa(msa_string)

            if not alignments:
                logging.debug(f"No alignments found in MSA for outlier detection")
                return sampled_ids, set()

            # Map MSA IDs (seq0, seq1, ...) to actual read IDs
            # SPOA writes sequences in the order they were provided
            sampled_list = list(sampled_ids)
            msa_to_read_id = {f"seq{i}": sampled_list[i] for i in range(len(sampled_list))}

            keep_ids = set()
            outlier_ids = set()

            consensus_length = len(consensus_seq)
            if consensus_length == 0:
                return sampled_ids, set()

            for alignment in alignments:
                error_rate = alignment.edit_distance / consensus_length

                # Get actual read ID
                actual_read_id = msa_to_read_id.get(alignment.read_id, alignment.read_id)

                if error_rate <= threshold:
                    keep_ids.add(actual_read_id)
                else:
                    outlier_ids.add(actual_read_id)

            return keep_ids, outlier_ids

        except Exception as e:
            logging.warning(f"Failed to identify outlier reads: {e}")
            return sampled_ids, set()

    def calculate_read_variance(self, msa_string: str, consensus_seq: str) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """Calculate per-read variance (edit distance to consensus) from MSA.

        Variance measures the difference between individual reads and the consensus sequence.
        In well-clustered data with accurate consensus, this approximates the sequencing error rate.
        High variance may indicate either high sequencing error OR imperfect clustering/consensus.

        Args:
            msa_string: MSA in FASTA format (string) generated by SPOA
            consensus_seq: Ungapped consensus sequence

        Returns:
            Tuple of (mean_variance, median_variance, p95_variance) as fractions (0.0-1.0)
            Returns (None, None, None) if MSA cannot be parsed or has no valid reads
        """
        if not msa_string or not consensus_seq:
            return None, None, None

        try:
            # Extract alignments from MSA string directly (no temp file needed)
            alignments, msa_consensus, _ = extract_alignments_from_msa(msa_string)

            if not alignments:
                logging.debug(f"No alignments found in MSA")
                return None, None, None

            # Calculate variance (edit distance / consensus length) for each read
            consensus_length = len(consensus_seq)
            if consensus_length == 0:
                return None, None, None

            variances = []
            for alignment in alignments:
                variance = alignment.edit_distance / consensus_length
                variances.append(variance)

            if not variances:
                return None, None, None

            # Calculate statistics
            mean_var = np.mean(variances)
            median_var = np.median(variances)
            p95_var = np.percentile(variances, 95)

            return mean_var, median_var, p95_var

        except Exception as e:
            logging.warning(f"Failed to calculate read variance from MSA: {e}")
            return None, None, None

    def detect_variant_positions(self, msa_string: str, consensus_seq: str,
                                error_threshold: float, min_alt_freq: float = 0.20) -> List[Dict]:
        """Detect variant positions in MSA for future phasing.

        Identifies positions with significant alternative alleles that may indicate
        unphased biological variation requiring cluster subdivision.

        Args:
            msa_string: MSA in FASTA format from SPOA
            consensus_seq: Ungapped consensus sequence
            error_threshold: Error rate threshold for variant detection (e.g., P95)
            min_alt_freq: Minimum alternative allele frequency (default: 0.20)

        Returns:
            List of variant position dictionaries with metadata, or empty list if none found
        """
        if not msa_string or not consensus_seq:
            return []

        try:
            # Extract alignments from MSA string directly (no temp file needed)
            alignments, msa_consensus, msa_to_consensus_pos = extract_alignments_from_msa(msa_string)

            if not alignments:
                logging.debug(f"No alignments found in MSA for variant detection")
                return []

            # Analyze positional variation
            position_stats = analyze_positional_variation(
                alignments,
                consensus_seq,
                error_threshold,
                msa_string
            )

            # Identify variant positions
            variant_positions = []
            for pos_stat in position_stats:
                is_variant, variant_bases, reason = is_variant_position_with_composition(
                    pos_stat,
                    error_threshold,
                    min_alt_freq=min_alt_freq,
                    confidence_sigma=3.0
                )

                if is_variant:
                    variant_positions.append({
                        'msa_position': pos_stat.msa_position,
                        'consensus_position': pos_stat.consensus_position,
                        'coverage': pos_stat.coverage,
                        'variant_bases': variant_bases,
                        'base_composition': pos_stat.base_composition,
                        'error_rate': pos_stat.error_rate,
                        'reason': reason
                    })

            return variant_positions

        except Exception as e:
            logging.warning(f"Failed to detect variant positions: {e}")
            return []

    def phase_reads_by_variants(
        self,
        msa_string: str,
        consensus_seq: str,
        cluster_read_ids: Set[str],
        variant_positions: List[Dict]
    ) -> List[Tuple[str, Set[str]]]:
        """Phase reads into haplotypes based on variant positions.

        Splits a cluster into sub-clusters where each sub-cluster represents reads
        sharing the same alleles at all variant positions. Every unique combination
        of alleles observed in the data becomes a separate haplotype.

        Args:
            msa_string: MSA in FASTA format from SPOA
            consensus_seq: Ungapped consensus sequence
            cluster_read_ids: Set of read IDs in this cluster (for mapping MSA IDs)
            variant_positions: List of variant position dicts from detect_variant_positions()

        Returns:
            List of (allele_combo_string, read_id_set) tuples
            e.g., [("C-T-A", {id1, id2}), ("T-C-A", {id3, id4})]

            Allele combination format: bases separated by hyphens, in MSA position order
            Gap positions are represented as "-" (the gap character itself)
        """
        if not msa_string or not variant_positions:
            # No variants to phase - return single group with all reads
            return [(None, cluster_read_ids)]

        try:
            # Parse MSA to get aligned sequences
            from io import StringIO
            from Bio import SeqIO

            msa_handle = StringIO(msa_string)
            records = list(SeqIO.parse(msa_handle, 'fasta'))

            if not records:
                logging.warning("No sequences found in MSA for phasing")
                return [(None, cluster_read_ids)]

            # Separate consensus from read sequences
            read_sequences = {}
            for record in records:
                if 'Consensus' not in record.description and 'Consensus' not in record.id:
                    read_sequences[record.id] = str(record.seq).upper()

            if not read_sequences:
                logging.warning("No read sequences found in MSA for phasing")
                return [(None, cluster_read_ids)]

            # Map MSA IDs (seq0, seq1, ...) to actual read IDs
            # SPOA writes sequences in the order they were provided
            cluster_read_list = list(cluster_read_ids)
            msa_to_read_id = {f"seq{i}": cluster_read_list[i] for i in range(len(cluster_read_list))}

            # Extract MSA positions for variants (sorted by position)
            variant_msa_positions = sorted([v['msa_position'] for v in variant_positions])

            # For each read, extract alleles at variant positions
            read_to_alleles = {}
            for msa_id, aligned_seq in read_sequences.items():
                # Get actual read ID
                actual_read_id = msa_to_read_id.get(msa_id, msa_id)

                # Extract allele at each variant position
                alleles = []
                for msa_pos in variant_msa_positions:
                    if msa_pos < len(aligned_seq):
                        allele = aligned_seq[msa_pos]
                        alleles.append(allele)
                    else:
                        # Read doesn't cover this position - treat as gap
                        alleles.append('-')

                # Create allele combination string
                allele_combo = '-'.join(alleles)
                read_to_alleles[actual_read_id] = allele_combo

            # Group reads by allele combination
            combo_to_reads = defaultdict(set)
            for read_id, allele_combo in read_to_alleles.items():
                combo_to_reads[allele_combo].add(read_id)

            # Convert to list of tuples
            result = [(combo, reads) for combo, reads in combo_to_reads.items()]

            # Sort by size (largest first) for consistency
            result.sort(key=lambda x: len(x[1]), reverse=True)

            logging.debug(f"Phased {len(cluster_read_ids)} reads into {len(result)} haplotypes "
                         f"based on {len(variant_positions)} variant positions")

            return result

        except Exception as e:
            logging.warning(f"Failed to phase reads by variants: {e}")
            # On error, return all reads as single group
            return [(None, cluster_read_ids)]

    def load_primers(self, primer_file: str) -> None:
        """Load primers from FASTA file with position awareness."""
        # Store primers in separate lists by position
        self.forward_primers = []
        self.reverse_primers = []
        self.forward_primers_rc = []  # RC of forward primers
        self.reverse_primers_rc = []  # RC of reverse primers
        
        # For backward compatibility with trim_primers
        self.primers = []  # Will be populated with all primers for existing code
        
        try:
            primer_count = {'forward': 0, 'reverse': 0, 'unknown': 0}
            
            for record in SeqIO.parse(primer_file, "fasta"):
                sequence = str(record.seq)
                sequence_rc = str(reverse_complement(sequence))
                
                # Parse position from header
                if "position=forward" in record.description:
                    self.forward_primers.append((record.id, sequence))
                    self.forward_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    primer_count['forward'] += 1
                elif "position=reverse" in record.description:
                    self.reverse_primers.append((record.id, sequence))
                    self.reverse_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    primer_count['reverse'] += 1
                else:
                    # For primers without position info, add to both lists
                    logging.warning(f"Primer {record.id} has no position annotation, treating as bidirectional")
                    self.forward_primers.append((record.id, sequence))
                    self.forward_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    self.reverse_primers.append((record.id, sequence))
                    self.reverse_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    primer_count['unknown'] += 1
                
                # Maintain backward compatibility
                self.primers.append((record.id, sequence))
                self.primers.append((f"{record.id}_RC", sequence_rc))

            total_primers = sum(primer_count.values())
            if total_primers == 0:
                logging.warning("No primers were loaded. Primer trimming will be disabled.")
            else:
                logging.info(f"Loaded {total_primers} primers: {primer_count['forward']} forward, "
                           f"{primer_count['reverse']} reverse, {primer_count['unknown']} unknown")
        except Exception as e:
            logging.error(f"Error loading primers: {str(e)}")
            raise

    def orient_sequences(self) -> set:
        """Normalize sequence orientations based on primer matches.
        
        Scoring system:
        - +1 point if a forward primer is found at the expected position
        - +1 point if a reverse primer is found at the expected position
        - Maximum score: 2 (both primers found)
        
        Decision logic:
        - If one orientation scores >0 and the other scores 0: use the non-zero orientation
        - If both score 0 or both score >0: keep original orientation (ambiguous/failed)
        """
        if not hasattr(self, 'forward_primers') or not hasattr(self, 'reverse_primers'):
            logging.warning("No positioned primers loaded, skipping orientation")
            return set()
            
        if len(self.forward_primers) == 0 and len(self.reverse_primers) == 0:
            logging.warning("No positioned primers available, skipping orientation")
            return set()
        
        logging.info("Starting sequence orientation based on primer positions...")
        
        oriented_count = 0
        already_correct = 0
        failed_count = 0
        failed_sequences = set()  # Track which sequences failed orientation
        
        # Process each sequence
        for seq_id in tqdm(self.sequences, desc="Orienting sequences"):
            sequence = self.sequences[seq_id]
            
            # Test both orientations (scores will be 0, 1, or 2)
            forward_score = self._score_orientation(sequence, "forward")
            reverse_score = self._score_orientation(sequence, "reverse")
            
            # Decision logic
            if forward_score > 0 and reverse_score == 0:
                # Clear forward orientation
                already_correct += 1
                logging.debug(f"Kept {seq_id} as-is: forward_score={forward_score}, reverse_score={reverse_score}")
            elif reverse_score > 0 and forward_score == 0:
                # Clear reverse orientation - needs to be flipped
                self.sequences[seq_id] = str(reverse_complement(sequence))
                
                # Also update the record if it exists
                if seq_id in self.records:
                    record = self.records[seq_id]
                    record.seq = reverse_complement(record.seq)
                    # Reverse quality scores too if they exist
                    if 'phred_quality' in record.letter_annotations:
                        record.letter_annotations['phred_quality'] = \
                            record.letter_annotations['phred_quality'][::-1]
                
                oriented_count += 1
                logging.debug(f"Reoriented {seq_id}: forward_score={forward_score}, reverse_score={reverse_score}")
            else:
                # Both zero (no primers) or both non-zero (ambiguous) - orientation failed
                failed_count += 1
                failed_sequences.add(seq_id)  # Track this sequence as failed
                if forward_score == 0 and reverse_score == 0:
                    logging.debug(f"No primer matches for {seq_id}")
                else:
                    logging.debug(f"Ambiguous orientation for {seq_id}: forward_score={forward_score}, reverse_score={reverse_score}")
        
        logging.info(f"Orientation complete: {already_correct} kept as-is, "
                    f"{oriented_count} reverse-complemented, {failed_count} orientation failed")
        
        # Return set of failed sequence IDs for potential filtering
        return failed_sequences
    
    def _score_orientation(self, sequence: str, orientation: str) -> int:
        """Score how well primers match in the given orientation.
        
        Simple binary scoring:
        - +1 if a forward primer is found at the expected position
        - +1 if a reverse primer is found at the expected position
        
        Args:
            sequence: The sequence to test
            orientation: Either "forward" or "reverse"
            
        Returns:
            Score from 0-2 (integer)
        """
        score = 0
        
        if orientation == "forward":
            # Forward orientation:
            # - Check for forward primers at 5' end (as-is)
            # - Check for RC of reverse primers at 3' end
            if self._has_primer_match(sequence, self.forward_primers, "start"):
                score += 1
            if self._has_primer_match(sequence, self.reverse_primers_rc, "end"):
                score += 1
        else:
            # Reverse orientation:
            # - Check for reverse primers at 5' end (as-is)
            # - Check for RC of forward primers at 3' end
            if self._has_primer_match(sequence, self.reverse_primers, "start"):
                score += 1
            if self._has_primer_match(sequence, self.forward_primers_rc, "end"):
                score += 1
        
        return score
    
    def _has_primer_match(self, sequence: str, primers: List[Tuple[str, str]], end: str) -> bool:
        """Check if any primer matches at the specified end of sequence.
        
        Args:
            sequence: The sequence to search in
            primers: List of (name, sequence) tuples to search for
            end: Either "start" or "end"
            
        Returns:
            True if any primer has a good match, False otherwise
        """
        if not primers or not sequence:
            return False
        
        # Determine search region
        max_primer_len = max(len(p[1]) for p in primers) if primers else 50
        if end == "start":
            search_region = sequence[:min(max_primer_len * 2, len(sequence))]
        else:
            search_region = sequence[-min(max_primer_len * 2, len(sequence)):]
        
        for primer_name, primer_seq in primers:
            # Allow up to 25% errors
            k = len(primer_seq) // 4
            
            # Use edlib to find best match
            result = edlib.align(primer_seq, search_region, task="distance", mode="HW", k=k)
            
            if result["editDistance"] != -1:
                # Consider it a match if identity is >= 75%
                identity = 1.0 - (result["editDistance"] / len(primer_seq))
                if identity >= 0.75:
                    logging.debug(f"Found {primer_name} at {end} with identity {identity:.2%} "
                                f"(edit_dist={result['editDistance']}, len={len(primer_seq)})")
                    return True
        
        return False
    
    def trim_primers(self, sequence: str) -> Tuple[str, List[str]]:
        """
        Trim primers from start and end of sequence.
        Returns (trimmed_sequence, [found_primers])
        """
        if not hasattr(self, 'primers') or not self.primers:
            return sequence, []

        found_primers = []
        trimmed_seq = sequence

        logging.debug(f"Starting primer trimming on sequence of length {len(sequence)}")

        # Look for primer at 5' end
        best_start_dist = float('inf')
        best_start_primer = None
        best_start_end = None

        for primer_name, primer_seq in self.primers:
            k = len(primer_seq) // 4  # Allow ~25% errors
            search_region = sequence[:len(primer_seq) * 2]

            result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)
            logging.debug(f"Testing {primer_name} at 5' end (k={k})")

            if result["editDistance"] != -1:
                dist = result["editDistance"]
                if dist < best_start_dist:
                    best_start_dist = dist
                    best_start_primer = primer_name
                    best_start_end = result["locations"][0][1] + 1

        if best_start_primer:
            logging.debug(f"Found {best_start_primer} at 5' end (k={best_start_dist})")
            found_primers.append(f"5'-{best_start_primer}")
            trimmed_seq = trimmed_seq[best_start_end:]

        # Look for primer at 3' end
        best_end_dist = float('inf')
        best_end_primer = None
        best_end_start = None

        for primer_name, primer_seq in self.primers:
            k = len(primer_seq) // 4  # Allow ~25% errors
            search_region = sequence[-len(primer_seq) * 2:]

            result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)

            if result["editDistance"] != -1:
                dist = result["editDistance"]
                if dist < best_end_dist:
                    best_end_dist = dist
                    best_end_primer = primer_name
                    base_pos = len(trimmed_seq) - len(search_region)
                    best_end_start = base_pos + result["locations"][0][0]

        if best_end_primer:
            logging.debug(f"Found {best_end_primer} at 3' end (k={best_end_dist})")
            found_primers.append(f"3'-{best_end_primer}")
            trimmed_seq = trimmed_seq[:best_end_start]

        return trimmed_seq, found_primers

    def calculate_consensus_distance(self, seq1: str, seq2: str, require_merge_compatible: bool = False) -> int:
        """Calculate distance between two consensus sequences using adjusted identity.
        
        Uses custom adjustment parameters that enable only homopolymer normalization:
        - Homopolymer differences (e.g., AAA vs AAAAA) are treated as identical
        - Regular substitutions count as mismatches
        - Non-homopolymer indels optionally prevent merging
        
        Args:
            seq1: First consensus sequence
            seq2: Second consensus sequence
            require_merge_compatible: If True, return -1 when sequences have variations
                                     that cannot be represented in IUPAC consensus (indels)
        
        Returns:
            Distance between sequences (substitutions only), or -1 if require_merge_compatible=True
            and sequences contain non-homopolymer indels
        """
        if not seq1 or not seq2:
            return max(len(seq1), len(seq2))

        # Get alignment from edlib (uses global NW alignment by default)
        result = edlib.align(seq1, seq2, task="path")
        if result["editDistance"] == -1:
            # Alignment failed, return maximum possible distance
            return max(len(seq1), len(seq2))
            
        # Get nice alignment for adjusted identity scoring
        alignment = edlib.getNiceAlignment(result, seq1, seq2)
        if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
            # Fall back to edit distance if alignment extraction fails
            return result["editDistance"]
            
        # Configure custom adjustment parameters for homopolymer normalization only
        custom_params = AdjustmentParams(
            normalize_homopolymers=True,    # Enable homopolymer normalization
            handle_iupac_overlap=False,     # Disable IUPAC overlap handling
            normalize_indels=False,         # Disable indel normalization
            end_skip_distance=0,            # Disable end trimming
            max_repeat_motif_length=2       # Keep default motif length
        )
        
        # Create custom scoring format to distinguish indels from substitutions
        from adjusted_identity import ScoringFormat
        custom_format = ScoringFormat(
            match='|',
            substitution='X',     # Distinct code for substitutions
            indel_start='I',      # Distinct code for indels
            indel_extension='-',
            homopolymer_extension='=',
            end_trimmed='.'
        )
        
        # Calculate adjusted identity with custom format
        score_result = score_alignment(
            alignment['query_aligned'], 
            alignment['target_aligned'],
            adjustment_params=custom_params,
            scoring_format=custom_format
        )
        
        # Check for merge compatibility if requested
        # Both non-homopolymer indels ('I') and terminal overhangs ('.') prevent merging
        if require_merge_compatible:
            if 'I' in score_result.score_aligned:
                logging.debug(f"Non-homopolymer indel detected, sequences not merge-compatible")
                return -1  # Signal that merging should not occur
            if '.' in score_result.score_aligned:
                logging.debug(f"Terminal overhang detected, sequences not merge-compatible")
                return -1  # Signal that merging should not occur
        
        # Count only substitutions (not homopolymer adjustments or indels)
        # Note: mismatches includes both substitutions and non-homopolymer indels
        # For accurate distance when indels are present, we use the mismatches count
        distance = score_result.mismatches
        
        # Log details about the variations found
        substitutions = score_result.score_aligned.count('X')
        indels = score_result.score_aligned.count('I')
        homopolymers = score_result.score_aligned.count('=')
        
        logging.debug(f"Consensus distance: {distance} total mismatches "
                     f"({substitutions} substitutions, {indels} indels, "
                     f"{homopolymers} homopolymer adjustments)")
        
        return distance

    def are_homopolymer_equivalent(self, seq1: str, seq2: str) -> bool:
        """Check if two sequences are equivalent when considering only homopolymer differences.

        Uses adjusted-identity scoring with global alignment. Terminal overhangs (marked as '.')
        and non-homopolymer indels (marked as 'I') prevent merging, ensuring truncated sequences
        don't merge with full-length sequences.
        """
        if not seq1 or not seq2:
            return seq1 == seq2

        # Use calculate_consensus_distance with merge compatibility check
        # Global alignment ensures terminal gaps are counted as indels
        # Returns: -1 (non-homopolymer indels), 0 (homopolymer-equivalent), >0 (substitutions)
        # Only distance == 0 means truly homopolymer-equivalent
        distance = self.calculate_consensus_distance(seq1, seq2, require_merge_compatible=True)
        return distance == 0

    def parse_mcl_output(self, mcl_output_file: str) -> List[Set[str]]:
        """Parse MCL output file into clusters of original sequence IDs."""
        clusters = []
        with open(mcl_output_file) as f:
            for line in f:
                # Each line is a tab-separated list of cluster members
                short_ids = line.strip().split('\t')
                # Map short IDs back to original sequence IDs
                cluster = {self.id_map[short_id] for short_id in short_ids}
                clusters.append(cluster)
        return clusters



def main():
    parser = argparse.ArgumentParser(
        description="MCL-based clustering of nanopore amplicon reads"
    )
    parser.add_argument("input_file", help="Input FASTQ file")
    parser.add_argument("--augment-input", help="Additional FASTQ/FASTA file with sequences recovered after primary demultiplexing (e.g., from specimine)")
    parser.add_argument("--algorithm", type=str, default="graph", choices=["graph", "greedy"],
                        help="Clustering algorithm to use (default: graph)")
    parser.add_argument("--min-identity", type=float, default=0.85,
                        help="Minimum sequence identity threshold (default: 0.85)")
    parser.add_argument("--inflation", type=float, default=4.0,
                        help="MCL inflation parameter (default: 4.0)")
    parser.add_argument("--min-size", type=int, default=5,
                        help="Minimum cluster size (default: 5, 0 to disable)")
    parser.add_argument("--min-cluster-ratio", type=float, default=0.2,
                        help="Minimum size ratio between a cluster and the largest cluster (default: 0.2, 0 to disable)")
    parser.add_argument("--max-sample-size", type=int, default=500,
                        help="Maximum cluster size for consensus (default: 500)")
    parser.add_argument("--mean-error-rate", type=float, default=None,
                        help="Estimated mean read-to-consensus error rate (e.g., 0.013 for 1.3%%). "
                             "When set, outlier reads exceeding threshold (mean × 3.0) will be removed "
                             "and consensus regenerated. Leave unset to disable outlier removal.")
    parser.add_argument("--presample", type=int, default=1000,
                        help="Presample size for initial reads (default: 1000, 0 to disable)")
    parser.add_argument("--k-nearest-neighbors", type=int, default=5,
                        help="Number of nearest neighbors for graph construction (default: 5)")
    parser.add_argument("--primers", help="FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)")
    parser.add_argument("-O", "--output-dir", default="clusters",
                        help="Output directory for all files (default: clusters)")
    parser.add_argument("--disable-homopolymer-equivalence", action="store_true",
                        help="Disable homopolymer equivalence in cluster merging (only merge identical sequences)")
    parser.add_argument("--orient-mode", choices=["skip", "keep-all", "filter-failed"], default="skip",
                        help="Sequence orientation mode: skip (default, no orientation), keep-all (orient but keep failed), or filter-failed (orient and remove failed)")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    parser.add_argument("--version", action="version",
                        version=f"Speconsense {__version__}",
                        help="Show program's version number and exit")

    args = parser.parse_args()

    # Setup standard logging
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format=log_format
    )

    sample = os.path.splitext(os.path.basename(args.input_file))[0]
    clusterer = SpecimenClusterer(
        min_identity=args.min_identity,
        inflation=args.inflation,
        min_size=args.min_size,
        min_cluster_ratio=args.min_cluster_ratio,
        max_sample_size=args.max_sample_size,
        presample_size=args.presample,
        k_nearest_neighbors=args.k_nearest_neighbors,
        sample_name=sample,
        disable_homopolymer_equivalence=args.disable_homopolymer_equivalence,
        output_dir=args.output_dir,
        mean_error_rate=args.mean_error_rate
    )

    # Log outlier removal configuration if enabled
    if args.mean_error_rate is not None:
        logging.info(f"Outlier removal enabled: mean_error_rate={args.mean_error_rate*100:.2f}%, "
                    f"threshold={clusterer.outlier_threshold*100:.2f}%")

    # Set additional attributes for metadata
    clusterer.input_file = os.path.abspath(args.input_file)
    clusterer.augment_input = os.path.abspath(args.augment_input) if args.augment_input else None
    clusterer.algorithm = args.algorithm
    clusterer.orient_mode = args.orient_mode

    # Read primary sequences
    logging.info(f"Reading sequences from {args.input_file}")
    format = "fasta" if args.input_file.endswith(".fasta") else "fastq"
    records = list(SeqIO.parse(args.input_file, format))
    logging.info(f"Loaded {len(records)} primary sequences")

    # Load augmented sequences if specified
    augment_records = None
    if args.augment_input:
        # Check if augment input file exists
        if not os.path.exists(args.augment_input):
            logging.error(f"Augment input file not found: {args.augment_input}")
            sys.exit(1)
            
        logging.info(f"Reading augmented sequences from {args.augment_input}")
        
        # Auto-detect format like main input
        augment_format = "fasta" if args.augment_input.endswith(".fasta") else "fastq"
        
        try:
            augment_records = list(SeqIO.parse(args.augment_input, augment_format))
            logging.info(f"Loaded {len(augment_records)} augmented sequences")
            
            if len(augment_records) == 0:
                logging.warning(f"No sequences found in augment input file: {args.augment_input}")
            
            # Add dummy quality scores to FASTA sequences so they can be written as FASTQ later
            if augment_format == "fasta":
                for record in augment_records:
                    if not hasattr(record, 'letter_annotations') or 'phred_quality' not in record.letter_annotations:
                        # Add dummy quality scores (quality 30 = '?' in FASTQ)
                        record.letter_annotations = {'phred_quality': [30] * len(record.seq)}
                logging.debug(f"Added quality scores to {len(augment_records)} FASTA sequences for downstream compatibility")
                
        except Exception as e:
            logging.error(f"Failed to read augment input file '{args.augment_input}': {e}")
            sys.exit(1)

    # Add sequences to clusterer (both primary and augmented)
    clusterer.add_sequences(records, augment_records)

    if args.primers:
        clusterer.primers_file = os.path.abspath(args.primers)
        clusterer.load_primers(args.primers)
    else:
        # Look for primers.fasta in the same directory as the input file
        input_dir = os.path.dirname(os.path.abspath(args.input_file))
        auto_primer_path = os.path.join(input_dir, "primers.fasta")

        if os.path.exists(auto_primer_path):
            logging.info(f"Found primers.fasta in input directory: {auto_primer_path}")
            clusterer.primers_file = os.path.abspath(auto_primer_path)
            clusterer.load_primers(auto_primer_path)
        else:
            logging.warning("No primer file specified and primers.fasta not found in input directory. Primer trimming will be disabled.")
            clusterer.primers_file = None
    
    # Handle sequence orientation based on mode
    if args.orient_mode != "skip":
        if hasattr(clusterer, 'forward_primers') and hasattr(clusterer, 'reverse_primers'):
            failed_sequences = clusterer.orient_sequences()

            # Filter failed sequences if requested
            if args.orient_mode == "filter-failed" and failed_sequences:
                logging.info(f"Filtering out {len(failed_sequences)} sequences with failed orientation")

                # Remove failed sequences from clusterer
                for seq_id in failed_sequences:
                    del clusterer.sequences[seq_id]
                    del clusterer.records[seq_id]

                remaining = len(clusterer.sequences)
                logging.info(f"Continuing with {remaining} successfully oriented sequences")
        else:
            logging.warning(f"--orient-mode={args.orient_mode} specified but no primers with position information loaded")

    # Write metadata file for use by post-processing tools
    clusterer.write_metadata()

    clusterer.cluster(algorithm=args.algorithm)
    print()

if __name__ == "__main__":
    main()

