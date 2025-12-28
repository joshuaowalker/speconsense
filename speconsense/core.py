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
from typing import List, Set, Tuple, Optional, Dict, Any, NamedTuple
from datetime import datetime

import edlib
from adjusted_identity import score_alignment, AdjustmentParams
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

try:
    from speconsense import __version__
except ImportError:
    # Fallback for when running as a script directly (e.g., in tests)
    __version__ = "dev"

# Import MSA analysis functions and data structures from msa module
from speconsense.msa import (
    # Data structures
    ErrorPosition,
    ReadAlignment,
    PositionStats,
    MSAResult,
    IUPAC_CODES,
    # Functions
    parse_score_aligned_for_errors,
    extract_alignments_from_msa,
    analyze_positional_variation,
    is_variant_position_with_composition,
    call_iupac_ambiguities,
    group_reads_by_single_position,
    filter_qualifying_haplotypes,
    calculate_within_cluster_error,
)
from speconsense.scalability import (
    VsearchCandidateFinder,
    ScalablePairwiseOperation,
    ScalabilityConfig,
)


class SpecimenClusterer:
    def __init__(self, min_identity: float = 0.9,
                 inflation: float = 4.0,
                 min_size: int = 5,
                 min_cluster_ratio: float = 0.2,
                 max_sample_size: int = 500,
                 presample_size: int = 1000,
                 k_nearest_neighbors: int = 20,
                 sample_name: str = "sample",
                 disable_homopolymer_equivalence: bool = False,
                 disable_cluster_merging: bool = False,
                 output_dir: str = "clusters",
                 outlier_identity_threshold: Optional[float] = None,
                 enable_secondpass_phasing: bool = True,
                 min_variant_frequency: float = 0.20,
                 min_variant_count: int = 5,
                 min_ambiguity_frequency: float = 0.10,
                 min_ambiguity_count: int = 3,
                 enable_iupac_calling: bool = True,
                 enable_scalability: Optional[int] = None,
                 max_threads: int = 1,
                 early_filter: bool = True,
                 collect_discards: bool = False):
        self.min_identity = min_identity
        self.inflation = inflation
        self.min_size = min_size
        self.min_cluster_ratio = min_cluster_ratio
        self.max_sample_size = max_sample_size
        self.presample_size = presample_size
        self.k_nearest_neighbors = k_nearest_neighbors
        self.sample_name = sample_name
        self.disable_homopolymer_equivalence = disable_homopolymer_equivalence
        self.disable_cluster_merging = disable_cluster_merging
        self.output_dir = output_dir

        # Auto-calculate outlier identity threshold if not provided
        # Logic: min_identity accounts for 2×error (read-to-read comparison)
        # outlier_identity_threshold accounts for 1×error (read-to-consensus)
        # Therefore: outlier_identity_threshold = (1 + min_identity) / 2
        if outlier_identity_threshold is None:
            self.outlier_identity_threshold = (1.0 + min_identity) / 2.0
        else:
            self.outlier_identity_threshold = outlier_identity_threshold

        self.enable_secondpass_phasing = enable_secondpass_phasing
        self.min_variant_frequency = min_variant_frequency
        self.min_variant_count = min_variant_count
        self.min_ambiguity_frequency = min_ambiguity_frequency
        self.min_ambiguity_count = min_ambiguity_count
        self.enable_iupac_calling = enable_iupac_calling
        self.enable_scalability = enable_scalability
        self.max_threads = max_threads
        self.early_filter = early_filter
        self.collect_discards = collect_discards
        self.discarded_read_ids: Set[str] = set()  # Track all discarded reads (outliers + filtered)

        # Initialize scalability configuration
        # enable_scalability: None=disabled, 0=always, N=threshold
        self.scalability_config = ScalabilityConfig(
            enabled=enable_scalability is not None,
            activation_threshold=enable_scalability if enable_scalability is not None else 0,
            max_threads=max_threads
        )
        self._candidate_finder = None
        if enable_scalability is not None:
            self._candidate_finder = VsearchCandidateFinder(num_threads=max_threads)
            if not self._candidate_finder.is_available:
                logging.warning("Scalability enabled but vsearch not found. Falling back to brute-force.")
                self._candidate_finder = None

        self.sequences = {}  # id -> sequence string
        self.records = {}  # id -> SeqRecord object
        self.id_map = {}  # short_id -> original_id
        self.rev_id_map = {}  # original_id -> short_id

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
                "disable_cluster_merging": self.disable_cluster_merging,
                "outlier_identity_threshold": self.outlier_identity_threshold,
                "enable_secondpass_phasing": self.enable_secondpass_phasing,
                "min_variant_frequency": self.min_variant_frequency,
                "min_variant_count": self.min_variant_count,
                "min_ambiguity_frequency": self.min_ambiguity_frequency,
                "min_ambiguity_count": self.min_ambiguity_count,
                "enable_iupac_calling": self.enable_iupac_calling,
                "enable_scalability": self.enable_scalability,
                "max_threads": self.max_threads,
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

        logging.debug(f"Wrote run metadata to {metadata_file}")

    def write_phasing_stats(self, initial_clusters_count: int, after_prephasing_merge_count: int,
                           subclusters_count: int, merged_count: int, final_count: int,
                           clusters_with_ambiguities: int = 0,
                           total_ambiguity_positions: int = 0) -> None:
        """Write phasing statistics to JSON file after clustering completes.

        Args:
            initial_clusters_count: Number of clusters from initial clustering
            after_prephasing_merge_count: Number of clusters after pre-phasing merge
            subclusters_count: Number of sub-clusters after phasing
            merged_count: Number of clusters after post-phasing merge
            final_count: Number of final clusters after filtering
            clusters_with_ambiguities: Number of clusters with at least one ambiguity code
            total_ambiguity_positions: Total number of ambiguity positions across all clusters
        """
        phasing_stats = {
            "phasing_enabled": self.enable_secondpass_phasing,
            "initial_clusters": initial_clusters_count,
            "after_prephasing_merge": after_prephasing_merge_count,
            "phased_subclusters": subclusters_count,
            "after_postphasing_merge": merged_count,
            "after_filtering": final_count,
            "prephasing_clusters_merged": after_prephasing_merge_count < initial_clusters_count,
            "clusters_split": subclusters_count > after_prephasing_merge_count,
            "postphasing_clusters_merged": merged_count < subclusters_count,
            "net_change": final_count - initial_clusters_count,
            "ambiguity_calling_enabled": self.enable_iupac_calling,
            "clusters_with_ambiguities": clusters_with_ambiguities,
            "total_ambiguity_positions": total_ambiguity_positions
        }

        # Write phasing stats to separate JSON file
        stats_file = os.path.join(self.debug_dir, f"{self.sample_name}-phasing_stats.json")
        with open(stats_file, 'w') as f:
            json.dump(phasing_stats, f, indent=2)

        logging.debug(f"Wrote phasing statistics to {stats_file}")

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

        # Log recommendation for large datasets
        if len(self.sequences) > self.scalability_config.recommendation_threshold and self.enable_scalability is None:
            logging.info(f"Large dataset detected ({len(self.sequences)} sequences). "
                         "Consider using --enable-scalability for faster processing.")

    def _get_scalable_operation(self) -> ScalablePairwiseOperation:
        """Get configured scalable operation for K-NN computation."""
        return ScalablePairwiseOperation(
            candidate_finder=self._candidate_finder,
            scoring_function=self.calculate_similarity,
            config=self.scalability_config
        )

    def write_mcl_input(self, output_file: str) -> None:
        """Write similarity matrix in MCL input format using k-nearest neighbors approach."""
        self._create_id_mapping()

        n = len(self.sequences)
        k = min(self.k_nearest_neighbors, n - 1)  # Connect to at most k neighbors

        # Use scalable operation to compute K-NN edges
        operation = self._get_scalable_operation()
        knn_edges = operation.compute_top_k_neighbors(
            sequences=self.sequences,
            k=k,
            min_identity=self.min_identity,
            output_dir=self.output_dir,
            min_edges_per_node=3
        )

        # Write edges to MCL input file
        with open(output_file, 'w') as f:
            for id1, neighbors in knn_edges.items():
                short_id1 = self.rev_id_map[id1]
                for id2, sim in neighbors:
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
            "-te", str(self.max_threads),  # Number of threads
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

    def merge_similar_clusters(self, clusters: List[Dict], phase_name: str = "Post-phasing") -> List[Dict]:
        """
        Merge clusters whose consensus sequences are identical or homopolymer-equivalent.
        Preserves provenance metadata through the merging process.

        This function is used for both pre-phasing merge (combining initial clusters before
        variant detection) and post-phasing merge (combining subclusters after phasing).

        Note: Primer trimming is performed before comparison to ensure clusters that differ
        only in primer regions are properly merged. Trimmed consensuses are used only for
        comparison and are discarded after merging.

        Args:
            clusters: List of cluster dictionaries with 'read_ids' and provenance fields
            phase_name: Name of the merge phase for logging (e.g., "Pre-phasing", "Post-phasing")

        Returns:
            List of merged cluster dictionaries with combined provenance
        """
        if not clusters:
            return []

        # Sort clusters by size, largest first
        clusters = sorted(clusters, key=lambda c: len(c['read_ids']), reverse=True)

        # Generate a consensus sequence for each cluster
        logging.debug(f"{phase_name} merge: Generating consensus sequences...")
        consensuses = []
        cluster_to_consensus = {}  # Map from cluster index to its consensus

        # First pass: prepare sampled sequences and handle single-read clusters
        clusters_needing_spoa = []  # (cluster_index, sampled_seqs)

        for i, cluster_dict in enumerate(clusters):
            cluster_reads = cluster_dict['read_ids']

            # Skip empty clusters
            if not cluster_reads:
                logging.warning(f"Cluster {i} is empty, skipping")
                continue

            # Single-read clusters don't need SPOA - use the read directly
            if len(cluster_reads) == 1:
                seq_id = next(iter(cluster_reads))
                consensus = self.sequences[seq_id]
                # Trim primers before comparison
                if hasattr(self, 'primers'):
                    consensus, _ = self.trim_primers(consensus)
                consensuses.append(consensus)
                cluster_to_consensus[i] = consensus
                continue

            # Sample from larger clusters to speed up consensus generation
            if len(cluster_reads) > self.max_sample_size:
                # Sample by quality
                qualities = []
                for seq_id in cluster_reads:
                    record = self.records[seq_id]
                    mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                    qualities.append((mean_quality, seq_id))

                # Sort by quality (descending), then by read ID (ascending) for deterministic tie-breaking
                sampled_ids = [seq_id for _, seq_id in
                               sorted(qualities, key=lambda x: (-x[0], x[1]))[:self.max_sample_size]]
                sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sampled_ids}
            else:
                # Sort all reads by quality for optimal SPOA ordering
                qualities = []
                for seq_id in cluster_reads:
                    record = self.records[seq_id]
                    mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                    qualities.append((mean_quality, seq_id))
                sorted_ids = [seq_id for _, seq_id in
                              sorted(qualities, key=lambda x: (-x[0], x[1]))]
                sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sorted_ids}

            clusters_needing_spoa.append((i, sampled_seqs))

        # Run SPOA for multi-read clusters
        if clusters_needing_spoa:
            if self.max_threads > 1 and len(clusters_needing_spoa) > 10:
                # Parallel SPOA execution
                from concurrent.futures import ThreadPoolExecutor

                def run_spoa_for_cluster(args):
                    cluster_idx, sampled_seqs = args
                    return cluster_idx, self.run_spoa(sampled_seqs)

                with ThreadPoolExecutor(max_workers=self.max_threads) as executor:
                    from tqdm import tqdm
                    results = list(tqdm(
                        executor.map(run_spoa_for_cluster, clusters_needing_spoa),
                        total=len(clusters_needing_spoa),
                        desc=f"{phase_name} consensus generation"
                    ))

                for cluster_idx, result in results:
                    if result is None:
                        logging.warning(f"Cluster {cluster_idx} produced no consensus, skipping")
                        continue
                    consensus = result.consensus
                    if hasattr(self, 'primers'):
                        consensus, _ = self.trim_primers(consensus)
                    consensuses.append(consensus)
                    cluster_to_consensus[cluster_idx] = consensus
            else:
                # Sequential SPOA execution
                for cluster_idx, sampled_seqs in clusters_needing_spoa:
                    result = self.run_spoa(sampled_seqs)
                    if result is None:
                        logging.warning(f"Cluster {cluster_idx} produced no consensus, skipping")
                        continue
                    consensus = result.consensus
                    if hasattr(self, 'primers'):
                        consensus, _ = self.trim_primers(consensus)
                    consensuses.append(consensus)
                    cluster_to_consensus[cluster_idx] = consensus

        consensus_to_clusters = defaultdict(list)

        if self.disable_homopolymer_equivalence:
            # Only merge exactly identical sequences
            for i, consensus in enumerate(consensuses):
                consensus_to_clusters[consensus].append(i)
        else:
            # Group by homopolymer-equivalent sequences
            # Use scalable method when enabled and there are many clusters
            use_scalable = (
                self.enable_scalability is not None and
                self._candidate_finder is not None and
                self._candidate_finder.is_available and
                len(cluster_to_consensus) > 50
            )

            if use_scalable:
                # Map cluster indices to string IDs for scalability module
                str_to_index = {str(i): i for i in cluster_to_consensus.keys()}
                consensus_seq_dict = {str(i): seq for i, seq in cluster_to_consensus.items()}

                # Use scalable equivalence grouping
                operation = self._get_scalable_operation()
                equivalence_groups = operation.compute_equivalence_groups(
                    sequences=consensus_seq_dict,
                    equivalence_fn=self.are_homopolymer_equivalent,
                    output_dir=self.output_dir,
                    min_candidate_identity=0.95
                )

                # Convert groups back to indices and populate consensus_to_clusters
                for group in equivalence_groups:
                    if group:
                        representative = group[0]
                        repr_consensus = consensus_seq_dict[representative]
                        for str_id in group:
                            consensus_to_clusters[repr_consensus].append(str_to_index[str_id])
            else:
                # Original O(n²) approach for small sets
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

        # Determine merge type for logging
        merge_type = "identical" if self.disable_homopolymer_equivalence else "homopolymer-equivalent"

        # Handle clusters with equivalent consensus sequences
        for equivalent_clusters in consensus_to_clusters.values():
            if len(equivalent_clusters) > 1:
                # Merge clusters with equivalent consensus
                merged_read_ids = set()
                merged_from_list = []

                # Check if we're merging phased subclusters from the same initial cluster
                initial_clusters_involved = set()
                phased_subclusters_merged = []

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

                    # Track if phased subclusters are being merged
                    if clusters[idx].get('allele_combo') is not None:
                        phased_subclusters_merged.append(cluster_info)
                        initial_clusters_involved.add(clusters[idx]['initial_cluster_num'])

                # Log if we're merging phased subclusters that came from the same initial cluster
                # This can happen when SPOA consensus generation doesn't preserve variant differences
                # that were detected during phasing (e.g., due to homopolymer normalization differences)
                if len(phased_subclusters_merged) > 1 and len(initial_clusters_involved) == 1:
                    initial_cluster = list(initial_clusters_involved)[0]
                    logging.debug(
                        f"Merging {len(phased_subclusters_merged)} phased subclusters from initial cluster {initial_cluster} "
                        f"back together (consensus sequences are {merge_type})"
                    )
                    for info in phased_subclusters_merged:
                        logging.debug(f"  Subcluster: allele_combo='{info['allele_combo']}', size={info['size']}")

                # Create merged cluster with provenance
                merged_cluster = {
                    'read_ids': merged_read_ids,
                    'initial_cluster_num': None,  # Multiple sources
                    'allele_combo': None,  # Multiple alleles merged
                    'merged_from': merged_from_list  # Track merge provenance
                }
                merged.append(merged_cluster)

        # Add remaining unmerged clusters
        for i, cluster_dict in enumerate(clusters):
            if i not in merged_indices:
                merged.append(cluster_dict)

        if len(merged) < len(clusters):
            logging.info(f"{phase_name} merge: Combined {len(clusters)} clusters into {len(merged)} "
                        f"({len(clusters) - len(merged)} merged due to {merge_type} consensus)")
        else:
            logging.info(f"{phase_name} merge: No clusters merged (no {merge_type} consensus found)")

        return merged

    def _find_root(self, merged_to: List[int], i: int) -> int:
        """Find the root index of a merged cluster using path compression."""
        if merged_to[i] != i:
            merged_to[i] = self._find_root(merged_to, merged_to[i])
        return merged_to[i]

    def write_cluster_files(self, cluster_num: int, cluster: Set[str],
                            consensus: str, trimmed_consensus: Optional[str] = None,
                            found_primers: Optional[List[str]] = None,
                            rid: Optional[float] = None,
                            rid_min: Optional[float] = None,
                            actual_size: Optional[int] = None,
                            consensus_fasta_handle = None,
                            sampled_ids: Optional[Set[str]] = None,
                            msa: Optional[str] = None,
                            sorted_cluster_ids: Optional[List[str]] = None,
                            sorted_sampled_ids: Optional[List[str]] = None,
                            iupac_count: int = 0) -> None:
        """Write cluster files: reads FASTQ, MSA, and consensus FASTA.

        Read identity metrics measure internal cluster consistency (not accuracy vs. ground truth):
        - rid: Mean read identity - measures average agreement between reads and consensus
        - rid_min: Minimum read identity - captures worst-case outlier reads

        High identity values indicate homogeneous clusters with consistent reads.
        Low values may indicate heterogeneity, outliers, or poor consensus (especially at low RiC).
        """
        cluster_size = len(cluster)
        ric_size = min(actual_size or cluster_size, self.max_sample_size)

        # Create info string with size first
        info_parts = [f"size={cluster_size}", f"ric={ric_size}"]

        # Add read identity metrics (as percentages for readability)
        if rid is not None:
            info_parts.append(f"rid={rid*100:.1f}")
        if rid_min is not None:
            info_parts.append(f"rid_min={rid_min*100:.1f}")

        if found_primers:
            info_parts.append(f"primers={','.join(found_primers)}")
        if iupac_count > 0:
            info_parts.append(f"ambig={iupac_count}")
        info_str = " ".join(info_parts)

        # Write reads FASTQ to debug directory with new naming convention
        # Use sorted order (by quality descending) if available, matching MSA order
        reads_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-reads.fastq")
        with open(reads_file, 'w') as f:
            read_ids_to_write = sorted_cluster_ids if sorted_cluster_ids is not None else cluster
            for seq_id in read_ids_to_write:
                SeqIO.write(self.records[seq_id], f, "fastq")

        # Write sampled reads FASTQ (only sequences used for consensus generation)
        # Use sorted order (by quality descending) if available, matching MSA order
        if sampled_ids is not None or sorted_sampled_ids is not None:
            sampled_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-sampled.fastq")
            with open(sampled_file, 'w') as f:
                sampled_to_write = sorted_sampled_ids if sorted_sampled_ids is not None else sampled_ids
                for seq_id in sampled_to_write:
                    SeqIO.write(self.records[seq_id], f, "fastq")

        # Write MSA (multiple sequence alignment) to debug directory
        if msa is not None:
            msa_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-msa.fasta")
            with open(msa_file, 'w') as f:
                f.write(msa)

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

        # Sort for deterministic order
        seq_ids = sorted(self.sequences.keys())
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

        # Sort for deterministic iteration (important for tie-breaking)
        for seq_id in sorted(available_ids):
            # Get all sequences that align with current sequence
            members = {other_id for other_id in self.alignments.get(seq_id, {})
                       if other_id in available_ids}

            # Use > (not >=) so first alphabetically wins ties
            if len(members) > best_count:
                best_count = len(members)
                best_center = seq_id
                best_members = members

        if best_center is None:
            # No alignments found, create singleton cluster with smallest ID
            singleton_id = min(available_ids)
            return singleton_id, {singleton_id}

        best_members.add(best_center)  # Include center in cluster
        return best_center, best_members


    # ========================================================================
    # Clustering Phase Helper Methods
    # ========================================================================

    def _run_initial_clustering(self, temp_dir: str, algorithm: str) -> List[Set[str]]:
        """Phase 1: Run initial clustering algorithm.
        
        Args:
            temp_dir: Temporary directory for intermediate files
            algorithm: 'graph' for MCL or 'greedy' for greedy clustering
            
        Returns:
            List of clusters (sets of read IDs), sorted by size (largest first)
        """
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
        return initial_clusters

    def _run_prephasing_merge(self, initial_clusters: List[Set[str]]) -> List[Set[str]]:
        """Phase 2: Merge initial clusters with HP-equivalent consensus.
        
        Maximizes read depth for variant detection in the phasing phase.
        
        Args:
            initial_clusters: List of initial clusters from Phase 1
            
        Returns:
            List of merged clusters (sets of read IDs)
        """
        if self.disable_cluster_merging:
            logging.info("Cluster merging disabled, skipping pre-phasing merge")
            return initial_clusters
        
        # Convert initial clusters to dict format for merge_similar_clusters
        initial_cluster_dicts = [
            {'read_ids': cluster, 'initial_cluster_num': i, 'allele_combo': None}
            for i, cluster in enumerate(initial_clusters, 1)
        ]
        merged_dicts = self.merge_similar_clusters(initial_cluster_dicts, phase_name="Pre-phasing")
        # Extract back to sets for Phase 3
        return [d['read_ids'] for d in merged_dicts]

    def _apply_early_filter(self, clusters: List[Set[str]]) -> Tuple[List[Set[str]], List[Set[str]]]:
        """Apply early size filtering after pre-phasing merge.

        Uses the same logic as _run_size_filtering() but operates before
        variant phasing to avoid expensive processing of small clusters.

        Args:
            clusters: List of merged clusters from Phase 2

        Returns:
            Tuple of (clusters_to_process, filtered_clusters)
        """
        if not self.early_filter:
            return clusters, []

        # Get size of each cluster for filtering
        cluster_sizes = [(c, len(c)) for c in clusters]

        # Find largest cluster size for ratio filtering
        if not cluster_sizes:
            return [], []
        largest_size = max(size for _, size in cluster_sizes)

        keep_clusters = []
        filtered_clusters = []

        for cluster, size in cluster_sizes:
            # Apply min_size filter
            if size < self.min_size:
                filtered_clusters.append(cluster)
                continue

            # Apply min_cluster_ratio filter
            if self.min_cluster_ratio > 0 and size / largest_size < self.min_cluster_ratio:
                filtered_clusters.append(cluster)
                continue

            keep_clusters.append(cluster)

        if filtered_clusters:
            # Collect discarded read IDs
            discarded_count = 0
            for cluster in filtered_clusters:
                self.discarded_read_ids.update(cluster)
                discarded_count += len(cluster)

            logging.info(f"Early filter: {len(filtered_clusters)} clusters ({discarded_count} reads) "
                        f"below threshold, {len(keep_clusters)} clusters proceeding to phasing")

        return keep_clusters, filtered_clusters

    def _process_single_cluster(self, initial_idx: int, cluster: Set[str]) -> Tuple[List[Dict], Set[str]]:
        """Process a single cluster for outlier removal and phasing.

        This method is designed to be thread-safe for parallel execution.
        It does not mutate any shared state directly.

        Args:
            initial_idx: The 1-based index of this cluster
            cluster: Set of read IDs in this cluster

        Returns:
            Tuple of (subclusters, discarded_read_ids) where:
            - subclusters: List of subcluster dicts with read_ids, initial_cluster_num, allele_combo
            - discarded_read_ids: Set of outlier read IDs that were removed
        """
        subclusters = []
        discarded_ids = set()

        # Sort all reads by quality for optimal SPOA ordering
        qualities = []
        for seq_id in cluster:
            record = self.records[seq_id]
            mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
            qualities.append((mean_quality, seq_id))
        sorted_ids = [seq_id for _, seq_id in
                      sorted(qualities, key=lambda x: (-x[0], x[1]))]
        cluster_ids = cluster

        # Generate consensus and MSA
        cluster_seqs = {seq_id: self.sequences[seq_id] for seq_id in sorted_ids}
        result = self.run_spoa(cluster_seqs)

        if result is None:
            logging.warning(f"Initial cluster {initial_idx}: Failed to generate consensus, skipping")
            return subclusters, discarded_ids

        consensus = result.consensus
        msa = result.msa_string
        alignments = result.alignments
        msa_to_consensus_pos = result.msa_to_consensus_pos

        # Calculate per-read identity from MSA
        _, _ = self.calculate_read_identity(alignments, consensus)

        # Optional: Remove outlier reads and regenerate consensus
        if self.outlier_identity_threshold is not None:
            keep_ids, outlier_ids = self.identify_outlier_reads(
                alignments, consensus, cluster_ids, self.outlier_identity_threshold
            )

            if outlier_ids:
                # Special case: 2-read cluster with exactly 1 outlier
                if len(cluster) == 2 and len(outlier_ids) == 1:
                    logging.info(f"Initial cluster {initial_idx}: Split 2-read cluster due to 1 outlier")
                    for read_id in cluster:
                        subcluster = {
                            'read_ids': {read_id},
                            'initial_cluster_num': initial_idx,
                            'allele_combo': 'single-read-split'
                        }
                        subclusters.append(subcluster)
                    return subclusters, discarded_ids

                logging.info(f"Initial cluster {initial_idx}: Removing {len(outlier_ids)}/{len(cluster_ids)} outlier reads, "
                           f"regenerating consensus")

                # Track discarded outlier reads (returned, not mutated)
                discarded_ids.update(outlier_ids)

                # Update cluster_ids to exclude outliers
                cluster_ids = keep_ids
                cluster = cluster - outlier_ids

                # Regenerate consensus with filtered reads
                qualities_filtered = []
                for seq_id in cluster_ids:
                    record = self.records[seq_id]
                    mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                    qualities_filtered.append((mean_quality, seq_id))
                sorted_ids_filtered = [seq_id for _, seq_id in
                                      sorted(qualities_filtered, key=lambda x: (-x[0], x[1]))]

                cluster_seqs = {seq_id: self.sequences[seq_id] for seq_id in sorted_ids_filtered}
                result = self.run_spoa(cluster_seqs)

                if result is not None:
                    consensus = result.consensus
                    msa = result.msa_string
                    alignments = result.alignments
                    msa_to_consensus_pos = result.msa_to_consensus_pos
                    _, _ = self.calculate_read_identity(alignments, consensus)

        # Detect variant positions (if second-pass phasing enabled)
        variant_positions = []
        if consensus and alignments and self.enable_secondpass_phasing:
            variant_positions = self.detect_variant_positions(
                alignments, consensus, msa_to_consensus_pos
            )

            if variant_positions:
                logging.info(f"Initial cluster {initial_idx}: Detected {len(variant_positions)} variant positions")

        # Phase reads into haplotypes
        phased_haplotypes = self.phase_reads_by_variants(
            msa, consensus, cluster, variant_positions, alignments
        )

        # Store each haplotype as a sub-cluster with provenance
        for haplotype_idx, (allele_combo, haplotype_reads) in enumerate(phased_haplotypes):
            subcluster = {
                'read_ids': haplotype_reads,
                'initial_cluster_num': initial_idx,
                'allele_combo': allele_combo
            }
            subclusters.append(subcluster)

        return subclusters, discarded_ids

    def _run_variant_phasing(self, merged_clusters: List[Set[str]]) -> List[Dict]:
        """Phase 3: Detect variants and phase reads into haplotypes.

        For each merged cluster:
        1. Sample reads if needed
        2. Generate consensus and MSA
        3. Optionally remove outliers
        4. Detect variant positions
        5. Phase reads by their alleles at variant positions

        Args:
            merged_clusters: List of merged clusters from Phase 2

        Returns:
            List of subclusters with provenance info (dicts with read_ids,
            initial_cluster_num, allele_combo)
        """
        all_subclusters = []
        all_discarded = set()
        logging.debug("Processing clusters for variant detection and phasing...")

        indexed_clusters = list(enumerate(merged_clusters, 1))

        if self.max_threads > 1 and len(merged_clusters) > 10:
            # Parallel processing
            from concurrent.futures import ThreadPoolExecutor

            def process_cluster(args):
                initial_idx, cluster = args
                return self._process_single_cluster(initial_idx, cluster)

            with ThreadPoolExecutor(max_workers=self.max_threads) as executor:
                from tqdm import tqdm
                results = list(tqdm(
                    executor.map(process_cluster, indexed_clusters),
                    total=len(indexed_clusters),
                    desc="Processing clusters"
                ))

            # Collect results maintaining order
            for subclusters, discarded_ids in results:
                all_subclusters.extend(subclusters)
                all_discarded.update(discarded_ids)
        else:
            # Sequential processing
            for initial_idx, cluster in indexed_clusters:
                subclusters, discarded_ids = self._process_single_cluster(initial_idx, cluster)
                all_subclusters.extend(subclusters)
                all_discarded.update(discarded_ids)

        # Update shared state after all processing complete
        self.discarded_read_ids.update(all_discarded)

        logging.info(f"After phasing, created {len(all_subclusters)} sub-clusters from {len(merged_clusters)} merged clusters")
        return all_subclusters

    def _run_postphasing_merge(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4: Merge subclusters with HP-equivalent consensus.
        
        Args:
            subclusters: List of subclusters from Phase 3
            
        Returns:
            List of merged subclusters
        """
        if self.disable_cluster_merging:
            logging.info("Cluster merging disabled, skipping post-phasing merge")
            return subclusters
        
        return self.merge_similar_clusters(subclusters, phase_name="Post-phasing")

    def _run_size_filtering(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 5: Filter clusters by size and ratio thresholds.
        
        Args:
            subclusters: List of subclusters from Phase 4
            
        Returns:
            List of filtered clusters, sorted by size (largest first)
        """
        # Filter by absolute size
        large_clusters = [c for c in subclusters if len(c['read_ids']) >= self.min_size]

        if len(large_clusters) < len(subclusters):
            filtered_count = len(subclusters) - len(large_clusters)
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

        if total_sequences > 0:
            logging.info(f"Final: {len(large_clusters)} clusters covering {sequences_covered} sequences "
                        f"({sequences_covered / total_sequences:.1%} of total)")
        else:
            logging.info(f"Final: {len(large_clusters)} clusters (no sequences to cluster)")
        
        return large_clusters

    def _write_cluster_outputs(self, clusters: List[Dict], output_file: str) -> Tuple[int, int]:
        """Phase 6: Generate final consensus and write output files.
        
        Args:
            clusters: List of filtered clusters from Phase 5
            output_file: Path to the output FASTA file
            
        Returns:
            Tuple of (clusters_with_ambiguities, total_ambiguity_positions)
        """
        total_ambiguity_positions = 0
        clusters_with_ambiguities = 0

        with open(output_file, 'w') as consensus_fasta_handle:
            for final_idx, cluster_dict in enumerate(clusters, 1):
                cluster = cluster_dict['read_ids']
                actual_size = len(cluster)

                # Sort all cluster reads by quality for consistent output ordering
                qualities = []
                for seq_id in cluster:
                    record = self.records[seq_id]
                    mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                    qualities.append((mean_quality, seq_id))
                # Sort by quality (descending), then by read ID (ascending) for deterministic tie-breaking
                sorted_cluster_ids = [seq_id for _, seq_id in
                                      sorted(qualities, key=lambda x: (-x[0], x[1]))]

                # Sample sequences for final consensus generation if needed
                if len(cluster) > self.max_sample_size:
                    logging.info(f"Cluster {final_idx}: Sampling {self.max_sample_size} from {len(cluster)} reads for final consensus")
                    sorted_sampled_ids = sorted_cluster_ids[:self.max_sample_size]
                    sampled_ids = set(sorted_sampled_ids)  # Keep set for membership testing
                else:
                    # All reads used for consensus
                    sorted_sampled_ids = sorted_cluster_ids
                    sampled_ids = cluster

                # Generate final consensus and MSA
                # Pass sequences to SPOA in quality-descending order
                sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sorted_sampled_ids}
                result = self.run_spoa(sampled_seqs)

                # Calculate final identity metrics
                rid, rid_min = None, None
                consensus = None
                msa = None

                if result is not None:
                    consensus = result.consensus
                    msa = result.msa_string
                    alignments = result.alignments
                    rid, rid_min = self.calculate_read_identity(alignments, consensus)

                # Note: No outlier removal or variant detection in output phase
                # (already done in detection phase)

                if consensus:
                    # Apply IUPAC ambiguity calling for unphased variant positions
                    iupac_count = 0
                    if self.enable_iupac_calling and result is not None:
                        consensus, iupac_count, iupac_details = call_iupac_ambiguities(
                            consensus=consensus,
                            alignments=result.alignments,
                            msa_to_consensus_pos=result.msa_to_consensus_pos,
                            min_variant_frequency=self.min_ambiguity_frequency,
                            min_variant_count=self.min_ambiguity_count
                        )
                        if iupac_count > 0:
                            logging.debug(f"Cluster {final_idx}: Called {iupac_count} IUPAC ambiguity position(s)")
                            total_ambiguity_positions += iupac_count
                            clusters_with_ambiguities += 1

                    # Perform primer trimming (on potentially modified consensus)
                    trimmed_consensus = None
                    found_primers = None
                    if hasattr(self, 'primers'):
                        trimmed_consensus, found_primers = self.trim_primers(consensus)

                    # Write output files
                    self.write_cluster_files(
                        cluster_num=final_idx,
                        cluster=cluster,
                        consensus=consensus,
                        trimmed_consensus=trimmed_consensus,
                        found_primers=found_primers,
                        rid=rid,
                        rid_min=rid_min,
                        actual_size=actual_size,
                        consensus_fasta_handle=consensus_fasta_handle,
                        sampled_ids=sampled_ids,
                        msa=msa,
                        sorted_cluster_ids=sorted_cluster_ids,
                        sorted_sampled_ids=sorted_sampled_ids,
                        iupac_count=iupac_count
                    )

        return clusters_with_ambiguities, total_ambiguity_positions

    def _write_discarded_reads(self) -> None:
        """Write discarded reads to a FASTQ file for inspection.

        Discards include:
        - Outlier reads removed during variant phasing
        - Reads from clusters filtered out by early filtering

        Output: cluster_debug/{sample_name}-discards.fastq
        """
        if not self.discarded_read_ids:
            return

        discards_file = os.path.join(self.debug_dir, f"{self.sample_name}-discards.fastq")
        with open(discards_file, 'w') as f:
            for seq_id in sorted(self.discarded_read_ids):
                if seq_id in self.records:
                    SeqIO.write(self.records[seq_id], f, "fastq")

        logging.info(f"Wrote {len(self.discarded_read_ids)} discarded reads to {discards_file}")

    def cluster(self, algorithm: str = "graph") -> None:
        """Perform complete clustering process with variant phasing and write output files.

        Pipeline:
            1. Initial clustering (MCL or greedy)
            2. Pre-phasing merge (combine HP-equivalent initial clusters)
            2b. Early filtering (optional, skip small clusters before expensive phasing)
            3. Variant detection + phasing (split clusters by haplotype)
            4. Post-phasing merge (combine HP-equivalent subclusters)
            5. Filtering (size and ratio thresholds)
            6. Output generation
            7. Write discarded reads (optional)

        Args:
            algorithm: Clustering algorithm to use ('graph' for MCL or 'greedy')
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            # Phase 1: Initial clustering
            initial_clusters = self._run_initial_clustering(temp_dir, algorithm)

            # Phase 2: Pre-phasing merge
            merged_clusters = self._run_prephasing_merge(initial_clusters)

            # Phase 2b: Early filtering (optional)
            clusters_to_phase, early_filtered = self._apply_early_filter(merged_clusters)

            # Phase 3: Variant detection + phasing
            all_subclusters = self._run_variant_phasing(clusters_to_phase)

            # Phase 4: Post-phasing merge
            merged_subclusters = self._run_postphasing_merge(all_subclusters)

            # Phase 5: Size filtering
            filtered_clusters = self._run_size_filtering(merged_subclusters)

            # Phase 6: Output generation
            consensus_output_file = os.path.join(self.output_dir, f"{self.sample_name}-all.fasta")
            clusters_with_ambiguities, total_ambiguity_positions = self._write_cluster_outputs(
                filtered_clusters, consensus_output_file
            )

            # Phase 7: Write discarded reads (optional)
            if self.collect_discards and self.discarded_read_ids:
                self._write_discarded_reads()

            # Write phasing statistics
            self.write_phasing_stats(
                initial_clusters_count=len(initial_clusters),
                after_prephasing_merge_count=len(merged_clusters),
                subclusters_count=len(all_subclusters),
                merged_count=len(merged_subclusters),
                final_count=len(filtered_clusters),
                clusters_with_ambiguities=clusters_with_ambiguities,
                total_ambiguity_positions=total_ambiguity_positions
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

    def run_spoa(self, read_sequences: Dict[str, str]) -> Optional[MSAResult]:
        """Run SPOA to generate consensus sequence and MSA.

        Args:
            read_sequences: Dictionary mapping read IDs to sequence strings

        Returns:
            MSAResult containing consensus sequence, raw MSA string, and parsed records,
            or None if SPOA fails or input is empty
        """
        if not read_sequences:
            return None

        try:
            # Create temporary input file with actual read IDs
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
                for read_id, seq in read_sequences.items():
                    f.write(f">{read_id}\n{seq}\n")
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

            # Parse MSA into ReadAlignment objects
            # Homopolymer normalization is enabled unless explicitly disabled
            enable_normalization = not self.disable_homopolymer_equivalence
            alignments, consensus, msa_to_consensus_pos = extract_alignments_from_msa(
                msa_output,
                enable_homopolymer_normalization=enable_normalization
            )

            if not consensus:
                raise RuntimeError("SPOA did not generate consensus sequence")

            return MSAResult(
                consensus=consensus,
                msa_string=msa_output,
                alignments=alignments,
                msa_to_consensus_pos=msa_to_consensus_pos
            )

        except subprocess.CalledProcessError as e:
            logging.error(f"SPOA failed with return code {e.returncode}")
            logging.error(f"Command: {' '.join(cmd)}")
            logging.error(f"Stderr: {e.stderr}")
            return None

        except RuntimeError:
            # Re-raise RuntimeError (e.g., consensus generation failure) as hard error
            raise

        except (FileNotFoundError, OSError) as e:
            # FileNotFoundError: SPOA not in PATH
            # OSError: Other OS-level errors (e.g., permission denied)
            logging.error(f"Error running SPOA: {str(e)}")
            return None

    def identify_outlier_reads(self, alignments: List[ReadAlignment], consensus_seq: str,
                              sampled_ids: Set[str], threshold: float) -> Tuple[Set[str], Set[str]]:
        """Identify outlier reads below identity threshold using homopolymer-normalized metrics.

        Uses normalized_edit_distance which excludes homopolymer length differences when
        homopolymer normalization is enabled (default behavior).

        Args:
            alignments: List of ReadAlignment objects from MSA
            consensus_seq: Ungapped consensus sequence
            sampled_ids: Set of read IDs that were sampled
            threshold: Identity threshold (e.g., 0.95 for 95% identity)
                      Reads with identity < threshold are considered outliers

        Returns:
            Tuple of (keep_ids, outlier_ids) where both are sets of read IDs
        """
        if not alignments or not consensus_seq:
            return sampled_ids, set()

        try:
            keep_ids = set()
            outlier_ids = set()

            consensus_length = len(consensus_seq)
            if consensus_length == 0:
                return sampled_ids, set()

            for alignment in alignments:
                # Calculate read identity using normalized edit distance
                # (excludes homopolymer length differences when normalization is enabled)
                error_rate = alignment.normalized_edit_distance / consensus_length
                identity = 1.0 - error_rate

                if identity >= threshold:
                    keep_ids.add(alignment.read_id)
                else:
                    outlier_ids.add(alignment.read_id)

            return keep_ids, outlier_ids

        except (ZeroDivisionError, AttributeError, TypeError) as e:
            # ZeroDivisionError: consensus_length is 0
            # AttributeError: alignment missing expected attributes
            # TypeError: unexpected None values
            logging.warning(f"Failed to identify outlier reads: {e}")
            return sampled_ids, set()

    def calculate_read_identity(self, alignments: List[ReadAlignment], consensus_seq: str) -> Tuple[Optional[float], Optional[float]]:
        """Calculate read identity metrics from MSA alignments using homopolymer-normalized metrics.

        Read identity measures how well individual reads agree with the consensus sequence.
        This is an internal consistency metric - it does NOT measure accuracy against ground truth.
        High values indicate homogeneous clustering and/or low read error rates.
        Low values may indicate heterogeneous clusters, outliers, or poor consensus quality (esp. at low RiC).

        Uses normalized_edit_distance which excludes homopolymer length differences when
        homopolymer normalization is enabled (default behavior).

        Args:
            alignments: List of ReadAlignment objects from MSA
            consensus_seq: Ungapped consensus sequence

        Returns:
            Tuple of (mean_rid, min_rid) where:
            - mean_rid: Mean read identity across all reads (0.0-1.0)
            - min_rid: Minimum read identity (worst-case read) (0.0-1.0)
            Returns (None, None) if no valid alignments or consensus
        """
        if not alignments or not consensus_seq:
            return None, None

        try:
            # Calculate read identity (1 - error_rate) for each read using normalized metrics
            consensus_length = len(consensus_seq)
            if consensus_length == 0:
                return None, None

            identities = []
            for alignment in alignments:
                # Use normalized edit distance (excludes homopolymer differences)
                error_rate = alignment.normalized_edit_distance / consensus_length
                identity = 1.0 - error_rate
                identities.append(identity)

            if not identities:
                return None, None

            # Calculate per-read statistics
            mean_rid = np.mean(identities)
            min_rid = np.min(identities)

            logging.debug(f"Calculated identity metrics: rid={mean_rid:.3f}, rid_min={min_rid:.3f}")

            return mean_rid, min_rid

        except (ZeroDivisionError, ValueError, TypeError) as e:
            # ZeroDivisionError: consensus_length is 0
            # ValueError: np.mean/np.min on invalid data
            # TypeError: unexpected None values in alignment data
            logging.warning(f"Failed to calculate read identity from MSA: {e}")
            return None, None

    def detect_variant_positions(self, alignments: List[ReadAlignment], consensus_seq: str,
                                msa_to_consensus_pos: Dict[int, Optional[int]]) -> List[Dict]:
        """Detect variant positions in MSA for future phasing.

        Identifies positions with significant alternative alleles that may indicate
        unphased biological variation requiring cluster subdivision.

        Args:
            alignments: List of ReadAlignment objects from MSA
            consensus_seq: Ungapped consensus sequence
            msa_to_consensus_pos: Mapping from MSA position to consensus position

        Returns:
            List of variant position dictionaries with metadata, or empty list if none found
        """
        if not alignments or not consensus_seq:
            return []

        try:
            # Extract consensus_aligned from first alignment (they all have the same length)
            if not alignments:
                logging.debug(f"No alignments found in MSA for variant detection")
                return []

            # Reconstruct consensus_aligned from consensus_seq and msa_to_consensus_pos
            msa_length = len(alignments[0].aligned_sequence)
            consensus_aligned = []
            for msa_pos in range(msa_length):
                cons_pos = msa_to_consensus_pos.get(msa_pos)
                if cons_pos is not None:
                    consensus_aligned.append(consensus_seq[cons_pos])
                else:
                    consensus_aligned.append('-')
            consensus_aligned = ''.join(consensus_aligned)

            # Analyze positional variation (error_threshold parameter unused but kept for compatibility)
            position_stats = analyze_positional_variation(alignments, consensus_aligned, msa_to_consensus_pos)

            # Identify variant positions
            variant_positions = []
            high_error_not_variants = []  # Track positions with high error but not detected as variants

            for pos_stat in position_stats:
                is_variant, variant_bases, reason = is_variant_position_with_composition(
                    pos_stat,
                    min_variant_frequency=self.min_variant_frequency,
                    min_variant_count=self.min_variant_count
                )

                if is_variant:
                    variant_positions.append({
                        'msa_position': pos_stat.msa_position,
                        'consensus_position': pos_stat.consensus_position,
                        'coverage': pos_stat.coverage,
                        'variant_bases': variant_bases,
                        'base_composition': pos_stat.base_composition,
                        'homopolymer_composition': pos_stat.homopolymer_composition,
                        'error_rate': pos_stat.error_rate,
                        'reason': reason
                    })
                elif pos_stat.error_rate >= self.min_variant_frequency:
                    # High error but didn't meet variant criteria (scattered errors)
                    high_error_not_variants.append({
                        'msa_position': pos_stat.msa_position,
                        'error_rate': pos_stat.error_rate,
                        'base_composition': pos_stat.base_composition,
                        'homopolymer_composition': pos_stat.homopolymer_composition,
                        'reason': reason
                    })

            # Debug logging for variant detection results
            if variant_positions:
                logging.debug(f"Variant detection: Found {len(variant_positions)} variant positions "
                            f"(thresholds: freq >= {self.min_variant_frequency*100:.1f}%, "
                            f"count >= {self.min_variant_count})")
                for var in variant_positions:
                    cons_pos_str = f"consensus:{var['consensus_position']}" if var['consensus_position'] is not None else "insertion"
                    logging.debug(f"  MSA pos {var['msa_position']} ({cons_pos_str}): "
                                f"error={var['error_rate']*100:.1f}%, variants={var['variant_bases']}, "
                                f"raw={var['base_composition']}, hp={var['homopolymer_composition']}")
            else:
                logging.debug(f"Variant detection: No variant positions detected "
                            f"(thresholds: freq >= {self.min_variant_frequency*100:.1f}%, "
                            f"count >= {self.min_variant_count})")

            if high_error_not_variants:
                logging.debug(f"Note: {len(high_error_not_variants)} positions have error >= {self.min_variant_frequency*100:.1f}% "
                            f"but no single alternative allele meets both thresholds (scattered errors)")
                for pos_info in high_error_not_variants[:3]:  # Show first 3 examples
                    logging.debug(f"  MSA pos {pos_info['msa_position']}: error={pos_info['error_rate']*100:.1f}%, "
                                f"raw={pos_info['base_composition']}, hp={pos_info['homopolymer_composition']}")

            return variant_positions

        except Exception as e:
            logging.warning(f"Failed to detect variant positions: {e}")
            return []

    def recursive_phase_cluster(
        self,
        read_ids: Set[str],
        read_sequences: Dict[str, str],
        path: List[str],
        depth: int = 0
    ) -> Tuple[List[Tuple[List[str], str, Set[str]]], Set[str]]:
        """Recursively phase a cluster, regenerating MSA at each level.

        At each recursion level:
        1. Generate new MSA for the current subset of reads
        2. Detect variant positions in the new MSA
        3. Find the best single-position split (minimum within-cluster error)
        4. Recurse on qualifying subclusters, collecting deferred reads

        This allows "scattered errors" to become detectable variants once
        correlated reads are grouped together through earlier phasing.

        Args:
            read_ids: Set of read IDs to phase at this level
            read_sequences: Dict mapping read_id -> sequence string
            path: List of alleles representing the path to this node
            depth: Current recursion depth (for logging)

        Returns:
            Tuple of (leaf_haplotypes, deferred_reads):
            - leaf_haplotypes: List of (path, consensus_seq, read_ids) for leaf nodes
            - deferred_reads: Set of read_ids that didn't qualify at any level
        """
        indent = "  " * depth
        total_reads = len(read_ids)
        path_str = '-'.join(path) if path else '(root)'

        # Base case: cluster too small to split
        if total_reads < self.min_variant_count * 2:
            logging.debug(f"{indent}Leaf (too small): path={path_str}, reads={total_reads}")
            # Generate consensus for this leaf
            leaf_seqs = {rid: read_sequences[rid] for rid in read_ids}
            result = self.run_spoa(leaf_seqs)
            consensus = result.consensus if result else ""
            return [(path, consensus, read_ids)], set()

        # Generate MSA for this cluster
        cluster_seqs = {rid: read_sequences[rid] for rid in read_ids}
        result = self.run_spoa(cluster_seqs)

        if result is None:
            logging.debug(f"{indent}Leaf (SPOA failed): path={path_str}, reads={total_reads}")
            return [(path, "", read_ids)], set()

        consensus = result.consensus
        alignments = result.alignments
        msa_to_consensus_pos = result.msa_to_consensus_pos

        # Detect variant positions in this new MSA
        variant_positions = self.detect_variant_positions(alignments, consensus, msa_to_consensus_pos)

        if not variant_positions:
            logging.debug(f"{indent}Leaf (no variants): path={path_str}, reads={total_reads}")
            return [(path, consensus, read_ids)], set()

        logging.debug(f"{indent}Level {depth}: path={path_str}, reads={total_reads}, variants={len(variant_positions)}")

        # Extract alleles at variant positions for each read
        # Parse MSA to get consensus_aligned for normalized allele extraction
        from io import StringIO
        msa_handle = StringIO(result.msa_string)
        records = list(SeqIO.parse(msa_handle, 'fasta'))

        consensus_aligned = None
        for record in records:
            if 'Consensus' in record.description or 'Consensus' in record.id:
                consensus_aligned = str(record.seq).upper()
                break

        if not consensus_aligned:
            logging.debug(f"{indent}Leaf (no consensus in MSA): path={path_str}")
            return [(path, consensus, read_ids)], set()

        # Build mapping from read_id to alignment
        read_to_alignment = {a.read_id: a for a in alignments}

        # Extract normalized alleles at variant positions
        variant_msa_positions = sorted([v['msa_position'] for v in variant_positions])
        read_to_position_alleles = {}

        for read_id in read_ids:
            alignment = read_to_alignment.get(read_id)
            if not alignment:
                continue

            aligned_seq = alignment.aligned_sequence
            score_aligned = alignment.score_aligned

            position_alleles = {}
            for msa_pos in variant_msa_positions:
                if msa_pos < len(aligned_seq):
                    allele = aligned_seq[msa_pos]
                    # Use normalized allele if homopolymer-equivalent
                    if score_aligned and msa_pos < len(score_aligned):
                        if score_aligned[msa_pos] == '=':
                            allele = consensus_aligned[msa_pos]
                    position_alleles[msa_pos] = allele
                else:
                    position_alleles[msa_pos] = '-'

            read_to_position_alleles[read_id] = position_alleles

        # Find the best single-position split
        best_pos = None
        best_error = float('inf')
        best_qualifying = None
        best_non_qualifying = None
        all_positions = set(variant_msa_positions)

        for pos in variant_msa_positions:
            # Group reads by allele at this position
            allele_groups = group_reads_by_single_position(
                read_to_position_alleles, pos, set(read_to_position_alleles.keys())
            )

            # Filter to qualifying haplotypes
            qualifying, non_qualifying = filter_qualifying_haplotypes(
                allele_groups, total_reads, self.min_variant_count, self.min_variant_frequency
            )

            # Need at least 2 qualifying groups to consider this split
            if len(qualifying) < 2:
                continue

            # Calculate within-cluster error for this split
            error = calculate_within_cluster_error(
                qualifying,
                read_to_position_alleles,
                {pos},
                all_positions
            )

            if error < best_error:
                best_error = error
                best_pos = pos
                best_qualifying = qualifying
                best_non_qualifying = non_qualifying

        # If no valid split found, this is a leaf node
        if best_pos is None:
            logging.debug(f"{indent}Leaf (no valid split): path={path_str}, reads={total_reads}")
            return [(path, consensus, read_ids)], set()

        # Log the split decision
        qual_summary = ', '.join(f"{a}:{len(r)}" for a, r in sorted(best_qualifying.items()))
        logging.debug(f"{indent}Split at pos {best_pos}: error={best_error:.4f}, qualifying=[{qual_summary}]")

        # Collect deferred reads from non-qualifying groups at this level
        all_deferred = set()
        for allele, reads in best_non_qualifying.items():
            all_deferred.update(reads)
            if reads:
                logging.debug(f"{indent}  Deferring {len(reads)} reads from allele '{allele}'")

        # Recurse on each qualifying sub-cluster
        all_leaves = []

        for allele, sub_read_ids in sorted(best_qualifying.items()):
            new_path = path + [allele]
            sub_leaves, sub_deferred = self.recursive_phase_cluster(
                sub_read_ids,
                read_sequences,
                new_path,
                depth + 1
            )
            all_leaves.extend(sub_leaves)
            all_deferred.update(sub_deferred)

        return all_leaves, all_deferred

    def phase_reads_by_variants(
        self,
        msa_string: str,
        consensus_seq: str,
        cluster_read_ids: Set[str],
        variant_positions: List[Dict],
        alignments: Optional[List[ReadAlignment]] = None
    ) -> List[Tuple[str, Set[str]]]:
        """Phase reads into haplotypes with recursive MSA regeneration.

        Splits a cluster into sub-clusters by recursively:
        1. Generating new MSA for each subcluster
        2. Rediscovering variant positions (previously "scattered errors" may become variants)
        3. Splitting on the best single position
        4. Reassigning deferred reads to nearest leaf consensus

        Args:
            msa_string: Initial MSA in FASTA format from SPOA (used for initial variant detection)
            consensus_seq: Ungapped consensus sequence (unused, kept for API compatibility)
            cluster_read_ids: Set of read IDs in this cluster
            variant_positions: List of variant position dicts from initial detect_variant_positions()
            alignments: Optional pre-parsed alignments (unused, kept for API compatibility)

        Returns:
            List of (allele_combo_string, read_id_set) tuples
            e.g., [("C-T-A", {id1, id2}), ("T-C-A", {id3, id4})]

            If cluster should not be split (0-1 qualifying haplotypes), returns:
            [(None, cluster_read_ids)]
        """
        if not variant_positions:
            # No variants to phase - return single group with all reads
            return [(None, cluster_read_ids)]

        try:
            # Build read_sequences dict for the cluster
            read_sequences = {rid: self.sequences[rid] for rid in cluster_read_ids if rid in self.sequences}

            if not read_sequences:
                logging.warning("No sequences found for cluster reads")
                return [(None, cluster_read_ids)]

            logging.info(f"Recursive phasing with MSA regeneration: {len(variant_positions)} initial variants, {len(read_sequences)} reads")

            # Run recursive phasing with MSA regeneration
            leaves, deferred = self.recursive_phase_cluster(
                set(read_sequences.keys()),
                read_sequences,
                path=[],
                depth=0
            )

            # Check if we got any splits
            if len(leaves) <= 1 and not deferred:
                if leaves:
                    path, consensus, reads = leaves[0]
                    if not path:  # Empty path means no splits happened
                        return [(None, cluster_read_ids)]
                return [(None, cluster_read_ids)]

            # Handle case where we have only 1 leaf but with deferred reads
            if len(leaves) == 1:
                return [(None, cluster_read_ids)]

            logging.info(f"Recursive phasing: {len(leaves)} leaf haplotypes, {len(deferred)} deferred reads")

            # Reassign deferred reads to nearest leaf haplotype by consensus alignment
            if deferred:
                logging.debug(f"Reassigning {len(deferred)} deferred reads to nearest leaf consensus")

                # For each deferred read, align to each leaf consensus and pick best match
                leaf_reads_updated = {tuple(path): set(reads) for path, consensus, reads in leaves}
                leaf_consensuses = {tuple(path): consensus for path, consensus, reads in leaves}

                for read_id in deferred:
                    if read_id not in read_sequences:
                        continue

                    read_seq = read_sequences[read_id]
                    min_distance = float('inf')
                    nearest_path = None

                    for path_tuple, consensus in leaf_consensuses.items():
                        if not consensus:
                            continue
                        result = edlib.align(read_seq, consensus)
                        distance = result['editDistance']
                        if distance < min_distance:
                            min_distance = distance
                            nearest_path = path_tuple

                    if nearest_path is not None:
                        leaf_reads_updated[nearest_path].add(read_id)

                # Update leaves with reassigned reads
                leaves = [(list(path), leaf_consensuses[path], reads)
                          for path, reads in leaf_reads_updated.items()]

            # Convert paths to allele combo strings and build result
            result = []
            sampled_read_count = len(read_sequences)

            for path, consensus, reads in leaves:
                combo = '-'.join(path) if path else 'unsplit'
                result.append((combo, reads))

            # Sort by size (largest first)
            result.sort(key=lambda x: len(x[1]), reverse=True)

            # Log results
            logging.info(f"Phasing decision: SPLITTING cluster into {len(result)} haplotypes")
            for i, (combo, reads) in enumerate(result, 1):
                logging.debug(f"  Haplotype {i}: '{combo}' with {len(reads)} reads ({len(reads)/sampled_read_count*100:.1f}%)")

            return result

        except Exception as e:
            logging.warning(f"Failed to phase reads by variants: {e}")
            import traceback
            logging.debug(traceback.format_exc())
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
                logging.debug(f"Loaded {total_primers} primers: {primer_count['forward']} forward, "
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

        # logging.debug(f"Starting primer trimming on sequence of length {len(sequence)}")

        # Look for primer at 5' end
        best_start_dist = float('inf')
        best_start_primer = None
        best_start_end = None

        for primer_name, primer_seq in self.primers:
            k = len(primer_seq) // 4  # Allow ~25% errors
            search_region = sequence[:len(primer_seq) * 2]

            result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)
            # logging.debug(f"Testing {primer_name} at 5' end (k={k})")

            if result["editDistance"] != -1:
                dist = result["editDistance"]
                if dist < best_start_dist:
                    best_start_dist = dist
                    best_start_primer = primer_name
                    best_start_end = result["locations"][0][1] + 1

        if best_start_primer:
            # logging.debug(f"Found {best_start_primer} at 5' end (k={best_start_dist})")
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
            # logging.debug(f"Found {best_end_primer} at 3' end (k={best_end_dist})")
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
        # Use max_repeat_motif_length=1 to be consistent with variant detection
        # (extract_alignments_from_msa also uses length=1)
        custom_params = AdjustmentParams(
            normalize_homopolymers=True,    # Enable homopolymer normalization
            handle_iupac_overlap=False,     # Disable IUPAC overlap handling
            normalize_indels=False,         # Disable indel normalization
            end_skip_distance=0,            # Disable end trimming
            max_repeat_motif_length=1       # Single-base repeats only (consistent with variant detection)
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
                # logging.debug(f"Non-homopolymer indel detected, sequences not merge-compatible")
                return -1  # Signal that merging should not occur
            if '.' in score_result.score_aligned:
                # logging.debug(f"Terminal overhang detected, sequences not merge-compatible")
                return -1  # Signal that merging should not occur
        
        # Count only substitutions (not homopolymer adjustments or indels)
        # Note: mismatches includes both substitutions and non-homopolymer indels
        # For accurate distance when indels are present, we use the mismatches count
        distance = score_result.mismatches
        
        # Log details about the variations found
        substitutions = score_result.score_aligned.count('X')
        indels = score_result.score_aligned.count('I')
        homopolymers = score_result.score_aligned.count('=')
        
        # logging.debug(f"Consensus distance: {distance} total mismatches "
        #              f"({substitutions} substitutions, {indels} indels, "
        #              f"{homopolymers} homopolymer adjustments)")
        
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
    parser.add_argument("--min-identity", type=float, default=0.9,
                        help="Minimum sequence identity threshold for clustering (default: 0.9)")
    parser.add_argument("--inflation", type=float, default=4.0,
                        help="MCL inflation parameter (default: 4.0)")
    parser.add_argument("--min-size", type=int, default=5,
                        help="Minimum cluster size (default: 5, 0 to disable)")
    parser.add_argument("--min-cluster-ratio", type=float, default=0.2,
                        help="Minimum size ratio between a cluster and the largest cluster (default: 0.2, 0 to disable)")
    parser.add_argument("--max-sample-size", type=int, default=500,
                        help="Maximum cluster size for consensus (default: 500)")
    parser.add_argument("--outlier-identity", type=float, default=None,
                        help="Minimum read-to-consensus identity to keep a read (default: auto). "
                             "Reads below this threshold are removed as outliers before final "
                             "consensus generation. Auto-calculated as (1 + min_identity) / 2. "
                             "This threshold is typically higher than --min-identity because "
                             "the consensus is error-corrected through averaging.")
    parser.add_argument("--disable-position-phasing", action="store_true",
                        help="Disable position-based variant phasing (enabled by default). "
                             "MCL graph clustering already separates most variants; this "
                             "second pass analyzes MSA positions to phase remaining variants.")
    parser.add_argument("--min-variant-frequency", type=float, default=0.20,
                        help="Minimum alternative allele frequency to call variant (default: 0.20 for 20%%)")
    parser.add_argument("--min-variant-count", type=int, default=5,
                        help="Minimum alternative allele read count to call variant (default: 5)")
    parser.add_argument("--min-ambiguity-frequency", type=float, default=0.10,
                        help="Minimum alternative allele frequency for IUPAC ambiguity calling (default: 0.10 for 10%%)")
    parser.add_argument("--min-ambiguity-count", type=int, default=3,
                        help="Minimum alternative allele read count for IUPAC ambiguity calling (default: 3)")
    parser.add_argument("--presample", type=int, default=1000,
                        help="Presample size for initial reads (default: 1000, 0 to disable)")
    parser.add_argument("--k-nearest-neighbors", type=int, default=5,
                        help="Number of nearest neighbors for graph construction (default: 5)")
    parser.add_argument("--enable-scalability", nargs='?', const=0, default=None, type=int,
                        metavar="THRESHOLD",
                        help="Enable scalable mode for large datasets (requires vsearch). "
                             "Optional: minimum sequence count to activate (default: 0 = always).")
    parser.add_argument("--threads", type=int, default=1, metavar="N",
                        help="Max threads for internal parallelism (vsearch, SPOA). "
                             "Default: 1. Use higher values for single large jobs.")
    parser.add_argument("--disable-early-filter", action="store_true",
                        help="Disable early filtering; process all clusters through variant phasing (default: early filter enabled)")
    parser.add_argument("--collect-discards", action="store_true",
                        help="Write discarded reads (outliers and filtered clusters) to cluster_debug/{sample}-discards.fastq")
    parser.add_argument("--primers", help="FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)")
    parser.add_argument("-O", "--output-dir", default="clusters",
                        help="Output directory for all files (default: clusters)")
    parser.add_argument("--disable-homopolymer-equivalence", action="store_true",
                        help="Disable homopolymer equivalence in cluster merging (only merge identical sequences)")
    parser.add_argument("--disable-cluster-merging", action="store_true",
                        help="Disable merging of clusters with identical consensus sequences")
    parser.add_argument("--disable-ambiguity-calling", action="store_true",
                        help="Disable IUPAC ambiguity code calling for unphased variant positions")
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
        disable_cluster_merging=args.disable_cluster_merging,
        output_dir=args.output_dir,
        outlier_identity_threshold=args.outlier_identity,
        enable_secondpass_phasing=not args.disable_position_phasing,
        min_variant_frequency=args.min_variant_frequency,
        min_variant_count=args.min_variant_count,
        min_ambiguity_frequency=args.min_ambiguity_frequency,
        min_ambiguity_count=args.min_ambiguity_count,
        enable_iupac_calling=not args.disable_ambiguity_calling,
        enable_scalability=args.enable_scalability,
        max_threads=args.threads,
        early_filter=not args.disable_early_filter,
        collect_discards=args.collect_discards
    )

    # Log configuration
    if args.outlier_identity is not None:
        logging.info(f"Outlier removal enabled: outlier_identity={args.outlier_identity*100:.1f}% (user-specified)")
    else:
        # Auto-calculated threshold
        auto_threshold = (1.0 + args.min_identity) / 2.0
        logging.info(f"Outlier removal enabled: outlier_identity={auto_threshold*100:.1f}% (auto-calculated from min_identity={args.min_identity*100:.1f}%)")

    if not args.disable_position_phasing:
        logging.info(f"Position-based variant phasing enabled: min_freq={args.min_variant_frequency:.0%}, "
                    f"min_count={args.min_variant_count}")

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

    if len(records) == 0:
        logging.warning("No sequences found in input file. Nothing to cluster.")
        sys.exit(0)

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
            logging.debug(f"Found primers.fasta in input directory: {auto_primer_path}")
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

