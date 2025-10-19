#!/usr/bin/env python3

from collections import defaultdict
import argparse
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

import edlib
from adjusted_identity import score_alignment, AdjustmentParams
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

__version__ = "0.3.3"

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
                 disable_stability: bool = False,
                 sample_name: str = "sample",
                 disable_homopolymer_equivalence: bool = False,
                 output_dir: str = "clusters"):
        self.min_identity = min_identity
        self.inflation = inflation
        self.min_size = min_size
        self.min_cluster_ratio = min_cluster_ratio
        self.max_sample_size = max_sample_size
        self.presample_size = presample_size
        self.k_nearest_neighbors = k_nearest_neighbors
        self.disable_stability = disable_stability
        self.sample_name = sample_name
        self.disable_homopolymer_equivalence = disable_homopolymer_equivalence
        self.output_dir = output_dir
        self.sequences = {}  # id -> sequence string
        self.records = {}  # id -> SeqRecord object
        self.id_map = {}  # short_id -> original_id
        self.rev_id_map = {}  # original_id -> short_id

        # Create output directory and debug subdirectory
        os.makedirs(self.output_dir, exist_ok=True)
        self.debug_dir = os.path.join(self.output_dir, "cluster_debug")
        os.makedirs(self.debug_dir, exist_ok=True)

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

    def assess_cluster_stability(self, cluster: Set[str],
                                 num_trials: int = 100,
                                 sample_size: int = 20) -> Tuple[str, Optional[float], Optional[float]]:
        """Assess cluster stability through subsampling."""
        # If stability assessment is disabled, return consensus with no stability metrics
        if self.disable_stability:
            full_seqs = [self.sequences[seq_id] for seq_id in cluster]
            full_consensus = self.run_spoa(full_seqs)
            return full_consensus, None, None

        # If cluster is too small for proper sampling, return just the consensus
        if len(cluster) < sample_size + 1:
            return self.run_spoa([self.sequences[seq_id] for seq_id in cluster]), None, None

        # Generate full cluster consensus
        full_seqs = [self.sequences[seq_id] for seq_id in cluster]
        full_consensus = self.run_spoa(full_seqs)
        if not full_consensus:
            return "", None, None

        # Determine the number of trials to run
        # If the number of possible combinations is small, run each combination once
        # otherwise run the specified number of random trials
        total_combinations = math.comb(len(cluster), sample_size)
        if total_combinations <= num_trials:
            # Enumerate all combinations
            if total_combinations <= 10000:  # Practical limit for combinations
                all_ids = list(cluster)
                combinations = []

                def generate_combinations(start, current_combination):
                    if len(current_combination) == sample_size:
                        combinations.append(current_combination.copy())
                        return
                    for i in range(start, len(all_ids)):
                        current_combination.append(all_ids[i])
                        generate_combinations(i + 1, current_combination)
                        current_combination.pop()

                generate_combinations(0, [])
                actual_trials = total_combinations
                trial_mode = "exhaustive"
            else:
                # Too many combinations, revert to random sampling
                combinations = None
                actual_trials = num_trials
                trial_mode = "random"
        else:
            combinations = None
            actual_trials = num_trials
            trial_mode = "random"

        # Run stability trials and collect differences
        differences = []
        with tqdm(total=actual_trials, desc="Assessing cluster stability") as pbar:
            if combinations:  # Use pre-generated combinations
                for combo in combinations:
                    sampled_seqs = [self.sequences[seq_id] for seq_id in combo]
                    trial_consensus = self.run_spoa(sampled_seqs)

                    if trial_consensus:
                        diff = self.calculate_consensus_distance(full_consensus, trial_consensus)
                        differences.append(diff)
                    pbar.update(1)
            else:  # Use random sampling
                for _ in range(actual_trials):
                    sampled_ids = random.sample(list(cluster), sample_size)
                    sampled_seqs = [self.sequences[seq_id] for seq_id in sampled_ids]
                    trial_consensus = self.run_spoa(sampled_seqs)

                    if trial_consensus:
                        diff = self.calculate_consensus_distance(full_consensus, trial_consensus)
                        differences.append(diff)
                    pbar.update(1)

        if differences:
            median_diff = statistics.median(differences)
            percentile_95 = np.percentile(differences, 95)
            logging.info(f"Stability assessment completed: median_diff={median_diff:.1f}, "
                         f"p95_diff={percentile_95:.1f} ({len(differences)} {trial_mode} trials)")
        else:
            median_diff = None
            percentile_95 = None
            logging.warning("No valid trials completed for stability assessment")

        return full_consensus, median_diff, percentile_95

    def merge_similar_clusters(self, clusters: List[Set[str]]) -> List[Set[str]]:
        """
        Merge clusters whose consensus sequences are identical or homopolymer-equivalent.
        Only counts unique consensus sequences when tracking merge information.

        Note: Primer trimming is performed before comparison to ensure clusters that differ
        only in primer regions are properly merged. Trimmed consensuses are used only for
        comparison and are discarded after merging.

        Args:
            clusters: List of clusters, where each cluster is a set of sequence IDs

        Returns:
            List of merged clusters
        """
        if not clusters:
            return []

        # Sort clusters by size, largest first
        clusters = sorted(clusters, key=len, reverse=True)

        # Generate a consensus sequence for each cluster
        logging.info("Generating consensus sequences for cluster merging...")
        consensuses = []
        cluster_to_consensus = {}  # Map from cluster index to its consensus

        with tqdm(total=len(clusters), desc="Generating cluster consensuses") as pbar:
            for i, cluster in enumerate(clusters):
                # Sample from larger clusters to speed up consensus generation
                if len(cluster) > self.max_sample_size:
                    # Sample by quality
                    qualities = []
                    for seq_id in cluster:
                        record = self.records[seq_id]
                        mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                        qualities.append((mean_quality, seq_id))

                    sampled_ids = [seq_id for _, seq_id in
                                   sorted(qualities, reverse=True)[:self.max_sample_size]]
                    sampled_seqs = [self.sequences[seq_id] for seq_id in sampled_ids]
                else:
                    sampled_seqs = [self.sequences[seq_id] for seq_id in cluster]

                consensus = self.run_spoa(sampled_seqs)

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
                new_cluster = set()
                for idx in equivalent_clusters:
                    new_cluster.update(clusters[idx])
                    merged_indices.add(idx)


                merged.append(new_cluster)

        # Add remaining unmerged clusters
        for i, cluster in enumerate(clusters):
            if i not in merged_indices:
                merged.append(cluster)

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
                            median_diff: Optional[float] = None,
                            p95_diff: Optional[float] = None,
                            actual_size: Optional[int] = None,
                            consensus_fasta_handle = None,
                            sampled_ids: Optional[Set[str]] = None) -> None:
        """Write cluster files: reads FASTQ and consensus FASTA/FASTG."""
        cluster_size = len(cluster)
        ric_size = min(actual_size or cluster_size, self.max_sample_size)

        # Create info string with size first
        info_parts = [f"size={cluster_size}", f"ric={ric_size}"]

        if median_diff is not None:
            info_parts.append(f"median_diff={median_diff:.1f}")
        if p95_diff is not None:
            info_parts.append(f"p95_diff={p95_diff:.1f}")
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
        """Perform complete clustering process and write output files.

        Args:
            algorithm: Clustering algorithm to use ('graph' for MCL or 'greedy')
        """
        # Create temporary directory outside the clustering method to allow
        # flexibility in choosing different clustering algorithms
        with tempfile.TemporaryDirectory() as temp_dir:
            # Select clustering algorithm
            if algorithm == "graph":
                try:
                    clusters = self.run_mcl_clustering(temp_dir)
                except (subprocess.SubprocessError, FileNotFoundError) as e:
                    logging.error(f"MCL clustering failed: {str(e)}")
                    logging.error("You may need to install MCL: https://micans.org/mcl/")
                    logging.error("Falling back to greedy clustering algorithm...")
                    clusters = self.run_greedy_clustering(temp_dir)
            elif algorithm == "greedy":
                clusters = self.run_greedy_clustering(temp_dir)
            else:
                raise ValueError(f"Unknown clustering algorithm: {algorithm}")

            # Sort by size (largest first)
            clusters.sort(key=lambda c: len(c), reverse=True)

            # Merge clusters with identical/homopolymer-equivalent consensus sequences
            # This happens BEFORE size filtering so that small clusters with identical
            # consensus can merge to meet the size threshold
            merged_clusters = self.merge_similar_clusters(clusters)

            # Re-sort by size after merging
            merged_clusters.sort(key=lambda c: len(c), reverse=True)

            # Filter clusters by size AFTER merging
            large_clusters = [c for c in merged_clusters if len(c) >= self.min_size]

            cluster_sizes = [len(c) for c in large_clusters]
            cluster_sizes_str = ', '.join(str(s) for s in cluster_sizes[:10])
            if len(cluster_sizes) > 10:
                cluster_sizes_str += f", ... ({len(cluster_sizes) - 10} more)"

            logging.info(
                f"After merging, found {len(large_clusters)} clusters meeting size threshold (>={self.min_size}): {cluster_sizes_str}")

            total_sequences = len(self.sequences)

            sequences_covered = 0
            clusters_to_output = []

            largest_cluster_size = len(large_clusters[0]) if large_clusters else 0

            skipped_clusters = []
            for cluster in large_clusters:
                cluster_size = len(cluster)

                # Apply min_cluster_ratio filter if specified and not the first cluster
                if self.min_cluster_ratio > 0 and largest_cluster_size > 0:
                    size_ratio = cluster_size / largest_cluster_size

                    # Skip this cluster if it's too small relative to the largest cluster
                    if size_ratio < self.min_cluster_ratio:
                        skipped_clusters.append((cluster_size, size_ratio))
                        continue

                clusters_to_output.append(cluster)
                sequences_covered += cluster_size

            # Log summary of skipped clusters if any were skipped
            if skipped_clusters:
                skipped_count = len(skipped_clusters)
                skipped_sizes = [size for size, _ in skipped_clusters]
                min_ratio = min([ratio for _, ratio in skipped_clusters])
                max_ratio = max([ratio for _, ratio in skipped_clusters])

                logging.info(
                    f"Skipped {skipped_count} clusters with sizes {min(skipped_sizes)}-{max(skipped_sizes)} "
                    f"(ratios to largest: {min_ratio:.3f}-{max_ratio:.3f} < {self.min_cluster_ratio})")

            logging.info(f"Outputting {len(clusters_to_output)} clusters covering {sequences_covered} sequences "
                         f"({sequences_covered / total_sequences:.1%} of total)")

            # Create the main consensus output file
            consensus_output_file = os.path.join(self.output_dir, f"{self.sample_name}-all.fasta")
            with open(consensus_output_file, 'w') as consensus_fasta_handle:
                for i, cluster in enumerate(clusters_to_output, 1):
                    actual_size = len(cluster)

                    # Check if this cluster was formed by merging
                    cluster_idx = large_clusters.index(cluster)

                    # Sample sequences for consensus generation if needed
                    if len(cluster) > self.max_sample_size:
                        logging.info(f"Sampling {self.max_sample_size} sequences from cluster {i} (size: {len(cluster)})")
                        qualities = []
                        for seq_id in cluster:
                            record = self.records[seq_id]
                            mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                            qualities.append((mean_quality, seq_id))
                        sampled_ids = {seq_id for _, seq_id in
                                       sorted(qualities, reverse=True)[:self.max_sample_size]}
                    else:
                        sampled_ids = cluster

                    # Generate consensus
                    consensus, median_diff, p95_diff = self.assess_cluster_stability(
                        sampled_ids,
                        num_trials=getattr(self, 'num_trials', 100),
                        sample_size=getattr(self, 'stability_sample', 20)
                    )

                    # For merged homopolymer-equivalent clusters, the consensus is already correct
                    # (we use the consensus from the larger cluster during merging)

                    if consensus:
                        # Perform primer trimming
                        trimmed_consensus = None
                        found_primers = None
                        if hasattr(self, 'primers'):
                            trimmed_consensus, found_primers = self.trim_primers(consensus)

                        self.write_cluster_files(
                            cluster_num=i,
                            cluster=cluster,
                            consensus=consensus,
                            trimmed_consensus=trimmed_consensus,
                            found_primers=found_primers,
                            median_diff=median_diff,
                            p95_diff=p95_diff,
                            actual_size=actual_size,
                            consensus_fasta_handle=consensus_fasta_handle,
                            sampled_ids=sampled_ids
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

    def run_spoa(self, sequences: List[str]) -> Optional[str]:
        """Run SPOA to generate consensus sequence."""
        if not sequences:
            return None

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
                "-r", "0",  # Result mode 0: consensus only
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

            # Parse SPOA output for consensus sequence
            consensus = None
            for line in result.stdout.split('\n'):
                if not line.startswith('>') and line.strip():
                    consensus = line.strip()
                    break

            if not consensus:
                logging.warning("SPOA did not generate consensus sequence")
                return sequences[0]  # Fall back to first sequence

            return consensus

        except subprocess.CalledProcessError as e:
            logging.error(f"SPOA failed with return code {e.returncode}")
            logging.error(f"Command: {' '.join(cmd)}")
            logging.error(f"Stderr: {e.stderr}")
            return sequences[0]  # Fall back to first sequence

        except Exception as e:
            logging.error(f"Error running SPOA: {str(e)}")
            return sequences[0]  # Fall back to first sequence

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
            
        # Get alignment from edlib
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
        if require_merge_compatible and 'I' in score_result.score_aligned:
            logging.debug(f"Non-homopolymer indel detected, sequences not merge-compatible")
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
        """Check if two sequences are equivalent when considering only homopolymer differences."""
        if not seq1 or not seq2:
            return seq1 == seq2
            
        # Use calculate_consensus_distance with merge compatibility check
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
    parser.add_argument("--presample", type=int, default=1000,
                        help="Presample size for initial reads (default: 1000, 0 to disable)")
    parser.add_argument("--k-nearest-neighbors", type=int, default=5,
                        help="Number of nearest neighbors for graph construction (default: 5)")
    parser.add_argument("--primers", help="FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)")
    parser.add_argument("-O", "--output-dir", default="clusters",
                        help="Output directory for all files (default: clusters)")
    parser.add_argument("--stability-trials", type=int, default=100,
                        help="Number of sampling trials to assess stability (default: 100)")
    parser.add_argument("--stability-sample", type=int, default=20,
                        help="Size of stability samples (default: 20)")
    parser.add_argument("--disable-stability", action="store_true",
                        help="Disable stability assessment")
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
        disable_stability=args.disable_stability,
        sample_name=sample,
        disable_homopolymer_equivalence=args.disable_homopolymer_equivalence,
        output_dir=args.output_dir
    )

    clusterer.num_trials = args.stability_trials
    clusterer.stability_sample = args.stability_sample

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
        clusterer.load_primers(args.primers)
    else:
        # Look for primers.fasta in the same directory as the input file
        input_dir = os.path.dirname(os.path.abspath(args.input_file))
        auto_primer_path = os.path.join(input_dir, "primers.fasta")
        
        if os.path.exists(auto_primer_path):
            logging.info(f"Found primers.fasta in input directory: {auto_primer_path}")
            clusterer.load_primers(auto_primer_path)
        else:
            logging.warning("No primer file specified and primers.fasta not found in input directory. Primer trimming will be disabled.")
    
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

    clusterer.cluster(algorithm=args.algorithm)
    print()

if __name__ == "__main__":
    main()

