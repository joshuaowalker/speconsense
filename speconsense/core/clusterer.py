"""Main SpecimenClusterer class for clustering and consensus generation."""

from collections import defaultdict
import json
import logging
import os
import statistics
import subprocess
import tempfile
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple

import edlib
from adjusted_identity import score_alignment, AdjustmentParams, ScoringFormat
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

try:
    from speconsense import __version__
except ImportError:
    __version__ = "dev"

from speconsense.msa import ReadAlignment
from speconsense.scalability import (
    VsearchCandidateFinder,
    ScalablePairwiseOperation,
    ScalabilityConfig,
)

from .workers import (
    ClusterProcessingConfig,
    ConsensusGenerationConfig,
    _run_spoa_worker,
    _process_cluster_worker,
    _generate_cluster_consensus_worker,
    _trim_primers_standalone,
    _phase_reads_by_variants_standalone,
)


class SpecimenClusterer:
    def __init__(self, min_identity: float = 0.9,
                 inflation: float = 4.0,
                 min_size: int = 5,
                 min_cluster_ratio: float = 0.2,
                 max_sample_size: int = 100,
                 presample_size: int = 1000,
                 k_nearest_neighbors: int = 20,
                 sample_name: str = "sample",
                 disable_homopolymer_equivalence: bool = False,
                 disable_cluster_merging: bool = False,
                 output_dir: str = "clusters",
                 outlier_identity_threshold: Optional[float] = None,
                 enable_secondpass_phasing: bool = True,
                 min_variant_frequency: float = 0.10,
                 min_variant_count: int = 5,
                 min_ambiguity_frequency: float = 0.10,
                 min_ambiguity_count: int = 3,
                 enable_iupac_calling: bool = True,
                 scale_threshold: int = 1001,
                 max_threads: int = 1,
                 early_filter: bool = False,
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
        self.scale_threshold = scale_threshold
        self.max_threads = max_threads
        self.early_filter = early_filter
        self.collect_discards = collect_discards
        self.discarded_read_ids: Set[str] = set()  # Track all discarded reads (outliers + filtered)

        # Initialize scalability configuration
        # scale_threshold: 0=disabled, N>0=enabled for datasets >= N sequences
        self.scalability_config = ScalabilityConfig(
            enabled=scale_threshold > 0,
            activation_threshold=scale_threshold,
            max_threads=max_threads
        )
        self._candidate_finder = None
        if scale_threshold > 0:
            finder = VsearchCandidateFinder(num_threads=max_threads)
            if finder.is_available:
                self._candidate_finder = finder

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
                "scale_threshold": self.scale_threshold,
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

        # Log scalability mode status for large datasets
        if len(self.sequences) >= self.scale_threshold and self.scale_threshold > 0:
            if self._candidate_finder is not None:
                logging.info(f"Scalability mode active for {len(self.sequences)} sequences (threshold: {self.scale_threshold})")
            else:
                logging.warning(f"Dataset has {len(self.sequences)} sequences (>= threshold {self.scale_threshold}) "
                               "but vsearch not found. Using brute-force.")

    def _get_scalable_operation(self) -> ScalablePairwiseOperation:
        """Get a ScalablePairwiseOperation for pairwise comparisons."""
        # Wrap calculate_similarity to match expected signature (seq1, seq2, id1, id2)
        # IDs are unused in core.py - only needed for primer-aware scoring in summarize.py
        return ScalablePairwiseOperation(
            candidate_finder=self._candidate_finder,
            scoring_function=lambda seq1, seq2, id1, id2: self.calculate_similarity(seq1, seq2),
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
            for id1, neighbors in sorted(knn_edges.items()):
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
                logging.debug(f"Cluster {i} is empty, skipping")
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
                # Parallel SPOA execution using ProcessPoolExecutor
                from concurrent.futures import ProcessPoolExecutor

                # Prepare work packages with config
                work_packages = [
                    (cluster_idx, sampled_seqs, self.disable_homopolymer_equivalence)
                    for cluster_idx, sampled_seqs in clusters_needing_spoa
                ]

                with ProcessPoolExecutor(max_workers=self.max_threads) as executor:
                    results = list(tqdm(
                        executor.map(_run_spoa_worker, work_packages),
                        total=len(work_packages),
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
                # Sequential SPOA execution using same worker function as parallel path
                for cluster_idx, sampled_seqs in clusters_needing_spoa:
                    _, result = _run_spoa_worker((cluster_idx, sampled_seqs, self.disable_homopolymer_equivalence))
                    if result is None:
                        logging.warning(f"Cluster {cluster_idx} produced no consensus, skipping")
                        continue
                    consensus = result.consensus
                    if hasattr(self, 'primers'):
                        consensus, _ = self.trim_primers(consensus)
                    consensuses.append(consensus)
                    cluster_to_consensus[cluster_idx] = consensus

        consensus_to_clusters = defaultdict(list)

        # IMPORTANT: Use cluster_to_consensus.items() instead of enumerate(consensuses)
        # because the feature branch processes single-read and multi-read clusters separately,
        # which changes the order in the consensuses list. The cluster_to_consensus dict
        # maintains the correct mapping from cluster index to consensus.

        if self.disable_homopolymer_equivalence:
            # Only merge exactly identical sequences
            # Sort by cluster index to match main branch iteration order
            for cluster_idx, consensus in sorted(cluster_to_consensus.items()):
                consensus_to_clusters[consensus].append(cluster_idx)
        else:
            # Group by homopolymer-equivalent sequences
            # Use scalable method when enabled and there are many clusters
            use_scalable = (
                self.scale_threshold > 0 and
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
                # Sort by cluster index to match main branch iteration order
                # (main branch iterates via enumerate(consensuses) where consensuses list
                # order matches clusters order; our dict may have different insertion order
                # due to single-read vs multi-read separation)
                for cluster_idx, consensus in sorted(cluster_to_consensus.items()):
                    # Find if this consensus is homopolymer-equivalent to any existing group
                    found_group = False
                    for existing_consensus in consensus_to_clusters.keys():
                        if self.are_homopolymer_equivalent(consensus, existing_consensus):
                            consensus_to_clusters[existing_consensus].append(cluster_idx)
                            found_group = True
                            break

                    if not found_group:
                        consensus_to_clusters[consensus].append(cluster_idx)

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
        split_count = 0
        logging.debug("Processing clusters for variant detection and phasing...")

        indexed_clusters = list(enumerate(merged_clusters, 1))

        # Create config object for workers (used by both parallel and sequential paths)
        config = ClusterProcessingConfig(
            outlier_identity_threshold=self.outlier_identity_threshold,
            enable_secondpass_phasing=self.enable_secondpass_phasing,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            min_variant_frequency=self.min_variant_frequency,
            min_variant_count=self.min_variant_count
        )

        # Build work packages with per-cluster data
        work_packages = []
        for initial_idx, cluster in indexed_clusters:
            cluster_seqs = {sid: self.sequences[sid] for sid in cluster}
            cluster_quals = {
                sid: statistics.mean(self.records[sid].letter_annotations["phred_quality"])
                for sid in cluster
            }
            work_packages.append((initial_idx, cluster, cluster_seqs, cluster_quals, config))

        if self.max_threads > 1 and len(merged_clusters) > 10:
            # Parallel processing with ProcessPoolExecutor for true parallelism
            from concurrent.futures import ProcessPoolExecutor

            with ProcessPoolExecutor(max_workers=self.max_threads) as executor:
                from tqdm import tqdm
                results = list(tqdm(
                    executor.map(_process_cluster_worker, work_packages),
                    total=len(work_packages),
                    desc="Processing clusters"
                ))

            # Collect results maintaining order
            for subclusters, discarded_ids in results:
                if len(subclusters) > 1:
                    split_count += 1
                all_subclusters.extend(subclusters)
                all_discarded.update(discarded_ids)
        else:
            # Sequential processing using same worker function as parallel path
            for work_package in work_packages:
                subclusters, discarded_ids = _process_cluster_worker(work_package)
                if len(subclusters) > 1:
                    split_count += 1
                all_subclusters.extend(subclusters)
                all_discarded.update(discarded_ids)

        # Update shared state after all processing complete
        self.discarded_read_ids.update(all_discarded)

        split_info = f" ({split_count} split)" if split_count > 0 else ""
        logging.info(f"After phasing, created {len(all_subclusters)} sub-clusters from {len(merged_clusters)} merged clusters{split_info}")
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
        small_clusters = [c for c in subclusters if len(c['read_ids']) < self.min_size]

        if small_clusters:
            filtered_count = len(small_clusters)
            logging.info(f"Filtered {filtered_count} clusters below minimum size ({self.min_size})")
            # Track discarded reads from size-filtered clusters
            for cluster in small_clusters:
                self.discarded_read_ids.update(cluster['read_ids'])

        # Filter by relative size ratio
        if large_clusters and self.min_cluster_ratio > 0:
            largest_size = max(len(c['read_ids']) for c in large_clusters)
            before_ratio_filter = len(large_clusters)
            passing_ratio = [c for c in large_clusters
                            if len(c['read_ids']) / largest_size >= self.min_cluster_ratio]
            failing_ratio = [c for c in large_clusters
                            if len(c['read_ids']) / largest_size < self.min_cluster_ratio]

            if failing_ratio:
                filtered_count = len(failing_ratio)
                logging.info(f"Filtered {filtered_count} clusters below minimum ratio ({self.min_cluster_ratio})")
                # Track discarded reads from ratio-filtered clusters
                for cluster in failing_ratio:
                    self.discarded_read_ids.update(cluster['read_ids'])

            large_clusters = passing_ratio

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

        # Create config for consensus generation workers
        primers = getattr(self, 'primers', None)
        config = ConsensusGenerationConfig(
            max_sample_size=self.max_sample_size,
            enable_iupac_calling=self.enable_iupac_calling,
            min_ambiguity_frequency=self.min_ambiguity_frequency,
            min_ambiguity_count=self.min_ambiguity_count,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            primers=primers
        )

        # Build work packages for each cluster
        work_packages = []
        for final_idx, cluster_dict in enumerate(clusters, 1):
            cluster = cluster_dict['read_ids']
            # Pre-compute quality means for each read
            qualities = {}
            for seq_id in cluster:
                record = self.records[seq_id]
                qualities[seq_id] = statistics.mean(record.letter_annotations["phred_quality"])
            # Extract sequences for this cluster
            sequences = {seq_id: self.sequences[seq_id] for seq_id in cluster}
            work_packages.append((final_idx, cluster, sequences, qualities, config))

        # Run consensus generation (parallel or sequential based on settings)
        if self.max_threads > 1 and len(clusters) > 4:
            # Parallel execution with ProcessPoolExecutor
            from concurrent.futures import ProcessPoolExecutor

            with ProcessPoolExecutor(max_workers=self.max_threads) as executor:
                results = list(tqdm(
                    executor.map(_generate_cluster_consensus_worker, work_packages),
                    total=len(work_packages),
                    desc="Final consensus generation"
                ))
        else:
            # Sequential execution using same worker function
            results = []
            for work_package in work_packages:
                result = _generate_cluster_consensus_worker(work_package)
                results.append(result)

        # Sort results by final_idx to ensure correct order
        results.sort(key=lambda r: r['final_idx'])

        # Write output files sequentially (I/O bound, must preserve order)
        with open(output_file, 'w') as consensus_fasta_handle:
            for result in results:
                final_idx = result['final_idx']
                cluster = result['cluster']
                actual_size = result['actual_size']

                # Log sampling info for large clusters
                if len(cluster) > self.max_sample_size:
                    logging.debug(f"Cluster {final_idx}: Sampling {self.max_sample_size} from {len(cluster)} reads for final consensus")

                consensus = result['consensus']
                iupac_count = result['iupac_count']

                if consensus:
                    if iupac_count > 0:
                        logging.debug(f"Cluster {final_idx}: Called {iupac_count} IUPAC ambiguity position(s)")
                        total_ambiguity_positions += iupac_count
                        clusters_with_ambiguities += 1

                    # Write output files
                    self.write_cluster_files(
                        cluster_num=final_idx,
                        cluster=cluster,
                        consensus=consensus,
                        trimmed_consensus=result['trimmed_consensus'],
                        found_primers=result['found_primers'],
                        rid=result['rid'],
                        rid_min=result['rid_min'],
                        actual_size=actual_size,
                        consensus_fasta_handle=consensus_fasta_handle,
                        sampled_ids=result['sampled_ids'],
                        msa=result['msa'],
                        sorted_cluster_ids=result['sorted_cluster_ids'],
                        sorted_sampled_ids=result['sorted_sampled_ids'],
                        iupac_count=iupac_count
                    )

        return clusters_with_ambiguities, total_ambiguity_positions

    def _write_discarded_reads(self) -> None:
        """Write discarded reads to a FASTQ file for inspection.

        Discards include:
        - Outlier reads removed during variant phasing
        - Reads from clusters filtered out by early filtering (Phase 2b)
        - Reads from clusters filtered out by size/ratio thresholds (Phase 5)
        - Reads filtered during orientation (when --orient-mode filter-failed)

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

    def phase_reads_by_variants(
        self,
        msa_string: str,
        consensus_seq: str,
        cluster_read_ids: Set[str],
        variant_positions: List[Dict],
        alignments: Optional[List[ReadAlignment]] = None
    ) -> List[Tuple[str, Set[str]]]:
        """Phase reads into haplotypes. Wrapper around standalone function.

        This method is provided for backward compatibility and testing.
        Internal processing uses _phase_reads_by_variants_standalone directly.
        """
        if not variant_positions:
            return [(None, cluster_read_ids)]

        # Build sequences dict from self.sequences
        read_sequences = {rid: self.sequences[rid] for rid in cluster_read_ids if rid in self.sequences}

        if not read_sequences:
            logging.warning("No sequences found for cluster reads")
            return [(None, cluster_read_ids)]

        config = ClusterProcessingConfig(
            outlier_identity_threshold=self.outlier_identity_threshold,
            enable_secondpass_phasing=self.enable_secondpass_phasing,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            min_variant_frequency=self.min_variant_frequency,
            min_variant_count=self.min_variant_count
        )

        return _phase_reads_by_variants_standalone(
            cluster_read_ids, self.sequences, variant_positions, config
        )

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
        """Trim primers from start and end of sequence. Wrapper around standalone function."""
        primers = getattr(self, 'primers', None)
        return _trim_primers_standalone(sequence, primers)

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
