"""Main SpecimenClusterer class for clustering and consensus generation."""

from collections import defaultdict
import copy
import json
import logging
import os
import statistics
import subprocess
import tempfile
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple

import edlib
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

try:
    from speconsense import __version__
except ImportError:
    __version__ = "dev"

from speconsense.msa import ReadAlignment, analyze_positional_variation, has_no_majority, call_iupac_ambiguities
from speconsense.context import classify_pairwise_differences
from speconsense.distances import count_variant_differences
from speconsense.qctx import DORADO_V5_0, get_qctx
from speconsense.significance import (
    compute_cer_factor,
    compute_critical_error_rate,
    compute_per_position_qstar,
    format_cer_details,
    is_cer_reportable,
)
from speconsense.scalability import (
    VsearchCandidateFinder,
    ScalablePairwiseOperation,
    ScalabilityConfig,
)

from .workers import (
    ClusterProcessingConfig,
    ConsensusGenerationConfig,
    _run_spoa_worker,
    _run_spoa_for_cluster_worker,
    _process_cluster_worker,
    _generate_cluster_consensus_worker,
    _trim_primers_standalone,
    _phase_reads_by_variants_standalone,
    _detect_variant_positions_standalone,
)


def _extract_consensus_aligned_from_msa(msa_string: str) -> Optional[str]:
    """Extract SPOA's consensus aligned sequence from MSA FASTA output."""
    from io import StringIO
    msa_handle = StringIO(msa_string)
    records = list(SeqIO.parse(msa_handle, 'fasta'))
    for record in records:
        if 'Consensus' in record.description or 'Consensus' in record.id:
            return str(record.seq).upper()
    return None


def _find_variant_columns(alignment_by_id: Dict, consensus_spoa_ids: Dict[str, int],
                          spoa_consensus_aligned: str) -> List[int]:
    """Find MSA columns where any two cluster consensuses differ.

    Args:
        alignment_by_id: {id: ReadAlignment} for all MSA entries
        consensus_spoa_ids: {spoa_id: cluster_index} for cluster consensus entries
        spoa_consensus_aligned: SPOA's own consensus aligned sequence (for HP normalization)

    Returns:
        Sorted list of MSA column indices with variant differences.
    """
    msa_length = len(spoa_consensus_aligned)

    # Extract HP-normalized base for each cluster consensus at each position
    cluster_bases: Dict[int, Dict[int, str]] = {}  # cluster_idx -> {msa_pos -> base}
    for spoa_id, cluster_idx in consensus_spoa_ids.items():
        aln = alignment_by_id.get(spoa_id)
        if not aln:
            continue
        bases = {}
        for pos in range(msa_length):
            base = aln.aligned_sequence[pos] if pos < len(aln.aligned_sequence) else '-'
            if aln.score_aligned and pos < len(aln.score_aligned):
                if aln.score_aligned[pos] == '=':
                    base = spoa_consensus_aligned[pos]
            bases[pos] = base
        cluster_bases[cluster_idx] = bases

    if len(cluster_bases) < 2:
        return []

    # Find positions where any two consensus bases differ
    cluster_indices = list(cluster_bases.keys())
    variant_columns = []
    for pos in range(msa_length):
        bases_at_pos = set()
        for idx in cluster_indices:
            b = cluster_bases[idx].get(pos, '-')
            if b != '-':
                bases_at_pos.add(b)
        if len(bases_at_pos) > 1:
            variant_columns.append(pos)

    return variant_columns


def _extract_bases_at_columns(aln, columns: List[int], spoa_consensus_aligned: str) -> Dict[int, str]:
    """Extract HP-normalized bases at specific MSA columns from an alignment."""
    bases = {}
    for col in columns:
        base = aln.aligned_sequence[col] if col < len(aln.aligned_sequence) else '-'
        if aln.score_aligned and col < len(aln.score_aligned):
            if aln.score_aligned[col] == '=':
                base = spoa_consensus_aligned[col] if col < len(spoa_consensus_aligned) else '-'
        bases[col] = base
    return bases


def _hp_normalized_pairwise_compare(seq1: str, seq2: str, min_hp_length: int = 6) -> Tuple[float, bool]:
    """Compare two sequences with HP-normalized global alignment.

    Runs edlib global alignment, builds HP context from both aligned sequences,
    and scores each position. Within qualifying HP runs (length >= min_hp_length),
    insertions/deletions of the HP base are treated as matches.

    Args:
        seq1, seq2: Ungapped sequences to compare
        min_hp_length: Minimum HP run length to normalize

    Returns:
        (distance, is_hp_equivalent) where:
        - distance: 1.0 - (matches / scored_positions), or 1.0 on failure
        - is_hp_equivalent: True if all differences are HP length differences
    """
    from speconsense.msa import _build_hp_context

    if not seq1 or not seq2:
        return 1.0, False
    if seq1 == seq2:
        return 0.0, True

    result = edlib.align(seq1, seq2, task="path")
    if result["editDistance"] == -1:
        return 1.0, False

    alignment = edlib.getNiceAlignment(result, seq1, seq2)
    if not alignment:
        return 1.0, False

    s1 = alignment['query_aligned']
    s2 = alignment['target_aligned']
    n = len(s1)

    hp1 = _build_hp_context(s1, min_hp_length)
    hp2 = _build_hp_context(s2, min_hp_length)

    matches = 0
    mismatches = 0
    hp_equiv = 0

    for i in range(n):
        b1, b2 = s1[i], s2[i]
        if b1 == '-' and b2 == '-':
            continue  # not a scored position
        if b1 == b2:
            matches += 1
        elif hp1[i] is not None and b1 == hp1[i] and b2 == '-':
            hp_equiv += 1  # deletion of HP base in seq2
        elif hp2[i] is not None and b2 == hp2[i] and b1 == '-':
            hp_equiv += 1  # deletion of HP base in seq1
        else:
            mismatches += 1

    scored = matches + hp_equiv + mismatches
    if scored == 0:
        return 0.0, True

    distance = mismatches / scored
    return distance, mismatches == 0


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
                 collect_discards: bool = False,
                 assumed_error_rate: float = 0.015,
                 significance_level: float = 1e-5,
                 group_identity: float = 0.85):
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
        self.assumed_error_rate = assumed_error_rate
        self.significance_level = significance_level
        self.group_identity = group_identity
        # Context-aware CER configuration. The factor threshold controls
        # whether a candidate cluster passes pairwise CER validation:
        # candidates pass when their per-position CER factor (q*/q_ctx) is
        # >= cer_factor_threshold. Setting assumed_error_rate <= 0 disables
        # the gate entirely; this is preserved as the disable sentinel.
        self.qctx_table = DORADO_V5_0
        self.cer_factor_threshold = 0.0 if assumed_error_rate <= 0 else 1.0
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
        """Write run metadata to JSON file for use by post-processing tools.

        Includes per-cluster CER reproduction data when clustering has
        completed (self.final_cluster_dicts populated). The schema follows the
        CER-in-practice paper Appendix A: parameters identify the q_ctx table
        and CER thresholds, identity_groups describe pairwise comparison
        membership, and variants contain N, M, K, L, context tags, and
        per-position q_ctx values sufficient to recompute cer_factor and
        cer_pstar under revised assumptions without re-running the pipeline.
        """
        metadata = {
            "schema_version": "2.0",
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
                "assumed_error_rate": self.assumed_error_rate,
                "significance_level": self.significance_level,
                "group_identity": self.group_identity,
                "cer_factor_threshold": self.cer_factor_threshold,
                "qctx_table": "dorado-v5.0",
            },
            "input_file": self.input_file,
            "augment_input": self.augment_input,
            "total_input_reads": self.total_input_reads,
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

        # Per-cluster CER reproduction data (populated after cluster() completes)
        cluster_dicts = getattr(self, 'final_cluster_dicts', None)
        if cluster_dicts is not None:
            metadata["identity_groups"] = self._build_identity_group_summary(cluster_dicts)
            metadata["variants"] = [
                self._build_variant_record(c) for c in cluster_dicts
            ]

        # Write metadata file
        metadata_file = os.path.join(self.debug_dir, f"{self.sample_name}-metadata.json")
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logging.debug(f"Wrote run metadata to {metadata_file}")

    @staticmethod
    def _build_identity_group_summary(cluster_dicts: List[Dict]) -> List[Dict]:
        """Aggregate per-cluster identity_group_id tags into group records."""
        groups: Dict[str, Dict] = {}
        for c in cluster_dicts:
            gid = c.get('identity_group_id')
            if gid is None:
                continue
            if gid not in groups:
                groups[gid] = {
                    'group_id': gid,
                    'group_N': c.get('cer_group_N'),
                    'members': [],
                }
            groups[gid]['members'].append(c.get('cluster_id') or str(c.get('cluster_idx')))
        return [groups[k] for k in sorted(groups)]

    @staticmethod
    def _build_variant_record(cluster_dict: Dict) -> Dict:
        """Serialize one cluster's CER reproduction data."""
        details = cluster_dict.get('cer_details') or {}
        return {
            'cluster_id': cluster_dict.get('cluster_id'),
            'identity_group': cluster_dict.get('identity_group_id'),
            'cer_status': cluster_dict.get('cer_status'),
            'M': len(cluster_dict.get('read_ids', [])),
            'N': cluster_dict.get('cer_group_N'),
            'K': details.get('K') if details else None,
            'context_tags': details.get('tags') if details else None,
            'q_ctx_per_position': details.get('q_ctx') if details else None,
            'compared_against_idx': details.get('ref_idx') if details else None,
            'cer_factor': cluster_dict.get('cer_factor'),
            'cer_pstar': cluster_dict.get('cer_pstar'),
        }

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

        self.total_input_reads = len(self.sequences)

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
                        logging.warning(f"Cluster {cluster_idx} produced no consensus, discarding {len(clusters[cluster_idx]['read_ids'])} reads")
                        self.discarded_read_ids.update(clusters[cluster_idx]['read_ids'])
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
                        logging.warning(f"Cluster {cluster_idx} produced no consensus, discarding {len(clusters[cluster_idx]['read_ids'])} reads")
                        self.discarded_read_ids.update(clusters[cluster_idx]['read_ids'])
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
                            iupac_count: int = 0,
                            cer_factor: Optional[float] = None,
                            cer_details: Optional[str] = None) -> None:
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
        if cer_factor is not None:
            if cer_factor == float('inf'):
                info_parts.append("cer_factor=inf")
            else:
                info_parts.append(f"cer_factor={cer_factor:.3f}")
        if cer_details is not None:
            info_parts.append(f"cer_details={cer_details}")
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
        clustered_count = sum(len(c) for c in initial_clusters)
        logging.info(f"Initial clustering produced {len(initial_clusters)} clusters "
                     f"covering {clustered_count}/{self.total_input_reads} reads")
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
            min_variant_count=self.min_variant_count,
            total_specimen_reads=self.total_input_reads,
            assumed_error_rate=self.assumed_error_rate,
            significance_level=self.significance_level,
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

    def _filter_noisy_clusters(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4a: Remove unreliable small clusters with no-majority positions.

        Only applies to clusters below the phasing floor (< min_variant_count * 2),
        since larger clusters have enough reads for reliable consensus. Disbanded
        reads go to discards for recovery by the global discard reassignment pass.
        """
        phasing_floor = self.min_variant_count * 2
        result = []
        total_disbanded = 0
        disbanded_reads = 0

        for cluster_dict in subclusters:
            read_ids = cluster_dict['read_ids']

            if len(read_ids) >= phasing_floor:
                result.append(cluster_dict)
                continue

            # Generate MSA for this small cluster
            qualities = {
                sid: statistics.mean(self.records[sid].letter_annotations["phred_quality"])
                for sid in read_ids
            }
            sorted_ids = sorted(read_ids, key=lambda x: (-qualities.get(x, 0), x))
            cluster_seqs = {sid: self.sequences[sid] for sid in sorted_ids}

            msa_result = _run_spoa_for_cluster_worker(
                cluster_seqs, self.disable_homopolymer_equivalence)

            if msa_result is None:
                self.discarded_read_ids.update(read_ids)
                total_disbanded += 1
                disbanded_reads += len(read_ids)
                continue

            # Build consensus_aligned for positional analysis
            msa_length = max(msa_result.msa_to_consensus_pos.keys()) + 1
            consensus_aligned = []
            for pos in range(msa_length):
                cons_pos = msa_result.msa_to_consensus_pos.get(pos)
                if cons_pos is not None and cons_pos < len(msa_result.consensus):
                    consensus_aligned.append(msa_result.consensus[cons_pos])
                else:
                    consensus_aligned.append('-')
            consensus_aligned = ''.join(consensus_aligned)

            position_stats = analyze_positional_variation(
                msa_result.alignments, consensus_aligned, msa_result.msa_to_consensus_pos)

            no_majority_count = sum(
                1 for ps in position_stats
                if ps.consensus_position is not None and has_no_majority(ps)
            )

            if no_majority_count > 0:
                self.discarded_read_ids.update(read_ids)
                total_disbanded += 1
                disbanded_reads += len(read_ids)
                logging.debug(f"Noise filter: disbanded cluster(n={len(read_ids)}) "
                             f"with {no_majority_count} no-majority positions")
            else:
                result.append(cluster_dict)

        if total_disbanded > 0:
            logging.info(f"Noise filter: disbanded {total_disbanded} unreliable clusters "
                        f"({disbanded_reads} reads moved to discards)")

        return result

    def _run_read_reassignment(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4b: Reassign reads to best-matching clusters within identity groups.

        Uses a shared MSA per identity group to identify variant positions between
        cluster consensuses, then moves reads to the cluster they best match at
        those positions. Iterates until stable.
        """
        if len(subclusters) <= 1:
            return subclusters

        # Generate consensus for each cluster
        consensuses = {}
        for i, cluster_dict in enumerate(subclusters):
            consensuses[i] = self._generate_consensus_for_validation(cluster_dict['read_ids'])

        # Form identity groups
        identity_groups = self._form_identity_groups(subclusters, consensuses)

        # Reassign within each multi-cluster group
        total_reassigned = 0
        for group_indices in identity_groups.values():
            if len(group_indices) <= 1:
                continue
            reassigned = self._reassign_reads_in_group(subclusters, consensuses, group_indices)
            total_reassigned += reassigned

        if total_reassigned > 0:
            logging.info(f"Read reassignment: moved {total_reassigned} reads")
            # Remove empty clusters
            subclusters = [c for c in subclusters if c['read_ids']]
        else:
            logging.info("Read reassignment: no reads moved")

        return subclusters

    def _run_second_phasing_pass(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4b2: Re-phase subclusters after read reassignment.

        Reassignment can move reads into clusters where they create detectable
        variants. This pass splits any such clusters. No outlier removal —
        just variant detection and phasing.
        """
        config = ClusterProcessingConfig(
            outlier_identity_threshold=self.outlier_identity_threshold,
            enable_secondpass_phasing=self.enable_secondpass_phasing,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            min_variant_frequency=self.min_variant_frequency,
            min_variant_count=self.min_variant_count,
            total_specimen_reads=self.total_input_reads,
            assumed_error_rate=self.assumed_error_rate,
            significance_level=self.significance_level,
        )

        result = []
        split_count = 0

        for cluster_dict in subclusters:
            read_ids = cluster_dict['read_ids']
            if len(read_ids) < config.min_variant_count * 2:
                result.append(cluster_dict)
                continue

            # Build quality-sorted sequences
            qualities = {
                sid: statistics.mean(self.records[sid].letter_annotations["phred_quality"])
                for sid in read_ids
            }
            sorted_ids = sorted(read_ids, key=lambda x: (-qualities.get(x, 0), x))
            cluster_seqs = {sid: self.sequences[sid] for sid in sorted_ids}

            msa_result = _run_spoa_for_cluster_worker(cluster_seqs, self.disable_homopolymer_equivalence)
            if msa_result is None:
                result.append(cluster_dict)
                continue

            variant_positions = _detect_variant_positions_standalone(
                msa_result.alignments, msa_result.consensus, msa_result.msa_to_consensus_pos,
                config.min_variant_frequency, config.min_variant_count
            )

            n_discards = sum(1 for sid in read_ids if sid.startswith('d-'))
            cluster_label = f"cluster(n={len(read_ids)}, d={n_discards})"

            if not variant_positions:
                logging.debug(f"Second phasing: {cluster_label} — no variants detected")
                result.append(cluster_dict)
                continue

            logging.debug(f"Second phasing: {cluster_label} — {len(variant_positions)} variant positions detected")

            # Phase
            phased = _phase_reads_by_variants_standalone(
                read_ids, self.sequences, qualities, variant_positions, config
            )

            if len(phased) <= 1:
                logging.debug(f"Second phasing: {cluster_label} — variants detected but phasing produced {len(phased)} group(s)")
                result.append(cluster_dict)
                continue

            # Split succeeded — track any reads lost during phasing
            split_count += 1
            phased_reads = set()
            for allele_combo, haplotype_reads in phased:
                phased_reads.update(haplotype_reads)
                result.append({
                    'read_ids': haplotype_reads,
                    'initial_cluster_num': cluster_dict.get('initial_cluster_num'),
                    'allele_combo': allele_combo,
                })
            lost = read_ids - phased_reads
            if lost:
                self.discarded_read_ids.update(lost)

        if split_count > 0:
            logging.info(f"Second phasing pass: {split_count} clusters split, "
                        f"{len(result)} sub-clusters from {len(subclusters)}")
            # Remove empty clusters
            result = [c for c in result if c['read_ids']]
        else:
            logging.info("Second phasing pass: no clusters split")

        return result

    def _reassign_reads_in_group(self, subclusters: List[Dict],
                                 consensuses: Dict[int, Optional[str]],
                                 group_indices: List[int],
                                 max_iterations: int = 3) -> int:
        """Reassign reads to best-matching clusters within an identity group.

        Returns total number of reads moved.
        """
        total_reassigned = 0

        for iteration in range(max_iterations):
            # Collect all reads and consensus sequences in this group
            cluster_consensus_seqs: Dict[int, str] = {}
            read_to_cluster: Dict[str, int] = {}
            all_read_seqs: Dict[str, str] = {}

            for idx in group_indices:
                consensus = consensuses.get(idx)
                if not consensus:
                    continue
                cluster_consensus_seqs[idx] = consensus
                for rid in subclusters[idx]['read_ids']:
                    read_to_cluster[rid] = idx
                    all_read_seqs[rid] = self.sequences[rid]

            if len(cluster_consensus_seqs) <= 1:
                break

            # Build SPOA input: consensus sequences (special IDs sort first) + reads
            spoa_input: Dict[str, str] = {}
            consensus_spoa_ids: Dict[str, int] = {}
            for idx, cons in cluster_consensus_seqs.items():
                spoa_id = f"__cons_{idx}__"
                spoa_input[spoa_id] = cons
                consensus_spoa_ids[spoa_id] = idx

            # Subsample reads if group is large
            max_reads = 500
            if len(all_read_seqs) > max_reads:
                per_cluster_cap = max(10, max_reads // len(cluster_consensus_seqs))
                sampled_reads: Dict[str, str] = {}
                for idx in group_indices:
                    cluster_rids = sorted(subclusters[idx]['read_ids'],
                        key=lambda x: (-statistics.mean(
                            self.records[x].letter_annotations["phred_quality"]), x))
                    for rid in cluster_rids[:per_cluster_cap]:
                        sampled_reads[rid] = self.sequences[rid]
                spoa_input.update(sampled_reads)
                reads_in_msa = set(sampled_reads.keys())
            else:
                spoa_input.update(all_read_seqs)
                reads_in_msa = set(all_read_seqs.keys())

            # Run SPOA
            msa_result = _run_spoa_for_cluster_worker(spoa_input, self.disable_homopolymer_equivalence)
            if not msa_result:
                break

            # Extract SPOA's consensus aligned sequence for HP normalization
            spoa_consensus_aligned = _extract_consensus_aligned_from_msa(msa_result.msa_string)
            if not spoa_consensus_aligned:
                break

            # Build alignment lookup
            alignment_by_id = {a.read_id: a for a in msa_result.alignments}

            # Find variant columns
            variant_columns = _find_variant_columns(alignment_by_id, consensus_spoa_ids,
                                                    spoa_consensus_aligned)
            if not variant_columns:
                break

            logging.debug(f"Reassignment iteration {iteration + 1}: "
                         f"{len(cluster_consensus_seqs)} clusters, "
                         f"{len(reads_in_msa)} reads, "
                         f"{len(variant_columns)} variant columns")

            # Extract cluster consensus bases at variant columns
            cluster_bases: Dict[int, Dict[int, str]] = {}
            for spoa_id, idx in consensus_spoa_ids.items():
                aln = alignment_by_id.get(spoa_id)
                if aln:
                    cluster_bases[idx] = _extract_bases_at_columns(
                        aln, variant_columns, spoa_consensus_aligned)

            # Score each read and reassign
            iteration_moves = 0
            for rid in reads_in_msa:
                if rid not in read_to_cluster:
                    continue
                current_idx = read_to_cluster[rid]
                aln = alignment_by_id.get(rid)
                if not aln:
                    continue

                read_bases = _extract_bases_at_columns(aln, variant_columns, spoa_consensus_aligned)

                # Concordance: count matches at variant positions for each cluster
                best_idx = current_idx
                best_score = -1
                for idx, cbases in cluster_bases.items():
                    score = sum(1 for col in variant_columns
                                if read_bases.get(col) == cbases.get(col))
                    if score > best_score or (score == best_score and idx == current_idx):
                        best_score = score
                        best_idx = idx

                if best_idx != current_idx:
                    # Verify via edit distance: read must be at least as close
                    # to target consensus as to source consensus
                    read_seq = all_read_seqs[rid]
                    source_cons = cluster_consensus_seqs.get(current_idx, '')
                    target_cons = cluster_consensus_seqs.get(best_idx, '')
                    if source_cons and target_cons:
                        dist_source = edlib.align(read_seq, source_cons)['editDistance']
                        dist_target = edlib.align(read_seq, target_cons)['editDistance']
                        if dist_target > dist_source:
                            continue  # Read is closer to current cluster overall

                    subclusters[current_idx]['read_ids'].discard(rid)
                    subclusters[best_idx]['read_ids'].add(rid)
                    read_to_cluster[rid] = best_idx
                    iteration_moves += 1

            total_reassigned += iteration_moves

            if iteration_moves == 0:
                break

            logging.debug(f"Reassignment iteration {iteration + 1}: {iteration_moves} reads moved")

            # Regenerate consensus for changed clusters
            for idx in group_indices:
                if subclusters[idx]['read_ids']:
                    consensuses[idx] = self._generate_consensus_for_validation(
                        subclusters[idx]['read_ids'])
                else:
                    consensuses[idx] = None

        return total_reassigned

    def _run_discard_reassignment(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4b3: Reassign discarded reads to existing clusters.

        For each discarded read:
        1. Find best-matching cluster consensus via edlib (group selection)
        2. For multi-cluster groups, use SPOA MSA + variant concordance
        3. Confirm via edit distance

        Returns updated subclusters with recovered reads added.
        """
        if not self.discarded_read_ids or not subclusters:
            return subclusters

        # Generate consensus for each cluster
        consensuses: Dict[int, Optional[str]] = {}
        for i, cluster_dict in enumerate(subclusters):
            consensuses[i] = self._generate_consensus_for_validation(cluster_dict['read_ids'])

        valid_consensuses = {i: c for i, c in consensuses.items() if c}
        if not valid_consensuses:
            return subclusters

        # Form identity groups (reuse existing logic)
        identity_groups = self._form_identity_groups(subclusters, consensuses)
        # Build reverse map: cluster_index -> group_id
        cluster_to_group: Dict[int, int] = {}
        for group_id, indices in identity_groups.items():
            for idx in indices:
                cluster_to_group[idx] = group_id

        # Step 1: Screen each discard against all cluster consensuses
        # Group candidates by identity group
        group_candidates: Dict[int, List[Tuple[str, str]]] = defaultdict(list)  # group_id -> [(read_id, seq)]
        rejected = 0

        # Use relaxed threshold: group_identity (0.85) since we just need group membership
        max_distance = 1.0 - self.group_identity

        for rid in sorted(self.discarded_read_ids):
            read_seq = self.sequences.get(rid)
            if not read_seq:
                continue

            best_idx = None
            best_dist = float('inf')
            for idx, cons in valid_consensuses.items():
                dist = edlib.align(read_seq, cons)['editDistance']
                norm_dist = dist / max(len(read_seq), len(cons))
                if norm_dist < best_dist:
                    best_dist = norm_dist
                    best_idx = idx

            if best_idx is None or best_dist > max_distance:
                rejected += 1
                continue

            group_id = cluster_to_group[best_idx]
            group_candidates[group_id].append((rid, read_seq))

        total_candidates = sum(len(v) for v in group_candidates.values())
        if total_candidates == 0:
            logging.info(f"Discard reassignment: no candidates (all {rejected} reads below threshold)")
            return subclusters

        def reassign_read(rid, target_idx):
            """Move a discarded read to a cluster with 'd-' prefix for traceability."""
            new_rid = f"d-{rid}"
            self.sequences[new_rid] = self.sequences[rid]
            if rid in self.records:
                rec = copy.copy(self.records[rid])
                rec.id = new_rid
                self.records[new_rid] = rec
            subclusters[target_idx]['read_ids'].add(new_rid)
            self.discarded_read_ids.discard(rid)

        # Step 2 & 3: For each group with candidates, run SPOA + concordance
        total_assigned = 0

        for group_id, candidates in group_candidates.items():
            group_indices = identity_groups[group_id]

            # Single-cluster group: assign directly via edit distance
            if len(group_indices) == 1:
                target_idx = group_indices[0]
                for rid, read_seq in candidates:
                    reassign_read(rid, target_idx)
                    total_assigned += 1
                continue

            # Multi-cluster group: SPOA + variant concordance
            cluster_consensus_seqs: Dict[int, str] = {}
            for idx in group_indices:
                cons = consensuses.get(idx)
                if cons:
                    cluster_consensus_seqs[idx] = cons

            if len(cluster_consensus_seqs) <= 1:
                # Only one valid consensus — assign all to it
                target_idx = next(iter(cluster_consensus_seqs))
                for rid, read_seq in candidates:
                    reassign_read(rid, target_idx)
                    total_assigned += 1
                continue

            # Build SPOA input: consensus sequences + candidate reads
            spoa_input: Dict[str, str] = {}
            consensus_spoa_ids: Dict[str, int] = {}
            for idx, cons in cluster_consensus_seqs.items():
                spoa_id = f"__cons_{idx}__"
                spoa_input[spoa_id] = cons
                consensus_spoa_ids[spoa_id] = idx

            # Subsample candidates if too many
            max_candidates = 500
            if len(candidates) > max_candidates:
                # Sort by quality descending
                candidates_with_qual = []
                for rid, seq in candidates:
                    rec = self.records.get(rid)
                    if rec:
                        qual = statistics.mean(rec.letter_annotations["phred_quality"])
                    else:
                        qual = 0.0
                    candidates_with_qual.append((qual, rid, seq))
                candidates_with_qual.sort(reverse=True)
                candidates = [(rid, seq) for _, rid, seq in candidates_with_qual[:max_candidates]]

            for rid, seq in candidates:
                spoa_input[rid] = seq

            # Run SPOA
            msa_result = _run_spoa_for_cluster_worker(spoa_input, self.disable_homopolymer_equivalence)
            if not msa_result:
                continue

            spoa_consensus_aligned = _extract_consensus_aligned_from_msa(msa_result.msa_string)
            if not spoa_consensus_aligned:
                continue

            alignment_by_id = {a.read_id: a for a in msa_result.alignments}

            # Find variant columns between cluster consensuses
            variant_columns = _find_variant_columns(alignment_by_id, consensus_spoa_ids,
                                                    spoa_consensus_aligned)
            if not variant_columns:
                # No variant positions — assign to closest by edit distance
                for rid, read_seq in candidates:
                    best_idx = min(cluster_consensus_seqs.keys(),
                                   key=lambda i: edlib.align(read_seq, cluster_consensus_seqs[i])['editDistance'])
                    reassign_read(rid, best_idx)
                    total_assigned += 1
                continue

            # Extract cluster consensus bases at variant columns
            cluster_bases: Dict[int, Dict[int, str]] = {}
            for spoa_id, idx in consensus_spoa_ids.items():
                aln = alignment_by_id.get(spoa_id)
                if aln:
                    cluster_bases[idx] = _extract_bases_at_columns(
                        aln, variant_columns, spoa_consensus_aligned)

            # Score each candidate and assign
            for rid, read_seq in candidates:
                aln = alignment_by_id.get(rid)
                if not aln:
                    continue

                read_bases = _extract_bases_at_columns(aln, variant_columns, spoa_consensus_aligned)

                # Concordance: count matches at variant positions for each cluster
                best_idx = None
                best_score = -1
                for idx, cbases in cluster_bases.items():
                    score = sum(1 for col in variant_columns
                                if read_bases.get(col) == cbases.get(col))
                    if score > best_score:
                        best_score = score
                        best_idx = idx

                if best_idx is None:
                    continue

                # Confirm via edit distance: read must be within group threshold
                target_cons = cluster_consensus_seqs.get(best_idx, '')
                if target_cons:
                    dist = edlib.align(read_seq, target_cons)['editDistance']
                    norm_dist = dist / max(len(read_seq), len(target_cons))
                    if norm_dist > max_distance:
                        continue

                reassign_read(rid, best_idx)
                total_assigned += 1

        if total_assigned > 0:
            logging.info(f"Discard reassignment: recovered {total_assigned} of "
                        f"{total_candidates} candidates ({rejected} reads below threshold)")
        else:
            logging.info(f"Discard reassignment: no reads recovered "
                        f"({total_candidates} candidates, {rejected} below threshold)")

        return subclusters

    def _run_cer_validation(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4c: Validate clusters via pairwise CER testing.

        Groups clusters by adjusted identity, then within each group:
        - The largest cluster (anchor) auto-passes
        - Each non-anchor is tested pairwise against validated clusters
        - CER = min(pairwise p*) using K = variant differences
        - Clusters failing CER (p* < assumed_error_rate) are marked 'ns'
        - All clusters are kept in output with CER annotations

        Annotates each cluster dict with 'cer_status' and 'cer_value'.
        """
        if self.assumed_error_rate <= 0:
            return subclusters

        if len(subclusters) <= 1:
            if subclusters:
                subclusters[0]['cer_status'] = 'anchor'
                subclusters[0]['cer_factor'] = None
                subclusters[0]['cer_pstar'] = None
                subclusters[0]['cer_details'] = None
                subclusters[0]['identity_group_id'] = 'g0'
                subclusters[0]['cer_group_N'] = len(subclusters[0].get('read_ids', []))
            return subclusters

        # Generate consensus for each cluster with ambiguity calling
        # so CER ignores positions where reads disagree
        consensuses = {}
        for i, cluster_dict in enumerate(subclusters):
            consensuses[i] = self._generate_consensus_for_validation(
                cluster_dict['read_ids'], apply_ambiguity_calling=True)

        # Form identity groups
        identity_groups = self._form_identity_groups(subclusters, consensuses)

        # Tag each subcluster with its identity_group_id and group_N
        # for downstream metadata reporting.
        for group_index, (_, group_indices) in enumerate(sorted(identity_groups.items())):
            group_id = f"g{group_index}"
            group_N = sum(len(subclusters[i].get('read_ids', [])) for i in group_indices)
            for i in group_indices:
                subclusters[i]['identity_group_id'] = group_id
                subclusters[i]['cer_group_N'] = group_N

        # Validate within each group
        validated_count = 0
        ns_count = 0
        for group_indices in identity_groups.values():
            v, n = self._validate_identity_group(subclusters, consensuses, group_indices)
            validated_count += v
            ns_count += n

        logging.info(f"CER validation: {validated_count} passed, {ns_count} not significant "
                     f"({len(identity_groups)} identity group(s))")

        return subclusters

    def _generate_consensus_for_validation(self, read_ids: Set[str],
                                            apply_ambiguity_calling: bool = False) -> Optional[str]:
        """Generate consensus sequence via SPOA.

        Args:
            read_ids: Set of read IDs to include
            apply_ambiguity_calling: If True, apply IUPAC ambiguity calling so
                positions where reads disagree are marked with IUPAC codes.
                Used for CER validation where ambiguous positions should be
                treated as 'no signal' for inter-cluster distinction.
                Must NOT be used when the consensus will be input to SPOA
                (e.g., reassignment), since SPOA doesn't understand IUPAC.
        """
        if not read_ids:
            return None
        if len(read_ids) == 1:
            rid = next(iter(read_ids))
            return self.sequences.get(rid)

        # Sort by quality descending, sample if needed
        sorted_ids = sorted(read_ids,
            key=lambda x: (-statistics.mean(self.records[x].letter_annotations["phred_quality"]), x))
        if len(sorted_ids) > self.max_sample_size:
            sorted_ids = sorted_ids[:self.max_sample_size]

        seqs = {sid: self.sequences[sid] for sid in sorted_ids}
        result = _run_spoa_for_cluster_worker(seqs, self.disable_homopolymer_equivalence)
        if not result or not result.consensus:
            return None

        if apply_ambiguity_calling and self.enable_iupac_calling:
            consensus, _, _ = call_iupac_ambiguities(
                consensus=result.consensus,
                alignments=result.alignments,
                msa_to_consensus_pos=result.msa_to_consensus_pos,
                min_variant_frequency=self.min_ambiguity_frequency,
                min_variant_count=self.min_ambiguity_count
            )
            return consensus

        return result.consensus

    def _form_identity_groups(self, subclusters: List[Dict],
                              consensuses: Dict[int, Optional[str]]) -> Dict[int, List[int]]:
        """Group clusters by pairwise adjusted identity using union-find."""
        n = len(subclusters)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        threshold = 1.0 - self.group_identity
        for i in range(n):
            for j in range(i + 1, n):
                ci, cj = consensuses.get(i), consensuses.get(j)
                if ci and cj:
                    dist, _ = _hp_normalized_pairwise_compare(ci, cj)
                    if dist <= threshold:
                        union(i, j)

        groups: Dict[int, List[int]] = defaultdict(list)
        for i in range(n):
            groups[find(i)].append(i)
        return dict(groups)

    def _validate_identity_group(self, subclusters: List[Dict],
                                 consensuses: Dict[int, Optional[str]],
                                 group_indices: List[int]) -> Tuple[int, int]:
        """Validate clusters within an identity group via context-aware pairwise CER.

        Each candidate cluster is compared pairwise against all already-validated
        clusters in the group. Each pairwise comparison classifies the differing
        positions, looks up empirical q_ctx values, and computes a CER factor
        (per-position multiplicative inflation needed for the variant to be
        plausible artifact). The candidate's reported factor is the *minimum*
        across all pairwise comparisons (the nearest plausible artifact source).

        A candidate passes when its minimum factor meets the configured factor
        threshold. Setting --assumed-error-rate <= 0 disables the CER gate.

        Returns (validated_count, ns_count).
        """
        # Sort by size descending
        sorted_indices = sorted(group_indices,
            key=lambda i: len(subclusters[i]['read_ids']), reverse=True)

        # N = total reads in this identity group
        group_N = sum(len(subclusters[i]['read_ids']) for i in sorted_indices)

        # Anchor: largest cluster auto-passes
        anchor_idx = sorted_indices[0]
        subclusters[anchor_idx]['cer_status'] = 'anchor'
        subclusters[anchor_idx]['cer_factor'] = None
        subclusters[anchor_idx]['cer_pstar'] = None
        subclusters[anchor_idx]['cer_details'] = None

        if len(sorted_indices) == 1:
            return 1, 0

        validated = [anchor_idx]
        validated_count = 1
        ns_count = 0
        factor_threshold = self.cer_factor_threshold

        for candidate_idx in sorted_indices[1:]:
            candidate_consensus = consensuses.get(candidate_idx)
            candidate_M = len(subclusters[candidate_idx]['read_ids'])

            if not candidate_consensus:
                subclusters[candidate_idx]['cer_status'] = 'ns'
                subclusters[candidate_idx]['cer_factor'] = 0.0
                subclusters[candidate_idx]['cer_pstar'] = None
                subclusters[candidate_idx]['cer_details'] = None
                ns_count += 1
                continue

            best = self._compare_candidate_against_validated(
                candidate_consensus=candidate_consensus,
                candidate_M=candidate_M,
                group_N=group_N,
                validated=validated,
                consensuses=consensuses,
            )

            if best is None:
                # No valid pairwise comparisons (all K=0 or alignment failed
                # or all positions had unsupported q_ctx). Auto-pass.
                subclusters[candidate_idx]['cer_status'] = 'pass'
                subclusters[candidate_idx]['cer_factor'] = None
                subclusters[candidate_idx]['cer_pstar'] = None
                subclusters[candidate_idx]['cer_details'] = None
                validated.append(candidate_idx)
                validated_count += 1
                continue

            min_factor, pstar, details = best
            subclusters[candidate_idx]['cer_factor'] = min_factor
            subclusters[candidate_idx]['cer_pstar'] = pstar
            subclusters[candidate_idx]['cer_details'] = details

            if factor_threshold <= 0 or min_factor >= factor_threshold:
                subclusters[candidate_idx]['cer_status'] = 'pass'
                validated.append(candidate_idx)
                validated_count += 1
            else:
                subclusters[candidate_idx]['cer_status'] = 'ns'
                ns_count += 1
                logging.debug(
                    f"CER ns: cluster with {candidate_M} reads, "
                    f"factor={min_factor:.3f} < {factor_threshold}"
                )

        return validated_count, ns_count

    def _compare_candidate_against_validated(
        self,
        candidate_consensus: str,
        candidate_M: int,
        group_N: int,
        validated: List[int],
        consensuses: Dict[int, Optional[str]],
    ) -> Optional[Tuple[float, Optional[float], dict]]:
        """Compute the worst-case (minimum) CER factor for a candidate.

        Walks each validated reference consensus, classifies the pairwise
        differences with context-aware tags, looks up q_ctx for each, and
        computes the per-position CER factor for that pair. The pair that
        produces the smallest factor (the nearest plausible artifact source)
        determines the candidate's reported metrics.

        Returns:
            (min_factor, per_position_pstar, details_dict) for the worst pair,
            or None if no pair produced a valid evaluation. The details dict
            contains the K, context tags, q_ctx values, and reference index
            for the determining pair, suitable for both FASTA annotation and
            metadata reproducibility.
        """
        best: Optional[Tuple[float, Optional[float], dict]] = None

        for ref_idx in validated:
            ref_consensus = consensuses.get(ref_idx)
            if not ref_consensus:
                continue

            tags = classify_pairwise_differences(candidate_consensus, ref_consensus)
            if tags is None or not tags:
                # K=0 (identical/IUPAC-compat) or alignment failure
                continue

            qctx_values: List[float] = []
            kept_tags = []
            for tag in tags:
                q = get_qctx(tag, table=self.qctx_table)
                if q is None:
                    continue
                qctx_values.append(q)
                kept_tags.append(tag)

            if not qctx_values:
                # All positions had unsupported q_ctx (e.g., HP L>=6 only).
                continue

            K = len(qctx_values)
            L = max(len(candidate_consensus), len(ref_consensus))

            factor = compute_cer_factor(
                N=group_N, M=candidate_M, n_sites=L,
                q_ctx_per_position=qctx_values,
                alpha=self.significance_level,
            )
            if factor is None:
                continue

            if best is None or factor < best[0]:
                pstar = compute_per_position_qstar(
                    N=group_N, M=candidate_M, n_sites=L, K=K,
                    alpha=self.significance_level,
                )
                details = {
                    'K': K,
                    'tags': [t.to_string() for t in kept_tags],
                    'q_ctx': qctx_values,
                    'ref_idx': ref_idx,
                }
                best = (factor, pstar, details)

        return best

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
            # Exclude "ns" clusters from ratio denominator
            validated_sizes = [len(c['read_ids']) for c in large_clusters
                               if c.get('cer_status') != 'ns']
            if validated_sizes:
                largest_size = max(validated_sizes)
            else:
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

        total_sequences = self.total_input_reads
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
        cluster_dicts_by_idx = {}
        for final_idx, cluster_dict in enumerate(clusters, 1):
            cluster = cluster_dict['read_ids']
            cluster_dict['cluster_id'] = f"c{final_idx}"
            cluster_dicts_by_idx[final_idx] = cluster_dict
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

                    # Look up CER from validation pass
                    cluster_dict = cluster_dicts_by_idx.get(final_idx, {})
                    cluster_cer_factor = cluster_dict.get('cer_factor')
                    cluster_cer_details_dict = cluster_dict.get('cer_details')
                    cluster_cer_details = format_cer_details(
                        cluster_dict.get('cer_pstar'),
                        cluster_cer_details_dict,
                    )

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
                        iupac_count=iupac_count,
                        cer_factor=cluster_cer_factor,
                        cer_details=cluster_cer_details,
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
            4b. Read reassignment (concordance-based, within identity groups)
            4c. CER validation (pairwise significance testing within identity groups)
            5. Filtering (size and ratio thresholds)
            6. Output generation
            7. Write discarded reads (optional)

        Args:
            algorithm: Clustering algorithm to use ('graph' for MCL or 'greedy')
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            # Phase 1: Initial clustering
            initial_clusters = self._run_initial_clustering(temp_dir, algorithm)

            # Track reads not assigned to any cluster
            clustered_reads = set().union(*initial_clusters) if initial_clusters else set()
            unclustered = set(self.sequences.keys()) - clustered_reads
            if unclustered:
                self.discarded_read_ids.update(unclustered)

            # Hard floor: drop clusters with < 3 reads (no error correction possible)
            min_consensus_reads = 3
            small = [c for c in initial_clusters if len(c) < min_consensus_reads]
            if small:
                for c in small:
                    self.discarded_read_ids.update(c)
                initial_clusters = [c for c in initial_clusters if len(c) >= min_consensus_reads]
                logging.info(f"Dropped {len(small)} clusters with < {min_consensus_reads} reads "
                            f"({sum(len(c) for c in small)} reads)")

            # Phase 2: Pre-phasing merge
            merged_clusters = self._run_prephasing_merge(initial_clusters)

            # Phase 2b: Early filtering (optional)
            clusters_to_phase, early_filtered = self._apply_early_filter(merged_clusters)

            # Phase 3: Variant detection + phasing
            all_subclusters = self._run_variant_phasing(clusters_to_phase)

            # Phase 4: Post-phasing merge
            merged_subclusters = self._run_postphasing_merge(all_subclusters)

            # Phase 4a: Filter noisy small clusters
            cleaned_subclusters = self._filter_noisy_clusters(merged_subclusters)

            # Phase 4b: Read reassignment
            reassigned_subclusters = self._run_read_reassignment(cleaned_subclusters)

            # Phase 4b2: Discard reassignment (recover reads dropped earlier)
            discard_reassigned = self._run_discard_reassignment(reassigned_subclusters)

            # Phase 4b3: Second phasing pass (split variants introduced by reassignment/recovery)
            rephased_subclusters = self._run_second_phasing_pass(discard_reassigned)

            # Phase 4c: CER validation
            validated_subclusters = self._run_cer_validation(rephased_subclusters)

            # Phase 5: Size filtering
            filtered_clusters = self._run_size_filtering(validated_subclusters)

            # Capture for metadata serialization (write_metadata reads this).
            self.final_cluster_dicts = filtered_clusters

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
            min_variant_count=self.min_variant_count,
            total_specimen_reads=self.total_input_reads,
            assumed_error_rate=self.assumed_error_rate,
            significance_level=self.significance_level,
        )

        # Build qualities dict for consistent SPOA ordering
        qualities = {}
        for rid in cluster_read_ids:
            if rid in self.records:
                qualities[rid] = statistics.mean(self.records[rid].letter_annotations["phred_quality"])

        return _phase_reads_by_variants_standalone(
            cluster_read_ids, self.sequences, qualities, variant_positions, config
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

    def are_homopolymer_equivalent(self, seq1: str, seq2: str) -> bool:
        """Check if two sequences differ only in HP run lengths.

        Uses edlib global alignment with consensus-derived HP context.
        Only HP runs of length >= min_hp_length (default 6) are normalized.
        Any non-HP difference (substitution, non-HP indel, terminal overhang)
        means the sequences are NOT equivalent.
        """
        if not seq1 or not seq2:
            return seq1 == seq2

        _, is_equiv = _hp_normalized_pairwise_compare(seq1, seq2)
        return is_equiv

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
