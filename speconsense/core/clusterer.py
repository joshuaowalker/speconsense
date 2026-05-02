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
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

try:
    from speconsense import __version__
except ImportError:
    __version__ = "dev"

from speconsense.msa import ReadAlignment, analyze_positional_variation, has_no_majority, call_iupac_ambiguities, compute_cluster_err_factor
from speconsense.context import classify_pairwise_differences
from speconsense.distances import count_variant_differences
from speconsense.qctx import DORADO_V5_0, get_qctx, load_table as load_error_model
from speconsense.significance import (
    compute_cer_factor,
    compute_critical_error_rate,
    compute_per_position_qstar,
)
from speconsense.scalability import (
    VsearchCandidateFinder,
    ScalablePairwiseOperation,
    ScalabilityConfig,
)

from .workers import (
    ClusterProcessingConfig,
    ConsensusGenerationConfig,
    MSAResult,
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


def _count_independent_sites(seq: str) -> int:
    """Bonferroni site count for context-aware CER (cer_in_practice §2.3).

    Each maximal run of identical bases counts as one site regardless of run
    length: a non-HP position is a length-1 run (1 site), and an HP run of
    length L>=2 is also 1 site (not L), since an HP-length variant tests one
    hypothesis per run rather than one per base. Equivalent to counting base
    transitions plus one.
    """
    if not seq:
        return 0
    n = 1
    for i in range(1, len(seq)):
        if seq[i] != seq[i - 1]:
            n += 1
    return n


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
                 enable_secondpass_phasing: bool = True,
                 enable_read_reassignment: bool = True,
                 enable_discard_recovery: bool = True,
                 min_variant_frequency: float = 0.10,
                 min_variant_count: int = 3,
                 min_ambiguity_frequency: float = 0.10,
                 min_ambiguity_count: int = 3,
                 enable_iupac_calling: bool = True,
                 scale_threshold: int = 1001,
                 max_threads: int = 1,
                 collect_discards: bool = False,
                 significance_level: float = 1e-5,
                 group_identity: float = 0.85,
                 min_hp_length: int = 6,
                 error_model: str = "dorado-v5.0"):
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

        self.enable_secondpass_phasing = enable_secondpass_phasing
        self.enable_read_reassignment = enable_read_reassignment
        self.enable_discard_recovery = enable_discard_recovery
        self.min_variant_frequency = min_variant_frequency
        self.min_variant_count = min_variant_count
        self.min_ambiguity_frequency = min_ambiguity_frequency
        self.min_ambiguity_count = min_ambiguity_count
        self.enable_iupac_calling = enable_iupac_calling
        self.scale_threshold = scale_threshold
        self.max_threads = max_threads
        self.collect_discards = collect_discards
        self.significance_level = significance_level
        self.group_identity = group_identity
        # HP run length threshold for MSA variant detection. Runs of length
        # >= min_hp_length are treated as HP context (length variants in them
        # are suppressed at phasing time). Runs of length < min_hp_length
        # surface length variants as regular candidates, which the
        # context-aware CER framework scores via q_ctx lookup.
        self.min_hp_length = min_hp_length
        # Context-aware CER configuration. Every non-anchor candidate is
        # compared pairwise against the accumulated validated pool and
        # annotated with its per-position cer_factor. Summarize applies the
        # user-visible pass/ns decision via --min-cer-factor.
        self.error_model = error_model
        self.qctx_table = load_error_model(error_model)
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

    def _log_stage(self, description: str, clusters,
                   reads_in_clusters: Optional[int] = None) -> None:
        """Emit one INFO line: fixed-width running cluster/read totals (left) + stage description (right).

        Accepts either a list of cluster dicts (with 'read_ids') or a list of
        read-id collections (sets/lists).
        """
        n_clusters = len(clusters)
        if reads_in_clusters is None:
            n_reads = sum(
                len(c['read_ids']) if isinstance(c, dict) else len(c)
                for c in clusters
            )
        else:
            n_reads = reads_in_clusters
        total = getattr(self, 'total_input_reads', 0) or 0
        pct = 100.0 * n_reads / total if total > 0 else 0.0
        width = max(len(str(total)), 1)
        logging.info(
            f"[{n_clusters:>3} clusters, {n_reads:>{width}}/{total} reads, "
            f"{pct:>5.1f}%]  {description}"
        )

    def write_metadata(self) -> None:
        """Write run metadata to JSON file for use by post-processing tools.

        Includes per-cluster CER reproduction data when clustering has
        completed (self.final_cluster_dicts populated). The schema follows the
        CER-in-practice paper Appendix A: parameters identify the error model
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
                "enable_secondpass_phasing": self.enable_secondpass_phasing,
                "enable_read_reassignment": self.enable_read_reassignment,
                "enable_discard_recovery": self.enable_discard_recovery,
                "min_variant_frequency": self.min_variant_frequency,
                "min_variant_count": self.min_variant_count,
                "min_ambiguity_frequency": self.min_ambiguity_frequency,
                "min_ambiguity_count": self.min_ambiguity_count,
                "enable_iupac_calling": self.enable_iupac_calling,
                "scale_threshold": self.scale_threshold,
                "max_threads": self.max_threads,
                "orient_mode": self.orient_mode,
                "significance_level": self.significance_level,
                "group_identity": self.group_identity,
                "error_model": self.error_model,
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
    def _assign_identity_ranks(cluster_dicts: List[Dict]) -> None:
        """Stamp identity_group_rank and identity_variant_rank on clusters.

        Groups are ranked by anchor (largest member) size descending. Variants
        within each group are ranked by size descending. Clusters without an
        ``identity_group_id`` are skipped (they won't emit gid/vid).
        """
        group_anchor_size: Dict[str, int] = {}
        for cluster_dict in cluster_dicts:
            gid = cluster_dict.get('identity_group_id')
            if gid is None:
                continue
            size = len(cluster_dict.get('read_ids', []))
            group_anchor_size[gid] = max(group_anchor_size.get(gid, 0), size)

        group_rank_map = {
            gid: rank for rank, (gid, _) in enumerate(
                sorted(group_anchor_size.items(), key=lambda kv: kv[1], reverse=True),
                start=1,
            )
        }
        for cluster_dict in cluster_dicts:
            gid = cluster_dict.get('identity_group_id')
            if gid is not None:
                cluster_dict['identity_group_rank'] = group_rank_map[gid]

        group_members: Dict[str, List[Dict]] = defaultdict(list)
        for cluster_dict in cluster_dicts:
            gid = cluster_dict.get('identity_group_id')
            if gid is not None:
                group_members[gid].append(cluster_dict)
        for members in group_members.values():
            members.sort(key=lambda c: len(c.get('read_ids', [])), reverse=True)
            for rank, cluster_dict in enumerate(members, start=1):
                cluster_dict['identity_variant_rank'] = rank

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
        ef_details = cluster_dict.get('err_factor_details') or {}
        return {
            'cluster_id': cluster_dict.get('cluster_id'),
            'identity_group': cluster_dict.get('identity_group_id'),
            'group_rank': cluster_dict.get('identity_group_rank'),
            'variant_rank': cluster_dict.get('identity_variant_rank'),
            'M': len(cluster_dict.get('read_ids', [])),
            'N': cluster_dict.get('cer_group_N'),
            'K': details.get('K') if details else None,
            'n_sites': details.get('n_sites') if details else None,
            'context_tags': details.get('tags') if details else None,
            'q_ctx_per_position': details.get('q_ctx') if details else None,
            'compared_against_idx': details.get('ref_idx') if details else None,
            'cer_factor': cluster_dict.get('cer_factor'),
            'cer_pstar': cluster_dict.get('cer_pstar'),
            'err_factor': cluster_dict.get('err_factor'),
            'err_factor_obs_sum': ef_details.get('obs_sum'),
            'err_factor_exp_sum': ef_details.get('exp_sum'),
            'err_factor_cols': ef_details.get('cols'),
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
            logging.info(f"Presampling {self.presample_size} reads from {len(all_records)} total "
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

            logging.info(f"Presampled {len(presampled)} reads "
                         f"({primary_to_include} primary, {len(presampled) - primary_to_include} augmented)")
            all_records = presampled

        # Add all selected records to internal storage
        for record in all_records:
            self.sequences[record.id] = str(record.seq)
            self.records[record.id] = record

        # Log scalability mode status for large datasets
        if len(self.sequences) >= self.scale_threshold and self.scale_threshold > 0:
            if self._candidate_finder is not None:
                logging.info(f"Scalability mode active for {len(self.sequences)} reads (threshold: {self.scale_threshold})")
            else:
                logging.warning(f"Dataset has {len(self.sequences)} sequences (>= threshold {self.scale_threshold}) "
                               "but vsearch not found. Using brute-force.")

        self.total_input_reads = len(self.sequences)

    def _get_scalable_operation(self, scoring_function=None) -> ScalablePairwiseOperation:
        """Get a ScalablePairwiseOperation for pairwise comparisons.

        Args:
            scoring_function: Optional override for the (seq1, seq2, id1, id2)
                -> similarity scoring callable. Defaults to
                self.calculate_similarity (the K-NN / equivalence-grouping
                metric). Pass an alternative when the call site needs a
                different distance semantics (e.g. HP-normalized identity
                grouping).
        """
        if scoring_function is None:
            scoring_function = lambda seq1, seq2, id1, id2: self.calculate_similarity(seq1, seq2)
        return ScalablePairwiseOperation(
            candidate_finder=self._candidate_finder,
            scoring_function=scoring_function,
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
                    (cluster_idx, sampled_seqs, self.disable_homopolymer_equivalence, self.min_hp_length)
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
                    _, result = _run_spoa_worker((cluster_idx, sampled_seqs, self.disable_homopolymer_equivalence, self.min_hp_length))
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
            self._log_stage(
                f"{phase_name}: combined {len(clusters)}→{len(merged)} ({len(clusters) - len(merged)} {merge_type})",
                merged,
            )
        else:
            logging.debug(f"{phase_name}: no merges ({merge_type} not found)")

        return merged

    def _find_root(self, merged_to: List[int], i: int) -> int:
        """Find the root index of a merged cluster using path compression."""
        if merged_to[i] != i:
            merged_to[i] = self._find_root(merged_to, merged_to[i])
        return merged_to[i]

    def write_cluster_files(self, cluster_name: str, cluster: Set[str],
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
                            err_factor: Optional[float] = None,
                            identity_group_rank: Optional[int] = None,
                            identity_variant_rank: Optional[int] = None) -> None:
        """Write cluster files: reads FASTQ, MSA, and consensus FASTA.

        ``cluster_name`` is the group.variant suffix used in filenames and the
        FASTA header (e.g., ``"1.v2"``). It matches the core-assigned identity
        group/variant rank so summarize can honor the same naming end-to-end.

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
        if err_factor is not None:
            info_parts.append(f"err_factor={err_factor:.3f}")
        if identity_group_rank is not None:
            info_parts.append(f"gid={identity_group_rank}")
        if identity_variant_rank is not None:
            info_parts.append(f"vid={identity_variant_rank}")
        info_str = " ".join(info_parts)

        # Write reads FASTQ to debug directory with new naming convention
        # Use sorted order (by quality descending) if available, matching MSA order
        reads_file = os.path.join(self.debug_dir, f"{self.sample_name}-{cluster_name}-RiC{ric_size}-reads.fastq")
        with open(reads_file, 'w') as f:
            read_ids_to_write = sorted_cluster_ids if sorted_cluster_ids is not None else cluster
            for seq_id in read_ids_to_write:
                SeqIO.write(self.records[seq_id], f, "fastq")

        # Write sampled reads FASTQ (only sequences used for consensus generation)
        # Use sorted order (by quality descending) if available, matching MSA order
        if sampled_ids is not None or sorted_sampled_ids is not None:
            sampled_file = os.path.join(self.debug_dir, f"{self.sample_name}-{cluster_name}-RiC{ric_size}-sampled.fastq")
            with open(sampled_file, 'w') as f:
                sampled_to_write = sorted_sampled_ids if sorted_sampled_ids is not None else sampled_ids
                for seq_id in sampled_to_write:
                    SeqIO.write(self.records[seq_id], f, "fastq")

        # Write MSA (multiple sequence alignment) to debug directory
        if msa is not None:
            msa_file = os.path.join(self.debug_dir, f"{self.sample_name}-{cluster_name}-RiC{ric_size}-msa.fasta")
            with open(msa_file, 'w') as f:
                f.write(msa)

        # Write untrimmed consensus to debug directory
        with open(os.path.join(self.debug_dir, f"{self.sample_name}-{cluster_name}-RiC{ric_size}-untrimmed.fasta"),
                  'w') as f:
            f.write(f">{self.sample_name}-{cluster_name} {info_str}\n")
            f.write(consensus + "\n")

        # Write consensus to main output file if handle is provided
        if consensus_fasta_handle:
            final_consensus = trimmed_consensus if trimmed_consensus else consensus
            consensus_fasta_handle.write(f">{self.sample_name}-{cluster_name} {info_str}\n")
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

        logging.debug(f"Running MCL algorithm with inflation {self.inflation}...")
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
        logging.debug("Running greedy clustering algorithm...")

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
        logging.debug("Calculating pairwise sequence similarities...")

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
        self._log_stage("Initial clustering", initial_clusters, reads_in_clusters=clustered_count)
        return initial_clusters

    def _run_prephasing_merge(self, initial_clusters: List[Set[str]]) -> List[Set[str]]:
        """Phase 2: Pre-phasing merge — combine HP-equivalent initial clusters.

        Maximizes read depth for variant detection in the phasing phase.

        Args:
            initial_clusters: List of initial clusters from Phase 1

        Returns:
            List of merged clusters (sets of read IDs)
        """
        if self.disable_cluster_merging:
            logging.debug("Cluster merging disabled, skipping cluster equivalence merge")
            return initial_clusters

        # Convert initial clusters to dict format for merge_similar_clusters
        initial_cluster_dicts = [
            {'read_ids': cluster, 'initial_cluster_num': i, 'allele_combo': None}
            for i, cluster in enumerate(initial_clusters, 1)
        ]
        merged_dicts = self.merge_similar_clusters(initial_cluster_dicts, phase_name="Cluster equivalence merge")
        # Extract back to sets for Phase 3
        return [d['read_ids'] for d in merged_dicts]

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
            enable_secondpass_phasing=self.enable_secondpass_phasing,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            min_variant_frequency=self.min_variant_frequency,
            min_variant_count=self.min_variant_count,
            significance_level=self.significance_level,
            min_hp_length=self.min_hp_length,
            max_sample_size=self.max_sample_size,
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

        if split_count > 0 or len(all_subclusters) != len(merged_clusters):
            self._log_stage(
                f"Variant phasing: split {split_count} cluster(s) → {len(all_subclusters)} sub-clusters",
                all_subclusters,
            )
        else:
            logging.debug("Variant phasing: no splits")
        return all_subclusters

    def _run_postphasing_merge(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 4: Post-phasing merge — combine HP-equivalent subclusters.

        Args:
            subclusters: List of subclusters from Phase 3

        Returns:
            List of merged subclusters
        """
        if self.disable_cluster_merging:
            logging.debug("Cluster merging disabled, skipping cluster equivalence merge (round 2)")
            return subclusters

        return self.merge_similar_clusters(subclusters, phase_name="Cluster equivalence merge (round 2)")

    def _filter_noisy_clusters(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 5: Noise filter — remove unreliable small clusters with no-majority positions.

        Only applies to clusters below the phasing floor (< min_variant_count * 2),
        since larger clusters have enough reads for reliable consensus. Disbanded
        reads go to discards for recovery by the global discard reassignment pass.

        SPOA dispatch parallels via ProcessPoolExecutor when max_threads > 1 and
        more than 10 small clusters are queued; otherwise the same worker is
        invoked sequentially (deterministic-equivalent fallback).
        """
        phasing_floor = self.min_variant_count * 2
        result: List[Optional[Dict]] = []
        total_disbanded = 0
        disbanded_reads = 0

        # Pass 1: partition. Large clusters pass through; small ones queue for SPOA.
        # `result` carries None placeholders at small-cluster slots so we can
        # stitch SPOA outcomes back in input order.
        small_queue: List[Tuple[int, Dict[str, str]]] = []
        small_meta: List[Tuple[int, Dict, Set[str]]] = []  # (result_idx, cluster_dict, read_ids)

        for cluster_dict in subclusters:
            read_ids = cluster_dict['read_ids']
            if len(read_ids) >= phasing_floor:
                result.append(cluster_dict)
                continue

            qualities = {
                sid: statistics.mean(self.records[sid].letter_annotations["phred_quality"])
                for sid in read_ids
            }
            sorted_ids = sorted(read_ids, key=lambda x: (-qualities.get(x, 0), x))
            cluster_seqs = {sid: self.sequences[sid] for sid in sorted_ids}

            placeholder_idx = len(result)
            result.append(None)
            queue_idx = len(small_queue)
            small_queue.append((queue_idx, cluster_seqs))
            small_meta.append((placeholder_idx, cluster_dict, read_ids))

        # Pass 2: dispatch SPOA in parallel when worthwhile.
        msa_results: Dict[int, Optional[MSAResult]] = {}
        if small_queue:
            work_packages = [
                (queue_idx, seqs, self.disable_homopolymer_equivalence, self.min_hp_length)
                for queue_idx, seqs in small_queue
            ]
            if self.max_threads > 1 and len(work_packages) > 10:
                from concurrent.futures import ProcessPoolExecutor
                with ProcessPoolExecutor(max_workers=self.max_threads) as executor:
                    spoa_results = list(tqdm(
                        executor.map(_run_spoa_worker, work_packages),
                        total=len(work_packages),
                        desc="Noise filter SPOA",
                    ))
            else:
                spoa_results = [_run_spoa_worker(pkg) for pkg in work_packages]

            for queue_idx, msa_result in spoa_results:
                msa_results[queue_idx] = msa_result

        # Pass 3: stitch results in input order. Positional analysis runs on the
        # parent process; only SPOA itself is parallel.
        for queue_idx, (placeholder_idx, cluster_dict, read_ids) in enumerate(small_meta):
            msa_result = msa_results.get(queue_idx)

            if msa_result is None:
                self.discarded_read_ids.update(read_ids)
                total_disbanded += 1
                disbanded_reads += len(read_ids)
                continue

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
                result[placeholder_idx] = cluster_dict

        # Drop None placeholders left by disbanded clusters.
        result = [c for c in result if c is not None]

        if total_disbanded > 0:
            self._log_stage(
                f"Noise filter: disbanded {total_disbanded} cluster(s) (-{disbanded_reads} reads)",
                result,
            )
        else:
            logging.debug("Noise filter: no clusters disbanded")

        return result

    def _run_read_reassignment(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 6: Read reassignment — move reads to best-matching clusters within identity groups.

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
            # Remove empty clusters
            subclusters = [c for c in subclusters if c['read_ids']]
            self._log_stage(
                f"Read reassignment: moved {total_reassigned} read(s)",
                subclusters,
            )
        else:
            logging.debug("Read reassignment: no reads moved")

        return subclusters

    def _run_second_phasing_pass(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 8: Second phasing pass — re-phase subclusters after read reassignment / discard recovery.

        Reassignment can move reads into clusters where they create detectable
        variants. This pass splits any such clusters. No outlier removal —
        just variant detection and phasing.
        """
        config = ClusterProcessingConfig(
            enable_secondpass_phasing=self.enable_secondpass_phasing,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            min_variant_frequency=self.min_variant_frequency,
            min_variant_count=self.min_variant_count,
            significance_level=self.significance_level,
            min_hp_length=self.min_hp_length,
            max_sample_size=self.max_sample_size,
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

            msa_result = _run_spoa_for_cluster_worker(cluster_seqs, self.disable_homopolymer_equivalence, self.min_hp_length)
            if msa_result is None:
                result.append(cluster_dict)
                continue

            # Sample = top-N-by-quality. Mirrors the final consensus sample so
            # a variant crossing threshold only on the sample still surfaces
            # here instead of silently becoming an IUPAC code later.
            sample_for_detect = None
            if config.max_sample_size is not None and len(sorted_ids) > config.max_sample_size:
                sample_for_detect = set(sorted_ids[:config.max_sample_size])

            variant_positions = _detect_variant_positions_standalone(
                msa_result.alignments, msa_result.consensus, msa_result.msa_to_consensus_pos,
                config.min_variant_frequency, config.min_variant_count,
                sample_read_ids=sample_for_detect,
            )

            cluster_label = f"cluster(n={len(read_ids)})"

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
            # Remove empty clusters
            result = [c for c in result if c['read_ids']]
            self._log_stage(
                f"Variant phasing (round 2): split {split_count} cluster(s) → {len(result)} sub-clusters",
                result,
            )
        else:
            logging.debug("Variant phasing (round 2): no splits")

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
            msa_result = _run_spoa_for_cluster_worker(spoa_input, self.disable_homopolymer_equivalence, self.min_hp_length)
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
        """Phase 7: Discard recovery — reassign previously-discarded reads to existing clusters.

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

        # Decide screening strategy. For large discard pools the dense
        # discards × clusters edlib scan dominates; replace it with a
        # vsearch top-K query, refined by edlib on the candidate slice
        # (same metric as the dense path).
        use_scalable_screen = (
            self.scalability_config.enabled
            and self._candidate_finder is not None
            and self._candidate_finder.is_available
            and len(self.discarded_read_ids) >= self.scalability_config.activation_threshold
        )

        if use_scalable_screen:
            # Index over cluster consensuses; query with discards.
            consensus_input = {f"__c_{i}__": c for i, c in valid_consensuses.items()}
            self._candidate_finder.build_index(consensus_input, self.output_dir)
            try:
                discard_seqs = {
                    rid: self.sequences[rid]
                    for rid in sorted(self.discarded_read_ids)
                    if rid in self.sequences
                }
                query_ids = list(discard_seqs.keys())
                relaxed_threshold = self.group_identity * self.scalability_config.relaxed_identity_factor
                top_k = 10
                vsearch_candidates = self._candidate_finder.find_candidates(
                    query_ids=query_ids,
                    sequences=discard_seqs,
                    min_identity=relaxed_threshold,
                    max_candidates=top_k,
                )

                logging.debug(
                    f"Discard screening (scalable): {len(query_ids)} discards, "
                    f"{len(valid_consensuses)} cluster consensuses, top_k={top_k}, "
                    f"dense edlib calls avoided={len(query_ids) * len(valid_consensuses)}"
                )

                cons_prefix_len = len("__c_")
                cons_suffix_len = len("__")
                for rid in query_ids:
                    seq = discard_seqs[rid]
                    cand_strs = vsearch_candidates.get(rid, [])
                    if not cand_strs:
                        rejected += 1
                        continue

                    best_idx = None
                    best_dist = float('inf')
                    for cand_str in cand_strs:
                        # Recover cluster index from "__c_<i>__" namespace.
                        i = int(cand_str[cons_prefix_len:-cons_suffix_len])
                        cons = valid_consensuses[i]
                        dist = edlib.align(seq, cons)['editDistance']
                        norm_dist = dist / max(len(seq), len(cons))
                        if norm_dist < best_dist:
                            best_dist = norm_dist
                            best_idx = i

                    if best_idx is None or best_dist > max_distance:
                        rejected += 1
                        continue

                    group_id = cluster_to_group[best_idx]
                    group_candidates[group_id].append((rid, seq))
            finally:
                self._candidate_finder.cleanup()
        else:
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
            logging.debug(f"Discard recovery: no candidates (all {rejected} reads below threshold)")
            return subclusters

        def reassign_read(rid, target_idx):
            subclusters[target_idx]['read_ids'].add(rid)
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
            msa_result = _run_spoa_for_cluster_worker(spoa_input, self.disable_homopolymer_equivalence, self.min_hp_length)
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
            self._log_stage(
                f"Discard recovery: recovered {total_assigned} of {total_candidates} candidate(s)",
                subclusters,
            )
        else:
            logging.debug(
                f"Discard recovery: no reads recovered "
                f"({total_candidates} candidates, {rejected} below threshold)"
            )

        return subclusters

    def _run_cer_validation(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 9: CER validation — annotate clusters with pairwise CER factors.

        Groups clusters by adjusted identity, then within each group:
        - The largest cluster (anchor) carries no pairwise comparison
        - Each non-anchor is tested pairwise against all larger clusters in
          the group; the minimum (worst-case) factor is reported
        - All clusters flow to output with cer_factor / cer_details annotations;
          summarize applies the user-visible pass/ns decision via --min-cer-factor
        """
        if len(subclusters) <= 1:
            if subclusters:
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

        annotated = 0
        for group_indices in identity_groups.values():
            annotated += self._validate_identity_group(subclusters, consensuses, group_indices)

        self._log_stage(
            f"Significance testing: {annotated} candidate(s), {len(identity_groups)} identity group(s)",
            subclusters,
        )

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
        result = _run_spoa_for_cluster_worker(seqs, self.disable_homopolymer_equivalence, self.min_hp_length)
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
        """Group clusters by pairwise adjusted identity using complete linkage.

        Every pair of clusters in a group must satisfy the identity threshold
        (complete linkage), preventing transitive chaining that can merge
        closely related but distinct variants in eDNA-style mixtures.

        For large cluster counts, the dense O(n²) distance matrix is replaced
        with a vsearch-derived sparse matrix: candidate pairs (those above the
        relaxed identity threshold) are scored under the same HP-normalized
        metric as the dense path; missing pairs default to 1.0 (max distance),
        which under complete linkage forces the corresponding clusters into
        separate groups — the desired semantics, since an absent pair is below
        the relaxed threshold and therefore below the user threshold.
        """
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform

        n = len(subclusters)
        if n == 0:
            return {}
        if n == 1:
            return {0: [0]}

        threshold = 1.0 - self.group_identity
        dist_matrix = [[0.0] * n for _ in range(n)]

        use_scalable = (
            self.scalability_config.enabled
            and self._candidate_finder is not None
            and self._candidate_finder.is_available
            and n > 50
            and n >= self.scalability_config.activation_threshold
        )

        if use_scalable:
            # Score function uses the SAME HP-normalized metric as the dense
            # path; vsearch identity is only used for candidate filtering.
            def hp_similarity(seq1, seq2, id1, id2):
                dist, _ = _hp_normalized_pairwise_compare(seq1, seq2, self.min_hp_length)
                return 1.0 - dist

            seq_dict = {str(i): c for i, c in consensuses.items() if c}
            op = self._get_scalable_operation(scoring_function=hp_similarity)
            sparse_distances = op.compute_distance_matrix(
                seq_dict, self.output_dir, min_identity=self.group_identity)

            logging.debug(
                f"Identity grouping (scalable): n={n}, candidate pairs={len(sparse_distances) // 2}, "
                f"dense pairs avoided={n * (n - 1) // 2}"
            )

            for i in range(n):
                for j in range(i + 1, n):
                    si, sj = str(i), str(j)
                    if si in seq_dict and sj in seq_dict:
                        d = sparse_distances.get((si, sj), 1.0)
                    else:
                        d = 1.0  # one or both consensuses None
                    dist_matrix[i][j] = d
                    dist_matrix[j][i] = d
        else:
            for i in range(n):
                for j in range(i + 1, n):
                    ci, cj = consensuses.get(i), consensuses.get(j)
                    if ci and cj:
                        dist, _ = _hp_normalized_pairwise_compare(ci, cj, self.min_hp_length)
                    else:
                        dist = 1.0
                    dist_matrix[i][j] = dist
                    dist_matrix[j][i] = dist

        condensed = squareform(dist_matrix, checks=False)
        Z = linkage(condensed, method='complete')
        labels = fcluster(Z, t=threshold, criterion='distance')

        groups: Dict[int, List[int]] = defaultdict(list)
        for idx, label in enumerate(labels):
            groups[int(label)].append(idx)
        return dict(groups)

    def _validate_identity_group(self, subclusters: List[Dict],
                                 consensuses: Dict[int, Optional[str]],
                                 group_indices: List[int]) -> int:
        """Annotate clusters within an identity group with pairwise CER factors.

        Each non-anchor candidate is compared pairwise against every larger
        cluster in the group. Each pairwise comparison classifies the differing
        positions, looks up empirical q_ctx values, and computes a CER factor
        (per-position multiplicative inflation needed for the variant to be
        plausible artifact). The candidate's reported factor is the *minimum*
        across all pairwise comparisons (the nearest plausible artifact source).

        For large identity groups (size > 50, scalability enabled), each
        candidate is compared against the top-K most similar larger peers
        rather than all larger peers. The minimum factor is dominated by the
        most-similar peer (smaller pairwise distance → smaller K positions →
        tighter Bonferroni → smaller factor), so top-K preserves the reported
        minimum with very high probability. min over a subset is always ≥ min
        over the full set, so the approximation only over-reports factors,
        never under-reports — i.e., approximation errors result in slightly
        less filtering, never over-filtering.

        Every candidate flows through to downstream phases. Summarize applies
        the user-visible pass/ns decision via --min-cer-factor.

        Returns the number of non-anchor candidates that received annotations.
        """
        # Sort by size descending; the largest cluster is the anchor (no
        # pairwise comparison). Subsequent candidates are compared against
        # all earlier ones.
        sorted_indices = sorted(group_indices,
            key=lambda i: len(subclusters[i]['read_ids']), reverse=True)

        group_N = sum(len(subclusters[i]['read_ids']) for i in sorted_indices)

        anchor_idx = sorted_indices[0]
        subclusters[anchor_idx]['cer_factor'] = None
        subclusters[anchor_idx]['cer_pstar'] = None
        subclusters[anchor_idx]['cer_details'] = None

        if len(sorted_indices) == 1:
            return 0

        # Decide whether to use within-group top-K. Threshold matches the
        # other vsearch-activation sites in the codebase (>50, e.g.
        # _form_identity_groups, merge_similar_clusters HP-equivalence).
        use_topk = (
            self.scalability_config.enabled
            and self._candidate_finder is not None
            and self._candidate_finder.is_available
            and len(sorted_indices) > 50
        )

        # In top-K mode, pre-compute each cluster's most-similar peers via a
        # single batched vsearch query over the group's consensuses. The
        # reference-pool filter inside the loop intersects this ranking with
        # the strict-larger-peer subset.
        #
        # K_query = 30 is 3× the K_use = 10 we actually consume, providing
        # headroom for the case where many of a candidate's most-similar
        # peers turn out to be smaller-than-self and therefore not in the
        # reference pool.
        K_use = 10
        K_query = 30
        topk_neighbors: Dict[int, List[int]] = {}
        if use_topk:
            group_seqs: Dict[str, str] = {}
            spoa_to_idx: Dict[str, int] = {}
            for idx in sorted_indices:
                cons = consensuses.get(idx)
                if cons:
                    spoa_id = f"__g_{idx}__"
                    group_seqs[spoa_id] = cons
                    spoa_to_idx[spoa_id] = idx

            finder = self._candidate_finder
            finder.build_index(group_seqs, self.output_dir)
            try:
                relaxed_threshold = self.group_identity * self.scalability_config.relaxed_identity_factor
                vsearch_candidates = finder.find_candidates(
                    query_ids=list(group_seqs.keys()),
                    sequences=group_seqs,
                    min_identity=relaxed_threshold,
                    max_candidates=K_query,
                )

                for spoa_id, idx in spoa_to_idx.items():
                    ranked: List[int] = []
                    for cand_str in vsearch_candidates.get(spoa_id, []):
                        cand_idx = spoa_to_idx.get(cand_str)
                        if cand_idx is not None and cand_idx != idx:
                            ranked.append(cand_idx)
                    topk_neighbors[idx] = ranked

                logging.debug(
                    f"CER within-group top-K: group_size={len(sorted_indices)}, "
                    f"K_query={K_query}, K_use={K_use}, "
                    f"dense pairs avoided≈{len(sorted_indices) * (len(sorted_indices) - 1) // 2 - len(sorted_indices) * K_use}"
                )
            finally:
                finder.cleanup()

        reference_pool = [anchor_idx]
        annotated = 0

        for candidate_idx in sorted_indices[1:]:
            candidate_consensus = consensuses.get(candidate_idx)
            candidate_M = len(subclusters[candidate_idx]['read_ids'])

            if not candidate_consensus:
                subclusters[candidate_idx]['cer_factor'] = 0.0
                subclusters[candidate_idx]['cer_pstar'] = None
                subclusters[candidate_idx]['cer_details'] = None
                reference_pool.append(candidate_idx)
                annotated += 1
                continue

            if use_topk:
                # Top-K most similar peers from the reference pool, preserving
                # vsearch's identity-descending order. Empty intersections
                # propagate to _compare_candidate_against_validated which
                # returns None — same handling as the "no valid pairwise" case.
                ref_set = set(reference_pool)
                validated_for_compare = [
                    i for i in topk_neighbors.get(candidate_idx, []) if i in ref_set
                ][:K_use]
            else:
                validated_for_compare = reference_pool

            best = self._compare_candidate_against_validated(
                candidate_consensus=candidate_consensus,
                candidate_M=candidate_M,
                group_N=group_N,
                validated=validated_for_compare,
                consensuses=consensuses,
            )

            if best is None:
                # No valid pairwise comparisons (K=0, alignment failed, or no
                # supported q_ctx). Leave factor unset.
                subclusters[candidate_idx]['cer_factor'] = None
                subclusters[candidate_idx]['cer_pstar'] = None
                subclusters[candidate_idx]['cer_details'] = None
            else:
                min_factor, pstar, details = best
                subclusters[candidate_idx]['cer_factor'] = min_factor
                subclusters[candidate_idx]['cer_pstar'] = pstar
                subclusters[candidate_idx]['cer_details'] = details

            reference_pool.append(candidate_idx)
            annotated += 1

        return annotated

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
            longer = (ref_consensus if len(ref_consensus) >= len(candidate_consensus)
                      else candidate_consensus)
            n_sites = _count_independent_sites(longer)

            factor = compute_cer_factor(
                N=group_N, M=candidate_M, n_sites=n_sites,
                q_ctx_per_position=qctx_values,
                alpha=self.significance_level,
            )
            if factor is None:
                continue

            if best is None or factor < best[0]:
                pstar = compute_per_position_qstar(
                    N=group_N, M=candidate_M, n_sites=n_sites, K=K,
                    alpha=self.significance_level,
                )
                details = {
                    'K': K,
                    'n_sites': n_sites,
                    'tags': [t.to_string() for t in kept_tags],
                    'q_ctx': qctx_values,
                    'ref_idx': ref_idx,
                }
                best = (factor, pstar, details)

        return best

    def _run_size_filtering(self, subclusters: List[Dict]) -> List[Dict]:
        """Phase 10: Size filtering — drop clusters by size and ratio thresholds.

        Args:
            subclusters: List of subclusters from Phase 9

        Returns:
            List of filtered clusters, sorted by size (largest first)
        """
        # Filter by absolute size
        large_clusters = [c for c in subclusters if len(c['read_ids']) >= self.min_size]
        small_clusters = [c for c in subclusters if len(c['read_ids']) < self.min_size]

        if small_clusters:
            filtered_count = len(small_clusters)
            # Track discarded reads from size-filtered clusters
            for cluster in small_clusters:
                self.discarded_read_ids.update(cluster['read_ids'])
            self._log_stage(
                f"Min-size filter: dropped {filtered_count} cluster(s) below {self.min_size}",
                large_clusters,
            )

        # Filter by relative size ratio
        if large_clusters and self.min_cluster_ratio > 0:
            largest_size = max(len(c['read_ids']) for c in large_clusters)
            passing_ratio = [c for c in large_clusters
                            if len(c['read_ids']) / largest_size >= self.min_cluster_ratio]
            failing_ratio = [c for c in large_clusters
                            if len(c['read_ids']) / largest_size < self.min_cluster_ratio]

            if failing_ratio:
                filtered_count = len(failing_ratio)
                # Track discarded reads from ratio-filtered clusters
                for cluster in failing_ratio:
                    self.discarded_read_ids.update(cluster['read_ids'])
                self._log_stage(
                    f"Min-ratio filter: dropped {filtered_count} cluster(s) below {self.min_cluster_ratio}",
                    passing_ratio,
                )

            large_clusters = passing_ratio

        # Sort by size and renumber as c1, c2, c3...
        large_clusters.sort(key=lambda c: len(c['read_ids']), reverse=True)

        return large_clusters

    def _write_cluster_outputs(self, clusters: List[Dict], output_file: str) -> Tuple[int, int]:
        """Phase 11: Output generation — generate final consensus and write output files.

        Args:
            clusters: List of filtered clusters from Phase 10
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

        # Compute identity group/variant ranks for the surviving clusters.
        # Groups are ranked by anchor (largest-member) size desc; variants
        # within each group are ranked by size desc. These drive gid=/vid=
        # header emission and the group_rank/variant_rank metadata fields,
        # giving summarize a stable end-to-end identity for each cluster.
        self._assign_identity_ranks(clusters)

        # Build work packages for each cluster. cluster_id is the core-assigned
        # group.variant designator (e.g., "1.v2"), matching the filename suffix
        # and summarize's output naming.
        work_packages = []
        cluster_dicts_by_idx = {}
        for final_idx, cluster_dict in enumerate(clusters, 1):
            cluster = cluster_dict['read_ids']
            group_rank = cluster_dict.get('identity_group_rank')
            variant_rank = cluster_dict.get('identity_variant_rank')
            if group_rank is not None and variant_rank is not None:
                cluster_id = f"{group_rank}.v{variant_rank}"
            else:
                # Defensive fallback for pipelines that bypass CER validation.
                cluster_id = f"c{final_idx}"
            cluster_dict['cluster_id'] = cluster_id
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

        # Aggregate counters for a single Outlier filter summary line at INFO
        outlier_dropped_clusters = 0
        outlier_trimmed_clusters = 0
        outlier_reads_moved = 0

        # Write output files sequentially (I/O bound, must preserve order)
        with open(output_file, 'w') as consensus_fasta_handle:
            for result in results:
                final_idx = result['final_idx']
                cluster = result['cluster']
                actual_size = result['actual_size']

                # Route any MAD-flagged outlier reads (and the full sampled set
                # if the cluster fell below min-reads after removal) to the
                # specimen-level discard pool. These reads will appear in the
                # {specimen}-discards.fastq and are not reassigned (Phase 6
                # reassignment has already completed).
                mad_outlier_ids = result.get('mad_outlier_ids') or set()
                if mad_outlier_ids:
                    self.discarded_read_ids.update(mad_outlier_ids)
                    outlier_reads_moved += len(mad_outlier_ids)
                    if result.get('dropped_by_min_reads'):
                        outlier_dropped_clusters += 1
                        logging.debug(
                            f"Cluster {final_idx}: dropped after MAD outlier "
                            f"removal left <3 reads ({len(mad_outlier_ids)} "
                            f"reads moved to discards)"
                        )
                    else:
                        outlier_trimmed_clusters += 1
                        logging.debug(
                            f"Cluster {final_idx}: removed "
                            f"{len(mad_outlier_ids)} outlier read(s) before "
                            f"final consensus"
                        )

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

                    # Compute err_factor from MSA + q_ctx (cluster-wide homogeneity).
                    err_factor = None
                    msa_string = result.get('msa')
                    if msa_string:
                        ef, obs_sum, exp_sum, ef_cols = compute_cluster_err_factor(
                            msa_string, self.qctx_table
                        )
                        err_factor = ef
                        cluster_dict['err_factor'] = ef
                        cluster_dict['err_factor_details'] = {
                            'obs_sum': obs_sum,
                            'exp_sum': exp_sum,
                            'cols': ef_cols,
                        }

                    # Write output files
                    self.write_cluster_files(
                        cluster_name=cluster_dict['cluster_id'],
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
                        err_factor=err_factor,
                        identity_group_rank=cluster_dict.get('identity_group_rank'),
                        identity_variant_rank=cluster_dict.get('identity_variant_rank'),
                    )

        # Build the surviving cluster list for end-of-pipeline summary lines.
        # A cluster survived if MAD didn't drop it; its post-outlier read count
        # is the assigned read_ids minus any MAD-flagged outliers.
        surviving_clusters = []
        surviving_reads = 0
        for result in results:
            if result.get('dropped_by_min_reads'):
                continue
            cluster_dict = cluster_dicts_by_idx.get(result['final_idx'], {})
            mad_outlier_ids = result.get('mad_outlier_ids') or set()
            surviving_clusters.append(cluster_dict)
            surviving_reads += len(cluster_dict.get('read_ids') or ()) - len(mad_outlier_ids)

        if outlier_dropped_clusters > 0 or outlier_trimmed_clusters > 0:
            parts = []
            if outlier_dropped_clusters:
                parts.append(f"dropped {outlier_dropped_clusters} cluster(s)")
            if outlier_trimmed_clusters:
                parts.append(f"trimmed {outlier_trimmed_clusters} cluster(s)")
            detail = ", ".join(parts)
            self._log_stage(
                f"Outlier filter: {detail} (-{outlier_reads_moved} reads)",
                surviving_clusters,
                reads_in_clusters=surviving_reads,
            )
        else:
            logging.debug("Outlier filter: no outliers removed")

        self._log_stage("Final", surviving_clusters, reads_in_clusters=surviving_reads)

        return clusters_with_ambiguities, total_ambiguity_positions

    def _write_discarded_reads(self) -> None:
        """Phase 12: Write discarded reads to a FASTQ file for inspection.

        Discards include:
        - Outlier reads removed during variant phasing
        - Reads from clusters dropped by the hard-floor filter or noise filter
        - Reads from clusters filtered out by size/ratio thresholds (Phase 10)
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
             3. Variant detection + phasing (split clusters by haplotype)
             4. Post-phasing merge (combine HP-equivalent subclusters)
             5. Noise filter (drop small clusters with no-majority columns)
             6. Read reassignment (concordance-based, within identity groups; optional)
             7. Discard recovery (re-admit dropped reads to clusters; optional, coupled to 6)
             8. Second phasing pass (split variants introduced by 6/7; coupled to phasing+6)
             9. CER validation (pairwise significance testing within identity groups)
            10. Size filtering (min size and cluster-ratio thresholds)
            11. Output generation (final consensus, FASTA writing)
            12. Write discarded reads (optional)

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
                small_reads = sum(len(c) for c in small)
                for c in small:
                    self.discarded_read_ids.update(c)
                initial_clusters = [c for c in initial_clusters if len(c) >= min_consensus_reads]
                self._log_stage(
                    f"Hard-floor filter: dropped {len(small)} cluster(s) below {min_consensus_reads} reads (-{small_reads} reads)",
                    initial_clusters,
                )

            # Phase 2: Pre-phasing merge
            merged_clusters = self._run_prephasing_merge(initial_clusters)

            # Phase 3: Variant detection + phasing
            all_subclusters = self._run_variant_phasing(merged_clusters)

            # Phase 4: Post-phasing merge
            merged_subclusters = self._run_postphasing_merge(all_subclusters)

            # Phase 5: Noise filter (drop small clusters with no-majority columns)
            cleaned_subclusters = self._filter_noisy_clusters(merged_subclusters)

            # Phase 6: Read reassignment (optional)
            if self.enable_read_reassignment:
                reassigned_subclusters = self._run_read_reassignment(cleaned_subclusters)
            else:
                reassigned_subclusters = cleaned_subclusters

            # Phase 7: Discard recovery (optional, requires Phase 6 enabled).
            # Discard recovery uses the same identity-group concordance logic
            # as read reassignment; running it without Phase 6 would silently
            # inject previously-discarded reads into clusters whose membership
            # the user asked to keep frozen. So Phase 7 is gated on Phase 6.
            if self.enable_discard_recovery and self.enable_read_reassignment:
                discard_reassigned = self._run_discard_reassignment(reassigned_subclusters)
            else:
                if self.enable_discard_recovery and not self.enable_read_reassignment:
                    logging.info("Discard recovery skipped: requires --enable-read-reassignment")
                discard_reassigned = reassigned_subclusters

            # Phase 8: Second phasing pass. Gated on the phasing flag AND on
            # Phase 6 having run — Phase 7 is coupled to Phase 6, so if 6 did
            # not run nothing has changed since Phase 3 and there is no new
            # work for Phase 8 to do.
            if self.enable_secondpass_phasing and self.enable_read_reassignment:
                rephased_subclusters = self._run_second_phasing_pass(discard_reassigned)
            else:
                rephased_subclusters = discard_reassigned

            # Phase 9: CER validation
            validated_subclusters = self._run_cer_validation(rephased_subclusters)

            # Phase 10: Size filtering
            filtered_clusters = self._run_size_filtering(validated_subclusters)

            # Capture for metadata serialization (write_metadata reads this).
            self.final_cluster_dicts = filtered_clusters

            # Phase 11: Output generation
            consensus_output_file = os.path.join(self.output_dir, f"{self.sample_name}-all.fasta")
            clusters_with_ambiguities, total_ambiguity_positions = self._write_cluster_outputs(
                filtered_clusters, consensus_output_file
            )

            # Phase 12: Write discarded reads (optional)
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
        """Raw edlib identity for initial read-to-read clustering.

        No homopolymer normalization: reads that differ only at HP length
        contribute to distance here and are pooled only by virtue of the
        loose --min-identity threshold (default 0.90). HP-aware comparisons
        live elsewhere: _hp_normalized_pairwise_compare for identity-group
        formation and post-clustering consensus merging, and the adjusted-
        identity path in speconsense.distances for downstream summarize work.
        """
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
            enable_secondpass_phasing=self.enable_secondpass_phasing,
            disable_homopolymer_equivalence=self.disable_homopolymer_equivalence,
            min_variant_frequency=self.min_variant_frequency,
            min_variant_count=self.min_variant_count,
            significance_level=self.significance_level,
            min_hp_length=self.min_hp_length,
            max_sample_size=self.max_sample_size,
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
        Only HP runs of length >= self.min_hp_length are normalized.
        Any non-HP difference (substitution, non-HP indel, terminal overhang)
        means the sequences are NOT equivalent.
        """
        if not seq1 or not seq2:
            return seq1 == seq2

        _, is_equiv = _hp_normalized_pairwise_compare(seq1, seq2, self.min_hp_length)
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
