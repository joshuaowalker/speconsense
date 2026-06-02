"""Command-line interface for speconsense-summarize.

Provides argument parsing, logging setup, and main entry point.
"""

import glob
import json
import os
import shutil
import sys
import argparse
import logging
import tempfile
from typing import List, Tuple, Dict, Optional
from collections import defaultdict

# Python 3.8 compatibility: BooleanOptionalAction was added in Python 3.9
if not hasattr(argparse, 'BooleanOptionalAction'):
    class BooleanOptionalAction(argparse.Action):
        def __init__(self, option_strings, dest, default=None, required=False, help=None):
            _option_strings = []
            for option_string in option_strings:
                _option_strings.append(option_string)
                if option_string.startswith('--'):
                    _option_strings.append('--no-' + option_string[2:])
            super().__init__(option_strings=_option_strings, dest=dest, nargs=0,
                           default=default, required=required, help=help)

        def __call__(self, parser, namespace, values, option_string=None):
            if option_string.startswith('--no-'):
                setattr(namespace, self.dest, False)
            else:
                setattr(namespace, self.dest, True)
    argparse.BooleanOptionalAction = BooleanOptionalAction

from tqdm import tqdm

try:
    from speconsense import __version__
except ImportError:
    # Fallback for when running as a script directly (e.g., in tests)
    __version__ = "dev"

from speconsense.profiles import (
    Profile,
    ProfileError,
    print_profiles_list,
)
from speconsense._help import install_advanced_help, add_advanced_argument
from speconsense.scalability import ScalabilityConfig
from speconsense.types import ConsensusInfo, OverlapMergeInfo

from .fields import parse_fasta_fields

from .io import (
    load_consensus_sequences,
    load_existing_specimen_outputs,
    build_fastq_lookup_table,
    write_specimen_data_files,
    write_ns_variant_files,
    write_lq_variant_files,
    write_output_files,
    load_metadata_from_json,
    strip_cluster_suffix,
)
from .clustering import (
    group_by_core_identity,
    merge_groups_by_anchor_overlap,
    select_variants,
)
from .merging import merge_group_with_msa, _refresh_cluster_metrics, _refresh_cer_factor
from .analysis import MAX_MSA_MERGE_VARIANTS
from .tree import write_specimen_variant_tree


# Merge effort configuration
MERGE_EFFORT_PRESETS = {
    'fast': 8,
    'balanced': 10,
    'thorough': 12,
}


def parse_merge_effort(spec: str) -> int:
    """Parse merge effort specification into numeric value.

    Args:
        spec: Preset name (fast, balanced, thorough) or numeric 6-14

    Returns:
        Effort level as integer

    Raises:
        ValueError: If spec is invalid
    """
    spec = spec.strip().lower()
    if spec in MERGE_EFFORT_PRESETS:
        return MERGE_EFFORT_PRESETS[spec]
    try:
        value = int(spec)
        if 6 <= value <= 14:
            return value
        raise ValueError(f"Numeric merge-effort must be 6-14, got {value}")
    except ValueError as e:
        if "invalid literal" in str(e):
            raise ValueError(
                f"Unknown merge-effort: '{spec}'. "
                f"Use preset (fast, balanced, thorough) or numeric 6-14"
            )
        raise


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Speconsense output with advanced variant handling.")
    install_advanced_help(parser)

    # Input/Output group
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("--source", type=str, default="clusters",
                          help="Source directory containing Speconsense output (default: clusters)")
    io_group.add_argument("--summary-dir", type=str, default="__Summary__",
                          help="Output directory for summary files (default: __Summary__)")
    io_group.add_argument("--specimen", type=str, default=None,
                          help="Process only this specimen. Loads only <specimen>-all.fasta from --source.")
    io_group.add_argument("--aggregate-only", action="store_true",
                          help="Skip processing. Generate aggregate summary from existing per-specimen outputs.")
    io_group.add_argument("--fasta-fields", type=str, default="default",
                          help="FASTA header fields to output. Can be: "
                               "(1) a preset name (default, minimal, qc, full, id-only), "
                               "(2) comma-separated field names (size, ric, length, rawric, "
                               "snp, rid, rid_min, primers, group, variant), or "
                               "(3) a combination of presets and fields (e.g., minimal,qc or "
                               "minimal,rid). Duplicates removed, order preserved "
                               "left to right. Default: default")

    # Filtering group
    filtering_group = parser.add_argument_group("Filtering")
    filtering_group.add_argument("--min-ric", type=int, default=3,
                                 help="Minimum Reads in Consensus (RiC) threshold (default: 3)")
    filtering_group.add_argument("--min-len", type=int, default=0,
                                 help="Minimum sequence length in bp (default: 0 = disabled)")
    filtering_group.add_argument("--max-len", type=int, default=0,
                                 help="Maximum sequence length in bp (default: 0 = disabled)")
    filtering_group.add_argument("--min-cer-factor", type=float, default=1.0,
                                 help="Minimum per-position CER factor for a variant to be kept "
                                      "as a primary output. Variants with cer_factor below this "
                                      "are routed to __Summary__/variants/ as .ns records. "
                                      "Variants with cer_factor=None (anchors, clusters without a "
                                      "valid pairwise comparison) always pass. Set to 0 to disable "
                                      "CER filtering. (default: 1.0)")
    filtering_group.add_argument("--max-err-factor", type=float, default=1.5,
                                 help="Maximum cluster err_factor (observed/q_ctx-expected "
                                      "disagreement ratio). Clusters above this threshold are "
                                      "routed to __Summary__/variants/ as .lq records. Variants "
                                      "with err_factor=None (legacy output) always pass. Set to 0 "
                                      "to disable err_factor filtering. (default: 1.5)")
    # Grouping group
    grouping_group = parser.add_argument_group("Grouping")
    grouping_group.add_argument("--group-identity", "--variant-group-identity",
                                dest="group_identity", type=float, default=0.85,
                                help="Anchor-to-anchor identity threshold for cross-primer overlap "
                                     "merging between core-assigned groups. Matches core's "
                                     "--group-identity default. (default: 0.85)")

    # Merging group
    merging_group = parser.add_argument_group("Merging")
    merging_group.add_argument("--merge-snp", action=argparse.BooleanOptionalAction, default=True,
                               help="Enable SNP-based merging (default: True, use --no-merge-snp to disable)")
    merging_group.add_argument("--merge-indel-length", type=int, default=0,
                               help="Maximum length of individual indels allowed in merging (default: 0 = disabled)")
    merging_group.add_argument("--merge-position-count", type=int, default=2,
                               help="Maximum total SNP+indel positions allowed in merging (default: 2)")
    merging_group.add_argument("--merge-min-size-ratio", type=float, default=0.1,
                               help="Minimum size ratio (smaller/larger) for merging clusters (default: 0.1, 0 to disable)")
    merging_group.add_argument("--min-merge-overlap", type=int, default=200,
                               help="Minimum overlap in bp for merging sequences of different lengths (default: 200, 0 to disable)")
    merging_group.add_argument("--merge-effort", type=str, default="balanced", metavar="LEVEL",
                               help="Merging effort level: fast (8), balanced (10), thorough (12), "
                                    "or numeric 6-14. Higher values allow larger batch sizes for "
                                    "exhaustive subset search. Default: balanced")
    merging_group.add_argument("--hp-normalization-length", type=int, default=6,
                               help="Minimum homopolymer run length at/above which HP length "
                                    "differences are blanket-normalized (treated as noise). Runs "
                                    "shorter than this are surfaced as real edits in both distance "
                                    "calculations and MSA merging. Matches core's "
                                    "--hp-normalization-length default. (default: 6)")

    # Backward compatibility: support old --snp-merge-limit parameter
    parser.add_argument("--snp-merge-limit", type=int, dest="_snp_merge_limit_deprecated",
                        help=argparse.SUPPRESS)  # Hidden but functional

    # Selection group
    selection_group = parser.add_argument_group("Selection")
    selection_group.add_argument("--select-max-groups", "--max-groups",
                                 dest="select_max_groups", type=int, default=-1,
                                 help="Maximum number of groups to output per specimen (default: -1 = all groups)")
    selection_group.add_argument("--select-max-variants", "--max-variants",
                                 dest="select_max_variants", type=int, default=-1,
                                 help="Maximum total variants to output per group (default: -1 = no limit, 0 also means no limit)")
    selection_group.add_argument("--select-strategy", "--variant-selection",
                                 dest="select_strategy", choices=["size", "diversity"], default="size",
                                 help="Variant selection strategy: size or diversity (default: size)")
    selection_group.add_argument("--select-min-size-ratio", type=float, default=0,
                                 help="Minimum size ratio (variant/largest) to include in output "
                                      "(default: 0 = disabled, e.g. 0.2 for 20%% cutoff)")
    selection_group.add_argument("--enable-full-consensus", action="store_true", default=False,
                                 help="Per identity group with >=2 selected variants, emit an "
                                      "additional ``-{gid}-full`` consensus built from a "
                                      "size-weighted, quality-sorted sample of the pre-merge core "
                                      "variants' reads. Intended as a query substrate for BLAST "
                                      "against legacy unphased references. Uses local SPOA "
                                      "alignment when the group spans multiple primer sets.")
    selection_group.add_argument("--disable-full-consensus", action="store_false",
                                 dest="enable_full_consensus",
                                 help="Override --enable-full-consensus or profile setting "
                                      "(e.g. when using the compressed profile but the -full "
                                      "artifact isn't wanted for this run).")

    # Performance group
    perf_group = parser.add_argument_group("Performance")
    perf_group.add_argument("--scale-threshold", type=int, default=1001,
                            help="Sequence count threshold for scalable mode in HAC clustering (requires vsearch). "
                                 "Set to 0 to disable. Default: 1001")
    perf_group.add_argument("--threads", type=int, default=0, metavar="N",
                            help="Max threads for internal parallelism. "
                                 "0=auto-detect (default), N>0 for explicit count.")

    # Advanced group (hidden from --help; view with --help-advanced).
    # Disables for integral merging behavior — primarily retained for
    # ablation studies, not user-facing tuning.
    advanced_group = parser.add_argument_group(
        "Advanced (pre-tuned — rarely needed)",
        "Hidden from --help; view with --help-advanced. Disable flags for "
        "integral merging phases — primarily retained for ablation studies.",
    )
    add_advanced_argument(advanced_group, "--disable-merging", action="store_true",
                          help="Disable all variant merging (skip MSA-based merge evaluation entirely)")
    add_advanced_argument(advanced_group, "--enable-merging", action="store_false", dest="disable_merging",
                          help="Override --disable-merging or profile setting")
    add_advanced_argument(advanced_group, "--disable-homopolymer-equivalence", action="store_true",
                          help="Disable homopolymer equivalence in merging (treat AAA vs AAAA as different)")
    add_advanced_argument(advanced_group, "--enable-homopolymer-equivalence", action="store_false",
                          dest="disable_homopolymer_equivalence",
                          help="Override --disable-homopolymer-equivalence or profile setting")

    # Version and profile options (default group)
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level")
    parser.add_argument("--version", action="version",
                        version=f"speconsense-summarize {__version__}",
                        help="Show program's version number and exit")
    parser.add_argument("-p", "--profile", metavar="NAME",
                        help="Load parameter profile (use --list-profiles to see available)")
    parser.add_argument("--list-profiles", action="store_true",
                        help="List available profiles and exit")

    # Handle --list-profiles early (before requiring other args)
    if '--list-profiles' in sys.argv:
        print_profiles_list('speconsense-summarize')
        sys.exit(0)

    # First pass: get profile name if specified
    pre_args, _ = parser.parse_known_args()

    # Track which arguments were explicitly provided on CLI
    explicit_args = set()
    for arg in sys.argv[1:]:
        if arg.startswith('--') and '=' in arg:
            explicit_args.add(arg.split('=')[0][2:].replace('-', '_'))
        elif arg.startswith('--'):
            explicit_args.add(arg[2:].replace('-', '_'))

    # Load and apply profile if specified
    loaded_profile = None
    if pre_args.profile:
        try:
            loaded_profile = Profile.load(pre_args.profile)
        except ProfileError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

        # Apply profile values to parser defaults (explicit CLI args will override)
        for key, value in loaded_profile.speconsense_summarize.items():
            attr_name = key.replace('-', '_')
            if attr_name not in explicit_args:
                parser.set_defaults(**{attr_name: value})

    args = parser.parse_args()

    # Store loaded profile for logging later
    args._loaded_profile = loaded_profile

    # Handle backward compatibility for deprecated parameters
    if args._snp_merge_limit_deprecated is not None:
        if '--snp-merge-limit' in sys.argv:
            logging.warning("--snp-merge-limit is deprecated, use --merge-position-count instead")
        args.merge_position_count = args._snp_merge_limit_deprecated

    # Validate mutual exclusion of --specimen and --aggregate-only
    if args.specimen and args.aggregate_only:
        parser.error("--specimen and --aggregate-only are mutually exclusive")

    if '--variant-group-identity' in sys.argv:
        logging.warning("--variant-group-identity is deprecated, use --group-identity instead")

    if '--max-variants' in sys.argv:
        logging.warning("--max-variants is deprecated, use --select-max-variants instead")

    if '--max-groups' in sys.argv:
        logging.warning("--max-groups is deprecated, use --select-max-groups instead")

    if '--variant-selection' in sys.argv:
        logging.warning("--variant-selection is deprecated, use --select-strategy instead")

    return args


def setup_logging(log_level: str, log_file: str = None):
    """Setup logging configuration with optional file output."""
    # Clear any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    # Set up root logger
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, log_level))
    logger.addHandler(console_handler)

    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return log_file

    return None


def process_single_specimen(file_consensuses: List[ConsensusInfo],
                           args,
                           ns_for_specimen: List[ConsensusInfo] = None,
                           lq_for_specimen: List[ConsensusInfo] = None,
                           fastq_lookup: Dict[str, List[str]] = None,
                           qctx_table: Dict[str, float] = None,
                           cer_alpha: float = 1e-5,
                           full_min_ambiguity_frequency: float = 0.10,
                           full_max_sample_size: int = 100,
                           specimen_global_size_total: Optional[int] = None,
                           ) -> Tuple[List[ConsensusInfo], Dict[str, List[str]], Dict, int, List[OverlapMergeInfo], Dict[str, List], List[ConsensusInfo], List[ConsensusInfo]]:
    """
    Process a single specimen file: bucket by core gid, MSA-merge within each
    bucket, conflate cross-primer groups, select variants, and emit final
    names that honor core's gid/vid except where cross-primer conflation moved
    a record between groups. Returns final consensus list, merge traceability,
    naming info (keyed by final gid), limited_count, overlap merge info,
    full-consensus reads, and the *annotated* ns/lq lists (with
    ``group_size_total`` / ``global_size_total`` populated for the
    frequency fields).

    ``ns_for_specimen`` and ``lq_for_specimen`` are the specimen's filtered
    records (CER-routed and err_factor-routed). They are not renamed, but their
    core-assigned vids contribute to the collision set when summarize allocates
    fresh vids for cross-primer-moved variants, and their ``size`` contributes
    to the conflation-aware ``group_size_total`` denominator. Callers should
    use the returned annotated lists when writing ``.ns`` / ``.lq`` outputs.

    ``specimen_global_size_total`` is the per-specimen denominator for
    ``global_frequency=`` (``total_input_reads`` from the metadata JSON).
    None when metadata is missing.
    """
    if not file_consensuses:
        return [], {}, {}, 0, [], {}, ns_for_specimen or [], lq_for_specimen or []

    ns_for_specimen = ns_for_specimen or []
    lq_for_specimen = lq_for_specimen or []

    file_name = os.path.basename(file_consensuses[0].file_path)
    logging.info(f"Processing specimen from file: {file_name}")

    # Phase 1: honor core-assigned group_rank (gid=) to form variant groups.
    # Core runs complete-linkage identity grouping on the clusters directly and
    # emits the resulting rank as gid=/vid= header fields. Summarize no longer
    # re-runs HAC inside a specimen; cross-primer overlap conflation is applied
    # as a separate step below.
    variant_groups = group_by_core_identity(file_consensuses)

    # Phase 1b: Cross-primer overlap conflation between core groups.
    # Core places ITS vs ITS2 (and other cross-locus pairs) into distinct
    # identity groups because global identity fails on length mismatch.
    # Summarize conflates them here via anchor-to-anchor overlap-aware distance
    # so the subsequent MSA-based merge sees them as one group. Supports the
    # primer-pool use case (multiple primers amplifying overlapping loci in
    # one tissue sample).
    if args.min_merge_overlap > 0 and len(variant_groups) > 1:
        variant_groups = merge_groups_by_anchor_overlap(
            variant_groups, args.min_merge_overlap, args.group_identity,
            hp_normalization_length=args.hp_normalization_length)

    # Phase 1c: Compute conflation-aware bucket totals for the
    # ``group_frequency=`` field. Each record's denominator is the sum of
    # ``size`` across every record (passed + ns + lq) in its post-conflation
    # bucket, so cross-primer-conflated variants are measured against the
    # merged-group total. Annotate every record's ``group_size_total`` and
    # ``global_size_total`` here so the field formatters can read them
    # directly without recomputing.
    core_to_bucket: Dict[int, int] = {}
    for final_gid, members in variant_groups.items():
        for member in members:
            if member.group_rank is not None:
                core_to_bucket[member.group_rank] = final_gid

    bucket_totals: Dict[int, int] = defaultdict(int)
    for record in file_consensuses:
        if record.group_rank is not None:
            bucket = core_to_bucket.get(record.group_rank, record.group_rank)
            bucket_totals[bucket] += record.size
    for record in ns_for_specimen:
        if record.group_rank is not None:
            bucket = core_to_bucket.get(record.group_rank, record.group_rank)
            bucket_totals[bucket] += record.size
    for record in lq_for_specimen:
        if record.group_rank is not None:
            bucket = core_to_bucket.get(record.group_rank, record.group_rank)
            bucket_totals[bucket] += record.size

    def _annotate(record: ConsensusInfo) -> ConsensusInfo:
        if record.group_rank is None:
            group_total = None
        else:
            bucket = core_to_bucket.get(record.group_rank, record.group_rank)
            group_total = bucket_totals.get(bucket)
        return record._replace(
            group_size_total=group_total,
            global_size_total=specimen_global_size_total,
        )

    # Re-bind annotated copies of each list. ``variant_groups`` is rebuilt
    # so downstream Phases see the annotated records.
    file_consensuses = [_annotate(r) for r in file_consensuses]
    ns_for_specimen = [_annotate(r) for r in ns_for_specimen]
    lq_for_specimen = [_annotate(r) for r in lq_for_specimen]
    variant_groups = {
        gid: [_annotate(m) for m in members]
        for gid, members in variant_groups.items()
    }

    # Filter to max groups if specified
    if args.select_max_groups > 0 and len(variant_groups) > args.select_max_groups:
        # Sort groups by size of largest member
        sorted_for_filtering = sorted(
            variant_groups.items(),
            key=lambda x: max(m.size for m in x[1]),
            reverse=True
        )
        # Keep only top N groups
        variant_groups = dict(sorted_for_filtering[:args.select_max_groups])
        logging.info(f"Filtered to top {args.select_max_groups} groups by size (from {len(sorted_for_filtering)} total groups)")

    # Phase 2: MSA-based merging within each group
    merged_groups = {}
    all_merge_traceability = {}
    total_limited_count = 0
    all_overlap_merges = []

    if args.disable_merging:
        # Skip merging entirely - pass variants through unchanged
        logging.info("Merging disabled - skipping MSA-based merge evaluation")
        for group_id, group_members in variant_groups.items():
            merged_groups[group_id] = group_members
    else:
        for group_id, group_members in variant_groups.items():
            merged, traceability, limited_count, overlap_merges = merge_group_with_msa(group_members, args)
            merged_groups[group_id] = merged
            all_merge_traceability.update(traceability)
            total_limited_count += limited_count
            all_overlap_merges.extend(overlap_merges)

    # Phase 2b: Recompute per-cluster metrics on merged records. For each
    # record that absorbed ≥2 contributors, re-derive:
    #   - rid, rid_min, err_factor from a SPOA MSA over the union of
    #     contributor reads (cluster-local, always run when fastq_lookup is
    #     available).
    #   - cer_factor from the post-merge peer landscape via
    #     ``classify_pairwise_differences`` + ``compute_cer_factor``. Same-
    #     primer merges recompute against larger peers in the same bucket;
    #     cross-primer merges set cer_factor=None (CER noise model doesn't
    #     apply to merged-locus candidates).
    # Inherited values from ``_build_merged_consensus_info`` are otherwise
    # stale because they describe the largest contributor's pre-merge state.
    if fastq_lookup is not None:
        original_consensus_lookup = {c.sample_name: c for c in file_consensuses}
        for group_id, group_members in merged_groups.items():
            for idx, record in enumerate(group_members):
                contributors_names = all_merge_traceability.get(record.sample_name)
                if not contributors_names or len(contributors_names) <= 1:
                    continue
                contributors = [
                    original_consensus_lookup[n]
                    for n in contributors_names
                    if n in original_consensus_lookup
                ]
                if len(contributors) <= 1:
                    continue
                refreshed = _refresh_cluster_metrics(
                    record,
                    contributors,
                    fastq_lookup=fastq_lookup,
                    qctx_table=qctx_table,
                    hp_normalization_length=args.hp_normalization_length,
                    disable_homopolymer_equivalence=args.disable_homopolymer_equivalence,
                )
                group_members[idx] = refreshed

        # cer_factor recompute: needs the post-rid/err_factor refresh to have
        # settled, and uses each bucket's full membership for peer lookup.
        for group_id, group_members in merged_groups.items():
            for idx, record in enumerate(group_members):
                contributors_names = all_merge_traceability.get(record.sample_name)
                if not contributors_names or len(contributors_names) <= 1:
                    continue
                contributors = [
                    original_consensus_lookup[n]
                    for n in contributors_names
                    if n in original_consensus_lookup
                ]
                if len(contributors) <= 1:
                    continue
                refreshed = _refresh_cer_factor(
                    record,
                    contributors,
                    bucket_members=group_members,
                    qctx_table=qctx_table,
                    alpha=cer_alpha,
                    hp_min_length=args.hp_normalization_length,
                )
                group_members[idx] = refreshed

    # Phase 3: Select representative variants and emit final names.
    #
    # Naming policy: preserve core's gid/vid verbatim except for variants moved
    # between groups by cross-primer overlap conflation in Phase 1b. The dict
    # key in ``merged_groups`` is the survivor's pre-conflation core gid (the
    # absorbing bucket). For each variant in that bucket:
    #   - If ``variant.group_rank == final_gid``, the variant originated in
    #     this group; emit ``(final_gid, variant.variant_rank)`` directly,
    #     letting MSA-merge and filter drops manifest as vid gaps.
    #   - Else the variant was moved here from an absorbed group; allocate a
    #     new vid above ``max(used_vids)``.
    #
    # ``used_vids`` for a given final_gid is seeded from every record core
    # emitted under that gid (passed + ns + lq), pre-MSA-merge, so moved
    # members never collide with names core ever produced — including vids
    # that became gaps via within-group merging or filtering. Each newly
    # allocated vid is added to ``used_vids`` so subsequent moved members in
    # the same pass (e.g., 3+ group conflation) cannot collide with each other.
    final_consensus = []
    naming_info = {}
    full_consensus_reads: Dict[str, List] = {}

    # Pre-index core-emitted vids by core gid. Pre-merge passed records plus
    # ns/lq records — covers everything core wrote to disk under each gid.
    core_vids_by_gid: Dict[int, set] = defaultdict(set)
    for record in file_consensuses:
        if record.group_rank is not None and record.variant_rank is not None:
            core_vids_by_gid[record.group_rank].add(record.variant_rank)
    for record in ns_for_specimen:
        if record.group_rank is not None and record.variant_rank is not None:
            core_vids_by_gid[record.group_rank].add(record.variant_rank)
    for record in lq_for_specimen:
        if record.group_rank is not None and record.variant_rank is not None:
            core_vids_by_gid[record.group_rank].add(record.variant_rank)

    # Iterate groups in anchor-size desc order for stable selection logging;
    # naming itself is gid-based so iteration order does not affect output.
    sorted_groups = sorted(merged_groups.items(),
                          key=lambda x: max(m.size for m in x[1]),
                          reverse=True)

    for final_gid, group_members in sorted_groups:
        # Apply select-min-size-ratio filter
        if args.select_min_size_ratio > 0 and len(group_members) > 1:
            largest_size = max(v.size for v in group_members)
            filtered = [v for v in group_members
                        if (v.size / largest_size) >= args.select_min_size_ratio]
            if len(filtered) < len(group_members):
                filtered_count = len(group_members) - len(filtered)
                logging.debug(f"Group {final_gid}: filtered out {filtered_count} "
                              f"variants with size ratio < {args.select_min_size_ratio} "
                              f"relative to largest (size={largest_size})")
                group_members = filtered

        # Select variants for this group
        selected_variants = select_variants(
            group_members, args.select_max_variants, args.select_strategy,
            group_number=final_gid,
            hp_normalization_length=args.hp_normalization_length)

        used_vids = set(core_vids_by_gid.get(final_gid, set()))

        group_naming = []
        for variant in selected_variants:
            specimen_base = strip_cluster_suffix(variant.sample_name)
            if variant.group_rank == final_gid and variant.variant_rank is not None:
                new_vid = variant.variant_rank
            else:
                new_vid = (max(used_vids) + 1) if used_vids else 1
                used_vids.add(new_vid)
            new_name = f"{specimen_base}-{final_gid}.v{new_vid}"

            renamed_variant = variant._replace(sample_name=new_name)
            final_consensus.append(renamed_variant)
            group_naming.append((variant.sample_name, new_name))

        naming_info[final_gid] = group_naming

        # Phase 3b: Optional -full group consensus. Source from
        # variant_groups[final_gid] (pre-MSA-merge, post-conflation core
        # variants — .ns/.lq already excluded). Independent of
        # merged_groups / select_variants output, so MSA merging never
        # collapses the pool the -full pipeline draws from.
        if getattr(args, 'enable_full_consensus', False) and len(selected_variants) >= 2:
            from .full_consensus import build_group_full_consensus
            pre_merge_variants = variant_groups.get(final_gid, [])
            specimen_base_for_full = strip_cluster_suffix(selected_variants[0].sample_name)
            full_result = build_group_full_consensus(
                final_gid=final_gid,
                group_members=pre_merge_variants,
                pass_track_count=len(selected_variants),
                fastq_lookup=fastq_lookup,
                min_ambiguity_frequency=full_min_ambiguity_frequency,
                max_sample_size=full_max_sample_size,
                specimen_base=specimen_base_for_full,
                primary_file_path=selected_variants[0].file_path,
                group_size_total=bucket_totals.get(final_gid),
                global_size_total=specimen_global_size_total,
            )
            if full_result is not None:
                full_record, sampled_reads = full_result
                final_consensus.append(full_record)
                full_consensus_reads[full_record.sample_name] = sampled_reads

    logging.info(f"Processed {file_name}: {len(final_consensus)} final variants across {len(merged_groups)} groups")

    return (final_consensus, all_merge_traceability, naming_info,
            total_limited_count, all_overlap_merges, full_consensus_reads,
            ns_for_specimen, lq_for_specimen)


def _resolve_specimen_cer_context(
    source_folder: str,
    file_path: str,
    cache: Dict[str, Optional[Dict[str, float]]],
) -> Tuple[Optional[Dict[str, float]], float]:
    """Load the q_ctx table and CER alpha that core used for this specimen.

    Reads ``parameters.error_model`` and ``parameters.significance_level``
    from the specimen's metadata JSON in ``cluster_debug/`` and loads the
    corresponding q_ctx table via ``speconsense.qctx.load_table``. Tables
    are cached by model name across specimens. Returns ``(qctx_table,
    alpha)`` where ``qctx_table`` may be None (metadata missing, model
    unloadable) — callers must treat None as "skip err_factor / cer_factor
    recompute". ``alpha`` falls back to the canonical 1e-5 default when
    metadata doesn't specify one.
    """
    DEFAULT_ALPHA = 1e-5
    specimen_base = os.path.basename(file_path)
    if specimen_base.endswith('-all.fasta'):
        specimen_base = specimen_base[:-len('-all.fasta')]
    metadata = load_metadata_from_json(source_folder, specimen_base)
    if not metadata:
        return None, DEFAULT_ALPHA
    params = metadata.get('parameters') or {}
    alpha = float(params.get('significance_level') or DEFAULT_ALPHA)
    model_name = params.get('error_model')
    if not model_name:
        return None, alpha
    if model_name in cache:
        return cache[model_name], alpha
    try:
        from speconsense.qctx import load_table
        table = load_table(model_name)
    except Exception as e:
        logging.debug(
            f"Could not load q_ctx model '{model_name}' for merge recompute: {e}"
        )
        table = None
    cache[model_name] = table
    return table, alpha


def _resolve_specimen_global_total(
    source_folder: str,
    file_path: str,
) -> Optional[int]:
    """Return the per-specimen denominator for ``global_frequency=``.

    Reads ``total_input_reads`` from the specimen's metadata JSON — the
    post-presample-cap count actually fed into clustering. Returns ``None``
    when the metadata is missing or doesn't carry that field; the
    ``global_frequency=`` formatter then silently omits the field on
    that specimen's records.
    """
    specimen_base = os.path.basename(file_path)
    if specimen_base.endswith('-all.fasta'):
        specimen_base = specimen_base[:-len('-all.fasta')]
    metadata = load_metadata_from_json(source_folder, specimen_base)
    if not metadata:
        return None
    value = metadata.get('total_input_reads')
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _resolve_full_consensus_params(
    source_folder: str,
    file_path: str,
) -> Tuple[float, int]:
    """Return ``(min_ambiguity_frequency, max_sample_size)`` for the
    ``-full`` builder, sourced from the specimen's metadata JSON. Falls
    back to the speconsense defaults when metadata is missing — older
    runs predate per-field serialization.
    """
    DEFAULT_MIN_AMBIG = 0.10
    DEFAULT_MAX_SAMPLE = 100
    specimen_base = os.path.basename(file_path)
    if specimen_base.endswith('-all.fasta'):
        specimen_base = specimen_base[:-len('-all.fasta')]
    metadata = load_metadata_from_json(source_folder, specimen_base)
    if not metadata:
        return DEFAULT_MIN_AMBIG, DEFAULT_MAX_SAMPLE
    params = metadata.get('parameters') or {}
    min_ambig = params.get('min_ambiguity_frequency')
    max_sample = params.get('max_sample_size')
    return (
        float(min_ambig) if min_ambig is not None else DEFAULT_MIN_AMBIG,
        int(max_sample) if max_sample is not None else DEFAULT_MAX_SAMPLE,
    )


def _clean_specimen_output(summary_dir: str, specimen_id: str) -> None:
    """Remove previous output files for a specimen before reprocessing.

    Cleans files matching {specimen_id}* from summary_dir and its subdirs
    to prevent stale files from accumulating when variant counts change.
    """
    dirs_to_clean = [
        summary_dir,
        os.path.join(summary_dir, 'FASTQ Files'),
        os.path.join(summary_dir, 'variants'),
        os.path.join(summary_dir, 'variants', 'FASTQ Files'),
        os.path.join(summary_dir, 'trees'),
    ]
    removed = 0
    for d in dirs_to_clean:
        if not os.path.exists(d):
            continue
        for f in glob.glob(os.path.join(d, f"{glob.escape(specimen_id)}*")):
            try:
                os.unlink(f)
                removed += 1
            except OSError as e:
                logging.warning(f"Could not remove stale file {f}: {e}")
    if removed:
        logging.info(f"Cleaned {removed} previous output files for {specimen_id}")


def _cleanup_log(log_path: str) -> None:
    """Clean up temporary log file."""
    try:
        os.unlink(log_path)
    except Exception as e:
        logging.debug(f"Could not clean up temporary log file: {e}")


def main():
    """Main function to process command line arguments and run the summarization."""
    args = parse_arguments()

    # Parse FASTA field specification early
    try:
        fasta_fields = parse_fasta_fields(args.fasta_fields)
    except ValueError as e:
        logging.error(f"Invalid --fasta-fields specification: {e}")
        sys.exit(1)

    # Parse merge effort specification
    try:
        args.merge_effort_value = parse_merge_effort(args.merge_effort)
    except ValueError as e:
        logging.error(f"Invalid --merge-effort: {e}")
        sys.exit(1)

    # Set up logging with temporary log file
    temp_log_file = tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False)
    temp_log_file.close()

    setup_logging(args.log_level, temp_log_file.name)

    logging.info(f"speconsense-summarize version {__version__}")
    if args._loaded_profile:
        logging.info(f"Using profile '{args._loaded_profile.name}': {args._loaded_profile.description}")
    logging.info(f"Command: speconsense-summarize {' '.join(sys.argv[1:])}")
    logging.info("")
    logging.info("Starting enhanced speconsense summarization")
    logging.info(f"Parameters:")
    logging.info(f"  --source: {args.source}")
    logging.info(f"  --summary-dir: {args.summary_dir}")
    logging.info(f"  --min-ric: {args.min_ric}")
    logging.info(f"  --min-len: {args.min_len}")
    logging.info(f"  --max-len: {args.max_len}")
    logging.info(f"  --fasta-fields: {args.fasta_fields}")
    logging.info(f"  --merge-snp: {args.merge_snp}")
    logging.info(f"  --merge-indel-length: {args.merge_indel_length}")
    logging.info(f"  --merge-position-count: {args.merge_position_count}")
    logging.info(f"  --merge-min-size-ratio: {args.merge_min_size_ratio}")
    logging.info(f"  --disable-homopolymer-equivalence: {args.disable_homopolymer_equivalence}")
    logging.info(f"  --min-merge-overlap: {args.min_merge_overlap}")
    logging.info(f"  --merge-effort: {args.merge_effort} ({args.merge_effort_value})")
    logging.info(f"  --hp-normalization-length: {args.hp_normalization_length}")
    logging.info(f"  --group-identity: {args.group_identity}")
    logging.info(f"  --select-max-variants: {args.select_max_variants}")
    logging.info(f"  --select-max-groups: {args.select_max_groups}")
    logging.info(f"  --select-strategy: {args.select_strategy}")
    logging.info(f"  --select-min-size-ratio: {args.select_min_size_ratio}")
    logging.info(f"  --log-level: {args.log_level}")
    logging.info("")

    # --- Aggregate-only mode ---
    if args.aggregate_only:
        logging.info("Aggregate-only mode: generating summary from existing per-specimen outputs")
        all_final_consensus = load_existing_specimen_outputs(args.summary_dir)
        if not all_final_consensus:
            logging.error("No existing specimen outputs found in summary directory")
            _cleanup_log(temp_log_file.name)
            return

        write_output_files(
            all_final_consensus,
            [],  # no raw consensuses available in aggregate mode
            args.summary_dir,
            temp_log_file.name,
            fasta_fields
        )

        logging.info(f"Aggregate summary complete: {len(all_final_consensus)} sequences")
        _cleanup_log(temp_log_file.name)
        return

    logging.info("Processing each specimen file independently to organize variants within specimens")

    # Load consensus sequences (optionally filtered to single specimen)
    consensus_list, ns_list, lq_list = load_consensus_sequences(
        args.source, args.min_ric, args.min_len, args.max_len,
        specimen_id=args.specimen,
        min_cer_factor=args.min_cer_factor,
        max_err_factor=args.max_err_factor,
    )
    if not consensus_list and not ns_list and not lq_list:
        logging.error("No consensus sequences found")
        _cleanup_log(temp_log_file.name)
        return

    # Group consensus sequences by input file (one file per specimen)
    file_groups = defaultdict(list)
    for cons in consensus_list:
        file_groups[cons.file_path].append(cons)
    ns_by_file = defaultdict(list)
    for cons in ns_list:
        ns_by_file[cons.file_path].append(cons)
    lq_by_file = defaultdict(list)
    for cons in lq_list:
        lq_by_file[cons.file_path].append(cons)

    # Create output directories before processing
    os.makedirs(args.summary_dir, exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'FASTQ Files'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'variants'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'variants', 'FASTQ Files'), exist_ok=True)

    # In single-specimen mode, clean previous output for this specimen
    if args.specimen:
        _clean_specimen_output(args.summary_dir, args.specimen)

    # Build lookup tables once before processing loop
    fastq_lookup = build_fastq_lookup_table(args.source)
    original_consensus_lookup = {cons.sample_name: cons for cons in consensus_list}

    # Cache of qctx tables keyed by error_model name (read per-specimen from
    # the metadata JSON). Used by the merge-time recompute of err_factor and
    # cer_factor. Populated lazily inside the per-specimen loop.
    qctx_table_cache: Dict[str, Optional[Dict[str, float]]] = {}

    # Process each specimen file independently
    all_final_consensus = []
    all_merge_traceability = {}
    all_naming_info = {}
    all_raw_consensuses = []  # Collect .raw files from all specimens
    all_overlap_merges = []  # Collect overlap merge info for quality reporting
    total_limited_merges = 0

    # Iterate the union of file paths from all three lists so specimens whose
    # variants are entirely routed to .ns/.lq still get their filtered-variant
    # files written. process_single_specimen and write_specimen_data_files
    # both no-op cleanly on an empty consensus list.
    all_file_paths = (
        set(file_groups.keys()) | set(ns_by_file.keys()) | set(lq_by_file.keys())
    )
    sorted_file_paths = sorted(all_file_paths)
    for file_path in tqdm(sorted_file_paths, desc="Processing specimens", unit="specimen"):
        file_consensuses = file_groups.get(file_path, [])
        specimen_ns = ns_by_file.get(file_path, [])
        specimen_lq = lq_by_file.get(file_path, [])

        # Resolve the q_ctx table + CER alpha for this specimen so merge-time
        # recompute of err_factor and cer_factor uses the same model and
        # significance level core did. Falls back silently if the metadata
        # JSON is missing or the model is unloadable — merge recompute then
        # skips err_factor/cer_factor.
        specimen_qctx, specimen_alpha = _resolve_specimen_cer_context(
            args.source, file_path, qctx_table_cache
        )

        # When --enable-full-consensus is on, also resolve the per-specimen
        # IUPAC frequency and read budget from the core run's metadata so
        # the -full builder reuses the exact thresholds the per-cluster
        # consensus already applied.
        if args.enable_full_consensus:
            full_min_ambig, full_max_sample = _resolve_full_consensus_params(args.source, file_path)
        else:
            full_min_ambig, full_max_sample = 0.10, 100

        # Per-specimen denominator for global_frequency=. Sourced from
        # parameters.total_input_reads in the metadata JSON (post-presample
        # count). None when metadata is missing → global_frequency= is
        # silently omitted from headers for that specimen.
        specimen_global_total = _resolve_specimen_global_total(args.source, file_path)

        # Process specimen. ns/lq slices feed the vid collision-avoidance
        # allocator for cross-primer-conflated records. fastq_lookup +
        # qctx_table feed the merge-time metric recompute (rid, rid_min,
        # err_factor, cer_factor) on records that absorbed ≥2 contributors.
        # full_* params drive the optional -{gid}-full group consensus.
        # specimen_global_size_total feeds the global_frequency= field.
        # ns/lq lists come back annotated with the new frequency
        # denominators so the .ns/.lq writers emit those fields too.
        (final_consensus, merge_traceability, naming_info,
         limited_count, overlap_merges, full_reads_by_name,
         specimen_ns, specimen_lq) = process_single_specimen(
            file_consensuses, args,
            ns_for_specimen=specimen_ns,
            lq_for_specimen=specimen_lq,
            fastq_lookup=fastq_lookup,
            qctx_table=specimen_qctx,
            cer_alpha=specimen_alpha,
            full_min_ambiguity_frequency=full_min_ambig,
            full_max_sample_size=full_max_sample,
            specimen_global_size_total=specimen_global_total,
        )

        # Write individual data files immediately
        specimen_raw_consensuses = write_specimen_data_files(
            final_consensus,
            merge_traceability,
            naming_info,
            args.summary_dir,
            os.path.join(args.summary_dir, 'FASTQ Files'),
            fastq_lookup,
            original_consensus_lookup,
            fasta_fields,
            full_reads_by_name=full_reads_by_name,
        )

        # Emit ns variants (CER-filtered) for this specimen
        if specimen_ns:
            write_ns_variant_files(
                specimen_ns, args.summary_dir, fastq_lookup, fasta_fields
            )

        # Emit lq variants (err_factor-filtered) for this specimen
        if specimen_lq:
            write_lq_variant_files(
                specimen_lq, args.summary_dir, fastq_lookup, fasta_fields
            )

        # Emit per-specimen variant tree (passed + ns + lq, grouped by core gid)
        specimen_id = os.path.basename(file_path)
        if specimen_id.endswith('-all.fasta'):
            specimen_id = specimen_id[:-len('-all.fasta')]
        write_specimen_variant_tree(
            specimen_id=specimen_id,
            passed=final_consensus,
            ns=specimen_ns,
            lq=specimen_lq,
            output_dir=os.path.join(args.summary_dir, 'trees'),
            hp_normalization_length=args.hp_normalization_length,
        )

        # Accumulate results for summary files
        all_final_consensus.extend(final_consensus)
        all_merge_traceability.update(merge_traceability)
        all_raw_consensuses.extend(specimen_raw_consensuses)
        all_overlap_merges.extend(overlap_merges)
        total_limited_merges += limited_count

        # Update naming info with unique keys per specimen
        file_name = os.path.basename(file_path)
        for group_id, group_naming in naming_info.items():
            unique_key = f"{file_name}_{group_id}"
            all_naming_info[unique_key] = group_naming

    # --- Single-specimen mode: print JSON to stdout, skip aggregate files ---
    if args.specimen:
        variants = []
        for cons in all_final_consensus:
            variants.append({
                "name": cons.sample_name,
                "ric": cons.ric,
                "size": cons.size,
                "length": len(cons.sequence),
            })
        result = {
            "specimen_id": args.specimen,
            "variant_count": len(variants),
            "variants": variants,
        }
        print(json.dumps(result))
        logging.info(f"Single-specimen mode complete: {len(variants)} variants for {args.specimen}")
        _cleanup_log(temp_log_file.name)
        return

    # Write summary files at end (after all processing)
    write_output_files(
        all_final_consensus,
        all_raw_consensuses,
        args.summary_dir,
        temp_log_file.name,
        fasta_fields
    )

    # Write quality report (deferred import to avoid circular dependency).
    # Pass pre-merge consensus_list (with full per-cluster err/cer metrics)
    # and the routed ns/lq lists so the report can surface filter decisions.
    from speconsense import quality_report
    quality_report.write_quality_report(
        all_final_consensus,
        all_raw_consensuses,
        args.summary_dir,
        args.source,
        all_overlap_merges,
        args.min_merge_overlap,
        consensus_list=consensus_list,
        ns_list=ns_list,
        lq_list=lq_list,
        min_cer_factor=args.min_cer_factor,
        max_err_factor=args.max_err_factor,
    )

    logging.info(f"Enhanced summarization completed successfully")
    logging.info(f"Final output: {len(all_final_consensus)} consensus sequences in {args.summary_dir}")

    # Report if any variant groups were potentially suboptimal due to size
    if total_limited_merges > 0:
        logging.info(f"Note: {total_limited_merges} variant group(s) had >{MAX_MSA_MERGE_VARIANTS} variants (results potentially suboptimal)")

    _cleanup_log(temp_log_file.name)


if __name__ == "__main__":
    main()
