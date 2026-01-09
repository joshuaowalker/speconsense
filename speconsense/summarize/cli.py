"""Command-line interface for speconsense-summarize.

Provides argument parsing, logging setup, and main entry point.
"""

import os
import sys
import argparse
import logging
import tempfile
from typing import List, Tuple, Dict
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
from speconsense.scalability import ScalabilityConfig
from speconsense.types import ConsensusInfo, OverlapMergeInfo

from .fields import parse_fasta_fields
from .io import (
    load_consensus_sequences,
    build_fastq_lookup_table,
    write_specimen_data_files,
    write_output_files,
)
from .clustering import perform_hac_clustering, select_variants
from .merging import merge_group_with_msa
from .analysis import MAX_MSA_MERGE_VARIANTS, MIN_MERGE_BATCH, MAX_MERGE_BATCH


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

    # Input/Output group
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("--source", type=str, default="clusters",
                          help="Source directory containing Speconsense output (default: clusters)")
    io_group.add_argument("--summary-dir", type=str, default="__Summary__",
                          help="Output directory for summary files (default: __Summary__)")
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

    # Grouping group
    grouping_group = parser.add_argument_group("Grouping")
    grouping_group.add_argument("--group-identity", "--variant-group-identity",
                                dest="group_identity", type=float, default=0.9,
                                help="Identity threshold for variant grouping using HAC (default: 0.9)")

    # Merging group
    merging_group = parser.add_argument_group("Merging")
    merging_group.add_argument("--disable-merging", action="store_true",
                               help="Disable all variant merging (skip MSA-based merge evaluation entirely)")
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
    merging_group.add_argument("--disable-homopolymer-equivalence", action="store_true",
                               help="Disable homopolymer equivalence in merging (treat AAA vs AAAA as different)")
    merging_group.add_argument("--merge-effort", type=str, default="balanced", metavar="LEVEL",
                               help="Merging effort level: fast (8), balanced (10), thorough (12), "
                                    "or numeric 6-14. Higher values allow larger batch sizes for "
                                    "exhaustive subset search. Default: balanced")

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

    # Performance group
    perf_group = parser.add_argument_group("Performance")
    perf_group.add_argument("--scale-threshold", type=int, default=1001,
                            help="Sequence count threshold for scalable mode in HAC clustering (requires vsearch). "
                                 "Set to 0 to disable. Default: 1001")
    perf_group.add_argument("--threads", type=int, default=0, metavar="N",
                            help="Max threads for internal parallelism. "
                                 "0=auto-detect (default), N>0 for explicit count.")

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
                           args) -> Tuple[List[ConsensusInfo], Dict[str, List[str]], Dict, int, List[OverlapMergeInfo]]:
    """
    Process a single specimen file: HAC cluster, MSA-based merge per group, and select final variants.
    Returns final consensus list, merge traceability, naming info, limited_count, and overlap merge info.

    Architecture (Phase 3):
    1. HAC clustering to separate variant groups (primary vs contaminants)
    2. MSA-based merging within each group
    3. Select representative variants per group
    """
    if not file_consensuses:
        return [], {}, {}, 0, []

    file_name = os.path.basename(file_consensuses[0].file_path)
    logging.info(f"Processing specimen from file: {file_name}")

    # Phase 1: HAC clustering to separate variant groups (moved before merging!)
    scale_threshold = getattr(args, 'scale_threshold', 1001)
    threads_arg = getattr(args, 'threads', 0)
    max_threads = threads_arg if threads_arg > 0 else os.cpu_count()
    scalability_config = None
    if scale_threshold > 0:
        scalability_config = ScalabilityConfig(
            enabled=True,
            activation_threshold=scale_threshold,
            max_threads=max_threads
        )

    variant_groups = perform_hac_clustering(
        file_consensuses, args.group_identity, min_overlap_bp=args.min_merge_overlap,
        scalability_config=scalability_config, output_dir=getattr(args, 'source', '.')
    )

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

    # Phase 3: Select representative variants for each group in this specimen
    final_consensus = []
    naming_info = {}

    # Sort variant groups by size of largest member (descending)
    sorted_groups = sorted(merged_groups.items(),
                          key=lambda x: max(m.size for m in x[1]),
                          reverse=True)

    for group_idx, (_, group_members) in enumerate(sorted_groups):
        final_group_name = group_idx + 1

        # Select variants for this group
        selected_variants = select_variants(group_members, args.select_max_variants, args.select_strategy, group_number=final_group_name)

        # Create naming for this group within this specimen
        group_naming = []

        for variant_idx, variant in enumerate(selected_variants):
            # All variants get .v suffix (primary is .v1, additional are .v2, .v3, etc.)
            # Use rsplit to split on the LAST '-c' (specimen names may contain '-c')
            specimen_base = variant.sample_name.rsplit('-c', 1)[0]
            new_name = f"{specimen_base}-{group_idx + 1}.v{variant_idx + 1}"

            # Use _replace to preserve all fields while updating sample_name
            renamed_variant = variant._replace(sample_name=new_name)

            final_consensus.append(renamed_variant)
            group_naming.append((variant.sample_name, new_name))

        naming_info[group_idx + 1] = group_naming

    logging.info(f"Processed {file_name}: {len(final_consensus)} final variants across {len(merged_groups)} groups")

    return final_consensus, all_merge_traceability, naming_info, total_limited_count, all_overlap_merges


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
    logging.info(f"  --group-identity: {args.group_identity}")
    logging.info(f"  --select-max-variants: {args.select_max_variants}")
    logging.info(f"  --select-max-groups: {args.select_max_groups}")
    logging.info(f"  --select-strategy: {args.select_strategy}")
    logging.info(f"  --log-level: {args.log_level}")
    logging.info("")
    logging.info("Processing each specimen file independently to organize variants within specimens")

    # Load all consensus sequences
    consensus_list = load_consensus_sequences(
        args.source, args.min_ric, args.min_len, args.max_len
    )
    if not consensus_list:
        logging.error("No consensus sequences found")
        return

    # Group consensus sequences by input file (one file per specimen)
    file_groups = defaultdict(list)
    for cons in consensus_list:
        file_groups[cons.file_path].append(cons)

    # Create output directories before processing
    os.makedirs(args.summary_dir, exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'FASTQ Files'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'variants'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'variants', 'FASTQ Files'), exist_ok=True)

    # Build lookup tables once before processing loop
    fastq_lookup = build_fastq_lookup_table(args.source)
    original_consensus_lookup = {cons.sample_name: cons for cons in consensus_list}

    # Process each specimen file independently
    all_final_consensus = []
    all_merge_traceability = {}
    all_naming_info = {}
    all_raw_consensuses = []  # Collect .raw files from all specimens
    all_overlap_merges = []  # Collect overlap merge info for quality reporting
    total_limited_merges = 0

    sorted_file_paths = sorted(file_groups.keys())
    for file_path in tqdm(sorted_file_paths, desc="Processing specimens", unit="specimen"):
        file_consensuses = file_groups[file_path]

        # Process specimen
        final_consensus, merge_traceability, naming_info, limited_count, overlap_merges = process_single_specimen(
            file_consensuses, args
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
            fasta_fields
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

    # Write summary files at end (after all processing)
    write_output_files(
        all_final_consensus,
        all_raw_consensuses,
        args.summary_dir,
        temp_log_file.name,
        fasta_fields
    )

    # Write quality report (deferred import to avoid circular dependency)
    from speconsense import quality_report
    quality_report.write_quality_report(
        all_final_consensus,
        all_raw_consensuses,
        args.summary_dir,
        args.source,
        all_overlap_merges,
        args.min_merge_overlap
    )

    logging.info(f"Enhanced summarization completed successfully")
    logging.info(f"Final output: {len(all_final_consensus)} consensus sequences in {args.summary_dir}")

    # Report if any variant groups were potentially suboptimal due to size
    if total_limited_merges > 0:
        logging.info(f"Note: {total_limited_merges} variant group(s) had >{MAX_MSA_MERGE_VARIANTS} variants (results potentially suboptimal)")

    # Clean up temporary log file
    try:
        os.unlink(temp_log_file.name)
    except Exception as e:
        logging.debug(f"Could not clean up temporary log file: {e}")


if __name__ == "__main__":
    main()
