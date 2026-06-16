"""CLI and entry point for speconsense core clustering tool."""

import argparse
import logging
import os
import sys

from Bio import SeqIO

try:
    from speconsense import __version__
except ImportError:
    __version__ = "dev"

from speconsense.profiles import (
    Profile,
    ProfileError,
    print_profiles_list,
    show_profile,
)
from speconsense._help import install_advanced_help, add_advanced_argument

from .clusterer import SpecimenClusterer


def main():
    parser = argparse.ArgumentParser(
        description="MCL-based clustering of nanopore amplicon reads"
    )
    install_advanced_help(parser)

    # Input/Output group
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("input_file", help="Input FASTQ file")
    io_group.add_argument("-O", "--output-dir", default="clusters",
                          help="Output directory for all files (default: clusters)")
    io_group.add_argument("--primers", help="FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)")

    # Clustering group
    clustering_group = parser.add_argument_group("Clustering")
    clustering_group.add_argument("--algorithm", type=str, default="graph", choices=["graph", "greedy"],
                                  help="Clustering algorithm to use (default: graph)")
    clustering_group.add_argument("--min-identity", type=float, default=0.9,
                                  help="Minimum sequence identity threshold for clustering (default: 0.9)")

    # Filtering group
    filtering_group = parser.add_argument_group("Filtering")
    filtering_group.add_argument("--min-size", type=int, default=3,
                                 help="Minimum cluster size (default: 3, 0 to disable)")
    filtering_group.add_argument("--max-sample-size", type=int, default=100,
                                 help="Maximum cluster size for consensus (default: 100)")

    # Variant Phasing & Validation group
    phasing_group = parser.add_argument_group("Variant Phasing & Validation")
    phasing_group.add_argument("--min-variant-frequency", type=float, default=0.10,
                               help="Minimum alternative allele frequency to call variant (default: 0.10 for 10%%)")
    phasing_group.add_argument("--min-variant-count", type=int, default=3,
                               help="Minimum alternative allele read count to call variant (default: 3)")
    phasing_group.add_argument("--significance-level", type=float, default=1e-5,
                               help="Significance level (alpha) for variant significance testing (default: 1e-5)")
    phasing_group.add_argument("--group-identity", type=float, default=0.85,
                               help="Minimum pairwise identity to group clusters for read reassignment, "
                                    "discard recovery, and CER validation. Grouping uses complete linkage: "
                                    "every pair within a group must meet this threshold. (default: 0.85)")
    phasing_group.add_argument("--hp-normalization-length", type=int, default=6,
                               help="Minimum homopolymer run length at/above which HP length variants "
                                    "are blanket-normalized (treated as noise). Runs of length >= this "
                                    "value have their length variants suppressed in MSA variant "
                                    "detection; runs shorter than this surface as candidates and are "
                                    "evaluated by context-aware CER. Default 6 matches the HP paper "
                                    "recommendation of CER-evaluating L <= 5 HP variants. (default: 6)")
    phasing_group.add_argument("--error-model", type=str, default="dorado-v5.0",
                               help="Per-basecaller error model used for context-aware variant "
                                    "validation. Either a shipped model name (use "
                                    "--list-error-models to see available), a user model in "
                                    "~/.config/speconsense/error_models/, or a filesystem path to "
                                    "a YAML file with a 'rates' mapping. (default: dorado-v5.0)")

    # Ambiguity Calling group
    ambiguity_group = parser.add_argument_group("Ambiguity Calling")
    ambiguity_group.add_argument("--min-ambiguity-frequency", type=float, default=0.10,
                                 help="Minimum alternative allele frequency for IUPAC ambiguity calling (default: 0.10 for 10%%)")
    ambiguity_group.add_argument("--min-ambiguity-count", type=int, default=3,
                                 help="Minimum alternative allele read count for IUPAC ambiguity calling (default: 3)")

    # Orientation group
    orient_group = parser.add_argument_group("Orientation")
    orient_group.add_argument("--orient-mode", choices=["none", "primer", "pyitsx", "pyitsx+primer"],
                              default="none",
                              help="Sequence orientation mode: none (default, no orientation), "
                                   "primer (orient via primer matching, discard failed), "
                                   "pyitsx (orient via ITS HMM profiles, discard failed/chimeric; "
                                   "requires pyitsx + ITSx), "
                                   "or pyitsx+primer (pyitsx with primer fallback for unrecognized reads)")
    orient_group.add_argument("--pyitsx-organism", default="F",
                              help="Organism group for pyitsx (default: F for Fungi). "
                                   "Used for --orient-mode pyitsx/pyitsx+primer and locus labeling in summarize. "
                                   "Common codes: F (Fungi), T (Tracheophyta), M (Metazoa)")

    # Performance group
    perf_group = parser.add_argument_group("Performance")
    perf_group.add_argument("--presample", type=int, default=1000,
                            help="Presample size for initial reads (default: 1000, 0 to disable)")
    perf_group.add_argument("--scale-threshold", type=int, default=1001,
                            help="Sequence count threshold for scalable mode (requires vsearch). "
                                 "Set to 0 to disable. Default: 1001")
    perf_group.add_argument("--threads", type=int, default=1, metavar="N",
                            help="Max threads for internal parallelism (vsearch, SPOA). "
                                 "0=auto-detect, default=1 (safe for parallel workflows).")

    # Debugging group
    debug_group = parser.add_argument_group("Debugging")
    debug_group.add_argument("--collect-discards", action="store_true",
                             help="Write discarded reads (outliers and filtered clusters) to cluster_debug/{sample}-discards.fastq")
    debug_group.add_argument("--no-collect-discards", action="store_false",
                             dest="collect_discards",
                             help="Override --collect-discards or profile setting")
    debug_group.add_argument("--log-level", default="INFO",
                             choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    # Advanced group (pre-tuned / experimental). Hidden from --help; users
    # discover via --help-advanced. Most of these existed to support paper-
    # side empirical analysis or expose pre-tuned algorithm internals.
    advanced_group = parser.add_argument_group(
        "Advanced (pre-tuned — rarely needed)",
        "Hidden from --help; view with --help-advanced. Pre-tuned algorithm "
        "internals, MAD knobs, and phase-disable flags retained primarily "
        "for empirical analysis.",
    )
    # MCL internals
    add_advanced_argument(advanced_group, "--inflation", type=float, default=4.0,
                          help="MCL inflation parameter (default: 4.0)")
    add_advanced_argument(advanced_group, "--k-nearest-neighbors", type=int, default=5,
                          help="Number of nearest neighbors for graph construction (default: 5)")
    # Phase disable flags (integral phases — disabling is for ablation studies)
    add_advanced_argument(advanced_group, "--disable-position-phasing", action="store_true",
                          help="Disable position-based variant phasing (enabled by default). "
                               "MCL graph clustering already separates most variants; this "
                               "second pass analyzes MSA positions to phase remaining variants.")
    add_advanced_argument(advanced_group, "--enable-position-phasing", action="store_false",
                          dest="disable_position_phasing",
                          help="Override --disable-position-phasing or profile setting")
    add_advanced_argument(advanced_group, "--disable-read-reassignment", action="store_true",
                          help="Disable post-phasing read reassignment within identity groups. "
                               "Reassignment moves reads between clusters based on consensus "
                               "concordance; disable to preserve original phasing-time membership.")
    add_advanced_argument(advanced_group, "--enable-read-reassignment", action="store_false",
                          dest="disable_read_reassignment",
                          help="Override --disable-read-reassignment or profile setting")
    add_advanced_argument(advanced_group, "--disable-discard-recovery", action="store_true",
                          help="Disable recovery of previously discarded reads into existing "
                               "clusters. Has no effect if --disable-read-reassignment is also "
                               "set (recovery requires reassignment).")
    add_advanced_argument(advanced_group, "--enable-discard-recovery", action="store_false",
                          dest="disable_discard_recovery",
                          help="Override --disable-discard-recovery or profile setting")
    add_advanced_argument(advanced_group, "--disable-second-phasing", action="store_true",
                          help="Disable the Phase 8 second phasing pass without disabling "
                               "the Phase 3 first pass. Has no effect if --disable-position-phasing "
                               "or --disable-read-reassignment is also set (Phase 8 is gated on both).")
    add_advanced_argument(advanced_group, "--enable-second-phasing", action="store_false",
                          dest="disable_second_phasing",
                          help="Override --disable-second-phasing or profile setting")
    add_advanced_argument(advanced_group, "--disable-noise-filter", action="store_true",
                          help="Disable the Phase 5 noise filter. Small clusters whose MSA "
                               "has no-majority columns are normally disbanded; this flag "
                               "preserves them and lets them flow through to final consensus.")
    add_advanced_argument(advanced_group, "--enable-noise-filter", action="store_false",
                          dest="disable_noise_filter",
                          help="Override --disable-noise-filter or profile setting")
    add_advanced_argument(advanced_group, "--disable-mad-outlier-removal", action="store_true",
                          help="Disable MAD-based outlier removal at final consensus generation. "
                               "Reads whose identity to the cluster consensus is more than n*MAD "
                               "from the median are normally dropped; this flag preserves them.")
    add_advanced_argument(advanced_group, "--enable-mad-outlier-removal", action="store_false",
                          dest="disable_mad_outlier_removal",
                          help="Override --disable-mad-outlier-removal or profile setting")
    add_advanced_argument(advanced_group, "--disable-ambiguity-calling", action="store_true",
                          help="Disable IUPAC ambiguity code calling for unphased variant positions")
    add_advanced_argument(advanced_group, "--enable-ambiguity-calling", action="store_false",
                          dest="disable_ambiguity_calling",
                          help="Override --disable-ambiguity-calling or profile setting")
    add_advanced_argument(advanced_group, "--disable-cluster-merging", action="store_true",
                          help="Disable merging of clusters with identical consensus sequences")
    add_advanced_argument(advanced_group, "--enable-cluster-merging", action="store_false",
                          dest="disable_cluster_merging",
                          help="Override --disable-cluster-merging or profile setting")
    add_advanced_argument(advanced_group, "--disable-homopolymer-equivalence", action="store_true",
                          help="Disable homopolymer equivalence in cluster merging (only merge identical sequences)")
    add_advanced_argument(advanced_group, "--enable-homopolymer-equivalence", action="store_false",
                          dest="disable_homopolymer_equivalence",
                          help="Override --disable-homopolymer-equivalence or profile setting")
    # MAD outlier removal tuning. Defaults mirror
    # ``speconsense.outliers.detect_rid_outliers``; once calibrated for the
    # release, end users shouldn't need to touch these.
    add_advanced_argument(advanced_group, "--mad-z-threshold", type=float, default=1.5,
                          help="Modified Z-score cutoff for the MAD rule. A read with "
                               "(0.6745 * (rid - median) / MAD) below -threshold is flagged. "
                               "Lower values flag more aggressively. The 0.8.1 default of 1.5 "
                               "is empirically tuned for speconsense's 3-10-read clusters; the "
                               "literature-standard 3.0 is too conservative for that regime "
                               "(see speconsense.outliers.detect_rid_outliers docstring). "
                               "(default: 1.5)")
    add_advanced_argument(advanced_group, "--mad-gap-factor", type=float, default=2.5,
                          help="Gap-rule multiplier. The worst read is flagged if "
                               "(r_second - r_worst) > gap_factor * (r_best - r_second). "
                               "Lower values flag more aggressively. (default: 2.5)")
    add_advanced_argument(advanced_group, "--mad-min-mad", type=float, default=0.002,
                          help="Floor for the MAD value to avoid divide-by-zero when most "
                               "reads have near-identical rid. Expressed on the 0..1 rid "
                               "scale (0.002 = 0.2pp). (default: 0.002)")
    add_advanced_argument(advanced_group, "--mad-min-drop-from-median", type=float, default=0.02,
                          help="Safety floor on the absolute drop below the cluster's median "
                               "rid. A statistically-unusual read is only flagged when its "
                               "rid is at least this far below the median. Expressed on the "
                               "0..1 rid scale (0.02 = 2pp). (default: 0.02)")

    # Version and profile options (default group)
    parser.add_argument("--version", action="version",
                        version=f"Speconsense {__version__}",
                        help="Show program's version number and exit")
    parser.add_argument("-p", "--profile", metavar="NAME",
                        help="Load parameter profile (use --list-profiles to see available)")
    parser.add_argument("--list-profiles", action="store_true",
                        help="List available profiles and exit")
    parser.add_argument("--show-profile", metavar="NAME",
                        help="Show contents of a named profile and exit")
    parser.add_argument("--list-error-models", action="store_true",
                        help="List available error models and exit")

    # Handle --list-profiles and --show-profile early (before requiring input_file)
    if '--list-profiles' in sys.argv:
        print_profiles_list('speconsense')
        sys.exit(0)

    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == '--show-profile' and i < len(sys.argv):
            show_profile(sys.argv[i + 1])
            sys.exit(0)
        if arg.startswith('--show-profile='):
            show_profile(arg.split('=', 1)[1])
            sys.exit(0)

    # Handle --list-error-models early (before requiring input_file)
    if '--list-error-models' in sys.argv:
        from speconsense.qctx import print_models_list
        print_models_list()
        sys.exit(0)

    # First pass: get profile name if specified
    # We need to detect which args were explicitly provided to not override them
    pre_args, _ = parser.parse_known_args()

    # Track which arguments were explicitly provided on CLI
    explicit_args = set()
    for arg in sys.argv[1:]:
        if arg.startswith('--') and '=' in arg:
            explicit_args.add(arg.split('=')[0][2:].replace('-', '_'))
        elif arg.startswith('--'):
            explicit_args.add(arg[2:].replace('-', '_'))
        elif arg.startswith('-') and len(arg) == 2:
            # Short option - would need to map to long name
            # For now, we skip this since profile args use long names
            pass

    # Load and apply profile if specified
    loaded_profile = None
    if pre_args.profile:
        try:
            loaded_profile = Profile.load(pre_args.profile)
        except ProfileError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

        # Apply profile values to parser defaults (explicit CLI args will override)
        for key, value in loaded_profile.speconsense.items():
            attr_name = key.replace('-', '_')
            if attr_name not in explicit_args:
                parser.set_defaults(**{attr_name: value})

    args = parser.parse_args()

    # Setup standard logging
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format=log_format
    )

    # Log profile usage after logging is configured
    if loaded_profile:
        logging.info(f"Using profile '{loaded_profile.name}': {loaded_profile.description}")

    # Resolve threads: 0 means auto-detect
    threads = args.threads if args.threads > 0 else os.cpu_count()

    sample = os.path.splitext(os.path.basename(args.input_file))[0]
    clusterer = SpecimenClusterer(
        min_identity=args.min_identity,
        inflation=args.inflation,
        min_size=args.min_size,
        max_sample_size=args.max_sample_size,
        presample_size=args.presample,
        k_nearest_neighbors=args.k_nearest_neighbors,
        sample_name=sample,
        disable_homopolymer_equivalence=args.disable_homopolymer_equivalence,
        disable_cluster_merging=args.disable_cluster_merging,
        output_dir=args.output_dir,
        enable_secondpass_phasing=not args.disable_position_phasing,
        enable_read_reassignment=not args.disable_read_reassignment,
        enable_discard_recovery=not args.disable_discard_recovery,
        enable_phase8=not args.disable_second_phasing,
        enable_noise_filter=not args.disable_noise_filter,
        enable_mad_outlier_removal=not args.disable_mad_outlier_removal,
        mad_z_threshold=args.mad_z_threshold,
        mad_gap_factor=args.mad_gap_factor,
        mad_min_mad=args.mad_min_mad,
        mad_min_drop_from_median=args.mad_min_drop_from_median,
        min_variant_frequency=args.min_variant_frequency,
        min_variant_count=args.min_variant_count,
        min_ambiguity_frequency=args.min_ambiguity_frequency,
        min_ambiguity_count=args.min_ambiguity_count,
        enable_iupac_calling=not args.disable_ambiguity_calling,
        scale_threshold=args.scale_threshold,
        max_threads=threads,
        collect_discards=args.collect_discards,
        significance_level=args.significance_level,
        group_identity=args.group_identity,
        min_hp_length=args.hp_normalization_length,
        error_model=args.error_model,
    )

    # Log configuration
    if not args.disable_position_phasing:
        logging.debug(f"Position-based variant phasing enabled: min_freq={args.min_variant_frequency:.0%}, "
                     f"min_count={args.min_variant_count}")

    # Set additional attributes for metadata
    clusterer.input_file = os.path.abspath(args.input_file)
    clusterer.algorithm = args.algorithm
    clusterer.orient_mode = args.orient_mode
    clusterer.pyitsx_organism = args.pyitsx_organism
    clusterer.profile_name = args.profile

    # Read primary sequences
    logging.info(f"Reading sequences from {args.input_file}")
    format = "fasta" if args.input_file.endswith(".fasta") else "fastq"
    records = list(SeqIO.parse(args.input_file, format))
    duplicate_ids = len(records) - len({r.id for r in records})
    if duplicate_ids:
        logging.info(f"Loaded {len(records)} primary reads ({len(records) - duplicate_ids} unique; {duplicate_ids} duplicate ID(s) will be deduplicated)")
    else:
        logging.info(f"Loaded {len(records)} primary reads")

    if len(records) == 0:
        logging.warning("No reads found in input file. Nothing to cluster.")
        sys.exit(0)

    clusterer.add_sequences(records)

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
    if args.orient_mode == "primer":
        if hasattr(clusterer, 'forward_primers') and hasattr(clusterer, 'reverse_primers'):
            failed_sequences = clusterer.orient_sequences()

            if failed_sequences:
                logging.info(f"Filtering out {len(failed_sequences)} reads with failed orientation")
                clusterer.discarded_read_ids.update(failed_sequences)
                for seq_id in failed_sequences:
                    del clusterer.sequences[seq_id]

                remaining = len(clusterer.sequences)
                logging.info(f"Continuing with {remaining} successfully oriented reads")
        else:
            logging.warning("--orient-mode=primer specified but no primers with position information loaded")

    elif args.orient_mode in ("pyitsx", "pyitsx+primer"):
        failed_ids, chimeric_ids = clusterer.orient_sequences_pyitsx(
            organism=args.pyitsx_organism,
        )

        # Primer fallback for pyitsx+primer mode
        # (chimeric reads are always discarded — that's a pyitsx-specific signal)
        if args.orient_mode == "pyitsx+primer" and failed_ids:
            has_primers = (hasattr(clusterer, 'forward_primers')
                           and hasattr(clusterer, 'reverse_primers'))
            if has_primers:
                primer_failed = clusterer.orient_sequences(seq_ids=failed_ids)
                primer_rescued = failed_ids - primer_failed
                if primer_rescued:
                    logging.info(f"Primer fallback rescued {len(primer_rescued)} of "
                                f"{len(failed_ids)} pyitsx-failed reads")
                failed_ids = primer_failed
            else:
                logging.warning("--orient-mode=pyitsx+primer specified but no primers "
                               "with position information loaded; skipping primer fallback")

        discard_ids = failed_ids | chimeric_ids
        if discard_ids:
            logging.info(f"Filtering out {len(discard_ids)} reads after orientation "
                        f"({len(failed_ids)} failed, {len(chimeric_ids)} chimeric)")
            clusterer.discarded_read_ids.update(discard_ids)
            for seq_id in discard_ids:
                del clusterer.sequences[seq_id]

        remaining = len(clusterer.sequences)
        logging.info(f"Continuing with {remaining} successfully oriented reads")

    clusterer.cluster(algorithm=args.algorithm)

    # Write metadata file (after clustering so per-cluster CER reproduction
    # data can be included for downstream tools).
    clusterer.write_metadata()
    print()

if __name__ == "__main__":
    main()
