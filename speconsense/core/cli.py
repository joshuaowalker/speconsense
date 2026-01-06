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
)

from .clusterer import SpecimenClusterer


def main():
    parser = argparse.ArgumentParser(
        description="MCL-based clustering of nanopore amplicon reads"
    )

    # Input/Output group
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("input_file", help="Input FASTQ file")
    io_group.add_argument("-O", "--output-dir", default="clusters",
                          help="Output directory for all files (default: clusters)")
    io_group.add_argument("--primers", help="FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)")
    io_group.add_argument("--augment-input", help="Additional FASTQ/FASTA file with sequences recovered after primary demultiplexing (e.g., from specimine)")

    # Clustering group
    clustering_group = parser.add_argument_group("Clustering")
    clustering_group.add_argument("--algorithm", type=str, default="graph", choices=["graph", "greedy"],
                                  help="Clustering algorithm to use (default: graph)")
    clustering_group.add_argument("--min-identity", type=float, default=0.9,
                                  help="Minimum sequence identity threshold for clustering (default: 0.9)")
    clustering_group.add_argument("--inflation", type=float, default=4.0,
                                  help="MCL inflation parameter (default: 4.0)")
    clustering_group.add_argument("--k-nearest-neighbors", type=int, default=5,
                                  help="Number of nearest neighbors for graph construction (default: 5)")

    # Filtering group
    filtering_group = parser.add_argument_group("Filtering")
    filtering_group.add_argument("--min-size", type=int, default=5,
                                 help="Minimum cluster size (default: 5, 0 to disable)")
    filtering_group.add_argument("--min-cluster-ratio", type=float, default=0.01,
                                 help="Minimum size ratio between a cluster and the largest cluster (default: 0.01, 0 to disable)")
    filtering_group.add_argument("--max-sample-size", type=int, default=100,
                                 help="Maximum cluster size for consensus (default: 100)")
    filtering_group.add_argument("--outlier-identity", type=float, default=None,
                                 help="Minimum read-to-consensus identity to keep a read (default: auto). "
                                      "Reads below this threshold are removed as outliers before final "
                                      "consensus generation. Auto-calculated as (1 + min_identity) / 2. "
                                      "This threshold is typically higher than --min-identity because "
                                      "the consensus is error-corrected through averaging.")

    # Variant Phasing group
    phasing_group = parser.add_argument_group("Variant Phasing")
    phasing_group.add_argument("--disable-position-phasing", action="store_true",
                               help="Disable position-based variant phasing (enabled by default). "
                                    "MCL graph clustering already separates most variants; this "
                                    "second pass analyzes MSA positions to phase remaining variants.")
    phasing_group.add_argument("--min-variant-frequency", type=float, default=0.10,
                               help="Minimum alternative allele frequency to call variant (default: 0.10 for 10%%)")
    phasing_group.add_argument("--min-variant-count", type=int, default=5,
                               help="Minimum alternative allele read count to call variant (default: 5)")

    # Ambiguity Calling group
    ambiguity_group = parser.add_argument_group("Ambiguity Calling")
    ambiguity_group.add_argument("--disable-ambiguity-calling", action="store_true",
                                 help="Disable IUPAC ambiguity code calling for unphased variant positions")
    ambiguity_group.add_argument("--min-ambiguity-frequency", type=float, default=0.10,
                                 help="Minimum alternative allele frequency for IUPAC ambiguity calling (default: 0.10 for 10%%)")
    ambiguity_group.add_argument("--min-ambiguity-count", type=int, default=3,
                                 help="Minimum alternative allele read count for IUPAC ambiguity calling (default: 3)")

    # Cluster Merging group
    merging_group = parser.add_argument_group("Cluster Merging")
    merging_group.add_argument("--disable-cluster-merging", action="store_true",
                               help="Disable merging of clusters with identical consensus sequences")
    merging_group.add_argument("--disable-homopolymer-equivalence", action="store_true",
                               help="Disable homopolymer equivalence in cluster merging (only merge identical sequences)")

    # Orientation group
    orient_group = parser.add_argument_group("Orientation")
    orient_group.add_argument("--orient-mode", choices=["skip", "keep-all", "filter-failed"], default="skip",
                              help="Sequence orientation mode: skip (default, no orientation), keep-all (orient but keep failed), or filter-failed (orient and remove failed)")

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
    perf_group.add_argument("--enable-early-filter", action="store_true",
                            help="Enable early filtering to skip small clusters before variant phasing (improves performance for large datasets)")

    # Debugging group
    debug_group = parser.add_argument_group("Debugging")
    debug_group.add_argument("--collect-discards", action="store_true",
                             help="Write discarded reads (outliers and filtered clusters) to cluster_debug/{sample}-discards.fastq")
    debug_group.add_argument("--log-level", default="INFO",
                             choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    # Version and profile options (default group)
    parser.add_argument("--version", action="version",
                        version=f"Speconsense {__version__}",
                        help="Show program's version number and exit")
    parser.add_argument("-p", "--profile", metavar="NAME",
                        help="Load parameter profile (use --list-profiles to see available)")
    parser.add_argument("--list-profiles", action="store_true",
                        help="List available profiles and exit")

    # Handle --list-profiles early (before requiring input_file)
    if '--list-profiles' in sys.argv:
        print_profiles_list('speconsense')
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
        scale_threshold=args.scale_threshold,
        max_threads=threads,
        early_filter=args.enable_early_filter,
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

                # Track as discarded and remove from clustering (but keep records for discards file)
                clusterer.discarded_read_ids.update(failed_sequences)
                for seq_id in failed_sequences:
                    del clusterer.sequences[seq_id]
                    # Keep records so they can be written to discards file

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
