#!/usr/bin/env python3
"""
Cluster refinement tool for Speconsense.

This tool analyzes error characteristics across all clusters from a speconsense run
and refines clusters by:
1. Identifying and removing outlier reads that don't belong to their clusters
2. Detecting unphased variants where clusters need to be subdivided

The tool uses global error statistics across all clusters to set thresholds,
making outlier detection more reliable than per-cluster analysis.
"""

import argparse
import logging
import os
import sys
from typing import List, Tuple
import numpy as np
from tqdm import tqdm

# Import analysis functions from analyze module
from speconsense.analyze import (
    load_consensus_sequences,
    find_cluster_files,
    find_msa_files,
    analyze_cluster,
    analyze_all_cluster_positions,
    ClusterInfo,
    ClusterStats,
    ReadAlignment
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="""
        Refine speconsense clusters by removing outlier reads and detecting unphased variants.

        This tool analyzes error characteristics across ALL clusters to establish global
        thresholds for outlier detection and variant phasing. It produces refined clusters
        that can be used as input to speconsense-summarize.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        'input_dir',
        help="Path to speconsense output directory"
    )

    parser.add_argument(
        '--output-dir', '-o',
        default=None,
        help="Directory for refined cluster output (default: {input_dir}_refined)"
    )

    parser.add_argument(
        '--min-ric',
        type=int,
        default=5,
        help="Minimum RiC threshold for analysis (default: 5)"
    )

    parser.add_argument(
        '--outlier-percentile',
        type=float,
        default=95.0,
        help="Percentile threshold for outlier detection (default: 95.0)"
    )

    parser.add_argument(
        '--variant-percentile',
        type=float,
        default=95.0,
        help="Percentile threshold for variant position detection (default: 95.0)"
    )

    parser.add_argument(
        '--use-sampled',
        action='store_true',
        help="Use sampled reads (*-sampled.fastq) instead of full cluster (*-reads.fastq)"
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help="Analyze and report refinements without writing output files"
    )

    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help="Logging level (default: INFO)"
    )

    return parser.parse_args()


def setup_logging(log_level: str):
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def compute_global_statistics(all_cluster_stats: List[ClusterStats]) -> dict:
    """
    Compute global error statistics across all clusters.

    Returns dictionary with global thresholds and distributions.
    """
    # Collect all error rates
    all_error_rates = [stats.mean_error_rate for stats in all_cluster_stats]

    # Collect all edit distances from cluster statistics
    all_mean_edit_distances = [stats.mean_edit_distance for stats in all_cluster_stats]

    global_stats = {
        'mean_error_rate': np.mean(all_error_rates),
        'median_error_rate': np.median(all_error_rates),
        'std_error_rate': np.std(all_error_rates),
        'p95_error_rate': np.percentile(all_error_rates, 95),
        'p99_error_rate': np.percentile(all_error_rates, 99),
        'mean_edit_distance': np.mean(all_mean_edit_distances),
        'p95_edit_distance': np.percentile(all_mean_edit_distances, 95),
    }

    return global_stats


def identify_outlier_reads(
    cluster_info: ClusterInfo,
    alignments: List[ReadAlignment],
    global_threshold: float
) -> Tuple[List[str], List[str]]:
    """
    Identify outlier reads that should be removed from the cluster.

    Args:
        cluster_info: Cluster metadata
        alignments: Read alignments for this cluster
        global_threshold: Global error rate threshold (e.g., p95)

    Returns:
        Tuple of (keep_read_ids, outlier_read_ids)
    """
    keep_reads = []
    outlier_reads = []

    consensus_length = len(cluster_info.consensus_seq)

    for alignment in alignments:
        error_rate = alignment.edit_distance / consensus_length if consensus_length > 0 else 0

        if error_rate <= global_threshold:
            keep_reads.append(alignment.read_id)
        else:
            outlier_reads.append(alignment.read_id)

    return keep_reads, outlier_reads


def main():
    """Main entry point for speconsense-refine."""
    args = parse_arguments()
    setup_logging(args.log_level)

    # Set output directory
    if args.output_dir is None:
        args.output_dir = args.input_dir.rstrip('/') + '_refined'

    logging.info("="*80)
    logging.info("Speconsense Cluster Refinement")
    logging.info("="*80)
    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Minimum RiC: {args.min_ric}")
    logging.info(f"Outlier percentile: {args.outlier_percentile}")
    logging.info(f"Variant percentile: {args.variant_percentile}")
    if args.dry_run:
        logging.info("DRY RUN MODE: No files will be written")
    logging.info("="*80)

    # Create output directory
    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    # Load consensus sequences
    logging.info("Loading consensus sequences...")
    consensus_map = load_consensus_sequences(args.input_dir)
    if not consensus_map:
        logging.error("No consensus sequences found")
        return 1
    logging.info(f"Found {len(consensus_map)} consensus sequences")

    # Find cluster FASTQ files
    logging.info("Finding cluster files...")
    cluster_files = find_cluster_files(args.input_dir, args.use_sampled)
    if not cluster_files:
        logging.error("No cluster files found")
        return 1
    logging.info(f"Found {len(cluster_files)} cluster files")

    # Find MSA files (required)
    logging.info("Finding MSA files...")
    msa_files = find_msa_files(args.input_dir)
    if not msa_files:
        logging.error("No MSA files found")
        logging.error("This tool requires MSA files generated by speconsense.")
        logging.error("Please re-run speconsense on your data to generate MSA files.")
        return 1
    logging.info(f"Found {len(msa_files)} MSA files")

    # Match consensus sequences with their reads files and MSA files
    clusters_to_analyze = []
    for key, reads_file in cluster_files.items():
        if key in consensus_map:
            msa_file = msa_files.get(key, None)
            if not msa_file:
                logging.warning(f"Skipping cluster {key} - no MSA file found")
                continue

            cluster_info = consensus_map[key]._replace(
                reads_file=reads_file,
                msa_file=msa_file
            )

            if cluster_info.ric >= args.min_ric:
                clusters_to_analyze.append(cluster_info)

    if not clusters_to_analyze:
        logging.error("No matching clusters found")
        return 1

    logging.info(f"Analyzing {len(clusters_to_analyze)} clusters (RiC >= {args.min_ric})")

    # STEP 1: Analyze all clusters to compute global statistics
    logging.info("\nSTEP 1: Computing global error statistics...")
    all_cluster_stats = []
    all_cluster_alignments = []

    for cluster_info in tqdm(clusters_to_analyze, desc="Analyzing clusters", unit="cluster"):
        stats, alignments = analyze_cluster(cluster_info)

        if stats:
            all_cluster_stats.append(stats)
            all_cluster_alignments.append((cluster_info, alignments))

    if not all_cluster_stats:
        logging.error("No clusters successfully analyzed")
        return 1

    # Compute global statistics
    global_stats = compute_global_statistics(all_cluster_stats)

    logging.info("\nGlobal Error Statistics:")
    logging.info(f"  Mean error rate: {global_stats['mean_error_rate']*100:.2f}%")
    logging.info(f"  Median error rate: {global_stats['median_error_rate']*100:.2f}%")
    logging.info(f"  Std dev error rate: {global_stats['std_error_rate']*100:.2f}%")
    logging.info(f"  P95 error rate: {global_stats['p95_error_rate']*100:.2f}%")
    logging.info(f"  P99 error rate: {global_stats['p99_error_rate']*100:.2f}%")

    # STEP 2: Identify outlier reads
    logging.info("\nSTEP 2: Identifying outlier reads...")
    outlier_threshold = global_stats['p95_error_rate']

    total_reads = 0
    total_outliers = 0
    clusters_with_outliers = 0

    for cluster_info, alignments in tqdm(all_cluster_alignments, desc="Finding outliers", unit="cluster"):
        keep_reads, outlier_reads = identify_outlier_reads(
            cluster_info,
            alignments,
            outlier_threshold
        )

        total_reads += len(alignments)
        total_outliers += len(outlier_reads)

        if outlier_reads:
            clusters_with_outliers += 1
            logging.debug(f"{cluster_info.sample_name}: {len(outlier_reads)}/{len(alignments)} outliers")

    logging.info(f"\nOutlier Detection Results:")
    logging.info(f"  Total reads: {total_reads:,}")
    logging.info(f"  Outliers identified: {total_outliers:,} ({total_outliers/total_reads*100:.1f}%)")
    logging.info(f"  Clusters with outliers: {clusters_with_outliers}/{len(all_cluster_alignments)}")

    # STEP 3: Positional analysis for variant detection
    logging.info("\nSTEP 3: Analyzing positional variation...")
    all_position_stats = analyze_all_cluster_positions(all_cluster_alignments)

    if all_position_stats:
        error_rates = [ps.error_rate * 100 for ps in all_position_stats]
        p95_positional = np.percentile(error_rates, args.variant_percentile)

        high_error_positions = [ps for ps in all_position_stats if ps.error_rate * 100 >= p95_positional]

        # Categorize by priority (exclude homopolymers and ends)
        high_priority = [
            ps for ps in high_error_positions
            if not ps.is_homopolymer
            and ps.position >= 20
            and ps.position < (ps.consensus_length - 20)
        ]

        logging.info(f"\nPositional Variation Results:")
        logging.info(f"  Total positions analyzed: {len(all_position_stats):,}")
        logging.info(f"  P{args.variant_percentile} positional error rate: {p95_positional:.2f}%")
        logging.info(f"  High-error positions: {len(high_error_positions):,}")
        logging.info(f"  High-priority (internal, non-HP): {len(high_priority)}")

        if high_priority:
            logging.info(f"\nTop 10 high-priority variant positions:")
            high_priority.sort(key=lambda x: x.error_rate, reverse=True)
            for i, ps in enumerate(high_priority[:10], 1):
                cluster_id = f"{ps.sample_name}-c{ps.cluster_num}"
                logging.info(f"    {i}. {cluster_id} pos {ps.position}: {ps.error_rate*100:.1f}% error")

    # TODO: STEP 4: Subdivide clusters based on high-priority variant positions
    # TODO: STEP 5: Write refined clusters to output directory

    if args.dry_run:
        logging.info("\nDRY RUN: Analysis complete, no files written")
    else:
        logging.info("\nTODO: Write refined cluster files")
        # Will implement cluster writing in future commits

    logging.info("\n" + "="*80)
    logging.info("Refinement analysis completed successfully!")
    logging.info("="*80)

    return 0


if __name__ == '__main__':
    sys.exit(main())
