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
import glob
import logging
import os
import sys
import subprocess
import tempfile
from typing import List, Tuple, Set, Optional, Dict
import numpy as np
from tqdm import tqdm
from Bio import SeqIO

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


def compute_global_statistics(
    all_cluster_stats: List[ClusterStats],
    all_cluster_alignments: List[Tuple[ClusterInfo, List[ReadAlignment]]]
) -> dict:
    """
    Compute global error statistics across all clusters.

    Computes statistics at two levels:
    1. Cluster-level: mean/median across cluster mean error rates
    2. Read-level: percentiles from ALL individual read error rates

    Args:
        all_cluster_stats: Cluster-level statistics
        all_cluster_alignments: Individual read alignments for all clusters

    Returns dictionary with global thresholds and distributions.
    """
    # Cluster-level statistics
    cluster_error_rates = [stats.mean_error_rate for stats in all_cluster_stats]
    cluster_edit_distances = [stats.mean_edit_distance for stats in all_cluster_stats]

    # Read-level statistics: collect ALL individual read error rates
    all_read_error_rates = []
    all_read_edit_distances = []

    for cluster_info, alignments in all_cluster_alignments:
        consensus_length = len(cluster_info.consensus_seq)
        if consensus_length == 0:
            continue

        for alignment in alignments:
            # Calculate per-read error rate
            error_rate = alignment.edit_distance / consensus_length
            all_read_error_rates.append(error_rate)
            all_read_edit_distances.append(alignment.edit_distance)

    global_stats = {
        # Cluster-level statistics
        'cluster_mean_error_rate': np.mean(cluster_error_rates),
        'cluster_median_error_rate': np.median(cluster_error_rates),
        'cluster_std_error_rate': np.std(cluster_error_rates),

        # Read-level statistics (for outlier detection)
        'read_mean_error_rate': np.mean(all_read_error_rates),
        'read_median_error_rate': np.median(all_read_error_rates),
        'read_std_error_rate': np.std(all_read_error_rates),
        'read_p95_error_rate': np.percentile(all_read_error_rates, 95),
        'read_p99_error_rate': np.percentile(all_read_error_rates, 99),

        # Edit distance statistics
        'mean_edit_distance': np.mean(cluster_edit_distances),
        'p95_edit_distance': np.percentile(cluster_edit_distances, 95),

        # Counts
        'num_clusters': len(all_cluster_stats),
        'num_reads': len(all_read_error_rates),
    }

    return global_stats


def map_msa_ids_to_read_ids(cluster_info: ClusterInfo, alignments: List[ReadAlignment]) -> Dict[str, str]:
    """
    Map MSA sequence IDs (seq0, seq1, etc.) to actual read IDs.

    The MSA uses generic seq0, seq1, ... IDs. We need to map these back to the
    actual read IDs from the sampled FASTQ file that was used to generate the MSA.

    Args:
        cluster_info: Cluster metadata
        alignments: Read alignments from MSA (with seq0, seq1, ... IDs)

    Returns:
        Dict mapping MSA ID (e.g., 'seq0') to actual read ID
    """
    # Find the sampled FASTQ file
    # Pattern: {sample_name}-RiC{N}-sampled.fastq
    cluster_debug_dir = os.path.dirname(cluster_info.reads_file)
    sampled_pattern = os.path.join(cluster_debug_dir, f"{cluster_info.sample_name}-RiC*-sampled.fastq")

    sampled_files = glob.glob(sampled_pattern)

    if not sampled_files:
        # No sampled file means all reads were used (cluster size <= max_sample_size)
        # In this case, use the reads file directly
        sampled_file = cluster_info.reads_file
    else:
        sampled_file = sampled_files[0]

    # Load sampled reads in order
    sampled_reads = list(SeqIO.parse(sampled_file, 'fastq'))

    # Create mapping: seq0 -> first read ID, seq1 -> second read ID, etc.
    msa_to_read_id = {}
    for i, read in enumerate(sampled_reads):
        msa_id = f"seq{i}"
        msa_to_read_id[msa_id] = read.id

    return msa_to_read_id


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
        Tuple of (keep_read_ids, outlier_read_ids) using actual read IDs
    """
    # Map MSA IDs (seq0, seq1, ...) to actual read IDs
    msa_to_read_id = map_msa_ids_to_read_ids(cluster_info, alignments)

    keep_reads = []
    outlier_reads = []

    consensus_length = len(cluster_info.consensus_seq)

    for alignment in alignments:
        error_rate = alignment.edit_distance / consensus_length if consensus_length > 0 else 0

        # Get actual read ID (not seqN)
        actual_read_id = msa_to_read_id.get(alignment.read_id, alignment.read_id)

        if error_rate <= global_threshold:
            keep_reads.append(actual_read_id)
        else:
            outlier_reads.append(actual_read_id)

    return keep_reads, outlier_reads


def run_spoa(sequences: List[str], max_sample_size: int = 100) -> Tuple[Optional[str], Optional[str]]:
    """
    Run SPOA to generate consensus sequence and MSA.

    Args:
        sequences: List of sequences to align
        max_sample_size: Maximum number of sequences to use

    Returns:
        Tuple of (consensus_sequence, msa_fasta) or (None, None) on failure
    """
    if not sequences:
        return None, None

    # Sample if too many sequences
    if len(sequences) > max_sample_size:
        import random
        sequences = random.sample(sequences, max_sample_size)

    # Write sequences to temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i}\n{seq}\n")
        temp_file = f.name

    try:
        # Run SPOA with MSA output
        result = subprocess.run(
            ['spoa', temp_file, '-r', '2'],  # -r 2: output MSA + consensus
            capture_output=True,
            text=True,
            check=True
        )

        # Parse output to extract consensus and MSA
        lines = result.stdout.strip().split('\n')
        if not lines:
            return None, None

        # Find consensus (last sequence with header containing "Consensus")
        consensus = None
        msa_lines = []
        current_seq = []
        current_header = None

        for line in lines:
            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    seq_str = ''.join(current_seq)
                    msa_lines.append(current_header)
                    msa_lines.append(seq_str)

                    if 'Consensus' in current_header:
                        consensus = seq_str.replace('-', '')  # Remove gaps from consensus

                # Start new sequence
                current_header = line
                current_seq = []
            else:
                current_seq.append(line.strip())

        # Don't forget last sequence
        if current_header is not None:
            seq_str = ''.join(current_seq)
            msa_lines.append(current_header)
            msa_lines.append(seq_str)

            if 'Consensus' in current_header:
                consensus = seq_str.replace('-', '')

        msa = '\n'.join(msa_lines) if msa_lines else None

        return consensus, msa

    except subprocess.CalledProcessError as e:
        logging.error(f"SPOA failed: {e}")
        return None, None
    except FileNotFoundError:
        logging.error("SPOA not found in PATH. Please install SPOA.")
        return None, None
    finally:
        os.unlink(temp_file)


def refine_cluster(
    cluster_info: ClusterInfo,
    outlier_read_ids: Set[str],
    output_dir: str,
    max_sample_size: int = 100
) -> Tuple[Optional[str], int, int]:
    """
    Refine a single cluster by removing outliers and regenerating consensus.

    Args:
        cluster_info: Original cluster information
        outlier_read_ids: Set of read IDs to remove
        output_dir: Directory for refined output
        max_sample_size: Maximum reads to use for consensus

    Returns:
        Tuple of (consensus_sequence, num_reads_kept, num_reads_removed)
    """
    # Load all reads from cluster file
    all_reads = list(SeqIO.parse(cluster_info.reads_file, 'fastq'))

    # Filter out outliers
    filtered_reads = [r for r in all_reads if r.id not in outlier_read_ids]

    num_kept = len(filtered_reads)
    num_removed = len(all_reads) - num_kept

    # Debug logging
    if num_removed > 0:
        logging.debug(f"{cluster_info.sample_name}: Removed {num_removed}/{len(all_reads)} reads (outlier_ids={len(outlier_read_ids)})")
    elif len(outlier_read_ids) > 0:
        logging.warning(f"{cluster_info.sample_name}: {len(outlier_read_ids)} outlier IDs provided but no reads removed! First outlier ID: {list(outlier_read_ids)[0] if outlier_read_ids else 'none'}, First read ID: {all_reads[0].id if all_reads else 'none'}")

    if num_kept == 0:
        logging.warning(f"Cluster {cluster_info.sample_name} has no reads after filtering!")
        return None, 0, num_removed

    # Create output directories
    cluster_debug_dir = os.path.join(output_dir, 'cluster_debug')
    os.makedirs(cluster_debug_dir, exist_ok=True)

    # Determine RiC (size after filtering)
    ric = min(num_kept, max_sample_size)

    # Write filtered reads file
    reads_filename = f"{cluster_info.sample_name}-RiC{ric}-reads.fastq"
    reads_file = os.path.join(cluster_debug_dir, reads_filename)
    with open(reads_file, 'w') as f:
        SeqIO.write(filtered_reads, f, 'fastq')

    # Sample reads if necessary
    if num_kept > max_sample_size:
        import random
        sampled_reads = random.sample(filtered_reads, max_sample_size)

        # Write sampled reads
        sampled_filename = f"{cluster_info.sample_name}-RiC{ric}-sampled.fastq"
        sampled_file = os.path.join(cluster_debug_dir, sampled_filename)
        with open(sampled_file, 'w') as f:
            SeqIO.write(sampled_reads, f, 'fastq')

        consensus_reads = sampled_reads
    else:
        consensus_reads = filtered_reads

    # Extract sequences for SPOA
    sequences = [str(r.seq) for r in consensus_reads]

    # Generate consensus and MSA
    consensus, msa = run_spoa(sequences, max_sample_size)

    if consensus is None:
        logging.warning(f"Failed to generate consensus for {cluster_info.sample_name}")
        return None, num_kept, num_removed

    # Write MSA file
    if msa:
        msa_filename = f"{cluster_info.sample_name}-RiC{ric}-msa.fasta"
        msa_file = os.path.join(cluster_debug_dir, msa_filename)
        with open(msa_file, 'w') as f:
            f.write(msa)

    # Write untrimmed consensus (we're not doing primer trimming in this version)
    untrimmed_filename = f"{cluster_info.sample_name}-RiC{ric}-untrimmed.fasta"
    untrimmed_file = os.path.join(cluster_debug_dir, untrimmed_filename)

    # Create header with metadata (no stability metrics)
    header = f">{cluster_info.sample_name} size={num_kept} ric={ric}"

    with open(untrimmed_file, 'w') as f:
        f.write(f"{header}\n{consensus}\n")

    return consensus, num_kept, num_removed


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
    global_stats = compute_global_statistics(all_cluster_stats, all_cluster_alignments)

    logging.info("\nGlobal Error Statistics:")
    logging.info(f"  Analyzed {global_stats['num_clusters']:,} clusters with {global_stats['num_reads']:,} total reads")
    logging.info(f"\n  Cluster-level statistics (mean of cluster means):")
    logging.info(f"    Mean error rate: {global_stats['cluster_mean_error_rate']*100:.2f}%")
    logging.info(f"    Median error rate: {global_stats['cluster_median_error_rate']*100:.2f}%")
    logging.info(f"    Std dev: {global_stats['cluster_std_error_rate']*100:.2f}%")
    logging.info(f"\n  Read-level statistics (all individual reads):")
    logging.info(f"    Mean error rate: {global_stats['read_mean_error_rate']*100:.2f}%")
    logging.info(f"    Median error rate: {global_stats['read_median_error_rate']*100:.2f}%")
    logging.info(f"    Std dev: {global_stats['read_std_error_rate']*100:.2f}%")
    logging.info(f"    P95: {global_stats['read_p95_error_rate']*100:.2f}%")
    logging.info(f"    P99: {global_stats['read_p99_error_rate']*100:.2f}%")

    # STEP 2: Identify outlier reads
    logging.info("\nSTEP 2: Identifying outlier reads...")
    outlier_threshold = global_stats['read_p95_error_rate']
    logging.info(f"  Using read-level P95 threshold: {outlier_threshold*100:.2f}%")

    total_reads = 0
    total_outliers = 0
    clusters_with_outliers = 0

    # Store outliers for each cluster
    cluster_outliers = {}  # Maps cluster_info to set of outlier read IDs

    for cluster_info, alignments in tqdm(all_cluster_alignments, desc="Finding outliers", unit="cluster"):
        keep_reads, outlier_reads = identify_outlier_reads(
            cluster_info,
            alignments,
            outlier_threshold
        )

        total_reads += len(alignments)
        total_outliers += len(outlier_reads)

        # Store outliers for this cluster
        cluster_outliers[cluster_info] = set(outlier_reads)

        if outlier_reads:
            clusters_with_outliers += 1
            logging.debug(f"{cluster_info.sample_name}: {len(outlier_reads)}/{len(alignments)} outliers, IDs={list(outlier_reads)[:3]}...")

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
    # (Will implement variant-based subdivision in future version)

    # STEP 5: Refine clusters by removing outliers and regenerating consensus
    if args.dry_run:
        logging.info("\nDRY RUN: Analysis complete, no files written")
    else:
        logging.info("\nSTEP 4: Refining clusters...")
        logging.info(f"  Removing outliers and regenerating consensus sequences")

        # Track refinement results
        refined_consensus_seqs = []  # List of (sample_name, consensus_seq, ric, size)
        total_reads_kept = 0
        total_reads_removed = 0
        clusters_failed = 0

        # Refine each cluster
        for cluster_info in tqdm(all_cluster_alignments, desc="Refining clusters", unit="cluster"):
            cluster_info_obj, alignments = cluster_info
            outlier_ids = cluster_outliers.get(cluster_info_obj, set())

            # Refine the cluster
            consensus, num_kept, num_removed = refine_cluster(
                cluster_info_obj,
                outlier_ids,
                args.output_dir,
                max_sample_size=100  # Match speconsense default
            )

            total_reads_kept += num_kept
            total_reads_removed += num_removed

            if consensus:
                refined_consensus_seqs.append((
                    cluster_info_obj.sample_name,
                    consensus,
                    min(num_kept, 100),  # RiC
                    num_kept  # Size
                ))
            else:
                clusters_failed += 1

        # Write combined consensus FASTA files (one per sample)
        logging.info(f"\n  Writing refined consensus sequences...")

        # Group by sample base name
        from collections import defaultdict
        samples = defaultdict(list)

        for sample_name, consensus, ric, size in refined_consensus_seqs:
            # Extract base sample name (before first '-c')
            base_name = sample_name.split('-c')[0] if '-c' in sample_name else sample_name
            samples[base_name].append((sample_name, consensus, ric, size))

        # Write one FASTA file per sample
        for base_name, consensus_list in samples.items():
            output_fasta = os.path.join(args.output_dir, f"{base_name}-all.fasta")

            with open(output_fasta, 'w') as f:
                for sample_name, consensus, ric, size in consensus_list:
                    header = f">{sample_name} size={size} ric={ric}"
                    f.write(f"{header}\n{consensus}\n")

        logging.info(f"\n  Refinement Results:")
        logging.info(f"    Clusters refined: {len(refined_consensus_seqs)}")
        logging.info(f"    Clusters failed: {clusters_failed}")
        logging.info(f"    Total reads kept: {total_reads_kept:,}")
        logging.info(f"    Total reads removed: {total_reads_removed:,}")
        logging.info(f"    Sample FASTA files written: {len(samples)}")

    logging.info("\n" + "="*80)
    logging.info("Refinement completed successfully!")
    logging.info("="*80)

    return 0


if __name__ == '__main__':
    sys.exit(main())
