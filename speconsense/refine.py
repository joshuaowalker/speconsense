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
import json
import logging
import os
import sys
import subprocess
import tempfile
from typing import List, Tuple, Set, Optional, Dict
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
import edlib

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
        '--min-alt-freq',
        type=float,
        default=0.20,
        help="Minimum alternative allele frequency for variant detection (default: 0.20)"
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


def read_metadata(input_dir: str, sample_name: str) -> Optional[Dict]:
    """
    Read metadata file for a sample.

    Args:
        input_dir: Speconsense output directory
        sample_name: Sample name (e.g., 'sample1-c1')

    Returns:
        Metadata dictionary or None if not found
    """
    # Extract base sample name (remove cluster suffix if present)
    # sample_name might be like "sample1-c1" but metadata is for "sample1"
    base_name = sample_name.split('-c')[0] if '-c' in sample_name else sample_name

    metadata_file = os.path.join(input_dir, 'cluster_debug', f'{base_name}-metadata.json')

    if not os.path.exists(metadata_file):
        logging.debug(f"Metadata file not found: {metadata_file}")
        return None

    try:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        logging.debug(f"Loaded metadata from {metadata_file}")
        return metadata
    except Exception as e:
        logging.warning(f"Failed to read metadata file {metadata_file}: {e}")
        return None


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


def trim_primers(sequence: str, primers: Dict[str, str]) -> Tuple[str, List[str]]:
    """
    Trim primers from start and end of sequence.

    Args:
        sequence: Sequence to trim
        primers: Dict mapping primer names to sequences

    Returns:
        Tuple of (trimmed_sequence, [found_primers])
    """
    if not primers:
        return sequence, []

    # Convert primer dict to list of tuples (matching core.py format)
    primer_list = list(primers.items())

    found_primers = []
    trimmed_seq = sequence

    # Look for primer at 5' end
    best_start_dist = float('inf')
    best_start_primer = None
    best_start_end = None

    for primer_name, primer_seq in primer_list:
        k = len(primer_seq) // 4  # Allow ~25% errors
        search_region = sequence[:len(primer_seq) * 2]

        result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)

        if result["editDistance"] != -1:
            dist = result["editDistance"]
            if dist < best_start_dist:
                best_start_dist = dist
                best_start_primer = primer_name
                best_start_end = result["locations"][0][1] + 1

    if best_start_primer:
        logging.debug(f"Found {best_start_primer} at 5' end (dist={best_start_dist})")
        found_primers.append(f"5'-{best_start_primer}")
        trimmed_seq = trimmed_seq[best_start_end:]

    # Look for primer at 3' end
    best_end_dist = float('inf')
    best_end_primer = None
    best_end_start = None

    for primer_name, primer_seq in primer_list:
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
        logging.debug(f"Found {best_end_primer} at 3' end (dist={best_end_dist})")
        found_primers.append(f"3'-{best_end_primer}")
        trimmed_seq = trimmed_seq[:best_end_start]

    return trimmed_seq, found_primers


def refine_cluster(
    cluster_info: ClusterInfo,
    outlier_read_ids: Set[str],
    output_dir: str,
    max_sample_size: int = 100,
    metadata: Optional[Dict] = None
) -> Tuple[Optional[str], int, int]:
    """
    Refine a single cluster by removing outliers and regenerating consensus.

    Args:
        cluster_info: Original cluster information
        outlier_read_ids: Set of read IDs to remove
        output_dir: Directory for refined output
        max_sample_size: Maximum reads to use for consensus (can be overridden by metadata)
        metadata: Optional metadata dict containing primers and other parameters

    Returns:
        Tuple of (consensus_sequence, num_reads_kept, num_reads_removed)
    """
    # Use max_sample_size from metadata if available
    if metadata and 'parameters' in metadata:
        params = metadata['parameters']
        if 'max_sample_size' in params and params['max_sample_size'] is not None:
            max_sample_size = params['max_sample_size']
            logging.debug(f"Using max_sample_size={max_sample_size} from metadata")
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

    # Write untrimmed consensus
    untrimmed_filename = f"{cluster_info.sample_name}-RiC{ric}-untrimmed.fasta"
    untrimmed_file = os.path.join(cluster_debug_dir, untrimmed_filename)

    # Create header with metadata (no stability metrics)
    header = f">{cluster_info.sample_name} size={num_kept} ric={ric}"

    with open(untrimmed_file, 'w') as f:
        f.write(f"{header}\n{consensus}\n")

    # Perform primer trimming if primers are available in metadata
    trimmed_consensus = consensus
    if metadata and 'primers' in metadata and metadata['primers']:
        primers = metadata['primers']
        trimmed_consensus, found_primers = trim_primers(consensus, primers)

        # Write trimmed consensus to main output file
        output_filename = f"{cluster_info.sample_name}.fasta"
        output_file = os.path.join(output_dir, output_filename)
        with open(output_file, 'a') as f:  # Append mode for multiple clusters per sample
            f.write(f"{header}\n{trimmed_consensus}\n")

        logging.debug(f"Trimmed {cluster_info.sample_name}: {len(consensus)} -> {len(trimmed_consensus)} bp, found primers: {found_primers}")
    else:
        # No primers available, just copy untrimmed consensus to output
        output_filename = f"{cluster_info.sample_name}.fasta"
        output_file = os.path.join(output_dir, output_filename)
        with open(output_file, 'a') as f:
            f.write(f"{header}\n{consensus}\n")

    return trimmed_consensus, num_kept, num_removed


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

    # Load metadata for each sample
    logging.info("Loading run metadata...")
    metadata_cache = {}  # Cache metadata by base sample name
    sample_names = set()
    for key in consensus_map.keys():
        sample_name, cluster_num = key
        # Extract base name (without cluster number)
        base_name = sample_name.split('-c')[0] if '-c' in sample_name else sample_name
        sample_names.add(base_name)

    for base_name in sample_names:
        metadata = read_metadata(args.input_dir, base_name)
        if metadata:
            metadata_cache[base_name] = metadata
            logging.debug(f"Loaded metadata for sample: {base_name}")
        else:
            logging.debug(f"No metadata found for sample: {base_name}")

    if metadata_cache:
        logging.info(f"Loaded metadata for {len(metadata_cache)} samples")
    else:
        logging.warning("No metadata files found - will use default parameters")

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

            # Get metadata for this cluster's sample
            base_name = cluster_info_obj.sample_name.split('-c')[0] if '-c' in cluster_info_obj.sample_name else cluster_info_obj.sample_name
            cluster_metadata = metadata_cache.get(base_name, None)

            # Refine the cluster
            consensus, num_kept, num_removed = refine_cluster(
                cluster_info_obj,
                outlier_ids,
                args.output_dir,
                max_sample_size=100,  # Default, will be overridden by metadata if available
                metadata=cluster_metadata
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

        # STEP 5: Variant position analysis on refined clusters
        logging.info("\nSTEP 5: Analyzing variant positions in refined clusters...")

        # Import additional functions
        from speconsense.analyze import (
            analyze_positional_variation,
            calculate_min_coverage_for_variant_detection,
            is_variant_position_with_composition,
            extract_alignments_from_msa
        )

        # Re-load refined MSA files and consensus sequences
        refined_msa_files = find_msa_files(args.output_dir)
        refined_consensus_map = load_consensus_sequences(args.output_dir)

        if not refined_msa_files or not refined_consensus_map:
            logging.warning("  No refined MSA files found - skipping variant analysis")
        else:
            # Calculate minimum coverage threshold
            min_coverage = calculate_min_coverage_for_variant_detection(
                min_alt_freq=args.min_alt_freq if hasattr(args, 'min_alt_freq') else 0.20,
                global_error_rate=global_stats['read_p95_error_rate'],
                confidence_sigma=3.0
            )

            logging.info(f"  Minimum coverage for variant detection: {min_coverage} reads")
            logging.info(f"  Minimum alternative allele frequency: {args.min_alt_freq if hasattr(args, 'min_alt_freq') else 0.20:.1%}")

            # Analyze each refined cluster
            all_variant_positions = []

            for key in refined_consensus_map.keys():
                if key not in refined_msa_files:
                    continue

                cluster_info = refined_consensus_map[key]
                msa_file = refined_msa_files[key]

                # Skip low-coverage clusters
                if cluster_info.ric < min_coverage:
                    continue

                # Extract alignments from MSA
                alignments, consensus = extract_alignments_from_msa(msa_file)

                if not alignments:
                    continue

                # Analyze positional variation with composition tracking
                position_stats = analyze_positional_variation(
                    alignments,
                    consensus,
                    global_stats['read_p95_error_rate'],
                    msa_file=msa_file
                )

                # Identify variant positions
                for pos_stat in position_stats:
                    is_variant, variant_bases, reason = is_variant_position_with_composition(
                        pos_stat,
                        global_stats['read_p95_error_rate'],
                        min_alt_freq=args.min_alt_freq if hasattr(args, 'min_alt_freq') else 0.20,
                        confidence_sigma=3.0
                    )

                    if is_variant:
                        all_variant_positions.append({
                            'sample_name': cluster_info.sample_name,
                            'cluster_num': cluster_info.cluster_num,
                            'position': pos_stat.position,
                            'coverage': pos_stat.coverage,
                            'variant_bases': variant_bases,
                            'base_composition': pos_stat.base_composition,
                            'error_rate': pos_stat.error_rate,
                            'is_homopolymer': pos_stat.is_homopolymer,
                            'reason': reason
                        })

            # Report variant positions
            logging.info(f"\n  Variant Position Detection Results:")
            logging.info(f"    Total variant positions found: {len(all_variant_positions)}")

            if all_variant_positions:
                # Filter for high-priority variants (internal, non-homopolymer)
                high_priority_variants = [
                    v for v in all_variant_positions
                    if not v['is_homopolymer']
                    and v['position'] >= 20
                ]

                logging.info(f"    High-priority variants (non-HP, internal): {len(high_priority_variants)}")

                if high_priority_variants:
                    # Sort by error rate (highest first)
                    high_priority_variants.sort(key=lambda x: x['error_rate'], reverse=True)

                    logging.info(f"\n  Top 20 high-priority variant positions:")
                    for i, v in enumerate(high_priority_variants[:20], 1):
                        cluster_id = f"{v['sample_name']}"
                        comp = v['base_composition']
                        sorted_comp = sorted(comp.items(), key=lambda x: x[1], reverse=True)
                        comp_str = ', '.join([f"{b}:{c}" for b, c in sorted_comp[:3]])

                        logging.info(
                            f"    {i}. {cluster_id} pos {v['position']}: "
                            f"{v['error_rate']*100:.1f}% error, "
                            f"cov={v['coverage']}, "
                            f"bases=[{comp_str}], "
                            f"variants={v['variant_bases']}"
                        )

    logging.info("\n" + "="*80)
    logging.info("Refinement completed successfully!")
    logging.info("="*80)

    return 0


if __name__ == '__main__':
    sys.exit(main())
