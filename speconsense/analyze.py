#!/usr/bin/env python3

"""
Empirical error analysis for Speconsense output.

This module compares raw reads against their consensus sequences to estimate
real-world Nanopore error rates and characteristics. It adapts the synthetic
testing framework to work with real data, where consensus sequences serve as
proxies for biological truth.
"""

import os
import re
import glob
import csv
import logging
from io import StringIO
from typing import List, Dict, Tuple, Optional, NamedTuple, Union
from pathlib import Path
from collections import defaultdict, Counter

from Bio import SeqIO
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import binomtest
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from tqdm import tqdm

# Import shared types and functions from core
from speconsense.core import (
    ErrorPosition,
    ReadAlignment,
    PositionStats,
    extract_alignments_from_msa,
    detect_homopolymer_run,
    calculate_min_coverage_for_variant_detection,
    is_variant_position_with_composition,
)


class ClusterInfo(NamedTuple):
    """Information about a cluster and its consensus."""
    sample_name: str
    cluster_num: int
    consensus_seq: str
    reads_file: str  # Path to FASTQ file with raw reads
    ric: int  # Reads in consensus (from header)
    size: int  # Total cluster size (from header)
    rid: Optional[float] = None  # Mean read identity
    rid_min: Optional[float] = None  # Minimum read identity
    msa_file: Optional[str] = None  # Path to MSA file (if available)


# Note: ErrorPosition, ReadAlignment, PositionStats now imported from core.py
# Wrapper functions below provide file I/O support for the CLI tool


class ClusterStats(NamedTuple):
    """Statistics for a single cluster."""
    sample_name: str
    cluster_num: int
    ric: int
    size: int
    actual_cluster_size: int  # Actual number of reads in the cluster file
    num_reads_analyzed: int
    mean_edit_distance: float
    median_edit_distance: float
    p95_edit_distance: float
    mean_error_rate: float
    median_error_rate: float
    insertion_rate: float
    deletion_rate: float
    substitution_rate: float
    consensus_length: int
    mean_read_length: float
    has_contiguous_errors: bool  # True if reads show contiguous high-error regions
    rid: Optional[float] = None  # Mean read identity from consensus header
    rid_min: Optional[float] = None  # Minimum read identity from consensus header


def parse_consensus_header(header: str) -> Tuple[str, int, int, int, Optional[float], Optional[float]]:
    """
    Parse a consensus FASTA header to extract metadata.

    Header format: >sample-c1 size=200 ric=100 [rid=99.3] [rid_min=97.5]

    Returns:
        Tuple of (sample_name, cluster_num, size, ric, rid, rid_min)
    """
    # Remove '>' if present
    header = header.lstrip('>')

    # Extract sample name (everything before first space)
    parts = header.split()
    sample_name = parts[0]

    # Extract cluster number from sample name
    match = re.search(r'-c(\d+)$', sample_name)
    if not match:
        raise ValueError(f"Cannot extract cluster number from: {sample_name}")
    cluster_num = int(match.group(1))

    # Parse key=value pairs
    size = ric = None
    rid = rid_min = None

    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            if key == 'size':
                size = int(value)
            elif key == 'ric':
                ric = int(value)
            elif key == 'rid':
                rid = float(value) / 100.0  # Convert percentage to fraction
            elif key == 'rid_min':
                rid_min = float(value) / 100.0  # Convert percentage to fraction

    if size is None or ric is None:
        raise ValueError(f"Missing size or ric in header: {header}")

    return sample_name, cluster_num, size, ric, rid, rid_min


def find_cluster_files(output_dir: str, use_sampled: bool = False) -> Dict[Tuple[str, int], str]:
    """
    Find all cluster FASTQ files in the output directory.

    Args:
        output_dir: Path to speconsense output directory
        use_sampled: If True, use *-sampled.fastq files; otherwise use *-reads.fastq

    Returns:
        Dict mapping (sample_base_name, cluster_num) to FASTQ file path
    """
    cluster_debug_dir = os.path.join(output_dir, 'cluster_debug')
    if not os.path.exists(cluster_debug_dir):
        logging.warning(f"cluster_debug directory not found: {cluster_debug_dir}")
        return {}

    # Pattern for cluster files: {sample}-c{N}-RiC{M}-{stage}.fastq
    suffix = 'sampled.fastq' if use_sampled else 'reads.fastq'
    pattern = os.path.join(cluster_debug_dir, f'*-c*-RiC*-{suffix}')

    cluster_files = {}
    for file_path in glob.glob(pattern):
        basename = os.path.basename(file_path)

        # Extract sample name and cluster number
        # Pattern: {sample}-c{num}-RiC{ric}-{stage}.fastq
        match = re.match(r'(.+)-c(\d+)-RiC\d+-.*\.fastq$', basename)
        if not match:
            logging.warning(f"Cannot parse cluster file name: {basename}")
            continue

        sample_name = match.group(1)
        cluster_num = int(match.group(2))

        cluster_files[(sample_name, cluster_num)] = file_path

    return cluster_files


def find_msa_files(output_dir: str) -> Dict[Tuple[str, int], str]:
    """
    Find all MSA FASTA files in the output directory.

    Args:
        output_dir: Path to speconsense output directory

    Returns:
        Dict mapping (sample_base_name, cluster_num) to MSA file path
    """
    cluster_debug_dir = os.path.join(output_dir, 'cluster_debug')
    if not os.path.exists(cluster_debug_dir):
        logging.warning(f"cluster_debug directory not found: {cluster_debug_dir}")
        return {}

    # Pattern for MSA files: {sample}-c{N}-RiC{M}-msa.fasta
    pattern = os.path.join(cluster_debug_dir, '*-c*-RiC*-msa.fasta')

    msa_files = {}
    for file_path in glob.glob(pattern):
        basename = os.path.basename(file_path)

        # Extract sample name and cluster number
        # Pattern: {sample}-c{num}-RiC{ric}-msa.fasta
        match = re.match(r'(.+)-c(\d+)-RiC\d+-msa\.fasta$', basename)
        if not match:
            logging.warning(f"Cannot parse MSA file name: {basename}")
            continue

        sample_name = match.group(1)
        cluster_num = int(match.group(2))

        msa_files[(sample_name, cluster_num)] = file_path

    return msa_files


def load_consensus_sequences(output_dir: str) -> Dict[Tuple[str, int], ClusterInfo]:
    """
    Load consensus sequences from speconsense output FASTA files.

    IMPORTANT: Loads UNTRIMMED consensus sequences from cluster_debug directory.
    This is critical for apples-to-apples comparison with untrimmed raw reads.
    Trimmed consensus would cause primer regions to appear as errors.

    Args:
        output_dir: Path to speconsense output directory

    Returns:
        Dict mapping (sample_base_name, cluster_num) to ClusterInfo
    """
    consensus_map = {}

    # Find all *-untrimmed.fasta files in cluster_debug directory
    # These contain consensus BEFORE primer trimming, matching the raw reads
    cluster_debug_dir = os.path.join(output_dir, 'cluster_debug')
    if not os.path.exists(cluster_debug_dir):
        logging.warning(f"cluster_debug directory not found: {cluster_debug_dir}")
        return {}

    pattern = os.path.join(cluster_debug_dir, '*-untrimmed.fasta')

    for fasta_file in glob.glob(pattern):
        try:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                sample_name, cluster_num, size, ric, rid, rid_min = parse_consensus_header(record.description)

                # Extract base sample name (remove -c{num} suffix)
                base_name = re.sub(r'-c\d+$', '', sample_name)

                consensus_map[(base_name, cluster_num)] = ClusterInfo(
                    sample_name=sample_name,
                    cluster_num=cluster_num,
                    consensus_seq=str(record.seq),
                    reads_file='',  # Will be filled in later
                    ric=ric,
                    size=size,
                    rid=rid,
                    rid_min=rid_min
                )
        except Exception as e:
            logging.warning(f"Error parsing {fasta_file}: {e}")
            continue

    return consensus_map


def extract_alignments_from_msa_file(msa_input: Union[str, Path]) -> Tuple[List[ReadAlignment], str, Dict[int, Optional[int]]]:
    """
    Extract read alignments from an MSA file or string.

    This is a wrapper around core.extract_alignments_from_msa() that adds file I/O support.

    Args:
        msa_input: Path to MSA FASTA file OR MSA content as string

    Returns:
        Tuple of:
        - list of ReadAlignment objects (with errors at MSA positions)
        - consensus sequence without gaps
        - mapping from MSA position to consensus position (None for insertion columns)
    """
    # Determine if input is MSA content (string with newlines) or file path
    if isinstance(msa_input, str) and '\n' in msa_input:
        # Input is MSA content as string - call core function directly
        from speconsense.core import extract_alignments_from_msa as core_extract
        return core_extract(msa_input)
    else:
        # Input is file path - read it first
        with open(msa_input, 'r') as f:
            msa_string = f.read()
        from speconsense.core import extract_alignments_from_msa as core_extract
        return core_extract(msa_string)




def detect_contiguous_error_regions(
    error_positions: List[ErrorPosition],
    consensus_length: int,
    msa_to_consensus_pos: Dict[int, Optional[int]],
    window_size: int = 50
) -> bool:
    """
    Detect if errors cluster into contiguous high-error regions.

    This helps identify read-end degradation and systematic error patterns.

    Args:
        error_positions: List of ErrorPosition objects (with MSA positions)
        consensus_length: Length of consensus sequence
        msa_to_consensus_pos: Mapping from MSA position to consensus position
        window_size: Size of sliding window for error density calculation

    Returns:
        True if contiguous high-error region detected
    """
    if len(error_positions) < 5:
        return False

    # Map MSA positions to consensus positions, filtering out insertions
    positions = []
    for ep in error_positions:
        cons_pos = msa_to_consensus_pos.get(ep.msa_position)
        if cons_pos is not None:
            positions.append(cons_pos)

    if len(positions) < 5:
        return False

    # Calculate error density in sliding windows
    error_counts = Counter(positions)

    # Check windows
    for start in range(0, consensus_length - window_size + 1, 10):
        end = start + window_size
        errors_in_window = sum(error_counts[pos] for pos in range(start, end))

        # Consider high-error if >20% of window has errors
        if errors_in_window / window_size > 0.2:
            return True

    return False


# detect_homopolymer_run is now imported from core
# Keeping this comment to mark where it was

def _placeholder_detect_homopolymer_run(seq: str, position: int, min_length: int = 4) -> bool:
    """
    Check if a position is within a homopolymer run.

    Args:
        seq: DNA sequence
        position: Position to check (0-indexed)
        min_length: Minimum length of homopolymer to detect

    Returns:
        True if position is within a homopolymer run of min_length or longer
    """
    if position >= len(seq):
        return False

    base = seq[position]

    # Find start of run
    start = position
    while start > 0 and seq[start - 1] == base:
        start -= 1

    # Find end of run
    end = position
    while end < len(seq) - 1 and seq[end + 1] == base:
        end += 1

    run_length = end - start + 1
    return run_length >= min_length


class ClusterPositionStats(NamedTuple):
    """Statistics for a position within a specific cluster."""
    sample_name: str
    cluster_num: int
    position: int
    consensus_length: int
    error_rate: float
    coverage: int
    error_count: int
    sub_count: int
    ins_count: int
    del_count: int
    is_homopolymer: bool
    nucleotide: str


def analyze_all_cluster_positions(
    cluster_alignments: List[Tuple[ClusterInfo, List[ReadAlignment]]]
) -> List[ClusterPositionStats]:
    """
    Analyze positional variation across all clusters.

    NOTE: This is legacy code that works in consensus space, not MSA space.
    For variant detection, use analyze_positional_variation instead.

    Args:
        cluster_alignments: List of (cluster_info, alignments) tuples

    Returns:
        List of ClusterPositionStats for every position in every cluster
    """
    all_position_stats = []

    for cluster_info, alignments in tqdm(cluster_alignments, desc="Analyzing positions", unit="cluster"):
        if not alignments:
            continue

        consensus_seq = cluster_info.consensus_seq
        consensus_len = len(consensus_seq)

        # Get MSA-to-consensus mapping from MSA file
        if not cluster_info.msa_file or not os.path.exists(cluster_info.msa_file):
            logging.warning(f"MSA file not found for {cluster_info.sample_name}, skipping positional analysis")
            continue

        try:
            # Parse MSA to get mapping
            _, _, msa_to_consensus_pos = extract_alignments_from_msa_file(cluster_info.msa_file)
        except Exception as e:
            logging.warning(f"Failed to parse MSA for {cluster_info.sample_name}: {e}")
            continue

        # Build error frequency matrix for this cluster (in consensus space)
        error_matrix = np.zeros((consensus_len, 4), dtype=int)

        for alignment in alignments:
            # Count coverage for all consensus positions
            for pos in range(consensus_len):
                error_matrix[pos, 3] += 1

            # Map MSA position errors to consensus positions
            for error_pos in alignment.error_positions:
                msa_pos = error_pos.msa_position
                cons_pos = msa_to_consensus_pos.get(msa_pos)

                # Only count errors that map to a consensus position (skip insertions)
                if cons_pos is not None and 0 <= cons_pos < consensus_len:
                    if error_pos.error_type == 'sub':
                        error_matrix[cons_pos, 0] += 1
                    elif error_pos.error_type == 'ins':
                        # Insertions don't really map to consensus positions, skip
                        pass
                    elif error_pos.error_type == 'del':
                        error_matrix[cons_pos, 2] += 1

        # Calculate statistics for each consensus position
        for pos in range(consensus_len):
            sub_count = error_matrix[pos, 0]
            ins_count = error_matrix[pos, 1]
            del_count = error_matrix[pos, 2]
            coverage = error_matrix[pos, 3]

            error_count = sub_count + ins_count + del_count
            error_rate = error_count / coverage if coverage > 0 else 0.0

            # Sequence context
            is_homopolymer = detect_homopolymer_run(consensus_seq, pos, min_length=4)
            nucleotide = consensus_seq[pos] if pos < len(consensus_seq) else 'N'

            all_position_stats.append(ClusterPositionStats(
                sample_name=cluster_info.sample_name,
                cluster_num=cluster_info.cluster_num,
                position=pos,
                consensus_length=consensus_len,
                error_rate=error_rate,
                coverage=coverage,
                error_count=error_count,
                sub_count=sub_count,
                ins_count=ins_count,
                del_count=del_count,
                is_homopolymer=is_homopolymer,
                nucleotide=nucleotide
            ))

    return all_position_stats


def analyze_cluster(cluster_info: ClusterInfo) -> Tuple[ClusterStats, List[ReadAlignment]]:
    """
    Analyze all reads in a cluster against their consensus using MSA.

    Args:
        cluster_info: Cluster information including MSA file and consensus

    Returns:
        Tuple of (cluster statistics, list of read alignments)
    """
    if not os.path.exists(cluster_info.reads_file):
        logging.warning(f"Reads file not found: {cluster_info.reads_file}")
        return None, []

    # Require MSA file
    if not cluster_info.msa_file or not os.path.exists(cluster_info.msa_file):
        logging.error(f"MSA file not found for {cluster_info.sample_name}: {cluster_info.msa_file}")
        logging.error(f"  This tool requires MSA files generated by speconsense.")
        logging.error(f"  Please re-run speconsense on your data to generate MSA files.")
        return None, []

    consensus_seq = cluster_info.consensus_seq
    total_reads_in_file = 0

    # Count total reads in cluster file
    for record in SeqIO.parse(cluster_info.reads_file, 'fastq'):
        total_reads_in_file += 1

    # Extract alignments from MSA
    logging.debug(f"Using MSA for {cluster_info.sample_name}: {cluster_info.msa_file}")
    alignments, consensus_from_msa, msa_to_consensus_pos = extract_alignments_from_msa_file(cluster_info.msa_file)

    # Verify consensus matches
    if consensus_from_msa and consensus_from_msa != consensus_seq:
        # Show detailed mismatch information
        logging.warning(f"Consensus mismatch for {cluster_info.sample_name}:")
        logging.warning(f"  MSA consensus:      length={len(consensus_from_msa)}, "
                      f"start={consensus_from_msa[:50]}..., end=...{consensus_from_msa[-50:]}")
        logging.warning(f"  Expected consensus: length={len(consensus_seq)}, "
                      f"start={consensus_seq[:50]}..., end=...{consensus_seq[-50:]}")

        # Check if one is a substring of the other (primer trimming issue?)
        if consensus_from_msa in consensus_seq:
            logging.warning(f"  → MSA consensus is a substring of expected (primers may have been trimmed from MSA)")
        elif consensus_seq in consensus_from_msa:
            logging.warning(f"  → Expected consensus is a substring of MSA (primers may not be trimmed from MSA)")

        # Use consensus from MSA as it's the one that was actually aligned
        consensus_seq = consensus_from_msa

    if not alignments:
        logging.warning(f"No alignments found in MSA for {cluster_info.sample_name}")
        return None, []

    # Calculate statistics
    edit_distances = [a.edit_distance for a in alignments]
    read_lengths = [a.read_length for a in alignments]
    error_rates = [a.edit_distance / len(consensus_seq) for a in alignments]

    total_insertions = sum(a.num_insertions for a in alignments)
    total_deletions = sum(a.num_deletions for a in alignments)
    total_substitutions = sum(a.num_substitutions for a in alignments)
    total_errors = total_insertions + total_deletions + total_substitutions

    # Check for contiguous error regions across all reads
    all_error_positions = []
    for a in alignments:
        all_error_positions.extend(a.error_positions)
    has_contiguous = detect_contiguous_error_regions(
        all_error_positions,
        len(consensus_seq),
        msa_to_consensus_pos
    )

    stats = ClusterStats(
        sample_name=cluster_info.sample_name,
        cluster_num=cluster_info.cluster_num,
        ric=cluster_info.ric,
        size=cluster_info.size,
        actual_cluster_size=total_reads_in_file,
        num_reads_analyzed=len(alignments),
        mean_edit_distance=np.mean(edit_distances),
        median_edit_distance=np.median(edit_distances),
        p95_edit_distance=np.percentile(edit_distances, 95),
        mean_error_rate=np.mean(error_rates),
        median_error_rate=np.median(error_rates),
        insertion_rate=total_insertions / total_errors if total_errors > 0 else 0.0,
        deletion_rate=total_deletions / total_errors if total_errors > 0 else 0.0,
        substitution_rate=total_substitutions / total_errors if total_errors > 0 else 0.0,
        consensus_length=len(consensus_seq),
        mean_read_length=np.mean(read_lengths),
        has_contiguous_errors=has_contiguous,
        rid=cluster_info.rid,
        rid_min=cluster_info.rid_min
    )

    return stats, alignments


def write_cluster_statistics(cluster_stats: List[ClusterStats], output_file: str):
    """Write per-cluster statistics to CSV file."""
    fieldnames = [
        'sample_name', 'cluster_num', 'ric', 'size', 'actual_cluster_size', 'num_reads_analyzed',
        'mean_edit_distance', 'median_edit_distance', 'p95_edit_distance',
        'mean_error_rate', 'median_error_rate',
        'insertion_rate', 'deletion_rate', 'substitution_rate',
        'consensus_length', 'mean_read_length', 'has_contiguous_errors',
        'rid', 'rid_min'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for stats in cluster_stats:
            writer.writerow(stats._asdict())

    logging.info(f"Wrote cluster statistics to: {output_file}")


def analyze_read_depth_convergence(cluster_stats: List[ClusterStats]) -> Dict:
    """
    Analyze how consensus quality converges with increasing cluster size.

    Uses power law decay curve fitting with fixed plateau:
    error_rate = plateau + A * cluster_size^(-alpha)

    The plateau is calculated as a weighted mean from high-coverage clusters.
    This approach is theoretically grounded: consensus error from majority voting
    typically decreases as a power law with read depth.

    Returns:
        Dict with convergence statistics and thresholds
    """
    # Extract cluster size vs error rate data
    cluster_sizes = np.array([s.actual_cluster_size for s in cluster_stats])
    error_rates = np.array([s.mean_error_rate for s in cluster_stats])

    # Calculate fixed plateau from high-coverage clusters (top 20% by cluster size)
    # Use weighted mean: larger clusters give more reliable error rate estimates
    high_depth_threshold = np.percentile(cluster_sizes, 80)
    high_depth_mask = cluster_sizes >= high_depth_threshold

    if np.sum(high_depth_mask) > 0:
        high_depth_sizes = cluster_sizes[high_depth_mask]
        high_depth_errors = error_rates[high_depth_mask]

        # Weighted mean: weight by cluster size (more reads = more reliable)
        plateau = np.sum(high_depth_errors * high_depth_sizes) / np.sum(high_depth_sizes)
        plateau_stderr = np.sqrt(np.sum(high_depth_sizes * (high_depth_errors - plateau)**2) / np.sum(high_depth_sizes))
        empirical_baseline_n = np.sum(high_depth_mask)
    else:
        # Fallback if insufficient data
        plateau = np.percentile(error_rates, 5)
        plateau_stderr = np.std(error_rates)
        empirical_baseline_n = 0

    # Define power law decay model with fixed plateau
    # error_rate = plateau + A * cluster_size^(-alpha)
    def power_law_decay(cluster_size, amplitude, alpha):
        return plateau + amplitude * np.power(cluster_size, -alpha)

    # Initial parameter guesses
    initial_error_est = np.percentile(error_rates, 95)
    initial_guess = [
        initial_error_est - plateau,  # amplitude A
        0.5  # alpha (0.5 would match sqrt scaling from statistics theory)
    ]

    # Bounds: amplitude > 0, alpha > 0
    bounds = (
        [0, 0],  # lower bounds
        [1, 5]  # upper bounds (alpha rarely > 5 in practice)
    )

    try:
        # Fit the curve
        params, covariance = curve_fit(
            power_law_decay,
            cluster_sizes,
            error_rates,
            p0=initial_guess,
            bounds=bounds,
            maxfev=10000
        )

        amplitude, alpha = params

        # Calculate standard errors from covariance matrix
        param_errors = np.sqrt(np.diag(covariance))
        amplitude_stderr = param_errors[0]
        alpha_stderr = param_errors[1]

        # Calculate convergence thresholds
        # 95% convergence: plateau + 0.05*A = plateau + A * cluster_size^(-alpha)
        # Solve: 0.05*A = A * cluster_size^(-alpha)
        # => 0.05 = cluster_size^(-alpha)
        # => cluster_size = (1/0.05)^(1/alpha) = 20^(1/alpha)

        if alpha > 0:
            cluster_size_90pct = np.power(1.0 / 0.10, 1.0 / alpha)  # 10^(1/alpha)
            cluster_size_95pct = np.power(1.0 / 0.05, 1.0 / alpha)  # 20^(1/alpha)
            cluster_size_99pct = np.power(1.0 / 0.01, 1.0 / alpha)  # 100^(1/alpha)
        else:
            cluster_size_90pct = cluster_size_95pct = cluster_size_99pct = np.inf

        curve_fit_success = True

    except Exception as e:
        logging.warning(f"Curve fitting failed: {e}. Using plateau estimate only.")
        amplitude = np.percentile(error_rates, 95) - plateau
        alpha = np.nan
        amplitude_stderr = np.nan
        alpha_stderr = np.nan
        cluster_size_90pct = cluster_size_95pct = cluster_size_99pct = np.nan
        curve_fit_success = False

    # Calculate improvement potential
    # How much does error rate improve from low to high depth?
    low_depth_threshold = np.percentile(cluster_sizes, 20)
    low_depth_mask = cluster_sizes <= low_depth_threshold

    if np.sum(low_depth_mask) > 0:
        low_depth_mean_error = np.mean(error_rates[low_depth_mask])
        improvement_from_depth = (low_depth_mean_error - plateau) / low_depth_mean_error
    else:
        low_depth_mean_error = np.nan
        improvement_from_depth = np.nan

    return {
        'curve_fit_success': curve_fit_success,
        'plateau_error_rate': plateau,
        'plateau_stderr': plateau_stderr,
        'amplitude': amplitude if curve_fit_success else np.nan,
        'amplitude_stderr': amplitude_stderr if curve_fit_success else np.nan,
        'alpha': alpha if curve_fit_success else np.nan,
        'alpha_stderr': alpha_stderr if curve_fit_success else np.nan,
        'convergence_cluster_size_90pct': cluster_size_90pct,
        'convergence_cluster_size_95pct': cluster_size_95pct,
        'convergence_cluster_size_99pct': cluster_size_99pct,
        'empirical_baseline': plateau,  # Now same as fitted plateau
        'empirical_baseline_std': plateau_stderr,
        'empirical_baseline_depth_threshold': high_depth_threshold,
        'empirical_baseline_n': empirical_baseline_n,
        'low_depth_mean_error': low_depth_mean_error,
        'improvement_from_depth': improvement_from_depth
    }


def write_dataset_summary_with_convergence(cluster_stats: List[ClusterStats], convergence: Dict, output_file: str):
    """Write dataset-wide summary statistics to text file."""
    if not cluster_stats:
        return

    # Aggregate statistics
    total_clusters = len(cluster_stats)
    total_reads = sum(s.num_reads_analyzed for s in cluster_stats)

    all_error_rates = [s.mean_error_rate for s in cluster_stats]
    all_insertions = [s.insertion_rate for s in cluster_stats]
    all_deletions = [s.deletion_rate for s in cluster_stats]
    all_substitutions = [s.substitution_rate for s in cluster_stats]

    clusters_with_contiguous = sum(1 for s in cluster_stats if s.has_contiguous_errors)

    with open(output_file, 'w') as f:
        f.write("=== Speconsense Empirical Error Analysis ===\n\n")

        f.write(f"Dataset Overview:\n")
        f.write(f"  Total clusters analyzed: {total_clusters}\n")
        f.write(f"  Total reads analyzed: {total_reads}\n")
        f.write(f"  Mean reads per cluster: {total_reads / total_clusters:.1f}\n\n")

        f.write(f"Error Rate Estimates:\n")
        f.write(f"  Mean error rate (all clusters): {np.mean(all_error_rates) * 100:.2f}%\n")
        f.write(f"  Median error rate: {np.median(all_error_rates) * 100:.2f}%\n")
        f.write(f"  95th percentile: {np.percentile(all_error_rates, 95) * 100:.2f}%\n")
        f.write(f"  5th percentile: {np.percentile(all_error_rates, 5) * 100:.2f}%\n\n")

        # Report true sequencing error rate from convergence analysis
        f.write(f"True Sequencing Error Rate (power law model):\n")
        f.write(f"  Plateau (weighted mean, top 20%): {convergence['plateau_error_rate'] * 100:.2f}% ± {convergence['plateau_stderr'] * 100:.2f}%\n")
        f.write(f"  Calculated from {convergence['empirical_baseline_n']} clusters with size ≥ {convergence['empirical_baseline_depth_threshold']:.0f}\n")

        if convergence['curve_fit_success'] and not np.isnan(convergence['alpha']):
            f.write(f"\n  Power law fit: error = plateau + A × (cluster_size)^(-α)\n")
            f.write(f"    Amplitude A: {convergence['amplitude'] * 100:.2f}%\n")
            f.write(f"    Exponent α: {convergence['alpha']:.2f}\n")
            if convergence['alpha'] < 0.3:
                f.write(f"    (α < 0.5: slower than sqrt scaling, harder to converge)\n")
            elif convergence['alpha'] > 0.7:
                f.write(f"    (α > 0.5: faster than sqrt scaling, converges easily)\n")
            else:
                f.write(f"    (α ≈ 0.5: matches theoretical sqrt scaling from statistics)\n")
        f.write(f"\n")

        f.write(f"Error Type Distribution:\n")
        f.write(f"  Mean insertion rate: {np.mean(all_insertions) * 100:.1f}%\n")
        f.write(f"  Mean deletion rate: {np.mean(all_deletions) * 100:.1f}%\n")
        f.write(f"  Mean substitution rate: {np.mean(all_substitutions) * 100:.1f}%\n\n")

        f.write(f"Positional Error Patterns:\n")
        f.write(f"  Clusters with contiguous high-error regions: {clusters_with_contiguous} ")
        f.write(f"({clusters_with_contiguous / total_clusters * 100:.1f}%)\n\n")

        # Write convergence analysis
        f.write(f"Consensus Quality Convergence:\n")
        if convergence['curve_fit_success']:
            if not np.isnan(convergence['convergence_cluster_size_95pct']):
                f.write(f"  Cluster size for 95% convergence: {convergence['convergence_cluster_size_95pct']:.0f} reads\n")
            if not np.isnan(convergence['convergence_cluster_size_90pct']):
                f.write(f"  Cluster size for 90% convergence: {convergence['convergence_cluster_size_90pct']:.0f} reads\n")

            # Provide recommendations based on thresholds
            size_95 = convergence['convergence_cluster_size_95pct']
            size_90 = convergence['convergence_cluster_size_90pct']

            if not np.isnan(size_95) and size_95 < 1000:
                f.write(f"\n  Recommendations:\n")
                f.write(f"    - Cluster size ≥ {size_95:.0f}: Recommended minimum (95% of optimal)\n")
                if not np.isnan(size_90) and size_90 < size_95:
                    f.write(f"    - Cluster size ≥ {size_90:.0f}: Acceptable (90% of optimal)\n")
                f.write(f"    - Diminishing returns beyond cluster size ≥ {size_95 * 2:.0f}\n")
        else:
            f.write(f"  Curve fitting failed - using empirical estimates\n")

        if not np.isnan(convergence['improvement_from_depth']):
            f.write(f"\n  Error rate improvement from consensus depth:\n")
            f.write(f"    Low depth (bottom 20%): {convergence['low_depth_mean_error'] * 100:.2f}%\n")
            f.write(f"    High depth (plateau): {convergence['plateau_error_rate'] * 100:.2f}%\n")
            f.write(f"    Improvement: {convergence['improvement_from_depth'] * 100:.1f}%\n")

        f.write(f"\n")
        f.write(f"Comparison to Synthetic Testing Expectations:\n")
        f.write(f"  Modern chemistry (R10.4.1): 1-2% error expected\n")
        f.write(f"  True error rate (plateau): {convergence['plateau_error_rate'] * 100:.2f}%\n")
        f.write(f"  Overall mean (all cluster sizes): {np.mean(all_error_rates) * 100:.2f}%\n")

        plateau = convergence['plateau_error_rate']
        if plateau < 0.02:
            f.write(f"  ✓ Consistent with modern chemistry\n")
        elif plateau < 0.05:
            f.write(f"  ✓ Good quality, slightly above modern expectations\n")
        else:
            f.write(f"  ⚠ Higher than expected, may indicate older chemistry or quality issues\n")

    logging.info(f"Wrote dataset summary to: {output_file}")


# CLI functions moved to refine.py
# This module now contains only analysis functions that can be imported
