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
import argparse
import logging
from typing import List, Dict, Tuple, Optional, NamedTuple
from collections import defaultdict, Counter

import edlib
from Bio import SeqIO
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import binomtest
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from tqdm import tqdm


class ClusterInfo(NamedTuple):
    """Information about a cluster and its consensus."""
    sample_name: str
    cluster_num: int
    consensus_seq: str
    reads_file: str  # Path to FASTQ file with raw reads
    ric: int  # Reads in consensus (from header)
    size: int  # Total cluster size (from header)
    p50_diff: Optional[float] = None
    p95_diff: Optional[float] = None


class ErrorPosition(NamedTuple):
    """An error at a specific position in the consensus."""
    position: int  # 0-indexed position in consensus
    error_type: str  # 'sub', 'ins', or 'del'


class ReadAlignment(NamedTuple):
    """Alignment result for a single read against consensus."""
    read_id: str
    read_length: int
    edit_distance: int
    num_insertions: int
    num_deletions: int
    num_substitutions: int
    error_positions: List[ErrorPosition]  # Detailed error information
    alignment_cigar: str


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
    p50_diff: Optional[float] = None
    p95_diff: Optional[float] = None


def parse_consensus_header(header: str) -> Tuple[str, int, int, Optional[float], Optional[float]]:
    """
    Parse a consensus FASTA header to extract metadata.

    Header format: >sample-c1 size=200 ric=100 [p50diff=0.0] [p95diff=1.0]

    Returns:
        Tuple of (sample_name, size, ric, p50_diff, p95_diff)
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
    p50_diff = p95_diff = None

    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            if key == 'size':
                size = int(value)
            elif key == 'ric':
                ric = int(value)
            elif key == 'p50diff':
                p50_diff = float(value)
            elif key == 'p95diff':
                p95_diff = float(value)

    if size is None or ric is None:
        raise ValueError(f"Missing size or ric in header: {header}")

    return sample_name, cluster_num, size, ric, p50_diff, p95_diff


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
                sample_name, cluster_num, size, ric, p50_diff, p95_diff = parse_consensus_header(record.description)

                # Extract base sample name (remove -c{num} suffix)
                base_name = re.sub(r'-c\d+$', '', sample_name)

                consensus_map[(base_name, cluster_num)] = ClusterInfo(
                    sample_name=sample_name,
                    cluster_num=cluster_num,
                    consensus_seq=str(record.seq),
                    reads_file='',  # Will be filled in later
                    ric=ric,
                    size=size,
                    p50_diff=p50_diff,
                    p95_diff=p95_diff
                )
        except Exception as e:
            logging.warning(f"Error parsing {fasta_file}: {e}")
            continue

    return consensus_map


def align_read_to_consensus(read_seq: str, consensus_seq: str, read_id: str) -> ReadAlignment:
    """
    Align a read to its consensus sequence and extract error characteristics.

    Args:
        read_seq: Read sequence
        consensus_seq: Consensus sequence
        read_id: Read identifier

    Returns:
        ReadAlignment with error statistics
    """
    # Perform alignment with edlib
    result = edlib.align(
        read_seq,
        consensus_seq,
        mode='NW',  # Global alignment
        task='path'  # Get full alignment path
    )

    if result['editDistance'] == -1:
        # Alignment failed
        logging.warning(f"Alignment failed for read {read_id}")
        return ReadAlignment(
            read_id=read_id,
            read_length=len(read_seq),
            edit_distance=-1,
            num_insertions=0,
            num_deletions=0,
            num_substitutions=0,
            error_positions=[],
            alignment_cigar=''
        )

    # Parse CIGAR string to count error types
    cigar = result.get('cigar', '')
    edit_distance = result['editDistance']

    # Count error types from CIGAR
    # CIGAR format: M (match/mismatch), I (insertion to ref), D (deletion from ref), X (mismatch), = (match)
    insertions = cigar.count('I') if cigar else 0
    deletions = cigar.count('D') if cigar else 0

    # For substitutions, we need to count X operations (or infer from distance - indels)
    # edlib uses: = (match), X (mismatch), I (insertion), D (deletion)
    substitutions = cigar.count('X') if 'X' in cigar else 0

    # If using simple CIGAR (only I/D/M), estimate substitutions
    if substitutions == 0 and cigar.count('M') > 0:
        # Estimate: total distance - insertions - deletions
        substitutions = max(0, edit_distance - insertions - deletions)

    # Extract error positions from alignment
    error_positions = extract_error_positions(result, consensus_seq)

    return ReadAlignment(
        read_id=read_id,
        read_length=len(read_seq),
        edit_distance=edit_distance,
        num_insertions=insertions,
        num_deletions=deletions,
        num_substitutions=substitutions,
        error_positions=error_positions,
        alignment_cigar=cigar
    )


def extract_error_positions(alignment_result: dict, consensus_seq: str) -> List[ErrorPosition]:
    """
    Extract positions in the consensus where errors occur, with error types.

    Args:
        alignment_result: edlib alignment result dictionary
        consensus_seq: Consensus sequence

    Returns:
        List of ErrorPosition objects with position and error type
    """
    cigar = alignment_result.get('cigar', '')
    if not cigar:
        return []

    error_positions = []
    consensus_pos = 0

    # Parse CIGAR operations
    i = 0
    while i < len(cigar):
        # Count consecutive digits
        num_str = ''
        while i < len(cigar) and cigar[i].isdigit():
            num_str += cigar[i]
            i += 1

        if i >= len(cigar):
            break

        count = int(num_str) if num_str else 1
        op = cigar[i]

        if op == '=' or op == 'M':
            # Match - advance position
            if op == 'M':
                # M could be match or mismatch; assume some are errors
                # For more precise analysis, would need to compare sequences directly
                pass
            consensus_pos += count
        elif op == 'X':
            # Substitution (mismatch)
            for j in range(count):
                error_positions.append(ErrorPosition(consensus_pos, 'sub'))
                consensus_pos += 1
        elif op == 'D':
            # Deletion from reference (gap in read)
            for j in range(count):
                error_positions.append(ErrorPosition(consensus_pos, 'del'))
                consensus_pos += 1
        elif op == 'I':
            # Insertion to reference (gap in consensus)
            # Attribute to current consensus position
            error_positions.append(ErrorPosition(consensus_pos, 'ins'))

        i += 1

    return error_positions


def detect_contiguous_error_regions(error_positions: List[ErrorPosition], consensus_length: int, window_size: int = 50) -> bool:
    """
    Detect if errors cluster into contiguous high-error regions.

    This helps identify read-end degradation and systematic error patterns.

    Args:
        error_positions: List of ErrorPosition objects
        consensus_length: Length of consensus sequence
        window_size: Size of sliding window for error density calculation

    Returns:
        True if contiguous high-error region detected
    """
    if len(error_positions) < 5:
        return False

    # Extract positions from ErrorPosition objects
    positions = [ep.position for ep in error_positions]

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


def detect_homopolymer_run(seq: str, position: int, min_length: int = 4) -> bool:
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


def calculate_local_gc(seq: str, position: int, window: int = 10) -> float:
    """
    Calculate GC content in a window around a position.

    Args:
        seq: DNA sequence
        position: Center position (0-indexed)
        window: Window size (±window/2 around position)

    Returns:
        GC content as a fraction (0.0 to 1.0)
    """
    half_window = window // 2
    start = max(0, position - half_window)
    end = min(len(seq), position + half_window + 1)

    if start >= end:
        return 0.0

    window_seq = seq[start:end]
    if len(window_seq) == 0:
        return 0.0

    gc_count = window_seq.count('G') + window_seq.count('C')
    return gc_count / len(window_seq)


def benjamini_hochberg_correction(p_values: List[float], alpha: float = 0.05) -> List[float]:
    """
    Apply Benjamini-Hochberg FDR correction to p-values.

    Args:
        p_values: List of p-values
        alpha: Significance level

    Returns:
        List of adjusted p-values (q-values)
    """
    n = len(p_values)
    if n == 0:
        return []

    # Create list of (index, p_value) and sort by p_value
    indexed_p = [(i, p) for i, p in enumerate(p_values)]
    indexed_p.sort(key=lambda x: x[1])

    # Calculate adjusted p-values
    adjusted = [0.0] * n
    for rank, (original_idx, p_value) in enumerate(indexed_p, 1):
        # BH adjustment: p * n / rank
        adjusted[original_idx] = min(1.0, p_value * n / rank)

    # Enforce monotonicity (optional but common)
    sorted_adjusted = sorted(enumerate(adjusted), key=lambda x: p_values[x[0]])
    for i in range(len(sorted_adjusted) - 2, -1, -1):
        idx_i = sorted_adjusted[i][0]
        idx_next = sorted_adjusted[i + 1][0]
        adjusted[idx_i] = min(adjusted[idx_i], adjusted[idx_next])

    return adjusted


class PositionStats(NamedTuple):
    """Statistics for a single position in the consensus."""
    position: int
    coverage: int
    error_count: int
    error_rate: float
    sub_count: int
    ins_count: int
    del_count: int
    p_value: float
    is_homopolymer: bool
    gc_content: float
    nucleotide: str


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


def analyze_positional_variation(
    alignments: List[ReadAlignment],
    consensus_seq: str,
    overall_error_rate: float
) -> List[PositionStats]:
    """
    Analyze error rates at each position in the consensus sequence.

    Identifies positions with significantly elevated error rates, which may
    indicate homopolymer regions, structural variants, or other error-prone
    sequence contexts.

    Args:
        alignments: List of read alignments
        consensus_seq: Consensus sequence
        overall_error_rate: Overall error rate for binomial test

    Returns:
        List of PositionStats for each consensus position
    """
    consensus_len = len(consensus_seq)

    # Build error frequency matrix
    # For each position: [sub_count, ins_count, del_count, total_coverage]
    error_matrix = np.zeros((consensus_len, 4), dtype=int)

    for alignment in alignments:
        # Count this read as coverage for all positions
        # (assuming full-length alignments; could be refined with alignment bounds)
        for pos in range(consensus_len):
            error_matrix[pos, 3] += 1  # coverage

        # Add errors at specific positions
        for error_pos in alignment.error_positions:
            pos = error_pos.position
            if 0 <= pos < consensus_len:
                if error_pos.error_type == 'sub':
                    error_matrix[pos, 0] += 1
                elif error_pos.error_type == 'ins':
                    error_matrix[pos, 1] += 1
                elif error_pos.error_type == 'del':
                    error_matrix[pos, 2] += 1

    # Calculate statistics for each position
    position_stats = []
    p_values = []

    for pos in range(consensus_len):
        sub_count = error_matrix[pos, 0]
        ins_count = error_matrix[pos, 1]
        del_count = error_matrix[pos, 2]
        coverage = error_matrix[pos, 3]

        error_count = sub_count + ins_count + del_count
        error_rate = error_count / coverage if coverage > 0 else 0.0

        # Binomial test: is this position's error rate significantly different from overall?
        # Using one-tailed test since we're interested in high error rates
        if coverage > 0:
            result = binomtest(error_count, coverage, overall_error_rate, alternative='greater')
            p_value = result.pvalue
        else:
            p_value = 1.0

        p_values.append(p_value)

        # Sequence context
        is_homopolymer = detect_homopolymer_run(consensus_seq, pos, min_length=4)
        gc_content = calculate_local_gc(consensus_seq, pos, window=10)
        nucleotide = consensus_seq[pos] if pos < len(consensus_seq) else 'N'

        position_stats.append(PositionStats(
            position=pos,
            coverage=coverage,
            error_count=error_count,
            error_rate=error_rate,
            sub_count=sub_count,
            ins_count=ins_count,
            del_count=del_count,
            p_value=p_value,
            is_homopolymer=is_homopolymer,
            gc_content=gc_content,
            nucleotide=nucleotide
        ))

    # Apply FDR correction
    adjusted_p_values = benjamini_hochberg_correction(p_values)

    # Update position_stats with adjusted p-values
    # (we'll add a q_value field later if needed, for now just use p_value)
    # For this implementation, we'll create a new list with adjusted p-values
    position_stats_with_fdr = []
    for i, stats in enumerate(position_stats):
        # Create a new PositionStats with q_value instead of p_value
        # Actually, PositionStats doesn't have q_value field, so we'll keep p_value
        # and interpret it as q_value after correction
        position_stats_with_fdr.append(PositionStats(
            position=stats.position,
            coverage=stats.coverage,
            error_count=stats.error_count,
            error_rate=stats.error_rate,
            sub_count=stats.sub_count,
            ins_count=stats.ins_count,
            del_count=stats.del_count,
            p_value=adjusted_p_values[i],  # Use adjusted p-value (q-value)
            is_homopolymer=stats.is_homopolymer,
            gc_content=stats.gc_content,
            nucleotide=stats.nucleotide
        ))

    return position_stats_with_fdr


class WindowStats(NamedTuple):
    """Statistics for a window in the consensus sequence."""
    window_id: int
    start: int
    end: int
    length: int
    error_rate: float
    error_count: int
    coverage: int
    sub_count: int
    ins_count: int
    del_count: int
    gc_content: float
    has_homopolymer: bool
    sequence: str


def analyze_windowed_consensus(
    alignments: List[ReadAlignment],
    consensus_seq: str,
    window_size: int = 50
) -> List[WindowStats]:
    """
    Analyze consensus sequence in fixed non-overlapping windows.

    Useful for identifying regional error patterns like read-end degradation.

    Args:
        alignments: List of read alignments
        consensus_seq: Consensus sequence
        window_size: Size of each window in bp

    Returns:
        List of WindowStats for each window
    """
    consensus_len = len(consensus_seq)
    num_windows = (consensus_len + window_size - 1) // window_size

    window_stats = []

    for w in range(num_windows):
        start = w * window_size
        end = min(start + window_size, consensus_len)
        window_len = end - start

        # Collect errors in this window from all reads
        window_errors = {'sub': 0, 'ins': 0, 'del': 0}
        total_coverage = len(alignments)  # Each read contributes to coverage

        for alignment in alignments:
            # Count errors in window [start, end)
            for error_pos in alignment.error_positions:
                pos = error_pos.position
                if start <= pos < end:
                    window_errors[error_pos.error_type] += 1

        # Calculate statistics
        total_errors = sum(window_errors.values())
        # Error rate per base per read
        error_rate = total_errors / (total_coverage * window_len) if (total_coverage > 0 and window_len > 0) else 0.0

        # Sequence context
        window_seq = consensus_seq[start:end]
        gc_count = window_seq.count('G') + window_seq.count('C')
        gc_content = gc_count / len(window_seq) if len(window_seq) > 0 else 0.0

        # Check for homopolymer runs in window
        has_homopolymer = False
        for pos in range(start, end):
            if detect_homopolymer_run(consensus_seq, pos, min_length=4):
                has_homopolymer = True
                break

        window_stats.append(WindowStats(
            window_id=w,
            start=start,
            end=end,
            length=window_len,
            error_rate=error_rate,
            error_count=total_errors,
            coverage=total_coverage,
            sub_count=window_errors['sub'],
            ins_count=window_errors['ins'],
            del_count=window_errors['del'],
            gc_content=gc_content,
            has_homopolymer=has_homopolymer,
            sequence=window_seq
        ))

    return window_stats


def analyze_all_cluster_positions(
    cluster_alignments: List[Tuple[ClusterInfo, List[ReadAlignment]]]
) -> List[ClusterPositionStats]:
    """
    Analyze positional variation across all clusters.

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

        # Build error frequency matrix for this cluster
        error_matrix = np.zeros((consensus_len, 4), dtype=int)

        for alignment in alignments:
            # Count coverage
            for pos in range(consensus_len):
                error_matrix[pos, 3] += 1

            # Add errors
            for error_pos in alignment.error_positions:
                pos = error_pos.position
                if 0 <= pos < consensus_len:
                    if error_pos.error_type == 'sub':
                        error_matrix[pos, 0] += 1
                    elif error_pos.error_type == 'ins':
                        error_matrix[pos, 1] += 1
                    elif error_pos.error_type == 'del':
                        error_matrix[pos, 2] += 1

        # Calculate statistics for each position in this cluster
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
    Analyze all reads in a cluster against their consensus.

    Args:
        cluster_info: Cluster information including reads file and consensus

    Returns:
        Tuple of (cluster statistics, list of read alignments)
    """
    if not os.path.exists(cluster_info.reads_file):
        logging.warning(f"Reads file not found: {cluster_info.reads_file}")
        return None, []

    consensus_seq = cluster_info.consensus_seq
    alignments = []
    total_reads_in_file = 0

    # Read and align all reads
    for record in SeqIO.parse(cluster_info.reads_file, 'fastq'):
        total_reads_in_file += 1
        alignment = align_read_to_consensus(
            str(record.seq),
            consensus_seq,
            record.id
        )
        if alignment.edit_distance >= 0:  # Skip failed alignments
            alignments.append(alignment)

    if not alignments:
        logging.warning(f"No successful alignments for {cluster_info.sample_name}")
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
        len(consensus_seq)
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
        p50_diff=cluster_info.p50_diff,
        p95_diff=cluster_info.p95_diff
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
        'p50_diff', 'p95_diff'
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


def plot_positional_variation(
    position_stats: List[PositionStats],
    consensus_seq: str,
    output_file: str,
    significance_threshold: float = 0.05
):
    """
    Create visualization of position-wise error rates in the consensus.

    Shows:
    - Error rate by position with error type breakdown
    - Homopolymer region highlighting
    - Statistically significant positions marked
    - Sequence context annotation

    Args:
        position_stats: List of position statistics
        consensus_seq: Consensus sequence
        output_file: Path to save the plot PNG
        significance_threshold: FDR threshold for marking significant positions
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), height_ratios=[3, 1])

    positions = [ps.position for ps in position_stats]
    error_rates = [ps.error_rate * 100 for ps in position_stats]  # Convert to percentage

    # Extract error type counts
    sub_rates = [(ps.sub_count / ps.coverage * 100) if ps.coverage > 0 else 0 for ps in position_stats]
    ins_rates = [(ps.ins_count / ps.coverage * 100) if ps.coverage > 0 else 0 for ps in position_stats]
    del_rates = [(ps.del_count / ps.coverage * 100) if ps.coverage > 0 else 0 for ps in position_stats]

    # Top panel: Stacked error rates by type
    ax1.fill_between(positions, 0, sub_rates, alpha=0.6, color='#e74c3c', label='Substitutions', step='mid')

    # Stack insertions on top of substitutions
    ins_bottom = sub_rates
    ins_top = [sub_rates[i] + ins_rates[i] for i in range(len(positions))]
    ax1.fill_between(positions, ins_bottom, ins_top, alpha=0.6, color='#3498db', label='Insertions', step='mid')

    # Stack deletions on top
    del_bottom = ins_top
    del_top = [ins_top[i] + del_rates[i] for i in range(len(positions))]
    ax1.fill_between(positions, del_bottom, del_top, alpha=0.6, color='#2ecc71', label='Deletions', step='mid')

    # Mark significant positions with red dots
    significant_positions = [ps.position for ps in position_stats if ps.p_value < significance_threshold]
    significant_error_rates = [ps.error_rate * 100 for ps in position_stats if ps.p_value < significance_threshold]

    if significant_positions:
        ax1.scatter(significant_positions, significant_error_rates,
                   color='red', s=30, marker='o', alpha=0.7,
                   label=f'Significant (FDR < {significance_threshold})', zorder=5)

    # Highlight homopolymer regions with vertical shading
    homopolymer_regions = []
    in_homopolymer = False
    start_pos = 0

    for ps in position_stats:
        if ps.is_homopolymer and not in_homopolymer:
            start_pos = ps.position
            in_homopolymer = True
        elif not ps.is_homopolymer and in_homopolymer:
            homopolymer_regions.append((start_pos, ps.position - 1))
            in_homopolymer = False

    if in_homopolymer:
        homopolymer_regions.append((start_pos, len(position_stats) - 1))

    for start, end in homopolymer_regions:
        ax1.axvspan(start, end, color='yellow', alpha=0.15, zorder=0)

    ax1.set_ylabel('Error Rate (%)', fontsize=11)
    ax1.set_title('Position-wise Error Rates with Error Type Breakdown', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, len(positions))

    # Bottom panel: Sequence context (homopolymer and GC content)
    gc_contents = [ps.gc_content * 100 for ps in position_stats]
    homopolymer_indicator = [1 if ps.is_homopolymer else 0 for ps in position_stats]

    ax2.plot(positions, gc_contents, color='purple', linewidth=1, alpha=0.7, label='Local GC%')
    ax2.fill_between(positions, 0, homopolymer_indicator, color='orange', alpha=0.4,
                     step='mid', label='Homopolymer (≥4bp)')

    ax2.set_xlabel('Position in Consensus Sequence', fontsize=11)
    ax2.set_ylabel('GC% / Homopoly', fontsize=10)
    ax2.set_ylim(-0.1, max(max(gc_contents) * 1.1, 1.2))
    ax2.set_xlim(0, len(positions))
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def plot_windowed_heatmap(
    window_stats: List[WindowStats],
    output_file: str
):
    """
    Create visualization of windowed error analysis.

    Shows:
    - Error rate by window as stacked bar chart
    - Error type breakdown
    - Window positions and sequence context

    Args:
        window_stats: List of window statistics
        output_file: Path to save the plot PNG
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), height_ratios=[3, 1])

    window_ids = [ws.window_id for ws in window_stats]
    window_positions = [ws.start for ws in window_stats]

    # Calculate error rates by type (percentage)
    sub_rates = [(ws.sub_count / (ws.coverage * ws.length) * 100) if (ws.coverage > 0 and ws.length > 0) else 0 for ws in window_stats]
    ins_rates = [(ws.ins_count / (ws.coverage * ws.length) * 100) if (ws.coverage > 0 and ws.length > 0) else 0 for ws in window_stats]
    del_rates = [(ws.del_count / (ws.coverage * ws.length) * 100) if (ws.coverage > 0 and ws.length > 0) else 0 for ws in window_stats]

    # Create stacked bar chart
    bar_width = 0.8
    ax1.bar(window_ids, sub_rates, bar_width, label='Substitutions', color='#e74c3c', alpha=0.8)
    ax1.bar(window_ids, ins_rates, bar_width, bottom=sub_rates, label='Insertions', color='#3498db', alpha=0.8)

    # Stack deletions on top
    bottom_for_del = [sub_rates[i] + ins_rates[i] for i in range(len(window_ids))]
    ax1.bar(window_ids, del_rates, bar_width, bottom=bottom_for_del, label='Deletions', color='#2ecc71', alpha=0.8)

    # Highlight windows with homopolymers
    homopoly_windows = [ws.window_id for ws in window_stats if ws.has_homopolymer]
    if homopoly_windows:
        for wid in homopoly_windows:
            ax1.axvline(wid, color='orange', linestyle='--', alpha=0.3, linewidth=1)

    ax1.set_xlabel('Window ID', fontsize=11)
    ax1.set_ylabel('Error Rate (%)', fontsize=11)
    ax1.set_title(f'Windowed Error Analysis (Window Size: {window_stats[0].length if window_stats else 50}bp)',
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')

    # Bottom panel: GC content and homopolymer indicator
    gc_contents = [ws.gc_content * 100 for ws in window_stats]
    homopoly_indicator = [1 if ws.has_homopolymer else 0 for ws in window_stats]

    ax2.bar(window_ids, gc_contents, bar_width, color='purple', alpha=0.6, label='GC%')
    ax2.scatter([ws.window_id for ws in window_stats if ws.has_homopolymer],
               [50] * len(homopoly_windows),
               color='orange', s=100, marker='v', alpha=0.8, label='Has Homopolymer', zorder=5)

    ax2.set_xlabel('Window ID (corresponds to position)', fontsize=11)
    ax2.set_ylabel('GC Content (%)', fontsize=10)
    ax2.set_ylim(0, 100)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')

    # Add position labels to help interpret window IDs
    if len(window_ids) <= 20:
        # Show all position labels
        ax2.set_xticks(window_ids)
        ax2.set_xticklabels([f'{ws.start}-{ws.end}' for ws in window_stats], rotation=45, ha='right', fontsize=8)
    else:
        # Show subset of labels
        step = max(1, len(window_ids) // 10)
        tick_indices = list(range(0, len(window_ids), step))
        ax2.set_xticks([window_ids[i] for i in tick_indices])
        ax2.set_xticklabels([f'{window_stats[i].start}' for i in tick_indices], fontsize=9)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def plot_positional_error_histogram(
    position_stats: List[ClusterPositionStats],
    output_file: str
):
    """
    Create histogram of positional error rates across all clusters.

    Shows the distribution of error rates at individual positions, which
    helps assess clustering quality (normal distribution is better) and
    identify outlier positions.

    Args:
        position_stats: List of ClusterPositionStats for all positions
        output_file: Path to save the plot PNG
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), height_ratios=[2, 1])

    # Extract error rates
    error_rates = np.array([ps.error_rate * 100 for ps in position_stats])  # Convert to percentage

    # Main histogram
    bins = np.linspace(0, min(20, np.percentile(error_rates, 99.5)), 50)  # Cap at 20% or 99.5th percentile
    counts, bin_edges, patches = ax1.hist(error_rates, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')

    # Add statistics
    mean_error = np.mean(error_rates)
    median_error = np.median(error_rates)
    std_error = np.std(error_rates)

    ax1.axvline(mean_error, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_error:.2f}%')
    ax1.axvline(median_error, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_error:.2f}%')

    ax1.set_xlabel('Error Rate (%)', fontsize=11)
    ax1.set_ylabel('Number of Positions', fontsize=11)
    ax1.set_title(f'Distribution of Positional Error Rates Across All Clusters (n={len(position_stats):,} positions)',
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')

    # Add text box with summary statistics
    textstr = f'Total positions: {len(position_stats):,}\n'
    textstr += f'Mean: {mean_error:.2f}%\n'
    textstr += f'Median: {median_error:.2f}%\n'
    textstr += f'Std Dev: {std_error:.2f}%\n'
    textstr += f'95th percentile: {np.percentile(error_rates, 95):.2f}%\n'
    textstr += f'99th percentile: {np.percentile(error_rates, 99):.2f}%'

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.98, 0.97, textstr, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    # Bottom panel: Log-scale histogram to show tail
    ax2.hist(error_rates, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')
    ax2.set_yscale('log')
    ax2.set_xlabel('Error Rate (%)', fontsize=11)
    ax2.set_ylabel('Count (log scale)', fontsize=11)
    ax2.set_title('Same Distribution (log scale) - Shows Outlier Tail', fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def plot_read_error_histogram(
    cluster_alignments: List[Tuple[ClusterInfo, List[ReadAlignment]]],
    output_file: str
):
    """
    Create histogram of read error rates across all clusters.

    Shows the distribution of error rates for individual reads, which
    helps identify outlier reads that may be misclassified or contaminants.

    Args:
        cluster_alignments: List of (cluster_info, alignments) tuples
        output_file: Path to save the plot PNG
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), height_ratios=[2, 1])

    # Extract error rates for all reads
    read_error_rates = []
    for cluster_info, alignments in cluster_alignments:
        for alignment in alignments:
            error_rate = alignment.edit_distance / alignment.read_length if alignment.read_length > 0 else 0.0
            read_error_rates.append(error_rate * 100)  # Convert to percentage

    error_rates = np.array(read_error_rates)

    # Main histogram
    bins = np.linspace(0, min(20, np.percentile(error_rates, 99.5)), 50)  # Cap at 20% or 99.5th percentile
    counts, bin_edges, patches = ax1.hist(error_rates, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')

    # Add statistics
    mean_error = np.mean(error_rates)
    median_error = np.median(error_rates)
    std_error = np.std(error_rates)

    ax1.axvline(mean_error, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_error:.2f}%')
    ax1.axvline(median_error, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_error:.2f}%')

    ax1.set_xlabel('Error Rate (%)', fontsize=11)
    ax1.set_ylabel('Number of Reads', fontsize=11)
    ax1.set_title(f'Distribution of Read Error Rates Across All Clusters (n={len(error_rates):,} reads)',
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')

    # Add text box with summary statistics
    textstr = f'Total reads: {len(error_rates):,}\n'
    textstr += f'Mean: {mean_error:.2f}%\n'
    textstr += f'Median: {median_error:.2f}%\n'
    textstr += f'Std Dev: {std_error:.2f}%\n'
    textstr += f'95th percentile: {np.percentile(error_rates, 95):.2f}%\n'
    textstr += f'99th percentile: {np.percentile(error_rates, 99):.2f}%'

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.98, 0.97, textstr, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    # Bottom panel: Log-scale histogram to show tail
    ax2.hist(error_rates, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')
    ax2.set_yscale('log')
    ax2.set_xlabel('Error Rate (%)', fontsize=11)
    ax2.set_ylabel('Count (log scale)', fontsize=11)
    ax2.set_title('Same Distribution (log scale) - Shows Outlier Tail', fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def plot_convergence_analysis(cluster_stats: List[ClusterStats], convergence: Dict, output_file: str):
    """
    Create visualization of consensus quality convergence vs cluster size.

    Shows:
    - Scatter plot of individual cluster error rates
    - Sliding window mean
    - Fitted exponential decay curve
    - Plateau reference line

    Args:
        cluster_stats: List of cluster statistics
        convergence: Convergence analysis results from analyze_read_depth_convergence
        output_file: Path to save the plot PNG
    """
    # Extract data
    cluster_sizes = np.array([s.actual_cluster_size for s in cluster_stats])
    error_rates = np.array([s.mean_error_rate for s in cluster_stats])

    # Sort by cluster size for plotting
    sort_idx = np.argsort(cluster_sizes)
    sizes_sorted = cluster_sizes[sort_idx]
    error_sorted = error_rates[sort_idx]

    # Calculate sliding window statistics with fixed window size
    max_size = np.max(cluster_sizes)
    min_size = np.min(cluster_sizes)

    # Fixed window size: ±30 reads (adjustable based on data range)
    # For very small ranges, scale down
    size_range = max_size - min_size
    half_window = min(30, size_range / 10)  # At least 10 windows across range

    # Create evaluation points
    eval_points = np.linspace(min_size, max_size, min(200, len(cluster_stats) // 5))

    window_means = []
    window_n = []

    for size_center in eval_points:
        # Find points in fixed-size window
        in_window = np.abs(cluster_sizes - size_center) <= half_window

        if np.sum(in_window) >= 3:  # Need at least 3 points for stats
            window_errors = error_rates[in_window]
            mean = np.mean(window_errors)
            window_means.append(mean)
            window_n.append(len(window_errors))
        else:
            window_means.append(np.nan)
            window_n.append(0)

    window_means = np.array(window_means)

    # Generate fitted curve if available
    if convergence['curve_fit_success']:
        size_curve = np.linspace(min_size, max_size, 500)
        plateau = convergence['plateau_error_rate']
        amplitude = convergence['amplitude']
        alpha = convergence['alpha']

        # Generate power law curve: error_rate = plateau + A * cluster_size^(-alpha)
        if not np.isnan(amplitude) and not np.isnan(alpha) and alpha > 0:
            fitted_curve = plateau + amplitude * np.power(size_curve, -alpha)
        else:
            fitted_curve = None
    else:
        size_curve = None
        fitted_curve = None

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot scatter of all points (semi-transparent)
    ax.scatter(cluster_sizes, error_rates * 100, alpha=0.15, s=20, c='gray', label='Individual clusters')

    # Plot sliding window mean
    valid_mask = ~np.isnan(window_means)
    if np.sum(valid_mask) > 0:
        ax.plot(eval_points[valid_mask], window_means[valid_mask] * 100,
                'b-', linewidth=2.5, label=f'Windowed mean (±{half_window:.0f} reads)', zorder=3)

    # Plot fitted curve
    if fitted_curve is not None:
        ax.plot(size_curve, fitted_curve * 100, 'r-', linewidth=2,
                label='Fitted power law decay', zorder=4, alpha=0.8)

    # Plot plateau line
    plateau = convergence['plateau_error_rate']
    ax.axhline(y=plateau * 100, color='green', linestyle='--', linewidth=1.5,
               label=f'Plateau (true error): {plateau*100:.2f}%', zorder=2)

    # Annotations
    ax.set_xlabel('Cluster Size (reads)', fontsize=12)
    ax.set_ylabel('Error Rate (%)', fontsize=12)
    ax.set_title('Consensus Quality Convergence Analysis', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3, zorder=0)

    # Add text box with key statistics
    textstr = f'Plateau (weighted mean): {convergence["plateau_error_rate"]*100:.2f}%'
    if convergence['curve_fit_success'] and not np.isnan(convergence['alpha']):
        textstr += f'\nPower law exponent α: {convergence["alpha"]:.2f}'
    if convergence['curve_fit_success'] and not np.isnan(convergence['convergence_cluster_size_95pct']):
        textstr += f'\n95% convergence: cluster size ≥ {convergence["convergence_cluster_size_95pct"]:.0f}'

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)

    # Set reasonable y-axis limits
    y_min = max(0, np.percentile(error_rates, 1) * 100 - 0.5)
    y_max = np.percentile(error_rates, 99) * 100 + 0.5
    ax.set_ylim(y_min, y_max)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    logging.info(f"Wrote convergence plot to: {output_file}")


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


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze empirical error characteristics by comparing raw reads to consensus sequences."
    )
    parser.add_argument(
        'output_dir',
        help="Path to speconsense output directory"
    )
    parser.add_argument(
        '--output-dir', '-o',
        dest='analysis_dir',
        default='analysis_results',
        help="Directory for analysis output files (default: analysis_results)"
    )
    parser.add_argument(
        '--use-sampled',
        action='store_true',
        help="Use sampled reads (*-sampled.fastq) instead of full cluster (*-reads.fastq)"
    )
    parser.add_argument(
        '--min-ric',
        type=int,
        default=0,
        help="Minimum RiC threshold for analysis (default: 0, analyze all clusters)"
    )
    parser.add_argument(
        '--per-read',
        action='store_true',
        help="Output detailed per-read alignment data (can be large)"
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


def main():
    """Main entry point for speconsense-analyze."""
    args = parse_arguments()
    setup_logging(args.log_level)

    logging.info("Starting Speconsense empirical error analysis")
    logging.info(f"Input directory: {args.output_dir}")
    logging.info(f"Output directory: {args.analysis_dir}")

    # Create output directory
    os.makedirs(args.analysis_dir, exist_ok=True)

    # Load consensus sequences
    logging.info("Loading consensus sequences...")
    consensus_map = load_consensus_sequences(args.output_dir)
    if not consensus_map:
        logging.error("No consensus sequences found")
        return 1
    logging.info(f"Found {len(consensus_map)} consensus sequences")

    # Find cluster FASTQ files
    logging.info("Finding cluster files...")
    cluster_files = find_cluster_files(args.output_dir, args.use_sampled)
    if not cluster_files:
        logging.error("No cluster files found")
        return 1
    logging.info(f"Found {len(cluster_files)} cluster files")

    # Match consensus sequences with their reads files
    clusters_to_analyze = []
    for key, reads_file in cluster_files.items():
        if key in consensus_map:
            cluster_info = consensus_map[key]._replace(reads_file=reads_file)

            # Apply RiC filter
            if cluster_info.ric >= args.min_ric:
                clusters_to_analyze.append(cluster_info)

    if not clusters_to_analyze:
        logging.error("No matching clusters found")
        return 1

    logging.info(f"Analyzing {len(clusters_to_analyze)} clusters (RiC >= {args.min_ric})")

    # Analyze each cluster
    all_cluster_stats = []
    all_read_alignments = []
    all_cluster_alignments = []  # Store (cluster_info, alignments) for positional analysis

    for cluster_info in tqdm(clusters_to_analyze, desc="Analyzing clusters", unit="cluster"):
        stats, alignments = analyze_cluster(cluster_info)

        if stats:
            all_cluster_stats.append(stats)
            all_cluster_alignments.append((cluster_info, alignments))

            if args.per_read:
                all_read_alignments.extend([
                    (cluster_info.sample_name, cluster_info.cluster_num, a)
                    for a in alignments
                ])

    if not all_cluster_stats:
        logging.error("No clusters successfully analyzed")
        return 1

    # Write outputs
    logging.info("Writing results...")

    # Cluster statistics CSV
    cluster_stats_file = os.path.join(args.analysis_dir, 'cluster_statistics.csv')
    write_cluster_statistics(all_cluster_stats, cluster_stats_file)

    # Perform convergence analysis (for both summary report and plot)
    convergence = analyze_read_depth_convergence(all_cluster_stats)

    # Dataset summary TXT
    summary_file = os.path.join(args.analysis_dir, 'analysis_report.txt')
    write_dataset_summary_with_convergence(all_cluster_stats, convergence, summary_file)

    # Convergence plot
    plot_file = os.path.join(args.analysis_dir, 'convergence_plot.png')
    plot_convergence_analysis(all_cluster_stats, convergence, plot_file)

    # Positional analysis across ALL clusters
    logging.info("Performing positional analysis across all clusters...")
    all_position_stats = analyze_all_cluster_positions(all_cluster_alignments)

    if all_position_stats:
        # Generate histogram
        histogram_file = os.path.join(args.analysis_dir, 'positional_error_histogram.png')
        plot_positional_error_histogram(all_position_stats, histogram_file)
        logging.info(f"Generated positional error histogram: {histogram_file}")

        # Calculate error rate statistics
        error_rates = [ps.error_rate * 100 for ps in all_position_stats]
        p95_threshold = np.percentile(error_rates, 95)

        # Filter positions exceeding p95 threshold
        high_error_positions = [ps for ps in all_position_stats if ps.error_rate * 100 >= p95_threshold]

        # Categorize positions by priority
        high_priority = []  # Internal, non-homopolymer
        medium_priority = []  # Homopolymer (but not at ends)
        low_priority = []  # At ends (first or last 20bp)

        for ps in high_error_positions:
            # Check if position is near either end
            is_at_end = ps.position < 20 or ps.position >= (ps.consensus_length - 20)

            if is_at_end:
                low_priority.append(ps)
            elif ps.is_homopolymer:
                medium_priority.append(ps)
            else:
                high_priority.append(ps)

        # Sort each category by error rate (descending)
        high_priority.sort(key=lambda x: x.error_rate, reverse=True)
        medium_priority.sort(key=lambda x: x.error_rate, reverse=True)
        low_priority.sort(key=lambda x: x.error_rate, reverse=True)

        # Write high-error positions to separate file
        high_error_file = os.path.join(args.analysis_dir, 'high_error_positions.txt')
        with open(high_error_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("HIGH-ERROR POSITIONS REPORT\n")
            f.write(f"Positions exceeding p95 error rate threshold ({p95_threshold:.2f}%)\n")
            f.write("="*80 + "\n\n")

            f.write(f"Total high-error positions: {len(high_error_positions):,}\n")
            f.write(f"  High priority (internal, non-homopolymer): {len(high_priority)}\n")
            f.write(f"  Medium priority (homopolymer regions): {len(medium_priority)}\n")
            f.write(f"  Low priority (sequence ends): {len(low_priority)}\n\n")

            # Helper function to write section
            def write_section(positions, title, description):
                if not positions:
                    return

                f.write("\n" + "="*80 + "\n")
                f.write(f"{title}\n")
                f.write(f"{description}\n")
                f.write("="*80 + "\n\n")
                f.write(f"Count: {len(positions)}\n\n")
                f.write(f"{'Rank':<6} {'Cluster':<40} {'Pos':<6} {'FromEnd':<8} {'Error%':<8} {'Size':<6} {'Sub':<5} {'Ins':<5} {'Del':<5} {'HP':<4} {'Base':<4}\n")
                f.write("-" * 125 + "\n")

                for rank, ps in enumerate(positions, 1):
                    cluster_id = f"{ps.sample_name}-c{ps.cluster_num}"
                    from_end = ps.consensus_length - ps.position
                    f.write(f"{rank:<6} {cluster_id:<40} {ps.position:<6} {from_end:<8} {ps.error_rate*100:>6.2f}% {ps.coverage:<6} {ps.sub_count:<5} {ps.ins_count:<5} {ps.del_count:<5} {'Y' if ps.is_homopolymer else 'N':<4} {ps.nucleotide:<4}\n")

            # Write sections in priority order
            write_section(
                high_priority,
                "SECTION 1: HIGH PRIORITY - Internal Non-Homopolymer Positions",
                "These positions are NOT at sequence ends and NOT in homopolymer regions.\nMost likely to indicate unphased variants or clustering issues."
            )

            write_section(
                medium_priority,
                "SECTION 2: MEDIUM PRIORITY - Homopolymer Region Positions",
                "These positions are in homopolymer regions (≥4bp) but not at sequence ends.\nMay indicate homopolymer length variation or sequencing errors."
            )

            write_section(
                low_priority,
                "SECTION 3: LOW PRIORITY - Sequence End Positions",
                "These positions are within first or last 20bp of sequences.\nLikely due to read-end degradation, a known ONT sequencing artifact."
            )

            # Footer notes
            f.write("\n" + "="*80 + "\n")
            f.write("NOTES\n")
            f.write("="*80 + "\n")
            f.write("HP = Homopolymer (≥4bp run)\n")
            f.write("Pos = Position in consensus sequence (0-indexed)\n")
            f.write("FromEnd = Distance from sequence tail (length - position)\n")
            f.write("Size = Cluster size (number of reads)\n")
            f.write("\nHigh-error positions may indicate:\n")
            f.write("  - Unphased variants (multiple haplotypes in same cluster)\n")
            f.write("  - Sequencing error hotspots\n")
            f.write("  - Homopolymer length ambiguity\n")
            f.write("  - Structural variation\n")
            f.write("  - Read-end degradation (for positions at sequence ends)\n")

        logging.info(f"Generated high-error positions report: {high_error_file}")

        # Add positional analysis summary to main report
        with open(summary_file, 'a') as f:
            f.write("\n" + "="*60 + "\n")
            f.write("=== Positional Error Analysis Across All Clusters ===\n")
            f.write("="*60 + "\n\n")

            # Overall statistics
            f.write(f"Dataset-wide Position Statistics:\n")
            f.write(f"  Total positions analyzed: {len(all_position_stats):,}\n")
            f.write(f"  From {len(all_cluster_alignments)} clusters\n\n")

            f.write(f"  Positional error rate distribution:\n")
            f.write(f"    Mean: {np.mean(error_rates):.2f}%\n")
            f.write(f"    Median: {np.median(error_rates):.2f}%\n")
            f.write(f"    Std Dev: {np.std(error_rates):.2f}%\n")
            f.write(f"    95th percentile: {p95_threshold:.2f}%\n")
            f.write(f"    99th percentile: {np.percentile(error_rates, 99):.2f}%\n")
            f.write(f"    Max: {np.max(error_rates):.2f}%\n\n")

            # Clustering quality assessment
            f.write(f"  Clustering Quality Assessment:\n")
            f.write(f"    Positions with error rate > 5%: {sum(1 for er in error_rates if er > 5):,} ({sum(1 for er in error_rates if er > 5)/len(error_rates)*100:.1f}%)\n")
            f.write(f"    Positions with error rate > 10%: {sum(1 for er in error_rates if er > 10):,} ({sum(1 for er in error_rates if er > 10)/len(error_rates)*100:.1f}%)\n")
            f.write(f"    (Lower percentages indicate better clustering)\n\n")

            # High-error position summary
            f.write(f"  High-Error Positions (exceeding p95 = {p95_threshold:.2f}%):\n")
            f.write(f"    Total: {len(high_error_positions):,}\n")
            f.write(f"    High priority (internal, non-homopolymer): {len(high_priority)}\n")
            f.write(f"    Medium priority (homopolymer regions): {len(medium_priority)}\n")
            f.write(f"    Low priority (sequence ends): {len(low_priority)}\n")
            f.write(f"    → See detailed report: {os.path.basename(high_error_file)}\n")

    # Read error histogram across ALL clusters
    logging.info("Generating read error histogram...")
    read_histogram_file = os.path.join(args.analysis_dir, 'read_error_histogram.png')
    plot_read_error_histogram(all_cluster_alignments, read_histogram_file)
    logging.info(f"Generated read error histogram: {read_histogram_file}")

    # Per-read data (if requested)
    if args.per_read and all_read_alignments:
        read_file = os.path.join(args.analysis_dir, 'read_distances.csv')
        with open(read_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'sample_name', 'cluster_num', 'read_id', 'read_length',
                'edit_distance', 'error_rate', 'insertions', 'deletions', 'substitutions'
            ])
            for sample_name, cluster_num, alignment in all_read_alignments:
                error_rate = alignment.edit_distance / alignment.read_length if alignment.read_length > 0 else 0
                writer.writerow([
                    sample_name, cluster_num, alignment.read_id, alignment.read_length,
                    alignment.edit_distance, f"{error_rate:.4f}",
                    alignment.num_insertions, alignment.num_deletions, alignment.num_substitutions
                ])
        logging.info(f"Wrote per-read data to: {read_file}")

    logging.info("Analysis completed successfully!")
    logging.info(f"Results written to: {args.analysis_dir}")

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
