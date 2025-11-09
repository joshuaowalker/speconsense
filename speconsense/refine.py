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
        '--mean-error-rate',
        type=float,
        default=0.013,
        help="Estimated mean read-to-consensus error rate (default: 0.013 for 1.3%%). "
             "Used to calculate outlier detection threshold (mean × 3.0 ≈ P95). "
             "Actual mean will be computed and compared to validate estimate."
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


def is_homopolymer_variant(variant: Dict, consensus_seq: str,
                           all_variants_by_cluster: Dict, msa_to_consensus_pos: Dict[int, Optional[int]],
                           min_alt_freq: float = 0.20) -> bool:
    """
    Determine if a variant is a homopolymer-related indel.

    A variant is considered homopolymer-related if:
    1. Only one significant alternative allele exists (100% homogeneous)
    2. At least one adjacent consensus position contains that same base

    Only alternative alleles with frequency >= min_alt_freq are considered when
    determining homogeneity. This prevents low-frequency sequencing errors from
    masking clear homopolymer length differences. Since we filter to significant
    alternatives, we require 100% homogeneity - multiple significant alternatives
    indicate true biological variation, not homopolymer length variation.

    For deletions: Check if the consensus base (being deleted) appears in adjacent positions
    For insertions: Check if the inserted base appears in adjacent consensus positions
    For substitutions: Check if the substituting base appears in adjacent positions

    Args:
        variant: Variant dictionary with base_composition, msa_position, consensus_position
        consensus_seq: Ungapped consensus sequence
        all_variants_by_cluster: Dict mapping (sample_name, cluster_num, msa_pos) to variant
        msa_to_consensus_pos: Mapping from MSA positions to consensus positions
        min_alt_freq: Minimum frequency threshold for significant alternative alleles (default 0.20)

    Returns:
        True if this is a homopolymer variant
    """
    base_comp = variant['base_composition']
    cons_pos = variant['consensus_position']
    msa_pos = variant['msa_position']
    coverage = variant.get('coverage', sum(base_comp.values()))

    # Get alternative alleles (excluding the consensus base)
    alt_bases = [(base, count) for base, count in base_comp.items()
                 if base != variant.get('consensus_nucleotide', '-')]

    if not alt_bases:
        return False

    # Filter to only significant alternative alleles (>= min_alt_freq)
    # This prevents low-frequency sequencing errors from affecting HP classification
    significant_alt_bases = [(base, count) for base, count in alt_bases
                            if count / coverage >= min_alt_freq]

    if not significant_alt_bases:
        return False

    # Sort by count and get dominant
    significant_alt_bases.sort(key=lambda x: x[1], reverse=True)
    dominant_base, dominant_count = significant_alt_bases[0]

    # Check if variant is homogeneous (100% same base among significant alternatives)
    # Since we've already filtered out low-frequency sequencing errors, multiple
    # significant alternatives indicate true biological variation, not HP length variation
    total_significant = sum(count for base, count in significant_alt_bases)
    if dominant_count / total_significant < 1.0:
        return False

    # Determine which base to look for in adjacent positions
    # For deletions (dominant_base == '-'), check if adjacent positions have the consensus base
    # For insertions/substitutions, check if adjacent positions have the dominant alternative base
    if dominant_base == '-':
        # This is a deletion - check if adjacent positions have the consensus nucleotide
        base_to_check = variant.get('consensus_nucleotide', None)
        if base_to_check is None or base_to_check == '-':
            return False
    else:
        # This is an insertion or substitution - check if adjacent positions have the dominant base
        base_to_check = dominant_base

    # Now check if adjacent consensus positions have this base
    # For insertions, we need to find adjacent consensus positions via MSA mapping
    if cons_pos is None:
        # This is an insertion - find adjacent consensus positions
        left_pos, right_pos = get_adjacent_consensus_positions(msa_pos, msa_to_consensus_pos)

        # Check if either adjacent position has the base we're looking for
        adjacent_has_base = False
        if left_pos is not None and left_pos < len(consensus_seq):
            if consensus_seq[left_pos] == base_to_check:
                adjacent_has_base = True
        if right_pos is not None and right_pos < len(consensus_seq):
            if consensus_seq[right_pos] == base_to_check:
                adjacent_has_base = True

        return adjacent_has_base
    else:
        # Regular position - check adjacent consensus positions
        adjacent_has_base = False

        # Check left neighbor
        if cons_pos > 0:
            if consensus_seq[cons_pos - 1] == base_to_check:
                adjacent_has_base = True

        # Check right neighbor
        if cons_pos + 1 < len(consensus_seq):
            if consensus_seq[cons_pos + 1] == base_to_check:
                adjacent_has_base = True

        return adjacent_has_base


def extract_variant_combinations(msa_file: str, variant_msa_positions: List[int]) -> Dict[Tuple, int]:
    """
    Extract variant allele combinations across multiple positions from an MSA.

    For each read in the MSA, extract the base at each variant position and
    count the frequency of each combination (haplotype).

    Args:
        msa_file: Path to MSA FASTA file
        variant_msa_positions: List of MSA positions to check (0-based)

    Returns:
        Dictionary mapping combinations (tuples of bases) to their counts
    """
    from collections import Counter
    from Bio import SeqIO

    if not os.path.exists(msa_file):
        return {}

    combinations = []

    try:
        for record in SeqIO.parse(msa_file, 'fasta'):
            # Skip the consensus sequence (appears at end, labeled "Consensus")
            if record.id.lower() == 'consensus':
                continue

            # Extract bases at variant positions
            seq_str = str(record.seq)
            combo = tuple(seq_str[pos] if pos < len(seq_str) else '-'
                         for pos in variant_msa_positions)
            combinations.append(combo)

        # Count combinations
        return dict(Counter(combinations))

    except Exception as e:
        logging.warning(f"Failed to extract variant combinations from {msa_file}: {e}")
        return {}


def get_adjacent_consensus_positions(msa_position: int,
                                     msa_to_consensus_pos: Dict[int, Optional[int]]) -> Tuple[Optional[int], Optional[int]]:
    """
    For an insertion column in the MSA, find the adjacent consensus positions.

    Args:
        msa_position: MSA position of the insertion
        msa_to_consensus_pos: Mapping from MSA positions to consensus positions

    Returns:
        Tuple of (left_consensus_pos, right_consensus_pos)
        Either may be None if at the start/end of the sequence
    """
    # Find previous consensus position
    left_pos = None
    for i in range(msa_position - 1, -1, -1):
        if i in msa_to_consensus_pos and msa_to_consensus_pos[i] is not None:
            left_pos = msa_to_consensus_pos[i]
            break

    # Find next consensus position
    right_pos = None
    max_msa_pos = max(msa_to_consensus_pos.keys())
    for i in range(msa_position + 1, max_msa_pos + 1):
        if i in msa_to_consensus_pos and msa_to_consensus_pos[i] is not None:
            right_pos = msa_to_consensus_pos[i]
            break

    return left_pos, right_pos


def write_cluster_metadata(output_dir: str, sample_name: str, cluster_num: int,
                           trim_start: int, trim_end: int) -> None:
    """
    Update metadata file with cluster-specific trimming information.

    Args:
        output_dir: Output directory where metadata will be written
        sample_name: Full sample name (e.g., 'sample1-c1')
        cluster_num: Cluster number
        trim_start: Start coordinate of trimmed region (0-based, inclusive)
        trim_end: End coordinate of trimmed region (0-based, exclusive)
    """
    # Extract base sample name
    base_name = sample_name.split('-c')[0] if '-c' in sample_name else sample_name

    metadata_file = os.path.join(output_dir, 'cluster_debug', f'{base_name}-metadata.json')

    # Load existing metadata or create new structure
    if os.path.exists(metadata_file):
        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
        except Exception as e:
            logging.warning(f"Failed to read existing metadata {metadata_file}: {e}")
            metadata = {}
    else:
        metadata = {}

    # Initialize clusters dict if not present
    if 'clusters' not in metadata:
        metadata['clusters'] = {}

    # Store cluster-specific trim coordinates
    cluster_key = f"c{cluster_num}"
    metadata['clusters'][cluster_key] = {
        'trim_start': trim_start,
        'trim_end': trim_end
    }

    # Write updated metadata
    try:
        os.makedirs(os.path.dirname(metadata_file), exist_ok=True)
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        logging.debug(f"Updated metadata for {sample_name}: trim_start={trim_start}, trim_end={trim_end}")
    except Exception as e:
        logging.warning(f"Failed to write metadata {metadata_file}: {e}")


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


def trim_primers(sequence: str, primers: Dict[str, str]) -> Tuple[str, List[str], int, int]:
    """
    Trim primers from start and end of sequence.

    Args:
        sequence: Sequence to trim
        primers: Dict mapping primer names to sequences

    Returns:
        Tuple of (trimmed_sequence, [found_primers], trim_start, trim_end)
        where trim_start and trim_end are 0-based coordinates in the original sequence:
        - trim_start: first position kept (inclusive)
        - trim_end: last position kept (exclusive, Python slice convention)
    """
    if not primers:
        return sequence, [], 0, len(sequence)

    # Convert primer dict to list of tuples (matching core.py format)
    primer_list = list(primers.items())

    found_primers = []
    trimmed_seq = sequence

    # Track trimming coordinates in original sequence
    trim_start = 0  # First position to keep
    trim_end = len(sequence)  # Last position to keep (exclusive)

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
        trim_start = best_start_end
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
        # Convert to original sequence coordinates
        trim_end = trim_start + best_end_start
        trimmed_seq = trimmed_seq[:best_end_start]

    return trimmed_seq, found_primers, trim_start, trim_end


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
    trim_start = 0
    trim_end = len(consensus)

    if metadata and 'primers' in metadata and metadata['primers']:
        primers = metadata['primers']
        trimmed_consensus, found_primers, trim_start, trim_end = trim_primers(consensus, primers)

        # Write trimmed consensus to main output file
        output_filename = f"{cluster_info.sample_name}.fasta"
        output_file = os.path.join(output_dir, output_filename)
        with open(output_file, 'a') as f:  # Append mode for multiple clusters per sample
            f.write(f"{header}\n{trimmed_consensus}\n")

        logging.debug(f"Trimmed {cluster_info.sample_name}: {len(consensus)} -> {len(trimmed_consensus)} bp, found primers: {found_primers}")

        # Persist trim coordinates to metadata
        write_cluster_metadata(output_dir, cluster_info.sample_name,
                              cluster_info.cluster_num, trim_start, trim_end)
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
    if args.dry_run:
        logging.info("DRY RUN MODE: No files will be written")
    logging.info("="*80)

    # Calculate error rate threshold from estimated mean
    OUTLIER_THRESHOLD_MULTIPLIER = 3.0
    error_rate_threshold = args.mean_error_rate * OUTLIER_THRESHOLD_MULTIPLIER

    logging.info("")
    logging.info("Error Rate Threshold Calculation:")
    logging.info(f"  Estimated mean error rate: {args.mean_error_rate*100:.2f}%")
    logging.info(f"  Multiplier (≈P95): {OUTLIER_THRESHOLD_MULTIPLIER}×")
    logging.info(f"  Calculated threshold: {error_rate_threshold*100:.2f}%")
    logging.info(f"  (Reads exceeding this threshold will be identified as outliers)")
    logging.info("")

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

    # Compare estimated vs actual mean error rate
    actual_mean = global_stats['read_mean_error_rate']
    estimated_mean = args.mean_error_rate
    mean_error_diff = actual_mean - estimated_mean
    mean_error_pct_diff = (mean_error_diff / estimated_mean * 100) if estimated_mean > 0 else 0

    logging.info(f"\n  Mean Error Rate Validation:")
    logging.info(f"    Estimated (--mean-error-rate): {estimated_mean*100:.2f}%")
    logging.info(f"    Actual (observed): {actual_mean*100:.2f}%")
    logging.info(f"    Difference: {mean_error_diff*100:+.2f}% ({mean_error_pct_diff:+.1f}%)")
    if abs(mean_error_pct_diff) < 20:
        logging.info(f"    ✓ Estimate is reasonable (within 20%)")
    else:
        logging.info(f"    ⚠ Estimate differs significantly - consider using --mean-error-rate {actual_mean:.4f}")

    # STEP 2: Identify outlier reads
    logging.info("\nSTEP 2: Identifying outlier reads...")
    outlier_threshold = error_rate_threshold
    logging.info(f"  Using calculated threshold: {outlier_threshold*100:.2f}%")

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

    # STEP 3: Refine clusters by removing outliers and regenerating consensus
    if args.dry_run:
        logging.info("\nDRY RUN: Analysis complete, no files written")
    else:
        logging.info("\nSTEP 3: Refining clusters...")
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

        # STEP 4: Variant position analysis on refined clusters
        logging.info("\nSTEP 4: Analyzing variant positions in refined clusters...")

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
                global_error_rate=error_rate_threshold,
                confidence_sigma=3.0
            )

            logging.info(f"  Minimum coverage for variant detection: {min_coverage} reads")
            logging.info(f"  Minimum alternative allele frequency: {args.min_alt_freq if hasattr(args, 'min_alt_freq') else 0.20:.1%}")

            # Analyze each refined cluster
            all_variant_positions = []

            for key in tqdm(refined_consensus_map.keys(), desc="Analyzing variants", unit="cluster"):
                if key not in refined_msa_files:
                    continue

                cluster_info = refined_consensus_map[key]
                msa_file = refined_msa_files[key]

                # Skip low-coverage clusters
                if cluster_info.ric < min_coverage:
                    continue

                # Extract alignments from MSA
                alignments, consensus, msa_to_consensus_pos = extract_alignments_from_msa(msa_file)

                if not alignments:
                    continue

                consensus_length = len(consensus)

                # Analyze positional variation with composition tracking
                position_stats = analyze_positional_variation(
                    alignments,
                    consensus,
                    error_rate_threshold,
                    msa_file
                )

                # Identify variant positions
                for pos_stat in position_stats:
                    is_variant, variant_bases, reason = is_variant_position_with_composition(
                        pos_stat,
                        error_rate_threshold,
                        min_alt_freq=args.min_alt_freq if hasattr(args, 'min_alt_freq') else 0.20,
                        confidence_sigma=3.0
                    )

                    if is_variant:
                        # Handle insertion positions - create position string showing adjacent positions
                        if pos_stat.consensus_position is None:
                            # This is an insertion - find adjacent consensus positions
                            left_pos, right_pos = get_adjacent_consensus_positions(
                                pos_stat.msa_position, msa_to_consensus_pos
                            )
                            if left_pos is not None and right_pos is not None:
                                consensus_pos_str = f"{left_pos}-{right_pos}"
                            elif left_pos is not None:
                                consensus_pos_str = f"{left_pos}-END"
                            elif right_pos is not None:
                                consensus_pos_str = f"START-{right_pos}"
                            else:
                                consensus_pos_str = "INSERTION"
                        else:
                            consensus_pos_str = str(pos_stat.consensus_position)

                        all_variant_positions.append({
                            'sample_name': cluster_info.sample_name,
                            'cluster_num': cluster_info.cluster_num,
                            'msa_position': pos_stat.msa_position,
                            'consensus_position': pos_stat.consensus_position,  # Keep numeric for filtering
                            'consensus_position_str': consensus_pos_str,  # For display
                            'consensus_length': consensus_length,
                            'consensus_seq': consensus,  # Needed for HP detection
                            'msa_to_consensus_pos': msa_to_consensus_pos,  # Needed for HP detection
                            'coverage': pos_stat.coverage,
                            'variant_bases': variant_bases,
                            'base_composition': pos_stat.base_composition,
                            'consensus_nucleotide': pos_stat.consensus_nucleotide,  # Needed for HP detection
                            'error_rate': pos_stat.error_rate,
                            'is_homopolymer': pos_stat.is_homopolymer,  # Legacy, will be replaced
                            'reason': reason
                        })

            # Filter variants based on primer trimming
            # Remove variants that fall outside the trimmed region
            trimmed_variants = []
            filtered_count = 0

            for v in all_variant_positions:
                # Try to load trimming coordinates from metadata
                cluster_key = f"c{v['cluster_num']}"
                sample_metadata = read_metadata(args.output_dir, v['sample_name'])

                trim_start = 0
                trim_end = v['consensus_length']

                if sample_metadata and 'clusters' in sample_metadata:
                    if cluster_key in sample_metadata['clusters']:
                        trim_start = sample_metadata['clusters'][cluster_key].get('trim_start', 0)
                        trim_end = sample_metadata['clusters'][cluster_key].get('trim_end', v['consensus_length'])

                # Check if variant is within trimmed region
                cons_pos = v['consensus_position']

                if cons_pos is None:
                    # Insertion - check if adjacent positions are within trimmed region
                    left_pos, right_pos = get_adjacent_consensus_positions(
                        v['msa_position'], v['msa_to_consensus_pos']
                    )
                    # Keep if at least one adjacent position is within trimmed region
                    keep = False
                    if left_pos is not None and trim_start <= left_pos < trim_end:
                        keep = True
                    if right_pos is not None and trim_start <= right_pos < trim_end:
                        keep = True

                    if keep:
                        trimmed_variants.append(v)
                    else:
                        filtered_count += 1
                else:
                    # Regular position - check if within trimmed region
                    if trim_start <= cons_pos < trim_end:
                        trimmed_variants.append(v)
                    else:
                        filtered_count += 1

            # Report variant positions
            logging.info(f"\n  Variant Position Detection Results:")
            logging.info(f"    Total variant positions found: {len(all_variant_positions)}")
            logging.info(f"    Filtered (outside trimmed region): {filtered_count}")
            logging.info(f"    Remaining after trimming filter: {len(trimmed_variants)}")

            # Replace all_variant_positions with filtered list
            all_variant_positions = trimmed_variants

            if all_variant_positions:
                # Count unique clusters with variants
                unique_clusters_all = len(set((v['sample_name'], v['cluster_num']) for v in all_variant_positions))
                logging.info(f"    Unique clusters with variants: {unique_clusters_all}")

                # Categorize variants by type (homopolymer vs non-homopolymer)
                # Use enhanced HP detection that checks base composition and adjacent positions

                hp_variants = []      # Homopolymer variants
                non_hp_variants = []  # Non-homopolymer variants

                # Create lookup dict for multi-position analysis (currently unused but prepared)
                variants_by_key = {}
                for v in all_variant_positions:
                    key = (v['sample_name'], v['cluster_num'], v['msa_position'])
                    variants_by_key[key] = v

                for v in all_variant_positions:
                    # Use enhanced HP detection
                    is_hp = is_homopolymer_variant(
                        v,
                        v['consensus_seq'],
                        variants_by_key,
                        v['msa_to_consensus_pos'],
                        min_alt_freq=args.min_alt_freq if hasattr(args, 'min_alt_freq') else 0.20
                    )

                    if is_hp:
                        hp_variants.append(v)
                    else:
                        non_hp_variants.append(v)

                # Count unique clusters for each category
                unique_clusters_hp = len(set((v['sample_name'], v['cluster_num']) for v in hp_variants))
                unique_clusters_non_hp = len(set((v['sample_name'], v['cluster_num']) for v in non_hp_variants))

                logging.info(f"\n  Breakdown by category:")
                logging.info(f"    Non-homopolymer variants: {len(non_hp_variants)} positions in {unique_clusters_non_hp} clusters")
                logging.info(f"    Homopolymer variants: {len(hp_variants)} positions in {unique_clusters_hp} clusters")

                # Generate cluster ranking report
                if non_hp_variants:
                    # Group non-HP variants by cluster
                    from collections import defaultdict
                    clusters_with_variants = defaultdict(list)

                    for v in non_hp_variants:
                        cluster_key = (v['sample_name'], v['cluster_num'])
                        clusters_with_variants[cluster_key].append(v)

                    # Calculate variation score for each cluster
                    cluster_scores = []

                    for cluster_key, variants in clusters_with_variants.items():
                        # For each variant, get the maximum alternative allele frequency
                        total_score = 0.0

                        for v in variants:
                            base_comp = v['base_composition']
                            cons_nucleotide = v.get('consensus_nucleotide', '-')

                            # Get alternative allele counts (excluding consensus base)
                            alt_counts = {base: count for base, count in base_comp.items()
                                        if base != cons_nucleotide}

                            if alt_counts and v['coverage'] > 0:
                                # Find maximum alternative allele frequency
                                max_alt_count = max(alt_counts.values())
                                max_alt_freq = (max_alt_count / v['coverage']) * 100  # As percentage

                                total_score += max_alt_freq

                        cluster_scores.append({
                            'sample_name': cluster_key[0],
                            'cluster_num': cluster_key[1],
                            'total_score': total_score,
                            'num_variants': len(variants),
                            'variants': variants
                        })

                    # Sort by total score (descending)
                    cluster_scores.sort(key=lambda x: x['total_score'], reverse=True)

                    # Write detailed report to file
                    report_file = os.path.join(args.output_dir, 'variant_clusters_ranked.txt')
                    with open(report_file, 'w') as f:
                        f.write("Cluster Ranking by Non-Homopolymer Variation\n")
                        f.write("=" * 80 + "\n\n")

                        for rank, cluster in enumerate(cluster_scores, 1):
                            f.write(f"Rank {rank}: {cluster['sample_name']}\n")
                            f.write(f"Total variation score: {cluster['total_score']:.1f}%\n")
                            f.write(f"Non-HP variant positions: {cluster['num_variants']}\n")
                            f.write(f"Variant details:\n")

                            for v in cluster['variants']:
                                base_comp = v['base_composition']
                                cons_nucleotide = v.get('consensus_nucleotide', '-')
                                sorted_comp = sorted(base_comp.items(), key=lambda x: x[1], reverse=True)
                                comp_str = ', '.join([f"{b}:{c}" for b, c in sorted_comp[:4]])

                                # Get max alt frequency
                                alt_counts = {base: count for base, count in base_comp.items()
                                            if base != cons_nucleotide}
                                if alt_counts and v['coverage'] > 0:
                                    max_alt_count = max(alt_counts.values())
                                    max_alt_freq = (max_alt_count / v['coverage']) * 100
                                else:
                                    max_alt_freq = 0.0

                                pos_str = f"MSA:{v['msa_position']}/Cons:{v['consensus_position_str']}"
                                f.write(f"  - {pos_str}: {max_alt_freq:.1f}% (cov={v['coverage']}, bases=[{comp_str}])\n")

                            # Add variant combinations section for multiple variant positions
                            if cluster['num_variants'] > 1:
                                # Find MSA file for this cluster
                                sample_name = cluster['sample_name']
                                cluster_debug_dir = os.path.join(args.input_dir, 'cluster_debug')

                                # MSA file has RiC suffix - use glob to find it
                                msa_pattern = os.path.join(cluster_debug_dir, f"{sample_name}-RiC*-msa.fasta")
                                msa_files = glob.glob(msa_pattern)

                                if msa_files:
                                    msa_file = msa_files[0]  # Take first match
                                    # Get MSA positions from variants (sorted by position)
                                    sorted_variants = sorted(cluster['variants'], key=lambda x: x['msa_position'])
                                    variant_msa_positions = [v['msa_position'] for v in sorted_variants]

                                    # Extract combinations
                                    combinations = extract_variant_combinations(msa_file, variant_msa_positions)

                                    if combinations:
                                        f.write(f"Variant combinations:\n")

                                        # Sort by count (descending)
                                        sorted_combos = sorted(combinations.items(), key=lambda x: x[1], reverse=True)

                                        for combo, count in sorted_combos:
                                            # Format combination as (base1, base2, ...)
                                            combo_str = '(' + ', '.join(combo) + ')'
                                            f.write(f"  - {combo_str}: {count} reads\n")

                            f.write("\n")

                    logging.info(f"\n  Cluster Ranking:")
                    logging.info(f"    Full ranking written to: {report_file}")
                    logging.info(f"    Top 5 clusters by variation score:")

                    for rank, cluster in enumerate(cluster_scores[:5], 1):
                        logging.info(
                            f"      {rank}. {cluster['sample_name']}: "
                            f"{cluster['total_score']:.1f}% ({cluster['num_variants']} non-HP variants)"
                        )

                if non_hp_variants:
                    # Sort by error rate (highest first)
                    non_hp_variants.sort(key=lambda x: x['error_rate'], reverse=True)

                    logging.info(f"\n  Top 20 non-homopolymer variant positions:")
                    for i, v in enumerate(non_hp_variants[:20], 1):
                        cluster_id = f"{v['sample_name']}"
                        comp = v['base_composition']
                        sorted_comp = sorted(comp.items(), key=lambda x: x[1], reverse=True)
                        comp_str = ', '.join([f"{b}:{c}" for b, c in sorted_comp[:3]])

                        # Display both MSA and consensus position for clarity
                        pos_str = f"MSA:{v['msa_position']}/Cons:{v['consensus_position_str']}"

                        logging.info(
                            f"    {i}. {cluster_id} {pos_str}: "
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
