"""MSA analysis and quality assessment for speconsense-summarize.

Provides functions for analyzing multiple sequence alignments, detecting outliers,
identifying indel events, and assessing cluster quality.
"""

import os
import re
import logging
import subprocess
import tempfile
from typing import List, Dict, Optional, Tuple, NamedTuple
from io import StringIO

import edlib
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from speconsense.types import ConsensusInfo
from speconsense.msa import (
    extract_alignments_from_msa,
    analyze_positional_variation,
)

from .iupac import IUPAC_EQUIV


# Maximum number of variants to evaluate for MSA-based merging (legacy constant)
# Batch size is now dynamically computed based on --merge-effort and group size.
# This constant is kept for backward compatibility and as the default MAX_MERGE_BATCH.
MAX_MSA_MERGE_VARIANTS = 8

# Merge effort batch size limits
MIN_MERGE_BATCH = 4
MAX_MERGE_BATCH = 8


def compute_merge_batch_size(group_size: int, effort: int) -> int:
    """Compute batch size for a group based on effort level.

    Uses formula: B = E + 1 - log2(V), clamped to [MIN_MERGE_BATCH, MAX_MERGE_BATCH]
    This keeps expected evaluations near 2^E per group.

    Args:
        group_size: Number of variants in the HAC group
        effort: Merge effort level (6-14, default 10)

    Returns:
        Batch size between MIN_MERGE_BATCH and MAX_MERGE_BATCH
    """
    import math

    if group_size <= 1:
        return 1

    log_v = int(math.log2(group_size))
    batch = effort + 1 - log_v

    return max(MIN_MERGE_BATCH, min(MAX_MERGE_BATCH, batch))


class ClusterQualityData(NamedTuple):
    """Quality metrics for a cluster (no visualization matrix)."""
    consensus_seq: str
    position_error_rates: List[float]  # Per-position error rates (0-1) in consensus space
    position_error_counts: List[int]  # Per-position error counts in consensus space
    read_identities: List[float]  # Per-read identity scores (0-1)
    position_stats: Optional[List] = None  # Detailed PositionStats for debugging (optional)


def identify_outliers(final_consensus: List, all_raw_consensuses: List, source_folder: str) -> Dict:
    """Identify sequences with low read identity using statistical outlier detection.

    Flags sequences with mean read identity (rid) below (mean - 2*std) for the dataset.
    This identifies the ~2.5% lowest values that may warrant review.

    Note: rid_min (minimum read identity) is not used because single outlier reads
    don't significantly impact consensus quality. Positional analysis better captures
    systematic issues like mixed clusters or variants.

    Args:
        final_consensus: List of final consensus sequences
        all_raw_consensuses: List of all raw consensus sequences (unused, kept for API compatibility)
        source_folder: Source directory (unused, kept for API compatibility)

    Returns:
        Dictionary with:
        {
            'statistical_outliers': List of (cons, rid),
            'no_issues': List of consensus sequences with good quality,
            'global_stats': {'mean_rid', 'std_rid', 'stat_threshold_rid'}
        }
    """
    # Calculate global statistics for all sequences with identity metrics
    all_rids = []

    for cons in final_consensus:
        if cons.rid is not None:
            all_rids.append(cons.rid)

    # Calculate mean and std for statistical outlier detection
    mean_rid = np.mean(all_rids) if all_rids else 1.0
    std_rid = np.std(all_rids) if len(all_rids) > 1 else 0.0

    # Threshold for statistical outliers (2 standard deviations below mean)
    stat_threshold_rid = mean_rid - 2 * std_rid

    # Categorize sequences
    statistical = []
    no_issues = []

    for cons in final_consensus:
        rid = cons.rid if cons.rid is not None else 1.0

        if rid < stat_threshold_rid:
            statistical.append((cons, rid))
        else:
            no_issues.append(cons)

    return {
        'statistical_outliers': statistical,
        'no_issues': no_issues,
        'global_stats': {
            'mean_rid': mean_rid,
            'std_rid': std_rid,
            'stat_threshold_rid': stat_threshold_rid
        }
    }


def analyze_positional_identity_outliers(
    consensus_info,
    source_folder: str,
    min_variant_frequency: float,
    min_variant_count: int
) -> Optional[Dict]:
    """Analyze positional error rates and identify high-error positions.

    Args:
        consensus_info: ConsensusInfo object for the sequence
        source_folder: Source directory containing cluster_debug folder
        min_variant_frequency: Global threshold for flagging positions (from metadata)
        min_variant_count: Minimum variant count for phasing (from metadata)

    Returns:
        Dictionary with positional analysis:
        {
            'num_outlier_positions': int,
            'mean_outlier_error_rate': float,  # Mean error rate across outlier positions only
            'total_nucleotide_errors': int,    # Sum of error counts at outlier positions
            'outlier_threshold': float,
            'outlier_positions': List of (position, error_rate, error_count) tuples
        }
        Returns None if MSA file not found or analysis fails

        Note: Error rates already exclude homopolymer length differences due to
        homopolymer normalization in analyze_positional_variation()
    """
    # Skip analysis for low-RiC sequences (insufficient data for meaningful statistics)
    # Need at least 2 * min_variant_count to confidently phase two variants
    min_ric_threshold = 2 * min_variant_count
    if consensus_info.ric < min_ric_threshold:
        logging.debug(f"Skipping positional analysis for {consensus_info.sample_name}: "
                     f"RiC {consensus_info.ric} < {min_ric_threshold}")
        return None

    # Construct path to MSA file
    debug_dir = os.path.join(source_folder, "cluster_debug")

    # Try to find the MSA file
    # MSA files use the original cluster naming (e.g., "specimen-c1")
    # not the summarized naming (e.g., "specimen-1.v1")
    msa_file = None

    # Extract specimen name and cluster ID
    # consensus_info.sample_name might be "specimen-1.v1" (summarized)
    # consensus_info.cluster_id should be "-c1" (original cluster)

    # Build the base name from specimen + cluster_id
    # If sample_name is "ONT01.23-...-1.v1" and cluster_id is "-c1"
    # we need to reconstruct "ONT01.23-...-c1"

    sample_name = consensus_info.sample_name
    cluster_id = consensus_info.cluster_id

    # Remove any HAC group/variant suffix from sample_name to get specimen base
    # Pattern: "-\d+\.v\d+" (e.g., "-1.v1")
    specimen_base = re.sub(r'-\d+\.v\d+$', '', sample_name)

    # Reconstruct original cluster name
    original_cluster_name = f"{specimen_base}{cluster_id}"

    # Look for the MSA file with correct extension
    msa_fasta = os.path.join(debug_dir, f"{original_cluster_name}-RiC{consensus_info.ric}-msa.fasta")
    if os.path.exists(msa_fasta):
        msa_file = msa_fasta

    if not msa_file:
        logging.debug(f"No MSA file found for {original_cluster_name}")
        return None

    # Analyze cluster quality using core.py's positional analysis
    quality_data = analyze_cluster_quality(msa_file, consensus_info.sequence)

    if not quality_data or not quality_data.position_error_rates:
        logging.debug(f"Failed to analyze cluster quality for {original_cluster_name}")
        return None

    position_error_rates = quality_data.position_error_rates
    position_error_counts = quality_data.position_error_counts
    position_stats = quality_data.position_stats

    # Use global min_variant_frequency as threshold
    # Positions above this could be undetected/unphased variants
    threshold = min_variant_frequency
    outlier_positions = [
        (i, rate, count)
        for i, (rate, count) in enumerate(zip(position_error_rates, position_error_counts))
        if rate > threshold
    ]

    # Build detailed outlier info including base composition
    outlier_details = []
    if position_stats:
        for i, rate, count in outlier_positions:
            if i < len(position_stats):
                ps = position_stats[i]
                outlier_details.append({
                    'consensus_position': ps.consensus_position,
                    'msa_position': ps.msa_position,
                    'error_rate': rate,
                    'error_count': count,
                    'coverage': ps.coverage,
                    'consensus_nucleotide': ps.consensus_nucleotide,
                    'base_composition': dict(ps.base_composition),
                    'homopolymer_composition': dict(ps.homopolymer_composition) if ps.homopolymer_composition else {},
                    'sub_count': ps.sub_count,
                    'ins_count': ps.ins_count,
                    'del_count': ps.del_count,
                })

    # Calculate statistics for outlier positions only
    if outlier_positions:
        mean_outlier_error = np.mean([rate for _, rate, _ in outlier_positions])
        total_nucleotide_errors = sum(count for _, _, count in outlier_positions)
    else:
        mean_outlier_error = 0.0
        total_nucleotide_errors = 0

    return {
        'num_outlier_positions': len(outlier_positions),
        'mean_outlier_error_rate': mean_outlier_error,
        'total_nucleotide_errors': total_nucleotide_errors,
        'outlier_threshold': threshold,
        'outlier_positions': outlier_positions,
        'outlier_details': outlier_details,
        'consensus_seq': quality_data.consensus_seq,
        'ric': consensus_info.ric,
    }


def run_spoa_msa(sequences: List[str], alignment_mode: int = 1) -> List:
    """
    Run SPOA to create multiple sequence alignment.

    Args:
        sequences: List of DNA sequence strings
        alignment_mode: SPOA alignment mode:
            0 = local (Smith-Waterman) - best for overlap merging
            1 = global (Needleman-Wunsch) - default, for same-length sequences
            2 = semi-global - alternative for overlap merging

    Returns:
        List of SeqRecord objects with aligned sequences (including gaps)
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
        try:
            # Write sequences to temporary file
            records = [
                SeqRecord(Seq(seq), id=f"seq{i}", description="")
                for i, seq in enumerate(sequences)
            ]
            SeqIO.write(records, temp_input, "fasta")
            temp_input.flush()

            # Run SPOA with alignment output (-r 2) and specified alignment mode
            result = subprocess.run(
                ['spoa', temp_input.name, '-r', '2', '-l', str(alignment_mode)],
                capture_output=True,
                text=True,
                check=True
            )

            # Parse aligned sequences from SPOA output
            aligned_sequences = []
            lines = result.stdout.strip().split('\n')
            current_id = None
            current_seq = []

            for line in lines:
                if line.startswith('>'):
                    if current_id is not None:
                        # Skip consensus sequence (usually last)
                        if not current_id.startswith('Consensus'):
                            aligned_sequences.append(SeqRecord(
                                Seq(''.join(current_seq)),
                                id=current_id,
                                description=""
                            ))
                    current_id = line[1:]
                    current_seq = []
                elif line.strip():
                    current_seq.append(line.strip())

            # Add last sequence (if not consensus)
            if current_id is not None and not current_id.startswith('Consensus'):
                aligned_sequences.append(SeqRecord(
                    Seq(''.join(current_seq)),
                    id=current_id,
                    description=""
                ))

            return aligned_sequences

        finally:
            if os.path.exists(temp_input.name):
                os.unlink(temp_input.name)


def identify_indel_events(aligned_seqs: List, alignment_length: int) -> List[Tuple[int, int]]:
    """
    Identify consecutive runs of indel columns (events).

    An indel event is a maximal consecutive run of columns containing gaps.
    Each event represents a single biological insertion or deletion.

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        alignment_length: Length of the alignment

    Returns:
        List of (start_col, end_col) tuples, where end_col is inclusive
    """
    events = []
    in_event = False
    start_col = None

    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        has_gap = '-' in column
        has_bases = any(c != '-' for c in column)

        # Indel column: mix of gaps and bases
        if has_gap and has_bases:
            if not in_event:
                # Start new event
                in_event = True
                start_col = col_idx
        else:
            # Not an indel column (either all gaps or all bases)
            if in_event:
                # End current event
                events.append((start_col, col_idx - 1))
                in_event = False

    # Handle event that extends to end of alignment
    if in_event:
        events.append((start_col, alignment_length - 1))

    return events


def is_homopolymer_event(aligned_seqs: List, start_col: int, end_col: int) -> bool:
    """
    Classify a complete indel event as homopolymer or structural.

    An event is homopolymer if:
    1. All bases in the event region (across all sequences, all columns) are identical
    2. At least one flanking solid column has all sequences showing the same base

    This matches adjusted-identity semantics where AAA ~ AAAA.

    Examples:
        Homopolymer:  ATAAA--GC vs ATAAAAGC  (event has all A's, flanked by A)
        Structural:   ATAA-GC vs ATG-AGC     (event has A, flanked by A vs G)
        Structural:   ATC--GC vs ATCATGC     (event has A and T - not homopolymer)

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        start_col: First column of the indel event (inclusive)
        end_col: Last column of the indel event (inclusive)

    Returns:
        True if homopolymer event, False if structural
    """
    # Extract all bases from the event region (excluding gaps)
    bases_in_event = set()
    for col_idx in range(start_col, end_col + 1):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        bases_in_event.update(c for c in column if c != '-')

    # Must have exactly one base type across the entire event
    if len(bases_in_event) != 1:
        return False

    event_base = list(bases_in_event)[0]
    alignment_length = len(aligned_seqs[0].seq)

    # Check flanking columns for matching homopolymer context
    # A valid flanking column must:
    # 1. Not be an indel column (all sequences have bases, no gaps)
    # 2. All bases match the event base

    # Check left flank
    if start_col > 0:
        left_col = start_col - 1
        left_column = [str(seq.seq[left_col]) for seq in aligned_seqs]
        left_bases = set(c for c in left_column if c != '-')
        left_has_gap = '-' in left_column

        if not left_has_gap and left_bases == {event_base}:
            return True

    # Check right flank
    if end_col < alignment_length - 1:
        right_col = end_col + 1
        right_column = [str(seq.seq[right_col]) for seq in aligned_seqs]
        right_bases = set(c for c in right_column if c != '-')
        right_has_gap = '-' in right_column

        if not right_has_gap and right_bases == {event_base}:
            return True

    # No valid homopolymer flanking found
    return False


def analyze_msa_columns(aligned_seqs: List) -> dict:
    """
    Analyze aligned sequences to count SNPs and indels.

    Distinguishes between structural indels (real insertions/deletions) and
    homopolymer indels (length differences in homopolymer runs like AAA vs AAAA).

    Uses event-based classification: consecutive indel columns are grouped into
    events, and each complete event is classified as homopolymer or structural.

    Important: All gaps (including terminal gaps) count as variant positions
    since variants within a group share the same primers.

    Returns dict with:
        'snp_count': number of positions with >1 non-gap base
        'structural_indel_count': number of structural indel events
        'structural_indel_length': length of longest structural indel event
        'homopolymer_indel_count': number of homopolymer indel events
        'homopolymer_indel_length': length of longest homopolymer indel event
        'indel_count': total indel events (for backward compatibility)
        'max_indel_length': max indel event length (for backward compatibility)
    """
    alignment_length = len(aligned_seqs[0].seq)

    # Step 1: Count SNPs
    snp_count = 0
    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        unique_bases = set(c for c in column if c != '-')
        has_gap = '-' in column

        # SNP position: multiple different bases with NO gaps
        # Columns with gaps are indels, not SNPs
        if len(unique_bases) > 1 and not has_gap:
            snp_count += 1

    # Step 2: Identify indel events (consecutive runs of indel columns)
    indel_events = identify_indel_events(aligned_seqs, alignment_length)

    # Step 3: Classify each event as homopolymer or structural
    structural_events = []
    homopolymer_events = []

    for start_col, end_col in indel_events:
        if is_homopolymer_event(aligned_seqs, start_col, end_col):
            homopolymer_events.append((start_col, end_col))
        else:
            structural_events.append((start_col, end_col))

    # Step 4: Calculate statistics
    # Count is number of events (not columns)
    structural_indel_count = len(structural_events)
    homopolymer_indel_count = len(homopolymer_events)

    # Length is the size of the longest event
    structural_indel_length = max((end - start + 1 for start, end in structural_events), default=0)
    homopolymer_indel_length = max((end - start + 1 for start, end in homopolymer_events), default=0)

    # Backward compatibility: total events and max length
    total_indel_count = structural_indel_count + homopolymer_indel_count
    max_indel_length = max(structural_indel_length, homopolymer_indel_length)

    return {
        'snp_count': snp_count,
        'structural_indel_count': structural_indel_count,
        'structural_indel_length': structural_indel_length,
        'homopolymer_indel_count': homopolymer_indel_count,
        'homopolymer_indel_length': homopolymer_indel_length,
        'indel_count': total_indel_count,  # Backward compatibility
        'max_indel_length': max_indel_length  # Backward compatibility
    }


def analyze_msa_columns_overlap_aware(aligned_seqs: List, min_overlap_bp: int,
                                       original_lengths: List[int]) -> dict:
    """
    Analyze MSA columns, distinguishing terminal gaps from structural indels.

    Terminal gaps (from length differences at sequence ends) are NOT counted
    as structural indels when sequences have sufficient overlap in their
    shared region. This enables merging sequences from primer pools with
    different endpoints.

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        min_overlap_bp: Minimum overlap required (0 to disable overlap mode)
        original_lengths: Original ungapped sequence lengths

    Returns dict with:
        'snp_count': SNPs in overlap region
        'structural_indel_count': Structural indels in overlap region only
        'structural_indel_length': Length of longest structural indel
        'homopolymer_indel_count': Homopolymer indels (anywhere)
        'homopolymer_indel_length': Length of longest homopolymer indel
        'terminal_gap_columns': Number of terminal gap columns (not counted as structural)
        'overlap_bp': Size of overlap region in base pairs
        'prefix_bp': Extension before overlap region (for logging)
        'suffix_bp': Extension after overlap region (for logging)
        'content_regions': List of (start, end) tuples per sequence (for span logging)
        'indel_count': Total events (backward compatibility)
        'max_indel_length': Max event length (backward compatibility)
    """
    alignment_length = len(aligned_seqs[0].seq)

    # Step 1: Find content region for each sequence (first non-gap to last non-gap)
    content_regions = []  # List of (start, end) tuples
    for seq in aligned_seqs:
        seq_str = str(seq.seq)
        # Find first and last non-gap positions
        first_base = next((i for i, c in enumerate(seq_str) if c != '-'), 0)
        last_base = alignment_length - 1 - next(
            (i for i, c in enumerate(reversed(seq_str)) if c != '-'), 0
        )
        content_regions.append((first_base, last_base))

    # Step 2: Calculate overlap region (intersection of all content regions)
    overlap_start = max(start for start, _ in content_regions)
    overlap_end = min(end for _, end in content_regions)

    # Calculate union region (for prefix/suffix extension reporting)
    union_start = min(start for start, _ in content_regions)
    union_end = max(end for _, end in content_regions)
    prefix_bp = overlap_start - union_start
    suffix_bp = union_end - overlap_end

    # Calculate actual overlap in base pairs (count only columns where all have bases)
    overlap_bp = 0
    if overlap_end >= overlap_start:
        for col_idx in range(overlap_start, overlap_end + 1):
            column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
            if all(c != '-' for c in column):
                overlap_bp += 1

    # Determine effective threshold for containment cases
    shorter_len = min(original_lengths)
    effective_threshold = min(min_overlap_bp, shorter_len)

    # Step 3: Count SNPs only within overlap region
    snp_count = 0
    for col_idx in range(overlap_start, overlap_end + 1):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        unique_bases = set(c for c in column if c != '-')
        has_gap = '-' in column

        # SNP position: multiple different bases with NO gaps
        if len(unique_bases) > 1 and not has_gap:
            snp_count += 1

    # Step 4: Identify indel events, but only count those within overlap region
    indel_events = identify_indel_events(aligned_seqs, alignment_length)

    # Step 5: Classify each event and determine if it's in overlap region
    structural_events = []
    homopolymer_events = []
    terminal_gap_columns = 0

    for start_col, end_col in indel_events:
        # Check if this event is entirely within the overlap region
        is_in_overlap = (start_col >= overlap_start and end_col <= overlap_end)

        # Check if this is a terminal gap event (at the boundary of a content region)
        is_terminal = False
        for seq_start, seq_end in content_regions:
            # Terminal if event is adjacent to or outside a sequence's content region
            if end_col < seq_start or start_col > seq_end:
                is_terminal = True
                break
            # Also terminal if event is at the very edge of content
            if start_col == seq_start or end_col == seq_end:
                # Check if the gaps in this event are from this sequence's terminal
                for col_idx in range(start_col, end_col + 1):
                    column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
                    for i, (s, e) in enumerate(content_regions):
                        if col_idx < s or col_idx > e:
                            if column[i] == '-':
                                is_terminal = True
                                break
                    if is_terminal:
                        break

        if is_terminal and overlap_bp >= effective_threshold:
            # Terminal gap from length difference - don't count as structural
            terminal_gap_columns += (end_col - start_col + 1)
        elif is_homopolymer_event(aligned_seqs, start_col, end_col):
            homopolymer_events.append((start_col, end_col))
        else:
            # Only count as structural if within overlap region
            if is_in_overlap:
                structural_events.append((start_col, end_col))
            else:
                # Outside overlap - this is a terminal gap
                terminal_gap_columns += (end_col - start_col + 1)

    # Step 6: Calculate statistics
    structural_indel_count = len(structural_events)
    homopolymer_indel_count = len(homopolymer_events)

    structural_indel_length = max((end - start + 1 for start, end in structural_events), default=0)
    homopolymer_indel_length = max((end - start + 1 for start, end in homopolymer_events), default=0)

    # Backward compatibility
    total_indel_count = structural_indel_count + homopolymer_indel_count
    max_indel_length = max(structural_indel_length, homopolymer_indel_length)

    return {
        'snp_count': snp_count,
        'structural_indel_count': structural_indel_count,
        'structural_indel_length': structural_indel_length,
        'homopolymer_indel_count': homopolymer_indel_count,
        'homopolymer_indel_length': homopolymer_indel_length,
        'terminal_gap_columns': terminal_gap_columns,
        'overlap_bp': overlap_bp,
        'prefix_bp': prefix_bp,
        'suffix_bp': suffix_bp,
        'content_regions': content_regions,
        'indel_count': total_indel_count,  # Backward compatibility
        'max_indel_length': max_indel_length  # Backward compatibility
    }


def analyze_cluster_quality(
    msa_file: str,
    consensus_seq: str,
    max_reads: Optional[int] = None
) -> Optional[ClusterQualityData]:
    """
    Analyze cluster quality using core.py's analyze_positional_variation().

    Uses the canonical positional analysis from core.py to ensure consistent
    treatment of homopolymer length differences across the pipeline.

    Args:
        msa_file: Path to MSA FASTA file
        consensus_seq: Ungapped consensus sequence
        max_reads: Maximum reads to include (for downsampling large clusters)

    Returns:
        ClusterQualityData with position error rates and read identities, or None if failed
    """
    if not os.path.exists(msa_file):
        logging.debug(f"MSA file not found: {msa_file}")
        return None

    # Load MSA file content
    try:
        with open(msa_file, 'r') as f:
            msa_string = f.read()
    except Exception as e:
        logging.debug(f"Failed to read MSA file {msa_file}: {e}")
        return None

    # Extract alignments from MSA using core.py function with homopolymer normalization
    # This returns ReadAlignment objects with score_aligned field
    alignments, msa_consensus, msa_to_consensus_pos = extract_alignments_from_msa(
        msa_string,
        enable_homopolymer_normalization=True
    )

    if not alignments:
        logging.debug(f"No alignments found in MSA: {msa_file}")
        return None

    # Verify consensus matches: the passed-in consensus_seq may be trimmed (shorter) with IUPAC codes
    # The MSA consensus is untrimmed (longer) without IUPAC codes
    # Use edlib in HW mode to check if trimmed consensus is contained within MSA consensus
    if msa_consensus and msa_consensus != consensus_seq:
        # Use edlib HW mode (semi-global) to find consensus_seq within msa_consensus
        # This handles primer trimming (length difference) and IUPAC codes (via equivalencies)
        result = edlib.align(consensus_seq, msa_consensus, mode="HW", task="distance",
                             additionalEqualities=IUPAC_EQUIV)
        edit_distance = result["editDistance"]
        if edit_distance > 0:  # Any edits indicate a real mismatch
            logging.warning(f"Consensus mismatch in MSA file: {msa_file}")
            logging.warning(f"  MSA length: {len(msa_consensus)}, consensus length: {len(consensus_seq)}, edit distance: {edit_distance}")

    # Use the passed-in consensus (with IUPAC codes) as authoritative for quality analysis
    # This reflects the actual output sequence
    consensus_length = len(consensus_seq)

    if consensus_length == 0:
        logging.debug(f"Empty consensus sequence: {msa_file}")
        return None

    # Get consensus aligned sequence by parsing MSA string
    msa_handle = StringIO(msa_string)
    records = list(SeqIO.parse(msa_handle, 'fasta'))
    consensus_aligned = None
    for record in records:
        if 'Consensus' in record.description or 'Consensus' in record.id:
            consensus_aligned = str(record.seq).upper()
            break

    if consensus_aligned is None:
        logging.debug(f"No consensus found in MSA: {msa_file}")
        return None

    # Downsample reads if needed
    if max_reads and len(alignments) > max_reads:
        # Sort by read identity (using normalized edit distance) and take worst reads first, then best
        # This gives us a representative sample showing the quality range
        read_identities_temp = []
        for alignment in alignments:
            # Use normalized edit distance for identity calculation
            identity = 1.0 - (alignment.normalized_edit_distance / consensus_length) if consensus_length > 0 else 0.0
            read_identities_temp.append((identity, alignment))

        # Sort by identity
        read_identities_temp.sort(key=lambda x: x[0])

        # Take worst half and best half
        n_worst = max_reads // 2
        n_best = max_reads - n_worst
        sampled = read_identities_temp[:n_worst] + read_identities_temp[-n_best:]

        alignments = [alignment for _, alignment in sampled]
        logging.debug(f"Downsampled {len(read_identities_temp)} reads to {len(alignments)} for analysis")

    # Use core.py's canonical positional analysis
    position_stats = analyze_positional_variation(alignments, consensus_aligned, msa_to_consensus_pos)

    # Extract position error rates and counts for consensus positions only (skip insertion columns)
    consensus_position_stats = [ps for ps in position_stats if ps.consensus_position is not None]
    # Sort by consensus position to ensure correct order
    consensus_position_stats.sort(key=lambda ps: ps.consensus_position)
    position_error_rates = [ps.error_rate for ps in consensus_position_stats]
    position_error_counts = [ps.error_count for ps in consensus_position_stats]

    # Calculate per-read identities from alignments
    read_identities = []
    for alignment in alignments:
        # Use normalized edit distance for identity calculation
        identity = 1.0 - (alignment.normalized_edit_distance / consensus_length) if consensus_length > 0 else 0.0
        read_identities.append(identity)

    return ClusterQualityData(
        consensus_seq=consensus_seq,
        position_error_rates=position_error_rates,
        position_error_counts=position_error_counts,
        read_identities=read_identities,
        position_stats=consensus_position_stats
    )
