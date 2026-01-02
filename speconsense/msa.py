#!/usr/bin/env python3
"""
MSA (Multiple Sequence Alignment) Analysis Module for Speconsense.

This module contains functions and data structures for analyzing MSA output from SPOA,
including:
- Homopolymer-normalized error detection
- Positional variation analysis
- Variant position detection and phasing support
- IUPAC ambiguity code generation

These functions were extracted from core.py to improve code organization and testability.
"""

from collections import defaultdict
import logging
from typing import List, Set, Tuple, Optional, Dict, NamedTuple

import edlib
from adjusted_identity import score_alignment, AdjustmentParams
import numpy as np
from Bio import SeqIO


# IUPAC nucleotide ambiguity codes mapping
# Maps sets of nucleotides to their corresponding IUPAC code
IUPAC_CODES = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N',
}


class ErrorPosition(NamedTuple):
    """An error at a specific position in the MSA."""
    msa_position: int  # 0-indexed position in MSA alignment
    error_type: str  # 'sub', 'ins', or 'del'


class ReadAlignment(NamedTuple):
    """Alignment result for a single read against consensus."""
    read_id: str
    aligned_sequence: str  # Gapped sequence from MSA
    read_length: int

    # Raw metrics (count all differences including homopolymer length)
    edit_distance: int
    num_insertions: int
    num_deletions: int
    num_substitutions: int
    error_positions: List[ErrorPosition]  # Detailed error information

    # Homopolymer-normalized metrics (exclude homopolymer extensions)
    normalized_edit_distance: int  # Edit distance excluding homopolymer length differences
    normalized_error_positions: List[ErrorPosition]  # Only non-homopolymer errors
    score_aligned: str  # Scoring string from adjusted-identity ('|'=match, '='=homopolymer, ' '=error)


class PositionStats(NamedTuple):
    """Statistics for a single position in the MSA."""
    msa_position: int  # Position in MSA (0-indexed)
    consensus_position: Optional[int]  # Position in consensus (None for insertion columns)
    coverage: int
    error_count: int
    error_rate: float
    sub_count: int
    ins_count: int
    del_count: int
    consensus_nucleotide: str  # Base in consensus at this MSA position (or '-' for insertion)
    base_composition: Dict[str, int]  # Raw base counts: {A: 50, C: 3, G: 45, T: 2, '-': 0}
    homopolymer_composition: Dict[str, int]  # HP extension counts: {A: 5, G: 2} (base and count)


class MSAResult(NamedTuple):
    """Result from SPOA multiple sequence alignment.

    Attributes:
        consensus: Ungapped consensus sequence
        msa_string: Raw MSA in FASTA format (for file writing)
        alignments: Parsed read alignments with gapped sequences
        msa_to_consensus_pos: Mapping from MSA position to consensus position
    """
    consensus: str
    msa_string: str
    alignments: List[ReadAlignment]
    msa_to_consensus_pos: Dict[int, Optional[int]]


# ============================================================================
# MSA Analysis Functions
# ============================================================================

def parse_score_aligned_for_errors(
    score_aligned: str,
    read_aligned: str,
    consensus_aligned: str
) -> List[ErrorPosition]:
    """
    Parse score_aligned string to extract non-homopolymer errors.

    The score_aligned string from adjusted-identity uses these codes:
    - '|' : Exact match (not an error)
    - '=' : Ambiguous match or homopolymer extension (not counted as error)
    - ' ' : Substitution or indel (IS an error)
    - '.' : End-trimmed position (not counted)

    Args:
        score_aligned: Scoring string from adjusted-identity
        read_aligned: Aligned read sequence with gaps
        consensus_aligned: Aligned consensus sequence with gaps

    Returns:
        List of ErrorPosition for positions marked as errors (excluding homopolymer extensions)
    """
    normalized_errors = []

    for msa_pos, (score_char, read_base, cons_base) in enumerate(
        zip(score_aligned, read_aligned, consensus_aligned)
    ):
        # Skip matches and homopolymer extensions
        if score_char in ('|', '=', '.'):
            continue

        # This is a real error (substitution or indel) - classify it
        if read_base == '-' and cons_base != '-':
            error_type = 'del'
        elif read_base != '-' and cons_base == '-':
            error_type = 'ins'
        elif read_base != cons_base:
            error_type = 'sub'
        else:
            # Both are gaps or identical - should not happen if score_char indicates error
            continue

        normalized_errors.append(ErrorPosition(msa_pos, error_type))

    return normalized_errors


def extract_alignments_from_msa(
    msa_string: str,
    enable_homopolymer_normalization: bool = True
) -> Tuple[List[ReadAlignment], str, Dict[int, Optional[int]]]:
    """
    Extract read alignments from an MSA string with optional homopolymer normalization.

    The MSA contains aligned sequences where the consensus has header containing "Consensus".
    This function compares each read to the consensus at each aligned position.

    Error classification (raw metrics):
    - Both '-': Not an error (read doesn't cover this position)
    - Read '-', consensus base: Deletion (missing base in read)
    - Read base, consensus '-': Insertion (extra base in read)
    - Different bases: Substitution
    - Same base: Match (not an error)

    When enable_homopolymer_normalization=True, also computes normalized metrics that
    exclude homopolymer length differences using adjusted-identity library.

    IMPORTANT: Errors are reported at MSA positions, not consensus positions.
    This avoids ambiguity when multiple insertion columns map to the same consensus position.

    Args:
        msa_string: MSA content in FASTA format
        enable_homopolymer_normalization: If True, compute homopolymer-normalized metrics

    Returns:
        Tuple of:
        - list of ReadAlignment objects (with both raw and normalized metrics)
        - consensus sequence without gaps
        - mapping from MSA position to consensus position (None for insertion columns)
    """
    from io import StringIO

    # Define adjustment parameters for homopolymer normalization
    # Only normalize homopolymers (single-base repeats), no other adjustments
    HOMOPOLYMER_ADJUSTMENT_PARAMS = AdjustmentParams(
        normalize_homopolymers=True,
        handle_iupac_overlap=False,
        normalize_indels=False,
        end_skip_distance=0,
        max_repeat_motif_length=1  # Single-base repeats only
    )

    # Parse MSA
    msa_handle = StringIO(msa_string)
    records = list(SeqIO.parse(msa_handle, 'fasta'))

    if not records:
        logging.warning("No sequences found in MSA string")
        return [], "", {}

    # Find consensus sequence
    consensus_record = None
    read_records = []

    for record in records:
        if 'Consensus' in record.description or 'Consensus' in record.id:
            consensus_record = record
        else:
            read_records.append(record)

    if consensus_record is None:
        logging.warning("No consensus sequence found in MSA string")
        return [], "", {}

    consensus_aligned = str(consensus_record.seq).upper()
    msa_length = len(consensus_aligned)

    # Build mapping from MSA position to consensus position (excluding gaps)
    # For insertion columns (consensus has '-'), maps to None
    msa_to_consensus_pos = {}
    consensus_pos = 0
    for msa_pos in range(msa_length):
        if consensus_aligned[msa_pos] != '-':
            msa_to_consensus_pos[msa_pos] = consensus_pos
            consensus_pos += 1
        else:
            # Insertion column - no consensus position
            msa_to_consensus_pos[msa_pos] = None

    # Get consensus without gaps for return value
    consensus_ungapped = consensus_aligned.replace('-', '')

    # Process each read
    alignments = []

    for read_record in read_records:
        read_aligned = str(read_record.seq).upper()

        if len(read_aligned) != msa_length:
            logging.warning(f"Read {read_record.id} length mismatch with MSA length")
            continue

        # Compare read to consensus at each position
        error_positions = []
        num_insertions = 0
        num_deletions = 0
        num_substitutions = 0

        for msa_pos in range(msa_length):
            read_base = read_aligned[msa_pos]
            cons_base = consensus_aligned[msa_pos]

            # Skip if both are gaps (read doesn't cover this position)
            if read_base == '-' and cons_base == '-':
                continue

            # Classify error type and record at MSA position
            if read_base == '-' and cons_base != '-':
                # Deletion (missing base in read)
                error_positions.append(ErrorPosition(msa_pos, 'del'))
                num_deletions += 1
            elif read_base != '-' and cons_base == '-':
                # Insertion (extra base in read)
                error_positions.append(ErrorPosition(msa_pos, 'ins'))
                num_insertions += 1
            elif read_base != cons_base:
                # Substitution (different bases)
                error_positions.append(ErrorPosition(msa_pos, 'sub'))
                num_substitutions += 1
            # else: match, no error

        # Calculate edit distance and read length
        edit_distance = num_insertions + num_deletions + num_substitutions
        read_length = len(read_aligned.replace('-', ''))  # Length without gaps

        # Compute homopolymer-normalized metrics if enabled
        if enable_homopolymer_normalization:
            try:
                # Use adjusted-identity to get homopolymer-normalized scoring
                # IMPORTANT: seq1=read, seq2=consensus. The score_aligned visualization
                # is asymmetric and shows HP extensions from seq1's (the READ's) perspective.
                # This is what we want since we're identifying which READ bases are extensions.
                result = score_alignment(
                    read_aligned,      # seq1 - the read
                    consensus_aligned, # seq2 - the consensus
                    HOMOPOLYMER_ADJUSTMENT_PARAMS
                )

                # Parse score_aligned string to extract normalized errors
                normalized_error_positions = parse_score_aligned_for_errors(
                    result.score_aligned,
                    read_aligned,
                    consensus_aligned
                )

                normalized_edit_distance = result.mismatches
                score_aligned_str = result.score_aligned

            except Exception as e:
                # If normalization fails, fall back to raw metrics
                logging.warning(f"Homopolymer normalization failed for read {read_record.id}: {e}")
                normalized_edit_distance = edit_distance
                normalized_error_positions = error_positions
                score_aligned_str = ""
        else:
            # Homopolymer normalization disabled - use raw metrics
            normalized_edit_distance = edit_distance
            normalized_error_positions = error_positions
            score_aligned_str = ""

        # Create alignment object with both raw and normalized metrics
        alignment = ReadAlignment(
            read_id=read_record.id,
            aligned_sequence=read_aligned,  # Store gapped sequence
            read_length=read_length,
            # Raw metrics
            edit_distance=edit_distance,
            num_insertions=num_insertions,
            num_deletions=num_deletions,
            num_substitutions=num_substitutions,
            error_positions=error_positions,
            # Normalized metrics
            normalized_edit_distance=normalized_edit_distance,
            normalized_error_positions=normalized_error_positions,
            score_aligned=score_aligned_str
        )
        alignments.append(alignment)

    return alignments, consensus_ungapped, msa_to_consensus_pos


def analyze_positional_variation(alignments: List[ReadAlignment], consensus_aligned: str,
                                 msa_to_consensus_pos: Dict[int, Optional[int]]) -> List[PositionStats]:
    """
    Analyze error rates at each position in the MSA with homopolymer tracking.

    Uses normalized error positions and base composition to identify true variants
    while tracking homopolymer length differences separately. For each position:
    - base_composition: Raw counts of each base observed
    - homopolymer_composition: Counts of bases that are homopolymer extensions (score_aligned='=')

    Downstream variant detection uses effective counts (raw - HP) to identify true
    biological variants while ignoring diversity due solely to homopolymer variation.

    IMPORTANT: All analysis is performed in MSA space (not consensus space).
    This correctly handles insertion columns where multiple MSA positions
    don't correspond to any consensus position.

    Args:
        alignments: List of read alignments (with normalized metrics)
        consensus_aligned: Consensus sequence (gapped, from MSA)
        msa_to_consensus_pos: Mapping from MSA position to consensus position

    Returns:
        List of PositionStats for each MSA position with normalized base composition
    """
    msa_length = len(consensus_aligned)

    # Build error frequency matrix in MSA space
    # For each MSA position: [sub_count, ins_count, del_count, total_coverage]
    error_matrix = np.zeros((msa_length, 4), dtype=int)

    # Build base composition matrix in MSA space (raw counts)
    base_composition_matrix = [
        {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
        for _ in range(msa_length)
    ]

    # Build homopolymer composition matrix in MSA space
    # Tracks bases that are homopolymer extensions (score_aligned='=')
    homopolymer_composition_matrix = [
        {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
        for _ in range(msa_length)
    ]

    # Process alignments to count errors at MSA positions
    for read_idx, alignment in enumerate(alignments):
        # Count this read as coverage for all MSA positions
        # Note: alignments span the full MSA
        for msa_pos in range(msa_length):
            error_matrix[msa_pos, 3] += 1  # coverage

        # Add errors at specific MSA positions using normalized errors
        # (excludes homopolymer extensions)
        for error_pos in alignment.normalized_error_positions:
            msa_pos = error_pos.msa_position
            if 0 <= msa_pos < msa_length:
                if error_pos.error_type == 'sub':
                    error_matrix[msa_pos, 0] += 1
                elif error_pos.error_type == 'ins':
                    error_matrix[msa_pos, 1] += 1
                elif error_pos.error_type == 'del':
                    error_matrix[msa_pos, 2] += 1

        # Extract base composition from aligned sequence with homopolymer normalization
        read_aligned = alignment.aligned_sequence
        if len(read_aligned) != msa_length:
            continue

        # Track what base each read has at each MSA position
        # Raw base composition plus separate tracking of homopolymer extensions
        for msa_pos in range(msa_length):
            read_base = read_aligned[msa_pos]

            # Track raw base composition
            if read_base in ['A', 'C', 'G', 'T', '-']:
                base_composition_matrix[msa_pos][read_base] += 1
            else:
                # Treat N or other ambiguous as gap
                base_composition_matrix[msa_pos]['-'] += 1

            # Additionally track if this is a homopolymer extension
            # NOTE: score_aligned is from the READ's perspective (seq1), which is what we want
            # since we're asking whether this particular READ base is an HP extension
            if alignment.score_aligned and msa_pos < len(alignment.score_aligned):
                if alignment.score_aligned[msa_pos] == '=':
                    # Homopolymer extension - track separately
                    if read_base in ['A', 'C', 'G', 'T', '-']:
                        homopolymer_composition_matrix[msa_pos][read_base] += 1
                    else:
                        homopolymer_composition_matrix[msa_pos]['-'] += 1

    # Calculate statistics for each MSA position
    position_stats = []

    for msa_pos in range(msa_length):
        sub_count = error_matrix[msa_pos, 0]
        ins_count = error_matrix[msa_pos, 1]
        del_count = error_matrix[msa_pos, 2]
        coverage = error_matrix[msa_pos, 3]

        # Total error events
        error_count = sub_count + ins_count + del_count
        error_rate = error_count / coverage if coverage > 0 else 0.0

        # Get consensus position (None for insertion columns)
        cons_pos = msa_to_consensus_pos[msa_pos]

        # Get consensus nucleotide at this MSA position
        cons_nucleotide = consensus_aligned[msa_pos]

        # Get base composition for this MSA position (raw counts)
        base_comp = base_composition_matrix[msa_pos].copy()

        # Get homopolymer extension composition for this MSA position
        hp_comp = homopolymer_composition_matrix[msa_pos].copy()

        position_stats.append(PositionStats(
            msa_position=msa_pos,
            consensus_position=cons_pos,
            coverage=coverage,
            error_count=error_count,
            error_rate=error_rate,
            sub_count=sub_count,
            ins_count=ins_count,
            del_count=del_count,
            consensus_nucleotide=cons_nucleotide,
            base_composition=base_comp,
            homopolymer_composition=hp_comp
        ))

    return position_stats


def is_variant_position_with_composition(
    position_stats: PositionStats,
    min_variant_frequency: float = 0.10,
    min_variant_count: int = 5
) -> Tuple[bool, List[str], str]:
    """
    Identify variant positions using simple frequency and count thresholds.

    This function determines if a position shows systematic variation (true biological
    variant) rather than scattered sequencing errors. Homopolymer extensions are
    excluded from consideration - diversity due solely to homopolymer length variation
    is not considered a true variant.

    Criteria for variant detection:
    1. At least one alternative allele must have frequency >= min_variant_frequency
    2. That allele must have count >= min_variant_count
    3. Counts are adjusted by subtracting homopolymer extension counts

    Args:
        position_stats: Position statistics including base composition
        min_variant_frequency: Minimum alternative allele frequency (default: 0.10 for 10%)
        min_variant_count: Minimum alternative allele read count (default: 5 reads)

    Returns:
        Tuple of (is_variant, variant_bases, reason)
        - is_variant: True if this position requires cluster separation
        - variant_bases: List of alternative bases meeting criteria (e.g., ['G', 'T'])
        - reason: Explanation of decision for logging/debugging
    """
    n = position_stats.coverage
    base_composition = position_stats.base_composition
    hp_composition = position_stats.homopolymer_composition

    # Check we have composition data
    if not base_composition or sum(base_composition.values()) == 0:
        return False, [], "No composition data available"

    # Calculate effective counts by subtracting homopolymer extensions
    # This excludes diversity that's purely due to HP length variation
    effective_composition = {}
    for base in ['A', 'C', 'G', 'T', '-']:
        raw_count = base_composition.get(base, 0)
        hp_count = hp_composition.get(base, 0) if hp_composition else 0
        effective_count = raw_count - hp_count
        if effective_count > 0:
            effective_composition[base] = effective_count

    # Check we have effective composition data after HP adjustment
    if not effective_composition or sum(effective_composition.values()) == 0:
        return False, [], "No composition data after HP adjustment"

    effective_total = sum(effective_composition.values())

    sorted_bases = sorted(
        effective_composition.items(),
        key=lambda x: x[1],
        reverse=True
    )

    if len(sorted_bases) < 2:
        return False, [], "No alternative alleles observed (after HP adjustment)"

    # Check each alternative allele (skip consensus base at index 0)
    variant_bases = []
    variant_details = []

    for base, count in sorted_bases[1:]:
        freq = count / effective_total if effective_total > 0 else 0

        # Must meet both frequency and count thresholds
        if freq >= min_variant_frequency and count >= min_variant_count:
            variant_bases.append(base)
            variant_details.append(f"{base}:{count}/{effective_total}({freq:.1%})")

    if variant_bases:
        return True, variant_bases, f"Variant alleles: {', '.join(variant_details)}"

    # Debug: Check if this would be a variant WITHOUT HP normalization
    # This helps identify cases where HP adjustment incorrectly eliminates variants
    raw_total = sum(base_composition.get(b, 0) for b in ['A', 'C', 'G', 'T', '-'])
    raw_sorted = sorted(
        [(b, base_composition.get(b, 0)) for b in ['A', 'C', 'G', 'T', '-'] if base_composition.get(b, 0) > 0],
        key=lambda x: x[1],
        reverse=True
    )
    if len(raw_sorted) >= 2:
        for base, count in raw_sorted[1:]:
            freq = count / raw_total if raw_total > 0 else 0
            if freq >= min_variant_frequency and count >= min_variant_count:
                # Would be variant without HP normalization!
                logging.debug(
                    f"HP normalization eliminated variant at MSA pos {position_stats.msa_position}: "
                    f"raw {base}:{count}/{raw_total}({freq:.1%}) meets threshold, "
                    f"but effective composition={effective_composition}, "
                    f"raw={base_composition}, hp={hp_composition}"
                )
                break

    return False, [], "No variants detected (after HP adjustment)"


def call_iupac_ambiguities(
    consensus: str,
    alignments: List['ReadAlignment'],
    msa_to_consensus_pos: Dict[int, Optional[int]],
    min_variant_frequency: float = 0.10,
    min_variant_count: int = 5
) -> Tuple[str, int, List[Dict]]:
    """
    Replace consensus bases at variant positions with IUPAC ambiguity codes.

    Analyzes positional variation in the MSA and identifies positions where
    significant variation remains after phasing. At these positions, the
    consensus base is replaced with the appropriate IUPAC code representing
    all variant alleles that meet the threshold criteria.

    Uses the same thresholds as phasing to ensure consistency. Homopolymer
    length variation is excluded (only true nucleotide variants are considered).

    Args:
        consensus: Ungapped consensus sequence from SPOA
        alignments: List of ReadAlignment objects from MSA
        msa_to_consensus_pos: Mapping from MSA position to consensus position
        min_variant_frequency: Minimum alternative allele frequency (default: 0.10)
        min_variant_count: Minimum alternative allele read count (default: 5)

    Returns:
        Tuple of:
        - Modified consensus sequence with IUPAC codes at variant positions
        - Count of IUPAC positions introduced
        - List of dicts with details about each IUPAC position:
          {
              'consensus_position': int,
              'original_base': str,
              'iupac_code': str,
              'variant_bases': List[str],
              'base_composition': Dict[str, int]
          }
    """
    if not consensus or not alignments:
        return consensus, 0, []

    # Reconstruct consensus_aligned from consensus and msa_to_consensus_pos
    # (same pattern as detect_variant_positions)
    msa_length = max(msa_to_consensus_pos.keys()) + 1 if msa_to_consensus_pos else 0
    if msa_length == 0:
        return consensus, 0, []

    consensus_aligned = []
    for msa_pos in range(msa_length):
        cons_pos = msa_to_consensus_pos.get(msa_pos)
        if cons_pos is not None and cons_pos < len(consensus):
            consensus_aligned.append(consensus[cons_pos])
        else:
            consensus_aligned.append('-')
    consensus_aligned_str = ''.join(consensus_aligned)

    # Analyze positional variation
    position_stats = analyze_positional_variation(alignments, consensus_aligned_str, msa_to_consensus_pos)

    # Build list of positions to replace
    iupac_positions = []

    for pos_stat in position_stats:
        # Skip insertion columns (no consensus position)
        if pos_stat.consensus_position is None:
            continue

        # Check if this position has significant variation
        is_variant, variant_bases, reason = is_variant_position_with_composition(
            pos_stat, min_variant_frequency, min_variant_count
        )

        if not is_variant:
            continue

        # Filter out gaps from variant bases (we can only represent nucleotide ambiguities)
        nucleotide_variants = [b for b in variant_bases if b in 'ACGT']

        if not nucleotide_variants:
            # Only gaps met the threshold - skip this position
            continue

        # Get the consensus base at this position
        cons_pos = pos_stat.consensus_position
        consensus_base = consensus[cons_pos] if cons_pos < len(consensus) else None

        if consensus_base is None or consensus_base not in 'ACGT':
            continue

        # Build set of all significant bases (consensus + variants)
        all_bases = set(nucleotide_variants)
        all_bases.add(consensus_base)

        # Look up IUPAC code
        iupac_code = IUPAC_CODES.get(frozenset(all_bases), 'N')

        # Only record if we actually need an ambiguity code (more than one base)
        if len(all_bases) > 1:
            iupac_positions.append({
                'consensus_position': cons_pos,
                'original_base': consensus_base,
                'iupac_code': iupac_code,
                'variant_bases': nucleotide_variants,
                'base_composition': pos_stat.base_composition
            })

    if not iupac_positions:
        return consensus, 0, []

    # Build modified consensus
    consensus_list = list(consensus)
    for pos_info in iupac_positions:
        cons_pos = pos_info['consensus_position']
        consensus_list[cons_pos] = pos_info['iupac_code']

    modified_consensus = ''.join(consensus_list)

    return modified_consensus, len(iupac_positions), iupac_positions


def calculate_within_cluster_error(
    haplotype_groups: Dict[str, Set[str]],
    read_alleles: Dict[str, Dict[int, str]],
    phasing_positions: Set[int],
    all_variant_positions: Set[int]
) -> float:
    """Calculate within-cluster error for a given haplotype grouping.

    Measures the average variation at ALL variant positions within each haplotype,
    including positions used for phasing. This ensures fair comparison across
    different candidate position sets and captures heterogeneity introduced by
    reassignment of non-qualifying haplotypes.

    Lower error indicates more homogeneous clusters.

    Args:
        haplotype_groups: Dict mapping allele_combo -> set of read_ids
        read_alleles: Dict mapping read_id -> {msa_position -> allele}
        phasing_positions: Set of MSA positions used for phasing (kept for API compatibility)
        all_variant_positions: Set of all variant MSA positions (error measured at all of these)

    Returns:
        Weighted average error rate across haplotypes (0.0 = perfect, 1.0 = maximum error)
    """
    # Measure error at ALL variant positions, not just non-phased ones.
    # This ensures fair comparison across candidate position sets and captures
    # heterogeneity introduced by reassignment at phasing positions.
    measured_positions = all_variant_positions

    if not measured_positions or not haplotype_groups:
        return 0.0

    total_weighted_error = 0.0
    total_reads = 0

    for combo, read_ids in haplotype_groups.items():
        if not read_ids:
            continue

        haplotype_error = 0.0
        positions_counted = 0

        for pos in measured_positions:
            # Count alleles at this position for reads in this haplotype
            allele_counts = defaultdict(int)
            for read_id in read_ids:
                allele = read_alleles.get(read_id, {}).get(pos, '-')
                allele_counts[allele] += 1

            if not allele_counts:
                continue

            # Find consensus (most common) allele
            total_at_pos = sum(allele_counts.values())
            max_count = max(allele_counts.values())

            # Error rate = fraction of reads NOT matching consensus
            error_at_pos = (total_at_pos - max_count) / total_at_pos
            haplotype_error += error_at_pos
            positions_counted += 1

        # Average error across all variant positions for this haplotype
        if positions_counted > 0:
            mean_haplotype_error = haplotype_error / positions_counted
            total_weighted_error += mean_haplotype_error * len(read_ids)
            total_reads += len(read_ids)

    if total_reads == 0:
        return 0.0

    return total_weighted_error / total_reads


def filter_qualifying_haplotypes(
    combo_to_reads: Dict[str, Set[str]],
    total_reads: int,
    min_count: int,
    min_frequency: float
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Filter haplotypes to those meeting count and frequency thresholds.

    Args:
        combo_to_reads: Dict mapping allele_combo -> set of read_ids
        total_reads: Total number of reads for frequency calculation
        min_count: Minimum read count threshold
        min_frequency: Minimum frequency threshold (0.0 to 1.0)

    Returns:
        Tuple of (qualifying_combos, non_qualifying_combos)
    """
    qualifying = {}
    non_qualifying = {}
    for combo, reads in combo_to_reads.items():
        count = len(reads)
        freq = count / total_reads if total_reads > 0 else 0
        if count >= min_count and freq >= min_frequency:
            qualifying[combo] = reads
        else:
            non_qualifying[combo] = reads
    return qualifying, non_qualifying


def group_reads_by_single_position(
    read_alleles: Dict[str, Dict[int, str]],
    position: int,
    read_ids: Set[str]
) -> Dict[str, Set[str]]:
    """Group a subset of reads by their allele at a single position.

    Args:
        read_alleles: Dict mapping read_id -> {msa_position -> allele}
        position: MSA position to group by
        read_ids: Subset of read IDs to consider

    Returns:
        Dict mapping allele -> set of read_ids
    """
    allele_to_reads = defaultdict(set)
    for read_id in read_ids:
        allele = read_alleles.get(read_id, {}).get(position, '-')
        allele_to_reads[allele].add(read_id)
    return dict(allele_to_reads)


