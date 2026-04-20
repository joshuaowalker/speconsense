"""IUPAC ambiguity code handling and distance calculations.

Provides utilities for working with IUPAC nucleotide codes and calculating
adjusted identity distances between sequences with homopolymer normalization.
"""

from typing import List, Optional

import edlib
from adjusted_identity import align_and_score

from speconsense.msa import IUPAC_CODES

# Re-export shared distance functions and constants from speconsense.distances
from speconsense.distances import (  # noqa: F401
    IUPAC_EQUIV,
    STANDARD_ADJUSTMENT_PARAMS,
    MIN_COVERAGE_THRESHOLD,
    MAX_ADJUSTMENT_RATIO,
    expand_iupac_code,
    bases_match_with_iupac,
    _count_alignment_differences,
    count_variant_differences,
    calculate_adjusted_identity_distance,
)


def primers_are_same(p1: Optional[List[str]], p2: Optional[List[str]]) -> bool:
    """Check if two primer annotations indicate the same amplicon.

    Used to determine whether overlap-aware merging should be allowed.
    When primers match, sequences should have the same amplicon length,
    so length differences indicate chimeras rather than primer pool variation.

    Returns True (use global distance, no overlap merging) when:
    - Either is None or empty (conservative: unknown = assume same)
    - Both have identical primer sets

    Returns False (allow overlap merging) when primers differ.

    Args:
        p1: Primer annotation from first sequence (e.g., ['ITS1', 'ITS4'])
        p2: Primer annotation from second sequence

    Returns:
        True if primers are same or unknown (use global distance)
        False if primers differ (allow overlap-aware distance)
    """
    if not p1 or not p2:
        return True  # Conservative: missing info -> assume same
    return set(p1) == set(p2)


def merge_bases_to_iupac(bases: set) -> str:
    """
    Merge a set of bases (which may include IUPAC codes) into a single IUPAC code.

    Expands any existing IUPAC codes to their constituent nucleotides,
    takes the union, and returns the appropriate IUPAC code.

    Examples:
        {'C', 'Y'} -> 'Y'  (Y=CT, so C+Y = CT = Y)
        {'A', 'R'} -> 'R'  (R=AG, so A+R = AG = R)
        {'C', 'R'} -> 'V'  (R=AG, so C+R = ACG = V)
    """
    # Expand all bases to their constituent nucleotides
    all_nucleotides = set()
    for base in bases:
        all_nucleotides.update(expand_iupac_code(base))

    # Look up the IUPAC code for the combined set
    return IUPAC_CODES.get(frozenset(all_nucleotides), 'N')


def create_variant_summary(primary_seq: str, variant_seq: str) -> str:
    """
    Compare a variant sequence to the primary sequence and create a summary string
    describing the differences. Returns a summary like:
    "3 substitutions, 1 single-nt indel, 1 short (<= 3nt) indel, 2 long indels"
    """
    if not primary_seq or not variant_seq:
        return "sequences empty - cannot compare"

    if primary_seq == variant_seq:
        return "identical sequences"

    try:
        result = edlib.align(primary_seq, variant_seq, task="path", additionalEqualities=IUPAC_EQUIV)
        if result["editDistance"] == -1:
            return "alignment failed"

        alignment = edlib.getNiceAlignment(result, primary_seq, variant_seq)
        if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
            return f"alignment parsing failed - edit distance {result['editDistance']}"

        substitutions, single_nt_indels, short_indels, long_indels = _count_alignment_differences(
            alignment['query_aligned'], alignment['target_aligned']
        )

        parts = []
        if substitutions > 0:
            parts.append(f"{substitutions} substitution{'s' if substitutions != 1 else ''}")
        if single_nt_indels > 0:
            parts.append(f"{single_nt_indels} single-nt indel{'s' if single_nt_indels != 1 else ''}")
        if short_indels > 0:
            parts.append(f"{short_indels} short (<= 3nt) indel{'s' if short_indels != 1 else ''}")
        if long_indels > 0:
            parts.append(f"{long_indels} long indel{'s' if long_indels != 1 else ''}")

        if not parts:
            return "identical sequences (IUPAC-compatible)"

        return ", ".join(parts)

    except Exception as e:
        return f"comparison failed: {str(e)}"


def calculate_overlap_aware_distance(seq1: str, seq2: str, min_overlap_bp: int) -> float:
    """
    Calculate distance that accounts for partial overlaps between sequences.

    When sequences have sufficient overlap with good identity, returns the
    overlap-region distance. Otherwise falls back to global distance.

    For containment cases where one sequence is shorter than min_overlap_bp,
    uses the shorter sequence length as the effective threshold.

    Args:
        seq1, seq2: DNA sequences (may have different lengths)
        min_overlap_bp: Minimum overlap required in base pairs

    Returns:
        Distance (0.0 to 1.0) based on overlap region if sufficient,
        otherwise global distance from calculate_adjusted_identity_distance()

    Note: This function calculates distance purely based on sequence content.
    Primer-based filtering (to prevent chimera grouping) is applied at the
    caller level via ``primers_are_same`` — see
    ``merge_groups_by_anchor_overlap``.
    """
    if not seq1 or not seq2:
        return 1.0  # Maximum distance

    if seq1 == seq2:
        return 0.0

    # Use align_and_score which handles bidirectional alignment internally
    result = align_and_score(seq1, seq2, STANDARD_ADJUSTMENT_PARAMS)

    # Calculate overlap in base pairs
    # Coverage is fraction of each sequence used in alignment
    len1, len2 = len(seq1), len(seq2)
    shorter_len = min(len1, len2)

    # Overlap is the minimum of the two coverages times the respective lengths
    # For containment, the shorter sequence should be fully covered
    overlap_bp = int(min(result.seq1_coverage * len1, result.seq2_coverage * len2))

    # Effective threshold: for containment cases, allow merge if short sequence is fully covered
    effective_threshold = min(min_overlap_bp, shorter_len)

    if overlap_bp >= effective_threshold:
        # Sufficient overlap - use overlap identity for distance
        return 1.0 - result.identity
    else:
        # Insufficient overlap - fall back to global distance
        # This will typically be high due to terminal gaps
        return calculate_adjusted_identity_distance(seq1, seq2)
