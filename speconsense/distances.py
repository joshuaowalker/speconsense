"""Shared distance and difference calculations for sequence comparison.

Provides IUPAC-aware alignment utilities, variant difference counting,
and adjusted identity distance calculations used by both core and
summarize subpackages.
"""

from typing import Tuple

import edlib
from adjusted_identity import score_alignment, AdjustmentParams

from speconsense.msa import IUPAC_CODES


# IUPAC equivalencies for edlib alignment
# This allows edlib to treat IUPAC ambiguity codes as matching their constituent bases
IUPAC_EQUIV = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"),
               ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
               ("W", "A"), ("W", "T"), ("M", "A"), ("M", "C"),
               ("S", "C"), ("S", "G"), ("K", "G"), ("K", "T"),
               ("B", "C"), ("B", "G"), ("B", "T"),
               ("D", "A"), ("D", "G"), ("D", "T"),
               ("H", "A"), ("H", "C"), ("H", "T"),
               ("V", "A"), ("V", "C"), ("V", "G"), ]

# Standard adjustment parameters for consistent sequence comparison
# Used by both substitution distance calculation and adjusted identity distance
STANDARD_ADJUSTMENT_PARAMS = AdjustmentParams(
    normalize_homopolymers=True,    # Enable homopolymer normalization
    handle_iupac_overlap=False,     # Disable IUPAC overlap - use standard IUPAC semantics (Y!=M)
    normalize_indels=False,         # Disable indel normalization
    end_skip_distance=0,            # No end trimming - sequences must match end-to-end
    max_repeat_motif_length=1       # Single-base repeats for homopolymer normalization
)

# Safeguards for unreliable adjusted identity on length-mismatched sequences.
# When terminal gap exclusion inflates identity, fall back to raw edlib identity.
MIN_COVERAGE_THRESHOLD = 0.5   # min(seq1_coverage, seq2_coverage) floor
MAX_ADJUSTMENT_RATIO = 1.5     # max adjusted_identity / raw_identity


def expand_iupac_code(base: str) -> set:
    """
    Expand an IUPAC code to its constituent nucleotides.
    Returns a set of nucleotides that the code represents.
    """
    iupac_expansion = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'},
    }

    return iupac_expansion.get(base.upper(), {'N'})


def bases_match_with_iupac(base1: str, base2: str) -> bool:
    """
    Check if two bases match, considering IUPAC ambiguity codes.
    Two bases match if their IUPAC expansions have any nucleotides in common.
    """
    if base1 == base2:
        return True

    # Handle gap characters
    if base1 == '-' or base2 == '-':
        return base1 == base2

    # Expand IUPAC codes and check for overlap
    expansion1 = expand_iupac_code(base1)
    expansion2 = expand_iupac_code(base2)

    # Bases match if their expansions have any nucleotides in common
    return bool(expansion1.intersection(expansion2))


def _count_alignment_differences(query_aligned: str, target_aligned: str) -> Tuple[int, int, int, int]:
    """Count differences in an IUPAC-aware alignment.

    Returns (substitutions, single_nt_indels, short_indels, long_indels)
    where each contiguous indel counts as one event regardless of length.
    """
    substitutions = 0
    single_nt_indels = 0  # Single nucleotide indels
    short_indels = 0      # 2-3 nt indels
    long_indels = 0       # 4+ nt indels

    i = 0
    while i < len(query_aligned):
        query_char = query_aligned[i]
        target_char = target_aligned[i]

        if not bases_match_with_iupac(query_char, target_char):
            if query_char == '-' or target_char == '-':
                indel_length = 1
                j = i + 1
                while j < len(query_aligned) and (query_aligned[j] == '-' or target_aligned[j] == '-'):
                    indel_length += 1
                    j += 1

                if indel_length == 1:
                    single_nt_indels += 1
                elif indel_length <= 3:
                    short_indels += 1
                else:
                    long_indels += 1

                i = j
                continue
            else:
                substitutions += 1

        i += 1

    return substitutions, single_nt_indels, short_indels, long_indels


def count_variant_differences(seq1: str, seq2: str) -> int:
    """Count the number of variant-level differences between two sequences.

    Uses IUPAC-aware alignment to count substitutions and indel events
    (each contiguous indel counts as one difference regardless of length).
    This count is used as K in multi-position CER evaluation.

    Returns 0 if sequences are identical or IUPAC-compatible.
    Returns -1 if alignment fails.
    """
    if not seq1 or not seq2:
        return -1

    if seq1 == seq2:
        return 0

    try:
        result = edlib.align(seq1, seq2, task="path", additionalEqualities=IUPAC_EQUIV)
        if result["editDistance"] == -1:
            return -1

        alignment = edlib.getNiceAlignment(result, seq1, seq2)
        if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
            return -1

        subs, single_indels, short_indels, long_indels = _count_alignment_differences(
            alignment['query_aligned'], alignment['target_aligned']
        )
        return subs + single_indels + short_indels + long_indels

    except Exception:
        return -1


def calculate_adjusted_identity_distance(seq1: str, seq2: str) -> float:
    """Calculate adjusted identity distance between two sequences.

    Uses homopolymer-normalized adjusted identity with safeguards against
    unreliable results from length-mismatched sequences.

    Returns distance from 0.0 (identical) to 1.0 (maximally different).
    """
    if not seq1 or not seq2:
        return 1.0  # Maximum distance

    if seq1 == seq2:
        return 0.0

    # Get alignment from edlib with IUPAC awareness
    result = edlib.align(seq1, seq2, task="path", additionalEqualities=IUPAC_EQUIV)
    if result["editDistance"] == -1:
        return 1.0

    # Raw identity from edlib (global, no adjustments)
    raw_identity = 1.0 - (result["editDistance"] / max(len(seq1), len(seq2)))

    # Get nice alignment for adjusted identity scoring
    alignment = edlib.getNiceAlignment(result, seq1, seq2)
    if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
        return 1.0

    # Calculate adjusted identity using standard parameters
    score_result = score_alignment(
        alignment['query_aligned'],
        alignment['target_aligned'],
        adjustment_params=STANDARD_ADJUSTMENT_PARAMS
    )
    adjusted_identity = score_result.identity

    # Guard: coverage too low → adjusted identity unreliable
    # (e.g. 220bp vs 660bp: terminal gap exclusion scores only the overlap)
    min_coverage = min(score_result.seq1_coverage, score_result.seq2_coverage)
    if min_coverage < MIN_COVERAGE_THRESHOLD:
        return 1.0 - raw_identity

    # Guard: adjustment inflated identity beyond reasonable bounds
    # (legitimate HP normalization raises identity by ~1-5%, not 3x)
    if raw_identity > 0 and adjusted_identity > raw_identity * MAX_ADJUSTMENT_RATIO:
        return 1.0 - raw_identity

    # Convert adjusted identity to distance
    return 1.0 - adjusted_identity
