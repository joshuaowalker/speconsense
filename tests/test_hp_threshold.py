#!/usr/bin/env python3
"""Tests for the --hp-normalization-length threshold plumbing.

Covers both the distance-side (adjusted-identity factory) and the MSA-side
(`is_homopolymer_event`) halves of the summarize pipeline. Confirms MIN-of-runs
semantics: an HP length difference is treated as noise (normalized) only when
the shorter per-sequence run length is at or above the threshold; otherwise it
counts as a real edit / structural indel.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from speconsense.distances import build_adjustment_params, STANDARD_ADJUSTMENT_PARAMS
from speconsense.summarize import (
    analyze_msa_columns,
    calculate_adjusted_identity_distance,
    calculate_overlap_aware_distance,
    is_homopolymer_event,
    run_spoa_msa,
)


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def test_build_adjustment_params_default_matches_standard():
    params = build_adjustment_params()
    assert params.normalize_homopolymers is True
    assert params.max_repeat_motif_length == 1
    assert params.hp_normalize_min_length == 1
    assert params.hp_normalize_min_length == STANDARD_ADJUSTMENT_PARAMS.hp_normalize_min_length


def test_build_adjustment_params_threshold_propagates():
    params = build_adjustment_params(hp_normalization_length=6)
    assert params.hp_normalize_min_length == 6
    # Other settings unchanged
    assert params.normalize_homopolymers is True
    assert params.max_repeat_motif_length == 1


# ---------------------------------------------------------------------------
# Distance side
# ---------------------------------------------------------------------------


def _flank_left(n: int) -> str:
    return "GCGC"[:max(0, n)] + "GCGCGCGCGC"[: max(0, 6 - n)]


# Use stable, non-repetitive flanking context so the HP run is the only
# disagreement.
LEFT = "ATGCATGCATGC"
RIGHT = "TTACGTACGTAC"


def _hp_seq(n: int, base: str = "A") -> str:
    return LEFT + base * n + RIGHT


def test_calculate_adjusted_identity_distance_default_normalizes_short_hp():
    # AAA vs AAAA collapses to identity 1.0 under the legacy default (1).
    d = calculate_adjusted_identity_distance(_hp_seq(3), _hp_seq(4))
    assert d == 0.0


def test_calculate_adjusted_identity_distance_threshold_demotes_short_hp():
    # min(L1, L2) = 3 < 6 => treated as edit, distance > 0
    d = calculate_adjusted_identity_distance(
        _hp_seq(3), _hp_seq(4), hp_normalization_length=6
    )
    assert d > 0.0


def test_calculate_adjusted_identity_distance_threshold_keeps_long_hp_normalized():
    # min(L1, L2) = 6 >= 6 => still normalized
    d = calculate_adjusted_identity_distance(
        _hp_seq(6), _hp_seq(8), hp_normalization_length=6
    )
    assert d == 0.0


def test_calculate_adjusted_identity_distance_boundary_straddle_demotes():
    # min(5, 7) = 5 < 6 => demoted (matches MIN-of-runs semantics)
    d = calculate_adjusted_identity_distance(
        _hp_seq(5), _hp_seq(7), hp_normalization_length=6
    )
    assert d > 0.0


def test_calculate_overlap_aware_distance_threshold_demotes_short_hp():
    s1 = _hp_seq(3)
    s2 = _hp_seq(4)
    d_default = calculate_overlap_aware_distance(s1, s2, min_overlap_bp=10)
    d_strict = calculate_overlap_aware_distance(
        s1, s2, min_overlap_bp=10, hp_normalization_length=6
    )
    assert d_default == 0.0
    assert d_strict > d_default


# ---------------------------------------------------------------------------
# MSA side: is_homopolymer_event threshold
# ---------------------------------------------------------------------------


def _msa_with_hp(short_n: int, long_n: int, base: str = "A"):
    """Build a 2-sequence MSA where one row has ``short_n`` of ``base`` and
    the other has ``long_n``, both flanked by stable non-HP context. Returns
    (aligned_seqs, start_col, end_col) for the indel span.
    """
    seqs = [_hp_seq(short_n, base), _hp_seq(long_n, base)]
    aligned = run_spoa_msa(seqs, alignment_mode=1)
    # Find the indel event span in the aligned strings
    s1 = str(aligned[0].seq)
    s2 = str(aligned[1].seq)
    span = None
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if c1 == '-' or c2 == '-':
            if span is None:
                span = [i, i]
            else:
                span[1] = i
        elif span is not None:
            break
    assert span is not None, "expected an indel column from MSA"
    return aligned, span[0], span[1]


def test_is_homopolymer_event_default_classifies_short_run_as_hp():
    aligned, s, e = _msa_with_hp(3, 4)
    assert is_homopolymer_event(aligned, s, e) is True


def test_is_homopolymer_event_threshold_demotes_short_run():
    aligned, s, e = _msa_with_hp(3, 4)
    assert is_homopolymer_event(aligned, s, e, min_hp_length=6) is False


def test_is_homopolymer_event_threshold_keeps_long_run_as_hp():
    aligned, s, e = _msa_with_hp(8, 10)
    assert is_homopolymer_event(aligned, s, e, min_hp_length=6) is True


def test_is_homopolymer_event_threshold_boundary_min_below_demotes():
    # min(L1, L2) = 5 < 6 => structural
    aligned, s, e = _msa_with_hp(5, 7)
    assert is_homopolymer_event(aligned, s, e, min_hp_length=6) is False


def test_is_homopolymer_event_threshold_boundary_min_inclusive():
    # min(L1, L2) = 6 >= 6 => homopolymer
    aligned, s, e = _msa_with_hp(6, 7)
    assert is_homopolymer_event(aligned, s, e, min_hp_length=6) is True


# ---------------------------------------------------------------------------
# MSA side: analyze_msa_columns flows the threshold through
# ---------------------------------------------------------------------------


def test_analyze_msa_columns_default_classifies_as_homopolymer():
    aligned, _, _ = _msa_with_hp(3, 4)
    stats = analyze_msa_columns(aligned)
    assert stats['homopolymer_indel_count'] == 1
    assert stats['structural_indel_count'] == 0


def test_analyze_msa_columns_threshold_demotes_short_run():
    aligned, _, _ = _msa_with_hp(3, 4)
    stats = analyze_msa_columns(aligned, min_hp_length=6)
    assert stats['homopolymer_indel_count'] == 0
    assert stats['structural_indel_count'] == 1


def test_analyze_msa_columns_threshold_passes_long_run():
    aligned, _, _ = _msa_with_hp(8, 10)
    stats = analyze_msa_columns(aligned, min_hp_length=6)
    assert stats['homopolymer_indel_count'] == 1
    assert stats['structural_indel_count'] == 0
