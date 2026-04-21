"""Tests for the context-aware CER functions in speconsense.significance."""

import math

import pytest

from speconsense.significance import (
    compute_cer_factor,
    compute_critical_error_rate,
    compute_joint_critical_q,
    compute_log_joint_critical_q,
    compute_per_position_qstar,
)


# ---------------------------------------------------------------------------
# compute_joint_critical_q
# ---------------------------------------------------------------------------


def test_joint_q_typical_case():
    # N=1000, M=100, n_sites=700, K=1, alpha=1e-5
    # Should give a small but positive q*
    q = compute_joint_critical_q(N=1000, M=100, n_sites=700, K=1, alpha=1e-5)
    assert q is not None
    assert 0 < q < 1
    # Sanity: at this many supporting reads, q* should be a few percent
    assert 0.01 < q < 0.10


def test_joint_q_consistent_with_uniform_compute_critical_error_rate():
    # compute_critical_error_rate returns p* such that q = p*/3 satisfies the
    # binomial equation. compute_joint_critical_q at K=1 with the same n_sites
    # should give that same q.
    p_star = compute_critical_error_rate(N=1000, M=100, L=700, alpha=1e-5, K=1)
    q_uniform = p_star / 3.0
    q_direct = compute_joint_critical_q(N=1000, M=100, n_sites=700, K=1, alpha=1e-5)
    assert q_direct is not None
    assert math.isclose(q_uniform, q_direct, rel_tol=1e-4)


def test_joint_q_returns_none_for_invalid_inputs():
    assert compute_joint_critical_q(N=0, M=10, n_sites=700) is None
    assert compute_joint_critical_q(N=100, M=0, n_sites=700) is None
    assert compute_joint_critical_q(N=100, M=200, n_sites=700) is None    # M > N
    assert compute_joint_critical_q(N=100, M=10, n_sites=0) is None
    assert compute_joint_critical_q(N=100, M=10, n_sites=700, K=701) is None


def test_joint_q_increases_with_more_M():
    # Higher M means the variant is more strongly supported, so the artifact
    # explanation requires a HIGHER per-read error probability — q* is larger.
    q_low_M = compute_joint_critical_q(N=1000, M=20, n_sites=700, alpha=1e-5)
    q_high_M = compute_joint_critical_q(N=1000, M=100, n_sites=700, alpha=1e-5)
    assert q_high_M > q_low_M


# ---------------------------------------------------------------------------
# compute_per_position_qstar
# ---------------------------------------------------------------------------


def test_per_position_qstar_k1_equals_joint():
    joint = compute_joint_critical_q(N=1000, M=50, n_sites=700, K=1, alpha=1e-5)
    per_pos = compute_per_position_qstar(N=1000, M=50, n_sites=700, K=1, alpha=1e-5)
    assert math.isclose(per_pos, joint)


def test_per_position_qstar_k2_takes_root():
    joint = compute_joint_critical_q(N=1000, M=10, n_sites=700, K=2, alpha=1e-5)
    per_pos = compute_per_position_qstar(N=1000, M=10, n_sites=700, K=2, alpha=1e-5)
    assert math.isclose(per_pos, joint ** 0.5, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# compute_cer_factor
# ---------------------------------------------------------------------------


def test_factor_k1_basic():
    # K=1, q_ctx=0.0059 (non-HP sub), strong variant
    factor = compute_cer_factor(
        N=1000, M=100, n_sites=700,
        q_ctx_per_position=[0.0059], alpha=1e-5,
    )
    assert factor is not None
    # At M=100/N=1000 the per-position q* is around 5-6%, vs q_ctx=0.59% -> factor ~10
    assert 5.0 < factor < 20.0


def test_factor_equals_qstar_over_qctx_at_k1():
    q_ctx = 0.0074
    qstar = compute_per_position_qstar(N=1000, M=50, n_sites=700, K=1, alpha=1e-5)
    factor = compute_cer_factor(
        N=1000, M=50, n_sites=700,
        q_ctx_per_position=[q_ctx], alpha=1e-5,
    )
    assert math.isclose(factor, qstar / q_ctx, rel_tol=1e-9)


def test_factor_k2_uniform_qctx():
    # At K=2 with uniform q_ctx, factor should still be in per-position units.
    q_ctx = 0.0059
    factor = compute_cer_factor(
        N=1000, M=10, n_sites=700,
        q_ctx_per_position=[q_ctx, q_ctx], alpha=1e-5,
    )
    assert factor is not None
    # Sanity: M=10 at K=2 is significant per the RFC (min M=6 at K=2).
    assert factor > 1.0


def test_factor_k2_mixed_contexts():
    # Mixed contexts: one non-HP sub and one HP-l3 length change.
    q_non_hp = 0.0059
    q_hp_l3 = 0.0097
    factor_mixed = compute_cer_factor(
        N=1000, M=10, n_sites=700,
        q_ctx_per_position=[q_non_hp, q_hp_l3], alpha=1e-5,
    )
    factor_uniform = compute_cer_factor(
        N=1000, M=10, n_sites=700,
        q_ctx_per_position=[q_non_hp, q_non_hp], alpha=1e-5,
    )
    assert factor_mixed is not None
    assert factor_uniform is not None
    # Higher denominator (q_hp_l3) -> smaller factor
    assert factor_mixed < factor_uniform


def test_factor_invalid_inputs():
    assert compute_cer_factor(N=1000, M=50, n_sites=700, q_ctx_per_position=[]) is None
    assert compute_cer_factor(N=1000, M=50, n_sites=700, q_ctx_per_position=[0.0]) is None
    assert compute_cer_factor(N=1000, M=50, n_sites=700, q_ctx_per_position=[None]) is None
    assert compute_cer_factor(N=1000, M=50, n_sites=700, q_ctx_per_position=[1.0]) is None


def test_factor_decreases_under_inflated_qctx():
    # Same N/M but a higher q_ctx (e.g., post-update from new chemistry table)
    # should produce a smaller factor.
    factor_low_q = compute_cer_factor(
        N=1000, M=50, n_sites=700,
        q_ctx_per_position=[0.0059], alpha=1e-5,
    )
    factor_high_q = compute_cer_factor(
        N=1000, M=50, n_sites=700,
        q_ctx_per_position=[0.015], alpha=1e-5,
    )
    assert factor_high_q < factor_low_q


def test_factor_rescaling_works():
    # If q_ctx is scaled by a factor r, the resulting CER factor scales by 1/r.
    # This is the property that makes stored output reinterpretable when the
    # q_ctx table is updated.
    base_q = 0.0059
    r = 1.5
    f1 = compute_cer_factor(N=1000, M=50, n_sites=700,
                            q_ctx_per_position=[base_q], alpha=1e-5)
    f2 = compute_cer_factor(N=1000, M=50, n_sites=700,
                            q_ctx_per_position=[base_q * r], alpha=1e-5)
    assert math.isclose(f2, f1 / r, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# High-K regression tests — log-space brentq fix
#
# Pre-fix: linear-space brentq with eps=1e-15 returned 0.0 whenever the true
# q* fell below the lower search bound (typical for K >= ~16). compute_cer_factor
# then wrote 0.000 to the FASTA header, and --min-cer-factor 1.0 routed the
# record to .ns. These tests pin the behavior end-to-end for the v21 case
# observed in /Users/josh/mm/data/ont98/sigtest/0420-v3 (ONT01.01-A01 v21).
# ---------------------------------------------------------------------------


def test_joint_q_high_k_returns_subepsilon_value():
    # ONT01.01-A01 v21: K=19, M=3, N=970, n_sites=759.
    # Hand-derived q* via the small-q asymptotic (M·log q ≈ log α − log_corr − log C(N,M))
    # is ~1.24e-17, well below the old 1e-15 lower bound.
    q = compute_joint_critical_q(N=970, M=3, n_sites=759, K=19, alpha=1e-5)
    assert q is not None
    assert 1e-18 < q < 1e-16


def test_log_joint_critical_q_returns_negative_log():
    # Same v21 case, exposed via the log-space wrapper for callers that need
    # to avoid exp/log round-trip precision loss.
    log_q = compute_log_joint_critical_q(N=970, M=3, n_sites=759, K=19, alpha=1e-5)
    assert log_q is not None
    # log(1e-18) ≈ -41.4; log(1e-16) ≈ -36.8. v21 sits near -38.9.
    assert -41.5 < log_q < -36.8


def test_cer_factor_high_k_v21_regression():
    # Exact q_ctx values from v21 metadata: 14 non-hp-sub (0.0059), 2 non-hp-indel
    # (0.0108), one each of hp-l5-G-del1 (0.0113), hp-l2-T-del2 (0.0083),
    # hp-l4-T-del1 (0.0099). K=19.
    q_ctx = ([0.0059] * 14 + [0.0108, 0.0108]
             + [0.0113, 0.0083, 0.0099])
    factor = compute_cer_factor(N=970, M=3, n_sites=759,
                                q_ctx_per_position=q_ctx, alpha=1e-5)
    assert factor is not None
    # Hand-derived target ≈ 18.9; allow ±5% for solver tolerance.
    assert 17.5 < factor < 20.5


def test_cer_factor_extreme_k_no_underflow():
    # K=100, q_ctx=0.008 each — actual_joint = 0.008^100 ≈ 1.6e-210, near the
    # float64 underflow boundary. Pre-fix the linear-space product was OK at
    # this K but log-space removes the hazard entirely.
    q_ctx = [0.008] * 100
    factor = compute_cer_factor(N=1000, M=3, n_sites=700,
                                q_ctx_per_position=q_ctx, alpha=1e-5)
    assert factor is not None
    assert factor > 1.0   # 100 distinct mutations in 3 reads → highly significant


def test_compute_critical_error_rate_high_k_returns_finite():
    # Parallel fix in the uniform-error variant. Pre-fix this returned 0.0 for
    # the same reason as the joint solver. We expect a small but positive p*.
    p_star = compute_critical_error_rate(N=970, M=3, L=759, alpha=1e-5, K=19)
    assert p_star > 0.0
    assert p_star < 1.0   # well below the signal destruction threshold


def test_factor_rescaling_holds_at_high_k():
    # The K-th-root scaling property must survive the log-space refactor.
    base = [0.006] * 19
    scaled = [q * 1.5 for q in base]
    f1 = compute_cer_factor(N=970, M=3, n_sites=759,
                            q_ctx_per_position=base, alpha=1e-5)
    f2 = compute_cer_factor(N=970, M=3, n_sites=759,
                            q_ctx_per_position=scaled, alpha=1e-5)
    assert math.isclose(f2, f1 / 1.5, rel_tol=1e-6)
