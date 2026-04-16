"""Tests for the context-aware CER functions in speconsense.significance."""

import math

import pytest

from speconsense.significance import (
    compute_cer_factor,
    compute_critical_error_rate,
    compute_joint_critical_q,
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
