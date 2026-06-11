"""Statistical significance testing for variant phasing.

Implements the critical error rate (p*) framework from the variant significance
paper. Uses the binomial survival function with Bonferroni correction to
determine whether an observed variant count M out of N total reads at L positions
is statistically distinguishable from sequencing error.

Error model: uniform (q = p/3), where p is total per-position error rate.

Numerical note: brentq solvers operate in log-space (parameter ``log_q`` or
``log_p``) so q* values below ~1e-15 — common at high K — are recoverable
rather than silently clipped to zero. When ``scipy.stats.binom.logsf``
underflows to ``-inf`` (q below ~1e-200), the closed-form asymptotic
``log P(>=M) ≈ log C(N,M) + M·log(q)`` takes over; the two formulas agree
to many digits across the overlap region, so the brentq objective remains
continuous.
"""

import math
from typing import Optional, Sequence

from scipy.stats import binom
from scipy.optimize import brentq

# Log-space brentq bounds. exp(-700) ≈ 1e-304, just above float64 underflow;
# the upper bound stops just below log(1) = 0.
_LOG_Q_MIN = -700.0
_LOG_Q_MAX = -1e-15
_LOG3 = math.log(3.0)

# Default boundary for the per-position CER factor below which summarize routes a
# variant to the ``.ns`` track. cer_factor < 1.0 means the pairwise difference is
# within single-position significance (i.e. plausibly a basecaller artifact of a
# larger peer). This is the natural floor, not merely a UI default — core's vid
# tiering (SpecimenClusterer._assign_identity_ranks) and summarize's
# ``--min-cer-factor`` default both reference it so they cannot drift.
DEFAULT_MIN_CER_FACTOR = 1.0


def _log_binom_coef(N: int, M: int) -> float:
    """log C(N, M) via lgamma (handles large N robustly)."""
    return (math.lgamma(N + 1) - math.lgamma(M + 1)
            - math.lgamma(N - M + 1))


def _log_p_at_least(M: int, N: int, log_q: float, log_binom_NM: float) -> float:
    """log P(Binom(N, q) >= M) at q = exp(log_q).

    Falls back to the small-q asymptotic ``log C(N,M) + M·log(q)`` when
    ``binom.logsf`` underflows. Asymptotic is exact in the limit Nq → 0
    and we only invoke it where the exact form has lost all precision.
    """
    if log_q <= _LOG_Q_MIN:
        q = 0.0
    else:
        q = math.exp(log_q)
    if q > 0.0:
        logsf = binom.logsf(M - 1, N, q)
        if math.isfinite(logsf):
            return logsf
    return log_binom_NM + M * log_q


def _solve_log_qstar(N: int, M: int, n_sites: int, K: int,
                     alpha: float) -> Optional[float]:
    """Solve C(n_sites, K) · P(Binom(N, q*) >= M) = alpha for log(q*).

    Returns:
        log(q*) in (_LOG_Q_MIN, _LOG_Q_MAX), or
        ``None`` if the equation has no solution for any q < 1
        (M too large relative to N), or
        ``float('-inf')`` if the test passes even at q ≈ 0
        (essentially unreachable with the asymptotic dominating).
    """
    if M <= 0 or N <= 0 or n_sites <= 0 or K <= 0:
        return None
    if M > N or K > n_sites:
        return None

    log_alpha = math.log(alpha)
    log_correction = math.log(math.comb(n_sites, K))
    log_binom_NM = _log_binom_coef(N, M)

    def obj(log_q: float) -> float:
        return (log_correction
                + _log_p_at_least(M, N, log_q, log_binom_NM)
                - log_alpha)

    try:
        val_low = obj(_LOG_Q_MIN)
        val_high = obj(_LOG_Q_MAX)
    except (ValueError, OverflowError):
        return None

    if val_high < 0:
        return None    # M too large: even q≈1 doesn't satisfy
    if val_low > 0:
        return float('-inf')   # essentially unreachable with the asymptotic

    try:
        return brentq(obj, _LOG_Q_MIN, _LOG_Q_MAX, xtol=1e-10)
    except (ValueError, RuntimeError):
        return None


def _solve_log_pstar_uniform(N: int, M: int, L: int, K: int,
                             alpha: float) -> Optional[float]:
    """Solve C(L,K) · P(Binom(N, (p*/3)^K) >= M) = alpha for log(p*).

    Uniform error model variant of ``_solve_log_qstar``. Returns log(p*) where
    p* is the total per-position error rate (q = p/3 per alternative base).
    Search domain is log(p) ∈ (_LOG_Q_MIN + log 3, log 3); the upper bound is
    log(3) since q = p/3 must remain ≤ 1.
    """
    if M <= 0 or N <= 0 or L <= 0 or K <= 0:
        return None
    if M > N or K > L:
        return None

    log_alpha = math.log(alpha)
    log_correction = math.log(math.comb(L, K))
    log_binom_NM = _log_binom_coef(N, M)

    log_p_min = _LOG_Q_MIN + _LOG3
    log_p_max = _LOG3 - 1e-15   # p just below 3, so q = p/3 just below 1

    def obj(log_p: float) -> float:
        # q_joint = (p/3)^K  →  log_q_joint = K · (log_p - log 3)
        log_q_joint = K * (log_p - _LOG3)
        return (log_correction
                + _log_p_at_least(M, N, log_q_joint, log_binom_NM)
                - log_alpha)

    try:
        val_low = obj(log_p_min)
        val_high = obj(log_p_max)
    except (ValueError, OverflowError):
        return None

    if val_high < 0:
        return None
    if val_low > 0:
        return float('-inf')

    try:
        return brentq(obj, log_p_min, log_p_max, xtol=1e-10)
    except (ValueError, RuntimeError):
        return None


def compute_critical_error_rate(N: int, M: int, L: int, alpha: float = 1e-5,
                                K: int = 1) -> float:
    """Compute p* under uniform error model (q = p/3).

    For K=1: L * P(Binom(N, p*/3) >= M) = alpha
    For K>1: C(L,K) * P(Binom(N, (p*/3)^K) >= M) = alpha

    The K>1 model requires the SAME M reads to error at ALL K positions.
    Under independent errors across positions, the per-read probability of
    matching the variant at all K positions is q^K.

    Args:
        N: Total specimen reads (denominator)
        M: Variant read count
        L: Amplicon length (number of positions for Bonferroni correction)
        alpha: Significance level
        K: Number of variant positions (default 1)

    Returns:
        p* as total error rate (uncapped), or 0.0 if M <= 0
    """
    log_pstar = _solve_log_pstar_uniform(N, M, L, K, alpha)
    if log_pstar is None:
        # Bug-for-bug compat with prior contract: None paths returned
        # float('inf') (val_high < 0) or 0.0 (everything else).
        # Distinguish using a quick endpoint probe.
        if M <= 0 or N <= 0 or L <= 0 or K <= 0 or M > N or K > L:
            return 0.0
        return float('inf')
    if log_pstar == float('-inf'):
        return 0.0
    if log_pstar <= _LOG_Q_MIN + _LOG3:
        return 0.0
    return math.exp(log_pstar)


# ---------------------------------------------------------------------------
# Context-aware CER (unified equation) — see docs/cer_in_practice/
# ---------------------------------------------------------------------------


def compute_joint_critical_q(N: int, M: int, n_sites: int,
                             K: int = 1, alpha: float = 1e-5) -> Optional[float]:
    """Solve C(n_sites, K) * P(Binom(N, q*) >= M) = alpha for joint q*.

    The unified CER equation (Walker 2026c §2). Unlike compute_critical_error_rate,
    this returns the joint per-read error probability directly, with no
    uniform-model q = p/3 conversion. The result is in (0, 1).

    Args:
        N: Read count denominator (typically the identity-group sum).
        M: Variant read count.
        n_sites: Number of independent test sites for Bonferroni correction.
            Typically (non-HP positions) + (L>=2 HP runs); approximately L for
            short HP-run densities.
        K: Number of variant positions distinguishing candidate from reference.
        alpha: Significance level.

    Returns:
        joint q* in (0, 1), or None if no solution exists in that range
        (M too large relative to N for any q < 1 to satisfy).
        Returns 0.0 if the test passes even at q ≈ 0, or if the solver's
        log_q* underflows when exponentiated; callers needing precision in
        that regime should use ``compute_log_joint_critical_q`` instead.
    """
    log_qstar = _solve_log_qstar(N, M, n_sites, K, alpha)
    if log_qstar is None:
        return None
    if log_qstar == float('-inf'):
        return 0.0
    if log_qstar <= _LOG_Q_MIN:
        return 0.0
    return math.exp(log_qstar)


def compute_log_joint_critical_q(N: int, M: int, n_sites: int,
                                 K: int = 1,
                                 alpha: float = 1e-5) -> Optional[float]:
    """log(joint q*) — direct log-space companion to compute_joint_critical_q.

    Use this when downstream arithmetic (per-position root, factor ratio)
    would lose precision under the exp/log round-trip — i.e., whenever K is
    large enough that q* falls below ~1e-15.

    Returns:
        log(q*) in (-700, 0), or
        ``None`` if M is too large for any q < 1 to satisfy, or
        ``float('-inf')`` if the test passes even at q ≈ 0.
    """
    return _solve_log_qstar(N, M, n_sites, K, alpha)


def compute_per_position_qstar(N: int, M: int, n_sites: int,
                               K: int = 1, alpha: float = 1e-5) -> Optional[float]:
    """K-th root of joint q*. Per-position critical rate under uniform error.

    Useful for reporting alongside cer_factor: a per-position rate the reader
    can compare against platform expectations.
    """
    log_joint = _solve_log_qstar(N, M, n_sites, K, alpha)
    if log_joint is None:
        return None
    if log_joint == float('-inf'):
        return 0.0
    log_per_position = log_joint / K
    if log_per_position <= _LOG_Q_MIN:
        return 0.0
    return math.exp(log_per_position)


def compute_cer_factor(N: int, M: int, n_sites: int,
                       q_ctx_per_position: Sequence[float],
                       alpha: float = 1e-5) -> Optional[float]:
    """Compute the CER factor for a context-aware pairwise comparison.

    The factor is the per-position multiplicative inflation that the empirical
    error rates would need to undergo to make the variant plausible as artifact:

        factor = (joint_q* / product(q_ctx_i))^(1/K)

    For K=1 this reduces to factor = q* / q_ctx, matching the paper's headline
    definition. For K>1 the K-th root keeps the factor in per-position units,
    so a factor of 4 always means "each position would need to error at 4x its
    empirical rate" regardless of K.

    Args:
        N: Read count denominator.
        M: Variant read count.
        n_sites: Number of independent test sites for Bonferroni correction.
        q_ctx_per_position: Sequence of per-position empirical error rates,
            one per differing position (length K).
        alpha: Significance level.

    Returns:
        The per-position CER factor (dimensionless), or None if any q_ctx is
        invalid or the equation has no solution. Returns float('inf') when
        the joint q* solver indicates the variant is so strongly supported
        that no error inflation can explain it (joint q* > 1 is impossible).
    """
    if not q_ctx_per_position:
        return None
    K = len(q_ctx_per_position)

    log_actual_joint = 0.0
    for q in q_ctx_per_position:
        if q is None or q <= 0.0 or q >= 1.0:
            return None
        log_actual_joint += math.log(q)

    log_joint_qstar = _solve_log_qstar(N, M, n_sites, K, alpha)
    if log_joint_qstar is None:
        return float('inf')   # M too large for any q < 1 — strongest possible signal
    if log_joint_qstar == float('-inf'):
        return 0.0            # essentially unreachable; preserved for defensiveness

    log_factor = (log_joint_qstar - log_actual_joint) / K
    try:
        return math.exp(log_factor)
    except OverflowError:
        return float('inf')
