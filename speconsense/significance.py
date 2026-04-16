"""Statistical significance testing for variant phasing.

Implements the critical error rate (p*) framework from the variant significance
paper. Uses the binomial survival function with Bonferroni correction to
determine whether an observed variant count M out of N total reads at L positions
is statistically distinguishable from sequencing error.

Error model: uniform (q = p/3), where p is total per-position error rate.

Interpretability constraint: Under the uniform model, the implied reference
population (reads carrying the true sequence) at error rate p* is N(1-p*).
At p* = 0.75 (the signal destruction threshold), the reference population
equals the variant count — each of the 4 bases is equally likely. Above this,
the "true sequence" has fewer supporting reads than the variant, making the
artifact hypothesis incoherent. CER is therefore only reported when p* < 0.75.
"""

import math
from typing import List, Optional, Sequence, Tuple

from scipy.stats import binom
from scipy.optimize import brentq

# Signal destruction threshold: at p* = 0.75, each of the 4 bases has equal
# probability (25%) under the uniform error model. Above this, the implied
# reference population is smaller than the variant count.
SIGNAL_DESTRUCTION_THRESHOLD = 0.75


def compute_critical_error_rate(N: int, M: int, L: int, alpha: float = 1e-5,
                                K: int = 1) -> float:
    """Compute p* under uniform error model (q = p/3).

    For K=1: L * P(Binom(N, p*/3) >= M) = alpha
    For K>1: C(L,K) * P(Binom(N, (p*/3)^K) >= M) = alpha

    The K>1 model requires the SAME M reads to error at ALL K positions.
    Under independent errors across positions, the per-read probability of
    matching the variant at all K positions is q^K.

    Returns the raw p* without capping. Values at or above
    SIGNAL_DESTRUCTION_THRESHOLD (0.75) indicate the variant has more support
    than the implied reference population under the error model. Callers
    should check p* < SIGNAL_DESTRUCTION_THRESHOLD before reporting.

    Args:
        N: Total specimen reads (denominator)
        M: Variant read count
        L: Amplicon length (number of positions for Bonferroni correction)
        alpha: Significance level
        K: Number of variant positions (default 1)

    Returns:
        p* as total error rate (uncapped), or 0.0 if M <= 0
    """
    if M <= 0 or N <= 0 or L <= 0 or K <= 0:
        return 0.0
    if M > N:
        return 0.0
    if K > L:
        return 0.0

    # We solve: log(C(L,K)) + logsf(M-1, N, (p/3)^K) - log(alpha) = 0
    # logsf(M-1, N, q) = log(P(X >= M)) where X ~ Binom(N, q)
    log_alpha = math.log(alpha)
    log_correction = math.log(math.comb(L, K))

    def objective(p):
        q = p / 3.0
        q_joint = q ** K
        # binom.logsf(M-1, N, q_joint) = log(P(X >= M))
        return log_correction + binom.logsf(M - 1, N, q_joint) - log_alpha

    # Search bounds: p in (epsilon, 3.0) since q = p/3 <= 1.0
    eps = 1e-15

    # Check if solution exists within bounds
    try:
        val_low = objective(eps)
        val_high = objective(3.0)
    except (ValueError, OverflowError):
        return 0.0

    # If objective is negative at upper bound, no solution (M too large relative to N)
    if val_high < 0:
        return float('inf')

    # If objective is positive at lower bound, M is too small
    if val_low > 0:
        return 0.0

    try:
        p_star = brentq(objective, eps, 3.0, xtol=1e-10)
        return p_star
    except (ValueError, RuntimeError):
        return 0.0


def is_cer_reportable(p_star: float) -> bool:
    """Check whether a computed p* value is meaningful for reporting.

    Under the uniform error model, at p* = 0.75 (the signal destruction
    threshold), each of the 4 bases is equally likely — the reference
    population carrying the true sequence equals the variant count. Above
    this, the artifact hypothesis is incoherent because the "true sequence"
    has fewer supporting reads than the variant.

    p* values below the threshold are physically meaningful and can be
    reported as CER annotations. p* values at or above the threshold
    indicate the variant is too well-supported for the error framework
    to apply.

    Args:
        p_star: Computed critical error rate

    Returns:
        True if p* is below the signal destruction threshold (meaningful),
        False otherwise.
    """
    return p_star < SIGNAL_DESTRUCTION_THRESHOLD


def compute_minimum_M(N: int, L: int, alpha: float = 1e-5,
                      assumed_error_rate: float = 0.015,
                      K: int = 1) -> int:
    """Compute minimum M for significance at assumed error rate.

    Binary search over M to find smallest M where p* >= assumed_error_rate.

    Args:
        N: Total specimen reads
        L: Amplicon length
        alpha: Significance level
        assumed_error_rate: Assumed per-position error rate
        K: Number of variant positions (default 1)

    Returns:
        Minimum M for significance
    """
    if N <= 0 or L <= 0 or assumed_error_rate <= 0:
        return 1

    # Binary search: find smallest M where p* >= assumed_error_rate
    lo, hi = 1, N

    # Check if even M=N is insufficient
    p_star_max = compute_critical_error_rate(N, N, L, alpha, K)
    if p_star_max < assumed_error_rate:
        return N + 1  # Impossible to achieve significance

    while lo < hi:
        mid = (lo + hi) // 2
        p_star = compute_critical_error_rate(N, mid, L, alpha, K)
        if p_star >= assumed_error_rate:
            hi = mid
        else:
            lo = mid + 1

    return lo


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
        Returns 0.0 if M is small enough that even q=0 satisfies (highly
        significant).
    """
    if M <= 0 or N <= 0 or n_sites <= 0 or K <= 0:
        return None
    if M > N:
        return None
    if K > n_sites:
        return None

    log_alpha = math.log(alpha)
    log_correction = math.log(math.comb(n_sites, K))

    def objective(q: float) -> float:
        return log_correction + binom.logsf(M - 1, N, q) - log_alpha

    eps = 1e-15
    upper = 1.0 - eps
    try:
        val_low = objective(eps)
        val_high = objective(upper)
    except (ValueError, OverflowError):
        return None

    if val_high < 0:
        return None    # even q=1 doesn't bring the prob above alpha (M too big)
    if val_low > 0:
        return 0.0     # q=0 already satisfies (M too small to need significance)

    try:
        return brentq(objective, eps, upper, xtol=1e-12)
    except (ValueError, RuntimeError):
        return None


def compute_per_position_qstar(N: int, M: int, n_sites: int,
                               K: int = 1, alpha: float = 1e-5) -> Optional[float]:
    """K-th root of joint q*. Per-position critical rate under uniform error.

    Useful for reporting alongside cer_factor: a per-position rate the reader
    can compare against platform expectations.
    """
    joint = compute_joint_critical_q(N, M, n_sites, K, alpha)
    if joint is None:
        return None
    if joint <= 0:
        return 0.0
    return joint ** (1.0 / K)


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

    actual_joint = 1.0
    for q in q_ctx_per_position:
        if q is None or q <= 0:
            return None
        actual_joint *= q
    if actual_joint <= 0 or actual_joint >= 1.0:
        return None

    joint_qstar = compute_joint_critical_q(N, M, n_sites, K, alpha)
    if joint_qstar is None:
        return float('inf')   # M too large for any q < 1 — strongest possible signal
    if joint_qstar <= 0:
        return 0.0            # M too small — variant fails CER even with no error

    return (joint_qstar / actual_joint) ** (1.0 / K)


def format_cer_details(pstar: Optional[float], details: Optional[dict]) -> Optional[str]:
    """Format the cer_details FASTA annotation.

    Renders a structured semicolon-separated string per the CER-in-practice
    paper §3.5. For K=1: ``p*=...;ctx=...;q=...``. For K>1, K is included and
    per-position values are joined with ``+``: ``p*=...;K=2;ctx=a+b;q=x+y``.

    Args:
        pstar: Per-position critical rate (output of compute_per_position_qstar).
            May be None when the candidate has no valid pairwise comparison.
        details: Dict with keys ``K``, ``tags`` (list[str]), ``q_ctx``
            (list[float]), and ``ref_idx`` (int). May be None.

    Returns:
        The formatted string, or None when either input is missing.
    """
    if pstar is None or not details:
        return None
    K = details.get('K') or len(details.get('tags', []))
    tags = details.get('tags', [])
    qs = details.get('q_ctx', [])
    if not tags or not qs:
        return None

    parts = [f"p*={pstar:.4f}"]
    if K and K > 1:
        parts.append(f"K={K}")
    parts.append("ctx=" + "+".join(tags))
    parts.append("q=" + "+".join(f"{q:.4f}" for q in qs))
    return ";".join(parts)


def is_variant_significant(M: int, N: int, L: int,
                           assumed_error_rate: float = 0.015,
                           alpha: float = 1e-5,
                           K: int = 1) -> Tuple[bool, float]:
    """Test whether a variant is statistically significant.

    A variant is significant when p* >= assumed_error_rate, meaning the
    observed pattern would require an error rate at least as high as the
    assumed rate to be explained as artifact.

    Always computes p* and uses it for the significance decision, even when
    the value is above the signal destruction threshold. Callers that need
    to display p* should check is_cer_reportable() to determine whether the
    value is physically meaningful.

    Args:
        M: Variant read count (smallest haplotype at split position)
        N: Total specimen reads
        L: Amplicon length (consensus length)
        assumed_error_rate: Assumed per-position error rate
        alpha: Significance level
        K: Number of variant positions (default 1)

    Returns:
        Tuple of (is_significant, p_star).
        When assumed_error_rate <= 0, returns (True, p_star) (feature disabled).
    """
    p_star = compute_critical_error_rate(N, M, L, alpha, K)

    if assumed_error_rate <= 0:
        return True, p_star

    return p_star >= assumed_error_rate, p_star
