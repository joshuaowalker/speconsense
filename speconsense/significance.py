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
from typing import Optional, Tuple

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
