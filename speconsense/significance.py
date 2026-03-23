"""Statistical significance testing for variant phasing.

Implements the critical error rate (p*) framework from the variant significance
paper. Uses the binomial survival function with Bonferroni correction to
determine whether an observed variant count M out of N total reads at L positions
is statistically distinguishable from sequencing error.

Error model: uniform (q = p/3), where p is total per-position error rate.
"""

import math
from typing import Tuple

from scipy.stats import binom
from scipy.optimize import brentq


def compute_critical_error_rate(N: int, M: int, L: int, alpha: float = 1e-5) -> float:
    """Compute p* under uniform error model (q = p/3).

    Solves: L * P(Binom(N, p*/3) >= M) = alpha for p*.

    Args:
        N: Total specimen reads (denominator)
        M: Variant read count
        L: Amplicon length (number of positions for Bonferroni correction)
        alpha: Significance level

    Returns:
        p* as total error rate, or 0.0 if M <= 0
    """
    if M <= 0 or N <= 0 or L <= 0:
        return 0.0
    if M > N:
        return 0.0

    # We solve: log(L) + logsf(M-1, N, p/3) - log(alpha) = 0
    # logsf(M-1, N, q) = log(P(X >= M)) where X ~ Binom(N, q)
    log_alpha = math.log(alpha)
    log_L = math.log(L)

    def objective(p):
        q = p / 3.0
        # binom.logsf(M-1, N, q) = log(P(X >= M))
        return log_L + binom.logsf(M - 1, N, q) - log_alpha

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
        return 1.0  # p* exceeds maximum — obviously not error

    # If objective is positive at lower bound, M is too small
    if val_low > 0:
        return 0.0

    try:
        p_star = brentq(objective, eps, 3.0, xtol=1e-10)
        return min(p_star, 1.0)
    except (ValueError, RuntimeError):
        return 0.0


def compute_minimum_M(N: int, L: int, alpha: float = 1e-5,
                      assumed_error_rate: float = 0.02) -> int:
    """Compute minimum M for significance at assumed error rate.

    Binary search over M to find smallest M where p* >= assumed_error_rate.

    Args:
        N: Total specimen reads
        L: Amplicon length
        alpha: Significance level
        assumed_error_rate: Assumed per-position error rate

    Returns:
        Minimum M for significance
    """
    if N <= 0 or L <= 0 or assumed_error_rate <= 0:
        return 1

    # Binary search: find smallest M where p* >= assumed_error_rate
    lo, hi = 1, N

    # Check if even M=N is insufficient
    p_star_max = compute_critical_error_rate(N, N, L, alpha)
    if p_star_max < assumed_error_rate:
        return N + 1  # Impossible to achieve significance

    while lo < hi:
        mid = (lo + hi) // 2
        p_star = compute_critical_error_rate(N, mid, L, alpha)
        if p_star >= assumed_error_rate:
            hi = mid
        else:
            lo = mid + 1

    return lo


def is_variant_significant(M: int, N: int, L: int,
                           assumed_error_rate: float = 0.02,
                           alpha: float = 1e-5) -> Tuple[bool, float]:
    """Test whether a variant is statistically significant.

    A variant is significant when p* >= assumed_error_rate, meaning the
    observed pattern would require an error rate at least as high as the
    assumed rate to be explained as artifact.

    Args:
        M: Variant read count (smallest haplotype at split position)
        N: Total specimen reads
        L: Amplicon length (consensus length)
        assumed_error_rate: Assumed per-position error rate
        alpha: Significance level

    Returns:
        Tuple of (is_significant, p_star).
        When assumed_error_rate <= 0, returns (True, p_star) (feature disabled).
    """
    p_star = compute_critical_error_rate(N, M, L, alpha)

    if assumed_error_rate <= 0:
        return True, p_star

    return p_star >= assumed_error_rate, p_star
