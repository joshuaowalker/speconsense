#!/usr/bin/env python3
"""Compute critical error rates and generate figures for variant significance paper.

For a variant defined by M of N reads sharing an alternative allele at K positions
in an amplicon of length L, the critical error rate p* is the per-position nanopore
error rate at which the observed pattern would become plausible under the null
hypothesis of independent sequencing errors.

Two error models are computed:
  - Conservative (q = p): all errors at a position produce the same alt base
  - Uniform (q = p/3): errors spread equally across three alternative bases
"""

import math
import numpy as np
from scipy.stats import binom
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path

matplotlib.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

FIGURES_DIR = Path(__file__).parent / 'figures'


def critical_error_rate(N, M, L=700, alpha=0.05, K=1, uniform=False):
    """Solve for the critical error rate p*.

    Finds p* such that the Bonferroni-corrected probability of observing
    M or more reads with the same alternative allele at K positions equals alpha.

    For K=1: L * P(Binom(N, q) >= M) = alpha
    For K>1: C(L,K) * [P(Binom(N, q) >= M)]^K = alpha

    where q = p (conservative, uniform=False) or q = p/3 (uniform=True).

    Returns p* as a per-position error rate.
    """
    if K == 1:
        correction = L
    else:
        correction = math.comb(L, K)

    # Per-position survival function target
    if K > 1:
        per_position_target = (alpha / correction) ** (1.0 / K)
    else:
        per_position_target = alpha / correction

    log_target = math.log(per_position_target)

    divisor = 3.0 if uniform else 1.0

    def equation(p):
        q = p / divisor
        if q >= 1.0:
            q = 1.0 - 1e-15
        return binom.logsf(M - 1, N, q) - log_target

    p_upper = 1.0 if not uniform else 3.0
    p_star = brentq(equation, 1e-15, p_upper - 1e-15, xtol=1e-12)
    return p_star


def effective_M(N, min_count=5, min_freq=0.10):
    """Minimum M satisfying both count and frequency thresholds (default params)."""
    return max(min_count, math.ceil(min_freq * N))


def minimum_M(N, p_assumed, L=700, alpha=0.05, uniform=False):
    """Minimum variant read count M for significance at assumed error rate.

    Returns the smallest M such that observing M of N reads with the same
    alternative allele is significant at level alpha (Bonferroni-corrected
    over L positions) when the true per-position error rate is p_assumed.

    Returns None if no M <= N is sufficient.
    """
    q = p_assumed / 3.0 if uniform else p_assumed
    target = alpha / L
    # Find M where P(Binom(N, q) >= M) <= target
    M_raw = binom.isf(target, N, q)
    M_cand = int(math.ceil(M_raw))
    # Verify (isf can be imprecise at discrete boundaries)
    while M_cand > 0 and binom.sf(M_cand - 1, N, q) > target:
        M_cand += 1
    # Check feasibility
    if M_cand > N:
        return None
    return M_cand


def generate_figure1():
    """Critical error rate vs N for K=1 (main figure).

    Shows p* vs M curves at multiple significance levels.
    """
    M_values = np.arange(5, 201)
    N = 1000

    alphas = [0.05, 1e-3, 1e-5]
    alpha_labels = [
        r'$\alpha = 0.05$',
        r'$\alpha = 10^{-3}$',
        r'$\alpha = 10^{-5}$',
    ]
    colors_cons = ['#1f77b4', '#ff7f0e', '#2ca02c']

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Conservative (q=p) curves, solid
    for a, label, color in zip(alphas, alpha_labels, colors_cons):
        curve = []
        for M in M_values:
            curve.append(critical_error_rate(N, M, alpha=a))
        ax.plot(M_values, curve, '-', color=color, linewidth=1.8,
                label=f'{label} (q = p)')

    # Uniform (q=p/3) curves, dashed, for α=0.05 and α=10⁻⁵
    for a, label, color in zip(
        [0.05, 1e-5],
        [r'$\alpha = 0.05$', r'$\alpha = 10^{-5}$'],
        ['#1f77b4', '#2ca02c'],
    ):
        curve = []
        for M in M_values:
            curve.append(critical_error_rate(N, M, alpha=a, uniform=True))
        ax.plot(M_values, curve, '--', color=color, linewidth=1.2,
                label=f'{label} (q = p/3)', alpha=0.7)

    ax.axhline(y=0.05, color='red', linestyle=':', linewidth=1, alpha=0.7)
    ax.text(195, 0.053, '5%', fontsize=8, color='red', alpha=0.7, ha='right')
    ax.axhline(y=0.01, color='darkred', linestyle=':', linewidth=1, alpha=0.7)
    ax.text(195, 0.0107, '1%', fontsize=8, color='darkred', alpha=0.7, ha='right')

    ax.set_yscale('log')
    ax.set_xlabel('M (variant read count)')
    ax.set_ylabel(r'Critical error rate $p^*$')
    ax.set_title('Critical error rate vs variant support (N = 1000, K = 1)')
    ax.set_xlim(5, 200)
    ax.set_ylim(0.001, 1.0)
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    fig.savefig(FIGURES_DIR / 'fig1_critical_rate_vs_M.pdf')
    plt.close(fig)
    print("Figure 1 saved.")


def generate_figure2():
    """Critical error rate vs M for multiple N values (K=1)."""
    M_values = np.arange(5, 201)
    N_values_fixed = [50, 100, 200, 500, 1000]
    colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728', '#9467bd']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

    # Left panel: conservative (q=p)
    for N_fixed, color in zip(N_values_fixed, colors):
        curve = []
        for M in M_values:
            if M > N_fixed:
                curve.append(np.nan)
            else:
                curve.append(critical_error_rate(N_fixed, M))
        ax1.plot(M_values, curve, '-', color=color, linewidth=1.8,
                 label=f'N = {N_fixed}')

    ax1.axhline(y=0.05, color='red', linestyle=':', linewidth=1, alpha=0.7)
    ax1.axhline(y=0.01, color='darkred', linestyle=':', linewidth=1, alpha=0.7)
    ax1.set_yscale('log')
    ax1.set_xlabel('M (variant read count)')
    ax1.set_ylabel(r'Critical error rate $p^*$')
    ax1.set_title(r'Conservative ($q = p$)')
    ax1.set_xlim(5, 200)
    ax1.set_ylim(0.001, 1.0)
    ax1.legend(loc='lower right', fontsize=8)
    ax1.grid(True, alpha=0.3, which='both')

    # Right panel: uniform (q=p/3)
    for N_fixed, color in zip(N_values_fixed, colors):
        curve = []
        for M in M_values:
            if M > N_fixed:
                curve.append(np.nan)
            else:
                curve.append(critical_error_rate(N_fixed, M, uniform=True))
        ax2.plot(M_values, curve, '-', color=color, linewidth=1.8,
                 label=f'N = {N_fixed}')

    ax2.axhline(y=0.05, color='red', linestyle=':', linewidth=1, alpha=0.7)
    ax2.axhline(y=0.01, color='darkred', linestyle=':', linewidth=1, alpha=0.7)
    ax2.set_yscale('log')
    ax2.set_xlabel('M (variant read count)')
    ax2.set_title(r'Uniform ($q = p/3$)')
    ax2.set_xlim(5, 200)
    ax2.legend(loc='lower right', fontsize=8)
    ax2.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / 'fig2_critical_rate_by_N.pdf')
    plt.close(fig)
    print("Figure 2 saved.")


def print_key_values():
    """Print key numerical results for use in the paper text."""
    print("\n=== Key numerical values (conservative q=p) ===\n")

    cases = [
        (1000, 100, 1, "N=1000, M=100, K=1"),
        (100, 10, 1, "N=100, M=10, K=1"),
        (50, 5, 1, "N=50, M=5, K=1"),
        (10000, 1000, 1, "N=10000, M=1000, K=1"),
        (1000, 100, 2, "N=1000, M=100, K=2"),
        (1000, 5, 1, "N=1000, M=5, K=1"),
    ]

    for N, M, K, desc in cases:
        p_cons = critical_error_rate(N, M, K=K, uniform=False)
        p_unif = critical_error_rate(N, M, K=K, uniform=True)
        print(f"  {desc}: p*(q=p) = {p_cons*100:.2f}%, p*(q=p/3) = {p_unif*100:.2f}%")

    print("\n=== Effect of significance level (N=1000, M=100, K=1) ===\n")
    for alpha in [0.05, 1e-3, 1e-5]:
        p_cons = critical_error_rate(1000, 100, alpha=alpha)
        p_unif = critical_error_rate(1000, 100, alpha=alpha, uniform=True)
        print(f"  alpha={alpha:.0e}: p*(q=p) = {p_cons*100:.2f}%, p*(q=p/3) = {p_unif*100:.2f}%")

    print("\n=== Boundary cases at alpha=1e-5 (run-corrected) ===\n")
    boundary_cases = [
        (50, 5, 1, "N=50, M=5"),
        (100, 10, 1, "N=100, M=10"),
        (1000, 100, 1, "N=1000, M=100"),
    ]
    for N, M, K, desc in boundary_cases:
        p_cons = critical_error_rate(N, M, K=K, alpha=1e-5)
        p_unif = critical_error_rate(N, M, K=K, alpha=1e-5, uniform=True)
        print(f"  {desc}: p*(q=p) = {p_cons*100:.2f}%, p*(q=p/3) = {p_unif*100:.2f}%")

    print()


def print_minimum_M_table():
    """Print minimum M values for Table 2 (assumed error rate framing)."""
    N_values = [100, 200, 500, 1000]
    p_values = [0.01, 0.03, 0.05]
    alphas = [0.05, 1e-3, 1e-5]

    for model_name, uniform in [('Conservative (q=p)', False), ('Uniform (q=p/3)', True)]:
        print(f"\n=== Minimum M — {model_name} ===\n")
        # Header
        header = f"{'N':>6}"
        for p in p_values:
            for a in alphas:
                header += f"  {p*100:.0f}%/a={a:.0e}"
            header += "  |"
        print(header)
        print("-" * len(header))

        for N in N_values:
            row = f"{N:>6}"
            for p in p_values:
                for a in alphas:
                    M = minimum_M(N, p, alpha=a, uniform=uniform)
                    row += f"  {M if M is not None else '>N':>13}"
                row += "  |"
            print(row)
    print()


if __name__ == '__main__':
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    generate_figure1()
    generate_figure2()
    print_key_values()
    print_minimum_M_table()
