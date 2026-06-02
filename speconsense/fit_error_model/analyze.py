"""Aggregation and biological-variant filtering for approach-1 HP observations.

Two filters work together to exclude real biological hp-only variants from
the error-rate estimate (HP paper §5):

1. **Outlier filter**: a one-sided binomial test asks whether a given
   position's deviation count exceeds what its (base, hp_length) pool's
   expected error rate would predict. Bonferroni correction over the total
   number of tested positions controls the false-positive rate.
2. **Bimodal filter**: positions whose second-most-common called length
   carries a fraction substantially above the context's expected error
   rate are flagged as likely heterozygous / multi-template signal.

After flagging, the same aggregation is rerun on the surviving positions
to produce the final per-context error rates that feed the YAML output.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

import numpy as np
from scipy import stats as scipy_stats

from .extract import HPObservation


@dataclass
class ContextResult:
    """Aggregated stats for one (base, hp_length) context."""
    n_runs: int = 0
    n_obs: int = 0
    all_called: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int16))
    frac_correct: float = 0.0
    frac_minus1: float = 0.0
    frac_plus1: float = 0.0
    frac_minus2plus: float = 0.0
    frac_plus2plus: float = 0.0
    distribution: Dict[int, int] = field(default_factory=dict)
    mean_called: float = 0.0
    std_called: float = 0.0


def aggregate_by_context(
    observations: List[HPObservation],
) -> Dict[Tuple[str, int], ContextResult]:
    """Group observations by ``(base, mode_length)`` and compute error stats.

    Returns a dict keyed by ``(base, hp_length)``.
    """
    groups: Dict[Tuple[str, int], List[np.ndarray]] = defaultdict(list)
    for obs in observations:
        groups[(obs.base, obs.mode_length)].append(obs.called_lengths)

    results: Dict[Tuple[str, int], ContextResult] = {}
    for (base, mode_len), arrays in groups.items():
        all_called = np.concatenate(arrays) if arrays else np.array([], dtype=np.int16)
        deltas = all_called - mode_len
        n = len(deltas)
        if n == 0:
            continue
        unique, counts = np.unique(deltas, return_counts=True)
        results[(base, mode_len)] = ContextResult(
            n_runs=len(arrays),
            n_obs=n,
            all_called=all_called,
            frac_correct=float(np.mean(deltas == 0)),
            frac_minus1=float(np.mean(deltas == -1)),
            frac_plus1=float(np.mean(deltas == 1)),
            frac_minus2plus=float(np.mean(deltas <= -2)),
            frac_plus2plus=float(np.mean(deltas >= 2)),
            distribution={int(u): int(c) for u, c in zip(unique, counts)},
            mean_called=float(np.mean(all_called)),
            std_called=float(np.std(all_called)),
        )
    return results


def detect_outliers(
    observations: List[HPObservation],
    context_results: Dict[Tuple[str, int], ContextResult],
    alpha: float = 0.01,
) -> List[Tuple[str, int]]:
    """Flag (specimen, consensus_pos) tuples whose error count is improbable
    under the per-context pooled rate.

    One-sided binomial survival probability, Bonferroni-corrected over the
    number of positions tested.
    """
    n_tests = len(observations)
    if n_tests == 0:
        return []
    threshold = alpha / n_tests
    flagged: List[Tuple[str, int]] = []
    for obs in observations:
        key = (obs.base, obs.mode_length)
        ctx = context_results.get(key)
        if ctx is None:
            continue
        p_err = 1.0 - ctx.frac_correct
        if p_err <= 0:
            continue
        n = obs.n_reads
        k = int(np.sum(obs.called_lengths != obs.mode_length))
        # P(X >= k) under Binomial(n, p_err)
        pval = scipy_stats.binom.sf(k - 1, n, p_err)
        if pval < threshold:
            flagged.append((obs.specimen, obs.consensus_pos))
    return flagged


def detect_bimodal(
    observations: List[HPObservation],
    context_results: Dict[Tuple[str, int], ContextResult],
    min_fraction: float = 0.10,
    rate_ratio: float = 2.0,
    min_reads: int = 20,
) -> List[Tuple[str, int]]:
    """Flag positions whose second-most-common length frac exceeds both
    ``min_fraction`` AND ``rate_ratio * expected_error_rate``.

    Mirrors HP paper §5.2: bimodal signal at frequencies far above the
    sequencing noise floor is the signature of a real second template.
    """
    flagged: List[Tuple[str, int]] = []
    for obs in observations:
        called = obs.called_lengths
        if len(called) < min_reads:
            continue
        unique, counts = np.unique(called, return_counts=True)
        if len(counts) < 2:
            continue
        order = np.argsort(-counts)
        second_frac = counts[order[1]] / len(called)
        key = (obs.base, obs.mode_length)
        ctx = context_results.get(key)
        if ctx is None:
            continue
        expected = 1.0 - ctx.frac_correct
        if expected <= 0:
            continue
        if second_frac > rate_ratio * expected and second_frac > min_fraction:
            flagged.append((obs.specimen, obs.consensus_pos))
    return flagged


def filtered_observations(
    observations: List[HPObservation],
    flagged: Set[Tuple[str, int]],
) -> List[HPObservation]:
    """Return observations whose (specimen, consensus_pos) is not in ``flagged``."""
    if not flagged:
        return list(observations)
    return [o for o in observations if (o.specimen, o.consensus_pos) not in flagged]


def run_filter_pass(
    observations: List[HPObservation],
    alpha: float = 0.01,
    bimodal_min_frac: float = 0.10,
    bimodal_rate_ratio: float = 2.0,
) -> Tuple[Dict[Tuple[str, int], ContextResult], Dict[Tuple[str, int], ContextResult], int]:
    """Aggregate → flag (outlier + bimodal) → re-aggregate.

    Returns ``(raw_results, filtered_results, n_excluded_runs)``.
    """
    raw = aggregate_by_context(observations)
    outliers = detect_outliers(observations, raw, alpha=alpha)
    bimodal = detect_bimodal(
        observations, raw,
        min_fraction=bimodal_min_frac,
        rate_ratio=bimodal_rate_ratio,
    )
    flagged: Set[Tuple[str, int]] = set(outliers) | set(bimodal)
    kept = filtered_observations(observations, flagged)
    filtered = aggregate_by_context(kept)
    return raw, filtered, len(observations) - len(kept)
