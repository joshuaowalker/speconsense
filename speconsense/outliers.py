"""Outlier read detection based on per-read identity distributions.

Provides ``detect_rid_outliers`` for flagging reads whose identity to their
cluster consensus is anomalously low relative to their clustermates. Used at
multiple pipeline stages to remove reads that are dragging down cluster
quality but are not caught by absolute-threshold filters.

Two complementary rules run in parallel; a read is flagged if either
triggers:

1. MAD-based modified Z-score. Principled under a normal-like distribution,
   robust to the single-outlier case via use of median and median-absolute-
   deviation rather than mean/stdev. Weaker at small N where only a few
   abs-deviations enter the MAD median.

2. Gap rule on the minimum. Flags the worst read when the gap between it
   and the next-worst is disproportionately large compared to the spread
   among the rest. Works at all N, especially useful at N=3,4 where MAD's
   median-of-deviations is pulled toward zero by the median's own zero
   deviation.

Only low-rid outliers are flagged; high-rid values are never considered
outliers.
"""

from typing import List, Sequence


def _median(sorted_values: Sequence[float]) -> float:
    n = len(sorted_values)
    if n == 0:
        raise ValueError("median of empty sequence")
    mid = n // 2
    if n % 2 == 1:
        return sorted_values[mid]
    return (sorted_values[mid - 1] + sorted_values[mid]) / 2.0


def detect_rid_outliers(
    rids: Sequence[float],
    modified_z_threshold: float = 3.0,
    gap_factor: float = 2.5,
    min_mad: float = 0.002,
    min_drop_from_median: float = 0.02,
) -> List[int]:
    """Return indices of reads flagged as low-rid outliers.

    Args:
        rids: Per-read identity values (fractions in [0, 1] or percentages;
            both work since the method is scale-invariant, but
            ``min_mad`` and ``min_drop_from_median`` defaults assume the
            [0, 1] scale — pass adjusted values for percentages). Order
            is preserved in returned indices.
        modified_z_threshold: Threshold for the MAD-based modified Z-score.
            A read with ``modified_Z < -threshold`` is flagged. The 3.5
            value commonly cited in the literature is tuned for larger-N
            outlier detection; for speconsense clusters of 3-10 reads with
            a single borderline outlier, 3.0 gives better recall without
            obviously over-flagging clean clusters.
        gap_factor: Threshold for the gap rule. The worst read is flagged
            if ``(r_second - r_worst) > gap_factor * (r_best - r_second)``.
        min_mad: Floor for MAD to avoid divide-by-zero when most reads are
            identical. Expressed in the same units as ``rids`` (default
            0.002 = 0.2pp for fractional rids).
        min_drop_from_median: Safety threshold on absolute magnitude. A
            read is not flagged unless its rid is at least this far below
            the cluster median. Prevents flagging reads that are
            statistically "far" in a cluster of near-identical reads but
            not meaningfully so. Default 0.02 = 2pp, roughly the expected
            per-read ONT noise at Q20+.

    Returns:
        Sorted list of flagged indices (into the input ``rids`` sequence).
    """
    n = len(rids)
    if n < 3:
        # At N=1,2 neither rule is meaningful.
        return []

    sorted_rids = sorted(rids)
    median = _median(sorted_rids)

    flagged = set()

    # --- Rule 1: MAD-based modified Z-score ---
    abs_devs = sorted(abs(r - median) for r in rids)
    mad = _median(abs_devs)
    mad_eff = max(mad, min_mad)
    for i, r in enumerate(rids):
        modified_z = 0.6745 * (r - median) / mad_eff
        if modified_z < -modified_z_threshold:
            flagged.add(i)

    # --- Rule 2: Gap rule on the worst read ---
    sorted_indices = sorted(range(n), key=lambda i: rids[i])
    r_worst_idx = sorted_indices[0]
    r_worst = rids[r_worst_idx]
    r_second = rids[sorted_indices[1]]
    r_best = rids[sorted_indices[-1]]

    gap_bottom = r_second - r_worst
    spread_top = r_best - r_second
    if spread_top > 1e-9:
        if gap_bottom > gap_factor * spread_top:
            flagged.add(r_worst_idx)
    else:
        # Degenerate: all other reads are identical. Flag worst if any
        # appreciable gap exists at all (>0.5pp by default scale).
        if gap_bottom > min_mad * 5:
            flagged.add(r_worst_idx)

    # --- Safety: require meaningful absolute drop from median ---
    # Both rules are scale-invariant and fire on statistically-unusual shape
    # regardless of magnitude. In a cluster of near-identical reads, the
    # "worst" might be only a fraction of a percent below the others —
    # statistically distinct but biologically uninteresting. Require the
    # flagged read's rid to be at least ``min_drop_from_median`` below the
    # cluster median before accepting the flag.
    return sorted(
        i for i in flagged
        if (median - rids[i]) >= min_drop_from_median
    )
