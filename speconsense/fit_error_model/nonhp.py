"""Approach 2: all-cluster pooled non-HP error rates.

Iterates every ``cluster_debug/*-msa.fasta`` with ``RiC >= nonhp_min_ric``
and computes per-position errors against each cluster's own consensus,
excluding positions inside HP runs of length > 1. The pooled rates match
the operational distribution CER evaluates against in production, where
non-anchor candidate clusters of all sizes are validated against their
peers (HP paper §8.4).

This module is independent of selection: it uses every cluster's MSA,
not just qualifying primary anchors.
"""

from __future__ import annotations

import glob
import logging
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional, Tuple

from .extract import (
    NonHPCounts,
    build_position_mappings,
    compute_non_hp_error_counts,
    hp_positions_of_runs,
    identify_hp_runs,
    parse_msa_lightweight,
)

_RIC_RE = re.compile(r'-RiC(\d+)-')


@dataclass
class NonHPBin:
    """Pooled non-HP counts plus cluster/read totals for one RiC bin."""
    label: str
    n_clusters: int = 0
    n_reads: int = 0
    counts: NonHPCounts = None  # type: ignore[assignment]

    def __post_init__(self):
        if self.counts is None:
            self.counts = NonHPCounts()


def _process_msa(msa_path: str) -> Optional[Tuple[NonHPCounts, int]]:
    """Return ``(counts, n_reads)`` for one MSA, or None on failure."""
    try:
        consensus_aligned, consensus_ungapped, read_sequences = parse_msa_lightweight(msa_path)
    except Exception as e:  # noqa: BLE001
        logging.debug(f"approach-2 skip {msa_path}: {e}")
        return None
    if not read_sequences:
        return None
    _ = build_position_mappings(consensus_aligned)
    hp_runs = identify_hp_runs(consensus_ungapped)
    hp_pos = hp_positions_of_runs(hp_runs)
    counts = compute_non_hp_error_counts(read_sequences, consensus_aligned, hp_pos)
    return counts, len(read_sequences)


def _ric_from_path(path: str) -> Optional[int]:
    m = _RIC_RE.search(path)
    return int(m.group(1)) if m else None


def _bin_label(ric: int, edges: list) -> Optional[str]:
    for low, high in zip(edges[:-1], edges[1:]):
        if low <= ric < high:
            return f"{low}-{high if high < 10**8 else 'inf'}"
    return None


def aggregate_nonhp(
    input_dir: str,
    nonhp_min_ric: int = 5,
    threads: int = 1,
    progress: bool = True,
) -> Tuple[NonHPCounts, int, int, dict]:
    """Walk every cluster MSA and pool non-HP error counts.

    Returns ``(pooled, n_clusters, n_skipped, by_bin)`` where ``by_bin`` is
    a dict mapping bin label (str) to a ``NonHPBin``. Bin edges mirror the
    paper's stratification: 5-20, 20-50, 50-100, 100-200, 200-500,
    500-1000, 1000-inf. The bins are diagnostic only — the YAML output
    uses the pooled totals.
    """
    debug = os.path.join(input_dir, "cluster_debug")
    msa_paths = sorted(glob.glob(os.path.join(debug, "*-msa.fasta")))

    edges = [nonhp_min_ric, 20, 50, 100, 200, 500, 1000, 10**9]
    edges = sorted(set(e for e in edges if e >= nonhp_min_ric))

    by_bin: dict = {f"{lo}-{hi if hi < 10**8 else 'inf'}": NonHPBin(label=f"{lo}-{hi if hi < 10**8 else 'inf'}")
                     for lo, hi in zip(edges[:-1], edges[1:])}
    pooled = NonHPCounts()
    n_clusters = 0
    n_skipped = 0

    try:
        from tqdm import tqdm
    except ImportError:
        def tqdm(it, **_kw):
            return it

    # Pre-filter by RiC from path (cheap)
    eligible = [(p, r) for p in msa_paths for r in [_ric_from_path(p)] if r is not None and r >= nonhp_min_ric]

    def accumulate(counts: NonHPCounts, n_reads: int, ric: int):
        nonlocal n_clusters
        n_clusters += 1
        pooled.matches += counts.matches
        pooled.substitutions += counts.substitutions
        pooled.deletions += counts.deletions
        pooled.insertions += counts.insertions
        pooled.total_positions += counts.total_positions
        label = _bin_label(ric, edges)
        if label is not None:
            bin_ = by_bin[label]
            bin_.n_clusters += 1
            bin_.n_reads += n_reads
            bin_.counts.matches += counts.matches
            bin_.counts.substitutions += counts.substitutions
            bin_.counts.deletions += counts.deletions
            bin_.counts.insertions += counts.insertions
            bin_.counts.total_positions += counts.total_positions

    if threads <= 1:
        for path, ric in tqdm(eligible, desc="Approach 2 (non-HP)", disable=not progress):
            out = _process_msa(path)
            if out is None:
                n_skipped += 1
                continue
            counts, n_reads = out
            accumulate(counts, n_reads, ric)
    else:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            futures = {ex.submit(_process_msa, path): ric for path, ric in eligible}
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc="Approach 2 (non-HP)", disable=not progress):
                out = fut.result()
                if out is None:
                    n_skipped += 1
                    continue
                counts, n_reads = out
                accumulate(counts, n_reads, futures[fut])

    return pooled, n_clusters, n_skipped, by_bin


def derive_nonhp_rates(pooled: NonHPCounts) -> Tuple[float, float]:
    """Return ``(non_hp_sub, non_hp_indel)`` from pooled counts.

    Matches the shipped-model construction (HP paper §8.5 Table 16, see
    ``dorado-v5.0.yaml`` comments): denominator is
    ``matches + substitutions + deletions`` (consensus-column observations,
    insertions are tallied separately at insertion columns); ``non-hp-sub``
    is the substitution rate; ``non-hp-indel`` is ``(deletions +
    insertions) / denominator``. Including the insertion count makes the
    indel rate inflate slightly above the column-level deletion rate alone,
    which is intentional in the shipped scheme.
    """
    denom = pooled.matches + pooled.substitutions + pooled.deletions
    if denom <= 0:
        return 0.0, 0.0
    sub_rate = pooled.substitutions / denom
    indel_rate = (pooled.deletions + pooled.insertions) / denom
    return sub_rate, indel_rate
