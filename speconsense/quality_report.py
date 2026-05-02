"""Quality report generation for speconsense-summarize.

The report is action-oriented: it surfaces variants that warrant inspection
or rescue rather than presenting raw statistics. Five sections:

1. **Executive summary** — specimen and cluster counts, filter routing, RiC
   and yield distributions. Calls out specimens that had input reads but
   produced no clusters.
2. **Passed variants worth inspecting** — variants that survived the CER
   and err_factor filters but stand out vs the run's own distribution
   (mean ± 2σ). Two sub-tables: high err_factor and high ambiguity.
3. **Possible rescues** — low-yield specimens (bottom-quartile yield) that
   have filtered variants close to the threshold (.ns with cer_factor > 0.9
   or .lq with err_factor < 1.7).
4. **Run-wide parameter signals** — aggregate patterns suggesting the
   user might want to retune defaults.
5. **Pipeline activity** — neutral run-wide stats (filter routing,
   merge counts, factor distributions).
"""

import glob
import json
import logging
import math
import os
import re
import statistics
from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable, Dict, List, Optional, Tuple

from speconsense.types import ConsensusInfo, OverlapMergeInfo


# ---------------------------------------------------------------------------
# q_ctx calibration tuning constants
# ---------------------------------------------------------------------------
# Pooled obs/exp ratio is a directional check; thresholds are weakly held.
# Asymmetry: underestimation is the more common failure mode, since residual
# heterogeneity inflates obs_sum and newer basecallers tend to be cleaner.
MIN_CALIBRATION_COLS = 5000
CALIBRATION_LOW = 0.70
CALIBRATION_HIGH = 1.40


# ---------------------------------------------------------------------------
# Per-cluster wrapper
# ---------------------------------------------------------------------------


@dataclass
class _Cluster:
    """View of a single source-level cluster for reporting purposes."""
    info: ConsensusInfo
    state: str           # "passed", "ns", or "lq"
    specimen: str
    short: str           # specimen-relative tail, e.g. "1.v3"

    @property
    def err_factor(self) -> Optional[float]:
        return self.info.err_factor

    @property
    def err_factor_obs_sum(self) -> Optional[float]:
        return self.info.err_factor_obs_sum

    @property
    def err_factor_exp_sum(self) -> Optional[float]:
        return self.info.err_factor_exp_sum

    @property
    def err_factor_cols(self) -> Optional[int]:
        return self.info.err_factor_cols

    @property
    def cer_factor(self) -> Optional[float]:
        return self.info.cer_factor

    @property
    def rid(self) -> Optional[float]:
        return self.info.rid

    @property
    def ric(self) -> int:
        return self.info.ric or 0

    @property
    def ambig_count(self) -> int:
        snp = self.info.snp_count or 0
        seq_ambig = sum(1 for c in self.info.sequence if c.upper() not in "ACGT-")
        return snp + seq_ambig


@dataclass
class _Specimen:
    name: str
    total_input_reads: Optional[int] = None
    clusters: List[_Cluster] = field(default_factory=list)

    @property
    def reads_in_passed(self) -> int:
        # Use size (actual cluster size), not ric (capped at --max-sample-size).
        # Yield should reflect how much of the raw input survived clustering,
        # not how many reads were sampled for the consensus.
        return sum(c.info.size or 0 for c in self.clusters if c.state == "passed")

    @property
    def yield_pct(self) -> Optional[float]:
        if not self.total_input_reads:
            return None
        return 100.0 * self.reads_in_passed / self.total_input_reads


def _specimen_from_name(cluster_name: str) -> str:
    return re.sub(r"-\d+\.v\d+(\.raw\d+)?$", "", cluster_name)


def _short_id(cluster_name: str) -> str:
    m = re.search(r"-(\d+\.v\d+(?:\.raw\d+)?)$", cluster_name)
    return m.group(1) if m else cluster_name


def _build_specimens(
    consensus_list: List[ConsensusInfo],
    ns_list: List[ConsensusInfo],
    lq_list: List[ConsensusInfo],
    source_folder: str,
) -> Tuple[Dict[str, _Specimen], List[Tuple[str, int]]]:
    """Group source-level clusters by specimen and pull total_input_reads.

    Returns (specimens, zero_cluster_specimens) where zero_cluster_specimens
    is a list of (specimen_name, input_reads) for specimens that had input
    reads but produced no source clusters.
    """
    specimens: Dict[str, _Specimen] = {}

    def _add(info: ConsensusInfo, state: str) -> None:
        spec_name = _specimen_from_name(info.sample_name)
        sp = specimens.setdefault(spec_name, _Specimen(spec_name))
        sp.clusters.append(_Cluster(
            info=info, state=state, specimen=spec_name,
            short=_short_id(info.sample_name),
        ))

    for info in consensus_list:
        _add(info, "passed")
    for info in ns_list:
        _add(info, "ns")
    for info in lq_list:
        _add(info, "lq")

    # Pull total_input_reads from per-specimen metadata JSONs. Also catch
    # specimens that had input but produced zero clusters.
    zero_cluster_specimens: List[Tuple[str, int]] = []
    debug_dir = os.path.join(source_folder, "cluster_debug")
    if os.path.isdir(debug_dir):
        for path in sorted(glob.glob(os.path.join(debug_dir, "*-metadata.json"))):
            base = os.path.basename(path).replace("-metadata.json", "")
            try:
                with open(path) as f:
                    data = json.load(f)
            except (OSError, json.JSONDecodeError):
                continue
            total = data.get("total_input_reads")
            if base in specimens:
                specimens[base].total_input_reads = total
            elif total and total > 0:
                zero_cluster_specimens.append((base, int(total)))

    return specimens, zero_cluster_specimens


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------


def _mean_std(values: List[float]) -> Tuple[float, float]:
    if len(values) < 2:
        return (values[0] if values else 0.0, 0.0)
    return statistics.mean(values), statistics.pstdev(values)


def _percentile(values: List[float], p: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    k = (len(s) - 1) * p
    f, c = math.floor(k), math.ceil(k)
    if f == c:
        return s[int(k)]
    return s[f] * (c - k) + s[c] * (k - f)


def _select_outliers_high(
    clusters: List[_Cluster],
    key: Callable[[_Cluster], Optional[float]],
    *,
    k: float = 2.0,
    floor: Optional[float] = None,
    cap: int = 20,
) -> Optional[Tuple[float, float, float, List[_Cluster], int]]:
    vals = [v for v in (key(c) for c in clusters) if v is not None]
    if not vals:
        return None
    mean, std = _mean_std(vals)
    threshold = mean + k * std
    if floor is not None and threshold < floor:
        threshold = floor
    selected = [c for c in clusters if (key(c) is not None and key(c) > threshold)]
    selected.sort(key=lambda c: key(c) or 0.0, reverse=True)
    return mean, std, threshold, selected[:cap], len(selected)


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------


def _section(title: str) -> str:
    return f"\n{'=' * 80}\n{title}\n{'=' * 80}\n"


def _format_pct(x: Optional[float]) -> str:
    """Format a 0..1 fraction as a right-aligned percentage."""
    return f"{x:>6.1%}" if x is not None else "     —"


def _format_factor(x: Optional[float]) -> str:
    return f"{x:>5.2f}" if x is not None else "    —"


def _render_header(source_folder: str, summary_folder: str,
                   min_cer_factor: float, max_err_factor: float,
                   min_merge_overlap: int) -> str:
    return (
        "=" * 80 + "\n"
        "QUALITY REPORT — speconsense-summarize\n"
        + "=" * 80 + "\n"
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"Source:    {source_folder}\n"
        f"Summary:   {summary_folder}\n"
        f"Filters:   --min-cer-factor={min_cer_factor}  "
        f"--max-err-factor={max_err_factor}  "
        f"--min-merge-overlap={min_merge_overlap}\n"
    )


def _render_executive_summary(
    specimens: Dict[str, _Specimen],
    zero_cluster_specimens: List[Tuple[str, int]],
    counts: Dict[str, int],
) -> str:
    yields = [sp.yield_pct for sp in specimens.values() if sp.yield_pct is not None]
    rics = [c.ric for sp in specimens.values() for c in sp.clusters if c.state == "passed"]
    variants_per_specimen = [
        sum(1 for c in sp.clusters if c.state == "passed")
        for sp in specimens.values()
    ]
    n_zero = len(zero_cluster_specimens)
    zero_reads = sum(t for _, t in zero_cluster_specimens)

    spec_line = f"Specimens:         {len(specimens)} with clusters"
    if n_zero:
        spec_line += (
            f"  ({n_zero} additional specimens had input reads but no clusters,"
            f" {zero_reads} reads total)"
        )

    lines = [
        spec_line,
        f"Total clusters:    {counts['total']}  "
        f"({counts['passed']} passed | {counts['ns']} routed to .ns | "
        f"{counts['lq']} routed to .lq)",
    ]
    if rics:
        lines.append(
            f"RiC (passed):      median {int(_percentile(rics, 0.5))}, "
            f"p25 {int(_percentile(rics, 0.25))}, "
            f"p75 {int(_percentile(rics, 0.75))}, max {max(rics)}"
        )
    if variants_per_specimen:
        lines.append(
            f"Variants/specimen: median {int(_percentile(variants_per_specimen, 0.5))}, "
            f"max {max(variants_per_specimen)}"
        )
    if yields:
        lines.append(
            f"Yield (passed reads / input): "
            f"median {_percentile(yields, 0.5):.1f}%, "
            f"p25 {_percentile(yields, 0.25):.1f}%, "
            f"min {min(yields):.1f}%"
        )
    return _section("EXECUTIVE SUMMARY") + "\n".join(lines) + "\n"


def _render_outlier_table(
    label: str,
    look_for: str,
    result,
    fmt_threshold: Callable[[float], str],
) -> str:
    if result is None:
        return ""
    mean, std, threshold, selected, total = result
    if not selected:
        return (
            f"\n{label}\n  None — all values within normal range "
            f"(run mean ± 2σ: {fmt_threshold(mean)} ± {fmt_threshold(std)}).\n"
        )

    body = [
        f"\n{label}",
        f"  Run mean ± 2σ: {fmt_threshold(mean)} ± {fmt_threshold(std)}  "
        f"(flagging values beyond {fmt_threshold(threshold)}; {total} total)",
        "",
        f"  {'Specimen':<42}  {'Variant':<8}  {'RiC':>4}  "
        f"{'err':>5}  {'rid':>6}  {'ambig':>5}",
    ]
    for c in selected:
        body.append(
            f"  {c.specimen:<42}  {c.short:<8}  {c.ric:>4}  "
            f"{_format_factor(c.err_factor)}  {_format_pct(c.rid)}  "
            f"{c.ambig_count:>5}"
        )
    if total > len(selected):
        body.append(f"  ... ({total - len(selected)} additional outliers not shown)")
    body.append(f"\n  Look for: {look_for}")
    return "\n".join(body) + "\n"


def _render_passed_outliers(specimens: Dict[str, _Specimen]) -> str:
    passed = [c for sp in specimens.values() for c in sp.clusters if c.state == "passed"]

    body = _section("PASSED VARIANTS WORTH INSPECTING")
    body += (
        "These survived CER and err_factor filters but stand out vs the run's own\n"
        "distribution. Inspect MSAs in clusters/cluster_debug/ to judge each case.\n"
    )

    body += _render_outlier_table(
        "Unusually high err_factor (cluster-wide read disagreement)",
        "heterogeneous reads in the MSA — long indels, mixed haplotype signal,"
        " or a true variant that would be split if rerun with stricter --max-err-factor.",
        _select_outliers_high(passed, lambda c: c.err_factor, k=2.0),
        lambda v: f"{v:.2f}",
    )

    body += _render_outlier_table(
        "Unusually high ambiguity (snp + IUPAC in consensus sequence)",
        "over-aggressive cluster-equivalence merging — the .raw* component"
        " sequences may be distinct variants that got conflated. Compare the"
        " .raw component sequences in clusters/cluster_debug/.",
        _select_outliers_high(
            passed,
            lambda c: c.ambig_count if c.ambig_count > 0 else None,
            k=2.0, floor=2,
        ),
        lambda v: f"{int(v):>3}",
    )

    return body


def _render_low_yield_rescues(
    specimens: Dict[str, _Specimen], max_specimens: int = 30,
) -> str:
    yields = [(sp, sp.yield_pct) for sp in specimens.values() if sp.yield_pct is not None]
    if not yields:
        return ""

    yield_vals = [y for _, y in yields]
    p25_yield = _percentile(yield_vals, 0.25)

    low_yield = sorted([(sp, y) for sp, y in yields if y < p25_yield],
                       key=lambda t: t[1])

    rows = []
    for sp, yp in low_yield:
        candidates = []
        for c in sp.clusters:
            if c.state == "ns" and c.cer_factor is not None and c.cer_factor > 0.9:
                candidates.append((c, f"cer={c.cer_factor:.2f}"))
            elif c.state == "lq" and c.err_factor is not None and c.err_factor < 1.7:
                candidates.append((c, f"err={c.err_factor:.2f}"))
        if candidates:
            candidates.sort(
                key=lambda t: -(t[0].cer_factor or 0)
                if t[0].state == "ns" else (t[0].err_factor or 0)
            )
            rows.append((sp, yp, candidates))

    body = _section("POSSIBLE RESCUES (low yield + borderline filtered)")
    body += (
        f"Low yield = below run p25 ({p25_yield:.1f}% of input reads in passing clusters).\n"
        "Borderline filter = .ns with cer_factor > 0.9, or .lq with err_factor < 1.7.\n"
    )

    if not rows:
        body += "\n  None — no low-yield specimens have borderline-filtered variants.\n"
        return body

    shown = rows[:max_specimens]
    body += "\n  " + f"{'Specimen':<42}  {'Yield':>6}  Borderline filtered\n"
    body += "  " + "-" * 88 + "\n"
    for sp, yp, candidates in shown:
        first = True
        for c, reason in candidates:
            spec_col = sp.name if first else ""
            yield_col = f"{yp:>5.1f}%" if first else ""
            tag = f".{c.state}"
            body += (
                f"  {spec_col:<42}  {yield_col:>6}  "
                f"{c.short:<8} {tag}  {reason}  RiC={c.ric}\n"
            )
            first = False
    if len(rows) > max_specimens:
        body += (
            f"  ... ({len(rows) - max_specimens} additional low-yield specimens "
            "not shown)\n"
        )
    body += (
        "\n  To inspect: __Summary__/variants/{specimen}-{variant}.{ns,lq}-RiC*.fasta\n"
        "  To rescue all .ns: rerun summarize with --min-cer-factor 0.\n"
        "  To rescue all .lq: rerun summarize with --max-err-factor 0.\n"
    )
    return body


def _render_overlap_merges(
    overlap_merges: List[OverlapMergeInfo], min_merge_overlap: int,
) -> str:
    # Only include merges that extended beyond full overlap (prefix or suffix > 0).
    true_merges = [m for m in overlap_merges if m.prefix_bp > 0 or m.suffix_bp > 0]
    if not true_merges:
        return ""

    body = _section("OVERLAP MERGE ANALYSIS")
    body += (
        "Cross-primer overlap merges combined variants from different primer pairs\n"
        "(e.g. ITS1 + ITS2). Inspect details below if any overlap is close to the\n"
        "--min-merge-overlap threshold or if length ratios are skewed.\n\n"
    )

    by_specimen: Dict[str, List[OverlapMergeInfo]] = {}
    for m in true_merges:
        by_specimen.setdefault(m.specimen, []).append(m)

    body += f"{len(by_specimen)} specimen(s) had overlap merges:\n\n"
    for specimen in sorted(by_specimen):
        merges = by_specimen[specimen]
        max_iter = max(m.iteration for m in merges)
        suffix = ", iterative" if max_iter > 1 else ""
        body += f"{specimen} ({len(merges)} merge(s){suffix}):\n"
        for m in sorted(merges, key=lambda x: x.iteration):
            iter_prefix = f"  Round {m.iteration}: " if max_iter > 1 else "  "
            input_parts = [
                f"{cluster.rsplit('-', 1)[-1] if '-' in cluster else cluster}"
                f" ({length}bp, RiC={ric})"
                for cluster, length, ric
                in zip(m.input_clusters, m.input_lengths, m.input_rics)
            ]
            body += (
                f"{iter_prefix}Merged: {' + '.join(input_parts)} "
                f"-> {m.output_length}bp\n"
            )
            shorter_len = min(m.input_lengths) if m.input_lengths else 0
            overlap_pct = (m.overlap_bp / shorter_len * 100) if shorter_len > 0 else 0
            body += (
                f"    Overlap: {m.overlap_bp}bp ({overlap_pct:.0f}% of shorter sequence)\n"
                f"    Extensions: prefix={m.prefix_bp}bp, suffix={m.suffix_bp}bp\n"
            )
        body += "\n"

    # Surface only the thin-overlap and skew warnings — these are the
    # "look at this" cases, not the routine ones.
    warnings = []
    for m in true_merges:
        if m.overlap_bp < min_merge_overlap * 1.1:
            shorter_len = min(m.input_lengths) if m.input_lengths else 0
            if shorter_len and m.overlap_bp < shorter_len:
                warnings.append(
                    f"{m.specimen}: thin overlap "
                    f"({m.overlap_bp}bp, threshold={min_merge_overlap}bp)"
                )
        if m.input_lengths:
            max_len = max(m.input_lengths)
            min_len = min(m.input_lengths)
            if max_len > min_len * 3:
                warnings.append(
                    f"{m.specimen}: skewed lengths "
                    f"({max_len}bp / {min_len}bp = {max_len/min_len:.1f}x)"
                )
    if warnings:
        body += "Attention:\n"
        for w in warnings:
            body += f"  * {w}\n"
    return body


def _render_parameter_signals(
    specimens: Dict[str, _Specimen], counts: Dict[str, int],
    max_err_factor: float,
) -> str:
    body = _section("RUN-WIDE PARAMETER SIGNALS")
    body += (
        "Aggregate patterns that may suggest tuning. Each item is a hint, not a\n"
        "directive — interpretation depends on what you expect from your samples.\n"
    )

    passed = [c for sp in specimens.values() for c in sp.clusters if c.state == "passed"]
    n_passed = len(passed)
    if n_passed == 0:
        return body + "\n  No passing clusters; signals not computable.\n"

    notes = []

    high_err = [
        c for c in passed
        if c.err_factor is not None and 1.30 < c.err_factor <= max_err_factor
    ]
    if high_err:
        pct = 100.0 * len(high_err) / n_passed
        notes.append(
            f"  ▸ {len(high_err)} ({pct:.1f}%) passing clusters have err_factor in "
            f"(1.30, {max_err_factor:.2f}] —\n"
            f"    just under the .lq threshold. If samples are expected to be clean,\n"
            f"    --max-err-factor 1.3 would route these to .lq instead."
        )

    border_cer = [
        c for c in passed
        if c.cer_factor is not None and 1.0 < c.cer_factor <= 1.10
    ]
    if border_cer:
        pct = 100.0 * len(border_cer) / n_passed
        notes.append(
            f"  ▸ {len(border_cer)} ({pct:.1f}%) passing clusters have cer_factor in "
            f"(1.00, 1.10] —\n"
            f"    weakly above the .ns threshold. --min-cer-factor 1.1 would route "
            f"these\n    to .ns. Trade-off: false positives vs false negatives in your"
            f" downstream use."
        )

    if counts["total"]:
        ns_pct = 100.0 * counts["ns"] / counts["total"]
        lq_pct = 100.0 * counts["lq"] / counts["total"]
        if ns_pct > 50 or lq_pct > 20:
            notes.append(
                f"  ▸ Filter routing is unusually high "
                f"(.ns: {ns_pct:.1f}%, .lq: {lq_pct:.1f}%).\n"
                f"    Consider whether --min-size or --min-ric should be raised so "
                f"weak\n    clusters drop out earlier in core."
            )

    if not notes:
        body += "\n  No notable patterns at default thresholds.\n"
    else:
        body += "\n" + "\n\n".join(notes) + "\n"

    return body


def _render_pipeline_activity(
    specimens: Dict[str, _Specimen], counts: Dict[str, int],
    overlap_merges: List[OverlapMergeInfo], min_merge_overlap: int,
    final_consensus: List[ConsensusInfo],
) -> str:
    body = _section("PIPELINE ACTIVITY (run-wide)")
    total = counts["total"]
    if total == 0:
        return body + "\n  No clusters.\n"

    rics_passed = [
        c.ric for sp in specimens.values()
        for c in sp.clusters if c.state == "passed"
    ]
    err_passed = [
        c.err_factor for sp in specimens.values()
        for c in sp.clusters if c.state == "passed" and c.err_factor is not None
    ]
    cer_passed = [
        c.cer_factor for sp in specimens.values()
        for c in sp.clusters if c.state == "passed" and c.cer_factor is not None
    ]
    # Merged variants are post-merge outputs that carry snp > 0; this lives on
    # final_consensus, not on the pre-merge source clusters.
    merged = sum(
        1 for cons in final_consensus
        if cons.snp_count is not None and cons.snp_count > 0
    )

    lines = [
        f"  Filter routing:    "
        f"{counts['passed']} passed ({100*counts['passed']/total:.1f}%), "
        f"{counts['ns']} .ns ({100*counts['ns']/total:.1f}%), "
        f"{counts['lq']} .lq ({100*counts['lq']/total:.1f}%)",
        f"  Merged variants:   {merged} (carry snp > 0)",
        f"  Cross-primer overlap merges: {len(overlap_merges)} "
        f"(threshold --min-merge-overlap={min_merge_overlap})",
    ]
    if rics_passed:
        lines.append(
            f"  RiC distribution:  min {min(rics_passed)}, "
            f"p25 {int(_percentile(rics_passed, 0.25))}, "
            f"median {int(_percentile(rics_passed, 0.5))}, "
            f"p75 {int(_percentile(rics_passed, 0.75))}, "
            f"max {max(rics_passed)}"
        )
    if err_passed:
        lines.append(
            f"  err_factor (passed): "
            f"median {_percentile(err_passed, 0.5):.2f}, "
            f"p95 {_percentile(err_passed, 0.95):.2f}, "
            f"max {max(err_passed):.2f}"
        )
    if cer_passed:
        lines.append(
            f"  cer_factor (passed): "
            f"median {_percentile(cer_passed, 0.5):.2f}, "
            f"p5 {_percentile(cer_passed, 0.05):.2f}, "
            f"min {min(cer_passed):.2f}"
        )

    return body + "\n".join(lines) + "\n"


def _render_qctx_calibration(specimens: Dict[str, _Specimen]) -> str:
    """Pooled-ratio calibration check across all clusters.

    Aggregates obs_sum / exp_sum from per-cluster err_factor details (loaded
    from per-specimen metadata JSON). A pooled ratio far from 1.0 is a
    directional signal that the bundled q_ctx model diverges from the
    basecaller's actual error behavior.

    The headline number pools all states (passed + .ns + .lq) to match the
    "all-cluster pooled" methodology used to derive the shipped q_ctx tables
    (see cer_in_practice §8.4 / Appendix B). Excluding .lq would top-truncate
    the right tail by construction and bias the estimator low regardless of
    model quality. The per-state breakdown below is more diagnostic than the
    headline alone — divergence between buckets is itself informative.
    """
    body = _section("q_ctx CALIBRATION CHECK (run-pooled)")

    bucket_order = ("passed", "ns", "lq")
    bucket_label = {"passed": "passed", "ns": ".ns", "lq": ".lq"}
    buckets: Dict[str, Dict[str, float]] = {
        s: {"obs": 0.0, "exp": 0.0, "cols": 0, "n": 0} for s in bucket_order
    }

    for sp in specimens.values():
        for c in sp.clusters:
            if c.state not in buckets:
                continue
            obs = c.err_factor_obs_sum
            exp = c.err_factor_exp_sum
            cols = c.err_factor_cols
            if obs is None or exp is None or cols is None:
                continue
            b = buckets[c.state]
            b["obs"] += obs
            b["exp"] += exp
            b["cols"] += cols
            b["n"] += 1

    obs_total = sum(b["obs"] for b in buckets.values())
    exp_total = sum(b["exp"] for b in buckets.values())
    cols_total = sum(b["cols"] for b in buckets.values())
    n_clusters = sum(b["n"] for b in buckets.values())

    if cols_total < MIN_CALIBRATION_COLS or exp_total <= 0 or n_clusters == 0:
        return (
            body
            + f"\n  Insufficient data for calibration check "
            f"({cols_total} cols across {n_clusters} clusters; "
            f"need >= {MIN_CALIBRATION_COLS} cols).\n"
        )

    pooled = obs_total / exp_total
    lines = [
        f"  Pooled obs/exp:    {pooled:.2f} "
        f"(all states: {n_clusters} clusters, {cols_total} consensus columns)",
    ]
    if pooled < CALIBRATION_LOW:
        lines.append(
            f"  WARNING: q_ctx model appears to overestimate per-column "
            f"disagreement; pooled ratio {pooled:.2f} < {CALIBRATION_LOW:.2f} "
            f"suggests recorded error rates exceed observed."
        )
    elif pooled > CALIBRATION_HIGH:
        lines.append(
            f"  WARNING: q_ctx model appears to underestimate per-column "
            f"disagreement; pooled ratio {pooled:.2f} > {CALIBRATION_HIGH:.2f} "
            f"suggests cluster heterogeneity beyond MAD reach or a "
            f"basecaller/error-model mismatch."
        )
    else:
        lines.append(
            f"  Within expected range "
            f"[{CALIBRATION_LOW:.2f}, {CALIBRATION_HIGH:.2f}]; "
            f"q_ctx model appears reasonably calibrated for this run."
        )

    lines.append("")
    lines.append("  Per-state breakdown:")
    lines.append(f"    {'state':<8} {'n':>6} {'cols':>10}  obs/exp")
    for state in bucket_order:
        b = buckets[state]
        if b["n"] == 0:
            continue
        sub_pooled = b["obs"] / b["exp"] if b["exp"] > 0 else float("nan")
        lines.append(
            f"    {bucket_label[state]:<8} "
            f"{b['n']:>6} {b['cols']:>10}  {sub_pooled:.2f}"
        )

    return body + "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_quality_report(
    final_consensus: List[ConsensusInfo],
    all_raw_consensuses: List[Tuple[ConsensusInfo, str]],
    summary_folder: str,
    source_folder: str,
    overlap_merges: Optional[List[OverlapMergeInfo]] = None,
    min_merge_overlap: int = 200,
    *,
    consensus_list: Optional[List[ConsensusInfo]] = None,
    ns_list: Optional[List[ConsensusInfo]] = None,
    lq_list: Optional[List[ConsensusInfo]] = None,
    min_cer_factor: float = 1.0,
    max_err_factor: float = 1.5,
) -> None:
    """Write the action-oriented quality report.

    Args:
        final_consensus: Post-merge passing variants. Currently unused by the
            new report (kept for backward compatibility); the pre-merge
            ``consensus_list`` carries per-source-cluster metrics.
        all_raw_consensuses: Pre-merge .raw components. Currently unused.
        summary_folder: Output directory for the report.
        source_folder: Source directory containing -all.fasta files and
            cluster_debug/ (with -metadata.json files).
        overlap_merges: Cross-primer overlap merge events from summarize.
        min_merge_overlap: Threshold used for cross-primer overlap merging.
        consensus_list: Pre-merge passing source clusters with full metrics.
            Required for the new report; if omitted, the report falls back
            to ``final_consensus`` (which lacks per-cluster err/cer factors
            for merged variants).
        ns_list: Source clusters routed to .ns by the CER filter.
        lq_list: Source clusters routed to .lq by the err_factor filter.
        min_cer_factor: Threshold used for .ns routing (for header context).
        max_err_factor: Threshold used for .lq routing (for header context).
    """
    if overlap_merges is None:
        overlap_merges = []

    # Use pre-merge consensus_list when supplied; falling back to
    # final_consensus is lossy (merged variants have no err/cer factors)
    # but keeps the function callable with legacy arguments.
    passing = consensus_list if consensus_list is not None else final_consensus
    ns_list = ns_list or []
    lq_list = lq_list or []

    specimens, zero_cluster_specimens = _build_specimens(
        passing, ns_list, lq_list, source_folder
    )
    counts = {
        "passed": len(passing),
        "ns": len(ns_list),
        "lq": len(lq_list),
        "total": len(passing) + len(ns_list) + len(lq_list),
    }

    report = (
        _render_header(source_folder, summary_folder, min_cer_factor,
                       max_err_factor, min_merge_overlap)
        + _render_executive_summary(specimens, zero_cluster_specimens, counts)
        + _render_passed_outliers(specimens)
        + _render_low_yield_rescues(specimens)
        + _render_overlap_merges(overlap_merges, min_merge_overlap)
        + _render_parameter_signals(specimens, counts, max_err_factor)
        + _render_qctx_calibration(specimens)
        + _render_pipeline_activity(specimens, counts, overlap_merges,
                                     min_merge_overlap, final_consensus)
    )

    quality_report_path = os.path.join(summary_folder, "quality_report.txt")
    with open(quality_report_path, "w") as f:
        f.write(report)
    logging.info(f"Quality report written to: {quality_report_path}")
