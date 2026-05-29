"""Command-line interface for speconsense-fit-error-model.

Productizes the offline q_ctx re-estimation procedure (HP paper §8 /
CER paper §4.2 Phase 1 deployment regime) into a user-facing tool that
emits a YAML error model deposited at
``~/.config/speconsense/error_models/{name}.yaml``.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Dict, List, Tuple

from speconsense import qctx
from speconsense._help import add_advanced_argument, install_advanced_help

try:
    from speconsense import __version__
except ImportError:  # pragma: no cover
    __version__ = "dev"

from . import analyze, compare, extract, nonhp, output, selection


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        prog="speconsense-fit-error-model",
        description=(
            "Estimate a custom q_ctx error model from a finished speconsense "
            "output tree and deposit it in ~/.config/speconsense/error_models/. "
            "Productizes the re-estimation procedure documented in Walker 2026b "
            "(HP paper) §8."
        ),
    )
    install_advanced_help(parser)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument(
        "input_dir",
        help="Speconsense output directory (the parent of cluster_debug/).",
    )
    parser.add_argument(
        "--name",
        required=True,
        help=(
            "Output model name. Written to "
            "~/.config/speconsense/error_models/{name}.yaml. The model becomes "
            "selectable via 'speconsense --error-model {name}' on subsequent "
            "runs."
        ),
    )

    sel = parser.add_argument_group("Selection (paper defaults)")
    sel.add_argument(
        "--min-ric", type=int, default=200,
        help=(
            "Minimum RiC for an approach-1 primary anchor to qualify. Paper "
            "default 200 matches the v5.0 retune; relax to 50 for "
            "lower-depth datasets (paper's ont37 setting). (default: 200)"
        ),
    )
    sel.add_argument(
        "--max-err-factor", type=float, default=1.0,
        help=(
            "Maximum primary-anchor err_factor for approach-1 qualification. "
            "Acts as a single-template plausibility gate (HP paper §8.2). "
            "Set 0 to disable. (default: 1.0)"
        ),
    )
    sel.add_argument(
        "--nonhp-min-ric", type=int, default=5,
        help=(
            "Minimum RiC for a cluster to contribute to approach-2 non-HP "
            "rate pooling. (default: 5)"
        ),
    )

    meta = parser.add_argument_group("Model metadata")
    meta.add_argument(
        "--chemistry", default=None,
        help="Chemistry string written to YAML frontmatter (default: inherit "
             "from source model, or 'unknown')."
    )
    meta.add_argument(
        "--basecaller", default=None,
        help="Basecaller string written to YAML frontmatter (default: inherit "
             "from source model, or 'unknown')."
    )
    meta.add_argument(
        "--dataset", default=None,
        help="Dataset identifier written to YAML frontmatter (default: input "
             "dir's basename)."
    )
    meta.add_argument(
        "--source", default=None,
        help="Source/provenance string for YAML frontmatter (default: "
             "'speconsense-fit-error-model {version} ({input-dir})')."
    )
    meta.add_argument(
        "--description", default=None,
        help="Free-text description block prepended to the YAML as a comment."
    )

    cmp_grp = parser.add_argument_group("Comparison")
    cmp_grp.add_argument(
        "--compare-against", default=None,
        help=(
            "Error model to diff against (name or path). When omitted, uses "
            "the most common 'parameters.error_model' across the input's "
            "metadata JSONs, or falls back to "
            f"{qctx.DEFAULT_MODEL_NAME!r}."
        ),
    )

    out = parser.add_argument_group("Output")
    out.add_argument(
        "--debug-dir", default=None,
        help=(
            "If set, write per-context error rate / distribution TSVs and "
            "qualifying-specimen list into this directory. Mirrors the "
            "layout of the original ~/mm/analysis/hp_error_rate output/."
        ),
    )
    out.add_argument(
        "--force", action="store_true",
        help="Overwrite an existing user error model with the same name.",
    )
    out.add_argument(
        "--dry-run", action="store_true",
        help="Compute and print the YAML to stdout instead of writing it.",
    )
    out.add_argument(
        "--threads", type=int, default=1, metavar="N",
        help="Worker processes for parallel approach-1 / approach-2 MSA "
             "scans. (default: 1)",
    )
    out.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )

    adv = parser.add_argument_group(
        "Advanced (pre-tuned — rarely needed)",
        "Hidden from --help; view with --help-advanced.",
    )
    add_advanced_argument(
        adv, "--alpha", type=float, default=0.01,
        help="Bonferroni-corrected binomial significance for the per-position "
             "outlier filter (HP paper §5.1). (default: 0.01)",
    )
    add_advanced_argument(
        adv, "--bimodal-min-frac", type=float, default=0.10,
        help="Minimum second-allele fraction for a position to be flagged as "
             "bimodal (HP paper §5.2). (default: 0.10)",
    )
    add_advanced_argument(
        adv, "--bimodal-rate-ratio", type=float, default=2.0,
        help="Second-allele must exceed N x the context's expected error rate "
             "for the bimodal filter to fire. (default: 2.0)",
    )

    return parser.parse_args(argv)


def setup_logging(log_level: str) -> None:
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    root = logging.getLogger()
    root.setLevel(getattr(logging, log_level))
    root.addHandler(console)


def _build_rates(
    filtered_hp: Dict[Tuple[str, int], analyze.ContextResult],
    pooled_nonhp: extract.NonHPCounts,
) -> Tuple[Dict[str, float], List[str]]:
    """Knit approach-1 and approach-2 results into the canonical YAML keys.

    Per the HP paper §3.5 / §8.5, ``hp-l{N}`` is the *implied per-position
    error rate* ``p`` such that the observed run-level fraction-correct
    satisfies ``(1 - p)^N = frac_correct``. Solving:
    ``p = 1 - frac_correct ** (1 / N)``. This is the form CER and
    ``compute_cluster_err_factor`` consume; it normalizes out the trivial
    effect of run length on aggregate accuracy and lets q_ctx be compared
    on a per-column basis. ``hp-l{N}`` pools across bases at length N
    (``A``, ``C``, ``G``, ``T``) via read-weighted fraction-correct, then
    applies the per-position transform.

    ``non-hp-sub`` and ``non-hp-indel`` come straight from approach 2's
    pooled counts (no per-position transform needed; those are already
    single-position events).

    Returns ``(rates, missing_keys)``.
    """
    rates: Dict[str, float] = {}
    missing: List[str] = []

    sub_rate, indel_rate = nonhp.derive_nonhp_rates(pooled_nonhp)
    if pooled_nonhp.matches + pooled_nonhp.substitutions + pooled_nonhp.deletions == 0:
        missing.extend(["non-hp-sub", "non-hp-indel"])
    else:
        rates["non-hp-sub"] = sub_rate
        rates["non-hp-indel"] = indel_rate

    for length in range(1, 6):
        total_obs = 0
        weighted_correct = 0.0
        for base in ("A", "C", "G", "T"):
            ctx = filtered_hp.get((base, length))
            if ctx is None or ctx.n_obs == 0:
                continue
            total_obs += ctx.n_obs
            weighted_correct += ctx.frac_correct * ctx.n_obs
        if total_obs == 0:
            missing.append(f"hp-l{length}")
            continue
        pooled_correct = weighted_correct / total_obs
        # Implied per-position error rate: p such that (1-p)^L = pooled_correct.
        if pooled_correct <= 0:
            rates[f"hp-l{length}"] = 1.0
        else:
            rates[f"hp-l{length}"] = 1.0 - pooled_correct ** (1.0 / length)

    return rates, missing


def main(argv=None) -> int:
    args = parse_arguments(argv)
    setup_logging(args.log_level)

    input_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(input_dir):
        logging.error(f"Input directory not found: {input_dir}")
        return 2

    # Step 1-2: index + select
    logging.info(f"Indexing specimens under {input_dir}")
    records = selection.index_specimens(input_dir)
    logging.info(f"Found {len(records)} specimens with parseable metadata")
    if not records:
        logging.error("No qualifying metadata found; cannot proceed.")
        return 2

    effective_max_ef = args.max_err_factor if args.max_err_factor > 0 else float("inf")
    qualifying = selection.select_qualifying(records, args.min_ric, effective_max_ef)
    logging.info(
        f"Approach-1 qualifying specimens: {len(qualifying)} "
        f"(min_ric={args.min_ric}, max_err_factor={args.max_err_factor})"
    )
    if not qualifying:
        logging.error(
            "No specimens met the selection criterion. Try lowering --min-ric "
            "or raising --max-err-factor. The ont37 paper retune used "
            "--min-ric 50 for similar reasons."
        )
        return 2

    # Step 3-4: approach 1 extraction + raw aggregation
    observations, calibration_nonhp, n_skipped_a1 = extract.extract_all(
        qualifying, input_dir, threads=args.threads,
    )
    logging.info(
        f"Approach 1: {len(observations)} HP-run observations from "
        f"{len(qualifying) - n_skipped_a1} specimens "
        f"({n_skipped_a1} skipped, missing MSA)"
    )
    if not observations:
        logging.error("No HP observations extracted. Aborting.")
        return 3

    # Step 5: filter biological variants
    raw_hp, filtered_hp, n_excluded = analyze.run_filter_pass(
        observations,
        alpha=args.alpha,
        bimodal_min_frac=args.bimodal_min_frac,
        bimodal_rate_ratio=args.bimodal_rate_ratio,
    )
    logging.info(
        f"Biological-variant filter: excluded {n_excluded} HP runs "
        f"({100*n_excluded/len(observations):.2f}%); "
        f"{len(observations) - n_excluded} remaining"
    )

    # Step 6: approach 2 non-HP rates
    pooled_nonhp, n_clusters_a2, n_skipped_a2, nonhp_by_bin = nonhp.aggregate_nonhp(
        input_dir, nonhp_min_ric=args.nonhp_min_ric, threads=args.threads,
    )
    logging.info(
        f"Approach 2: pooled non-HP rates across {n_clusters_a2} clusters "
        f"(RiC>={args.nonhp_min_ric}; {n_skipped_a2} skipped)"
    )

    # Step 7: synthesize final rates
    rates, missing = _build_rates(filtered_hp, pooled_nonhp)
    if missing:
        logging.error(
            f"Insufficient data for keys: {', '.join(missing)}. "
            "Cannot emit a complete model — try a lower --min-ric or "
            "--nonhp-min-ric for more coverage."
        )
        return 3

    # Step 8: compare against source
    dominant = selection.dominant_error_model(records)
    source_name, source_rates = compare.resolve_source_model(
        args.compare_against, dominant,
    )
    if args.compare_against:
        logging.info(f"Comparing against (explicit): {source_name}")
    elif dominant:
        logging.info(f"Comparing against (run dominant): {source_name}")
    else:
        logging.info(f"Comparing against (fallback): {source_name}")

    # Step 9: write YAML
    chemistry = args.chemistry or _meta_from_source(source_name, "chemistry") or "unknown"
    basecaller = args.basecaller or _meta_from_source(source_name, "basecaller") or "unknown"
    dataset = args.dataset or os.path.basename(input_dir.rstrip("/")) or "unknown"
    source_str = args.source or (
        f"speconsense-fit-error-model {__version__} ({input_dir})"
    )

    try:
        out_path, yaml_text = output.write_model(
            args.name, rates,
            chemistry=chemistry,
            basecaller=basecaller,
            source=source_str,
            dataset=dataset,
            description=args.description,
            force=args.force,
            dry_run=args.dry_run,
        )
    except FileExistsError as e:
        logging.error(str(e))
        return 4

    if args.dry_run:
        print("# Dry run — would write to:")
        print(f"# {out_path}")
        print()
        print(yaml_text)
    else:
        logging.info(f"Wrote model: {out_path}")

    # Step 10: final summary + comparison
    print("")
    compare.print_diff(source_name, source_rates, rates)
    print("")
    print(f"Qualifying specimens:      {len(qualifying)}")
    print(f"HP run observations:       {len(observations)} ({n_excluded} filtered)")
    nonhp_obs = pooled_nonhp.matches + pooled_nonhp.substitutions + pooled_nonhp.deletions
    print(f"Non-HP position obs (A2):  {nonhp_obs:,} across {n_clusters_a2} clusters")
    if not args.dry_run:
        print(f"Output:                    {out_path}")
        print(f"Use with:                  speconsense --error-model {args.name} ...")
    if args.debug_dir:
        debug_path = output.write_debug_outputs(
            args.debug_dir, qualifying, raw_hp, filtered_hp,
            pooled_nonhp, nonhp_by_bin,
        )
        logging.info(f"Debug TSVs written to {debug_path}")

    return 0


def _meta_from_source(name: str, field: str) -> str:
    """Best-effort inheritance of chemistry/basecaller from the source YAML."""
    try:
        path = qctx.get_user_path(name) or qctx.get_bundled_path(name)
    except Exception:  # noqa: BLE001
        return ""
    if path is None:
        return ""
    try:
        meta = qctx._load_metadata(path)  # internal helper, returns dict
    except Exception:  # noqa: BLE001
        return ""
    return meta.get(field, "")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
