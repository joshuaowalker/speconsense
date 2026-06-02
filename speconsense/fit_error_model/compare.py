"""Source-model loading and diff rendering against the new rate table.

The "source model" is the model the run was clustered under (the most
common ``parameters.error_model`` value in the run's metadata JSONs), or
an explicit ``--compare-against`` override. The diff highlights the eight
canonical keys with absolute values, ratio (new/source), and an arrow
indicator. Ratios outside [0.5, 2.0] trigger a loud warning since they
indicate the dataset is materially different from what the source model
assumes.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

from speconsense import qctx


CANONICAL_KEYS: Tuple[str, ...] = (
    "non-hp-sub",
    "non-hp-indel",
    "hp-l1",
    "hp-l2",
    "hp-l3",
    "hp-l4",
    "hp-l5",
)

WARN_RATIO_LOW = 0.5
WARN_RATIO_HIGH = 2.0


def resolve_source_model(
    explicit: Optional[str],
    dominant: Optional[str],
    fallback: str = qctx.DEFAULT_MODEL_NAME,
) -> Tuple[str, Dict[str, float]]:
    """Pick a comparison model and load its rates.

    Precedence: ``--compare-against`` override → most common
    ``parameters.error_model`` from metadata → ``DEFAULT_MODEL_NAME``.

    Returns ``(source_name, rates_dict)``. Raises ``FileNotFoundError`` if
    the chosen model cannot be resolved (propagated from
    ``qctx.load_table``).
    """
    if explicit:
        return explicit, qctx.load_table(explicit)
    if dominant:
        return dominant, qctx.load_table(dominant)
    return fallback, qctx.load_table(fallback)


def render_diff(
    source_name: str,
    source_rates: Dict[str, float],
    new_rates: Dict[str, float],
) -> Tuple[str, List[str]]:
    """Build the diff table as a string. Returns ``(text, warnings)``.

    ``warnings`` is the list of keys whose ratio falls outside the
    [WARN_RATIO_LOW, WARN_RATIO_HIGH] band — callers may print them as a
    secondary block.
    """
    lines = []
    lines.append(f"Comparison against source model: {source_name}")
    lines.append("")
    header = f"  {'key':<14}  {'source':>10}  {'new':>10}  {'ratio':>8}   delta"
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    warnings: List[str] = []
    for key in CANONICAL_KEYS:
        src = source_rates.get(key)
        new = new_rates.get(key)
        if src is None or new is None:
            lines.append(f"  {key:<14}  {'(missing)':>10}  {'(missing)':>10}")
            continue
        ratio = new / src if src > 0 else float("inf")
        delta = new - src
        if ratio < WARN_RATIO_LOW or ratio > WARN_RATIO_HIGH:
            arrow = "  <<"
            warnings.append(key)
        elif ratio > 1.05:
            arrow = "  ^"
        elif ratio < 0.95:
            arrow = "  v"
        else:
            arrow = ""
        lines.append(
            f"  {key:<14}  {src:>10.4f}  {new:>10.4f}  {ratio:>8.2f}   {delta:+.4f}{arrow}"
        )
    return "\n".join(lines), warnings


def print_diff(
    source_name: str,
    source_rates: Dict[str, float],
    new_rates: Dict[str, float],
) -> None:
    text, warnings = render_diff(source_name, source_rates, new_rates)
    print(text)
    if warnings:
        print("")
        print(
            f"  WARNING: {len(warnings)} key(s) outside [{WARN_RATIO_LOW:.1f}, "
            f"{WARN_RATIO_HIGH:.1f}]x — dataset differs materially from "
            f"{source_name!r}: {', '.join(warnings)}"
        )
        logging.warning(
            f"Source-model comparison flagged keys: {', '.join(warnings)}. "
            "This is expected when comparing across basecaller versions; "
            "investigate if both runs use the same basecaller."
        )
