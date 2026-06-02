"""Per-specimen variant tree text rendering.

Produces a human-readable hierarchical view of a specimen's variants for
validator review. Each identity group becomes a tree where every non-anchor
variant branches from its closest larger-size peer (highest pairwise
adjusted-identity), with a one-line diff summary versus that parent. Passed,
.ns (CER-filtered), and .lq (err_factor-filtered) variants all appear together
so validators can judge filter calls in context.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import edlib

from speconsense.types import ConsensusInfo
from speconsense.distances import (
    IUPAC_EQUIV,
    calculate_adjusted_identity_distance,
)

from .iupac import create_variant_summary
from .io import strip_cluster_suffix


@dataclass
class _Node:
    info: ConsensusInfo
    status: str  # "passed" | "ns" | "lq"
    short_name: str  # display suffix (e.g. "-1.v2") relative to specimen base
    parent_idx: Optional[int] = None
    parent_identity: Optional[float] = None
    children: List[int] = field(default_factory=list)


def _short_name(sample_name: str, specimen_base: str, status: str) -> str:
    """Return the cluster suffix (e.g. '-1.v2', '-1.v4.ns'), or the full name as fallback.

    Filtered variants get the on-disk status marker appended ('-...v4.ns' /
    '-...v4.lq') so they're visually distinct from any passed variant that may
    happen to share the same vid after summarize's renaming.
    """
    if specimen_base and sample_name.startswith(specimen_base):
        suffix = sample_name[len(specimen_base):] or sample_name
    else:
        suffix = sample_name
    if status in ('ns', 'lq'):
        return f"{suffix}.{status}"
    return suffix


def _tag_status(passed: List[ConsensusInfo],
                ns: List[ConsensusInfo],
                lq: List[ConsensusInfo]) -> List[Tuple[ConsensusInfo, str]]:
    out: List[Tuple[ConsensusInfo, str]] = []
    for c in passed:
        out.append((c, 'passed'))
    for c in ns:
        out.append((c, 'ns'))
    for c in lq:
        out.append((c, 'lq'))
    return out


def _group_by_gid(items: List[Tuple[ConsensusInfo, str]]) -> Dict[int, List[Tuple[ConsensusInfo, str]]]:
    groups: Dict[int, List[Tuple[ConsensusInfo, str]]] = {}
    for c, status in items:
        gid = c.group_rank
        if gid is None:
            continue
        groups.setdefault(gid, []).append((c, status))
    return groups


def _build_tree(nodes: List[_Node], hp_normalization_length: int) -> None:
    """Populate parent_idx/parent_identity/children for each non-anchor node.

    Anchor (largest by size) is index 0 after the caller sorts; for every other
    node we pick the larger-size peer with the highest pairwise IUPAC-aware
    adjusted identity. Tie-break: larger size, then lower index (which already
    encodes size desc).
    """
    n = len(nodes)
    if n <= 1:
        return

    for i in range(1, n):
        cand = nodes[i].info
        best_j: Optional[int] = None
        best_id: float = -1.0
        # Larger-size peers come first since nodes are sorted by size desc.
        for j in range(0, i):
            peer = nodes[j].info
            if peer.size is None or cand.size is None:
                continue
            if peer.size <= cand.size and j != 0:
                # Equal sizes: only the strictly-larger-or-anchor case counts as
                # a valid parent. Anchor (j==0) is always a valid root parent.
                continue
            distance = calculate_adjusted_identity_distance(
                cand.sequence, peer.sequence,
                hp_normalization_length=hp_normalization_length,
            )
            identity = 1.0 - distance
            if identity > best_id:
                best_id = identity
                best_j = j

        if best_j is None:
            # Fall back to the anchor — every group must form a tree rooted at index 0.
            best_j = 0
            best_id = 1.0 - calculate_adjusted_identity_distance(
                cand.sequence, nodes[0].info.sequence,
                hp_normalization_length=hp_normalization_length,
            )

        nodes[i].parent_idx = best_j
        nodes[i].parent_identity = best_id
        nodes[best_j].children.append(i)


def _diff_vs_parent(child: ConsensusInfo, parent: ConsensusInfo,
                    parent_identity: Optional[float]) -> str:
    """One-line edit summary using edlib + create_variant_summary."""
    try:
        edit = edlib.align(parent.sequence, child.sequence,
                           task='distance', additionalEqualities=IUPAC_EQUIV)
        edits = edit.get('editDistance', -1)
    except Exception:
        edits = -1

    summary = create_variant_summary(parent.sequence, child.sequence)
    id_str = f"{parent_identity * 100:.1f}% id" if parent_identity is not None else "id n/a"
    edit_str = f"{edits} edit{'s' if edits != 1 else ''}" if edits >= 0 else "edits n/a"
    return f"{edit_str} ({id_str}) — {summary}"


def _format_metrics(c: ConsensusInfo, status: str) -> str:
    parts: List[str] = []
    parts.append(f"size={c.size}" if c.size is not None else "size=?")
    parts.append(f"ric={c.ric}" if c.ric is not None else "ric=?")
    parts.append(status)
    if c.cer_factor is not None:
        parts.append(f"cer={c.cer_factor:.2f}")
    if c.err_factor is not None:
        parts.append(f"err={c.err_factor:.2f}")
    if c.raw_ric and len(c.raw_ric) > 1:
        parts.append(f"merged({len(c.raw_ric)})")
    return "  ".join(parts)


def _render_group(gid: int, nodes: List[_Node], anchor_full_name: str) -> List[str]:
    counts = {'passed': 0, 'ns': 0, 'lq': 0}
    for node in nodes:
        counts[node.status] = counts.get(node.status, 0) + 1
    breakdown_parts = []
    for key in ('passed', 'ns', 'lq'):
        n = counts.get(key, 0)
        if n > 0:
            breakdown_parts.append(f"{n} {key}")
    breakdown = ", ".join(breakdown_parts) if breakdown_parts else "0"
    plural = "s" if len(nodes) != 1 else ""

    lines: List[str] = []
    lines.append("-" * 80)
    lines.append(f"Group {gid}  —  {len(nodes)} variant{plural}: {breakdown}")
    lines.append("-" * 80)

    if not nodes:
        return lines

    anchor = nodes[0]
    anchor_label = anchor_full_name
    if anchor.status in ('ns', 'lq'):
        anchor_label = f"{anchor_full_name}.{anchor.status}"
    lines.append(f"{anchor_label}  {_format_metrics(anchor.info, anchor.status)}  [anchor]")

    # Render tree DFS with box-drawing characters.
    def render(idx: int, prefix: str, is_last: bool, is_root: bool) -> None:
        if not is_root:
            connector = "└─ " if is_last else "├─ "
            node = nodes[idx]
            line = f"{prefix}{connector}{node.short_name}  {_format_metrics(node.info, node.status)}"
            lines.append(line)
            if node.parent_idx is not None:
                parent = nodes[node.parent_idx]
                cont = "   " if is_last else "│  "
                diff = _diff_vs_parent(node.info, parent.info, node.parent_identity)
                lines.append(f"{prefix}{cont}     vs {parent.short_name}: {diff}")
            new_prefix = prefix + ("   " if is_last else "│  ")
        else:
            new_prefix = prefix

        kids = nodes[idx].children
        # Stable child ordering: larger size first (already encoded by node index).
        kids_sorted = sorted(kids)
        for k_idx, child in enumerate(kids_sorted):
            render(child, new_prefix, k_idx == len(kids_sorted) - 1, is_root=False)

    render(0, "", True, is_root=True)
    return lines


def write_specimen_variant_tree(
    specimen_id: str,
    passed: List[ConsensusInfo],
    ns: List[ConsensusInfo],
    lq: List[ConsensusInfo],
    output_dir: str,
    hp_normalization_length: int = 6,
) -> None:
    """Write a tree-view text file for one specimen at output_dir/{specimen_id}.txt.

    No-ops when there are no eligible variants (no group_rank populated).
    """
    items = _tag_status(passed, ns, lq)
    if not items:
        return

    groups = _group_by_gid(items)
    if not groups:
        return

    # Determine specimen base for short-name rendering. Use the anchor of the
    # first group (any variant works since they all share the specimen prefix).
    any_variant = items[0][0]
    specimen_base = strip_cluster_suffix(any_variant.sample_name)

    total_passed = sum(1 for _, s in items if s == 'passed')
    total_ns = sum(1 for _, s in items if s == 'ns')
    total_lq = sum(1 for _, s in items if s == 'lq')

    header_lines: List[str] = []
    header_lines.append("=" * 80)
    header_lines.append(f"VARIANT TREE — {specimen_id}")
    header_lines.append("=" * 80)
    breakdown = []
    if total_passed:
        breakdown.append(f"{total_passed} passed")
    if total_ns:
        breakdown.append(f"{total_ns} .ns")
    if total_lq:
        breakdown.append(f"{total_lq} .lq")
    plural = "s" if len(items) != 1 else ""
    header_lines.append(f"{len(items)} variant{plural} total: {', '.join(breakdown) if breakdown else '0'}")
    header_lines.append(f"{len(groups)} identity group{'s' if len(groups) != 1 else ''}")
    header_lines.append("")
    header_lines.append("Children branch from the larger-size peer with highest pairwise identity.")
    header_lines.append("Diffs (\"vs ...\") are relative to the immediate parent, not the anchor.")
    header_lines.append("")

    body_lines: List[str] = []
    for gid in sorted(groups.keys()):
        members = groups[gid]
        # Sort by size desc, with stable tie-break on sample_name.
        members.sort(key=lambda t: (-(t[0].size or 0), t[0].sample_name))
        nodes: List[_Node] = [
            _Node(info=c, status=status,
                  short_name=_short_name(c.sample_name, specimen_base, status))
            for c, status in members
        ]
        _build_tree(nodes, hp_normalization_length=hp_normalization_length)
        anchor_full_name = nodes[0].info.sample_name
        body_lines.extend(_render_group(gid, nodes, anchor_full_name))
        body_lines.append("")

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{specimen_id}.txt")
    try:
        with open(out_path, 'w') as fh:
            fh.write("\n".join(header_lines + body_lines))
            if not (header_lines + body_lines)[-1].endswith("\n"):
                fh.write("\n")
    except OSError as e:
        logging.warning(f"Failed to write variant tree for {specimen_id}: {e}")
