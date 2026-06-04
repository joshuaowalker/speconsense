"""Grouping and variant selection for speconsense-summarize.

Core is authoritative for within-specimen identity grouping (emitted as
``gid=`` / ``vid=`` header fields). This module bridges that into summarize's
output pipeline:

- ``group_by_core_identity`` buckets consensuses by ``group_rank``.
- ``merge_groups_by_anchor_overlap`` conflates cross-primer groups via
  overlap-aware distance (primer-pool use case).
- ``select_variants`` picks representative variants per group for output.
"""

import itertools
import logging
from typing import List, Dict, Optional
from collections import defaultdict

from tqdm import tqdm

from speconsense.types import ConsensusInfo

from .iupac import (
    primers_are_same,
    calculate_overlap_aware_distance,
    create_variant_summary,
)


def merge_groups_by_anchor_overlap(
    groups: Dict[int, List[ConsensusInfo]],
    min_overlap_bp: int,
    group_identity: float,
    max_iterations: int = 10,
    hp_normalization_length: int = 1,
) -> Dict[int, List[ConsensusInfo]]:
    """Merge cross-primer core-assigned groups whose members overlap well.

    Core groups clusters within a specimen at the `--group-identity` threshold
    via complete linkage. Sequences from different primer pools (e.g., full
    ITS vs ITS2) naturally fall into distinct core groups because length
    mismatch tanks global identity. This function conflates such groups in
    summarize by comparing representative members via the overlap-aware
    distance and merging when any cross-primer pair passes threshold.

    Rules:
    - Only *cross-primer* member pairs are considered (same-primer pairs were
      already evaluated by core). A single cross-primer pair passing the
      overlap threshold triggers a group merge (single linkage over groups),
      which is what enables prefix → full → suffix chains to collapse.
    - The larger-anchor group absorbs the smaller. Absorbed-group members
      append to the survivor's member list; downstream positional renumbering
      handles re-ranking by size.
    - Iterates until no more merges fire.
    """
    if len(groups) < 2 or min_overlap_bp <= 0:
        return groups

    threshold = 1.0 - group_identity
    for _ in range(max_iterations):
        ranks = sorted(groups, key=lambda g: max(v.size for v in groups[g]),
                       reverse=True)
        merged_this_round = False
        for outer_idx, i in enumerate(ranks):
            for j in ranks[outer_idx + 1:]:
                if i not in groups or j not in groups:
                    continue
                merge_pair_dist = _best_cross_primer_overlap(
                    groups[i], groups[j], min_overlap_bp, hp_normalization_length)
                if merge_pair_dist is None or merge_pair_dist >= threshold:
                    continue
                anchor_i = max(groups[i], key=lambda v: v.size)
                anchor_j = max(groups[j], key=lambda v: v.size)
                if anchor_i.size >= anchor_j.size:
                    survivor, absorbed = i, j
                else:
                    survivor, absorbed = j, i
                groups[survivor].extend(groups[absorbed])
                del groups[absorbed]
                merged_this_round = True
                logging.info(
                    f"Cross-primer merge: group {absorbed} absorbed into "
                    f"group {survivor} (best cross-primer distance="
                    f"{merge_pair_dist:.4f})"
                )
                break
            if merged_this_round:
                break
        if not merged_this_round:
            break
    return groups


def _best_cross_primer_overlap(
    group_a: List[ConsensusInfo],
    group_b: List[ConsensusInfo],
    min_overlap_bp: int,
    hp_normalization_length: int = 1,
) -> Optional[float]:
    """Return the minimum overlap-aware distance across cross-primer pairs.

    Only pairs with different primer sets contribute. Returns ``None`` when
    every pair has identical primers (no cross-primer comparison possible).
    """
    best: Optional[float] = None
    for a in group_a:
        for b in group_b:
            if primers_are_same(a.primers, b.primers):
                continue
            dist = calculate_overlap_aware_distance(
                a.sequence, b.sequence, min_overlap_bp, hp_normalization_length)
            if best is None or dist < best:
                best = dist
    return best


def group_by_core_identity(
    consensus_list: List[ConsensusInfo],
) -> Dict[int, List[ConsensusInfo]]:
    """Bucket consensuses by core-assigned ``group_rank`` (gid=) from headers.

    Core is authoritative for within-specimen identity grouping. Every member
    must have ``group_rank`` populated; if any lack it (e.g., FASTA emitted by
    an older core version), raise ``ValueError`` pointing the user to rerun.
    Cross-primer overlap conflation is handled separately downstream.
    """
    if not consensus_list:
        return {}

    missing = [c.sample_name for c in consensus_list if c.group_rank is None]
    if missing:
        example = missing[0]
        raise ValueError(
            f"Input FASTA is missing the gid= header field required by the "
            f"current summarize pipeline (first example: {example!r}). "
            f"Regenerate the consensus with the current version of "
            f"speconsense core, which emits gid=/vid= identity ranks."
        )

    groups: Dict[int, List[ConsensusInfo]] = defaultdict(list)
    for consensus in consensus_list:
        groups[consensus.group_rank].append(consensus)
    return dict(groups)


def select_variants(group: List[ConsensusInfo],
                   max_variants: int,
                   variant_selection: str = "size",
                   group_number: int = None,
                   hp_normalization_length: int = 1) -> List[ConsensusInfo]:
    """
    Select variants from a group by size (largest first).
    max_variants of 0 or -1 means no limit (return all variants).

    Logs variant summaries for ALL variants in the group, including those
    that will be skipped in the final output.
    """
    sorted_group = sorted(group, key=lambda x: x.size, reverse=True)

    if not sorted_group:
        return []

    primary_variant = sorted_group[0]

    prefix = f"Group {group_number}: " if group_number is not None else ""

    if len(sorted_group) > 1:
        logging.info(f"{prefix}Primary: {primary_variant.sample_name} (size={primary_variant.size}, ric={primary_variant.ric})")

    if max_variants <= 0 or len(group) <= max_variants:
        selected = sorted_group
    else:
        selected = sorted_group[:max_variants]

    selected_secondary = selected[1:]
    for i, variant in enumerate(selected_secondary, 1):
        variant_summary = create_variant_summary(primary_variant.sequence, variant.sequence)
        logging.info(f"{prefix}Variant {i}: (size={variant.size}, ric={variant.ric}) - {variant_summary}")

    selected_names = {variant.sample_name for variant in selected}
    skipped_variants = [v for v in sorted_group[1:] if v.sample_name not in selected_names]

    for i, variant in enumerate(skipped_variants):
        original_position = next(j for j, v in enumerate(sorted_group) if v.sample_name == variant.sample_name)
        variant_summary = create_variant_summary(primary_variant.sequence, variant.sequence)
        logging.info(f"{prefix}Variant {original_position}: (size={variant.size}, ric={variant.ric}) - {variant_summary} - skipping")

    return selected
