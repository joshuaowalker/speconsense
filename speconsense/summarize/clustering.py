"""HAC clustering and variant selection for speconsense-summarize.

Provides hierarchical agglomerative clustering to separate specimens from variants
and variant selection strategies.
"""

import itertools
import logging
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict

from tqdm import tqdm

from speconsense.types import ConsensusInfo
from speconsense.scalability import (
    VsearchCandidateFinder,
    ScalablePairwiseOperation,
    ScalabilityConfig,
)

from .iupac import (
    primers_are_same,
    calculate_adjusted_identity_distance,
    calculate_overlap_aware_distance,
    create_variant_summary,
)


def perform_hac_clustering(consensus_list: List[ConsensusInfo],
                          variant_group_identity: float,
                          min_overlap_bp: int = 0,
                          scalability_config: Optional[ScalabilityConfig] = None,
                          output_dir: str = ".") -> Dict[int, List[ConsensusInfo]]:
    """
    Perform Hierarchical Agglomerative Clustering.
    Separates specimens from variants based on identity threshold.
    Returns groups of consensus sequences.

    Linkage strategy:
    - When min_overlap_bp > 0 (overlap mode): Uses SINGLE linkage, which groups
      sequences if ANY pair has sufficient overlap. This allows prefix+suffix+full
      scenarios where partials only overlap with the full sequence, not each other.
    - When min_overlap_bp == 0 (standard mode): Uses COMPLETE linkage, which
      requires ALL pairs to be within threshold. More conservative for same-length
      sequences.

    When min_overlap_bp > 0, also uses overlap-aware distance calculation that
    allows sequences of different lengths to be grouped together if they
    share sufficient overlap with good identity.
    """
    if len(consensus_list) <= 1:
        return {0: consensus_list}

    # Determine linkage strategy based on overlap mode
    use_single_linkage = min_overlap_bp > 0
    linkage_type = "single" if use_single_linkage else "complete"

    if min_overlap_bp > 0:
        logging.debug(f"Performing HAC clustering with {variant_group_identity} identity threshold "
                     f"({linkage_type} linkage, overlap-aware mode, min_overlap={min_overlap_bp}bp)")
    else:
        logging.debug(f"Performing HAC clustering with {variant_group_identity} identity threshold "
                     f"({linkage_type} linkage)")

    n = len(consensus_list)
    logging.debug(f"perform_hac_clustering: {n} sequences, threshold={variant_group_identity}")
    distance_threshold = 1.0 - variant_group_identity

    # Initialize each sequence as its own cluster
    clusters = [[i] for i in range(n)]

    # Build initial distance matrix between individual sequences
    seq_distances = {}

    # Use scalability if enabled and we have enough sequences
    use_scalable = (
        scalability_config is not None and
        scalability_config.enabled and
        n >= scalability_config.activation_threshold and
        n > 50
    )
    logging.debug(f"perform_hac_clustering: use_scalable={use_scalable}")

    if use_scalable:
        # Build sequence dict with index keys
        sequences = {str(i): consensus_list[i].sequence for i in range(n)}

        # Build primers lookup by ID for the scoring function
        primers_lookup = {str(i): consensus_list[i].primers for i in range(n)}

        # Create scoring function that returns similarity (1 - distance)
        # Use overlap-aware distance when min_overlap_bp > 0
        if min_overlap_bp > 0:
            def score_func(seq1: str, seq2: str, id1: str, id2: str) -> float:
                # Check if primers match - same primers require global distance
                p1, p2 = primers_lookup.get(id1), primers_lookup.get(id2)
                if primers_are_same(p1, p2):
                    # Same primers -> global distance (no overlap merging)
                    return 1.0 - calculate_adjusted_identity_distance(seq1, seq2)
                else:
                    # Different primers -> overlap-aware distance
                    return 1.0 - calculate_overlap_aware_distance(seq1, seq2, min_overlap_bp)
        else:
            def score_func(seq1: str, seq2: str, id1: str, id2: str) -> float:
                return 1.0 - calculate_adjusted_identity_distance(seq1, seq2)

        candidate_finder = VsearchCandidateFinder(num_threads=scalability_config.max_threads)
        if candidate_finder.is_available:
            try:
                operation = ScalablePairwiseOperation(
                    candidate_finder=candidate_finder,
                    scoring_function=score_func,
                    config=scalability_config
                )
                distances = operation.compute_distance_matrix(sequences, output_dir, variant_group_identity)

                # Convert to integer-keyed distances
                for (id1, id2), dist in distances.items():
                    i, j = int(id1), int(id2)
                    seq_distances[(i, j)] = dist
                    seq_distances[(j, i)] = dist
            finally:
                candidate_finder.cleanup()
        else:
            logging.warning("Scalability enabled but vsearch not available. Using brute-force.")
            use_scalable = False

    if not use_scalable:
        # Standard brute-force calculation
        for i, j in itertools.combinations(range(n), 2):
            if min_overlap_bp > 0:
                # Check if primers match - same primers require global distance
                p1, p2 = consensus_list[i].primers, consensus_list[j].primers
                if primers_are_same(p1, p2):
                    # Same primers -> global distance (no overlap merging)
                    dist = calculate_adjusted_identity_distance(
                        consensus_list[i].sequence,
                        consensus_list[j].sequence
                    )
                else:
                    # Different primers -> overlap-aware distance for primer pool scenarios
                    dist = calculate_overlap_aware_distance(
                        consensus_list[i].sequence,
                        consensus_list[j].sequence,
                        min_overlap_bp
                    )
            else:
                # Use standard global distance
                dist = calculate_adjusted_identity_distance(
                    consensus_list[i].sequence,
                    consensus_list[j].sequence
                )
            seq_distances[(i, j)] = dist
            seq_distances[(j, i)] = dist

    # Build sequence adjacency from computed distances (works for both paths)
    # Only include edges where distance < 1.0 (excludes failed alignments and non-candidates)
    seq_adjacency: Dict[int, Set[int]] = defaultdict(set)
    for (i, j), dist in seq_distances.items():
        if dist < 1.0 and i != j:
            seq_adjacency[i].add(j)
            seq_adjacency[j].add(i)

    logging.debug(f"Built adjacency: {len(seq_adjacency)} sequences with edges, "
                  f"{sum(len(v) for v in seq_adjacency.values()) // 2} unique edges")

    # Union-find helper functions
    parent: Dict[int, int] = {i: i for i in range(n)}

    def find(x: int) -> int:
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]

    def union(x: int, y: int) -> None:
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    if use_single_linkage:
        # Single linkage = connected components on edges where dist < threshold
        # This is O(n + E) instead of O(merges * E)
        logging.debug("Single linkage: computing connected components on threshold-filtered edges")

        edge_count = 0
        for (i, j), dist in seq_distances.items():
            if i < j and dist < distance_threshold:
                union(i, j)
                edge_count += 1

        logging.debug(f"Processed {edge_count} edges below threshold {distance_threshold:.3f}")

        # Collect groups
        component_groups: Dict[int, List[int]] = defaultdict(list)
        for i in range(n):
            component_groups[find(i)].append(i)

        clusters = list(component_groups.values())
        logging.info(f"Found {len(clusters)} sequence groups")

    else:
        # Complete linkage: partition by connected components first
        # Clusters from different components can never merge (missing edge = dist 1.0)
        logging.debug("Complete linkage: partitioning into connected components")

        for i in range(n):
            for j in seq_adjacency[i]:
                if i < j:
                    union(i, j)

        # Group sequences by component
        components: Dict[int, List[int]] = defaultdict(list)
        for i in range(n):
            components[find(i)].append(i)

        # Count singletons vs multi-sequence components
        singletons = sum(1 for c in components.values() if len(c) == 1)
        multi_seq = len(components) - singletons
        logging.info(f"Found {len(components)} sequence groups "
                     f"({singletons} single-sequence, {multi_seq} multi-sequence)")

        # Run HAC within each component
        clusters: List[List[int]] = []

        for component_seqs in tqdm(components.values(), desc="HAC per component"):
            if len(component_seqs) == 1:
                clusters.append(component_seqs)
                continue

            # Convert to set for O(1) membership lookup
            component_set = set(component_seqs)

            # Build local adjacency for this component
            local_adjacency: Dict[int, Set[int]] = defaultdict(set)
            for i in component_seqs:
                for j in seq_adjacency[i]:
                    if j in component_set:
                        local_adjacency[i].add(j)

            # Initialize clusters for this component
            seq_to_cluster: Dict[int, int] = {i: i for i in component_seqs}
            cluster_map: Dict[int, List[int]] = {i: [i] for i in component_seqs}

            def get_cluster_adjacency() -> Set[Tuple[int, int]]:
                adjacent_pairs: Set[Tuple[int, int]] = set()
                for seq_i in component_seqs:
                    cluster_i = seq_to_cluster[seq_i]
                    for seq_j in local_adjacency[seq_i]:
                        cluster_j = seq_to_cluster[seq_j]
                        if cluster_i != cluster_j:
                            pair = (min(cluster_i, cluster_j), max(cluster_i, cluster_j))
                            adjacent_pairs.add(pair)
                return adjacent_pairs

            def cluster_distance(cluster1: List[int], cluster2: List[int]) -> float:
                # Complete linkage: max distance, early exit on missing edge or threshold
                max_dist = 0.0
                for i in cluster1:
                    for j in cluster2:
                        if i == j:
                            continue
                        if j not in local_adjacency[i]:
                            return 1.0  # Missing edge = max distance
                        key = (i, j) if (i, j) in seq_distances else (j, i)
                        dist = seq_distances.get(key, 1.0)
                        if dist >= distance_threshold:
                            return 1.0  # Early exit
                        max_dist = max(max_dist, dist)
                return max_dist

            # HAC within component
            while len(cluster_map) > 1:
                adjacent_pairs = get_cluster_adjacency()
                if not adjacent_pairs:
                    break

                min_distance = float('inf')
                merge_pair = None

                for cluster_i, cluster_j in adjacent_pairs:
                    if cluster_i not in cluster_map or cluster_j not in cluster_map:
                        continue
                    dist = cluster_distance(cluster_map[cluster_i], cluster_map[cluster_j])
                    if dist < min_distance:
                        min_distance = dist
                        merge_pair = (cluster_i, cluster_j)

                if min_distance >= distance_threshold or merge_pair is None:
                    break

                ci, cj = merge_pair
                merged = cluster_map[ci] + cluster_map[cj]
                for seq_idx in cluster_map[cj]:
                    seq_to_cluster[seq_idx] = ci
                cluster_map[ci] = merged
                del cluster_map[cj]

            clusters.extend(cluster_map.values())

    # Convert clusters to groups of ConsensusInfo
    groups = {}
    for group_id, cluster_indices in enumerate(clusters):
        group_members = [consensus_list[idx] for idx in cluster_indices]
        groups[group_id] = group_members

    logging.debug(f"HAC clustering created {len(groups)} groups")
    for group_id, group_members in groups.items():
        member_names = [m.sample_name for m in group_members]
        # Convert group_id to final naming (group 0 -> 1, group 1 -> 2, etc.)
        final_group_name = group_id + 1
        logging.debug(f"Group {final_group_name}: {member_names}")

    return groups


def select_variants(group: List[ConsensusInfo],
                   max_variants: int,
                   variant_selection: str,
                   group_number: int = None) -> List[ConsensusInfo]:
    """
    Select variants from a group based on the specified strategy.
    Always includes the largest variant first.
    max_variants of 0 or -1 means no limit (return all variants).

    Logs variant summaries for ALL variants in the group, including those
    that will be skipped in the final output.

    Args:
        group: List of ConsensusInfo to select from
        max_variants: Maximum total variants per group (0 or -1 for no limit)
        variant_selection: Selection strategy ("size" or "diversity")
        group_number: Group number for logging prefix (optional)
    """
    # Sort by size, largest first
    sorted_group = sorted(group, key=lambda x: x.size, reverse=True)

    if not sorted_group:
        return []

    # The primary variant is always the largest
    primary_variant = sorted_group[0]

    # Build prefix for logging
    prefix = f"Group {group_number}: " if group_number is not None else ""

    # Only log Primary when there are other variants to compare against
    if len(sorted_group) > 1:
        logging.info(f"{prefix}Primary: {primary_variant.sample_name} (size={primary_variant.size}, ric={primary_variant.ric})")

    # Handle no limit case (0 or -1 means unlimited)
    if max_variants <= 0:
        selected = sorted_group
    elif len(group) <= max_variants:
        selected = sorted_group
    else:
        # Always include the largest (main) variant
        selected = [primary_variant]
        candidates = sorted_group[1:]

        if variant_selection == "size":
            # Select by size (max_variants - 1 because we already have primary)
            selected.extend(candidates[:max_variants - 1])
        else:  # diversity
            # Select by diversity (maximum distance from already selected)
            while len(selected) < max_variants and candidates:
                best_candidate = None
                best_min_distance = -1

                for candidate in candidates:
                    # Calculate minimum distance to all selected variants
                    min_distance = min(
                        calculate_adjusted_identity_distance(candidate.sequence, sel.sequence)
                        for sel in selected
                    )

                    if min_distance > best_min_distance:
                        best_min_distance = min_distance
                        best_candidate = candidate

                if best_candidate:
                    selected.append(best_candidate)
                    candidates.remove(best_candidate)

    # Now generate variant summaries, showing selected variants first in their final order
    # Then show skipped variants

    # Log selected variants first (excluding primary, which is already logged)
    selected_secondary = selected[1:]  # Exclude primary variant
    for i, variant in enumerate(selected_secondary, 1):
        variant_summary = create_variant_summary(primary_variant.sequence, variant.sequence)
        logging.info(f"{prefix}Variant {i}: (size={variant.size}, ric={variant.ric}) - {variant_summary}")

    # Log skipped variants
    selected_names = {variant.sample_name for variant in selected}
    skipped_variants = [v for v in sorted_group[1:] if v.sample_name not in selected_names]

    for i, variant in enumerate(skipped_variants):
        # Calculate what the variant number would have been in the original sorted order
        original_position = next(j for j, v in enumerate(sorted_group) if v.sample_name == variant.sample_name)
        variant_summary = create_variant_summary(primary_variant.sequence, variant.sequence)
        logging.info(f"{prefix}Variant {original_position}: (size={variant.size}, ric={variant.ric}) - {variant_summary} - skipping")

    return selected
