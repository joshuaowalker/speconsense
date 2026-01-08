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


def _complete_linkage_subset(
    indices: List[int],
    seq_distances: Dict[Tuple[int, int], float],
    distance_threshold: float,
    seq_adjacency: Dict[int, Set[int]]
) -> List[List[int]]:
    """Run complete linkage HAC on a subset of sequences.

    First partitions the subset into connected components (based on seq_adjacency),
    then runs HAC within each component. This matches the behavior of the original
    complete linkage code.

    Args:
        indices: List of sequence indices to cluster
        seq_distances: Precomputed distances between sequence pairs
        distance_threshold: Maximum distance for merging (1.0 - identity)
        seq_adjacency: Adjacency dict showing which sequences have edges

    Returns:
        List of clusters, where each cluster is a list of original indices
    """
    if len(indices) <= 1:
        return [indices]

    component_set = set(indices)

    # First, partition into connected components using union-find
    # This matches the original complete linkage behavior
    parent: Dict[int, int] = {i: i for i in indices}

    def find(x: int) -> int:
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x: int, y: int) -> None:
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in indices:
        for j in seq_adjacency.get(i, set()):
            if j in component_set and i < j:
                union(i, j)

    # Group indices by component
    components: Dict[int, List[int]] = defaultdict(list)
    for i in indices:
        components[find(i)].append(i)

    # Run HAC within each connected component
    all_clusters: List[List[int]] = []
    for component_indices in components.values():
        if len(component_indices) == 1:
            all_clusters.append(component_indices)
            continue

        # Run HAC on this component
        component_clusters = _run_hac_on_component(
            component_indices, seq_distances, distance_threshold, seq_adjacency
        )
        all_clusters.extend(component_clusters)

    return all_clusters


def _run_hac_on_component(
    indices: List[int],
    seq_distances: Dict[Tuple[int, int], float],
    distance_threshold: float,
    seq_adjacency: Dict[int, Set[int]]
) -> List[List[int]]:
    """Run HAC on a single connected component.

    This is the inner HAC loop, separated from component partitioning.
    """
    if len(indices) <= 1:
        return [indices]

    component_set = set(indices)

    # Build local adjacency for this subset
    local_adjacency: Dict[int, Set[int]] = defaultdict(set)
    for i in indices:
        for j in seq_adjacency.get(i, set()):
            if j in component_set:
                local_adjacency[i].add(j)

    # Initialize each sequence as its own cluster
    seq_to_cluster: Dict[int, int] = {i: i for i in indices}
    cluster_map: Dict[int, List[int]] = {i: [i] for i in indices}

    def get_cluster_adjacency() -> Set[Tuple[int, int]]:
        adjacent_pairs: Set[Tuple[int, int]] = set()
        for seq_i in indices:
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

    # HAC merging loop
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

    return list(cluster_map.values())


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
    - When min_overlap_bp > 0 (overlap mode): Uses HYBRID linkage:
      - Phase 1: COMPLETE linkage within each primer set (prevents chaining)
      - Phase 2: SINGLE linkage across primer sets (allows ITS1+full+ITS2 merging)
      This prevents sequences with the same primers from chaining through
      intermediates while still allowing different-primer sequences to merge
      via overlap regions.
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
    use_hybrid_linkage = min_overlap_bp > 0
    linkage_type = "hybrid" if use_hybrid_linkage else "complete"

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

    if use_hybrid_linkage:
        # HYBRID LINKAGE: Complete within primer sets, single across primer sets
        # This prevents chaining within same-primer sequences while allowing
        # different-primer sequences (e.g., ITS1 + full ITS + ITS2) to merge.
        logging.debug("Hybrid linkage: complete within primer sets, single across")

        # Phase 1: Group sequences by primer set
        primer_groups: Dict[Tuple[str, ...], List[int]] = defaultdict(list)
        for i, cons in enumerate(consensus_list):
            primer_key = tuple(sorted(cons.primers)) if cons.primers else ('_none_',)
            primer_groups[primer_key].append(i)

        logging.debug(f"Found {len(primer_groups)} distinct primer sets")

        # Run complete linkage HAC within each primer group
        primer_coherent_clusters: List[Tuple[Tuple[str, ...], List[int]]] = []

        # Log info about the work to be done
        max_group_size = max(len(indices) for indices in primer_groups.values())
        if max_group_size > 1000:
            logging.info(f"Running HAC on {len(primer_groups)} primer groups "
                        f"(largest has {max_group_size} sequences, this may take several minutes)")

        for primer_key, indices in primer_groups.items():
            if len(indices) == 1:
                primer_coherent_clusters.append((primer_key, indices))
            else:
                # Run complete linkage on this primer subset
                sub_clusters = _complete_linkage_subset(
                    indices, seq_distances, distance_threshold, seq_adjacency
                )
                for cluster in sub_clusters:
                    primer_coherent_clusters.append((primer_key, cluster))

        logging.debug(f"Phase 1 complete: {len(primer_coherent_clusters)} primer-coherent clusters")

        # Phase 2: Connect clusters with different primers using single linkage
        # BUT: prevent transitive chaining that would connect same-primer clusters
        # via different-primer intermediates
        n_clusters = len(primer_coherent_clusters)

        # Track which clusters are in each group (list of cluster indices per group)
        groups: List[Set[int]] = [set([i]) for i in range(n_clusters)]
        cluster_to_group: Dict[int, int] = {i: i for i in range(n_clusters)}

        def get_group_primers(group_idx: int) -> Dict[Tuple[str, ...], List[int]]:
            """Get all primer keys and their cluster indices in a group."""
            result: Dict[Tuple[str, ...], List[int]] = defaultdict(list)
            for cluster_idx in groups[group_idx]:
                primer_key = primer_coherent_clusters[cluster_idx][0]
                result[primer_key].append(cluster_idx)
            return result

        def can_merge_groups(group_a: int, group_b: int) -> bool:
            """Check if merging would violate complete linkage for same-primer clusters."""
            # Get all primer->clusters mappings for the merged group
            primers_a = get_group_primers(group_a)
            primers_b = get_group_primers(group_b)

            # Check each primer key that appears in both groups
            for primer_key in primers_a:
                if primer_key in primers_b:
                    # Same primer key in both groups - need complete linkage check
                    clusters_a = [primer_coherent_clusters[i][1] for i in primers_a[primer_key]]
                    clusters_b = [primer_coherent_clusters[i][1] for i in primers_b[primer_key]]

                    # All pairs between clusters_a and clusters_b must satisfy complete linkage
                    for ca in clusters_a:
                        for cb in clusters_b:
                            # Complete linkage: max distance must be < threshold
                            max_dist = 0.0
                            for si in ca:
                                for sj in cb:
                                    dist = seq_distances.get((si, sj), seq_distances.get((sj, si), 1.0))
                                    max_dist = max(max_dist, dist)
                                    if max_dist >= distance_threshold:
                                        return False  # Would violate complete linkage
                            if max_dist >= distance_threshold:
                                return False
            return True

        def merge_groups(group_a: int, group_b: int) -> None:
            """Merge group_b into group_a."""
            if group_a == group_b:
                return
            for cluster_idx in groups[group_b]:
                cluster_to_group[cluster_idx] = group_a
                groups[group_a].add(cluster_idx)
            groups[group_b] = set()

        cross_primer_edges = 0
        cross_primer_blocked = 0
        for i in range(n_clusters):
            for j in range(i + 1, n_clusters):
                primer_i, cluster_i = primer_coherent_clusters[i]
                primer_j, cluster_j = primer_coherent_clusters[j]

                # Skip same-primer pairs (already handled in phase 1)
                if primer_i == primer_j:
                    continue

                # Different primers: check single linkage distance
                min_dist = 1.0
                for si in cluster_i:
                    for sj in cluster_j:
                        dist = seq_distances.get((si, sj), seq_distances.get((sj, si), 1.0))
                        min_dist = min(min_dist, dist)
                        if min_dist < distance_threshold:
                            break
                    if min_dist < distance_threshold:
                        break

                if min_dist < distance_threshold:
                    group_i = cluster_to_group[i]
                    group_j = cluster_to_group[j]
                    if group_i != group_j:
                        # Check if merge would create invalid same-primer connections
                        if can_merge_groups(group_i, group_j):
                            merge_groups(group_i, group_j)
                            cross_primer_edges += 1
                        else:
                            cross_primer_blocked += 1

        logging.debug(f"Phase 2: {cross_primer_edges} cross-primer connections, "
                     f"{cross_primer_blocked} blocked by complete linkage constraint")

        # Collect final groups using the new group structure
        final_groups: Dict[int, List[int]] = defaultdict(list)
        for i, (_, cluster) in enumerate(primer_coherent_clusters):
            group_idx = cluster_to_group[i]
            final_groups[group_idx].extend(cluster)

        clusters = list(final_groups.values())
        logging.info(f"Found {len(clusters)} sequence groups (hybrid linkage)")

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
