"""Base classes and protocols for scalability features."""

from abc import ABC, abstractmethod
from typing import Callable, Dict, List, Set, Tuple, Optional, Protocol, runtime_checkable
import logging

from tqdm import tqdm

from .config import ScalabilityConfig


@runtime_checkable
class CandidateFinder(Protocol):
    """Protocol for fast approximate candidate finding.

    Implementations must be able to:
    1. Build an index from a set of sequences
    2. Find candidate matches for query sequences
    3. Clean up resources when done
    """

    @property
    def name(self) -> str:
        """Human-readable name of this backend."""
        ...

    @property
    def is_available(self) -> bool:
        """Check if this backend is available (e.g., tool installed)."""
        ...

    def build_index(self,
                    sequences: Dict[str, str],
                    output_dir: str) -> None:
        """Build search index from sequences.

        Args:
            sequences: Dict mapping sequence_id -> sequence_string
            output_dir: Directory for any cache/index files
        """
        ...

    def find_candidates(self,
                        query_ids: List[str],
                        sequences: Dict[str, str],
                        min_identity: float,
                        max_candidates: int) -> Dict[str, List[str]]:
        """Find candidate matches for query sequences.

        Args:
            query_ids: List of sequence IDs to query
            sequences: Dict mapping sequence_id -> sequence_string
            min_identity: Minimum identity threshold (0.0-1.0)
            max_candidates: Maximum candidates to return per query

        Returns:
            Dict mapping query_id -> list of candidate target_ids
        """
        ...

    def cleanup(self) -> None:
        """Clean up any temporary files or resources."""
        ...


class ScalablePairwiseOperation:
    """Generic scalable pairwise operation using candidate pre-filtering.

    This class encapsulates the two-stage pattern:
    1. Fast candidate finding using CandidateFinder (e.g., vsearch)
    2. Exact scoring using provided scoring function
    """

    def __init__(self,
                 candidate_finder: Optional[CandidateFinder],
                 scoring_function: Callable[[str, str, str, str], float],
                 config: ScalabilityConfig):
        """Initialize scalable pairwise operation.

        Args:
            candidate_finder: Backend for fast candidate finding (None = brute force only)
            scoring_function: Function(seq1, seq2, id1, id2) -> similarity_score (0.0-1.0)
            config: Scalability configuration
        """
        self.candidate_finder = candidate_finder
        self.scoring_function = scoring_function
        self.config = config

    def compute_top_k_neighbors(self,
                                 sequences: Dict[str, str],
                                 k: int,
                                 min_identity: float,
                                 output_dir: str,
                                 min_edges_per_node: int = 3) -> Dict[str, List[Tuple[str, float]]]:
        """Compute top-k nearest neighbors for all sequences.

        Args:
            sequences: Dict mapping sequence_id -> sequence_string
            k: Number of neighbors to find per sequence
            min_identity: Minimum identity threshold for neighbors
            output_dir: Directory for temporary files
            min_edges_per_node: Minimum edges to ensure connectivity

        Returns:
            Dict mapping sequence_id -> list of (neighbor_id, similarity) tuples
        """
        n = len(sequences)

        # Decide whether to use scalable or brute-force approach
        use_scalable = (
            self.config.enabled and
            self.candidate_finder is not None and
            self.candidate_finder.is_available and
            n >= self.config.activation_threshold
        )

        if use_scalable:
            return self._compute_knn_scalable(sequences, k, min_identity, output_dir, min_edges_per_node)
        else:
            return self._compute_knn_brute_force(sequences, k, min_identity, min_edges_per_node)

    def _compute_knn_scalable(self,
                               sequences: Dict[str, str],
                               k: int,
                               min_identity: float,
                               output_dir: str,
                               min_edges_per_node: int) -> Dict[str, List[Tuple[str, float]]]:
        """Two-stage scalable K-NN computation."""
        logging.debug(f"Using {self.candidate_finder.name}-based scalable K-NN computation")

        # Build index
        self.candidate_finder.build_index(sequences, output_dir)

        # Find candidates with oversampling and relaxed threshold
        candidate_count = k * self.config.oversampling_factor
        relaxed_threshold = min_identity * self.config.relaxed_identity_factor

        seq_ids = sorted(sequences.keys())
        candidates = self.candidate_finder.find_candidates(
            seq_ids, sequences, relaxed_threshold, candidate_count
        )

        # Refine with exact scoring
        results: Dict[str, List[Tuple[str, float]]] = {}

        with tqdm(total=len(seq_ids), desc="Refining K-NN with exact scoring") as pbar:
            for seq_id in seq_ids:
                seq_candidates = candidates.get(seq_id, [])

                # Score all candidates
                scored = []
                for cand_id in seq_candidates:
                    if cand_id != seq_id:
                        score = self.scoring_function(sequences[seq_id], sequences[cand_id], seq_id, cand_id)
                        scored.append((cand_id, score))

                # Sort by score descending
                scored.sort(key=lambda x: x[1], reverse=True)

                # Take top k meeting threshold
                top_k = [(cid, score) for cid, score in scored[:k] if score >= min_identity]

                # Ensure minimum connectivity
                if len(top_k) < min_edges_per_node and len(scored) >= min_edges_per_node:
                    for cid, score in scored[len(top_k):]:
                        if score >= min_identity * self.config.relaxed_identity_factor:
                            top_k.append((cid, score))
                            if len(top_k) >= min_edges_per_node:
                                break

                results[seq_id] = top_k
                pbar.update(1)

        return results

    def _compute_knn_brute_force(self,
                                  sequences: Dict[str, str],
                                  k: int,
                                  min_identity: float,
                                  min_edges_per_node: int) -> Dict[str, List[Tuple[str, float]]]:
        """Standard O(n^2) brute-force K-NN computation."""
        logging.debug("Using brute-force K-NN computation")

        seq_ids = sorted(sequences.keys())
        n = len(seq_ids)

        # Compute all pairwise similarities
        # IMPORTANT: This matches main branch's asymmetric dict structure exactly.
        # Main branch creates similarities[id1] = {} inside the loop, which overwrites
        # any entries added via setdefault(). The result is that similarities[id1]
        # only contains entries for id2 > id1 (lexically). This affects tie-breaking
        # when selecting top-k neighbors.
        similarities: Dict[str, Dict[str, float]] = {}

        total = (n * (n - 1)) // 2
        with tqdm(total=total, desc="Computing pairwise similarities") as pbar:
            for id1 in seq_ids:
                similarities[id1] = {}
                for id2 in seq_ids:
                    if id1 >= id2:  # Only calculate upper triangle (id2 > id1)
                        continue
                    score = self.scoring_function(sequences[id1], sequences[id2], id1, id2)
                    similarities[id1][id2] = score
                    similarities.setdefault(id2, {})[id1] = score  # Mirror for lookup
                    pbar.update(1)

        # Extract top-k for each sequence
        results: Dict[str, List[Tuple[str, float]]] = {}

        for seq_id in seq_ids:
            neighbors = sorted(
                [(nid, score) for nid, score in similarities[seq_id].items()],
                key=lambda x: x[1],
                reverse=True
            )

            top_k = [(nid, score) for nid, score in neighbors[:k] if score >= min_identity]

            # Ensure minimum connectivity
            if len(top_k) < min_edges_per_node and len(neighbors) >= min_edges_per_node:
                for nid, score in neighbors[k:]:
                    if score >= min_identity * 0.9:
                        top_k.append((nid, score))
                        if len(top_k) >= min_edges_per_node:
                            break

            results[seq_id] = top_k

        return results

    def compute_distance_matrix(self,
                                 sequences: Dict[str, str],
                                 output_dir: str,
                                 min_identity: float = 0.9) -> Dict[Tuple[str, str], float]:
        """Compute pairwise distance matrix (for HAC clustering).

        Args:
            sequences: Dict mapping sequence_id -> sequence_string
            output_dir: Directory for temporary files
            min_identity: Identity threshold for clustering (used to filter candidates)

        Returns:
            Dict mapping (id1, id2) -> distance, symmetric
        """
        n = len(sequences)
        logging.debug(f"compute_distance_matrix called with {n} sequences")

        # For small sets or when scalability disabled, use brute force
        use_scalable = (
            self.config.enabled and
            self.candidate_finder is not None and
            self.candidate_finder.is_available and
            n >= self.config.activation_threshold and
            n > 50  # Only worthwhile for larger sets
        )

        logging.debug(f"use_scalable={use_scalable} (enabled={self.config.enabled}, "
                      f"finder={self.candidate_finder is not None}, "
                      f"available={self.candidate_finder.is_available if self.candidate_finder else 'N/A'}, "
                      f"threshold={self.config.activation_threshold})")

        if use_scalable:
            return self._compute_distance_matrix_scalable(sequences, output_dir, min_identity)
        else:
            return self._compute_distance_matrix_brute_force(sequences)

    def _compute_distance_matrix_brute_force(self,
                                              sequences: Dict[str, str]) -> Dict[Tuple[str, str], float]:
        """Brute-force distance matrix computation."""
        seq_ids = sorted(sequences.keys())
        distances: Dict[Tuple[str, str], float] = {}

        total = (len(seq_ids) * (len(seq_ids) - 1)) // 2
        with tqdm(total=total, desc="Computing pairwise distances") as pbar:
            for i, id1 in enumerate(seq_ids):
                for id2 in seq_ids[i + 1:]:
                    score = self.scoring_function(sequences[id1], sequences[id2], id1, id2)
                    distance = 1.0 - score  # Convert similarity to distance
                    distances[(id1, id2)] = distance
                    distances[(id2, id1)] = distance
                    pbar.update(1)

        return distances

    def _compute_distance_matrix_scalable(self,
                                           sequences: Dict[str, str],
                                           output_dir: str,
                                           min_identity: float) -> Dict[Tuple[str, str], float]:
        """Scalable distance matrix using candidates to reduce comparisons."""
        logging.debug(f"Using {self.candidate_finder.name}-based scalable distance matrix")

        # Build index
        self.candidate_finder.build_index(sequences, output_dir)

        seq_ids = sorted(sequences.keys())
        n = len(seq_ids)

        # Use same safety factors as K-NN computation
        relaxed_threshold = min_identity * self.config.relaxed_identity_factor
        max_candidates = 500

        logging.debug(f"Finding candidates: identity>={relaxed_threshold:.2f}, max_candidates={max_candidates}")
        all_candidates = self.candidate_finder.find_candidates(
            seq_ids, sequences, relaxed_threshold, max_candidates
        )

        distances: Dict[Tuple[str, str], float] = {}
        computed_pairs: set = set()

        with tqdm(total=len(seq_ids), desc="Computing distances for candidates") as pbar:
            for id1 in seq_ids:
                for id2 in all_candidates.get(id1, []):
                    pair = (min(id1, id2), max(id1, id2))
                    if pair not in computed_pairs and id1 != id2:
                        score = self.scoring_function(sequences[id1], sequences[id2], id1, id2)
                        distance = 1.0 - score
                        distances[(id1, id2)] = distance
                        distances[(id2, id1)] = distance
                        computed_pairs.add(pair)
                pbar.update(1)

        # Return sparse matrix - missing pairs are treated as distance 1.0 by consumers
        logging.debug(f"Computed {len(computed_pairs)} distance pairs (sparse matrix)")

        return distances

    def compute_equivalence_groups(self,
                                    sequences: Dict[str, str],
                                    equivalence_fn: Callable[[str, str], bool],
                                    output_dir: str,
                                    min_candidate_identity: float = 0.95) -> List[List[str]]:
        """Compute groups of equivalent sequences using candidate pre-filtering.

        This is useful for merging clusters whose consensus sequences are
        identical or homopolymer-equivalent. Uses union-find for transitive grouping.

        Args:
            sequences: Dict mapping sequence_id -> sequence_string
            equivalence_fn: Function(seq1, seq2) -> bool for exact equivalence check
            output_dir: Directory for temporary files
            min_candidate_identity: Min identity threshold for candidates (default 0.95)

        Returns:
            List of groups, where each group is a list of equivalent sequence IDs
        """
        n = len(sequences)

        # For small sets or when scalability disabled, use brute force
        use_scalable = (
            self.config.enabled and
            self.candidate_finder is not None and
            self.candidate_finder.is_available and
            n >= self.config.activation_threshold and
            n > 50  # Only worthwhile for larger sets
        )

        if use_scalable:
            return self._compute_equivalence_groups_scalable(
                sequences, equivalence_fn, output_dir, min_candidate_identity
            )
        else:
            return self._compute_equivalence_groups_brute_force(sequences, equivalence_fn)

    def _compute_equivalence_groups_brute_force(self,
                                                 sequences: Dict[str, str],
                                                 equivalence_fn: Callable[[str, str], bool]) -> List[List[str]]:
        """Brute-force O(n²) equivalence group computation."""
        seq_ids = sorted(sequences.keys())

        # Union-find data structure
        parent = {sid: sid for sid in seq_ids}

        def find(x: str) -> str:
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x: str, y: str) -> None:
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        # Check all pairs
        total = (len(seq_ids) * (len(seq_ids) - 1)) // 2
        with tqdm(total=total, desc="Finding equivalent pairs") as pbar:
            for i, id1 in enumerate(seq_ids):
                for id2 in seq_ids[i + 1:]:
                    if equivalence_fn(sequences[id1], sequences[id2]):
                        union(id1, id2)
                    pbar.update(1)

        # Collect groups
        groups: Dict[str, List[str]] = {}
        for sid in seq_ids:
            root = find(sid)
            if root not in groups:
                groups[root] = []
            groups[root].append(sid)

        return list(groups.values())

    def _compute_equivalence_groups_scalable(self,
                                              sequences: Dict[str, str],
                                              equivalence_fn: Callable[[str, str], bool],
                                              output_dir: str,
                                              min_candidate_identity: float) -> List[List[str]]:
        """Scalable equivalence group computation using candidate pre-filtering."""
        logging.debug(f"Using {self.candidate_finder.name}-based equivalence grouping")

        # Build index
        self.candidate_finder.build_index(sequences, output_dir)

        seq_ids = sorted(sequences.keys())
        n = len(seq_ids)

        # Find candidates with high identity (likely equivalent sequences)
        # Use a reasonable max_candidates to limit work while still finding all equivalents
        max_candidates = min(n, 100)
        all_candidates = self.candidate_finder.find_candidates(
            seq_ids, sequences, min_candidate_identity, max_candidates
        )

        # Union-find data structure
        parent = {sid: sid for sid in seq_ids}

        def find(x: str) -> str:
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x: str, y: str) -> None:
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        # Check only candidate pairs
        checked_pairs: set = set()
        equivalent_count = 0

        with tqdm(total=len(seq_ids), desc="Checking candidate equivalences") as pbar:
            for id1 in seq_ids:
                for id2 in all_candidates.get(id1, []):
                    pair = (min(id1, id2), max(id1, id2))
                    if pair not in checked_pairs and id1 != id2:
                        if equivalence_fn(sequences[id1], sequences[id2]):
                            union(id1, id2)
                            equivalent_count += 1
                        checked_pairs.add(pair)
                pbar.update(1)

        logging.debug(f"Found {equivalent_count} equivalent pairs from {len(checked_pairs)} candidates")

        # Collect groups
        groups: Dict[str, List[str]] = {}
        for sid in seq_ids:
            root = find(sid)
            if root not in groups:
                groups[root] = []
            groups[root].append(sid)

        return list(groups.values())
