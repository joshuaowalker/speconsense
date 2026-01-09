"""Scalability features for large-scale sequence comparison.

This module provides abstractions for accelerating O(n^2) pairwise sequence
comparison operations using external tools like vsearch.

Example usage:
    from speconsense.scalability import (
        VsearchCandidateFinder,
        ScalablePairwiseOperation,
        ScalabilityConfig
    )

    config = ScalabilityConfig(enabled=True)
    finder = VsearchCandidateFinder()
    operation = ScalablePairwiseOperation(finder, scoring_function, config)

    neighbors = operation.compute_top_k_neighbors(sequences, k=20, min_identity=0.8)
"""

from .config import ScalabilityConfig
from .base import CandidateFinder, ScalablePairwiseOperation
from .vsearch import VsearchCandidateFinder

__all__ = [
    'ScalabilityConfig',
    'CandidateFinder',
    'ScalablePairwiseOperation',
    'VsearchCandidateFinder',
]
