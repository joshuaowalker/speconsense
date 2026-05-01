"""
Core subpackage for speconsense.

Provides clustering and consensus generation for Oxford Nanopore amplicon reads.
"""

# CLI and entry point
from .cli import main

# Main class
from .clusterer import SpecimenClusterer

# Config classes (for advanced usage)
from .workers import (
    ClusterProcessingConfig,
    ConsensusGenerationConfig,
)

__all__ = [
    # CLI
    "main",
    # Main class
    "SpecimenClusterer",
    # Config classes
    "ClusterProcessingConfig",
    "ConsensusGenerationConfig",
]
