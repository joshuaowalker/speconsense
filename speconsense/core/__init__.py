"""
Core subpackage for speconsense.

Provides clustering and consensus generation for Oxford Nanopore amplicon reads.
"""

# CLI and entry point
from .cli import main

# Main class
from .clusterer import SpecimenClusterer

# Worker functions and config classes (for advanced usage)
from .workers import (
    ClusterProcessingConfig,
    ConsensusGenerationConfig,
    _run_spoa_worker,
    _process_cluster_worker,
    _generate_cluster_consensus_worker,
    _trim_primers_standalone,
    _phase_reads_by_variants_standalone,
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
