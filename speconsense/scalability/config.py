"""Configuration for scalability features."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ScalabilityConfig:
    """Configuration for scalability features.

    Attributes:
        enabled: Whether scalability mode is active
        activation_threshold: Minimum sequence count to activate scalability
        max_threads: Max threads for internal parallelism (default: 1 for backward compatibility)
        backend: Which backend to use (default: 'vsearch')
        oversampling_factor: Multiplier for candidate count in K-NN (default: 10)
        relaxed_identity_factor: Factor to relax identity threshold for candidates (default: 0.9)
        batch_size: Number of sequences per batch for vsearch queries (default: 1000)
    """
    enabled: bool = False
    activation_threshold: int = 0
    max_threads: int = 1
    backend: str = 'vsearch'
    oversampling_factor: int = 10
    relaxed_identity_factor: float = 0.9
    batch_size: int = 1000

    @classmethod
    def from_args(cls, args) -> 'ScalabilityConfig':
        """Create config from command-line arguments.

        The scale_threshold arg controls scalability:
        - 0: disabled
        - N > 0: enabled for datasets >= N sequences (default: 1001)
        """
        threshold = getattr(args, 'scale_threshold', 1001)
        return cls(
            enabled=threshold > 0,
            activation_threshold=threshold,
            max_threads=getattr(args, 'threads', 1),
            backend=getattr(args, 'scalability_backend', 'vsearch'),
        )
