"""Configuration for scalability features."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ScalabilityConfig:
    """Configuration for scalability features.

    Attributes:
        enabled: Whether scalability mode is active
        activation_threshold: Minimum sequence count to activate scalability (0 = always)
        backend: Which backend to use (default: 'vsearch')
        oversampling_factor: Multiplier for candidate count in K-NN (default: 10)
        relaxed_identity_factor: Factor to relax identity threshold for candidates (default: 0.9)
        batch_size: Number of sequences per batch for vsearch queries (default: 1000)
        recommendation_threshold: Sequence count above which to recommend scalability (default: 2000)
    """
    enabled: bool = False
    activation_threshold: int = 0
    backend: str = 'vsearch'
    oversampling_factor: int = 10
    relaxed_identity_factor: float = 0.9
    batch_size: int = 1000
    recommendation_threshold: int = 2000

    @classmethod
    def from_args(cls, args) -> 'ScalabilityConfig':
        """Create config from command-line arguments.

        The enable_scalability arg can be:
        - None: disabled
        - 0: enabled for all sizes
        - N: enabled for sequences >= N
        """
        threshold = getattr(args, 'enable_scalability', None)
        return cls(
            enabled=threshold is not None,
            activation_threshold=threshold if threshold is not None else 0,
            backend=getattr(args, 'scalability_backend', 'vsearch'),
        )
