"""Configuration for scalability features."""

from dataclasses import dataclass


@dataclass
class ScalabilityConfig:
    """Configuration for scalability features.

    Attributes:
        enabled: Whether scalability mode is active
        backend: Which backend to use (default: 'vsearch')
        oversampling_factor: Multiplier for candidate count in K-NN (default: 10)
        relaxed_identity_factor: Factor to relax identity threshold for candidates (default: 0.9)
        batch_size: Number of sequences per batch for vsearch queries (default: 1000)
        recommendation_threshold: Sequence count above which to recommend scalability (default: 2000)
    """
    enabled: bool = False
    backend: str = 'vsearch'
    oversampling_factor: int = 10
    relaxed_identity_factor: float = 0.9
    batch_size: int = 1000
    recommendation_threshold: int = 2000

    @classmethod
    def from_args(cls, args) -> 'ScalabilityConfig':
        """Create config from command-line arguments."""
        return cls(
            enabled=getattr(args, 'enable_scalability', False),
            backend=getattr(args, 'scalability_backend', 'vsearch'),
        )
