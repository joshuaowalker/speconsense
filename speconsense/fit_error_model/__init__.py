"""Fit-error-model subpackage for speconsense.

Re-estimates a per-context error-rate (q_ctx) table from a finished
speconsense output tree and writes a user error model YAML to
``~/.config/speconsense/error_models/``. Productizes the offline
re-estimation procedure documented in the HP paper §8 and the CER paper
§4.2 ("Phase 1 / Phase 2 deployment regimes").

Entry point: ``speconsense-fit-error-model`` (see ``cli.main``).
"""

from .cli import main

__all__ = ["main"]
