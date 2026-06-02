"""Qualifying-specimen selection from speconsense metadata JSONs.

Walks a finished speconsense output tree's ``cluster_debug/*-metadata.json``
files, extracts each specimen's primary-anchor RiC and err_factor, and
applies the paper's selection criterion (``RiC >= min_ric`` AND
``err_factor < max_err_factor``) to produce the list of specimens whose
primary cluster is suitable for HP rate estimation (approach 1).

The same index is reused by approach 2 to identify the source error model
the run was clustered under (for the comparison report). Approach 2 itself
operates on every cluster MSA, not just qualifying primaries, so it does
not consume the qualifying list.
"""

from __future__ import annotations

import glob
import json
import logging
import os
from dataclasses import dataclass
from typing import List, Optional


@dataclass(frozen=True)
class SpecimenRecord:
    """One specimen's metadata-derived selection inputs."""

    name: str  # sample_name
    metadata_path: str
    error_model: Optional[str]  # parameters.error_model (None if absent)
    total_input_reads: Optional[int]
    primary_m: Optional[int]  # variants[0].M when group_rank=1,variant_rank=1
    primary_err_factor: Optional[float]
    schema_version: str

    @property
    def qualifies(self) -> bool:
        """True when all selection inputs are present and well-defined."""
        return (
            self.primary_m is not None
            and self.primary_err_factor is not None
        )


def _find_primary_anchor(variants: list) -> Optional[dict]:
    """Return the variant dict with group_rank=1 AND variant_rank=1, if any."""
    for v in variants:
        if v.get("group_rank") == 1 and v.get("variant_rank") == 1:
            return v
    return None


def load_specimen(metadata_path: str) -> SpecimenRecord:
    """Parse one metadata JSON into a SpecimenRecord.

    Hard-fails on schema_version < 2.0 because variant-level ``err_factor``
    and ``M`` were introduced in 2.0; without them we cannot make a
    selection decision. The caller's error message should guide the user to
    re-run with a current speconsense version.
    """
    with open(metadata_path) as f:
        data = json.load(f)

    schema = str(data.get("schema_version", "1.0"))
    if not schema.startswith("2."):
        raise ValueError(
            f"{metadata_path}: schema_version={schema!r} is not supported. "
            "speconsense-fit-error-model requires schema_version 2.0+ "
            "(variant-level err_factor and M). Re-cluster with a recent "
            "speconsense version to regenerate metadata."
        )

    params = data.get("parameters") or {}
    error_model = params.get("error_model") if isinstance(params, dict) else None
    total = data.get("total_input_reads")
    name = data.get("sample_name") or os.path.basename(metadata_path)

    variants = data.get("variants") or []
    primary = _find_primary_anchor(variants)
    if primary is None:
        primary_m = None
        primary_ef = None
    else:
        primary_m = primary.get("M")
        primary_ef = primary.get("err_factor")

    return SpecimenRecord(
        name=name,
        metadata_path=metadata_path,
        error_model=error_model,
        total_input_reads=total,
        primary_m=primary_m,
        primary_err_factor=primary_ef,
        schema_version=schema,
    )


def index_specimens(input_dir: str) -> List[SpecimenRecord]:
    """Load every ``cluster_debug/*-metadata.json`` under input_dir.

    Skips specimens whose metadata cannot be parsed (logging a warning).
    Sorted by sample_name for deterministic downstream ordering.
    """
    debug = os.path.join(input_dir, "cluster_debug")
    if not os.path.isdir(debug):
        raise FileNotFoundError(
            f"No cluster_debug/ directory under {input_dir!r}. "
            "speconsense-fit-error-model expects a finished speconsense output "
            "tree (the parent of cluster_debug/)."
        )

    paths = sorted(glob.glob(os.path.join(debug, "*-metadata.json")))
    records: List[SpecimenRecord] = []
    n_failed = 0
    for path in paths:
        try:
            records.append(load_specimen(path))
        except ValueError:
            # schema-version hard-fail propagates with full guidance.
            raise
        except Exception as e:  # noqa: BLE001
            logging.warning(f"Skipping {os.path.basename(path)}: {e}")
            n_failed += 1
    if n_failed:
        logging.warning(f"{n_failed} metadata JSONs failed to parse")
    return records


def select_qualifying(
    records: List[SpecimenRecord],
    min_ric: int,
    max_err_factor: float,
) -> List[SpecimenRecord]:
    """Apply the approach-1 selection criterion.

    A specimen qualifies when its primary-anchor cluster has
    ``M >= min_ric`` and ``err_factor < max_err_factor``. Mirrors the HP
    paper §8 selection (defaults 200 / 1.0 for ont98; relax to 50 for
    lower-depth ont37).
    """
    qualifying = []
    for r in records:
        if not r.qualifies:
            continue
        if r.primary_m < min_ric:  # type: ignore[operator]
            continue
        if r.primary_err_factor >= max_err_factor:  # type: ignore[operator]
            continue
        qualifying.append(r)
    return qualifying


def dominant_error_model(records: List[SpecimenRecord]) -> Optional[str]:
    """Return the most common ``parameters.error_model`` value across records.

    Used by the comparison step when ``--compare-against`` is unset. Returns
    None when no record reports an error_model.
    """
    counts: dict = {}
    for r in records:
        if r.error_model:
            counts[r.error_model] = counts.get(r.error_model, 0) + 1
    if not counts:
        return None
    return max(counts.items(), key=lambda kv: kv[1])[0]
