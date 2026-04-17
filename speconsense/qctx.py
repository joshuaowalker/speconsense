"""Context-specific error rate (q_ctx) tables for context-aware CER.

Loads per-basecaller error rate tables from YAML files in
``speconsense/qctx_tables/``. Users can select a shipped profile by name or
supply a path to a custom YAML table.

Each table is a flat mapping from a string lookup key (produced by
``qctx_lookup_key``) to a per-position error rate. The lookup keys are:

    non-hp-sub          Substitution at a non-HP position
    non-hp-indel        Indel at a non-HP position
    hp-l{N}             Length change in an HP run of reference length N

HP runs of length beyond those represented in the table return None from
``get_qctx``, signaling that callers should route to blanket HP normalization
rather than CER evaluation.

Phase 2 per-run estimation, described in the CER-in-practice paper §4.2, is
not yet implemented; use the shipped per-basecaller tables in the meantime.
"""

from pathlib import Path
from typing import Dict, Optional

import yaml

from speconsense.context import ContextClass, ContextTag


_TABLES_DIR = Path(__file__).parent / "qctx_tables"
_cache: Dict[str, Dict[str, float]] = {}


DEFAULT_TABLE_NAME = "dorado-v5.0"

# Maximum HP length representable in shipped tables. HP runs longer than this
# are handled via blanket normalization rather than CER evaluation, matching
# the HP paper §8.1 recommendation. Custom YAML tables may extend beyond this
# by adding hp-l{N} keys for larger N; this constant is for reference only.
MAX_HP_LENGTH = 5


# ---------------------------------------------------------------------------
# Table loading
# ---------------------------------------------------------------------------


def load_table(name_or_path: str) -> Dict[str, float]:
    """Load a q_ctx table by shipped name or filesystem path.

    Shipped names (``dorado-v5.0``, ``dorado-v3.5``) resolve to YAML files
    bundled with speconsense. Anything that looks like a path (contains '/' or
    ends in '.yaml' / '.yml') is treated as a user-supplied file.

    Returns the ``rates`` mapping from the YAML file. Raises FileNotFoundError
    if the table cannot be resolved and ValueError if the YAML is malformed.
    """
    if name_or_path in _cache:
        return _cache[name_or_path]

    path = _resolve_path(name_or_path)
    with open(path) as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict) or "rates" not in data:
        raise ValueError(f"Invalid q_ctx table at {path}: missing 'rates' mapping")
    rates = data["rates"]
    if not isinstance(rates, dict):
        raise ValueError(f"Invalid q_ctx table at {path}: 'rates' is not a mapping")

    _cache[name_or_path] = rates
    return rates


def _resolve_path(name_or_path: str) -> Path:
    """Resolve a name-or-path argument to an absolute YAML path."""
    if "/" in name_or_path or name_or_path.endswith((".yaml", ".yml")):
        path = Path(name_or_path)
        if not path.exists():
            raise FileNotFoundError(f"q_ctx table not found: {name_or_path}")
        return path

    shipped = _TABLES_DIR / f"{name_or_path}.yaml"
    if not shipped.exists():
        available = sorted(p.stem for p in _TABLES_DIR.glob("*.yaml"))
        raise FileNotFoundError(
            f"Unknown q_ctx table '{name_or_path}'. "
            f"Available shipped tables: {', '.join(available)}; "
            f"or pass a filesystem path to a custom YAML."
        )
    return shipped


def list_shipped_tables() -> list:
    """Return the names of all shipped q_ctx tables (for CLI help text etc.)."""
    return sorted(p.stem for p in _TABLES_DIR.glob("*.yaml"))


# Convenience handles for the two shipped defaults. Loaded lazily on access
# to avoid I/O at import time for code paths that supply their own tables.
def _shipped(name: str) -> Dict[str, float]:
    return load_table(name)


class _LazyTable:
    """Dict-like lazy loader for shipped tables (preserves the previous API)."""
    def __init__(self, name: str):
        self._name = name
        self._table: Optional[Dict[str, float]] = None

    def _load(self) -> Dict[str, float]:
        if self._table is None:
            self._table = _shipped(self._name)
        return self._table

    def __getitem__(self, key):
        return self._load()[key]

    def __iter__(self):
        return iter(self._load())

    def __len__(self):
        return len(self._load())

    def __contains__(self, key):
        return key in self._load()

    def get(self, key, default=None):
        return self._load().get(key, default)

    def keys(self):
        return self._load().keys()

    def items(self):
        return self._load().items()

    def values(self):
        return self._load().values()


DORADO_V5_0 = _LazyTable("dorado-v5.0")
DORADO_V3_5 = _LazyTable("dorado-v3.5")

# Backward-compat map for code that expects a {name: table} registry.
TABLES: Dict[str, _LazyTable] = {
    "dorado-v5.0": DORADO_V5_0,
    "dorado-v3.5": DORADO_V3_5,
}


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------


def qctx_lookup_key(tag: ContextTag) -> str:
    """Derive the q_ctx table lookup key for a context tag.

    HP length-change tags collapse to ``hp-l{length}`` because the shipped
    tables are pooled across base and direction. The classifier still records
    the base in ``ContextTag.base`` for output annotation; only the lookup is
    pooled.

    Raises:
        ValueError if the tag's cls is not recognized, or if an HP tag is
        missing the required length.
    """
    if tag.cls == ContextClass.NON_HP_SUB:
        return "non-hp-sub"
    if tag.cls == ContextClass.NON_HP_INDEL:
        return "non-hp-indel"
    if tag.cls == ContextClass.HP_LENGTH_CHANGE:
        if tag.length is None:
            raise ValueError("HP_LENGTH_CHANGE tag missing required length")
        return f"hp-l{tag.length}"
    raise ValueError(f"Unknown context class: {tag.cls}")


def get_qctx(tag: ContextTag, table: Optional[Dict[str, float]] = None) -> Optional[float]:
    """Look up q_ctx for a context tag against a q_ctx table.

    Returns the per-position error rate, or None if the tag's lookup key is
    absent from the table. Callers receiving None should treat the variant as
    ineligible for CER evaluation (route to blanket normalization for long
    HPs, or skip the pair).

    Args:
        tag: Context tag from speconsense.context.
        table: Optional q_ctx table dict. Defaults to the shipped
            ``DORADO_V5_0`` table.
    """
    if table is None:
        table = DORADO_V5_0
    try:
        key = qctx_lookup_key(tag)
    except ValueError:
        return None
    return table.get(key)


def is_supported(tag: ContextTag, table: Optional[Dict[str, float]] = None) -> bool:
    """Whether a context tag has a defined q_ctx in the table.

    HP length-change tags beyond the table's coverage return False, signaling
    that the variant should be handled via blanket HP normalization.
    """
    return get_qctx(tag, table) is not None
