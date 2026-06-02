"""Error models (context-specific error-rate tables) for context-aware CER.

Loads per-basecaller error models from YAML files. Resolution order:

    1. Filesystem path (contains '/' or ends in .yaml/.yml)
    2. User error models in ~/.config/speconsense/error_models/
    3. Bundled error models in speconsense/error_models/

Each model is a flat mapping from a string lookup key (produced by
``qctx_lookup_key``) to a per-position error rate (q_ctx). The lookup keys are:

    non-hp-sub          Substitution at a non-HP position
    non-hp-indel        Indel at a non-HP position
    hp-l{N}             Length change in an HP run of reference length N

HP runs of length beyond those represented in the model return None from
``get_qctx``, signaling that callers should route to blanket HP normalization
rather than CER evaluation.

Phase 2 per-run estimation, described in the CER-in-practice paper §4.2, is
not yet implemented; use the shipped per-basecaller models in the meantime.
"""

from pathlib import Path
from typing import Dict, List, Optional

import yaml

from speconsense.context import ContextClass, ContextTag


_BUNDLED_DIR = Path(__file__).parent / "error_models"
_USER_DIR = Path.home() / ".config" / "speconsense" / "error_models"
_cache: Dict[str, Dict[str, float]] = {}


DEFAULT_MODEL_NAME = "dorado-v5.0"

# Maximum HP length representable in shipped models. HP runs longer than this
# are handled via blanket normalization rather than CER evaluation, matching
# the HP paper §8.1 recommendation. Custom YAML models may extend beyond this
# by adding hp-l{N} keys for larger N; this constant is for reference only.
MAX_HP_LENGTH = 5


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------


def load_table(name_or_path: str) -> Dict[str, float]:
    """Load an error model by name or filesystem path.

    Resolution order:
      1. Anything containing '/' or ending in .yaml/.yml is treated as a path.
      2. A user model in ~/.config/speconsense/error_models/{name}.yaml.
      3. A bundled model in speconsense/error_models/{name}.yaml.

    Returns the ``rates`` mapping from the YAML file. Raises FileNotFoundError
    if the model cannot be resolved and ValueError if the YAML is malformed.
    """
    if name_or_path in _cache:
        return _cache[name_or_path]

    path = _resolve_path(name_or_path)
    with open(path) as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict) or "rates" not in data:
        raise ValueError(f"Invalid error model at {path}: missing 'rates' mapping")
    rates = data["rates"]
    if not isinstance(rates, dict):
        raise ValueError(f"Invalid error model at {path}: 'rates' is not a mapping")

    _cache[name_or_path] = rates
    return rates


def _resolve_path(name_or_path: str) -> Path:
    """Resolve a name-or-path argument to an absolute YAML path."""
    if "/" in name_or_path or name_or_path.endswith((".yaml", ".yml")):
        path = Path(name_or_path)
        if not path.exists():
            raise FileNotFoundError(f"Error model not found: {name_or_path}")
        return path

    user_path = _USER_DIR / f"{name_or_path}.yaml"
    if user_path.exists():
        return user_path

    bundled = _BUNDLED_DIR / f"{name_or_path}.yaml"
    if bundled.exists():
        return bundled

    available_bundled = sorted(p.stem for p in _BUNDLED_DIR.glob("*.yaml"))
    available_user = sorted(p.stem for p in _USER_DIR.glob("*.yaml")) if _USER_DIR.exists() else []
    hint = f"Available shipped models: {', '.join(available_bundled)}"
    if available_user:
        hint += f"; user models: {', '.join(available_user)}"
    raise FileNotFoundError(
        f"Unknown error model '{name_or_path}'. {hint}; "
        f"or pass a filesystem path to a custom YAML."
    )


def list_bundled_models() -> List[str]:
    """Return the names of all shipped error models."""
    return sorted(p.stem for p in _BUNDLED_DIR.glob("*.yaml"))


def list_user_models() -> List[str]:
    """Return the names of user error models in ~/.config/speconsense/error_models/."""
    if not _USER_DIR.exists():
        return []
    return sorted(p.stem for p in _USER_DIR.glob("*.yaml"))


def list_models() -> List[str]:
    """Return the union of bundled and user model names (user takes precedence)."""
    return sorted(set(list_bundled_models()) | set(list_user_models()))


def get_bundled_path(name: str) -> Optional[Path]:
    """Absolute path to the bundled YAML for ``name``, or None."""
    p = _BUNDLED_DIR / f"{name}.yaml"
    return p if p.exists() else None


def get_user_path(name: str) -> Optional[Path]:
    """Absolute path to the user YAML for ``name``, or None."""
    p = _USER_DIR / f"{name}.yaml"
    return p if p.exists() else None


def user_dir() -> Path:
    """Path to the user error-models directory (may not yet exist)."""
    return _USER_DIR


def _load_metadata(path: Path) -> Dict[str, str]:
    """Load top-level metadata fields (name, chemistry, basecaller, source,
    dataset) from an error-model YAML, ignoring the rates block. Returns an
    empty dict on any parse error.
    """
    try:
        with open(path) as f:
            data = yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError):
        return {}
    if not isinstance(data, dict):
        return {}
    meta: Dict[str, str] = {}
    for key in ("name", "chemistry", "basecaller", "source", "dataset"):
        val = data.get(key)
        if isinstance(val, (str, int, float)):
            meta[key] = str(val)
    return meta


def print_models_list() -> None:
    """Print the set of known error models (bundled + user) to stdout.

    Mirrors the UX of ``speconsense.profiles.print_profiles_list``: prints
    each model's source (bundled/user), description line built from YAML
    metadata, and a usage footer.
    """
    names = list_models()
    if not names:
        print("No error models found.")
        print(f"\nUser error models directory: {_USER_DIR}")
        return

    print("Available error models:\n")

    user_names = set(list_user_models())
    for name in names:
        source = "user" if name in user_names else "bundled"
        path = get_user_path(name) if source == "user" else get_bundled_path(name)
        meta = _load_metadata(path) if path else {}
        print(f"  {name} ({source})")
        basecaller = meta.get("basecaller")
        chemistry = meta.get("chemistry")
        if basecaller or chemistry:
            bits = [b for b in (basecaller, chemistry) if b]
            print(f"    {' / '.join(bits)}")
        if meta.get("source"):
            print(f"    Source: {meta['source']}")
        print()

    print("Usage: speconsense --error-model <name> [other options]")
    print(f"User error models directory: {_USER_DIR}")


# Convenience handles for the two shipped defaults. Loaded lazily on access
# to avoid I/O at import time for code paths that supply their own models.
def _shipped(name: str) -> Dict[str, float]:
    return load_table(name)


class _LazyTable:
    """Dict-like lazy loader for shipped models (preserves the previous API)."""
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

# Registry of shipped models keyed by name.
MODELS: Dict[str, _LazyTable] = {
    "dorado-v5.0": DORADO_V5_0,
    "dorado-v3.5": DORADO_V3_5,
}


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------


def qctx_lookup_key(tag: ContextTag) -> str:
    """Derive the q_ctx lookup key for a context tag.

    HP length-change tags collapse to ``hp-l{length}`` because the shipped
    models are pooled across base and direction. The classifier still records
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
    """Look up q_ctx for a context tag against an error model.

    Returns the per-position error rate, or None if the tag's lookup key is
    absent from the model. Callers receiving None should treat the variant as
    ineligible for CER evaluation (route to blanket normalization for long
    HPs, or skip the pair).

    Args:
        tag: Context tag from speconsense.context.
        table: Optional error model dict. Defaults to the shipped
            ``DORADO_V5_0`` model.
    """
    if table is None:
        table = DORADO_V5_0
    try:
        key = qctx_lookup_key(tag)
    except ValueError:
        return None
    return table.get(key)


def is_supported(tag: ContextTag, table: Optional[Dict[str, float]] = None) -> bool:
    """Whether a context tag has a defined q_ctx in the model.

    HP length-change tags beyond the model's coverage return False, signaling
    that the variant should be handled via blanket HP normalization.
    """
    return get_qctx(tag, table) is not None
