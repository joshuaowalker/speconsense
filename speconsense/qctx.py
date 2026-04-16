"""Context-specific error rate (q_ctx) tables for context-aware CER.

Provides empirical per-context error rate values from the homopolymer error
rate analysis (HP paper). Each table is a flat dict from a string lookup key
(produced by qctx_lookup_key) to a per-position error rate.

The Phase 1 implementation ships per-basecaller pooled tables. Phase 2 may
add per-run estimation by replacing or augmenting the shipped tables at
specimen ingestion time.

Conventions:
    * All values are filtered post-clustering per-position error rates,
      conditioned on speconsense's --min-identity 0.90 default.
    * HP length-change rates are pooled across base; per-base granularity
      may be added when supplementary per-(base, length) data is wired up.
    * HP runs of length >= MAX_HP_LENGTH are not represented; callers should
      treat them via blanket normalization (per HP paper §8.1).
"""

from typing import Dict, Optional

from speconsense.context import ContextClass, ContextTag


MAX_HP_LENGTH = 5


# ---------------------------------------------------------------------------
# Default table: Dorado v5.0 / R10.4.1 (derived from ont98 dataset)
# ---------------------------------------------------------------------------

# Filtered per-position rates from the HP error rate analysis
# (Walker 2026b). Values are pooled across bases at each HP length;
# per-base variation is at most ~1.4x in this regime and per-base
# refinement is deferred to a future table version.
DORADO_V5_0: Dict[str, float] = {
    # Non-HP positions
    "non-hp-sub":     0.0059,   # Walker 2026b Table 5 (ont98 substitution)
    "non-hp-indel":   0.0108,   # ont98 deletion + insertion
    # HP length-change rates by reference run length
    "hp-l1":          0.0067,   # Walker 2026b Table 4 (filtered)
    "hp-l2":          0.0083,
    "hp-l3":          0.0097,
    "hp-l4":          0.0099,
    "hp-l5":          0.0113,
}


# ---------------------------------------------------------------------------
# Optional table: Dorado v3.5 (older basecaller) for users on legacy data
# ---------------------------------------------------------------------------

# Filtered per-position rates from the ont37 dataset (Dorado SUP v3.5).
# Approximately 2.2x the v5.0 rates across HP lengths, consistent with the
# basecaller-version effect described in Walker 2026b.
DORADO_V3_5: Dict[str, float] = {
    "non-hp-sub":     0.0144,   # ont37 substitution rate
    "non-hp-indel":   0.0222,   # ont37 deletion + insertion
    "hp-l1":          0.0156,
    "hp-l2":          0.0184,
    "hp-l3":          0.0217,
    "hp-l4":          0.0227,
    "hp-l5":          0.0213,
}


TABLES: Dict[str, Dict[str, float]] = {
    "dorado-v5.0": DORADO_V5_0,
    "dorado-v3.5": DORADO_V3_5,
}

DEFAULT_TABLE_NAME = "dorado-v5.0"


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------


def qctx_lookup_key(tag: ContextTag) -> str:
    """Derive the q_ctx table lookup key for a context tag.

    HP length-change tags collapse to ``hp-l{length}`` because the shipped
    Phase 1 tables are pooled across base. The classifier still records the
    base in ContextTag.base for output annotation; only the lookup is pooled.

    Raises:
        ValueError if the tag's cls is not recognized.
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

    Returns the per-position error rate, or None if the tag is unsupported by
    the table (e.g., HP length >= MAX_HP_LENGTH + 1, or a missing entry).
    Callers receiving None should treat the variant as ineligible for CER
    evaluation (typically: route to blanket normalization for long HPs).

    Args:
        tag: Context tag from speconsense.context.
        table: Optional q_ctx table dict. Defaults to DORADO_V5_0.
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

    HP length-change tags with length > MAX_HP_LENGTH return False, signaling
    that the variant should be treated by blanket HP normalization rather
    than CER evaluation.
    """
    return get_qctx(tag, table) is not None
