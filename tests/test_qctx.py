"""Tests for q_ctx tables and lookup (speconsense.qctx)."""

import pytest

from speconsense.context import ContextClass, ContextTag
from speconsense.qctx import (
    DEFAULT_TABLE_NAME,
    DORADO_V3_5,
    DORADO_V5_0,
    MAX_HP_LENGTH,
    TABLES,
    get_qctx,
    is_supported,
    qctx_lookup_key,
)


# ---------------------------------------------------------------------------
# qctx_lookup_key
# ---------------------------------------------------------------------------


def test_lookup_key_non_hp_sub():
    tag = ContextTag(cls=ContextClass.NON_HP_SUB)
    assert qctx_lookup_key(tag) == "non-hp-sub"


def test_lookup_key_non_hp_indel():
    tag = ContextTag(cls=ContextClass.NON_HP_INDEL)
    assert qctx_lookup_key(tag) == "non-hp-indel"


def test_lookup_key_hp_length_change_pools_base_and_direction():
    # Same lookup key regardless of base or direction — Phase 1 tables are
    # pooled across base and direction.
    tag_a_del = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", length=5, direction=-1)
    tag_g_ins = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="G", length=5, direction=+2)
    assert qctx_lookup_key(tag_a_del) == "hp-l5"
    assert qctx_lookup_key(tag_g_ins) == "hp-l5"


def test_lookup_key_hp_length_required():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", direction=-1)
    with pytest.raises(ValueError):
        qctx_lookup_key(tag)


# ---------------------------------------------------------------------------
# get_qctx
# ---------------------------------------------------------------------------


def test_get_qctx_default_table_non_hp_sub():
    tag = ContextTag(cls=ContextClass.NON_HP_SUB)
    q = get_qctx(tag)
    assert q == DORADO_V5_0["non-hp-sub"]


def test_get_qctx_default_table_hp_length_change():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", length=3, direction=-1)
    q = get_qctx(tag)
    assert q == DORADO_V5_0["hp-l3"]


def test_get_qctx_unsupported_long_hp_returns_none():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", length=MAX_HP_LENGTH + 1, direction=-1)
    assert get_qctx(tag) is None


def test_get_qctx_explicit_v3_5_table():
    tag = ContextTag(cls=ContextClass.NON_HP_SUB)
    q = get_qctx(tag, table=DORADO_V3_5)
    assert q == DORADO_V3_5["non-hp-sub"]
    assert q != DORADO_V5_0["non-hp-sub"]   # the two tables differ


# ---------------------------------------------------------------------------
# is_supported
# ---------------------------------------------------------------------------


def test_is_supported_short_hp():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", length=3, direction=-1)
    assert is_supported(tag) is True


def test_is_supported_long_hp_returns_false():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", length=MAX_HP_LENGTH + 1, direction=-1)
    assert is_supported(tag) is False


def test_is_supported_non_hp():
    assert is_supported(ContextTag(cls=ContextClass.NON_HP_SUB)) is True
    assert is_supported(ContextTag(cls=ContextClass.NON_HP_INDEL)) is True


# ---------------------------------------------------------------------------
# Table sanity
# ---------------------------------------------------------------------------


def test_default_table_name_resolves():
    assert DEFAULT_TABLE_NAME in TABLES
    assert TABLES[DEFAULT_TABLE_NAME] is DORADO_V5_0


def test_v5_0_rates_increase_with_hp_length():
    rates = [DORADO_V5_0[f"hp-l{L}"] for L in range(1, MAX_HP_LENGTH + 1)]
    # Monotonic non-decreasing — empirical observation from the HP paper.
    assert all(rates[i] <= rates[i + 1] for i in range(len(rates) - 1))


def test_v3_5_rates_higher_than_v5_0_at_short_lengths():
    # The basecaller-version effect: v3.5 is roughly 2x v5.0 at short HPs.
    for L in range(1, 5):
        v3 = DORADO_V3_5[f"hp-l{L}"]
        v5 = DORADO_V5_0[f"hp-l{L}"]
        assert v3 > v5


def test_v5_0_non_hp_sub_matches_paper_expectation():
    # The HP paper §3.2 reports 0.59% post-clustering substitution rate.
    assert abs(DORADO_V5_0["non-hp-sub"] - 0.0059) < 1e-9


def test_all_tables_have_required_keys():
    required_keys = (
        ["non-hp-sub", "non-hp-indel"]
        + [f"hp-l{L}" for L in range(1, MAX_HP_LENGTH + 1)]
    )
    for name, table in TABLES.items():
        missing = [k for k in required_keys if k not in table]
        assert not missing, f"Table {name} missing keys: {missing}"
