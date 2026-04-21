"""Tests for error models and q_ctx lookup (speconsense.qctx)."""

import pytest

from speconsense.context import ContextClass, ContextTag
from speconsense.qctx import (
    DEFAULT_MODEL_NAME,
    DORADO_V3_5,
    DORADO_V5_0,
    MAX_HP_LENGTH,
    MODELS,
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


def test_default_model_name_resolves():
    assert DEFAULT_MODEL_NAME in MODELS
    assert MODELS[DEFAULT_MODEL_NAME] is DORADO_V5_0


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


def test_all_models_have_required_keys():
    required_keys = (
        ["non-hp-sub", "non-hp-indel"]
        + [f"hp-l{L}" for L in range(1, MAX_HP_LENGTH + 1)]
    )
    for name, table in MODELS.items():
        missing = [k for k in required_keys if k not in table]
        assert not missing, f"Model {name} missing keys: {missing}"


# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------


def test_list_bundled_models_contains_both_defaults():
    from speconsense.qctx import list_bundled_models
    shipped = list_bundled_models()
    assert "dorado-v5.0" in shipped
    assert "dorado-v3.5" in shipped


def test_load_table_by_name():
    from speconsense.qctx import load_table
    t = load_table("dorado-v5.0")
    assert isinstance(t, dict)
    assert "non-hp-sub" in t
    assert t["non-hp-sub"] == pytest.approx(0.0059)


def test_load_table_unknown_name_raises():
    from speconsense.qctx import load_table
    with pytest.raises(FileNotFoundError):
        load_table("does-not-exist")


def test_load_table_from_custom_path(tmp_path):
    from speconsense.qctx import load_table
    custom = tmp_path / "custom.yaml"
    custom.write_text(
        "name: custom\n"
        "rates:\n"
        "  non-hp-sub: 0.01\n"
        "  non-hp-indel: 0.02\n"
        "  hp-l1: 0.005\n"
        "  hp-l2: 0.007\n"
    )
    table = load_table(str(custom))
    assert table["non-hp-sub"] == pytest.approx(0.01)
    assert table["hp-l2"] == pytest.approx(0.007)
    # Keys not in the table return None via get_qctx
    assert table.get("hp-l5") is None


def test_load_table_malformed_raises(tmp_path):
    from speconsense.qctx import load_table
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: bad\n# no rates block\n")
    with pytest.raises(ValueError):
        load_table(str(bad))
