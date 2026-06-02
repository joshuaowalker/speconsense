"""Tests for the variant context classifier (speconsense.context)."""

from speconsense.context import (
    ContextClass,
    ContextTag,
    HPRun,
    build_hp_runs,
    classify_pairwise_differences,
)


# ---------------------------------------------------------------------------
# ContextTag serialization
# ---------------------------------------------------------------------------


def test_tag_string_non_hp_sub():
    assert ContextTag(cls=ContextClass.NON_HP_SUB).to_string() == "non-hp-sub"


def test_tag_string_non_hp_indel():
    assert ContextTag(cls=ContextClass.NON_HP_INDEL).to_string() == "non-hp-indel"


def test_tag_string_hp_length_change_deletion():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="A", length=5, direction=-1)
    assert tag.to_string() == "hp-l5-A-del1"


def test_tag_string_hp_length_change_insertion():
    tag = ContextTag(cls=ContextClass.HP_LENGTH_CHANGE, base="G", length=3, direction=2)
    assert tag.to_string() == "hp-l3-G-ins2"


# ---------------------------------------------------------------------------
# build_hp_runs
# ---------------------------------------------------------------------------


def test_build_hp_runs_no_hp():
    runs = build_hp_runs("ACGTACGT", min_length=2)
    assert all(r is None for r in runs)


def test_build_hp_runs_simple_run():
    runs = build_hp_runs("ACAAAAAGT", min_length=2)
    # Cols 2..6 are AAAAA (length 5)
    for col in range(2, 7):
        assert runs[col] is not None
        assert runs[col].base == "A"
        assert runs[col].length == 5
    assert runs[1] is None  # 'C' not part of run
    assert runs[7] is None  # 'G' not part of run


def test_build_hp_runs_min_length_filter():
    # AA at cols 2..3 — only qualifies when min_length <= 2
    runs2 = build_hp_runs("ACAAGT", min_length=2)
    assert runs2[2] is not None and runs2[2].length == 2
    runs3 = build_hp_runs("ACAAGT", min_length=3)
    assert all(r is None for r in runs3)


def test_build_hp_runs_with_internal_gap():
    # Aligned consensus with an internal gap inside the HP run.
    # AAAAA in MSA cols 0,1,2,4,5 with a gap at col 3.
    aligned = "AAA-AA"
    runs = build_hp_runs(aligned, min_length=2)
    # The gap at col 3 should also be tagged because it's inside the run.
    assert runs[3] is not None
    assert runs[3].base == "A"
    assert runs[3].length == 5


def test_build_hp_runs_extends_into_adjacent_gaps():
    # Aligned consensus with leading gap before an HP run.
    aligned = "-AAAA-"
    runs = build_hp_runs(aligned, min_length=2)
    # Both leading and trailing gaps should be tagged with the run.
    assert runs[0] is not None and runs[0].base == "A"
    assert runs[5] is not None and runs[5].base == "A"


# ---------------------------------------------------------------------------
# classify_pairwise_differences — basic cases
# ---------------------------------------------------------------------------


def test_identical_returns_empty():
    assert classify_pairwise_differences("ACGT", "ACGT") == []


def test_empty_returns_none():
    assert classify_pairwise_differences("", "ACGT") is None
    assert classify_pairwise_differences("ACGT", "") is None


def test_single_substitution_non_hp():
    # GCTACGT vs GCTACGT with one sub at position 3 (T->G)
    candidate = "GCTGCGT"
    reference = "GCTACGT"
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 1
    assert tags[0].cls == ContextClass.NON_HP_SUB


def test_two_substitutions_non_hp():
    candidate = "GCTGCAT"
    reference = "GCTACGT"
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 2
    assert all(t.cls == ContextClass.NON_HP_SUB for t in tags)


# ---------------------------------------------------------------------------
# HP length-change cases
# ---------------------------------------------------------------------------


def test_hp_deletion_one_base():
    # Reference has AAAAA (length 5), candidate is missing one A
    reference = "GCTAAAAAGCT"
    candidate = "GCTAAAAGCT"
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 1
    tag = tags[0]
    assert tag.cls == ContextClass.HP_LENGTH_CHANGE
    assert tag.base == "A"
    assert tag.length == 5
    assert tag.direction == -1


def test_hp_insertion_one_base():
    # Reference has AAA (length 3), candidate has one extra A
    reference = "GCTAAAGCT"
    candidate = "GCTAAAAGCT"
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 1
    tag = tags[0]
    assert tag.cls == ContextClass.HP_LENGTH_CHANGE
    assert tag.base == "A"
    assert tag.length == 3
    assert tag.direction == 1


def test_hp_deletion_two_bases():
    reference = "GCTAAAAAGCT"  # AAAAA
    candidate = "GCTAAAGCT"    # AAA (lost two)
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 1
    tag = tags[0]
    assert tag.cls == ContextClass.HP_LENGTH_CHANGE
    assert tag.base == "A"
    assert tag.length == 5
    assert tag.direction == -2


# ---------------------------------------------------------------------------
# Mixed and edge cases
# ---------------------------------------------------------------------------


def test_mixed_hp_change_and_substitution():
    # HP deletion in AAAAA (unambiguously placed inside the run) plus a
    # substitution far from the HP run. Edlib places the gap at the right
    # edge of the HP run and the sub at the trailing position.
    reference = "ACGTAAAAATCGCAT"  # AAAAA at cols 4..8, T at col 14
    candidate = "ACGTAAAATCGCAG"   # AAAA + ...G at end
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 2
    classes = sorted(t.cls.value for t in tags)
    assert "hp-length-change" in classes
    assert "non-hp-sub" in classes
    # Verify the HP tag has the expected metadata
    hp_tag = next(t for t in tags if t.cls == ContextClass.HP_LENGTH_CHANGE)
    assert hp_tag.base == "A"
    assert hp_tag.length == 5
    assert hp_tag.direction == -1


def test_iupac_match_is_not_a_difference():
    # Reference has Y at one position, candidate has C — IUPAC match.
    reference = "GCTYCGT"
    candidate = "GCTCCGT"
    tags = classify_pairwise_differences(candidate, reference)
    assert tags == []


def test_non_hp_indel_classification():
    # Indel at a non-HP position (single base flanked by different bases).
    reference = "GCTAGCT"
    candidate = "GCTGCT"  # missing the 'A'
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 1
    # 'A' at position 3 is not part of an HP run (length-1 run with min_length=2).
    assert tags[0].cls == ContextClass.NON_HP_INDEL


def test_substitution_within_hp_run_classified_as_non_hp_sub():
    # A substitution inside an A-run that changes one A to G — this is a
    # non-HP-sub by the paper's classification (HP-interior subs treated as
    # non-HP per HP paper §4.3).
    reference = "GCTAAAAAGCT"  # AAAAA
    candidate = "GCTAAGAAGCT"  # one A->G in the middle
    tags = classify_pairwise_differences(candidate, reference)
    assert tags is not None
    assert len(tags) == 1
    assert tags[0].cls == ContextClass.NON_HP_SUB
