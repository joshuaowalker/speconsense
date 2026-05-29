"""Tests for the speconsense-fit-error-model subpackage."""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from speconsense.fit_error_model import (
    analyze,
    compare,
    extract,
    nonhp,
    output,
    selection,
)
from speconsense.fit_error_model.cli import _build_rates


# ---------------------------------------------------------------------------
# Synthetic helpers
# ---------------------------------------------------------------------------


def _make_msa(consensus: str, reads: list[str]) -> str:
    """Render a SPOA-style MSA FASTA string."""
    lines = [">Consensus", consensus]
    for i, r in enumerate(reads):
        lines.append(f">read{i}")
        lines.append(r)
    return "\n".join(lines) + "\n"


def _write_msa(tmp_path: Path, specimen: str, ric: int,
                consensus: str, reads: list[str]) -> Path:
    debug = tmp_path / "cluster_debug"
    debug.mkdir(parents=True, exist_ok=True)
    msa_path = debug / f"{specimen}-1.v1-RiC{ric}-msa.fasta"
    msa_path.write_text(_make_msa(consensus, reads))
    return msa_path


def _write_metadata(
    tmp_path: Path,
    specimen: str,
    *,
    primary_m: int,
    primary_err_factor: float,
    error_model: str = "dorado-v5.0",
    schema_version: str = "2.0",
) -> Path:
    debug = tmp_path / "cluster_debug"
    debug.mkdir(parents=True, exist_ok=True)
    meta = {
        "schema_version": schema_version,
        "version": "test",
        "sample_name": specimen,
        "parameters": {"error_model": error_model},
        "total_input_reads": primary_m * 2,
        "variants": [
            {
                "cluster_id": "1.v1",
                "group_rank": 1,
                "variant_rank": 1,
                "M": primary_m,
                "N": primary_m,
                "err_factor": primary_err_factor,
            },
        ],
    }
    path = debug / f"{specimen}-metadata.json"
    path.write_text(json.dumps(meta))
    return path


# ---------------------------------------------------------------------------
# extract: HP-run identification and read extraction
# ---------------------------------------------------------------------------


def test_identify_hp_runs_includes_length_one():
    runs = extract.identify_hp_runs("ACGTAACCC")
    lengths = [(r["base"], r["length"]) for r in runs]
    # A, C, G, T, AA, CCC
    assert lengths == [("A", 1), ("C", 1), ("G", 1), ("T", 1), ("A", 2), ("C", 3)]


def test_extract_hp_length_clean():
    # consensus aligned: "ACAAAG" (A-run length 3 between C and G)
    # read aligned:     "ACAAAG"
    # left anchor (C) = col 1; right anchor (G) = col 5; region cols 2..4
    length, status = extract.extract_hp_length_from_read("ACAAAG", "A", 1, 5)
    assert (length, status) == (3, "clean")


def test_extract_hp_length_deletion():
    # read has one A missing: ACAAG
    length, status = extract.extract_hp_length_from_read("ACAA-G", "A", 1, 5)
    assert (length, status) == (2, "clean")


def test_extract_hp_length_complex():
    # read has a T inside the A run: ACATAG
    length, status = extract.extract_hp_length_from_read("ACATAG", "A", 1, 5)
    assert status == "complex"
    assert length == 2  # 2 A's, not counting the T


# ---------------------------------------------------------------------------
# extract: per-specimen pipeline
# ---------------------------------------------------------------------------


def test_process_specimen_round_trip(tmp_path):
    # Consensus has an L=3 A-run. 100 reads: 95 correct, 5 with a deletion.
    consensus = "GCAAATG"
    correct = ["GCAAATG"] * 95
    deleted = ["GCAA-TG"] * 5
    _write_msa(tmp_path, "spA", 100, consensus, correct + deleted)

    obs, non_hp = extract.process_specimen("spA", str(tmp_path / "cluster_debug"))
    assert obs is not None and non_hp is not None
    # Find the A-run observation
    a_runs = [o for o in obs if o.base == "A" and o.spoa_length == 3]
    assert len(a_runs) == 1
    a = a_runs[0]
    assert a.mode_length == 3  # mode of 95 threes + 5 twos = 3
    assert a.n_reads == 100
    n_correct = int(np.sum(a.called_lengths == 3))
    assert n_correct == 95


# ---------------------------------------------------------------------------
# analyze: aggregation + filtering
# ---------------------------------------------------------------------------


def test_aggregate_by_context_basic():
    # Two specimens, both A/L=3, each with frac_correct = 0.95
    obs = [
        extract.HPObservation(
            specimen=f"sp{i}", consensus_pos=2, base="A",
            spoa_length=3, mode_length=3,
            called_lengths=np.array([3] * 95 + [2] * 5, dtype=np.int16),
            n_reads=100, n_complex=0, n_excluded=0,
        )
        for i in range(2)
    ]
    results = analyze.aggregate_by_context(obs)
    ctx = results[("A", 3)]
    assert ctx.n_runs == 2
    assert ctx.n_obs == 200
    assert ctx.frac_correct == pytest.approx(0.95)
    assert ctx.frac_minus1 == pytest.approx(0.05)


def test_detect_outliers_flags_high_error():
    # Pool of 20 typical positions (5% error) + 1 outlier (50% error)
    typical = [
        extract.HPObservation(
            specimen=f"sp{i}", consensus_pos=i, base="A",
            spoa_length=3, mode_length=3,
            called_lengths=np.array([3] * 95 + [2] * 5, dtype=np.int16),
            n_reads=100, n_complex=0, n_excluded=0,
        )
        for i in range(20)
    ]
    outlier = extract.HPObservation(
        specimen="spOUT", consensus_pos=99, base="A",
        spoa_length=3, mode_length=3,
        called_lengths=np.array([3] * 50 + [2] * 50, dtype=np.int16),
        n_reads=100, n_complex=0, n_excluded=0,
    )
    obs = typical + [outlier]
    raw = analyze.aggregate_by_context(obs)
    flagged = analyze.detect_outliers(obs, raw, alpha=0.01)
    assert ("spOUT", 99) in flagged
    assert all(("sp", i) not in flagged for i in range(20))


def test_run_filter_pass_smoke():
    obs = [
        extract.HPObservation(
            specimen=f"sp{i}", consensus_pos=2, base="A",
            spoa_length=3, mode_length=3,
            called_lengths=np.array([3] * 95 + [2] * 5, dtype=np.int16),
            n_reads=100, n_complex=0, n_excluded=0,
        )
        for i in range(20)
    ]
    raw, filtered, n_excluded = analyze.run_filter_pass(obs)
    assert ("A", 3) in raw
    assert ("A", 3) in filtered
    assert n_excluded == 0  # All positions consistent


# ---------------------------------------------------------------------------
# selection: metadata-JSON index
# ---------------------------------------------------------------------------


def test_select_qualifying_threshold(tmp_path):
    _write_metadata(tmp_path, "pass1", primary_m=300, primary_err_factor=0.5)
    _write_metadata(tmp_path, "fail_ric", primary_m=50, primary_err_factor=0.5)
    _write_metadata(tmp_path, "fail_ef", primary_m=300, primary_err_factor=1.5)
    records = selection.index_specimens(str(tmp_path))
    qual = selection.select_qualifying(records, min_ric=200, max_err_factor=1.0)
    assert {r.name for r in qual} == {"pass1"}


def test_select_qualifying_low_ric_threshold(tmp_path):
    _write_metadata(tmp_path, "p1", primary_m=80, primary_err_factor=0.5)
    _write_metadata(tmp_path, "p2", primary_m=300, primary_err_factor=0.5)
    records = selection.index_specimens(str(tmp_path))
    qual = selection.select_qualifying(records, min_ric=50, max_err_factor=1.0)
    assert {r.name for r in qual} == {"p1", "p2"}


def test_select_qualifying_hard_fails_old_schema(tmp_path):
    _write_metadata(tmp_path, "old", primary_m=300, primary_err_factor=0.5,
                    schema_version="1.0")
    with pytest.raises(ValueError, match="schema_version"):
        selection.index_specimens(str(tmp_path))


def test_dominant_error_model(tmp_path):
    _write_metadata(tmp_path, "a", primary_m=300, primary_err_factor=0.5,
                    error_model="dorado-v5.0")
    _write_metadata(tmp_path, "b", primary_m=300, primary_err_factor=0.5,
                    error_model="dorado-v5.0")
    _write_metadata(tmp_path, "c", primary_m=300, primary_err_factor=0.5,
                    error_model="dorado-v3.5")
    records = selection.index_specimens(str(tmp_path))
    assert selection.dominant_error_model(records) == "dorado-v5.0"


# ---------------------------------------------------------------------------
# nonhp (approach 2) aggregation
# ---------------------------------------------------------------------------


def test_aggregate_nonhp_pools_clusters(tmp_path):
    # Two clusters, each with a single non-HP error column.
    # Consensus: "GCGAT" (all length-1 runs → all non-HP), 10 reads each
    # cluster 1: 9 correct + 1 with G→T substitution at col 0
    consensus1 = "GCGAT"
    reads1 = ["GCGAT"] * 9 + ["TCGAT"]
    _write_msa(tmp_path, "sp1", 10, consensus1, reads1)
    # cluster 2: 8 correct + 2 with a deletion at col 4
    consensus2 = "GCGAT"
    reads2 = ["GCGAT"] * 8 + ["GCGA-"] * 2
    msa2 = tmp_path / "cluster_debug" / "sp2-1.v1-RiC10-msa.fasta"
    msa2.write_text(_make_msa(consensus2, reads2))

    pooled, n_clusters, n_skipped, by_bin = nonhp.aggregate_nonhp(
        str(tmp_path), nonhp_min_ric=5, threads=1, progress=False,
    )
    assert n_clusters == 2
    assert pooled.substitutions == 1
    assert pooled.deletions == 2
    # Each cluster: 5 cols x 10 reads = 50 obs, minus the errors counted
    # separately. matches = 50 - 1 - 2 = wait, no: each cluster contributes
    # 5 cols x 10 reads = 50 observations. Subs (1) + del (2) + matches.
    # matches = 50+50 - 1 - 2 = 97
    assert pooled.matches == 97


def test_derive_nonhp_rates_shipped_scheme():
    pooled = extract.NonHPCounts(
        matches=994, substitutions=4, deletions=2, insertions=4,
        total_positions=100,
    )
    sub, indel = nonhp.derive_nonhp_rates(pooled)
    # denominator = matches + sub + del = 1000
    assert sub == pytest.approx(0.004)
    # indel = (del + ins) / denom = 6 / 1000
    assert indel == pytest.approx(0.006)


# ---------------------------------------------------------------------------
# output: YAML writer round-trip
# ---------------------------------------------------------------------------


def test_write_model_yaml_round_trip(tmp_path, monkeypatch):
    rates = {
        "non-hp-sub": 0.004,
        "non-hp-indel": 0.007,
        "hp-l1": 0.005, "hp-l2": 0.006, "hp-l3": 0.008,
        "hp-l4": 0.009, "hp-l5": 0.010,
    }
    fake_user_dir = tmp_path / "error_models"
    monkeypatch.setattr("speconsense.qctx._USER_DIR", fake_user_dir)
    monkeypatch.setattr("speconsense.fit_error_model.output.qctx.user_dir",
                        lambda: fake_user_dir)

    path, text = output.write_model(
        "test-model", rates,
        chemistry="R10.4.1",
        basecaller="Dorado SUP v5.0.0",
        source="unit test",
        dataset="synthetic",
    )
    assert path == fake_user_dir / "test-model.yaml"
    assert path.exists()

    # Load via qctx
    from speconsense import qctx
    loaded = qctx.load_table(str(path))
    assert loaded == rates


def test_write_model_refuses_overwrite(tmp_path, monkeypatch):
    rates = {k: 0.01 for k in ("non-hp-sub", "non-hp-indel",
                                "hp-l1", "hp-l2", "hp-l3", "hp-l4", "hp-l5")}
    fake_user_dir = tmp_path / "error_models"
    monkeypatch.setattr("speconsense.fit_error_model.output.qctx.user_dir",
                        lambda: fake_user_dir)
    output.write_model("x", rates, chemistry="c", basecaller="b",
                       source="s", dataset="d")
    with pytest.raises(FileExistsError):
        output.write_model("x", rates, chemistry="c", basecaller="b",
                           source="s", dataset="d")
    # With force, the second write succeeds
    output.write_model("x", rates, chemistry="c", basecaller="b",
                       source="s", dataset="d", force=True)


# ---------------------------------------------------------------------------
# compare: diff rendering
# ---------------------------------------------------------------------------


def test_render_diff_table_format():
    source = {
        "non-hp-sub": 0.004, "non-hp-indel": 0.007,
        "hp-l1": 0.005, "hp-l2": 0.006, "hp-l3": 0.008,
        "hp-l4": 0.009, "hp-l5": 0.010,
    }
    new = dict(source)
    text, warnings = compare.render_diff("dorado-v5.0", source, new)
    assert "dorado-v5.0" in text
    assert "non-hp-sub" in text
    assert warnings == []


def test_render_diff_flags_large_drift():
    source = {k: 0.01 for k in compare.CANONICAL_KEYS}
    new = {k: 0.001 for k in compare.CANONICAL_KEYS}  # 10x lower → ratio 0.1
    _, warnings = compare.render_diff("dorado-v5.0", source, new)
    assert set(warnings) == set(compare.CANONICAL_KEYS)


def test_build_rates_implies_per_position_transform():
    """hp-lN must be the implied per-position p satisfying (1-p)^N = frac_correct.

    Regression: an earlier draft computed hp-lN = 1 - frac_correct directly,
    which over-counted error at L>=2 by ~2x against the shipped table.
    """
    # Single (A, L=4) context with frac_correct = 0.9412.
    # Implied per-position: 1 - 0.9412^(1/4) ≈ 1 - 0.9850 ≈ 0.01505
    frac_correct = 0.9412
    n_obs = 10_000
    ctx = analyze.ContextResult(
        n_runs=100, n_obs=n_obs,
        all_called=np.array([], dtype=np.int16),  # not consumed by _build_rates
        frac_correct=frac_correct,
    )
    filtered_hp = {("A", 4): ctx}
    # Stub non-HP so the function returns its hp-l4 entry.
    nonhp_counts = extract.NonHPCounts(
        matches=995, substitutions=3, deletions=2, insertions=2, total_positions=100,
    )
    rates, missing = _build_rates(filtered_hp, nonhp_counts)
    assert "hp-l4" in rates
    expected_p = 1.0 - frac_correct ** (1.0 / 4)
    assert rates["hp-l4"] == pytest.approx(expected_p, abs=1e-6)
    # Other lengths should be missing because we provided no data for them.
    assert {"hp-l1", "hp-l2", "hp-l3", "hp-l5"} <= set(missing)


def test_resolve_source_model_precedence(monkeypatch):
    # Explicit override beats dominant
    monkeypatch.setattr(compare.qctx, "load_table",
                        lambda name: {"name": name})
    name, _ = compare.resolve_source_model("explicit", "dominant", "fallback")
    assert name == "explicit"
    name, _ = compare.resolve_source_model(None, "dominant", "fallback")
    assert name == "dominant"
    name, _ = compare.resolve_source_model(None, None, "fallback")
    assert name == "fallback"
