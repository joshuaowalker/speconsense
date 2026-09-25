"""quality_report.txt in incremental (--specimen + --aggregate-only) mode.

specimux-suite builds the summary dir one specimen at a time with
``--specimen`` and finishes with ``--aggregate-only``. The report must come
out the same as a single full-mode run, which requires the per-specimen
sidecars in --source/cluster_debug/ to carry the .filtered and overlap-merge
state across invocations.
"""

import glob
import logging
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from speconsense.quality_report import (
    REPORT_SIDECAR_SUFFIX,
    assemble_from_sidecars,
    find_source_dirs,
    load_report_sidecars,
    write_report_sidecar,
)
from speconsense.types import ConsensusInfo, OverlapMergeInfo


TEST_DATA_DIR = Path(__file__).parent / "data" / "overlap_test_clusters"


def _summarize(*args):
    result = subprocess.run(["speconsense-summarize", *args],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return result


def _report_body(path):
    """Report text minus the lines that legitimately differ between runs."""
    with open(path) as f:
        return [line for line in f
                if not line.startswith(("Generated:", "Summary:"))]


@pytest.fixture
def source(tmp_path):
    src = tmp_path / "clusters"
    shutil.copytree(TEST_DATA_DIR, src)
    return src


@pytest.fixture
def nested_source(tmp_path, source):
    """One core output directory per specimen, as specimux-suite lays it out:
    <root>/<id>/<id>-all.fasta — no top-level cluster_debug/."""
    root = tmp_path / "consensus"
    for fasta in source.glob("*-all.fasta"):
        sid = fasta.name[:-len("-all.fasta")]
        (root / sid / "cluster_debug").mkdir(parents=True)
        shutil.copy(fasta, root / sid / fasta.name)
    return root


def test_incremental_report_matches_full_mode_nested_layout(source, nested_source, tmp_path):
    full_dir = tmp_path / "full"
    inc_dir = tmp_path / "inc"
    _summarize("--source", str(source), "--summary-dir", str(full_dir))

    specimens = sorted(p.name for p in nested_source.iterdir())
    for sid in specimens:
        staging = tmp_path / "staging" / sid
        _summarize("--source", str(nested_source / sid),
                   "--summary-dir", str(staging), "--specimen", sid)
        shutil.copytree(staging, inc_dir, dirs_exist_ok=True)  # suite "publish"
    _summarize("--source", str(nested_source), "--summary-dir", str(inc_dir),
               "--aggregate-only")

    for sid in specimens:
        assert (nested_source / sid / "cluster_debug" /
                f"{sid}{REPORT_SIDECAR_SUFFIX}").exists()

    full = [l for l in _report_body(full_dir / "quality_report.txt")
            if not l.startswith("Source:")]
    inc = [l for l in _report_body(inc_dir / "quality_report.txt")
           if not l.startswith("Source:")]
    assert "OVERLAP MERGE ANALYSIS\n" in full
    assert inc == full


def test_stale_sidecar_is_ignored(source, tmp_path):
    sid = "test-prefix"
    out = tmp_path / "out"
    _summarize("--source", str(source), "--summary-dir", str(out), "--specimen", sid)

    # Core re-runs the specimen (rewriting -all.fasta) but summarize doesn't follow.
    fasta = source / f"{sid}-all.fasta"
    fasta.write_text(fasta.read_text() + "\n")

    result = _summarize("--source", str(source), "--summary-dir", str(out),
                        "--aggregate-only")
    assert "stale report sidecar" in result.stderr
    assert "OVERLAP MERGE ANALYSIS" not in (out / "quality_report.txt").read_text()


def test_incremental_report_matches_full_mode(source, tmp_path):
    full_dir = tmp_path / "full"
    inc_dir = tmp_path / "inc"
    _summarize("--source", str(source), "--summary-dir", str(full_dir))

    specimens = sorted(p.name[:-len("-all.fasta")]
                       for p in source.glob("*-all.fasta"))
    for sid in specimens:
        _summarize("--source", str(source), "--summary-dir", str(inc_dir),
                   "--specimen", sid)
    _summarize("--source", str(source), "--summary-dir", str(inc_dir),
               "--aggregate-only")

    full_report = _report_body(full_dir / "quality_report.txt")
    assert "OVERLAP MERGE ANALYSIS\n" in full_report
    assert _report_body(inc_dir / "quality_report.txt") == full_report

    # Sidecars live beside core's metadata, never in the summary deliverable.
    sidecars = sorted(os.path.basename(p) for p in
                      glob.glob(str(source / "cluster_debug" / f"*{REPORT_SIDECAR_SUFFIX}")))
    assert sidecars == sorted(f"{sid}{REPORT_SIDECAR_SUFFIX}" for sid in specimens)
    assert not list(inc_dir.rglob(f"*{REPORT_SIDECAR_SUFFIX}"))


def test_full_mode_writes_no_sidecars(source, tmp_path):
    _summarize("--source", str(source), "--summary-dir", str(tmp_path / "out"))
    assert not list(source.rglob(f"*{REPORT_SIDECAR_SUFFIX}"))


def test_aggregate_without_sidecars_still_writes_report(source, tmp_path):
    out = tmp_path / "out"
    _summarize("--source", str(source), "--summary-dir", str(out))
    os.unlink(out / "quality_report.txt")

    result = _summarize("--source", str(source), "--summary-dir", str(out),
                        "--aggregate-only")
    assert (out / "quality_report.txt").exists()
    assert "have no report sidecar" in result.stderr


def test_find_source_dirs_flat_and_nested(tmp_path):
    (tmp_path / "flat" / "cluster_debug").mkdir(parents=True)
    (tmp_path / "flat" / "a-all.fasta").write_text("")
    (tmp_path / "nested" / "b" / "cluster_debug").mkdir(parents=True)
    (tmp_path / "nested" / "c").mkdir(parents=True)
    (tmp_path / "nested" / "c" / "c-all.fasta").write_text("")
    (tmp_path / "nested" / "empty").mkdir()
    assert find_source_dirs(str(tmp_path / "flat")) == [str(tmp_path / "flat")]
    assert find_source_dirs(str(tmp_path / "nested")) == [
        str(tmp_path / "nested" / "b"), str(tmp_path / "nested" / "c")]


def _info(name, size=10):
    return ConsensusInfo(
        sample_name=name, cluster_id=name.rsplit("-", 1)[-1], sequence="ACGT",
        ric=size, size=size, file_path="x", snp_count=None, primers=None,
        raw_ric=None, rid=None, rid_min=None,
    )


PARAMS = {"min_ric": 3, "min_cer_factor": 1.0}


def _fake_source(tmp_path, sid="s1"):
    (tmp_path / f"{sid}-all.fasta").write_text(">x\nACGT\n")
    return str(tmp_path)


def test_sidecar_round_trip_joins_filtered_by_name(tmp_path, caplog):
    merge = OverlapMergeInfo(
        specimen="s1", iteration=2, input_clusters=["s1-1.v1", "s1-2.v1"],
        input_lengths=[600, 300], input_rics=[50, 20], overlap_bp=250,
        prefix_bp=0, suffix_bp=50, output_length=650,
    )
    src = _fake_source(tmp_path)
    write_report_sidecar(src, "s1", [_info("s1-1.v2")], [merge], PARAMS, "test")
    sidecars, stale = load_report_sidecars([src])
    assert list(sidecars) == ["s1"] and stale == []

    source_records = [_info("s1-1.v1"), _info("s1-1.v2", size=7)]
    with caplog.at_level(logging.WARNING):
        filtered, merges = assemble_from_sidecars(sidecars, source_records, PARAMS)
    assert [c.sample_name for c in filtered] == ["s1-1.v2"]
    assert filtered[0].size == 7  # metrics come from the fresh source load
    assert merges == [merge]
    assert not caplog.records


def test_sidecar_param_mismatch_and_missing_record_warn(tmp_path, caplog):
    src = _fake_source(tmp_path)
    write_report_sidecar(src, "s1", [_info("s1-1.v9")], [], PARAMS, "test")
    sidecars, _ = load_report_sidecars([src])
    with caplog.at_level(logging.WARNING):
        filtered, _ = assemble_from_sidecars(
            sidecars, [_info("s1-1.v1")], {**PARAMS, "min_cer_factor": 0.0})
    assert filtered == []
    messages = " ".join(r.getMessage() for r in caplog.records)
    assert "different filter parameters" in messages
    assert "not found in --source" in messages
