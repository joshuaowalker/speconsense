"""Smoke test for the rationalized core INFO log surface.

Asserts the new uniform stage→state format, that legacy/no-op strings have
been demoted out of INFO, and that ``Final`` is the last core INFO line
emitted (i.e. it lands after the outlier-filter step, reflecting actual
post-MAD totals).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Patterns that should NEVER appear at INFO level after the cleanup
LEGACY_INFO_FORBIDDEN = [
    "Position-based variant phasing enabled",
    "Pre-phasing merge",
    "Post-phasing merge",
    "After phasing, created",
    "Second phasing pass",
    "CER annotation",
    "Initial clustering produced",
]

# No-op messages that previously logged at INFO and should now be DEBUG-only
NOOP_INFO_FORBIDDEN = [
    "No clusters merged",
    "no reads moved",
    "no clusters split",
    "no candidates",
]

STATE_COLUMN_RE = re.compile(
    r"\[\s*\d+ clusters,\s+\d+/\d+ reads,\s+\d+\.\d%\]"
)


@pytest.fixture
def temp_dir():
    test_dir = tempfile.mkdtemp(prefix="speconsense_log_test_")
    original_dir = os.getcwd()
    os.chdir(test_dir)
    yield test_dir
    os.chdir(original_dir)
    shutil.rmtree(test_dir)


def _write_fastq(records, path):
    with open(path, "w") as f:
        SeqIO.write(records, f, "fastq")


def _make_records(seqs, prefix="r"):
    out = []
    for i, s in enumerate(seqs):
        out.append(SeqRecord(
            Seq(s), id=f"{prefix}{i}", description="",
            letter_annotations={"phred_quality": [30] * len(s)},
        ))
    return out


def _info_lines(stderr: str) -> list[str]:
    return [line for line in stderr.splitlines() if " - INFO - " in line]


def _run_speconsense(input_path: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "speconsense.core",
         input_path, "--algorithm", "greedy", "--min-size", "2",
         "--log-level", "INFO"],
        capture_output=True, text=True,
    )


def test_log_format_smoke(temp_dir):
    # Two well-separated cluster groups, enough reads each to survive
    # min-size and produce real consensus output.
    base_a = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    base_b = "TTTTTTTTAAAAAAAAGGGGGGGGCCCCCCCCTTTTTTTTAAAAAAAAGGGGGGGG"
    seqs = [base_a] * 5 + [base_b] * 5
    records = _make_records(seqs)
    _write_fastq(records, "input.fastq")

    result = _run_speconsense("input.fastq")
    assert result.returncode == 0, f"speconsense failed: {result.stderr}"

    stderr = result.stderr
    info_lines = _info_lines(stderr)
    assert info_lines, "Expected at least one INFO log line"

    for forbidden in LEGACY_INFO_FORBIDDEN:
        assert forbidden not in stderr, (
            f"Legacy INFO string '{forbidden}' still appears in log:\n{stderr}"
        )
    for forbidden in NOOP_INFO_FORBIDDEN:
        for line in info_lines:
            assert forbidden not in line, (
                f"No-op string '{forbidden}' should not be at INFO level:\n{line}"
            )

    # At least one stage line uses the new state column format
    assert any(STATE_COLUMN_RE.search(line) for line in info_lines), (
        f"No INFO line matched the new state column format:\n"
        + "\n".join(info_lines)
    )

    # 'Final' is the last core INFO line (the only line that may follow it
    # is the discards-write notice, which only fires when there are discards).
    final_indices = [i for i, line in enumerate(info_lines) if "Final" in line]
    assert final_indices, "Expected at least one 'Final' INFO line"
    final_idx = final_indices[-1]
    trailing = info_lines[final_idx + 1:]
    for line in trailing:
        assert "Wrote" in line and "discarded reads" in line, (
            f"Unexpected INFO line after 'Final':\n{line}"
        )

    # And the Final line itself uses the state column format
    final_line = info_lines[final_idx]
    assert STATE_COLUMN_RE.search(final_line), (
        f"'Final' line does not match state column format:\n{final_line}"
    )
