#!/usr/bin/env python3
"""
Integration tests for overlap merge feature using real synthetic test data.

These tests use pre-generated FASTA files representing various overlap scenarios:
- test-core-overlap: Tests specimen name with '-c' in it (bug fix verification)
- test-prefix: Containment case (shorter sequence is prefix of longer)
- test-prefix-suffix-full: 3-way merge case (prefix + suffix + full)
"""

import os
import tempfile
import subprocess
import pytest
from pathlib import Path


# Path to test data
TEST_DATA_DIR = Path(__file__).parent / "data" / "overlap_test_clusters"


def run_summarize(source_dir: str, output_dir: str, extra_args: list = None) -> subprocess.CompletedProcess:
    """Run speconsense-summarize and return the result."""
    cmd = [
        "speconsense-summarize",
        "--source", source_dir,
        "--summary-dir", output_dir,
    ]
    if extra_args:
        cmd.extend(extra_args)

    return subprocess.run(cmd, capture_output=True, text=True)


def count_fasta_records(fasta_path: str) -> int:
    """Count the number of records in a FASTA file."""
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def parse_fasta_headers(fasta_path: str) -> list:
    """Parse FASTA headers and return list of (name, attributes) tuples."""
    headers = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split()
                name = parts[0]
                attrs = {}
                for part in parts[1:]:
                    if "=" in part:
                        key, value = part.split("=", 1)
                        attrs[key] = value
                headers.append((name, attrs))
    return headers


class TestOverlapMergeIntegration:
    """Integration tests for overlap merge feature."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_specimen_name_with_dash_c_not_truncated(self, temp_output_dir):
        """
        Test that specimen names containing '-c' are not incorrectly truncated.

        Bug: 'test-core-overlap-c1' was being parsed as 'test-' instead of 'test-core-overlap'
        because split('-c')[0] was used instead of rsplit('-c', 1)[0].
        """
        result = run_summarize(str(TEST_DATA_DIR), temp_output_dir)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        summary_fasta = os.path.join(temp_output_dir, "summary.fasta")
        assert os.path.exists(summary_fasta), "summary.fasta not created"

        headers = parse_fasta_headers(summary_fasta)

        # Find the test-core-overlap specimen
        core_overlap_headers = [h for h in headers if h[0].startswith("test-core-overlap")]
        assert len(core_overlap_headers) >= 1, "No test-core-overlap variants found"

        # Verify the full specimen name is preserved (not truncated to 'test-')
        for name, _ in core_overlap_headers:
            assert name.startswith("test-core-overlap-"), f"Specimen name truncated: {name}"
            # Should NOT start with just 'test-' followed by digit
            assert not name.startswith("test-1"), f"Specimen name incorrectly truncated to 'test': {name}"

    def test_containment_merge_with_overlap_enabled(self, temp_output_dir):
        """
        Test that containment case (prefix inside full) merges correctly.

        test-prefix has:
        - c1: 248 reads, ~361bp prefix sequence
        - c2: 223 reads, ~630bp full sequence

        Should merge to 1 variant with RiC ~471 (248+223).
        """
        result = run_summarize(str(TEST_DATA_DIR), temp_output_dir)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        summary_fasta = os.path.join(temp_output_dir, "summary.fasta")
        headers = parse_fasta_headers(summary_fasta)

        # Find test-prefix variants
        prefix_headers = [h for h in headers if h[0].startswith("test-prefix-") and not h[0].startswith("test-prefix-suffix")]

        # Should have exactly 1 merged variant
        assert len(prefix_headers) == 1, f"Expected 1 merged variant, got {len(prefix_headers)}: {prefix_headers}"

        name, attrs = prefix_headers[0]

        # Verify merge occurred
        assert "rawric" in attrs, "No rawric field - merge may not have occurred"
        assert "+" in attrs["rawric"], f"rawric should show merged clusters: {attrs['rawric']}"

        # Total RiC should be sum of original clusters
        ric = int(attrs.get("ric", 0))
        assert ric == 471, f"Expected RiC=471 (248+223), got {ric}"

    def test_three_way_iterative_merge(self, temp_output_dir):
        """
        Test 3-way merge case: prefix + suffix + full sequences.

        test-prefix-suffix-full has:
        - c1: 167 reads, prefix sequence (~360bp)
        - c2: 166 reads, suffix sequence (~269bp)
        - c3: 163 reads, full sequence (~629bp)

        prefix and suffix don't overlap with each other, only with full.
        Should use iterative merging:
        1. First iteration: merge prefix+full
        2. Second iteration: merge result+suffix

        Final result: 1 variant with RiC ~496 (167+166+163).
        """
        result = run_summarize(str(TEST_DATA_DIR), temp_output_dir)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        summary_fasta = os.path.join(temp_output_dir, "summary.fasta")
        headers = parse_fasta_headers(summary_fasta)

        # Find test-prefix-suffix-full variants
        psf_headers = [h for h in headers if h[0].startswith("test-prefix-suffix-full")]

        # Should have exactly 1 merged variant (all 3 merged together)
        assert len(psf_headers) == 1, f"Expected 1 merged variant, got {len(psf_headers)}: {psf_headers}"

        name, attrs = psf_headers[0]

        # Verify merge occurred
        assert "rawric" in attrs, "No rawric field - merge may not have occurred"
        assert "+" in attrs["rawric"], f"rawric should show merged clusters: {attrs['rawric']}"

        # Total RiC should be sum of all three original clusters
        ric = int(attrs.get("ric", 0))
        assert ric == 496, f"Expected RiC=496 (167+166+163), got {ric}"

        # Verify iterative merge happened (rawric should show intermediate merge)
        # After iteration 1: c1(167)+c3(163)=330, after iteration 2: 330+c2(166)=496
        rawric = attrs["rawric"]
        # rawric could be "330+166" or similar showing the iterative merge
        parts = rawric.split("+")
        assert len(parts) == 2, f"Expected 2 parts in rawric showing iterative merge: {rawric}"

    def test_overlap_disabled_no_merge(self, temp_output_dir):
        """
        Test that with --min-merge-overlap 0, different-length sequences don't merge.

        This verifies the original behavior is preserved when overlap merging is disabled.
        """
        result = run_summarize(str(TEST_DATA_DIR), temp_output_dir, ["--min-merge-overlap", "0"])
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        summary_fasta = os.path.join(temp_output_dir, "summary.fasta")
        headers = parse_fasta_headers(summary_fasta)

        # test-prefix should have 2 variants (not merged)
        prefix_headers = [h for h in headers if h[0].startswith("test-prefix-") and not h[0].startswith("test-prefix-suffix")]
        assert len(prefix_headers) == 2, f"Expected 2 unmerged variants with overlap disabled, got {len(prefix_headers)}"

        # test-prefix-suffix-full should have 3 variants (not merged)
        psf_headers = [h for h in headers if h[0].startswith("test-prefix-suffix-full")]
        assert len(psf_headers) == 3, f"Expected 3 unmerged variants with overlap disabled, got {len(psf_headers)}"

    def test_core_overlap_merges_correctly(self, temp_output_dir):
        """
        Test that test-core-overlap (two overlapping sequences) merges correctly.

        test-core-overlap has:
        - c1: 250 reads, ~500bp sequence
        - c2: 250 reads, ~500bp sequence (overlapping with c1)

        Should merge to 1 variant with RiC=500 (250+250).
        """
        result = run_summarize(str(TEST_DATA_DIR), temp_output_dir)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        summary_fasta = os.path.join(temp_output_dir, "summary.fasta")
        headers = parse_fasta_headers(summary_fasta)

        # Find test-core-overlap variants
        core_headers = [h for h in headers if h[0].startswith("test-core-overlap")]

        # Should have exactly 1 merged variant
        assert len(core_headers) == 1, f"Expected 1 merged variant, got {len(core_headers)}"

        name, attrs = core_headers[0]
        ric = int(attrs.get("ric", 0))
        assert ric == 500, f"Expected RiC=500 (250+250), got {ric}"

    def test_output_files_created(self, temp_output_dir):
        """
        Test that all expected output files are created.
        """
        result = run_summarize(str(TEST_DATA_DIR), temp_output_dir)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Check standard output files exist
        expected_files = [
            "summary.fasta",
            "summary.txt",
            "quality_report.txt",
            "summarize_log.txt",
        ]

        for filename in expected_files:
            filepath = os.path.join(temp_output_dir, filename)
            assert os.path.exists(filepath), f"Expected file not created: {filename}"

        # Check per-specimen files exist
        specimen_files = [
            "test-core-overlap-1.v1-RiC500.fasta",
            "test-prefix-1.v1-RiC471.fasta",
            "test-prefix-suffix-full-1.v1-RiC496.fasta",
        ]

        for filename in specimen_files:
            filepath = os.path.join(temp_output_dir, filename)
            assert os.path.exists(filepath), f"Expected specimen file not created: {filename}"
