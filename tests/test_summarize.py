#!/usr/bin/env python3
"""
Tests for speconsense-summarize functionality.

Tests focus on merge behavior with different sequence variants.
"""

import tempfile
import os
import shutil
import subprocess
import re
from Bio import SeqIO


def test_merge_behavior_with_full_hac_context():
    """Test merge behavior using complete real specimen file with multiple clusters.

    This uses the ONT01.06-F01--iNat233404862-all.fasta file which contains
    9 clusters. This test demonstrates that merge decisions depend on the full HAC
    group context, not just pairwise comparisons.

    The file contains:
    - c1: main cluster (ric=500) - majority pattern
    - c2: second major cluster (ric=250) - majority pattern
    - c3: contamination (ric=9) - separate group
    - c4: variant (ric=6) - ends with TAG, structural variant
    - c5: variant (ric=6) - structural variant
    - c6: contamination (ric=4) - separate group
    - c7: variant (ric=3) - ends with TAA
    - c8: variant (ric=3) - majority pattern
    - c9: variant (ric=3) - ends with TAA, homopolymer variation from c7

    Expected behavior:
    - c1 + c2 merge into 1.v1 with rawric=500+250
    - c7 + c9 merge into 1.v4 with rawric=3+3 (both end with TAA, differ only by homopolymers)
    - c4 stays separate as 1.v3 (ends with TAG, not compatible with TAA sequences)
    - Total: 7 sequences across 3 HAC groups
    """
    # Use the test data file
    test_file = os.path.join(os.path.dirname(__file__), "data", "ONT01.06-F01--iNat233404862-all.fasta")

    # Skip test if file doesn't exist (e.g., on CI)
    if not os.path.exists(test_file):
        import pytest
        pytest.skip(f"Test file not found: {test_file}")

    # Create temporary directory for output
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # Copy the test file to our temp directory
        import shutil as shutil_module
        dest_file = os.path.join(source_dir, "ONT01.06-F01--iNat233404862-all.fasta")
        shutil_module.copy(test_file, dest_file)

        # Run speconsense-summarize with default parameters
        # Disable overlap merge (--min-merge-overlap 0) to test original behavior
        result = subprocess.run(
            [
                "speconsense-summarize",
                "--source", source_dir,
                "--summary-dir", summary_dir,
                "--min-ric", "3",  # Include c4, c7, c8, c9 (all have ric >= 3)
                "--min-merge-overlap", "0"  # Disable overlap merge for this test
            ],
            capture_output=True,
            text=True
        )

        # Check that the command succeeded
        assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

        # Read the main output FASTA file
        output_fasta = os.path.join(summary_dir, "summary.fasta")
        assert os.path.exists(output_fasta), \
            f"Expected output file not found: {output_fasta}"

        # Read all sequences from the output
        output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

        # Print diagnostic information
        print(f"\nOutput sequences: {len(output_sequences)}")
        for seq in output_sequences:
            print(f"  {seq.id}: {seq.description}")

        # Verify the expected number of output sequences (7 total for this specimen)
        assert len(output_sequences) == 7, \
            f"Expected 7 output sequences, got {len(output_sequences)}"

        # Check if c7 and c9 were merged by examining the output sequences
        # Look for sequences with rawric field indicating a merge
        c7_c9_merged = False
        merged_into = None

        for seq in output_sequences:
            # Check if this sequence has rawric=3+3 (c7 and c9)
            if 'rawric=' in seq.description:
                # Extract rawric values
                rawric_match = re.search(r'rawric=([\d+]+)', seq.description)
                if rawric_match:
                    rawric_str = rawric_match.group(1)
                    ric_values = [int(x) for x in rawric_str.split('+')]
                    # Check if both values are 3 (c7 and c9)
                    if ric_values == [3, 3]:
                        c7_c9_merged = True
                        merged_into = seq.id
                        break

        # Alternative check: Look at .raw files in variants directory
        # If c7 and c9 are in the same variant group, they were merged
        variants_dir = os.path.join(summary_dir, "variants")
        if os.path.exists(variants_dir) and not c7_c9_merged:
            specimen_raw_files = sorted([f for f in os.listdir(variants_dir)
                                         if f.startswith('ONT01.06-F01--iNat233404862') and '.raw' in f])

            # Group raw files by their variant (e.g., "1.v4")
            variant_groups = {}
            for raw_file in specimen_raw_files:
                # Extract variant identifier (e.g., "1.v4" from "...1.v4.raw1...")
                match = re.search(r'-(\d+\.v\d+)\.raw', raw_file)
                if match:
                    variant_id = match.group(1)
                    if variant_id not in variant_groups:
                        variant_groups[variant_id] = []
                    variant_groups[variant_id].append(raw_file)

            # Check each variant group for both c7 and c9 (both end with TAA)
            for variant_id, raw_files in variant_groups.items():
                has_c7_or_c9 = 0

                for raw_file in raw_files:
                    raw_path = os.path.join(variants_dir, raw_file)
                    raw_seqs = list(SeqIO.parse(raw_path, "fasta"))

                    for seq in raw_seqs:
                        seq_str = str(seq.seq)
                        # Both c7 and c9 end with TAA
                        if seq_str.endswith('GACCTCAAATCAGGTAGGACTACCCGCTGAACTTAA'):
                            has_c7_or_c9 += 1

                if has_c7_or_c9 >= 2:  # Found at least 2 sequences ending with TAA
                    c7_c9_merged = True
                    merged_into = variant_id
                    break

        # Key assertion: c7 and c9 SHOULD be merged when in HAC group context
        # This is because the multi-sequence alignment reveals their differences
        # are homopolymer variations (both end with TAA)
        assert c7_c9_merged, \
            f"c7 (ric=3) and c9 (ric=3) should be merged in HAC group context, " \
            f"but they were not merged"

        # Verify they merged into the expected variant group (1.v4)
        assert merged_into == '1.v4' or merged_into == 'ONT01.06-F01--iNat233404862-1.v4', \
            f"Expected c7 and c9 to merge into variant 1.v4, but merged into: {merged_into}"

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)


def test_merge_with_homopolymer_only_differences():
    """Test that sequences differing only in homopolymer lengths DO merge.

    This test verifies that sequences with identical structure but different
    homopolymer lengths will merge with the homopolymer-aware algorithm.
    """
    # Create temporary directory structure
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # Create two sequences that differ only in homopolymer length
        # Base sequence with A homopolymer of length 5
        seq1 = "ATCGAAAAATCGATCGATCGATCG"
        # Same sequence with A homopolymer of length 8
        seq2 = "ATCGAAAAAAATCGATCGATCGATCG"

        fasta_content = f""">test-seq1 size=10 ric=10 primers=test
{seq1}
>test-seq2 size=8 ric=8 primers=test
{seq2}
"""

        fasta_file = os.path.join(source_dir, "test-homopoly-all.fasta")
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)

        # Run speconsense-summarize with default parameters
        result = subprocess.run(
            [
                "speconsense-summarize",
                "--source", source_dir,
                "--summary-dir", summary_dir,
                "--min-ric", "3"
            ],
            capture_output=True,
            text=True
        )

        # Check that the command succeeded
        assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

        # Read the main output FASTA file (summary.fasta combines all specimens)
        output_fasta = os.path.join(summary_dir, "summary.fasta")
        assert os.path.exists(output_fasta), \
            f"Expected output file not found: {output_fasta}"

        # Count sequences in output
        output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

        # Should have 1 sequence (merged due to homopolymer equivalence)
        assert len(output_sequences) == 1, \
            f"Expected 1 merged sequence, but got {len(output_sequences)}"

        # The merged sequence should have the combined size
        # Check the header for size information
        header = output_sequences[0].description
        assert "size=" in header, "Expected size field in output header"

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)


def test_merge_bases_to_iupac_expands_existing_codes():
    """Test that merge_bases_to_iupac correctly expands existing IUPAC codes.

    This tests the fix for a bug where merging a base with an existing IUPAC
    code would produce 'N' instead of the correct expanded code.
    For example, C + Y should produce Y (since Y = C|T, and C is already in Y).
    """
    from speconsense.summarize import merge_bases_to_iupac

    # Test cases: (input_bases, expected_output)
    test_cases = [
        # Bug fix cases: existing IUPAC codes should be expanded
        ({'C', 'Y'}, 'Y'),   # C + Y(CT) = CT = Y
        ({'T', 'Y'}, 'Y'),   # T + Y(CT) = CT = Y
        ({'A', 'R'}, 'R'),   # A + R(AG) = AG = R
        ({'G', 'R'}, 'R'),   # G + R(AG) = AG = R
        ({'C', 'R'}, 'V'),   # C + R(AG) = ACG = V
        ({'T', 'R'}, 'D'),   # T + R(AG) = AGT = D
        ({'Y', 'R'}, 'N'),   # Y(CT) + R(AG) = ACGT = N

        # Standard cases: no existing IUPAC codes
        ({'A'}, 'A'),        # Single base stays the same
        ({'C', 'T'}, 'Y'),   # C + T = CT = Y
        ({'A', 'G'}, 'R'),   # A + G = AG = R
        ({'A', 'C', 'G', 'T'}, 'N'),  # All four = N

        # More complex IUPAC expansion cases
        ({'M', 'K'}, 'N'),   # M(AC) + K(GT) = ACGT = N
        ({'S', 'W'}, 'N'),   # S(GC) + W(AT) = ACGT = N
        ({'B', 'A'}, 'N'),   # B(CGT) + A = ACGT = N
        ({'V', 'T'}, 'N'),   # V(ACG) + T = ACGT = N
    ]

    for bases, expected in test_cases:
        result = merge_bases_to_iupac(bases)
        assert result == expected, \
            f"merge_bases_to_iupac({bases}) returned '{result}', expected '{expected}'"


class TestPrimersAreSame:
    """Tests for primers_are_same() function used in overlap merge constraint."""

    def test_same_primers_exact_match(self):
        """Same primers should return True (use global distance)."""
        from speconsense.summarize import primers_are_same
        assert primers_are_same(['ITS1', 'ITS4'], ['ITS1', 'ITS4']) is True
        assert primers_are_same(['fwd', 'rev'], ['fwd', 'rev']) is True

    def test_same_primers_different_order(self):
        """Same primers in different order should return True."""
        from speconsense.summarize import primers_are_same
        assert primers_are_same(['ITS4', 'ITS1'], ['ITS1', 'ITS4']) is True
        assert primers_are_same(['rev', 'fwd'], ['fwd', 'rev']) is True

    def test_different_primers(self):
        """Different primers should return False (allow overlap merge)."""
        from speconsense.summarize import primers_are_same
        assert primers_are_same(['ITS1', 'ITS4'], ['ITS1', 'ITS2']) is False
        assert primers_are_same(['fwd_a', 'rev_a'], ['fwd_b', 'rev_b']) is False

    def test_none_primers_conservative(self):
        """None primers should return True (conservative: unknown = same)."""
        from speconsense.summarize import primers_are_same
        assert primers_are_same(None, None) is True
        assert primers_are_same(None, ['ITS1', 'ITS4']) is True
        assert primers_are_same(['ITS1', 'ITS4'], None) is True

    def test_empty_list_conservative(self):
        """Empty list should return True (conservative: unknown = same)."""
        from speconsense.summarize import primers_are_same
        assert primers_are_same([], []) is True
        assert primers_are_same([], ['ITS1', 'ITS4']) is True
        assert primers_are_same(['ITS1', 'ITS4'], []) is True

    def test_single_primer_overlap(self):
        """Partial primer overlap should be treated as different."""
        from speconsense.summarize import primers_are_same
        # Different sets = different amplicons
        assert primers_are_same(['ITS1'], ['ITS1', 'ITS4']) is False
        assert primers_are_same(['ITS1', 'ITS4'], ['ITS4']) is False

    def test_single_primer_same(self):
        """Single primer that matches should return True."""
        from speconsense.summarize import primers_are_same
        assert primers_are_same(['ITS1'], ['ITS1']) is True


class TestMergeEffort:
    """Tests for --merge-effort parameter parsing and batch size computation."""

    def test_parse_presets(self):
        """Test preset name parsing."""
        from speconsense.summarize.cli import parse_merge_effort
        assert parse_merge_effort("fast") == 8
        assert parse_merge_effort("balanced") == 10
        assert parse_merge_effort("thorough") == 12

    def test_parse_presets_case_insensitive(self):
        """Test that presets are case-insensitive."""
        from speconsense.summarize.cli import parse_merge_effort
        assert parse_merge_effort("BALANCED") == 10
        assert parse_merge_effort("Fast") == 8
        assert parse_merge_effort("THOROUGH") == 12

    def test_parse_presets_whitespace(self):
        """Test that whitespace is stripped."""
        from speconsense.summarize.cli import parse_merge_effort
        assert parse_merge_effort("  balanced  ") == 10
        assert parse_merge_effort("\tfast\n") == 8

    def test_parse_numeric(self):
        """Test numeric value parsing."""
        from speconsense.summarize.cli import parse_merge_effort
        assert parse_merge_effort("6") == 6
        assert parse_merge_effort("10") == 10
        assert parse_merge_effort("14") == 14

    def test_parse_numeric_at_bounds(self):
        """Test numeric values at the valid boundaries."""
        from speconsense.summarize.cli import parse_merge_effort
        assert parse_merge_effort("6") == 6   # Minimum
        assert parse_merge_effort("14") == 14  # Maximum

    def test_parse_invalid_preset(self):
        """Test that invalid preset names raise ValueError."""
        import pytest
        from speconsense.summarize.cli import parse_merge_effort
        with pytest.raises(ValueError, match="Unknown merge-effort"):
            parse_merge_effort("invalid")
        with pytest.raises(ValueError, match="Unknown merge-effort"):
            parse_merge_effort("medium")

    def test_parse_numeric_below_minimum(self):
        """Test that values below minimum raise ValueError."""
        import pytest
        from speconsense.summarize.cli import parse_merge_effort
        with pytest.raises(ValueError, match="must be 6-14"):
            parse_merge_effort("5")
        with pytest.raises(ValueError, match="must be 6-14"):
            parse_merge_effort("0")

    def test_parse_numeric_above_maximum(self):
        """Test that values above maximum raise ValueError."""
        import pytest
        from speconsense.summarize.cli import parse_merge_effort
        with pytest.raises(ValueError, match="must be 6-14"):
            parse_merge_effort("15")
        with pytest.raises(ValueError, match="must be 6-14"):
            parse_merge_effort("20")

    def test_batch_size_balanced_small_groups(self):
        """Test batch size computation for balanced effort with small groups."""
        from speconsense.summarize.analysis import compute_merge_batch_size
        # E=10 (balanced): groups <= 8 should get batch=8
        assert compute_merge_batch_size(4, 10) == 8
        assert compute_merge_batch_size(8, 10) == 8

    def test_batch_size_balanced_medium_groups(self):
        """Test batch size computation for balanced effort with medium groups."""
        from speconsense.summarize.analysis import compute_merge_batch_size
        # E=10: batch decreases as group size increases
        assert compute_merge_batch_size(16, 10) == 7
        assert compute_merge_batch_size(32, 10) == 6
        assert compute_merge_batch_size(64, 10) == 5

    def test_batch_size_balanced_large_groups(self):
        """Test batch size computation for balanced effort with large groups."""
        from speconsense.summarize.analysis import compute_merge_batch_size
        # E=10: large groups hit MIN_BATCH=4
        assert compute_merge_batch_size(128, 10) == 4
        assert compute_merge_batch_size(256, 10) == 4
        assert compute_merge_batch_size(512, 10) == 4

    def test_batch_size_fast(self):
        """Test batch size computation for fast effort (E=8)."""
        from speconsense.summarize.analysis import compute_merge_batch_size
        assert compute_merge_batch_size(8, 8) == 6
        assert compute_merge_batch_size(16, 8) == 5
        assert compute_merge_batch_size(32, 8) == 4

    def test_batch_size_thorough(self):
        """Test batch size computation for thorough effort (E=12)."""
        from speconsense.summarize.analysis import compute_merge_batch_size
        assert compute_merge_batch_size(32, 12) == 8
        assert compute_merge_batch_size(64, 12) == 7
        assert compute_merge_batch_size(128, 12) == 6

    def test_batch_size_edge_cases(self):
        """Test batch size computation edge cases."""
        from speconsense.summarize.analysis import compute_merge_batch_size
        # Single variant returns 1
        assert compute_merge_batch_size(1, 10) == 1
        # Two variants with high effort -> clamped to MAX_BATCH=8
        assert compute_merge_batch_size(2, 10) == 8

    def test_batch_size_clamped_to_max(self):
        """Test that batch size is clamped to MAX_MERGE_BATCH=8."""
        from speconsense.summarize.analysis import compute_merge_batch_size, MAX_MERGE_BATCH
        # Very small group with high effort should still be clamped to 8
        assert compute_merge_batch_size(2, 14) == MAX_MERGE_BATCH
        assert compute_merge_batch_size(4, 14) == MAX_MERGE_BATCH

    def test_batch_size_clamped_to_min(self):
        """Test that batch size is clamped to MIN_MERGE_BATCH=4."""
        from speconsense.summarize.analysis import compute_merge_batch_size, MIN_MERGE_BATCH
        # Very large group should be clamped to 4
        assert compute_merge_batch_size(1000, 10) == MIN_MERGE_BATCH
        assert compute_merge_batch_size(10000, 6) == MIN_MERGE_BATCH
