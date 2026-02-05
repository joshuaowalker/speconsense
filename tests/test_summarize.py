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


class TestFullConsensus:
    """Tests for create_full_consensus_from_msa and --enable-full-consensus feature."""

    def test_full_consensus_indel_always_included(self):
        """Test that gaps never win - indel content is always included."""
        from speconsense.summarize.merging import create_full_consensus_from_msa
        from speconsense.types import ConsensusInfo
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # Sequence 1: ACGT (no insertion)
        # Sequence 2: AC-GT with an extra base -> aligned as AC[A]GT
        # In MSA: seq1 = "AC-GT", seq2 = "ACAGT"
        aligned = [
            SeqRecord(Seq("AC-GT"), id="s1"),
            SeqRecord(Seq("ACAGT"), id="s2"),
        ]
        variants = [
            ConsensusInfo("s1", "c1", "ACGT", ric=10, size=10, file_path="f"),
            ConsensusInfo("s2", "c2", "ACAGT", ric=5, size=5, file_path="f"),
        ]

        result = create_full_consensus_from_msa(aligned, variants)

        # Full consensus should include the inserted A (gaps never win)
        assert result.sequence == "ACAGT"
        assert result.ric == 15
        assert result.size == 15

    def test_full_consensus_snp_produces_iupac(self):
        """Test that SNP positions produce IUPAC ambiguity codes."""
        from speconsense.summarize.merging import create_full_consensus_from_msa
        from speconsense.types import ConsensusInfo
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # Two sequences differing at position 2: A vs G
        aligned = [
            SeqRecord(Seq("AACGT"), id="s1"),
            SeqRecord(Seq("AGCGT"), id="s2"),
        ]
        variants = [
            ConsensusInfo("s1", "c1", "AACGT", ric=10, size=10, file_path="f"),
            ConsensusInfo("s2", "c2", "AGCGT", ric=8, size=8, file_path="f"),
        ]

        result = create_full_consensus_from_msa(aligned, variants)

        # Position 2 should be R (A or G)
        assert result.sequence == "ARCGT"
        assert result.snp_count == 1

    def test_full_consensus_metadata_aggregation(self):
        """Test that metadata is correctly aggregated from all variants."""
        from speconsense.summarize.merging import create_full_consensus_from_msa
        from speconsense.types import ConsensusInfo
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        aligned = [
            SeqRecord(Seq("ACGT"), id="s1"),
            SeqRecord(Seq("ACGT"), id="s2"),
        ]
        variants = [
            ConsensusInfo("s1", "c1", "ACGT", ric=20, size=20, file_path="f",
                         primers=["ITS1", "ITS4"], rid=0.95, rid_min=0.90),
            ConsensusInfo("s2", "c2", "ACGT", ric=10, size=10, file_path="f",
                         primers=["ITS1", "ITS4"], rid=0.90, rid_min=0.85),
        ]

        result = create_full_consensus_from_msa(aligned, variants)

        assert result.ric == 30
        assert result.size == 30
        assert result.raw_ric == [20, 10]  # Sorted descending
        assert result.raw_len == [4, 4]
        # Largest variant (s1) metadata used for primers/rid
        assert result.primers == ["ITS1", "ITS4"]
        assert result.rid == 0.95


    def test_full_consensus_filters_small_variants(self):
        """Integration test: merge_min_size_ratio filters small variants from full consensus."""
        temp_dir = tempfile.mkdtemp()
        source_dir = os.path.join(temp_dir, "clusters")
        summary_dir = os.path.join(temp_dir, "__Summary__")
        os.makedirs(source_dir)

        try:
            # Two similar sequences (1 SNP at position 12: G vs A)
            # Very different sizes so the small one is filtered by merge_min_size_ratio
            seq_large = "ATCGATCGATCGATCGATCGATCG"  # G at position 12
            seq_small = "ATCGATCGATCAATCGATCGATCG"  # A at position 12

            fasta_content = f""">test-c1 size=100 ric=100 primers=test
{seq_large}
>test-c2 size=5 ric=5 primers=test
{seq_small}
"""
            fasta_file = os.path.join(source_dir, "test-all.fasta")
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)

            # merge-min-size-ratio 0.1 filters 5/100=0.05 from full consensus
            result = subprocess.run(
                [
                    "speconsense-summarize",
                    "--source", source_dir,
                    "--summary-dir", summary_dir,
                    "--min-ric", "3",
                    "--enable-full-consensus",
                    "--merge-min-size-ratio", "0.1",
                    "--disable-merging",
                    "--min-merge-overlap", "0",
                ],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

            output_fasta = os.path.join(summary_dir, "summary.fasta")
            output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

            full_seqs = [s for s in output_sequences if '.full' in s.id]
            assert len(full_seqs) == 1, f"Expected 1 .full sequence, got {len(full_seqs)}"

            # Small variant was filtered — .full should be the large variant only (no IUPAC)
            full_seq_str = str(full_seqs[0].seq)
            assert full_seq_str == seq_large, \
                f"Expected large variant sequence, got {full_seq_str}"

        finally:
            shutil.rmtree(temp_dir)

    def test_full_consensus_no_filter_when_disabled(self):
        """Integration test: merge_min_size_ratio=0 preserves all variants in full consensus."""
        temp_dir = tempfile.mkdtemp()
        source_dir = os.path.join(temp_dir, "clusters")
        summary_dir = os.path.join(temp_dir, "__Summary__")
        os.makedirs(source_dir)

        try:
            # Same sequences as above — 1 SNP at position 12 (G vs A)
            seq_large = "ATCGATCGATCGATCGATCGATCG"
            seq_small = "ATCGATCGATCAATCGATCGATCG"

            fasta_content = f""">test-c1 size=100 ric=100 primers=test
{seq_large}
>test-c2 size=5 ric=5 primers=test
{seq_small}
"""
            fasta_file = os.path.join(source_dir, "test-all.fasta")
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)

            # merge-min-size-ratio 0 disables filtering — both contribute to .full
            result = subprocess.run(
                [
                    "speconsense-summarize",
                    "--source", source_dir,
                    "--summary-dir", summary_dir,
                    "--min-ric", "3",
                    "--enable-full-consensus",
                    "--merge-min-size-ratio", "0",
                    "--disable-merging",
                    "--min-merge-overlap", "0",
                ],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

            output_fasta = os.path.join(summary_dir, "summary.fasta")
            output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

            full_seqs = [s for s in output_sequences if '.full' in s.id]
            assert len(full_seqs) == 1, f"Expected 1 .full sequence, got {len(full_seqs)}"

            # Both variants contribute — SNP position should be IUPAC R (A/G)
            full_seq_str = str(full_seqs[0].seq)
            assert "R" in full_seq_str, \
                f"Expected IUPAC R (A/G) in full consensus, got {full_seq_str}"

        finally:
            shutil.rmtree(temp_dir)


class TestFieldRegexFullConsensus:
    """Tests for GroupField and VariantField regex handling of .full names."""

    def test_group_field_matches_full(self):
        """GroupField should extract group number from .full names."""
        from speconsense.summarize.fields import GroupField
        from speconsense.types import ConsensusInfo

        field = GroupField()
        cons = ConsensusInfo("specimen-1.full", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "group=1"

    def test_group_field_matches_full_multidigit(self):
        """GroupField should handle multi-digit group numbers with .full."""
        from speconsense.summarize.fields import GroupField
        from speconsense.types import ConsensusInfo

        field = GroupField()
        cons = ConsensusInfo("specimen-12.full", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "group=12"

    def test_group_field_still_matches_variant(self):
        """GroupField should still match standard .v names."""
        from speconsense.summarize.fields import GroupField
        from speconsense.types import ConsensusInfo

        field = GroupField()
        cons = ConsensusInfo("specimen-1.v1", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "group=1"

    def test_variant_field_matches_full(self):
        """VariantField should return variant=full for .full names."""
        from speconsense.summarize.fields import VariantField
        from speconsense.types import ConsensusInfo

        field = VariantField()
        cons = ConsensusInfo("specimen-1.full", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "variant=full"

    def test_variant_field_still_matches_v(self):
        """VariantField should still match standard .v names."""
        from speconsense.summarize.fields import VariantField
        from speconsense.types import ConsensusInfo

        field = VariantField()
        cons = ConsensusInfo("specimen-1.v2", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "variant=v2"

    def test_variant_field_still_matches_raw(self):
        """VariantField should still match .raw names."""
        from speconsense.summarize.fields import VariantField
        from speconsense.types import ConsensusInfo

        field = VariantField()
        cons = ConsensusInfo("specimen-1.v1.raw2", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "variant=v1"


class TestSelectMinSizeRatio:
    """Tests for --select-min-size-ratio filtering."""

    def test_select_min_size_ratio_filters_small_variants(self):
        """Integration test: --select-min-size-ratio 0.1 filters out tiny variants."""
        temp_dir = tempfile.mkdtemp()
        source_dir = os.path.join(temp_dir, "clusters")
        summary_dir = os.path.join(temp_dir, "__Summary__")
        os.makedirs(source_dir)

        try:
            seq1 = "ATCGATCGATCGATCGATCGATCG"
            seq2 = "ATCGATCGATCAATCGATCGATCG"  # One SNP — different enough to not merge

            fasta_content = f""">test-c1 size=100 ric=100 primers=test
{seq1}
>test-c2 size=3 ric=3 primers=test
{seq2}
"""
            fasta_file = os.path.join(source_dir, "test-all.fasta")
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)

            result = subprocess.run(
                [
                    "speconsense-summarize",
                    "--source", source_dir,
                    "--summary-dir", summary_dir,
                    "--min-ric", "3",
                    "--select-min-size-ratio", "0.1",
                    "--disable-merging",
                ],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

            output_fasta = os.path.join(summary_dir, "summary.fasta")
            output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

            # Only the large variant should remain (3/100 = 0.03 < 0.1)
            assert len(output_sequences) == 1, \
                f"Expected 1 sequence after filtering, got {len(output_sequences)}"
            assert "size=100" in output_sequences[0].description

        finally:
            shutil.rmtree(temp_dir)

    def test_select_min_size_ratio_disabled_preserves_all(self):
        """Integration test: --select-min-size-ratio 0 preserves all variants."""
        temp_dir = tempfile.mkdtemp()
        source_dir = os.path.join(temp_dir, "clusters")
        summary_dir = os.path.join(temp_dir, "__Summary__")
        os.makedirs(source_dir)

        try:
            seq1 = "ATCGATCGATCGATCGATCGATCG"
            seq2 = "ATCGATCGATCAATCGATCGATCG"  # One SNP

            fasta_content = f""">test-c1 size=100 ric=100 primers=test
{seq1}
>test-c2 size=3 ric=3 primers=test
{seq2}
"""
            fasta_file = os.path.join(source_dir, "test-all.fasta")
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)

            result = subprocess.run(
                [
                    "speconsense-summarize",
                    "--source", source_dir,
                    "--summary-dir", summary_dir,
                    "--min-ric", "3",
                    "--select-min-size-ratio", "0",
                    "--disable-merging",
                ],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

            output_fasta = os.path.join(summary_dir, "summary.fasta")
            output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

            # Both variants should be preserved
            assert len(output_sequences) == 2, \
                f"Expected 2 sequences with ratio=0, got {len(output_sequences)}"

        finally:
            shutil.rmtree(temp_dir)


class TestFullConsensusIntegration:
    """Integration test for --enable-full-consensus."""

    def test_full_consensus_appears_in_summary(self):
        """Test that .full sequences appear in summary.fasta with --enable-full-consensus."""
        # Create temporary directory structure
        temp_dir = tempfile.mkdtemp()
        source_dir = os.path.join(temp_dir, "clusters")
        summary_dir = os.path.join(temp_dir, "__Summary__")
        os.makedirs(source_dir)

        try:
            # Create two sequences that will be in the same group
            seq1 = "ATCGATCGATCGATCGATCGATCG"
            seq2 = "ATCGATCGATCAATCGATCGATCG"  # One SNP at position 12

            fasta_content = f""">test-c1 size=10 ric=10 primers=test
{seq1}
>test-c2 size=8 ric=8 primers=test
{seq2}
"""
            fasta_file = os.path.join(source_dir, "test-all.fasta")
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)

            result = subprocess.run(
                [
                    "speconsense-summarize",
                    "--source", source_dir,
                    "--summary-dir", summary_dir,
                    "--min-ric", "3",
                    "--enable-full-consensus",
                    "--min-merge-overlap", "0",
                ],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

            output_fasta = os.path.join(summary_dir, "summary.fasta")
            assert os.path.exists(output_fasta)

            output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

            # Find the .full sequence
            full_seqs = [s for s in output_sequences if '.full' in s.id]
            assert len(full_seqs) == 1, f"Expected 1 .full sequence, got {len(full_seqs)}"

            full_seq = full_seqs[0]
            # .full should be at least as long as the longest input
            assert len(full_seq.seq) >= max(len(seq1), len(seq2))

            # Verify no FASTQ for .full
            fastq_dir = os.path.join(summary_dir, "FASTQ Files")
            if os.path.exists(fastq_dir):
                fastq_files = os.listdir(fastq_dir)
                full_fastqs = [f for f in fastq_files if '.full' in f]
                assert len(full_fastqs) == 0, f"Should be no FASTQ for .full, got: {full_fastqs}"

        finally:
            shutil.rmtree(temp_dir)
