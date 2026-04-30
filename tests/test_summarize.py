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
    9 clusters with core gid=1 (all same-primer, grouped into one identity group).

    The file contains:
    - c1: main cluster (ric=500) - majority pattern
    - c2: second major cluster (ric=250) - majority pattern
    - c3: contamination (ric=9)
    - c4: variant (ric=6) - ends with TAG, structural variant
    - c5: variant (ric=6) - structural variant
    - c6: contamination (ric=4)
    - c7: variant (ric=3) - ends with TAA
    - c8: variant (ric=3) - majority pattern
    - c9: variant (ric=3) - ends with TAA, homopolymer variation from c7

    Expected behavior after switching to core-provided grouping:
    - c1 + c2 merge into 1.v1 with rawric=500+250 (primary merge)
    - Contaminations and other variants stay separate
    - c3, c6 should not merge into the majority pattern
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
                "--min-merge-overlap", "0",  # Disable overlap merge for this test
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

        # Primary cluster pair (c1 + c2) must still merge: same primers,
        # homopolymer-equivalent consensus, dominating this specimen.
        c1_c2_merged = False
        for seq in output_sequences:
            if 'rawric=' in seq.description:
                rawric_match = re.search(r'rawric=([\d+]+)', seq.description)
                if rawric_match:
                    ric_values = [int(x) for x in rawric_match.group(1).split('+')]
                    if sorted(ric_values, reverse=True)[:2] == [500, 250]:
                        c1_c2_merged = True
                        break
        assert c1_c2_merged, \
            "c1 (ric=500) and c2 (ric=250) should merge into the primary variant"

        # The specimen should produce a compact output (not every cluster
        # appearing as its own variant). At least one merge must have occurred.
        assert len(output_sequences) < 9, \
            f"Expected some merging to reduce from 9 input clusters, got {len(output_sequences)}"

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)


def test_merge_with_homopolymer_only_differences():
    """Test that sequences differing only in LONG homopolymer lengths DO merge.

    Under the default --hp-normalization-length=6, only HP runs whose shortest
    side is >= 6 are blanket-normalized. Short-HP length diffs now count as
    real edits. This test exercises the long-HP regime where merging still
    collapses length variation.
    """
    # Create temporary directory structure
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # Create two sequences that differ only in homopolymer length.
        # Both runs are >= 6 so min(L1,L2)=8 passes the default threshold.
        seq1 = "ATCGAAAAAAAATCGATCGATCGATCG"        # 8 A's
        seq2 = "ATCGAAAAAAAAAATCGATCGATCGATCG"      # 10 A's

        fasta_content = f""">test-seq-1.v1 size=10 ric=10 primers=test gid=1 vid=1
{seq1}
>test-seq-1.v2 size=8 ric=8 primers=test gid=1 vid=2
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


def test_lq_emitted_when_specimen_has_no_passing_variants():
    """Regression: specimens whose every variant is routed to .ns/.lq must
    still get their filtered-variant files written.

    Earlier the per-specimen processing loop iterated only over file paths
    present in the passing list, so an all-filtered specimen was skipped
    entirely and its .ns/.lq output was lost.
    """
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # One cluster, err_factor above the default 1.5 threshold → routed to .lq.
        # No other variants in this specimen, so file_groups is empty for it.
        fasta_content = (
            ">single-cluster-1.v1 size=3 ric=3 primers=test "
            "err_factor=1.675 gid=1 vid=1\n"
            "ACGTACGTACGTACGTACGTACGTACGTACGT\n"
        )
        fasta_file = os.path.join(source_dir, "single-cluster-all.fasta")
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)

        result = subprocess.run(
            [
                "speconsense-summarize",
                "--source", source_dir,
                "--summary-dir", summary_dir,
                "--min-ric", "3",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

        # The .lq file must exist for this specimen even though no variants passed.
        lq_files = [
            f for f in os.listdir(os.path.join(summary_dir, "variants"))
            if f.startswith("single-cluster-1.v1.lq-RiC")
        ]
        assert lq_files, (
            "Expected single-cluster-1.v1.lq-RiC*.fasta in variants/, "
            f"but found: {os.listdir(os.path.join(summary_dir, 'variants'))}"
        )

    finally:
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


class TestFieldRegex:
    """Tests for GroupField and VariantField regex on standard variant names."""

    def test_group_field_matches_variant(self):
        """GroupField should match standard .v names."""
        from speconsense.summarize.fields import GroupField
        from speconsense.types import ConsensusInfo

        field = GroupField()
        cons = ConsensusInfo("specimen-1.v1", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "group=1"

    def test_variant_field_matches_v(self):
        """VariantField should match standard .v names."""
        from speconsense.summarize.fields import VariantField
        from speconsense.types import ConsensusInfo

        field = VariantField()
        cons = ConsensusInfo("specimen-1.v2", "c1", "ACGT", ric=10, size=10, file_path="f")
        assert field.format_value(cons) == "variant=v2"

    def test_variant_field_matches_raw(self):
        """VariantField should match .raw names."""
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

            fasta_content = f""">test-1.v1 size=100 ric=100 primers=test gid=1 vid=1
{seq1}
>test-1.v2 size=3 ric=3 primers=test gid=1 vid=2
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

            fasta_content = f""">test-1.v1 size=100 ric=100 primers=test gid=1 vid=1
{seq1}
>test-1.v2 size=3 ric=3 primers=test gid=1 vid=2
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


class TestGroupByCoreIdentity:
    """group_by_core_identity should bucket by gid and hard-fail on missing."""

    def _make(self, name, group_rank, variant_rank=1, size=10):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name, cluster_id=name.rsplit('-', 1)[-1],
            sequence="ACGT", ric=size, size=size, file_path="f",
            group_rank=group_rank, variant_rank=variant_rank,
        )

    def test_buckets_by_group_rank(self):
        from speconsense.summarize.clustering import group_by_core_identity
        members = [
            self._make("s-c1", 1, 1, 100),
            self._make("s-c2", 2, 1, 50),
            self._make("s-c3", 1, 2, 25),
        ]
        groups = group_by_core_identity(members)
        assert set(groups.keys()) == {1, 2}
        assert len(groups[1]) == 2
        assert len(groups[2]) == 1

    def test_hard_fail_on_missing_gid(self):
        import pytest
        from speconsense.summarize.clustering import group_by_core_identity
        # One member lacks group_rank
        members = [
            self._make("s-c1", 1, 1, 100),
            self._make("s-c2", None, None, 50),
        ]
        with pytest.raises(ValueError, match="gid="):
            group_by_core_identity(members)

    def test_empty_returns_empty(self):
        from speconsense.summarize.clustering import group_by_core_identity
        assert group_by_core_identity([]) == {}


class TestMergeGroupsByAnchorOverlap:
    """merge_groups_by_anchor_overlap should conflate cross-primer groups."""

    def _make(self, name, seq, size, primers, group_rank):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name, cluster_id="c1", sequence=seq,
            ric=size, size=size, file_path="f", primers=primers,
            group_rank=group_rank, variant_rank=1,
        )

    def test_merges_cross_primer_overlap(self):
        from speconsense.summarize.clustering import merge_groups_by_anchor_overlap
        shared = "ACGT" * 100  # 400bp
        groups = {
            1: [self._make("s-1", shared + "TAAA" * 50, 100, ["P1"], 1)],
            2: [self._make("s-2", "AAA" * 20 + shared, 50, ["P2"], 2)],
        }
        result = merge_groups_by_anchor_overlap(groups, min_overlap_bp=200,
                                                group_identity=0.85)
        assert len(result) == 1, f"Expected merge, got {len(result)} groups"
        survivor_rank = next(iter(result))
        assert survivor_rank == 1, "Larger group should win"
        assert len(result[1]) == 2

    def test_skips_same_primer_pairs(self):
        from speconsense.summarize.clustering import merge_groups_by_anchor_overlap
        shared = "ACGT" * 100
        groups = {
            1: [self._make("s-1", shared, 100, ["P1"], 1)],
            2: [self._make("s-2", shared, 50, ["P1"], 2)],
        }
        result = merge_groups_by_anchor_overlap(groups, min_overlap_bp=200,
                                                group_identity=0.85)
        assert len(result) == 2, "Same-primer groups must not merge cross-primer"

    def test_noop_when_overlap_disabled(self):
        from speconsense.summarize.clustering import merge_groups_by_anchor_overlap
        groups = {
            1: [self._make("s-1", "A" * 400, 100, ["P1"], 1)],
            2: [self._make("s-2", "A" * 400, 50, ["P2"], 2)],
        }
        result = merge_groups_by_anchor_overlap(groups, min_overlap_bp=0,
                                                group_identity=0.85)
        assert len(result) == 2, "min_overlap_bp=0 disables merging"


class TestParseConsensusHeaderGidVid:
    """parse_consensus_header should extract the core-emitted gid/vid tokens."""

    def test_gid_vid_present(self):
        from speconsense.summarize.io import parse_consensus_header
        header = ">specimen-2.v1 size=100 ric=50 gid=2 vid=1"
        result = parse_consensus_header(header)
        assert len(result) == 13
        group_rank = result[11]
        variant_rank = result[12]
        assert group_rank == 2
        assert variant_rank == 1

    def test_gid_vid_absent_returns_none(self):
        from speconsense.summarize.io import parse_consensus_header
        header = ">specimen-c1 size=100 ric=50"
        result = parse_consensus_header(header)
        assert result[11] is None
        assert result[12] is None

    def test_gid_vid_coexist_with_cer_and_err(self):
        from speconsense.summarize.io import parse_consensus_header
        header = (">specimen-c1 size=100 ric=50 cer_factor=3.200 "
                  "err_factor=1.100 gid=3 vid=2")
        result = parse_consensus_header(header)
        sample_name = result[0]
        cer_factor = result[9]
        err_factor = result[10]
        group_rank = result[11]
        variant_rank = result[12]
        assert sample_name == "specimen-c1"
        assert cer_factor == 3.2
        assert err_factor == 1.1
        assert group_rank == 3
        assert variant_rank == 2
