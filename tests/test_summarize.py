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


class TestProcessSingleSpecimenNaming:
    """process_single_specimen should preserve core gid/vid except across cross-primer conflation."""

    @staticmethod
    def _make(name, group_rank, variant_rank, size,
              sequence=None, primers=None):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name,
            cluster_id=f"{group_rank}.v{variant_rank}",
            sequence=sequence or ("ACGT" * 100),
            ric=size,
            size=size,
            file_path="/tmp/test-all.fasta",
            primers=primers or ["P1"],
            group_rank=group_rank,
            variant_rank=variant_rank,
        )

    @staticmethod
    def _args(**overrides):
        from types import SimpleNamespace
        defaults = dict(
            min_merge_overlap=0,
            group_identity=0.85,
            hp_normalization_length=6,
            select_max_groups=-1,
            select_min_size_ratio=0,
            select_max_variants=-1,
            disable_merging=True,
            enable_full_consensus=False,
            prune_group_frac=0.10,
            prune_group_abs=15,
        )
        defaults.update(overrides)
        return SimpleNamespace(**defaults)

    def test_no_conflation_preserves_core_names(self):
        from speconsense.summarize.cli import process_single_specimen
        members = [
            self._make("test-1.v1", 1, 1, 100),
            self._make("test-1.v2", 1, 2, 50),
            self._make("test-1.v3", 1, 3, 10),
        ]
        final, _, _, _, _, _, _, _ = process_single_specimen(members, self._args())
        names = {v.sample_name for v in final}
        assert names == {"test-1.v1", "test-1.v2", "test-1.v3"}

    def test_gaps_from_size_ratio_filter_preserve_remaining_vids(self):
        from speconsense.summarize.cli import process_single_specimen
        members = [
            self._make("test-1.v1", 1, 1, 100),
            self._make("test-1.v2", 1, 2, 80),
            self._make("test-1.v3", 1, 3, 5),  # filtered by ratio
        ]
        args = self._args(select_min_size_ratio=0.5)
        final, _, _, _, _, _, _, _ = process_single_specimen(members, args)
        names = sorted(v.sample_name for v in final)
        assert names == ["test-1.v1", "test-1.v2"]

    def test_two_group_conflation_allocates_new_vids_above_max(self):
        from speconsense.summarize.cli import process_single_specimen
        shared = "ACGT" * 100  # 400bp common
        members = [
            self._make("test-1.v1", 1, 1, 100, shared + "TAAA" * 50, ["P1"]),
            self._make("test-1.v2", 1, 2, 80, shared + "TAAA" * 50, ["P1"]),
            self._make("test-2.v1", 2, 1, 50, "AAA" * 20 + shared, ["P2"]),
        ]
        args = self._args(min_merge_overlap=200)
        final, _, _, _, _, _, _, _ = process_single_specimen(members, args)
        names = sorted(v.sample_name for v in final)
        # gid=1 keeps v1/v2 verbatim; gid=2's v1 is moved to v3 under gid=1
        assert names == ["test-1.v1", "test-1.v2", "test-1.v3"]

    def test_two_group_conflation_skips_ns_vid(self):
        from speconsense.summarize.cli import process_single_specimen
        shared = "ACGT" * 100
        members = [
            self._make("test-1.v1", 1, 1, 100, shared + "TAAA" * 50, ["P1"]),
            self._make("test-2.v1", 2, 1, 50, "AAA" * 20 + shared, ["P2"]),
        ]
        ns_records = [self._make("test-1.v2", 1, 2, 5)]  # vid=2 already used in gid=1
        args = self._args(min_merge_overlap=200)
        final, _, _, _, _, _, _, _ = process_single_specimen(
            members, args, ns_for_specimen=ns_records)
        names = sorted(v.sample_name for v in final)
        # moved record gets v3, skipping v2 occupied by ns
        assert names == ["test-1.v1", "test-1.v3"]

    def test_two_group_conflation_skips_lq_vid(self):
        from speconsense.summarize.cli import process_single_specimen
        shared = "ACGT" * 100
        members = [
            self._make("test-1.v1", 1, 1, 100, shared + "TAAA" * 50, ["P1"]),
            self._make("test-2.v1", 2, 1, 50, "AAA" * 20 + shared, ["P2"]),
        ]
        lq_records = [self._make("test-1.v4", 1, 4, 5)]  # lq occupies vid=4
        args = self._args(min_merge_overlap=200)
        final, _, _, _, _, _, _, _ = process_single_specimen(
            members, args, lq_for_specimen=lq_records)
        names = sorted(v.sample_name for v in final)
        # moved record gets v5 (max(used={1,4}) + 1), not v2/v3 (gaps in core)
        assert names == ["test-1.v1", "test-1.v5"]

    def test_three_group_conflation_assigns_unique_vids(self):
        from speconsense.summarize.cli import process_single_specimen
        shared = "ACGT" * 100  # common backbone all three pairs share
        members = [
            self._make("test-1.v1", 1, 1, 100, shared + "TAAA" * 50, ["P1"]),
            self._make("test-2.v1", 2, 1, 80, "AAA" * 20 + shared, ["P2"]),
            self._make("test-2.v2", 2, 2, 60, "AAA" * 20 + shared, ["P2"]),
            self._make("test-3.v1", 3, 1, 40, shared + "GGG" * 50, ["P3"]),
        ]
        args = self._args(min_merge_overlap=200)
        final, _, _, _, _, _, _, _ = process_single_specimen(members, args)
        names = sorted(v.sample_name for v in final)
        assert len(names) == 4
        # All emit under gid=1
        assert all(n.startswith("test-1.v") for n in names)
        # All vids unique, gid=1's original v1 preserved, moved members
        # allocated sequentially above max
        vids = sorted(int(n.rsplit(".v", 1)[1]) for n in names)
        assert vids == [1, 2, 3, 4]

    def test_absorbed_groups_ns_lq_do_not_pollute_survivor_namespace(self):
        from speconsense.summarize.cli import process_single_specimen
        shared = "ACGT" * 100
        members = [
            self._make("test-1.v1", 1, 1, 100, shared + "TAAA" * 50, ["P1"]),
            self._make("test-2.v1", 2, 1, 50, "AAA" * 20 + shared, ["P2"]),
        ]
        # ns/lq under gid=2 should NOT enter gid=1's collision set
        ns_records = [self._make("test-2.v2", 2, 2, 5)]
        lq_records = [self._make("test-2.v3", 2, 3, 5)]
        args = self._args(min_merge_overlap=200)
        final, _, _, _, _, _, _, _ = process_single_specimen(
            members, args, ns_for_specimen=ns_records,
            lq_for_specimen=lq_records)
        names = sorted(v.sample_name for v in final)
        # moved record gets v2 (next free under gid=1 — gid=2's ns/lq vids
        # remain in gid=2 on disk, not blocking gid=1)
        assert names == ["test-1.v1", "test-1.v2"]


class TestMergeFieldHandling:
    """_build_merged_consensus_info should aggregate fields correctly."""

    @staticmethod
    def _make(name, size, primers, cer_factor=None, err_factor=None,
              group_rank=1, variant_rank=1):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name, cluster_id=f"{group_rank}.v{variant_rank}",
            sequence="ACGT" * 50, ric=size, size=size,
            file_path="f", primers=primers,
            cer_factor=cer_factor, err_factor=err_factor,
            group_rank=group_rank, variant_rank=variant_rank,
        )

    def test_primers_union_across_contributors(self):
        from speconsense.summarize.merging import _build_merged_consensus_info
        variants = [
            self._make("a", 100, ["5'-ITS1F", "3'-ITS4_RC"]),
            self._make("b", 50, ["5'-ITS3", "3'-ITS4_RC"]),
        ]
        merged = _build_merged_consensus_info(list("ACGT"), 0, variants)
        # Union of primer names, sorted
        assert merged.primers == sorted({
            "5'-ITS1F", "5'-ITS3", "3'-ITS4_RC"
        })

    def test_primers_single_pair_preserves_pair(self):
        from speconsense.summarize.merging import _build_merged_consensus_info
        primers = ["5'-ITS1F", "3'-ITS4_RC"]
        variants = [
            self._make("a", 100, primers),
            self._make("b", 50, primers),
        ]
        merged = _build_merged_consensus_info(list("ACGT"), 0, variants)
        assert sorted(merged.primers) == sorted(primers)

    def test_inherits_cer_err_factor_pre_recompute(self):
        # _build_merged_consensus_info hands off to the recompute pass which
        # may overwrite these; here we verify the inheritance happens at all,
        # so the recompute pass has something to start from.
        from speconsense.summarize.merging import _build_merged_consensus_info
        variants = [
            self._make("a", 100, ["P1"], cer_factor=3.2, err_factor=1.05),
            self._make("b", 50, ["P1"], cer_factor=8.0, err_factor=2.0),
        ]
        merged = _build_merged_consensus_info(list("ACGT"), 0, variants)
        # Largest contributor wins for these inherited fields.
        assert merged.cer_factor == 3.2
        assert merged.err_factor == 1.05


class TestRefreshCerFactor:
    """_refresh_cer_factor should distinguish same-primer vs cross-primer merges."""

    @staticmethod
    def _make(name, size, sequence, primers, group_rank=1, variant_rank=1):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name, cluster_id=f"{group_rank}.v{variant_rank}",
            sequence=sequence, ric=size, size=size,
            file_path="f", primers=primers,
            group_rank=group_rank, variant_rank=variant_rank,
        )

    def test_cross_primer_merge_sets_cer_to_none(self):
        from speconsense.summarize.merging import _refresh_cer_factor
        seq = "ACGT" * 50
        merged = self._make("merged", 150, seq, ["P1", "P2", "P3"])
        contributors = [
            self._make("a", 100, seq, ["P1", "P2"]),
            self._make("b", 50, seq, ["P3"]),
        ]
        peer = self._make("peer", 500, seq, ["P1"])
        bucket = [peer, merged]
        # Even with qctx_table available, cross-primer ⇒ None
        result = _refresh_cer_factor(
            merged, contributors, bucket,
            qctx_table={"non-hp-sub": 0.006, "non-hp-indel": 0.011},
            alpha=1e-5,
            hp_min_length=6,
        )
        assert result.cer_factor is None

    def test_same_primer_no_qctx_returns_unchanged(self):
        from speconsense.summarize.merging import _refresh_cer_factor
        # When qctx_table is None, recompute can't run; inherited value
        # passes through untouched.
        seq = "ACGT" * 50
        merged = self._make("merged", 150, seq, ["P1"])._replace(cer_factor=3.0)
        contributors = [
            self._make("a", 100, seq, ["P1"]),
            self._make("b", 50, seq, ["P1"]),
        ]
        peer = self._make("peer", 500, seq, ["P1"])
        bucket = [peer, merged]
        result = _refresh_cer_factor(
            merged, contributors, bucket,
            qctx_table=None, alpha=1e-5, hp_min_length=6,
        )
        assert result.cer_factor == 3.0

    def test_same_primer_no_larger_peer_returns_none(self):
        from speconsense.summarize.merging import _refresh_cer_factor
        # Merged record is the anchor of its bucket; no peer to compare
        # against ⇒ cer_factor None.
        seq = "ACGT" * 50
        merged = self._make("merged", 500, seq, ["P1"])._replace(cer_factor=3.0)
        contributors = [
            self._make("a", 400, seq, ["P1"]),
            self._make("b", 100, seq, ["P1"]),
        ]
        smaller = self._make("smaller", 50, seq, ["P1"])
        bucket = [merged, smaller]
        result = _refresh_cer_factor(
            merged, contributors, bucket,
            qctx_table={"non-hp-sub": 0.006, "non-hp-indel": 0.011},
            alpha=1e-5,
            hp_min_length=6,
        )
        assert result.cer_factor is None

    def test_same_primer_with_peer_recomputes(self):
        from speconsense.summarize.merging import _refresh_cer_factor
        # Construct a peer that differs from merged at one position so
        # classify_pairwise_differences returns K=1 and the recompute runs.
        merged_seq = "ACGT" * 50
        peer_seq = "ACGT" * 50
        peer_seq = peer_seq[:10] + "T" + peer_seq[11:]  # single substitution
        merged = self._make("merged", 50, merged_seq, ["P1"])._replace(cer_factor=999.0)
        contributors = [
            self._make("a", 30, merged_seq, ["P1"]),
            self._make("b", 20, merged_seq, ["P1"]),
        ]
        peer = self._make("peer", 500, peer_seq, ["P1"])
        bucket = [peer, merged]
        result = _refresh_cer_factor(
            merged, contributors, bucket,
            qctx_table={
                "non-hp-sub": 0.006, "non-hp-indel": 0.011,
                "hp-l1": 0.006, "hp-l2": 0.008, "hp-l3": 0.010,
                "hp-l4": 0.011, "hp-l5": 0.011,
            },
            alpha=1e-5,
            hp_min_length=6,
        )
        # The inherited 999.0 should be replaced by a recomputed value.
        # Exact value depends on solver; just assert it changed and is finite.
        assert result.cer_factor != 999.0
        assert result.cer_factor is None or 0.0 < result.cer_factor < float('inf')


class TestRawFieldPropagation:
    """.raw record construction should preserve every field from the source cluster."""

    def test_raw_preserves_cer_err_gid_vid_and_err_sums(self):
        # _build .raw records via the same _replace call that
        # write_specimen_data_files uses, and verify the new field
        # propagation. Direct test against the construction logic (not the
        # full file-writing path) — sufficient since the fix is a one-line
        # _replace.
        from speconsense.types import ConsensusInfo
        raw_info = ConsensusInfo(
            sample_name="specimen-1.v3",
            cluster_id="1.v3",
            sequence="ACGT",
            ric=42,
            size=42,
            file_path="f",
            primers=["P1"],
            rid=0.985,
            rid_min=0.971,
            cer_factor=4.2,
            err_factor=1.15,
            err_factor_obs_sum=0.5,
            err_factor_exp_sum=0.43,
            err_factor_cols=120,
            group_rank=1,
            variant_rank=3,
        )
        # Match the construction in write_specimen_data_files
        raw = raw_info._replace(
            sample_name="specimen-1.v1.raw2",
            snp_count=None,
            raw_ric=None,
            raw_len=None,
            merge_indel_count=None,
        )
        # Merge-only fields reset
        assert raw.snp_count is None
        assert raw.raw_ric is None
        assert raw.raw_len is None
        assert raw.merge_indel_count is None
        # Source cluster fields preserved
        assert raw.cer_factor == 4.2
        assert raw.err_factor == 1.15
        assert raw.err_factor_obs_sum == 0.5
        assert raw.err_factor_exp_sum == 0.43
        assert raw.err_factor_cols == 120
        assert raw.group_rank == 1
        assert raw.variant_rank == 3
        assert raw.rid == 0.985
        assert raw.rid_min == 0.971
        assert raw.primers == ["P1"]


class TestFullConsensus:
    """Unit tests for --enable-full-consensus group consensus builder."""

    @staticmethod
    def _variant(name, size, primers=None):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name,
            cluster_id=name.rsplit("-", 1)[-1],
            sequence="ACGT" * 50,
            ric=size,
            size=size,
            file_path="/tmp/test-all.fasta",
            primers=primers or ["P1"],
            group_rank=1,
            variant_rank=1,
        )

    def test_gate_admits_strong_minor_signal(self):
        from speconsense.summarize.full_consensus import _gate_variants
        variants = [
            self._variant("a", 100),
            self._variant("b", 30),
            self._variant("c", 20),
        ]
        # Sorted desc: 100; 30 >= 0.10 * 100; 20 >= 0.10 * 130 ⇒ all three
        gated = _gate_variants(variants, 0.10)
        assert [v.size for v in gated] == [100, 30, 20]

    def test_gate_excludes_noise_tail(self):
        from speconsense.summarize.full_consensus import _gate_variants
        variants = [
            self._variant("a", 100),
            self._variant("b", 30),
            self._variant("c", 8),
            self._variant("d", 2),
        ]
        # 30 passes; 8 < 0.10 * 130 = 13 ⇒ stop
        gated = _gate_variants(variants, 0.10)
        assert [v.size for v in gated] == [100, 30]

    def test_gate_handles_single_variant(self):
        from speconsense.summarize.full_consensus import _gate_variants
        gated = _gate_variants([self._variant("a", 100)], 0.10)
        assert len(gated) == 1

    def test_allocation_sums_to_budget(self):
        from speconsense.summarize.full_consensus import _allocate_reads
        alloc = _allocate_reads([80, 20], 100)
        assert sum(alloc) == 100
        # Size-weighted
        assert alloc == [80, 20]

    def test_allocation_rounds_fairly(self):
        from speconsense.summarize.full_consensus import _allocate_reads
        # 100 budget across [33, 33, 34] sizes: raw=[33,33,34], floored sums to 100 — no remainder
        alloc = _allocate_reads([33, 33, 34], 100)
        assert sum(alloc) == 100
        # 70 budget across [50, 50, 50] (total=150 > budget): raw=[23.3, 23.3, 23.3],
        # floored=[23,23,23] sums to 69, remainder=1 → goes to one of them
        alloc = _allocate_reads([50, 50, 50], 70)
        assert sum(alloc) == 70
        assert sorted(alloc) == [23, 23, 24]

    def test_allocation_caps_at_available(self):
        from speconsense.summarize.full_consensus import _allocate_reads
        # Total available smaller than budget: return sizes verbatim
        alloc = _allocate_reads([10, 5], 100)
        assert alloc == [10, 5]

    def test_cross_primer_detection(self):
        from speconsense.summarize.full_consensus import _is_cross_primer
        same = [self._variant("a", 100, ["P1"]), self._variant("b", 50, ["P1"])]
        cross = [self._variant("a", 100, ["P1"]), self._variant("b", 50, ["P2"])]
        assert not _is_cross_primer(same)
        assert _is_cross_primer(cross)

    def test_consensus_majority_gap_omits_position(self):
        from speconsense.summarize.merging import build_full_consensus_from_msa
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # Column 0: 4 gap, 2 A → omitted
        # Column 1: 2 gap, 4 A → kept as A
        rows = [
            SeqRecord(Seq("-A"), id="r1"),
            SeqRecord(Seq("-A"), id="r2"),
            SeqRecord(Seq("-A"), id="r3"),
            SeqRecord(Seq("-A"), id="r4"),
            SeqRecord(Seq("A-"), id="r5"),
            SeqRecord(Seq("A-"), id="r6"),
        ]
        consensus, snp = build_full_consensus_from_msa(rows, 0.10)
        assert consensus == "A"
        assert snp == 0

    def test_consensus_iupac_at_threshold(self):
        from speconsense.summarize.merging import build_full_consensus_from_msa
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # 9 A, 1 T at col 0 → T is at 10%, meets threshold → IUPAC W
        rows = [SeqRecord(Seq("A"), id=f"r{i}") for i in range(9)]
        rows.append(SeqRecord(Seq("T"), id="rT"))
        consensus, snp = build_full_consensus_from_msa(rows, 0.10)
        assert consensus == "W"
        assert snp == 1

    def test_consensus_iupac_below_threshold(self):
        from speconsense.summarize.merging import build_full_consensus_from_msa
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # 19 A, 1 T at col 0 → T is at 5%, below 10% threshold → plain A
        rows = [SeqRecord(Seq("A"), id=f"r{i}") for i in range(19)]
        rows.append(SeqRecord(Seq("T"), id="rT"))
        consensus, snp = build_full_consensus_from_msa(rows, 0.10)
        assert consensus == "A"
        assert snp == 0

    def test_consensus_three_way_iupac(self):
        from speconsense.summarize.merging import build_full_consensus_from_msa
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # 5 A, 3 C, 2 G → all above 10%: IUPAC V (A/C/G)
        rows = [SeqRecord(Seq(b), id=f"r{i}") for i, b in enumerate("AAAAACCCGG")]
        consensus, snp = build_full_consensus_from_msa(rows, 0.10)
        assert consensus == "V"
        assert snp == 1

    def test_pass_track_guard_blocks_single_variant(self):
        from speconsense.summarize.full_consensus import build_group_full_consensus
        result = build_group_full_consensus(
            final_gid=1,
            group_members=[self._variant("a", 100)],
            pass_track_count=1,  # only one variant on pass track
            fastq_lookup={},
            min_ambiguity_frequency=0.10,
            max_sample_size=100,
            specimen_base="test",
            primary_file_path="/tmp/test-all.fasta",
        )
        assert result is None

    def test_single_gate_guard_blocks_emission(self):
        from speconsense.summarize.full_consensus import build_group_full_consensus
        # Two pass-track variants but the smaller doesn't clear the gate
        result = build_group_full_consensus(
            final_gid=1,
            group_members=[self._variant("a", 100), self._variant("b", 2)],
            pass_track_count=2,
            fastq_lookup={},
            min_ambiguity_frequency=0.10,
            max_sample_size=100,
            specimen_base="test",
            primary_file_path="/tmp/test-all.fasta",
        )
        assert result is None

    def test_no_fastq_lookup_returns_none(self):
        from speconsense.summarize.full_consensus import build_group_full_consensus
        # Gate would pass but no reads can be loaded
        result = build_group_full_consensus(
            final_gid=1,
            group_members=[self._variant("a", 100), self._variant("b", 50)],
            pass_track_count=2,
            fastq_lookup=None,
            min_ambiguity_frequency=0.10,
            max_sample_size=100,
            specimen_base="test",
            primary_file_path="/tmp/test-all.fasta",
        )
        assert result is None

    def test_strip_cluster_suffix_handles_full(self):
        from speconsense.summarize.io import strip_cluster_suffix
        assert strip_cluster_suffix("specimen-1-full") == "specimen"
        assert strip_cluster_suffix("specimen-1.v2") == "specimen"
        assert strip_cluster_suffix("specimen-1.v2.raw3") == "specimen"
        assert strip_cluster_suffix("specimen-1-full.raw1") == "specimen"
        assert strip_cluster_suffix("plain_name") == "plain_name"

    def test_end_to_end_emits_full_record(self):
        """Integration smoke: with real cluster_debug FASTQs and SPOA available,
        process_single_specimen should emit a -{gid}-full record."""
        import tempfile
        from speconsense.types import ConsensusInfo
        from speconsense.summarize.cli import process_single_specimen

        with tempfile.TemporaryDirectory() as tmpdir:
            # Synthesize cluster_debug FASTQ files for two variants of the same gid.
            # Both share a backbone; v2 has a single SNP at position 50.
            backbone = "ACGT" * 40  # 160bp
            v1_seq = backbone
            v2_seq = backbone[:50] + "C" + backbone[51:]  # G→C at position 50

            cluster_dir = os.path.join(tmpdir, "cluster_debug")
            os.makedirs(cluster_dir)

            def write_reads(path, seq, count, ric_label):
                with open(path, "w") as f:
                    for i in range(count):
                        f.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")

            v1_fastq = os.path.join(cluster_dir, f"specimen-1.v1-RiC60-reads.fastq")
            v2_fastq = os.path.join(cluster_dir, f"specimen-1.v2-RiC20-reads.fastq")
            write_reads(v1_fastq, v1_seq, 60, "RiC60")
            write_reads(v2_fastq, v2_seq, 20, "RiC20")

            fastq_lookup = {"specimen": [v1_fastq, v2_fastq]}

            members = [
                ConsensusInfo(
                    sample_name="specimen-1.v1",
                    cluster_id="1.v1",
                    sequence=v1_seq,
                    ric=60, size=60,
                    file_path=os.path.join(tmpdir, "specimen-all.fasta"),
                    primers=["P1"],
                    group_rank=1, variant_rank=1,
                ),
                ConsensusInfo(
                    sample_name="specimen-1.v2",
                    cluster_id="1.v2",
                    sequence=v2_seq,
                    ric=20, size=20,
                    file_path=os.path.join(tmpdir, "specimen-all.fasta"),
                    primers=["P1"],
                    group_rank=1, variant_rank=2,
                ),
            ]

            from types import SimpleNamespace
            args = SimpleNamespace(
                min_merge_overlap=0,
                group_identity=0.85,
                hp_normalization_length=6,
                select_max_groups=-1,
                select_min_size_ratio=0,
                select_max_variants=-1,
                disable_merging=True,
                enable_full_consensus=True,
                prune_group_frac=0.10,
                prune_group_abs=15,
            )
            final, _, _, _, _, full_reads, _, _ = process_single_specimen(
                members, args,
                fastq_lookup=fastq_lookup,
                full_min_ambiguity_frequency=0.10,
                full_max_sample_size=50,
            )

            full_records = [c for c in final if c.sample_name.endswith("-full")]
            assert len(full_records) == 1, f"expected 1 -full record, got {[c.sample_name for c in final]}"
            full = full_records[0]
            assert full.sample_name == "specimen-1-full"
            assert full.group_rank == 1
            assert full.variant_rank is None
            assert full.size == 80  # 60 + 20
            # Sampled reads accessible for FASTQ emission
            assert full.sample_name in full_reads
            assert len(full_reads[full.sample_name]) > 0
            # Sequence should incorporate IUPAC at the SNP position (G/C → S)
            assert "S" in full.sequence, f"expected S code in {full.sequence}"


class TestFrequencyFields:
    """Tests for group_frequency= and global_frequency= FASTA header fields."""

    @staticmethod
    def _make(name, gid, vid, size, sequence=None, primers=None):
        from speconsense.types import ConsensusInfo
        return ConsensusInfo(
            sample_name=name,
            cluster_id=f"{gid}.v{vid}",
            sequence=sequence or "ACGT" * 50,
            ric=size,
            size=size,
            file_path="/tmp/test-all.fasta",
            primers=primers or ["P1"],
            group_rank=gid,
            variant_rank=vid,
        )

    @staticmethod
    def _args(**overrides):
        from types import SimpleNamespace
        defaults = dict(
            min_merge_overlap=0,
            group_identity=0.85,
            hp_normalization_length=6,
            select_max_groups=-1,
            select_min_size_ratio=0,
            select_max_variants=-1,
            disable_merging=True,
            enable_full_consensus=False,
            prune_group_frac=0.10,
            prune_group_abs=15,
        )
        defaults.update(overrides)
        return SimpleNamespace(**defaults)

    def test_group_frequency_format(self):
        from speconsense.summarize.fields import GroupFrequencyField
        c = self._make("test-1.v1", 1, 1, 30)._replace(group_size_total=100)
        assert GroupFrequencyField().format_value(c) == "group_frequency=30.0"

    def test_group_frequency_none_when_denom_missing(self):
        from speconsense.summarize.fields import GroupFrequencyField
        c = self._make("test-1.v1", 1, 1, 30)._replace(group_size_total=None)
        assert GroupFrequencyField().format_value(c) is None

    def test_group_frequency_none_when_denom_zero(self):
        from speconsense.summarize.fields import GroupFrequencyField
        c = self._make("test-1.v1", 1, 1, 30)._replace(group_size_total=0)
        assert GroupFrequencyField().format_value(c) is None

    def test_global_frequency_format(self):
        from speconsense.summarize.fields import GlobalFrequencyField
        c = self._make("test-1.v1", 1, 1, 30)._replace(global_size_total=1000)
        assert GlobalFrequencyField().format_value(c) == "global_frequency=3.0"

    def test_global_frequency_none_when_denom_missing(self):
        from speconsense.summarize.fields import GlobalFrequencyField
        c = self._make("test-1.v1", 1, 1, 30)._replace(global_size_total=None)
        assert GlobalFrequencyField().format_value(c) is None

    def test_full_preset_includes_both_fields(self):
        from speconsense.summarize.fields import FASTA_FIELD_PRESETS
        assert 'group_frequency' in FASTA_FIELD_PRESETS['full']
        assert 'global_frequency' in FASTA_FIELD_PRESETS['full']

    def test_default_preset_excludes_both_fields(self):
        from speconsense.summarize.fields import FASTA_FIELD_PRESETS
        assert 'group_frequency' not in FASTA_FIELD_PRESETS['default']
        assert 'global_frequency' not in FASTA_FIELD_PRESETS['default']
        assert 'group_frequency' not in FASTA_FIELD_PRESETS['qc']
        assert 'global_frequency' not in FASTA_FIELD_PRESETS['qc']

    def test_process_single_specimen_annotates_passed_variants(self):
        """Passed variants get group_size_total summed across their group
        including .ns/.lq records, plus the specimen's global_size_total."""
        from speconsense.summarize.cli import process_single_specimen
        members = [
            self._make("test-1.v1", 1, 1, 100),
            self._make("test-1.v2", 1, 2, 50),
        ]
        ns_records = [self._make("test-1.v3", 1, 3, 10)]
        lq_records = [self._make("test-1.v4", 1, 4, 5)]
        # Group total = 100 + 50 + 10 + 5 = 165
        # Global total provided as 2000.
        result = process_single_specimen(
            members, self._args(),
            ns_for_specimen=ns_records,
            lq_for_specimen=lq_records,
            specimen_global_size_total=2000,
        )
        final, _, _, _, _, _, annotated_ns, annotated_lq = result
        for v in final:
            assert v.group_size_total == 165
            assert v.global_size_total == 2000
        for v in annotated_ns:
            assert v.group_size_total == 165
            assert v.global_size_total == 2000
        for v in annotated_lq:
            assert v.group_size_total == 165
            assert v.global_size_total == 2000

    def test_process_single_specimen_cross_primer_uses_conflated_total(self):
        """Cross-primer conflation: the denominator includes records from
        every absorbed group, plus the .ns/.lq records under any of those
        gids. Verifies the conflation-aware bucket-total computation."""
        from speconsense.summarize.cli import process_single_specimen
        shared = "ACGT" * 100  # 400bp common anchor
        members = [
            self._make("test-1.v1", 1, 1, 100, shared + "TAAA" * 50, ["P1"]),
            self._make("test-2.v1", 2, 1, 30, "AAA" * 20 + shared, ["P2"]),
            self._make("test-2.v2", 2, 2, 20, "AAA" * 20 + shared, ["P2"]),
        ]
        # .ns record under gid=2 — should still contribute to the bucket
        # total of the conflated (1,2) bucket because gid=2 absorbs into gid=1.
        ns_records = [self._make("test-2.v3", 2, 3, 5, "AAA" * 20 + shared, ["P2"])]
        args = self._args(min_merge_overlap=200)
        result = process_single_specimen(
            members, args,
            ns_for_specimen=ns_records,
            specimen_global_size_total=500,
        )
        final, _, _, _, _, _, annotated_ns, _ = result
        # Bucket total = 100 + 30 + 20 + 5 = 155 — across all conflated gids
        for v in final:
            assert v.group_size_total == 155, \
                f"{v.sample_name}: expected 155, got {v.group_size_total}"
            assert v.global_size_total == 500
        for v in annotated_ns:
            assert v.group_size_total == 155
            assert v.global_size_total == 500

    def test_process_single_specimen_no_global_total(self):
        """When specimen metadata is missing, global_size_total stays None
        and global_frequency= is suppressed at format time."""
        from speconsense.summarize.cli import process_single_specimen
        members = [self._make("test-1.v1", 1, 1, 100), self._make("test-1.v2", 1, 2, 50)]
        result = process_single_specimen(
            members, self._args(), specimen_global_size_total=None,
        )
        final, _, _, _, _, _, _, _ = result
        for v in final:
            assert v.group_size_total == 150
            assert v.global_size_total is None
