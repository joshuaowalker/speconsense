"""Tests for the critical error rate (CER) statistical significance module.

Verifies compute_critical_error_rate, compute_minimum_M, and is_variant_significant
against values from the variant significance paper (docs/variant_significance/).
Also tests CER parsing and filtering in the summarize pipeline.
"""

import os
import tempfile
import pytest
from speconsense.significance import (
    SIGNAL_DESTRUCTION_THRESHOLD,
    compute_critical_error_rate,
    compute_minimum_M,
    is_cer_reportable,
    is_variant_significant,
)
from speconsense.summarize.io import parse_consensus_header, load_consensus_sequences


# ==============================================================================
# Tests for compute_critical_error_rate() — verify against paper Table 1
# ==============================================================================

class TestComputeCriticalErrorRate:
    """Verify p* values against the variant significance paper."""

    def test_standard_operating_point_uniform(self):
        """N=1000, M=100, alpha=0.05 => p*≈20.25% (uniform model, Table 1)."""
        p_star = compute_critical_error_rate(N=1000, M=100, L=700, alpha=0.05)
        assert abs(p_star - 0.2025) < 0.005, f"Expected ~20.25%, got {p_star*100:.2f}%"

    def test_run_corrected_uniform(self):
        """N=1000, M=100, alpha=1e-5 => p*≈16.59% (uniform model, Table 1)."""
        p_star = compute_critical_error_rate(N=1000, M=100, L=700, alpha=1e-5)
        assert abs(p_star - 0.1659) < 0.005, f"Expected ~16.59%, got {p_star*100:.2f}%"

    def test_moderate_coverage_uniform(self):
        """N=100, M=10, alpha=0.05 => p*≈6.55% (uniform model, Table 1)."""
        p_star = compute_critical_error_rate(N=100, M=10, L=700, alpha=0.05)
        assert abs(p_star - 0.0655) < 0.005, f"Expected ~6.55%, got {p_star*100:.2f}%"

    def test_moderate_run_corrected_uniform(self):
        """N=100, M=10, alpha=1e-5 => p*≈2.50% (uniform model, Table 1)."""
        p_star = compute_critical_error_rate(N=100, M=10, L=700, alpha=1e-5)
        assert abs(p_star - 0.0250) < 0.005, f"Expected ~2.50%, got {p_star*100:.2f}%"

    def test_edge_M_zero(self):
        """M=0 should return 0.0."""
        assert compute_critical_error_rate(N=1000, M=0, L=700) == 0.0

    def test_edge_M_negative(self):
        """M<0 should return 0.0."""
        assert compute_critical_error_rate(N=1000, M=-1, L=700) == 0.0

    def test_edge_N_zero(self):
        """N=0 should return 0.0."""
        assert compute_critical_error_rate(N=0, M=10, L=700) == 0.0

    def test_edge_L_zero(self):
        """L=0 should return 0.0."""
        assert compute_critical_error_rate(N=1000, M=10, L=0) == 0.0

    def test_edge_M_equals_N(self):
        """M=N should give p* > 1.0 (uncapped, beyond coherent range)."""
        p_star = compute_critical_error_rate(N=100, M=100, L=700, alpha=0.05)
        assert p_star > 1.0  # Uncapped — CER not applicable at this M/N ratio

    def test_monotonicity_in_M(self):
        """p* should increase with M (more reads = stronger evidence)."""
        prev_p_star = 0.0
        for M in [5, 10, 20, 50, 100, 200]:
            p_star = compute_critical_error_rate(N=1000, M=M, L=700, alpha=0.05)
            assert p_star > prev_p_star, f"p* should increase with M, but {p_star} <= {prev_p_star} at M={M}"
            prev_p_star = p_star

    def test_monotonicity_in_N(self):
        """p* should decrease with N for fixed M (same count in larger pool)."""
        prev_p_star = float('inf')
        for N in [100, 200, 500, 1000]:
            p_star = compute_critical_error_rate(N=N, M=20, L=700, alpha=0.05)
            assert p_star < prev_p_star, f"p* should decrease with N, but {p_star} >= {prev_p_star} at N={N}"
            prev_p_star = p_star

    def test_M_one(self):
        """M=1 should give very low p* (single read is weak evidence)."""
        p_star = compute_critical_error_rate(N=1000, M=1, L=700, alpha=0.05)
        assert p_star < 0.01  # Should be well below 1%

    def test_no_cap_at_high_M(self):
        """p* should not be capped — raw solver value returned for diagnostics."""
        p_star = compute_critical_error_rate(N=1000, M=500, L=700, alpha=1e-5)
        assert p_star > 1.0  # Raw value exceeds 1.0


# ==============================================================================
# Tests for is_cer_reportable() — signal destruction threshold
# ==============================================================================

class TestIsCerReportable:
    """CER is meaningful only when p* < 0.75 (signal destruction threshold).

    At p* = 0.75, each of the 4 bases is equally likely under the uniform
    error model — the reference population equals the variant count.
    """

    def test_low_p_star(self):
        """p*=0.10 is well below threshold — reportable."""
        assert is_cer_reportable(0.10) is True

    def test_at_threshold(self):
        """p*=0.75 is at threshold — not reportable."""
        assert is_cer_reportable(0.75) is False

    def test_above_threshold(self):
        """p*=0.90 is above threshold — not reportable."""
        assert is_cer_reportable(0.90) is False

    def test_well_above_threshold(self):
        """p*=1.5 (uncapped solver output) — not reportable."""
        assert is_cer_reportable(1.5) is False

    def test_zero(self):
        """p*=0.0 is reportable."""
        assert is_cer_reportable(0.0) is True

    def test_just_below_threshold(self):
        """p*=0.749 — reportable."""
        assert is_cer_reportable(0.749) is True

    def test_large_N_boundary(self):
        """N=1000, M=329: p* crosses 0.75, should not be reportable."""
        p_star = compute_critical_error_rate(N=1000, M=329, L=700, alpha=1e-5)
        assert is_cer_reportable(p_star) is False

    def test_large_N_below_boundary(self):
        """N=1000, M=328: p* just below 0.75, should be reportable."""
        p_star = compute_critical_error_rate(N=1000, M=328, L=700, alpha=1e-5)
        assert is_cer_reportable(p_star) is True

    def test_small_N_all_reportable(self):
        """At N=10, p* never reaches 0.75 — all variants are reportable."""
        for M in range(1, 11):
            p_star = compute_critical_error_rate(N=10, M=M, L=700, alpha=1e-5)
            assert is_cer_reportable(p_star) is True, f"N=10, M={M}: p*={p_star}"


# ==============================================================================
# Tests for compute_minimum_M() — verify against paper Table 2
# ==============================================================================

class TestComputeMinimumM:
    """Verify minimum M thresholds against the variant significance paper Table 2.

    Table 2 values are for uniform model (q = p/3).
    """

    # Uniform model, p=1%
    def test_N1000_p1pct_alpha005(self):
        """N=1000, p=1%, alpha=0.05 => M_min=13 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=0.05, assumed_error_rate=0.01) == 13

    def test_N1000_p1pct_alpha1e3(self):
        """N=1000, p=1%, alpha=1e-3 => M_min=16 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-3, assumed_error_rate=0.01) == 16

    def test_N1000_p1pct_alpha1e5(self):
        """N=1000, p=1%, alpha=1e-5 => M_min=19 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.01) == 19

    # Uniform model, p=3%
    def test_N1000_p3pct_alpha005(self):
        """N=1000, p=3%, alpha=0.05 => M_min=25 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=0.05, assumed_error_rate=0.03) == 25

    def test_N1000_p3pct_alpha1e5(self):
        """N=1000, p=3%, alpha=1e-5 => M_min=33 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.03) == 33

    # Uniform model, p=5%
    def test_N1000_p5pct_alpha005(self):
        """N=1000, p=5%, alpha=0.05 => M_min=35 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=0.05, assumed_error_rate=0.05) == 35

    def test_N1000_p5pct_alpha1e5(self):
        """N=1000, p=5%, alpha=1e-5 => M_min=44 (Table 2, uniform)."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.05) == 44

    # Smaller N values
    def test_N100_p1pct_alpha005(self):
        """N=100, p=1%, alpha=0.05 => M_min=5 (Table 2, uniform)."""
        assert compute_minimum_M(N=100, L=700, alpha=0.05, assumed_error_rate=0.01) == 5

    def test_N100_p3pct_alpha1e5(self):
        """N=100, p=3%, alpha=1e-5 => M_min=11 (Table 2, uniform)."""
        assert compute_minimum_M(N=100, L=700, alpha=1e-5, assumed_error_rate=0.03) == 11

    def test_N500_p1pct_alpha1e5(self):
        """N=500, p=1%, alpha=1e-5 => M_min=14 (Table 2, uniform)."""
        assert compute_minimum_M(N=500, L=700, alpha=1e-5, assumed_error_rate=0.01) == 14

    def test_error_rate_zero_returns_one(self):
        """assumed_error_rate=0 should return 1 (feature disabled)."""
        assert compute_minimum_M(N=1000, L=700, alpha=0.05, assumed_error_rate=0) == 1


# ==============================================================================
# Tests for is_variant_significant()
# ==============================================================================

class TestIsVariantSignificant:
    """Test the combined significance decision function."""

    def test_well_supported_variant(self):
        """M=100 out of N=1000 should be significant at 2% error rate."""
        is_sig, p_star = is_variant_significant(M=100, N=1000, L=700,
                                                 assumed_error_rate=0.02)
        assert is_sig is True
        assert p_star > 0.02

    def test_weak_variant_suppressed(self):
        """M=5 out of N=1000 should NOT be significant at 2% error rate."""
        is_sig, p_star = is_variant_significant(M=5, N=1000, L=700,
                                                 assumed_error_rate=0.02)
        assert is_sig is False
        assert p_star < 0.02

    def test_disabled_when_error_rate_zero(self):
        """assumed_error_rate=0 should always return significant."""
        is_sig, p_star = is_variant_significant(M=1, N=1000, L=700,
                                                 assumed_error_rate=0)
        assert is_sig is True

    def test_disabled_when_error_rate_negative(self):
        """Negative assumed_error_rate should always return significant."""
        is_sig, p_star = is_variant_significant(M=1, N=1000, L=700,
                                                 assumed_error_rate=-1)
        assert is_sig is True

    def test_returns_p_star(self):
        """Should return p_star value even when feature is disabled."""
        _, p_star = is_variant_significant(M=50, N=1000, L=700,
                                            assumed_error_rate=0)
        assert p_star > 0

    def test_large_M_still_significant(self):
        """Large M: p* exceeds assumed error rate, so significant."""
        is_sig, p_star = is_variant_significant(M=300, N=1000, L=700,
                                                 assumed_error_rate=0.02)
        assert is_sig is True
        assert p_star > 0.02  # p* is computed (not None), just above threshold

    def test_large_M_p_star_above_threshold(self):
        """Large M: p* above signal destruction threshold but still returned."""
        is_sig, p_star = is_variant_significant(M=500, N=1000, L=700,
                                                 assumed_error_rate=0.02)
        assert is_sig is True
        assert p_star >= SIGNAL_DESTRUCTION_THRESHOLD  # Above 0.75

    def test_small_N_large_fraction_not_auto_passed(self):
        """N=10, M=3 (30%): must be tested, not auto-passed."""
        is_sig, p_star = is_variant_significant(M=3, N=10, L=700,
                                                 assumed_error_rate=0.02)
        assert is_sig is False  # p*=0.0015 < 0.02
        assert p_star < 0.02

    def test_small_N_significant_variant(self):
        """N=10, M=5 (50%): genuinely significant."""
        is_sig, p_star = is_variant_significant(M=5, N=10, L=700,
                                                 assumed_error_rate=0.02)
        assert is_sig is True
        assert p_star > 0.02

    def test_boundary_case(self):
        """Test at the boundary: minimum M for 2% at N=1000."""
        min_M = compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.02)

        # At min_M, should be significant
        is_sig, _ = is_variant_significant(M=min_M, N=1000, L=700,
                                            assumed_error_rate=0.02, alpha=1e-5)
        assert is_sig is True

        # At min_M - 1, should not be significant
        is_sig, _ = is_variant_significant(M=min_M - 1, N=1000, L=700,
                                            assumed_error_rate=0.02, alpha=1e-5)
        assert is_sig is False


# ==============================================================================
# Tests for CER header parsing
# ==============================================================================

class TestCERHeaderParsing:
    """Test parsing of cer= and cer.a= fields from FASTA headers."""

    def test_parse_cer_fields(self):
        """Parse cer and cer.a from header."""
        header = ">sample-c1 size=100 ric=50 cer=0.15 cer.a=1e-05"
        result = parse_consensus_header(header)
        assert result[9] == pytest.approx(0.15)  # cer
        assert result[10] == pytest.approx(1e-5)  # cer_alpha

    def test_parse_no_cer_fields(self):
        """Missing cer fields should return None."""
        header = ">sample-c1 size=100 ric=50"
        result = parse_consensus_header(header)
        assert result[9] is None  # cer
        assert result[10] is None  # cer_alpha

    def test_cer_not_confused_with_cer_alpha(self):
        """cer= regex should not match cer.a= prefix."""
        header = ">sample-c1 size=100 ric=50 cer.a=1e-05"
        result = parse_consensus_header(header)
        assert result[9] is None  # cer should be None (not matched from cer.a)
        assert result[10] == pytest.approx(1e-5)  # cer_alpha


# ==============================================================================
# Tests for CER value loading (no filtering at load time)
# ==============================================================================

class TestCERLoading:
    """Test that CER values are loaded but not filtered in load_consensus_sequences.

    CER filtering now happens after HAC grouping in process_single_specimen,
    not at load time. Load should preserve all variants regardless of CER.
    """

    @pytest.fixture
    def fasta_dir(self, tmp_path):
        """Create a temporary directory with test FASTA files."""
        fasta_content = (
            ">specimen-c1 size=200 ric=100 cer=0.15 cer.a=1e-05\n"
            "ACGTACGTACGT\n"
            ">specimen-c2 size=50 ric=25 cer=0.01 cer.a=1e-05\n"
            "ACGTACGTACGT\n"
            ">specimen-c3 size=30 ric=15\n"
            "ACGTACGTACGT\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        return tmp_path

    def test_all_variants_loaded(self, fasta_dir):
        """All variants should load regardless of CER value."""
        result = load_consensus_sequences(str(fasta_dir), min_ric=1)
        assert len(result) == 3

    def test_cer_values_preserved(self, fasta_dir):
        """CER values should be preserved in ConsensusInfo."""
        result = load_consensus_sequences(str(fasta_dir), min_ric=1)
        by_name = {c.sample_name: c for c in result}
        assert by_name["specimen-c1"].cer == pytest.approx(0.15)
        assert by_name["specimen-c1"].cer_alpha == pytest.approx(1e-5)
        assert by_name["specimen-c2"].cer == pytest.approx(0.01)
        assert by_name["specimen-c3"].cer is None
        assert by_name["specimen-c3"].cer_alpha is None


# ==============================================================================
# Tests for CER-based filtering in process_single_specimen (post-HAC)
# ==============================================================================

class TestCERFilteringPostHAC:
    """Test CER-based variant filtering after HAC grouping.

    CER filtering protects the primary (largest) variant in each group and
    only filters secondary variants with cer < assumed_error_rate. This
    prevents eliminating the dominant sequence when there's no stronger
    competing signal.
    """

    @pytest.fixture
    def fasta_dir_multi(self, tmp_path):
        """Create test FASTA with a primary and secondary variants.

        All sequences are identical so they land in a single HAC group.
        """
        seq = "ACGT" * 50  # 200bp identical sequence
        fasta_content = (
            f">specimen-c1 size=200 ric=100 cer=0.15 cer.a=1e-05\n"
            f"{seq}\n"
            f">specimen-c2 size=50 ric=25 cer=0.01 cer.a=1e-05\n"
            f"{seq}\n"
            f">specimen-c3 size=30 ric=15\n"
            f"{seq}\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        return tmp_path

    @pytest.fixture
    def fasta_dir_single(self, tmp_path):
        """Create test FASTA with a single low-CER variant."""
        seq = "ACGT" * 50
        fasta_content = (
            f">specimen-c1 size=5 ric=5 cer=0.019 cer.a=1e-05\n"
            f"{seq}\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        return tmp_path

    @pytest.fixture
    def fasta_dir_primary_low_cer(self, tmp_path):
        """Create test FASTA where the primary variant has low CER."""
        seq = "ACGT" * 50
        fasta_content = (
            f">specimen-c1 size=100 ric=50 cer=0.015 cer.a=1e-05\n"
            f"{seq}\n"
            f">specimen-c2 size=10 ric=10 cer=0.005 cer.a=1e-05\n"
            f"{seq}\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        return tmp_path

    def _run_summarize(self, fasta_dir, assumed_error_rate):
        """Run process_single_specimen with given CER filter rate."""
        from speconsense.summarize.io import load_consensus_sequences
        from speconsense.summarize.cli import process_single_specimen
        from argparse import Namespace

        consensuses = load_consensus_sequences(str(fasta_dir), min_ric=1)
        args = Namespace(
            group_identity=0.95,
            min_merge_overlap=0,
            select_max_groups=0,
            disable_merging=True,
            select_min_size_ratio=0,
            select_max_variants=0,
            select_strategy='size',
            scale_threshold=0,
            threads=1,
            source=str(fasta_dir),
            enable_full_consensus=False,
        )
        final, _, _, _, _ = process_single_specimen(
            consensuses, args, cer_filter_rate=assumed_error_rate
        )
        return final

    def test_secondary_low_cer_filtered(self, fasta_dir_multi):
        """Secondary variant with cer < threshold should be filtered."""
        result = self._run_summarize(fasta_dir_multi, assumed_error_rate=0.02)
        # process_single_specimen renames variants, so check by size
        sizes = sorted([c.size for c in result], reverse=True)
        assert 200 in sizes  # Primary (cer=0.15, passes)
        assert 50 not in sizes  # Secondary (cer=0.01, filtered)
        assert 30 in sizes  # No CER header, never filtered
        assert len(result) == 2

    def test_primary_never_filtered(self, fasta_dir_primary_low_cer):
        """Primary variant should never be filtered, even with low CER."""
        result = self._run_summarize(fasta_dir_primary_low_cer, assumed_error_rate=0.02)
        sizes = [c.size for c in result]
        assert 100 in sizes  # Primary protected despite cer=0.015

    def test_secondary_with_low_cer_filtered_but_primary_kept(self, fasta_dir_primary_low_cer):
        """Secondary filtered, primary kept even when both have low CER."""
        result = self._run_summarize(fasta_dir_primary_low_cer, assumed_error_rate=0.02)
        sizes = [c.size for c in result]
        assert 100 in sizes  # Primary (largest)
        assert 10 not in sizes  # Secondary (cer=0.005 < 0.02)
        assert len(result) == 1

    def test_single_variant_never_filtered(self, fasta_dir_single):
        """Single-variant specimen should never be filtered by CER."""
        result = self._run_summarize(fasta_dir_single, assumed_error_rate=0.02)
        assert len(result) == 1
        assert result[0].size == 5

    def test_cer_filter_disabled(self, fasta_dir_multi):
        """CER filter disabled when rate=0 — all variants pass."""
        result = self._run_summarize(fasta_dir_multi, assumed_error_rate=0)
        assert len(result) == 3

    def test_no_cer_header_never_filtered(self, fasta_dir_multi):
        """Variants without CER header are never filtered."""
        result = self._run_summarize(fasta_dir_multi, assumed_error_rate=0.50)
        sizes = [c.size for c in result]
        # c3 (size=30) has no CER, should survive even at high threshold
        assert 30 in sizes
