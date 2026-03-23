"""Tests for the critical error rate (CER) statistical significance module.

Verifies compute_critical_error_rate, compute_minimum_M, and is_variant_significant
against values from the variant significance paper (docs/variant_significance/).
Also tests CER parsing and filtering in the summarize pipeline.
"""

import os
import tempfile
import pytest
from speconsense.significance import (
    compute_critical_error_rate,
    compute_minimum_M,
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
        """M=N should give high p* (all reads support variant)."""
        p_star = compute_critical_error_rate(N=100, M=100, L=700, alpha=0.05)
        assert p_star == 1.0  # Capped at 1.0 (obviously not error)

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
# Tests for CER-based filtering in load_consensus_sequences
# ==============================================================================

class TestCERFiltering:
    """Test CER-based variant filtering in load_consensus_sequences."""

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

    def test_cer_filter_removes_low_cer(self, fasta_dir):
        """Variants with cer < assumed_error_rate should be filtered."""
        result = load_consensus_sequences(
            str(fasta_dir), min_ric=1, assumed_error_rate=0.02
        )
        # c1 (cer=0.15) passes, c2 (cer=0.01) filtered, c3 (no cer) passes
        names = [c.sample_name for c in result]
        assert "specimen-c1" in names
        assert "specimen-c2" not in names
        assert "specimen-c3" in names
        assert len(result) == 2

    def test_cer_filter_disabled_by_default(self, fasta_dir):
        """Default assumed_error_rate=0 should not filter anything."""
        result = load_consensus_sequences(str(fasta_dir), min_ric=1)
        assert len(result) == 3

    def test_no_cer_header_not_filtered(self, fasta_dir):
        """Variants without cer header should never be filtered by CER."""
        result = load_consensus_sequences(
            str(fasta_dir), min_ric=1, assumed_error_rate=0.10
        )
        # c1 (cer=0.15) passes, c2 (cer=0.01) filtered, c3 (no cer) passes
        names = [c.sample_name for c in result]
        assert "specimen-c3" in names

    def test_cer_values_preserved(self, fasta_dir):
        """CER values should be preserved in ConsensusInfo."""
        result = load_consensus_sequences(str(fasta_dir), min_ric=1)
        by_name = {c.sample_name: c for c in result}
        assert by_name["specimen-c1"].cer == pytest.approx(0.15)
        assert by_name["specimen-c1"].cer_alpha == pytest.approx(1e-5)
        assert by_name["specimen-c3"].cer is None
        assert by_name["specimen-c3"].cer_alpha is None
