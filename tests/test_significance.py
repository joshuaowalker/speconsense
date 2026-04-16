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
# Tests for K>1 multi-position CER
# ==============================================================================

class TestMultiPositionCER:
    """Verify K>1 p* values against the reference implementation
    (docs/variant_significance/compute_critical_rates.py).
    """

    def test_k1_backward_compatible(self):
        """K=1 should produce identical results to default (no K argument)."""
        p_default = compute_critical_error_rate(N=1000, M=100, L=700, alpha=0.05)
        p_k1 = compute_critical_error_rate(N=1000, M=100, L=700, alpha=0.05, K=1)
        assert p_default == p_k1

    def test_k2_standard_point(self):
        """N=1000, M=100, K=2, alpha=0.05 => p*≈72.57% (uniform)."""
        p_star = compute_critical_error_rate(N=1000, M=100, L=700, alpha=0.05, K=2)
        assert abs(p_star - 0.7257) < 0.005, f"Expected ~72.57%, got {p_star*100:.2f}%"

    def test_k3_standard_point(self):
        """N=1000, M=100, K=3, alpha=0.05 => p*≈112.18% (uniform, uncapped)."""
        p_star = compute_critical_error_rate(N=1000, M=100, L=700, alpha=0.05, K=3)
        assert abs(p_star - 1.1218) < 0.01, f"Expected ~112.18%, got {p_star*100:.2f}%"

    def test_k2_run_corrected(self):
        """N=1000, M=100, K=2, alpha=1e-5 => p*≈66.63% (uniform)."""
        p_star = compute_critical_error_rate(N=1000, M=100, L=700, alpha=1e-5, K=2)
        assert abs(p_star - 0.6663) < 0.005, f"Expected ~66.63%, got {p_star*100:.2f}%"

    def test_k2_small_N(self):
        """N=100, M=10, K=2, alpha=1e-5 => p*≈20.13% (uniform)."""
        p_star = compute_critical_error_rate(N=100, M=10, L=700, alpha=1e-5, K=2)
        assert abs(p_star - 0.2013) < 0.005, f"Expected ~20.13%, got {p_star*100:.2f}%"

    def test_monotonicity_in_K(self):
        """p* should increase with K (more positions = stronger evidence)."""
        p_k1 = compute_critical_error_rate(N=1000, M=50, L=700, alpha=1e-5, K=1)
        p_k2 = compute_critical_error_rate(N=1000, M=50, L=700, alpha=1e-5, K=2)
        p_k3 = compute_critical_error_rate(N=1000, M=50, L=700, alpha=1e-5, K=3)
        assert p_k1 < p_k2 < p_k3, f"Expected p* to increase with K: {p_k1}, {p_k2}, {p_k3}"

    def test_minimum_M_k1_paper_table3(self):
        """Table 3: N=1000, K=1, p=1.5%, alpha=1e-5 => M_min=23."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.015, K=1) == 23

    def test_minimum_M_k2_paper_table3(self):
        """Table 3: N=1000, K=2, p=1.5%, alpha=1e-5 => M_min=6."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.015, K=2) == 6

    def test_minimum_M_k3_paper_table3(self):
        """Table 3: N=1000, K=3, p=1.5%, alpha=1e-5 => M_min=4."""
        assert compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.015, K=3) == 4

    def test_is_variant_significant_k2(self):
        """M=6 at K=2 should be significant (passes at K=2 but not K=1)."""
        is_sig_k1, _ = is_variant_significant(M=6, N=1000, L=700,
                                               assumed_error_rate=0.015, K=1)
        is_sig_k2, _ = is_variant_significant(M=6, N=1000, L=700,
                                               assumed_error_rate=0.015, K=2)
        assert is_sig_k1 is False
        assert is_sig_k2 is True

    def test_edge_K_zero(self):
        """K=0 should return 0.0."""
        assert compute_critical_error_rate(N=1000, M=100, L=700, K=0) == 0.0

    def test_edge_K_exceeds_L(self):
        """K > L should return 0.0 (can't have more positions than amplicon length)."""
        assert compute_critical_error_rate(N=1000, M=100, L=700, K=701) == 0.0

    def test_minimum_M_decreases_with_K(self):
        """Minimum M should decrease with K at fixed N, p, alpha."""
        m1 = compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.015, K=1)
        m2 = compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.015, K=2)
        m3 = compute_minimum_M(N=1000, L=700, alpha=1e-5, assumed_error_rate=0.015, K=3)
        assert m1 > m2 > m3, f"Expected M to decrease: K=1:{m1}, K=2:{m2}, K=3:{m3}"


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
    """Test parsing of cer=, cer.a=, and cer.ns fields from FASTA headers."""

    def test_parse_cer_fields(self):
        """Parse cer and cer.a from header."""
        header = ">sample-c1 size=100 ric=50 cer=0.15 cer.a=1e-05"
        result = parse_consensus_header(header)
        assert result[9] == pytest.approx(0.15)  # cer
        assert result[10] == pytest.approx(1e-5)  # cer_alpha
        assert result[11] is False  # cer_ns

    def test_parse_no_cer_fields(self):
        """Missing cer fields should return None/False."""
        header = ">sample-c1 size=100 ric=50"
        result = parse_consensus_header(header)
        assert result[9] is None  # cer
        assert result[10] is None  # cer_alpha
        assert result[11] is False  # cer_ns

    def test_cer_not_confused_with_cer_alpha(self):
        """cer= regex should not match cer.a= prefix."""
        header = ">sample-c1 size=100 ric=50 cer.a=1e-05"
        result = parse_consensus_header(header)
        assert result[9] is None  # cer should be None (not matched from cer.a)
        assert result[10] == pytest.approx(1e-5)  # cer_alpha
        assert result[11] is False  # cer_ns

    def test_parse_cer_ns_tag(self):
        """cer.ns tag should be detected."""
        header = ">sample-c3 size=5 ric=5 cer=0.0012 cer.ns cer.a=1e-05"
        result = parse_consensus_header(header)
        assert result[9] == pytest.approx(0.0012)  # cer (numeric value still parsed)
        assert result[10] == pytest.approx(1e-5)  # cer_alpha
        assert result[11] is True  # cer_ns

    def test_parse_cer_anchor(self):
        """cer=anchor should not parse as numeric cer, and is not cer.ns."""
        header = ">sample-c1 size=200 ric=200 cer=anchor"
        result = parse_consensus_header(header)
        assert result[9] is None  # cer (anchor is not numeric)
        assert result[11] is False  # cer_ns


# ==============================================================================
# Tests for CER value loading (no filtering at load time)
# ==============================================================================

class TestCERLoading:
    """Test CER loading and cer.ns filtering in load_consensus_sequences.

    Variants marked cer.ns by core are filtered at load time.
    All other variants (numeric cer, cer=anchor, no cer) are loaded.
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

    def test_all_non_ns_variants_loaded(self, fasta_dir):
        """Variants without cer.ns should load."""
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

    def test_cer_ns_variant_filtered_at_load(self, tmp_path):
        """Variants with cer.ns should be excluded at load time."""
        fasta_content = (
            ">specimen-c1 size=200 ric=100 cer=anchor\n"
            "ACGTACGTACGT\n"
            ">specimen-c2 size=50 ric=25 cer=0.15 cer.a=1e-05\n"
            "ACGTACGTACGT\n"
            ">specimen-c3 size=5 ric=5 cer=0.0012 cer.ns cer.a=1e-05\n"
            "ACGTACGTACGT\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        result = load_consensus_sequences(str(tmp_path), min_ric=1)
        names = [c.sample_name for c in result]
        assert len(result) == 2
        assert "specimen-c1" in names  # anchor loads
        assert "specimen-c2" in names  # pass loads
        assert "specimen-c3" not in names  # ns filtered


# ==============================================================================
# Tests for CER ns filtering at load time
# ==============================================================================

class TestCERNsFilterAtLoad:
    """Test that cer.ns variants from core are filtered at load time.

    CER computation is done exclusively in the core pipeline. Summarize
    trusts core's annotations: variants marked cer.ns are excluded at load
    time; all other variants pass through.
    """

    @pytest.fixture
    def fasta_dir_mixed(self, tmp_path):
        """Create test FASTA with a mix of CER statuses.

        - c1: anchor (cer=anchor)
        - c2: pass (numeric cer)
        - c3: not significant (cer.ns)
        - c4: no CER annotations (e.g., from older core version)
        """
        fasta_content = (
            ">specimen-c1 size=200 ric=200 cer=anchor\n"
            "ACGTACGTACGT\n"
            ">specimen-c2 size=50 ric=50 cer=0.1234 cer.a=1e-05\n"
            "ACGTACGTACGT\n"
            ">specimen-c3 size=5 ric=5 cer=0.0012 cer.ns cer.a=1e-05\n"
            "ACGTACGTACGT\n"
            ">specimen-c4 size=30 ric=30\n"
            "ACGTACGTACGT\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        return tmp_path

    def test_ns_excluded_others_loaded(self, fasta_dir_mixed):
        """Only cer.ns variant should be excluded."""
        result = load_consensus_sequences(str(fasta_dir_mixed), min_ric=1)
        names = [c.sample_name for c in result]
        assert len(result) == 3
        assert "specimen-c1" in names
        assert "specimen-c2" in names
        assert "specimen-c3" not in names
        assert "specimen-c4" in names

    def test_anchor_loads_with_no_numeric_cer(self, fasta_dir_mixed):
        """cer=anchor should load with cer=None (not numeric)."""
        result = load_consensus_sequences(str(fasta_dir_mixed), min_ric=1)
        by_name = {c.sample_name: c for c in result}
        assert by_name["specimen-c1"].cer is None

    def test_pass_preserves_cer_value(self, fasta_dir_mixed):
        """Numeric cer value should be preserved on passing variants."""
        result = load_consensus_sequences(str(fasta_dir_mixed), min_ric=1)
        by_name = {c.sample_name: c for c in result}
        assert by_name["specimen-c2"].cer == pytest.approx(0.1234)
        assert by_name["specimen-c2"].cer_alpha == pytest.approx(1e-5)

    def test_no_annotation_loads_normally(self, fasta_dir_mixed):
        """Variants without any CER annotation should load."""
        result = load_consensus_sequences(str(fasta_dir_mixed), min_ric=1)
        by_name = {c.sample_name: c for c in result}
        assert by_name["specimen-c4"].cer is None
        assert by_name["specimen-c4"].cer_alpha is None

    def test_no_ns_variants_all_pass(self, tmp_path):
        """When no cer.ns tags exist, all variants load (backward compatible)."""
        fasta_content = (
            ">specimen-c1 size=200 ric=100\n"
            "ACGTACGTACGT\n"
            ">specimen-c2 size=50 ric=25\n"
            "ACGTACGTACGT\n"
        )
        fasta_file = tmp_path / "specimen-all.fasta"
        fasta_file.write_text(fasta_content)
        result = load_consensus_sequences(str(tmp_path), min_ric=1)
        assert len(result) == 2
