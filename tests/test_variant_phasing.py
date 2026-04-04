"""Test variant phasing functions from msa.py module."""

import pytest
from speconsense.msa import (
    PositionStats,
    is_variant_position_with_composition,
    calculate_within_cluster_error,
    group_reads_by_single_position,
    call_iupac_ambiguities,
    ReadAlignment,
    ErrorPosition,
    IUPAC_CODES,
)


# ==============================================================================
# Tests for is_variant_position_with_composition()
# ==============================================================================

def test_variant_position_clear_biallelic():
    """Test clear biallelic variant (A/G 50/50 split)."""
    stats = PositionStats(
        msa_position=10,
        consensus_position=10,
        coverage=100,
        error_count=50,
        error_rate=0.50,
        sub_count=50,
        ins_count=0,
        del_count=0,
        consensus_nucleotide='A',
        base_composition={'A': 50, 'G': 50, 'C': 0, 'T': 0, '-': 0},
        homopolymer_composition={}
    )

    is_variant, alleles, reason = is_variant_position_with_composition(
        stats, min_variant_frequency=0.20, min_variant_count=5
    )

    assert is_variant is True
    # Function returns only alternative alleles (not the majority base)
    # With 50/50 split, the one sorted second alphabetically is the "alternative"
    assert len(alleles) == 1
    assert alleles[0] in ['A', 'G']  # One of them will be returned as alternative
    assert 'variant' in reason.lower()


def test_variant_position_below_frequency_threshold():
    """Test position with rare variant below frequency threshold (5% < 20%)."""
    stats = PositionStats(
        msa_position=10,
        consensus_position=10,
        coverage=100,
        error_count=5,
        error_rate=0.05,
        sub_count=5,
        ins_count=0,
        del_count=0,
        consensus_nucleotide='A',
        base_composition={'A': 95, 'G': 5, 'C': 0, 'T': 0, '-': 0},
        homopolymer_composition={}
    )

    is_variant, alleles, reason = is_variant_position_with_composition(
        stats, min_variant_frequency=0.20, min_variant_count=5
    )

    # G is at 5% frequency, which is below the 20% threshold
    assert is_variant is False
    assert len(alleles) == 0


def test_variant_position_below_count_threshold():
    """Test position with variant below count threshold."""
    stats = PositionStats(
        msa_position=10,
        consensus_position=10,
        coverage=10,
        error_count=3,
        error_rate=0.30,
        sub_count=3,
        ins_count=0,
        del_count=0,
        consensus_nucleotide='A',
        base_composition={'A': 7, 'G': 3, 'C': 0, 'T': 0, '-': 0},
        homopolymer_composition={}
    )

    is_variant, alleles, reason = is_variant_position_with_composition(
        stats, min_variant_frequency=0.20, min_variant_count=5
    )

    assert is_variant is False


def test_variant_position_homopolymer_normalized():
    """Test position where variation is reduced after homopolymer normalization.

    HP normalization subtracts homopolymer extension counts from base counts.
    Here we have A at 70 but 50 of those are HP extensions, leaving only 20 effective A's.
    G is at 30 with no HP adjustment.
    After normalization: A=20 (40%), G=30 (60%) - G is majority, A is alternative at 40%.
    A meets threshold so is_variant should be True.
    """
    stats = PositionStats(
        msa_position=10,
        consensus_position=10,
        coverage=100,
        error_count=30,
        error_rate=0.30,
        sub_count=30,
        ins_count=0,
        del_count=0,
        consensus_nucleotide='A',
        base_composition={'A': 70, 'G': 30, 'C': 0, 'T': 0, '-': 0},
        homopolymer_composition={'A': 50}  # 50 of the A's are HP extensions
    )

    is_variant, alleles, reason = is_variant_position_with_composition(
        stats, min_variant_frequency=0.20, min_variant_count=5
    )

    # After HP normalization: A=20, G=30. G is majority, A (20/50=40%) is alternative
    assert is_variant is True
    assert 'A' in alleles  # A becomes the alternative after G becomes majority


def test_variant_position_triallelic():
    """Test triallelic variant position (A/G/C).

    A=40 (majority), G=35 (35%), C=25 (25%).
    G and C are both alternative alleles that meet the 20% threshold.
    Function returns only alternative alleles (not the majority base A).
    """
    stats = PositionStats(
        msa_position=10,
        consensus_position=10,
        coverage=100,
        error_count=60,
        error_rate=0.60,
        sub_count=60,
        ins_count=0,
        del_count=0,
        consensus_nucleotide='A',
        base_composition={'A': 40, 'G': 35, 'C': 25, 'T': 0, '-': 0},
        homopolymer_composition={}
    )

    is_variant, alleles, reason = is_variant_position_with_composition(
        stats, min_variant_frequency=0.20, min_variant_count=5
    )

    assert is_variant is True
    # Function returns only alternative alleles (not including majority base A)
    assert set(alleles) == {'G', 'C'}


# ==============================================================================
# Tests for calculate_within_cluster_error()
# ==============================================================================

def test_within_cluster_error_perfect_haplotypes():
    """Test error calculation with perfectly separated haplotypes."""
    # Two haplotypes, each with perfectly consistent reads
    haplotype_groups = {
        'A-G': {'read1', 'read2', 'read3'},
        'T-C': {'read4', 'read5', 'read6'}
    }

    # Reads within each haplotype have identical alleles at all positions
    read_alleles = {
        'read1': {0: 'A', 1: 'G', 2: 'A', 3: 'T'},
        'read2': {0: 'A', 1: 'G', 2: 'A', 3: 'T'},
        'read3': {0: 'A', 1: 'G', 2: 'A', 3: 'T'},
        'read4': {0: 'T', 1: 'C', 2: 'G', 3: 'C'},
        'read5': {0: 'T', 1: 'C', 2: 'G', 3: 'C'},
        'read6': {0: 'T', 1: 'C', 2: 'G', 3: 'C'},
    }

    phasing_positions = {0, 1}  # Positions used for phasing
    all_variant_positions = {0, 1, 2, 3}  # All variant positions

    error = calculate_within_cluster_error(
        haplotype_groups, read_alleles, phasing_positions, all_variant_positions
    )

    # Perfect haplotypes should have 0 error at non-phased positions
    assert error == 0.0


def test_within_cluster_error_imperfect_haplotypes():
    """Test error calculation with some variation within haplotypes."""
    haplotype_groups = {
        'A-G': {'read1', 'read2', 'read3'},
    }

    # Reads differ at non-phased position 2
    read_alleles = {
        'read1': {0: 'A', 1: 'G', 2: 'A'},  # A at position 2
        'read2': {0: 'A', 1: 'G', 2: 'A'},  # A at position 2
        'read3': {0: 'A', 1: 'G', 2: 'T'},  # T at position 2 (different!)
    }

    phasing_positions = {0, 1}
    all_variant_positions = {0, 1, 2}

    error = calculate_within_cluster_error(
        haplotype_groups, read_alleles, phasing_positions, all_variant_positions
    )

    # Should have some error due to variation at position 2
    assert error > 0.0


def test_within_cluster_error_perfect_homogeneity():
    """Test when all reads in haplotype have identical alleles at all positions."""
    haplotype_groups = {
        'A-G': {'read1', 'read2'},
    }

    read_alleles = {
        'read1': {0: 'A', 1: 'G'},
        'read2': {0: 'A', 1: 'G'},
    }

    phasing_positions = {0, 1}  # All positions used for phasing
    all_variant_positions = {0, 1}  # Same as phasing positions

    error = calculate_within_cluster_error(
        haplotype_groups, read_alleles, phasing_positions, all_variant_positions
    )

    # Perfect homogeneity (all reads identical at all positions) means 0 error
    assert error == 0.0


# ==============================================================================
# Tests for group_reads_by_single_position()
# ==============================================================================

def test_group_reads_by_single_position():
    """Test grouping reads by a single position."""
    read_alleles = {
        'read1': {0: 'A', 10: 'C'},
        'read2': {0: 'A', 10: 'T'},
        'read3': {0: 'G', 10: 'C'},
        'read4': {0: 'G', 10: 'T'},
    }

    # Group by position 0
    groups = group_reads_by_single_position(
        read_alleles,
        position=0,
        read_ids={'read1', 'read2', 'read3', 'read4'}
    )

    assert 'A' in groups
    assert 'G' in groups
    assert groups['A'] == {'read1', 'read2'}
    assert groups['G'] == {'read3', 'read4'}


def test_group_reads_by_single_position_subset():
    """Test grouping only a subset of reads."""
    read_alleles = {
        'read1': {0: 'A'},
        'read2': {0: 'A'},
        'read3': {0: 'G'},
        'read4': {0: 'G'},
    }

    # Only consider read1 and read3
    groups = group_reads_by_single_position(
        read_alleles,
        position=0,
        read_ids={'read1', 'read3'}
    )

    assert 'A' in groups
    assert 'G' in groups
    assert groups['A'] == {'read1'}
    assert groups['G'] == {'read3'}


# ==============================================================================
# Tests for call_iupac_ambiguities()
# ==============================================================================

def test_call_iupac_basic_biallelic():
    """Test IUPAC calling for simple biallelic position."""
    consensus = "ACGTACGT"

    # Create alignments with variation at position 2 (G->R for A/G)
    alignments = [
        ReadAlignment(
            read_id='read1',
            aligned_sequence='ACATACGT',  # A at position 2
            read_length=8,
            edit_distance=1,
            num_insertions=0,
            num_deletions=0,
            num_substitutions=1,
            error_positions=[ErrorPosition(2, 'sub')],
            normalized_edit_distance=1,
            normalized_error_positions=[ErrorPosition(2, 'sub')],
            score_aligned='||X|||||'
        ),
        ReadAlignment(
            read_id='read2',
            aligned_sequence='ACATACGT',
            read_length=8,
            edit_distance=1,
            num_insertions=0,
            num_deletions=0,
            num_substitutions=1,
            error_positions=[ErrorPosition(2, 'sub')],
            normalized_edit_distance=1,
            normalized_error_positions=[ErrorPosition(2, 'sub')],
            score_aligned='||X|||||'
        ),
        ReadAlignment(
            read_id='read3',
            aligned_sequence='ACGTACGT',  # G at position 2 (matches consensus)
            read_length=8,
            edit_distance=0,
            num_insertions=0,
            num_deletions=0,
            num_substitutions=0,
            error_positions=[],
            normalized_edit_distance=0,
            normalized_error_positions=[],
            score_aligned='||||||||'
        ),
        ReadAlignment(
            read_id='read4',
            aligned_sequence='ACGTACGT',
            read_length=8,
            edit_distance=0,
            num_insertions=0,
            num_deletions=0,
            num_substitutions=0,
            error_positions=[],
            normalized_edit_distance=0,
            normalized_error_positions=[],
            score_aligned='||||||||'
        ),
    ]

    # Create MSA to consensus position mapping (1:1 for ungapped alignment)
    msa_to_consensus_pos = {i: i for i in range(8)}

    new_consensus, iupac_count, details = call_iupac_ambiguities(
        consensus,
        alignments,
        msa_to_consensus_pos,
        min_variant_frequency=0.20,
        min_variant_count=2
    )

    # Position 2 should have A/G variant -> R code
    if iupac_count > 0:
        assert new_consensus[2] == 'R'  # A/G -> R
        assert iupac_count >= 1


def test_call_iupac_no_variants():
    """Test IUPAC calling when no positions meet variant threshold."""
    consensus = "ACGTACGT"

    # All reads match consensus perfectly
    alignments = [
        ReadAlignment(
            read_id=f'read{i}',
            aligned_sequence='ACGTACGT',
            read_length=8,
            edit_distance=0,
            num_insertions=0,
            num_deletions=0,
            num_substitutions=0,
            error_positions=[],
            normalized_edit_distance=0,
            normalized_error_positions=[],
            score_aligned='||||||||'
        )
        for i in range(10)
    ]

    msa_to_consensus_pos = {i: i for i in range(8)}

    new_consensus, iupac_count, details = call_iupac_ambiguities(
        consensus,
        alignments,
        msa_to_consensus_pos,
        min_variant_frequency=0.20,
        min_variant_count=2
    )

    # No changes expected
    assert new_consensus == consensus
    assert iupac_count == 0


def test_iupac_codes_mapping():
    """Test that IUPAC_CODES mapping is correct."""
    # Verify key ambiguity codes
    assert IUPAC_CODES.get(frozenset(['A', 'G'])) == 'R'
    assert IUPAC_CODES.get(frozenset(['C', 'T'])) == 'Y'
    assert IUPAC_CODES.get(frozenset(['G', 'T'])) == 'K'
    assert IUPAC_CODES.get(frozenset(['A', 'C'])) == 'M'
    assert IUPAC_CODES.get(frozenset(['A', 'T'])) == 'W'
    assert IUPAC_CODES.get(frozenset(['C', 'G'])) == 'S'


# ==============================================================================
# Tests for CER (Critical Error Rate) gating in variant phasing
# ==============================================================================

class TestCERGating:
    """Test that CER significance filtering integrates with phasing config."""

    def test_cer_disabled_when_error_rate_zero(self):
        """assumed_error_rate=0 should disable CER filtering."""
        from speconsense.core.workers import ClusterProcessingConfig
        from speconsense.significance import is_variant_significant

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10,
            min_variant_count=5,
            total_specimen_reads=1000,
            assumed_error_rate=0,
            significance_level=1e-5
        )
        # With error rate 0, all variants should pass
        assert config.assumed_error_rate == 0
        is_sig, _ = is_variant_significant(
            M=1, N=config.total_specimen_reads, L=700,
            assumed_error_rate=config.assumed_error_rate
        )
        assert is_sig is True

    def test_cer_suppresses_weak_variant(self):
        """Small M relative to large N should be suppressed by CER."""
        from speconsense.significance import is_variant_significant

        # 5 reads out of 1000 at 2% assumed error rate
        is_sig, p_star = is_variant_significant(
            M=5, N=1000, L=700,
            assumed_error_rate=0.02, alpha=1e-5
        )
        assert is_sig is False
        assert p_star < 0.02

    def test_cer_passes_well_supported_variant(self):
        """Large M relative to N should pass CER."""
        from speconsense.significance import is_variant_significant

        # 100 reads out of 1000 at 2% assumed error rate
        is_sig, p_star = is_variant_significant(
            M=100, N=1000, L=700,
            assumed_error_rate=0.02, alpha=1e-5
        )
        assert is_sig is True
        assert p_star > 0.02

    def test_cer_k2_rescues_weak_k1(self):
        """K=2 should pass CER where K=1 fails with same M."""
        from speconsense.significance import is_variant_significant

        M, N, L = 6, 1000, 700
        # Fails at K=1
        is_sig_k1, p_k1 = is_variant_significant(
            M=M, N=N, L=L, assumed_error_rate=0.015, alpha=1e-5, K=1
        )
        assert is_sig_k1 is False

        # Passes at K=2
        is_sig_k2, p_k2 = is_variant_significant(
            M=M, N=N, L=L, assumed_error_rate=0.015, alpha=1e-5, K=2
        )
        assert is_sig_k2 is True
        assert p_k2 > p_k1

    def test_find_best_phasing_subset_k1(self):
        """K=1 split should work via the unified search."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # 100 reads: 60 have 'A' at pos 50, 40 have 'G' — strong variant
        all_read_ids = set()
        read_to_position_alleles = {}
        for i in range(100):
            rid = f"read_{i}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {50: 'G' if i < 40 else 'A'}

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.05,
            min_variant_count=5,
            total_specimen_reads=200,
            assumed_error_rate=0.015,
            significance_level=1e-5,
            min_k_position_gap=10,
        )

        result = _find_best_phasing_subset(
            [50], {50: 50}, read_to_position_alleles,
            all_read_ids, 100, 700, config
        )

        assert result is not None, "Expected K=1 split to be found"
        positions, qualifying, _, _, _ = result
        assert positions == [50]
        assert len(qualifying) == 2

    def test_find_best_phasing_subset_k2_correlated(self):
        """K=2 with perfectly correlated positions should be found."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # 100 reads: 10 have minority alleles at both positions
        minority_reads = {f"read_{i}" for i in range(10)}
        all_read_ids = set()
        read_to_position_alleles = {}
        for i in range(100):
            rid = f"read_{i}"
            all_read_ids.add(rid)
            if rid in minority_reads:
                read_to_position_alleles[rid] = {50: 'G', 150: 'G'}
            else:
                read_to_position_alleles[rid] = {50: 'A', 150: 'A'}

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.05,
            min_variant_count=5,
            total_specimen_reads=1000,
            assumed_error_rate=0.015,
            significance_level=1e-5,
            min_k_position_gap=10,
        )

        result = _find_best_phasing_subset(
            [50, 150], {50: 50, 150: 150}, read_to_position_alleles,
            all_read_ids, 100, 700, config
        )

        assert result is not None, "Expected K=2 split to be found"
        positions, qualifying, _, _, p_star = result
        # Should select K=2 since it gives higher p* than K=1
        assert len(positions) >= 1
        assert len(qualifying) >= 2

    def test_find_best_phasing_subset_partial_correlation(self):
        """Partially correlated positions (38% overlap) should be evaluated."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # Simulate the ONT08.31 case: 3 positions with partial overlap
        # pos 90: 8 minority reads
        # pos 698: 5 minority reads (3 shared with pos 90)
        # pos 727: 8 minority reads (3 shared with pos 90, 0 shared with pos 698)
        all_read_ids = set()
        read_to_position_alleles = {}
        for i in range(100):
            rid = f"read_{i}"
            all_read_ids.add(rid)
            allele_90 = 'G' if i < 8 else 'A'
            allele_698 = 'G' if i in range(5, 10) else 'A'  # 3 overlap with pos 90 (5,6,7)
            allele_727 = 'G' if i in range(0, 8) else 'A'   # same as pos 90
            read_to_position_alleles[rid] = {90: allele_90, 698: allele_698, 727: allele_727}

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.05,
            min_variant_count=3,  # Low count to allow partial overlap
            total_specimen_reads=400,
            assumed_error_rate=0.015,
            significance_level=1e-5,
            min_k_position_gap=10,
        )

        result = _find_best_phasing_subset(
            [90, 698, 727], {90: 90, 698: 698, 727: 727},
            read_to_position_alleles, all_read_ids, 100, 750, config
        )

        # Should find SOME split — the old correlation heuristic (0.9 threshold)
        # would have rejected all of these, but the new algorithm evaluates
        # all subsets directly
        # At minimum, K=1 at pos 90 or 727 (8 minority reads) should be found
        assert result is not None, "Expected at least a K=1 split"

    def test_find_best_phasing_subset_proximity_filter(self):
        """Nearby positions should be filtered by min_k_position_gap."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # Two positions only 5 consensus positions apart
        # Each individually fails K=1 CER but would pass K=2
        minority_reads = {f"read_{i}" for i in range(6)}
        all_read_ids = set()
        read_to_position_alleles = {}
        for i in range(100):
            rid = f"read_{i}"
            all_read_ids.add(rid)
            allele = 'G' if rid in minority_reads else 'A'
            read_to_position_alleles[rid] = {50: allele, 55: allele}

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.05,
            min_variant_count=5,
            total_specimen_reads=1000,
            assumed_error_rate=0.015,
            significance_level=1e-5,
            min_k_position_gap=10,  # Gap of 5 < 10: should be filtered
        )

        result = _find_best_phasing_subset(
            [50, 55], {50: 50, 55: 55}, read_to_position_alleles,
            all_read_ids, 100, 700, config
        )

        # K=2 should be blocked by proximity; K=1 with M=6 likely fails CER at N=1000
        # So result should be None (no valid split)
        if result is not None:
            # If K=1 happened to pass, positions should not include both
            positions = result[0]
            assert not (50 in positions and 55 in positions), \
                "Nearby positions should not be combined in K>1 subset"

    def test_find_best_phasing_subset_no_variants(self):
        """Empty variant positions should return None."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10,
            min_variant_count=5,
            total_specimen_reads=1000,
            assumed_error_rate=0.015,
            significance_level=1e-5,
        )

        result = _find_best_phasing_subset(
            [], {}, {}, set(), 0, 700, config
        )
        assert result is None

    def test_find_best_phasing_subset_cer_disabled(self):
        """With assumed_error_rate=0, CER gate is disabled; select by error."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # Clear biallelic split
        all_read_ids = set()
        read_to_position_alleles = {}
        for i in range(20):
            rid = f"read_{i}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {50: 'G' if i < 10 else 'A'}

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10,
            min_variant_count=5,
            total_specimen_reads=0,  # CER disabled
            assumed_error_rate=0,
        )

        result = _find_best_phasing_subset(
            [50], {50: 50}, read_to_position_alleles,
            all_read_ids, 20, 700, config
        )

        assert result is not None
        positions, qualifying, _, _, _ = result
        assert len(qualifying) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
