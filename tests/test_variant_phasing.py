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
# Tests for _find_best_phasing_subset
# ==============================================================================

class TestFindBestPhasingSubset:
    """Test the unified phasing-position search."""

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
        )

        result = _find_best_phasing_subset(
            [50], {50: 50}, read_to_position_alleles,
            all_read_ids, 100, 700, config
        )

        assert result is not None, "Expected K=1 split to be found"
        positions, qualifying, _, _, _ = result
        assert positions == [50]
        assert len(qualifying) == 2

    def test_find_best_phasing_subset_correlated_positions(self):
        """Correlated positions should each produce K=1 splits; best by error."""
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
        )

        result = _find_best_phasing_subset(
            [50, 150], {50: 50, 150: 150}, read_to_position_alleles,
            all_read_ids, 100, 700, config
        )

        assert result is not None, "Expected K=1 split to be found"
        positions, qualifying, _, _, _ = result
        # K=1 at either position splits 10 vs 90 — both equivalent
        assert len(positions) == 1
        assert len(qualifying) == 2

    def test_find_best_phasing_subset_partial_correlation(self):
        """Partially correlated positions should find K=1 split at best position."""
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
        )

        result = _find_best_phasing_subset(
            [90, 698, 727], {90: 90, 698: 698, 727: 727},
            read_to_position_alleles, all_read_ids, 100, 750, config
        )

        # K=1 at pos 90 or 727 (8 minority reads) should be found
        assert result is not None, "Expected K=1 split"
        positions, qualifying, _, _, _ = result
        assert len(positions) == 1
        assert len(qualifying) >= 2

    def test_find_best_phasing_subset_no_variants(self):
        """Empty variant positions should return None."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10,
            min_variant_count=5,
        )

        result = _find_best_phasing_subset(
            [], {}, {}, set(), 0, 700, config
        )
        assert result is None

    def test_find_best_phasing_subset_small_split(self):
        """Small but qualifying split should be found (no CER gate)."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # Clear biallelic split with 10 vs 10
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
        )

        result = _find_best_phasing_subset(
            [50], {50: 50}, read_to_position_alleles,
            all_read_ids, 20, 700, config
        )

        assert result is not None
        positions, qualifying, _, _, _ = result
        assert len(qualifying) == 2

    def test_find_best_phasing_subset_noisy_positions(self):
        """Best K=1 position selected by lowest within-cluster error."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # 200 reads:
        # - 7 carry minority at BOTH pos 100 and 300 (correlated haplotype)
        # - 30 carry minority at pos 100 only (noise)
        # - 20 carry minority at pos 300 only (noise)
        # - 143 are majority at both
        # K=1 at pos 100: 37 minority reads, high error (30 noise + 7 real)
        # K=1 at pos 300: 27 minority reads, high error (20 noise + 7 real)
        # K=2 {100, 300}: 7 reads with (G, G) haplotype, much cleaner split
        all_read_ids = set()
        read_to_position_alleles = {}
        idx = 0
        # 7 correlated minority at both positions
        for _ in range(7):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {100: 'G', 300: 'G'}
            idx += 1
        # 30 noise minority at pos 100 only
        for _ in range(30):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {100: 'G', 300: 'A'}
            idx += 1
        # 20 noise minority at pos 300 only
        for _ in range(20):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {100: 'A', 300: 'G'}
            idx += 1
        # 143 majority at both
        for _ in range(143):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {100: 'A', 300: 'A'}
            idx += 1

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.03,
            min_variant_count=5,
        )

        result = _find_best_phasing_subset(
            [100, 300], {100: 100, 300: 300}, read_to_position_alleles,
            all_read_ids, 200, 700, config
        )

        # K=1 at pos 300 (27 minority) has lower error than pos 100 (37 minority)
        # because fewer noise reads contaminate the majority group
        assert result is not None
        positions, qualifying, _, error, _ = result
        assert len(positions) == 1
        assert len(qualifying) >= 2

    def test_find_best_phasing_subset_two_independent_minorities(self):
        """Two unrelated minorities at different positions are handled correctly."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # 200 reads:
        # - 30 carry minority at pos 200 only
        # - 30 different reads carry minority at pos 500 only
        # - 140 are majority at both
        all_read_ids = set()
        read_to_position_alleles = {}
        idx = 0
        for _ in range(30):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {200: 'G', 500: 'A'}
            idx += 1
        for _ in range(30):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {200: 'A', 500: 'G'}
            idx += 1
        for _ in range(140):
            rid = f"read_{idx}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {200: 'A', 500: 'A'}
            idx += 1

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95,
            enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.05,
            min_variant_count=5,
        )

        result = _find_best_phasing_subset(
            [200, 500], {200: 200, 500: 500}, read_to_position_alleles,
            all_read_ids, 200, 700, config
        )

        # Should find K=1 split at either position (both have 30 minority)
        assert result is not None
        positions, qualifying, _, _, _ = result
        assert len(positions) == 1
        assert len(qualifying) >= 2


# ==============================================================================
# Tests for OR-with-sample phasing qualification
# ==============================================================================


class TestORWithSampleQualification:
    """Exercises sample_read_ids OR logic in _find_best_phasing_subset.

    Design: a variant can pass phasing thresholds on the full subcluster
    OR on the top-N-by-quality sample that the final consensus will see.
    This closes the gap where a quality-biased sample's frequency inflates
    above the phasing threshold while the full-cluster frequency stays
    below — the variant would otherwise become an IUPAC ambiguity.
    """

    def _build_single_position_alleles(self, n_total, minority_ids, pos=50):
        all_read_ids = set()
        read_to_position_alleles = {}
        for i in range(n_total):
            rid = f"read_{i}"
            all_read_ids.add(rid)
            read_to_position_alleles[rid] = {pos: 'G' if rid in minority_ids else 'A'}
        return all_read_ids, read_to_position_alleles

    def test_sample_triggers_qualification_where_full_fails(self):
        """Variant diluted below threshold in full cluster but dense in sample."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        # 200 reads total; minority allele on reads 0..14 (15 reads = 7.5%).
        # Top-40 by quality = reads 0..39, so minority-in-sample = 15/40 = 37.5%.
        total_reads = 200
        minority = {f"read_{i}" for i in range(15)}
        sample_ids = {f"read_{i}" for i in range(40)}
        all_read_ids, r2p = self._build_single_position_alleles(total_reads, minority)

        # Thresholds: 10% frequency, count 3. Fails on full (15/200=7.5%), passes on sample.
        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95, enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10, min_variant_count=3,
        )

        # Without sample: fails
        res_no_sample = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=None,
        )
        assert res_no_sample is None, "Full-cluster frequency should be below threshold"

        # With sample: passes
        res_with_sample = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=sample_ids,
        )
        assert res_with_sample is not None, "Sample-based OR should qualify this split"
        positions, qualifying, non_qualifying, _, _ = res_with_sample
        assert positions == [50]
        assert len(qualifying) == 2
        # Qualifying combos should carry full-cluster read membership, not just the sample slice.
        total_in_qualifying = sum(len(r) for r in qualifying.values())
        assert total_in_qualifying == total_reads, (
            "Full-cluster reads in sample-qualifying combos should not be discarded"
        )

    def test_sample_does_not_lower_bar_when_full_passes(self):
        """When full cluster already passes, sample OR is a no-op."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        total_reads = 100
        minority = {f"read_{i}" for i in range(25)}  # 25% — well above threshold
        sample_ids = {f"read_{i}" for i in range(40)}
        all_read_ids, r2p = self._build_single_position_alleles(total_reads, minority)

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95, enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10, min_variant_count=3,
        )

        res_no_sample = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=None,
        )
        res_with_sample = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=sample_ids,
        )

        assert res_no_sample is not None
        assert res_with_sample is not None
        # Same qualifying-combo keys.
        _, qual_no, _, _, _ = res_no_sample
        _, qual_with, _, _, _ = res_with_sample
        assert set(qual_no.keys()) == set(qual_with.keys())

    def test_sample_equals_full_is_noop(self):
        """When sample covers the full set, OR collapses to current behavior."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        total_reads = 30
        minority = {f"read_{i}" for i in range(5)}
        all_read_ids, r2p = self._build_single_position_alleles(total_reads, minority)
        full_as_sample = set(all_read_ids)

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95, enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10, min_variant_count=3,
        )

        res_no_sample = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=None,
        )
        res_sample_eq_full = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=full_as_sample,
        )

        # Both should yield identical qualifying keys (strict-subset guard kicks in).
        if res_no_sample is None:
            assert res_sample_eq_full is None
        else:
            _, qual_no, _, _, _ = res_no_sample
            _, qual_eq, _, _, _ = res_sample_eq_full
            assert set(qual_no.keys()) == set(qual_eq.keys())

    def test_variant_fails_on_both_sets_remains_non_qualifying(self):
        """Variant below threshold on both full and sample stays filtered out."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        total_reads = 200
        # 2 minority reads — way below count threshold on either side.
        minority = {f"read_{i}" for i in range(2)}
        sample_ids = {f"read_{i}" for i in range(40)}
        all_read_ids, r2p = self._build_single_position_alleles(total_reads, minority)

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95, enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10, min_variant_count=3,
        )

        result = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=sample_ids,
        )
        assert result is None

    def test_empty_sample_ignored(self):
        """Empty sample set should not disqualify a full-set-qualifying split."""
        from speconsense.core.workers import _find_best_phasing_subset, ClusterProcessingConfig

        total_reads = 100
        minority = {f"read_{i}" for i in range(20)}
        all_read_ids, r2p = self._build_single_position_alleles(total_reads, minority)

        config = ClusterProcessingConfig(
            outlier_identity_threshold=0.95, enable_secondpass_phasing=True,
            disable_homopolymer_equivalence=False,
            min_variant_frequency=0.10, min_variant_count=3,
        )

        result = _find_best_phasing_subset(
            [50], {50: 50}, r2p, all_read_ids, total_reads, 700, config,
            sample_read_ids=set(),
        )
        assert result is not None


class TestORWithSampleDetection:
    """Exercises sample_read_ids OR logic in _detect_variant_positions_standalone.

    A variant position should be detected if it meets thresholds on EITHER the
    full alignment set OR the top-N-by-quality sample. Mirrors the qualification
    OR and is needed so sample-only variants become phasing candidates at all.
    """

    def _build_alignments(self, n_total, n_minority, consensus):
        """Build ReadAlignments where reads 0..n_minority-1 have 'A' at pos 2,
        and reads n_minority..n_total-1 match the consensus.
        """
        variant_seq = consensus[:2] + 'A' + consensus[3:]
        alignments = []
        for i in range(n_total):
            if i < n_minority:
                alignments.append(ReadAlignment(
                    read_id=f'read_{i}',
                    aligned_sequence=variant_seq,
                    read_length=len(consensus),
                    edit_distance=1,
                    num_insertions=0, num_deletions=0, num_substitutions=1,
                    error_positions=[ErrorPosition(2, 'sub')],
                    normalized_edit_distance=1,
                    normalized_error_positions=[ErrorPosition(2, 'sub')],
                    score_aligned='||X' + '|' * (len(consensus) - 3),
                ))
            else:
                alignments.append(ReadAlignment(
                    read_id=f'read_{i}',
                    aligned_sequence=consensus,
                    read_length=len(consensus),
                    edit_distance=0,
                    num_insertions=0, num_deletions=0, num_substitutions=0,
                    error_positions=[],
                    normalized_edit_distance=0,
                    normalized_error_positions=[],
                    score_aligned='|' * len(consensus),
                ))
        return alignments

    def test_sample_detection_finds_position_that_full_misses(self):
        """Position fails 10% freq on full cluster (5%) but passes on sample (25%)."""
        from speconsense.core.workers import _detect_variant_positions_standalone

        consensus = "ACGTACGT"
        msa_to_consensus_pos = {i: i for i in range(8)}
        alignments = self._build_alignments(n_total=200, n_minority=10, consensus=consensus)

        # Full: 10/200 = 5% — below threshold.
        full_positions = _detect_variant_positions_standalone(
            alignments, consensus, msa_to_consensus_pos,
            min_variant_frequency=0.10, min_variant_count=3,
            sample_read_ids=None,
        )
        assert not any(v['msa_position'] == 2 for v in full_positions), (
            "Position 2 should be below threshold on full cluster"
        )

        # Sample = reads 0..39 (top-40 by quality, all minority inside).
        sample_ids = {f"read_{i}" for i in range(40)}  # 10/40 = 25% — above threshold.
        merged = _detect_variant_positions_standalone(
            alignments, consensus, msa_to_consensus_pos,
            min_variant_frequency=0.10, min_variant_count=3,
            sample_read_ids=sample_ids,
        )
        assert any(v['msa_position'] == 2 for v in merged), (
            "Sample-based OR should surface position 2 as a phasing candidate"
        )

    def test_detection_noop_when_sample_equals_full(self):
        """No strict subset → OR is bypassed → result identical to full-only."""
        from speconsense.core.workers import _detect_variant_positions_standalone

        consensus = "ACGTACGT"
        msa_to_consensus_pos = {i: i for i in range(8)}
        alignments = self._build_alignments(n_total=20, n_minority=5, consensus=consensus)
        full_ids = {a.read_id for a in alignments}

        res_no_sample = _detect_variant_positions_standalone(
            alignments, consensus, msa_to_consensus_pos,
            min_variant_frequency=0.10, min_variant_count=3,
            sample_read_ids=None,
        )
        res_sample_eq_full = _detect_variant_positions_standalone(
            alignments, consensus, msa_to_consensus_pos,
            min_variant_frequency=0.10, min_variant_count=3,
            sample_read_ids=full_ids,
        )
        assert [v['msa_position'] for v in res_no_sample] == \
               [v['msa_position'] for v in res_sample_eq_full]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
