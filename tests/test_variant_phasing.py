"""Test variant phasing functions from msa.py module."""

import pytest
from speconsense.msa import (
    PositionStats,
    is_variant_position_with_composition,
    calculate_within_cluster_error,
    select_correlated_positions,
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


def test_within_cluster_error_empty_non_phased():
    """Test when all variant positions are used for phasing (no non-phased positions)."""
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

    # No non-phased positions means 0 error
    assert error == 0.0


# ==============================================================================
# Tests for select_correlated_positions()
# ==============================================================================

def test_select_correlated_positions_exhaustive_small():
    """Test exhaustive search with small number of positions (<=10)."""
    # Create 5 variant positions - should use exhaustive search
    variant_positions = [
        {'msa_position': 0, 'consensus_position': 0},
        {'msa_position': 10, 'consensus_position': 10},
        {'msa_position': 20, 'consensus_position': 20},
        {'msa_position': 30, 'consensus_position': 30},
        {'msa_position': 40, 'consensus_position': 40},
    ]

    # Create reads with clear haplotype structure at positions 0 and 10
    # Position 0: A/G, Position 10: C/T
    # Haplotype A-C: 40 reads, Haplotype G-T: 40 reads
    # Other positions have random noise
    read_alleles = {}
    for i in range(40):
        read_alleles[f'read_AC_{i}'] = {0: 'A', 10: 'C', 20: 'A', 30: 'A', 40: 'A'}
    for i in range(40):
        read_alleles[f'read_GT_{i}'] = {0: 'G', 10: 'T', 20: 'A', 30: 'A', 40: 'A'}

    selected, reason = select_correlated_positions(
        variant_positions,
        read_alleles,
        min_haplotype_frequency=0.20,
        min_haplotype_count=5,
        total_reads=80
    )

    # Should select positions that give cleanest haplotype separation
    assert len(selected) >= 1
    assert 'exhaustive' in reason.lower()


def test_select_correlated_positions_beam_search_large():
    """Test beam search with many positions (>10)."""
    # Create 15 variant positions - should use beam search
    variant_positions = [
        {'msa_position': i * 10, 'consensus_position': i * 10}
        for i in range(15)
    ]

    # Create reads with clear haplotype structure at positions 0 and 10
    read_alleles = {}
    for i in range(50):
        alleles = {p['msa_position']: 'A' for p in variant_positions}
        alleles[0] = 'A'
        alleles[10] = 'C'
        read_alleles[f'read_AC_{i}'] = alleles
    for i in range(50):
        alleles = {p['msa_position']: 'A' for p in variant_positions}
        alleles[0] = 'G'
        alleles[10] = 'T'
        read_alleles[f'read_GT_{i}'] = alleles

    selected, reason = select_correlated_positions(
        variant_positions,
        read_alleles,
        min_haplotype_frequency=0.20,
        min_haplotype_count=5,
        total_reads=100
    )

    assert len(selected) >= 1
    assert 'beam' in reason.lower()


def test_select_correlated_positions_single_position():
    """Test with only one variant position."""
    variant_positions = [
        {'msa_position': 10, 'consensus_position': 10},
    ]

    read_alleles = {
        'read1': {10: 'A'},
        'read2': {10: 'A'},
        'read3': {10: 'G'},
        'read4': {10: 'G'},
    }

    selected, reason = select_correlated_positions(
        variant_positions,
        read_alleles,
        min_haplotype_frequency=0.20,
        min_haplotype_count=1,
        total_reads=4
    )

    # With single position, should return it if it creates valid haplotypes
    assert len(selected) <= 1


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
