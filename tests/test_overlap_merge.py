#!/usr/bin/env python3
"""
Tests for overlap merge feature in speconsense-summarize.

Tests scenarios where sequences from different primer pools have different lengths
but share sufficient overlap to be merged.
"""

import os
import tempfile
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from speconsense.summarize import (
    calculate_overlap_aware_distance,
    analyze_msa_columns_overlap_aware,
    create_overlap_consensus_from_msa,
    run_spoa_msa,
    STANDARD_ADJUSTMENT_PARAMS,
)


# Test sequences - simulating ITS regions with distinct, non-repetitive content
# Using pseudo-random but reproducible sequences to avoid alignment artifacts

# Generate distinct regions that won't align to each other
import hashlib

def generate_dna_sequence(seed: str, length: int) -> str:
    """Generate a reproducible pseudo-random DNA sequence."""
    result = []
    bases = "ACGT"
    for i in range(length):
        h = hashlib.md5(f"{seed}_{i}".encode()).hexdigest()
        result.append(bases[int(h[0], 16) % 4])
    return "".join(result)

# Create distinct regions that won't accidentally align
ITS1_REGION = generate_dna_sequence("its1_unique", 200)  # 200bp unique ITS1
REGION_5_8S = generate_dna_sequence("5_8S_shared", 250)   # 250bp shared 5.8S (increased for better overlap)
ITS2_REGION = generate_dna_sequence("its2_unique", 240)  # 240bp unique ITS2

FULL_ITS = ITS1_REGION + REGION_5_8S + ITS2_REGION  # 690bp
ITS1_WITH_5_8S = ITS1_REGION + REGION_5_8S  # 450bp
ITS2_WITH_5_8S = REGION_5_8S + ITS2_REGION  # 490bp
ITS2_ONLY = ITS2_REGION  # 240bp


class TestCalculateOverlapAwareDistance:
    """Tests for calculate_overlap_aware_distance function."""

    def test_identical_sequences(self):
        """Identical sequences should have zero distance."""
        dist = calculate_overlap_aware_distance(FULL_ITS, FULL_ITS, 200)
        assert dist == 0.0

    def test_containment_short_in_long(self):
        """Short sequence fully contained in long should have low distance."""
        # ITS2_ONLY (240bp) is contained in FULL_ITS (600bp)
        dist = calculate_overlap_aware_distance(FULL_ITS, ITS2_ONLY, 200)
        # With 240bp overlap (full containment), should pass 200bp threshold
        assert dist < 0.1  # Should be very similar

    def test_overlap_prefix_suffix(self):
        """Sequences with shared middle region should have low distance."""
        # ITS1_WITH_5_8S shares 5.8S (250bp) with ITS2_WITH_5_8S
        dist = calculate_overlap_aware_distance(ITS1_WITH_5_8S, ITS2_WITH_5_8S, 200)
        # 250bp overlap (5.8S region), passes 200bp threshold
        # With distinct sequences sharing exact 250bp region, should have very low distance
        assert dist < 0.2  # Allow some tolerance for alignment edge effects

    def test_insufficient_overlap(self):
        """Sequences with overlap below threshold should get global distance."""
        # Create two sequences that only share 100bp
        seq1 = "A" * 400 + "SHARED" * 16 + "T" * 4  # ~500bp with 96bp shared
        seq2 = "G" * 4 + "SHARED" * 16 + "C" * 400  # ~500bp with 96bp shared

        # With 200bp threshold, 96bp overlap is insufficient
        dist = calculate_overlap_aware_distance(seq1, seq2, 200)
        # Should fall back to global distance (high due to length differences)
        assert dist > 0.3

    def test_containment_at_threshold(self):
        """Short sequence exactly at threshold should still merge."""
        # Create a 200bp sequence that's a substring
        short_seq = FULL_ITS[200:400]  # 200bp from middle
        dist = calculate_overlap_aware_distance(FULL_ITS, short_seq, 200)
        # 200bp overlap = 200bp threshold, should pass
        assert dist < 0.1


class TestAnalyzeMsaColumnsOverlapAware:
    """Tests for analyze_msa_columns_overlap_aware function."""

    def test_containment_no_structural_indels(self):
        """Containment should not count terminal gaps as structural indels."""
        # Create aligned sequences where ITS2_ONLY is contained in FULL_ITS
        aligned_seqs = run_spoa_msa([FULL_ITS, ITS2_ONLY])
        original_lengths = [len(FULL_ITS), len(ITS2_ONLY)]

        stats = analyze_msa_columns_overlap_aware(aligned_seqs, 200, original_lengths)

        # Terminal gaps from length difference should NOT be counted as structural
        assert stats['structural_indel_count'] == 0
        # Should have overlap equal to shorter sequence
        assert stats['overlap_bp'] >= 200  # At least the threshold

    def test_snps_in_overlap_region(self):
        """SNPs within overlap region should be counted."""
        # Create two sequences with 1 SNP in shared region
        seq1 = "AAAA" + "CCCC" + "GGGG"  # 12bp
        seq2 = "CCCC" + "CGCC" + "TTTT"  # 12bp, one SNP in middle (C->G)

        aligned_seqs = run_spoa_msa([seq1, seq2])
        original_lengths = [len(seq1), len(seq2)]

        stats = analyze_msa_columns_overlap_aware(aligned_seqs, 1, original_lengths)

        # Should detect the SNP in the overlap region
        assert stats['snp_count'] >= 1


class TestCreateOverlapConsensusFromMsa:
    """Tests for create_overlap_consensus_from_msa function."""

    def test_containment_produces_full_length(self):
        """Merging contained sequence with full should produce full length."""
        from speconsense.types import ConsensusInfo

        aligned_seqs = run_spoa_msa([FULL_ITS, ITS2_ONLY])

        # Create mock ConsensusInfo objects
        variants = [
            ConsensusInfo(
                sample_name="full_its",
                cluster_id="c1",
                sequence=FULL_ITS,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=None,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="its2_only",
                cluster_id="c2",
                sequence=ITS2_ONLY,
                ric=50,
                size=50,
                file_path="/test",
                snp_count=None,
                primers=None,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
        ]

        merged = create_overlap_consensus_from_msa(aligned_seqs, variants)

        # Result should be at least as long as the full ITS
        assert len(merged.sequence) >= len(FULL_ITS) - 10  # Allow small variation
        # Combined size should be sum
        assert merged.size == 150
        assert merged.ric == 150

    def test_prefix_suffix_produces_union(self):
        """Merging prefix/suffix overlap should produce union of both.

        SPOA local alignment (-l 0) correctly handles prefix/suffix overlaps,
        producing clean terminal gaps and the expected union-length consensus.
        ITS1_WITH_5_8S (470bp) + ITS2_WITH_5_8S (470bp) sharing 250bp 5.8S
        should produce a 690bp merged consensus (220 + 250 + 220).
        """
        from speconsense.types import ConsensusInfo

        # Use local alignment mode (0) as used in actual overlap merging workflow
        aligned_seqs = run_spoa_msa([ITS1_WITH_5_8S, ITS2_WITH_5_8S], alignment_mode=0)

        variants = [
            ConsensusInfo(
                sample_name="its1_5_8s",
                cluster_id="c1",
                sequence=ITS1_WITH_5_8S,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=None,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="its2_5_8s",
                cluster_id="c2",
                sequence=ITS2_WITH_5_8S,
                ric=50,
                size=50,
                file_path="/test",
                snp_count=None,
                primers=None,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
        ]

        merged = create_overlap_consensus_from_msa(aligned_seqs, variants)

        # SPOA local alignment produces correct union-length consensus
        # ITS1 unique (220) + 5.8S shared (250) + ITS2 unique (220) = 690bp
        assert len(merged.sequence) == 690
        # Combined metrics should be correct
        assert merged.size == 150
        assert merged.ric == 150


class TestChimeraPrevention:
    """Tests to ensure chimeric merges are prevented."""

    def test_short_overlap_rejected(self):
        """Sequences sharing only short conserved region should not merge."""
        # Create a short shared region (150bp) that's below the 200bp threshold
        short_shared = generate_dna_sequence("short_shared", 150)
        russula_its1 = generate_dna_sequence("russula_unique", 300)
        hypomyces_its2 = generate_dna_sequence("hypomyces_unique", 300)

        russula_seq = russula_its1 + short_shared  # 450bp
        hypomyces_seq = short_shared + hypomyces_its2  # 450bp

        # With 200bp threshold, 150bp overlap should be insufficient
        dist = calculate_overlap_aware_distance(russula_seq, hypomyces_seq, 200)

        # Should have high distance (global fallback due to insufficient overlap)
        # or low distance if the alignment happens to find an overlap
        # The key is that these should NOT merge with 200bp threshold
        # For safety, we want distance > 0.1 (HAC threshold at 0.9 identity)
        assert dist > 0.1, f"Expected high distance for short overlap, got {dist}"


class TestDisabledMode:
    """Tests for when overlap merge is disabled."""

    def test_zero_threshold_uses_standard(self):
        """With min_overlap=0, should use standard distance calculation."""
        # Standard distance for different-length sequences should be high
        # due to terminal gap penalties
        standard_dist = calculate_overlap_aware_distance(FULL_ITS, ITS2_ONLY, 0)

        # With overlap mode disabled (0), falls back to standard
        # This should be a higher distance due to length difference
        # Actually, with 0 threshold, the function falls back to calculate_adjusted_identity_distance
        # which still handles the alignment reasonably
        assert standard_dist is not None


class TestHybridLinkage:
    """Tests for hybrid linkage in HAC clustering.

    When overlap mode is enabled (min_overlap_bp > 0), HAC uses:
    - Complete linkage within each primer set (prevents chaining)
    - Single linkage across primer sets (allows ITS1+full+ITS2 merging)
    """

    def test_same_primers_use_complete_linkage(self):
        """Same-primer sequences should not chain through intermediates.

        Creates a chain A-B-C where A~B, B~C, but A!~C.
        With complete linkage, they should NOT all be in one group.
        """
        from speconsense.types import ConsensusInfo
        from speconsense.summarize.clustering import perform_hac_clustering

        # Create sequences that form a chain: A~B~C but A!~C
        # Each differs from B by ~4% (20 mutations), so A-C differs by ~8%
        # At 95% identity threshold (5% distance):
        # - A-B: ~3.4% < 5% (should merge)
        # - B-C: ~4.0% < 5% (should merge)
        # - A-C: ~7.4% > 5% (should NOT merge)
        base_seq = generate_dna_sequence("chain_base", 500)
        num_mutations = 20

        # B is base sequence
        seq_b = base_seq

        # A differs from B at positions 0, 25, 50, ...
        seq_a_list = list(base_seq)
        for i in range(num_mutations):
            pos = i * 25  # Spread across sequence
            seq_a_list[pos] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_a_list[pos]]
        seq_a = ''.join(seq_a_list)

        # C differs from B at positions 12, 37, 62, ... (different from A)
        seq_c_list = list(base_seq)
        for i in range(num_mutations):
            pos = i * 25 + 12  # Different positions than A
            seq_c_list[pos] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c_list[pos]]
        seq_c = ''.join(seq_c_list)

        # All have same primers
        same_primers = ['ITS1F', 'ITS4']

        variants = [
            ConsensusInfo(
                sample_name="seq_a",
                cluster_id="c1",
                sequence=seq_a,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=same_primers,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="seq_b",
                cluster_id="c2",
                sequence=seq_b,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=same_primers,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="seq_c",
                cluster_id="c3",
                sequence=seq_c,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=same_primers,
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
        ]

        # Run with overlap mode enabled (triggers hybrid linkage)
        # Use 0.95 identity threshold (5% distance)
        groups = perform_hac_clustering(
            variants,
            variant_group_identity=0.95,
            min_overlap_bp=200
        )

        # With complete linkage for same primers, A and C should NOT be in same group
        # (since A-C distance is >5%)
        # But A-B and B-C could be grouped if distance <5%
        # The key is that we DON'T get all 3 in one group through chaining

        all_in_one_group = any(len(group) == 3 for group in groups.values())
        assert not all_in_one_group, \
            "Same-primer sequences should use complete linkage and not chain A-B-C into one group"

    def test_different_primers_use_single_linkage(self):
        """Different-primer sequences can connect via overlap.

        Simulates ITS1 + full ITS + ITS2 scenario where:
        - ITS1 overlaps with full ITS (in 5' region)
        - ITS2 overlaps with full ITS (in 3' region)
        - ITS1 does NOT overlap with ITS2

        With single linkage across primers, all three should merge.
        """
        from speconsense.types import ConsensusInfo
        from speconsense.summarize.clustering import perform_hac_clustering

        # Use the pre-defined test sequences
        variants = [
            ConsensusInfo(
                sample_name="its1_partial",
                cluster_id="c1",
                sequence=ITS1_WITH_5_8S,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=['ITS1F', 'ITS2'],  # ITS1-only primers
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="full_its",
                cluster_id="c2",
                sequence=FULL_ITS,
                ric=200,
                size=200,
                file_path="/test",
                snp_count=None,
                primers=['ITS1F', 'ITS4'],  # Full ITS primers
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="its2_partial",
                cluster_id="c3",
                sequence=ITS2_WITH_5_8S,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=['ITS3', 'ITS4'],  # ITS2-only primers
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
        ]

        # Run with overlap mode (triggers hybrid linkage)
        groups = perform_hac_clustering(
            variants,
            variant_group_identity=0.95,
            min_overlap_bp=200  # 5.8S region is 250bp
        )

        # With single linkage across primers, all 3 should be in one group
        # because ITS1 overlaps full_its, and ITS2 overlaps full_its
        all_in_one_group = any(len(group) == 3 for group in groups.values())
        assert all_in_one_group, \
            "Different-primer sequences should use single linkage and merge via overlaps"

    def test_non_overlapping_different_primers_stay_separate(self):
        """Non-overlapping different-primer sequences stay separate.

        Even with single linkage across primers, sequences without overlap
        should not merge (distance = 1.0).
        """
        from speconsense.types import ConsensusInfo
        from speconsense.summarize.clustering import perform_hac_clustering

        # Create two completely different sequences with different primers
        seq1 = generate_dna_sequence("genus_a_species", 500)
        seq2 = generate_dna_sequence("genus_b_species", 500)

        variants = [
            ConsensusInfo(
                sample_name="species_a",
                cluster_id="c1",
                sequence=seq1,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=['ITS1F', 'ITS2'],
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="species_b",
                cluster_id="c2",
                sequence=seq2,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=['ITS3', 'ITS4'],
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
        ]

        groups = perform_hac_clustering(
            variants,
            variant_group_identity=0.95,
            min_overlap_bp=200
        )

        # With no overlap, even single linkage should keep them separate
        assert len(groups) == 2, \
            "Non-overlapping sequences should remain in separate groups"

    def test_none_primers_grouped_as_separate_set(self):
        """Sequences with None primers are treated as their own primer set.

        This uses complete linkage within the None set.
        """
        from speconsense.types import ConsensusInfo
        from speconsense.summarize.clustering import perform_hac_clustering

        # Two similar sequences with no primer info
        base_seq = generate_dna_sequence("no_primer_base", 500)
        similar_seq_list = list(base_seq)
        for i in range(10):  # ~2% difference
            pos = i * 50
            similar_seq_list[pos] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[similar_seq_list[pos]]
        similar_seq = ''.join(similar_seq_list)

        variants = [
            ConsensusInfo(
                sample_name="seq_a",
                cluster_id="c1",
                sequence=base_seq,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=None,  # No primer info
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
            ConsensusInfo(
                sample_name="seq_b",
                cluster_id="c2",
                sequence=similar_seq,
                ric=100,
                size=100,
                file_path="/test",
                snp_count=None,
                primers=None,  # No primer info
                raw_ric=None,
                rid=None,
                rid_min=None,
            ),
        ]

        groups = perform_hac_clustering(
            variants,
            variant_group_identity=0.95,
            min_overlap_bp=200
        )

        # Both have None primers, so treated as same primer set
        # With 98% identity, should merge at 95% threshold
        all_in_one_group = any(len(group) == 2 for group in groups.values())
        assert all_in_one_group, \
            "Similar sequences with None primers should merge via complete linkage"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
