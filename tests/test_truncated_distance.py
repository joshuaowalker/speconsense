"""Tests for truncated sequence handling in adjusted identity distance.

Verifies that length-mismatched sequences (e.g. 220bp vs 660bp) don't produce
artificially low distances that distort HAC grouping.
"""

import random
import unittest

from speconsense.summarize.iupac import calculate_adjusted_identity_distance


def generate_dna_sequence(seed_str: str, length: int) -> str:
    """Generate a deterministic DNA sequence from a seed string."""
    rng = random.Random(seed_str)
    return ''.join(rng.choice('ACGT') for _ in range(length))


class TestTruncatedDistance(unittest.TestCase):
    """Unit tests for calculate_adjusted_identity_distance with length mismatches."""

    def test_truncated_prefix_gets_high_distance(self):
        """A 220bp prefix vs a 660bp sequence should have high distance (~0.67)."""
        full_seq = generate_dna_sequence("full", 660)
        truncated = full_seq[:220]

        distance = calculate_adjusted_identity_distance(truncated, full_seq)

        # Previously returned ~0.0 due to terminal gap exclusion.
        # With the fix, falls back to raw edlib identity.
        # Raw identity: 220 matching / 660 total → ~0.33 identity → ~0.67 distance
        self.assertGreater(distance, 0.5, "Truncated prefix should have high distance")

    def test_truncated_with_and_without_snp_consistent(self):
        """Truncated sequences should give consistent distances regardless of SNPs.

        This is the core inconsistency bug: without the fix, a truncated prefix
        gets 0.0 distance but a truncated prefix with a single SNP gets ~0.6
        because edlib produces a scattered alignment.
        """
        full_seq = generate_dna_sequence("consistency", 660)
        truncated_exact = full_seq[:220]

        # Add a single SNP near the middle of the truncated region
        trunc_list = list(truncated_exact)
        trunc_list[110] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[trunc_list[110]]
        truncated_snp = ''.join(trunc_list)

        dist_exact = calculate_adjusted_identity_distance(truncated_exact, full_seq)
        dist_snp = calculate_adjusted_identity_distance(truncated_snp, full_seq)

        # Both should be high distance (coverage-based fallback)
        self.assertGreater(dist_exact, 0.5)
        self.assertGreater(dist_snp, 0.5)

        # And they should be close to each other (not 0.0 vs 0.6)
        # The SNP shifts edlib's raw alignment slightly, so allow ~0.1 delta
        self.assertAlmostEqual(dist_exact, dist_snp, delta=0.1,
                               msg="Truncated distances should be consistent with/without SNP")

    def test_identical_same_length_zero_distance(self):
        """Identical sequences should still return 0.0 distance."""
        seq = generate_dna_sequence("identical", 660)
        self.assertEqual(calculate_adjusted_identity_distance(seq, seq), 0.0)

    def test_similar_full_length_low_distance(self):
        """Similar full-length sequences should have low distance (unchanged behavior)."""
        seq1 = generate_dna_sequence("similar", 660)
        # Introduce ~1% difference (7 mutations)
        seq2_list = list(seq1)
        for i in range(7):
            pos = i * 90 + 45
            seq2_list[pos] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq2_list[pos]]
        seq2 = ''.join(seq2_list)

        distance = calculate_adjusted_identity_distance(seq1, seq2)
        self.assertLess(distance, 0.05, "Similar full-length sequences should have low distance")

    def test_slightly_different_lengths_low_distance(self):
        """Sequences with slight length difference (650 vs 660) should have low distance."""
        seq_660 = generate_dna_sequence("slight_len", 660)
        # Trim 10bp — well above 50% coverage threshold
        seq_650 = seq_660[:650]

        distance = calculate_adjusted_identity_distance(seq_650, seq_660)
        self.assertLess(distance, 0.05,
                        "Slightly different length sequences should have low distance")

    def test_homopolymer_differences_still_adjusted(self):
        """Homopolymer length differences should still benefit from HP normalization."""
        # Two same-length sequences differing only by homopolymer run lengths
        base = generate_dna_sequence("hp_base", 100)

        # Insert a 3-base homopolymer run in the middle
        seq1 = base[:50] + "AAA" + base[50:]
        seq2 = base[:50] + "AAAAA" + base[50:]

        distance = calculate_adjusted_identity_distance(seq1, seq2)
        # HP normalization should keep distance low
        self.assertLess(distance, 0.05,
                        "Homopolymer differences should still be adjusted")


class TestTruncatedHACGrouping(unittest.TestCase):
    """Integration test: truncated sequences should not distort HAC grouping."""

    def test_truncated_sequence_does_not_split_full_length_group(self):
        """Full-length sequences >99% identical should stay grouped despite a truncated outlier.

        Mirrors the real bug: 4 full-length ITS sequences (~660bp) at >99% identity
        plus 1 truncated sequence (220bp) and 1 full-length with SNP at a position
        that causes edlib to scatter when aligning against the truncated sequence.
        """
        from speconsense.types import ConsensusInfo
        from speconsense.summarize.clustering import perform_hac_clustering

        # Base full-length sequence
        base_seq = generate_dna_sequence("hac_base", 660)

        # c1: identical to base
        seq_c1 = base_seq

        # c2: 1 SNP at position 400 (>99% identical to c1)
        seq_c2_list = list(base_seq)
        seq_c2_list[400] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c2_list[400]]
        seq_c2 = ''.join(seq_c2_list)

        # c3: 2 SNPs (still >99% identical)
        seq_c3_list = list(base_seq)
        seq_c3_list[200] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c3_list[200]]
        seq_c3_list[500] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c3_list[500]]
        seq_c3 = ''.join(seq_c3_list)

        # c4: truncated to 220bp (the problematic one)
        seq_c4 = base_seq[:220]

        # c5: 1 SNP at position 100, full-length
        seq_c5_list = list(base_seq)
        seq_c5_list[100] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c5_list[100]]
        seq_c5 = ''.join(seq_c5_list)

        # c6: 3 SNPs (still >99% identical)
        seq_c6_list = list(base_seq)
        seq_c6_list[150] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c6_list[150]]
        seq_c6_list[350] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c6_list[350]]
        seq_c6_list[550] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[seq_c6_list[550]]
        seq_c6 = ''.join(seq_c6_list)

        def make_ci(cluster_id, sequence, ric=100):
            return ConsensusInfo(
                sample_name="test",
                cluster_id=cluster_id,
                sequence=sequence,
                ric=ric,
                size=ric,
                file_path="/test",
            )

        consensus_list = [
            make_ci("c1", seq_c1, 200),
            make_ci("c2", seq_c2, 150),
            make_ci("c3", seq_c3, 120),
            make_ci("c4", seq_c4, 50),   # truncated
            make_ci("c5", seq_c5, 100),
            make_ci("c6", seq_c6, 80),
        ]

        # 97% identity threshold (standard HAC grouping)
        groups = perform_hac_clustering(consensus_list, variant_group_identity=0.97)

        # Find which group each cluster ended up in
        cluster_to_group = {}
        for group_id, members in groups.items():
            for member in members:
                cluster_to_group[member.cluster_id] = group_id

        # All full-length sequences should be in the same group
        full_length_ids = ["c1", "c2", "c3", "c5", "c6"]
        full_length_groups = {cluster_to_group[cid] for cid in full_length_ids}
        self.assertEqual(len(full_length_groups), 1,
                         f"Full-length sequences should all be in one group, "
                         f"but got groups: {full_length_groups}")

        # The truncated sequence should be in a different group
        truncated_group = cluster_to_group["c4"]
        full_length_group = full_length_groups.pop()
        self.assertNotEqual(truncated_group, full_length_group,
                            "Truncated sequence should be in a separate group")


if __name__ == '__main__':
    unittest.main()
