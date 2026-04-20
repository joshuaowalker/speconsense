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


if __name__ == '__main__':
    unittest.main()
