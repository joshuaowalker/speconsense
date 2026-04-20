"""Unit tests for detect_rid_outliers."""

import pytest

from speconsense.outliers import detect_rid_outliers


class TestCleanClusters:
    """Clusters with no outliers should return []."""

    def test_perfectly_identical(self):
        assert detect_rid_outliers([0.99, 0.99, 0.99]) == []

    def test_small_spread(self):
        assert detect_rid_outliers([0.985, 0.99, 0.995]) == []

    def test_smooth_gradient(self):
        assert detect_rid_outliers([0.98, 0.985, 0.99, 0.995, 1.0]) == []

    def test_n_10_clean(self):
        rids = [0.985, 0.988, 0.990, 0.991, 0.992, 0.993, 0.995, 0.996, 0.998, 1.0]
        assert detect_rid_outliers(rids) == []


class TestSmallNOutlierDetection:
    """At N=3,4 the gap rule is the main tool."""

    def test_n3_clear_outlier(self):
        # [95, 99, 99.5] — one clear outlier at bottom
        assert detect_rid_outliers([0.95, 0.99, 0.995]) == [0]

    def test_n3_mild_outlier_boundary(self):
        # [95.0, 98.0, 98.5] — ratio is 3.0/0.5=6.0, > gap_factor=2.5; drop
        # from median (98.0) is 3pp, above min_drop_from_median=2pp.
        assert detect_rid_outliers([0.95, 0.98, 0.985]) == [0]

    def test_n3_no_outlier(self):
        # [98.5, 99.0, 99.5] — smooth spacing, ratio=1.0
        assert detect_rid_outliers([0.985, 0.99, 0.995]) == []

    def test_n4_one_outlier(self):
        # [94, 99, 99.2, 99.5]
        assert detect_rid_outliers([0.94, 0.99, 0.992, 0.995]) == [0]

    def test_n4_no_outlier(self):
        # Smooth distribution
        assert detect_rid_outliers([0.97, 0.98, 0.99, 0.995]) == []


class TestLargerNMadDetection:
    """At N>=5 MAD should reliably flag outliers."""

    def test_n6_one_outlier(self):
        # Mimics 235757835: [94.42, 98.47, 99.02, 99.03, 99.16, 99.72]
        rids = [0.9442, 0.9847, 0.9902, 0.9903, 0.9916, 0.9972]
        flagged = detect_rid_outliers(rids)
        assert 0 in flagged  # worst read flagged

    def test_n7_one_outlier(self):
        # Mimics 234818038: [97.00, 98.56, 98.93, 99.16, 99.52, 99.64, 99.76]
        rids = [0.9700, 0.9856, 0.9893, 0.9916, 0.9952, 0.9964, 0.9976]
        flagged = detect_rid_outliers(rids)
        assert 0 in flagged

    def test_n10_one_clear_outlier(self):
        # Mimics 187834361: worst at 88.6, cluster otherwise tight
        rids = [1.0, 0.9986, 0.9973, 0.9919, 0.9905, 0.9878, 0.9771,
                0.9769, 0.9289, 0.8860]
        flagged = detect_rid_outliers(rids)
        # Index 9 is the 0.886 outlier (by input position)
        assert 9 in flagged
        # 0.9289 may or may not be flagged (boundary case); that's fine

    def test_n10_no_outlier(self):
        # Tight cluster, even with slight variation
        rids = [0.99, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997,
                0.998, 0.999]
        assert detect_rid_outliers(rids) == []


class TestEdgeCases:
    """Degenerate and boundary cases."""

    def test_n_below_3_returns_empty(self):
        assert detect_rid_outliers([]) == []
        assert detect_rid_outliers([0.99]) == []
        assert detect_rid_outliers([0.99, 0.95]) == []

    def test_all_identical_reads(self):
        # All reads identical → no outliers, no divide-by-zero
        assert detect_rid_outliers([0.99, 0.99, 0.99, 0.99]) == []

    def test_only_one_read_differs(self):
        # Most reads identical, one slightly different
        rids = [0.99, 0.99, 0.99, 0.99, 0.985]
        # The gap rule fires: spread_top=0 (all four 0.99 identical),
        # gap_bottom=0.005 > 5*min_mad=0.010? No, 0.005 < 0.010.
        # So this specific case should NOT flag (small difference).
        assert detect_rid_outliers(rids) == []

    def test_identical_rest_with_clear_outlier(self):
        # Four reads at 0.99, one at 0.90 — clear outlier despite
        # degenerate spread_top = 0.
        rids = [0.99, 0.99, 0.99, 0.99, 0.90]
        flagged = detect_rid_outliers(rids)
        assert 4 in flagged

    def test_two_low_outliers_n5(self):
        # MAD should catch both: [0.85, 0.85, 0.99, 0.99, 0.995]
        rids = [0.85, 0.85, 0.99, 0.99, 0.995]
        flagged = detect_rid_outliers(rids)
        # Both low reads flagged by MAD
        assert 0 in flagged and 1 in flagged

    def test_indices_preserved_under_reordering(self):
        # Same values in different input order should flag the same rids
        rids1 = [0.90, 0.99, 0.995, 0.99, 0.995]
        flagged1 = detect_rid_outliers(rids1)
        assert flagged1 == [0]
        rids2 = [0.995, 0.99, 0.90, 0.99, 0.995]  # moved outlier to index 2
        flagged2 = detect_rid_outliers(rids2)
        assert flagged2 == [2]

    def test_never_flags_high_outlier(self):
        # One read extremely high, others clustered at 0.95 — high read
        # should not be flagged (we only care about low-rid outliers).
        rids = [0.95, 0.95, 0.95, 0.95, 1.00]
        flagged = detect_rid_outliers(rids)
        assert 4 not in flagged

    def test_percentages_equivalent_to_fractions(self):
        # Should work regardless of whether input is in [0,1] or [0,100]
        frac = detect_rid_outliers([0.95, 0.99, 0.995])
        pct = detect_rid_outliers([95.0, 99.0, 99.5], min_mad=0.2)
        assert frac == pct


class TestThresholdSensitivity:
    """Knobs work as documented."""

    def test_stricter_z_threshold_flags_fewer(self):
        # One mild outlier; should be flagged at default but not stricter
        rids = [0.95, 0.98, 0.99, 0.99, 0.995]
        default = detect_rid_outliers(rids)
        strict = detect_rid_outliers(rids, modified_z_threshold=10.0, gap_factor=100.0)
        assert len(default) >= len(strict)

    def test_min_drop_safety(self):
        # Drop from median = 1.5pp, below default 2pp safety.
        rids = [0.975, 0.99, 0.99, 0.99]
        assert detect_rid_outliers(rids) == []
        # Relaxing the safety to 0.5pp lets it flag.
        assert detect_rid_outliers(rids, min_drop_from_median=0.005) == [0]

    def test_smaller_gap_factor_flags_more(self):
        # Borderline gap case; drops chosen to clear the min_drop_from_median
        # safety (needs at least 0.02 below median).
        rids = [0.96, 0.99, 0.995]  # gap_bottom=0.03, spread_top=0.005, ratio=6.0
        # At modified_z_threshold=10.0 MAD rule won't fire; gap rule alone.
        # gap_factor=10.0 → not flagged (6.0 < 10.0);
        # gap_factor=2.5 → flagged (6.0 > 2.5).
        assert detect_rid_outliers(rids, gap_factor=10.0, modified_z_threshold=10.0) == []
        assert detect_rid_outliers(rids, gap_factor=2.5, modified_z_threshold=10.0) == [0]
