"""Unit tests for the metadata CER reproduction helpers in SpecimenClusterer."""

from speconsense.core.clusterer import SpecimenClusterer


def test_build_identity_group_summary_groups_clusters_by_id():
    cluster_dicts = [
        {'cluster_id': 'c1', 'identity_group_id': 'g0', 'cer_group_N': 950},
        {'cluster_id': 'c2', 'identity_group_id': 'g0', 'cer_group_N': 950},
        {'cluster_id': 'c3', 'identity_group_id': 'g1', 'cer_group_N': 50},
    ]
    summary = SpecimenClusterer._build_identity_group_summary(cluster_dicts)
    assert summary == [
        {'group_id': 'g0', 'group_N': 950, 'members': ['c1', 'c2']},
        {'group_id': 'g1', 'group_N': 50, 'members': ['c3']},
    ]


def test_build_identity_group_summary_skips_clusters_with_no_group():
    cluster_dicts = [
        {'cluster_id': 'c1'},   # No identity_group_id; produced before validation
        {'cluster_id': 'c2', 'identity_group_id': 'g0', 'cer_group_N': 100},
    ]
    summary = SpecimenClusterer._build_identity_group_summary(cluster_dicts)
    assert summary == [
        {'group_id': 'g0', 'group_N': 100, 'members': ['c2']},
    ]


def test_build_variant_record_anchor_no_cer_details():
    cluster_dict = {
        'cluster_id': 'c1',
        'identity_group_id': 'g0',
        'cer_group_N': 950,
        'read_ids': set(range(900)),
        'cer_factor': None,
        'cer_pstar': None,
        'cer_details': None,
    }
    record = SpecimenClusterer._build_variant_record(cluster_dict)
    assert record == {
        'cluster_id': 'c1',
        'identity_group': 'g0',
        'group_rank': None,
        'variant_rank': None,
        'M': 900,
        'N': 950,
        'K': None,
        'context_tags': None,
        'q_ctx_per_position': None,
        'compared_against_idx': None,
        'cer_factor': None,
        'cer_pstar': None,
        'err_factor': None,
        'err_factor_obs_sum': None,
        'err_factor_exp_sum': None,
        'err_factor_cols': None,
    }


def test_build_variant_record_with_cer_details():
    cluster_dict = {
        'cluster_id': 'c2',
        'identity_group_id': 'g0',
        'cer_group_N': 950,
        'read_ids': set(range(40)),
        'cer_factor': 4.2,
        'cer_pstar': 0.0248,
        'cer_details': {
            'K': 2,
            'tags': ['non-hp-sub', 'hp-l3-A-del1'],
            'q_ctx': [0.0059, 0.0097],
            'ref_idx': 0,
        },
    }
    record = SpecimenClusterer._build_variant_record(cluster_dict)
    assert record == {
        'cluster_id': 'c2',
        'identity_group': 'g0',
        'group_rank': None,
        'variant_rank': None,
        'M': 40,
        'N': 950,
        'K': 2,
        'context_tags': ['non-hp-sub', 'hp-l3-A-del1'],
        'q_ctx_per_position': [0.0059, 0.0097],
        'compared_against_idx': 0,
        'cer_factor': 4.2,
        'cer_pstar': 0.0248,
        'err_factor': None,
        'err_factor_obs_sum': None,
        'err_factor_exp_sum': None,
        'err_factor_cols': None,
    }


def test_assign_identity_ranks_orders_groups_and_variants_by_size():
    cluster_dicts = [
        {'cluster_id': 'c1', 'identity_group_id': 'g0', 'read_ids': set(range(40))},
        {'cluster_id': 'c2', 'identity_group_id': 'g1', 'read_ids': set(range(25))},
        {'cluster_id': 'c3', 'identity_group_id': 'g0', 'read_ids': set(range(10))},
        {'cluster_id': 'c4', 'identity_group_id': 'g1', 'read_ids': set(range(8))},
        {'cluster_id': 'c5', 'identity_group_id': 'g2', 'read_ids': set(range(50))},
    ]
    SpecimenClusterer._assign_identity_ranks(cluster_dicts)
    # g2 anchor=50 (largest), g0 anchor=40, g1 anchor=25 — so ranks 1,2,3.
    ranks = {c['cluster_id']: (c['identity_group_rank'], c['identity_variant_rank'])
             for c in cluster_dicts}
    assert ranks == {
        'c5': (1, 1),
        'c1': (2, 1),
        'c3': (2, 2),
        'c2': (3, 1),
        'c4': (3, 2),
    }


def test_assign_identity_ranks_skips_ungrouped_clusters():
    cluster_dicts = [
        {'cluster_id': 'c1', 'identity_group_id': 'g0', 'read_ids': set(range(10))},
        {'cluster_id': 'c2', 'read_ids': set(range(5))},  # no identity_group_id
    ]
    SpecimenClusterer._assign_identity_ranks(cluster_dicts)
    assert cluster_dicts[0]['identity_group_rank'] == 1
    assert cluster_dicts[0]['identity_variant_rank'] == 1
    assert 'identity_group_rank' not in cluster_dicts[1]
    assert 'identity_variant_rank' not in cluster_dicts[1]


def test_build_variant_record_low_factor_with_cer_details():
    cluster_dict = {
        'cluster_id': 'c3',
        'identity_group_id': 'g0',
        'cer_group_N': 950,
        'read_ids': set(range(5)),
        'cer_factor': 0.34,
        'cer_pstar': 0.0020,
        'cer_details': {
            'K': 1,
            'tags': ['non-hp-sub'],
            'q_ctx': [0.0059],
            'ref_idx': 0,
        },
    }
    record = SpecimenClusterer._build_variant_record(cluster_dict)
    assert record['cer_factor'] == 0.34
    assert record['M'] == 5
