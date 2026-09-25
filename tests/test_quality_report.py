"""Unit tests for quality_report rendering, focused on the q_ctx calibration section."""

from speconsense.types import ConsensusInfo
from speconsense.quality_report import (
    CALIBRATION_HIGH,
    CALIBRATION_LOW,
    MIN_CALIBRATION_COLS,
    _Cluster,
    _Specimen,
    _render_qctx_calibration,
)


def _make_specimen(name: str, clusters):
    sp = _Specimen(name=name)
    for state, obs, exp, cols in clusters:
        info = ConsensusInfo(
            sample_name=f"{name}-1.v1",
            cluster_id="1.v1",
            sequence="ACGT",
            ric=10,
            size=10,
            file_path="",
            err_factor=(obs / exp) if exp else None,
            err_factor_obs_sum=obs,
            err_factor_exp_sum=exp,
            err_factor_cols=cols,
        )
        sp.clusters.append(_Cluster(info=info, state=state, specimen=name, short="1.v1"))
    return sp


def _render_with(clusters):
    return _render_qctx_calibration({"sp": _make_specimen("sp", clusters)})


def test_calibration_within_range_no_warning():
    # Pooled = 1.0 across enough cols.
    out = _render_with([("passed", 100.0, 100.0, MIN_CALIBRATION_COLS)])
    assert "Pooled obs/exp:    1.00" in out
    assert "Within expected range" in out
    assert "WARNING" not in out


def test_calibration_low_warning():
    # Pooled = 0.50, below CALIBRATION_LOW (0.70).
    out = _render_with([("passed", 50.0, 100.0, MIN_CALIBRATION_COLS)])
    assert "Pooled obs/exp:    0.50" in out
    assert "WARNING" in out
    assert "overestimate" in out


def test_calibration_high_warning():
    # Pooled = 2.00, above CALIBRATION_HIGH (1.40).
    out = _render_with([("passed", 200.0, 100.0, MIN_CALIBRATION_COLS)])
    assert "Pooled obs/exp:    2.00" in out
    assert "WARNING" in out
    assert "underestimate" in out


def test_calibration_insufficient_data():
    # Below MIN_CALIBRATION_COLS — no warning, no pooled value.
    out = _render_with([("passed", 50.0, 100.0, 100)])
    assert "Insufficient data" in out
    assert "WARNING" not in out
    assert "Pooled obs/exp" not in out


def test_calibration_pools_all_states():
    # All three buckets are pooled together for the headline ratio.
    # Synthetic case mirroring the empirical observation that a clean passing
    # subset (top-truncated by .lq routing) would alone read low (0.69), but
    # adding back .ns and .lq recovers ~1.0.
    out = _render_with([
        ("passed", 69.0, 100.0, MIN_CALIBRATION_COLS // 3),     # 0.69 alone
        ("ns",     63.0, 100.0, MIN_CALIBRATION_COLS // 3),     # 0.63 alone
        ("lq",    331.0, 100.0, MIN_CALIBRATION_COLS // 3 + 2),  # 3.31 alone
    ])
    # Combined obs/exp = 463/300 ~= 1.54; ensure the headline reflects all
    # three buckets, not just "passed".
    assert "Pooled obs/exp:    1.54" in out
    assert "all states:" in out


def test_calibration_breakdown_shows_each_state():
    # Each non-empty state appears in the per-state breakdown with its own
    # pooled obs/exp value.
    out = _render_with([
        ("passed", 50.0, 100.0, MIN_CALIBRATION_COLS // 3),
        ("ns",    100.0, 100.0, MIN_CALIBRATION_COLS // 3),
        ("lq",    300.0, 100.0, MIN_CALIBRATION_COLS // 3 + 2),
    ])
    assert "Per-state breakdown:" in out
    # passed alone = 0.50, ns alone = 1.00, lq alone = 3.00
    assert "passed" in out and "0.50" in out
    assert ".ns" in out and "1.00" in out
    assert ".lq" in out and "3.00" in out


def test_calibration_breakdown_omits_empty_states():
    # If only one state has data, the others are not listed.
    out = _render_with([("passed", 100.0, 100.0, MIN_CALIBRATION_COLS)])
    assert "passed" in out
    assert ".ns" not in out
    assert ".lq" not in out


def test_calibration_skips_clusters_without_metadata():
    # Clusters whose metadata JSON wasn't loadable have None for obs/exp/cols
    # and must be silently excluded from the pooled estimator.
    sp = _Specimen(name="sp")
    info_with = ConsensusInfo(
        sample_name="sp-1.v1", cluster_id="1.v1", sequence="A",
        ric=1, size=1, file_path="",
        err_factor_obs_sum=100.0, err_factor_exp_sum=100.0,
        err_factor_cols=MIN_CALIBRATION_COLS,
    )
    info_without = ConsensusInfo(
        sample_name="sp-1.v2", cluster_id="1.v2", sequence="A",
        ric=1, size=1, file_path="",
    )
    sp.clusters.append(_Cluster(info=info_with, state="passed", specimen="sp", short="1.v1"))
    sp.clusters.append(_Cluster(info=info_without, state="passed", specimen="sp", short="1.v2"))
    out = _render_qctx_calibration({"sp": sp})
    # Pooled reflects only the cluster with metadata, n=1.
    assert "Pooled obs/exp:    1.00" in out
    assert "1 clusters" in out


def test_calibration_includes_filtered_state():
    # .filtered records should participate in the q_ctx calibration just
    # like .ns and .lq — they are real clusters with valid err_factor data.
    out = _render_with([
        ("passed", 100.0, 100.0, MIN_CALIBRATION_COLS // 2),
        ("filtered", 100.0, 100.0, MIN_CALIBRATION_COLS // 2 + 2),
    ])
    assert "Pooled obs/exp:    1.00" in out
    assert ".filtered" in out


def test_calibration_thresholds_are_module_constants():
    # Defaults are intentionally weakly-held; verify they exist and are sane.
    assert 0 < CALIBRATION_LOW < 1.0 < CALIBRATION_HIGH
    assert MIN_CALIBRATION_COLS > 0


# ---------------------------------------------------------------------------
# write_quality_report: state accounting
# ---------------------------------------------------------------------------

from speconsense.quality_report import write_quality_report  # noqa: E402


def _info(name, size=10, obs=10.0, exp=10.0, cols=1000):
    return ConsensusInfo(
        sample_name=name, cluster_id=name.rsplit("-", 1)[-1], sequence="ACGT",
        ric=size, size=size, file_path="",
        err_factor=obs / exp,
        err_factor_obs_sum=obs, err_factor_exp_sum=exp, err_factor_cols=cols,
    )


def _write(tmp_path, **lists):
    write_quality_report([], [], str(tmp_path), str(tmp_path), **lists)
    return (tmp_path / "quality_report.txt").read_text()


def test_filtered_clusters_are_not_also_counted_as_passed(tmp_path):
    # Selection drops sp-1.v2 after load-time routing put it on the pass
    # list; pruning moves the .ns record sp-2.v1 to .filtered too. Each
    # cluster must be counted once, in its final state.
    passed = [_info("sp-1.v1"), _info("sp-1.v2")]
    ns = [_info("sp-2.v1")]
    filtered = [_info("sp-1.v2"), _info("sp-2.v1")]
    out = _write(tmp_path, consensus_list=passed, ns_list=ns,
                 filtered_list=filtered)
    assert ("Total clusters:    3  (1 passed | 0 routed to .ns | "
            "0 routed to .lq | 2 routed to .filtered)") in out
    assert "Variants/specimen: median 1" in out


def test_filtered_dedupe_keeps_calibration_sums_single_counted(tmp_path):
    cols = MIN_CALIBRATION_COLS
    passed = [_info("sp-1.v1", obs=100.0, exp=100.0, cols=cols),
              _info("sp-1.v2", obs=300.0, exp=100.0, cols=cols)]
    filtered = [_info("sp-1.v2", obs=300.0, exp=100.0, cols=cols)]
    out = _write(tmp_path, consensus_list=passed, filtered_list=filtered)
    # (100 + 300) / (100 + 100) = 2.00; double counting would give 3.50.
    assert "Pooled obs/exp:    2.00 (all states: 2 clusters" in out


def test_chimera_state_reported_when_present(tmp_path):
    cols = MIN_CALIBRATION_COLS
    out = _write(tmp_path,
                 consensus_list=[_info("sp-1.v1", cols=cols)],
                 chimera_list=[_info("sp-1.v3", cols=cols)])
    assert ("(1 passed | 0 routed to .ns | 0 routed to .lq | "
            "1 routed to .chimera | 0 routed to .filtered)") in out
    assert "1 .chimera (50.0%)" in out
    assert "    .chimera " in out  # calibration per-state row


def test_chimera_state_absent_by_default(tmp_path):
    out = _write(tmp_path, consensus_list=[_info("sp-1.v1")])
    assert ".chimera" not in out
