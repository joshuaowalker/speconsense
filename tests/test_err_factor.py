"""Unit tests for compute_cluster_err_factor."""

from speconsense.msa import compute_cluster_err_factor
from speconsense.qctx import load_table


QCTX = load_table("dorado-v5.0")


def _msa(reads, consensus):
    """Build a minimal MSA FASTA string."""
    parts = []
    for name, seq in reads:
        parts.append(f">{name}\n{seq}")
    parts.append(f">Consensus\n{consensus}")
    return "\n".join(parts) + "\n"


def test_err_factor_none_on_empty_msa():
    ef, obs, exp, cols = compute_cluster_err_factor("", QCTX)
    assert ef is None
    assert cols == 0


def test_err_factor_near_zero_for_perfectly_matching_reads():
    consensus = "ACGTACGTACGTACGTACGT"
    reads = [("r1", consensus), ("r2", consensus), ("r3", consensus)]
    ef, obs, exp, cols = compute_cluster_err_factor(_msa(reads, consensus), QCTX)
    assert ef is not None
    assert obs == 0.0
    assert exp > 0
    assert cols == len(consensus)
    assert ef == 0.0


def test_err_factor_elevated_for_heterogeneous_reads():
    # Consensus is all A's of length 30 (no HP run intent; just to test)
    # Use mixed bases so there's no HP
    consensus = "ACGTACGTACGTACGTACGTACGTACGTAC"
    # One read perfect, two reads with many mismatches
    reads = [
        ("r1", consensus),
        ("r2", "ATGTACCTACGTACGTACGTATGTACGTAC"),  # 3 substitutions
        ("r3", "ACGTACGAACGTACGTACATACGTACGTAC"),  # 2 substitutions
    ]
    ef, obs, exp, cols = compute_cluster_err_factor(_msa(reads, consensus), QCTX)
    assert ef is not None
    assert ef > 1.0  # Observed errors > q_ctx expected for this noisy cluster


def test_err_factor_handles_hp_runs():
    # Cluster with an HP run; one read has a deletion in the HP
    consensus = "ACGTAAAAAGTACGT"  # length 15, 5-A HP run
    reads = [
        ("r1", consensus),
        ("r2", "ACGTAAAA-GTACGT"),   # HP length deletion
        ("r3", consensus),
    ]
    ef, obs, exp, cols = compute_cluster_err_factor(_msa(reads, consensus), QCTX)
    assert ef is not None
    # Only one mismatch across 3 reads at one of 15 cols; HP deletion.
    assert 0.0 < ef < 3.0


def test_err_factor_obs_exp_sums_reconcile():
    consensus = "ACGT" * 10
    reads = [("r1", consensus), ("r2", consensus.replace("A", "T", 1)), ("r3", consensus)]
    ef, obs, exp, cols = compute_cluster_err_factor(_msa(reads, consensus), QCTX)
    # 1 mismatch in 3 reads at 1 column -> obs_sum = 1/3
    assert abs(obs - 1/3) < 1e-9
    assert cols == len(consensus)
    assert exp > 0
    assert abs(ef - obs / exp) < 1e-9
