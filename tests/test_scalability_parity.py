"""Parity tests for the scalability extensions added to the post-phasing pipeline.

The dense ("brute-force") and sparse ("vsearch-backed") code paths must produce
identical or near-identical results on the same input. These tests construct
controlled fixtures and exercise both paths to detect regressions in the
sparse code paths.

Tests are skipped when vsearch is not available, since the sparse paths
depend on it.
"""

import shutil
import statistics
import string
import tempfile
from typing import Dict, List, Set

import pytest

VSEARCH_AVAILABLE = shutil.which("vsearch") is not None
pytestmark = pytest.mark.skipif(not VSEARCH_AVAILABLE,
                                reason="vsearch not on PATH")


def _rng_seq(seed: int, length: int = 600, alphabet: str = "ACGT") -> str:
    """Deterministic pseudo-random sequence. Avoids low-complexity repeats that
    confuse vsearch's kmer heuristic."""
    import random
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(length))


def _mutate(base: str, seed: int, n_subs: int = 5, alphabet: str = "ACGT") -> str:
    """Deterministic per-position substitution."""
    out = list(base)
    for k in range(n_subs):
        pos = (seed * 1009 + k * 31) % len(out)
        new_base = alphabet[(seed + k * 7) % 4]
        if out[pos] != new_base:
            out[pos] = new_base
    return "".join(out)


def _build_consensuses(n_per_group: int = 20) -> List[str]:
    """Build cluster consensuses falling into 3 distinct identity groups.

    Each group has a 600bp pseudo-random prototype; members are derived by 4
    substitutions per member (~99% identity within a group, ~25% identity
    between groups since prototypes are independent random sequences).
    """
    proto_A = _rng_seq(seed=11)
    proto_B = _rng_seq(seed=22)
    proto_C = _rng_seq(seed=33)

    consensuses: List[str] = []
    for proto in (proto_A, proto_B, proto_C):
        for k in range(n_per_group):
            consensuses.append(_mutate(proto, seed=k + 1, n_subs=4))
    return consensuses


def _make_clusterer(scale_threshold: int, max_threads: int = 1):
    from speconsense.core import SpecimenClusterer
    output_dir = tempfile.mkdtemp(prefix="spec_par_")
    return SpecimenClusterer(
        output_dir=output_dir,
        scale_threshold=scale_threshold,
        max_threads=max_threads,
        group_identity=0.85,
        min_hp_length=6,
    )


def _run_form_identity_groups(consensuses: List[str], scale_threshold: int) -> Dict[int, frozenset]:
    """Run _form_identity_groups via a clusterer with the given scale_threshold.

    Returns membership as {sorted_repr_id -> frozenset(member_indices)} so the
    arbitrary scipy label numbering doesn't trip the comparison.
    """
    clusterer = _make_clusterer(scale_threshold=scale_threshold)
    subclusters = [{'read_ids': set([f'r{i}'])} for i in range(len(consensuses))]
    consensus_dict = {i: c for i, c in enumerate(consensuses)}
    raw_groups = clusterer._form_identity_groups(subclusters, consensus_dict)
    return {min(idxs): frozenset(idxs) for idxs in raw_groups.values()}


# ---------- Phase 4a parallelization parity ----------

def test_filter_noisy_clusters_parallel_matches_sequential():
    """Phase 4a: ProcessPoolExecutor path produces same disband decisions as
    the sequential path for the same input."""
    from speconsense.core import SpecimenClusterer

    # Build 30 small clusters each with 4 reads (below phasing_floor=6 default).
    # Half are concordant (same sequence), half are noisy (random reads).
    n_clusters = 30
    base = ("ACGTACGTACGT" * 50)  # 600bp

    def build_clusterer(threads: int):
        c = SpecimenClusterer(
            output_dir=tempfile.mkdtemp(prefix="spec_4a_"),
            scale_threshold=0,
            max_threads=threads,
            min_variant_count=3,  # phasing_floor = 6
        )
        return c

    def populate(c):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        subclusters = []
        rid = 0
        for k in range(n_clusters):
            read_ids: Set[str] = set()
            if k % 2 == 0:
                # Concordant cluster: 4 identical reads.
                seq = _mutate(base, seed=k + 1, n_subs=2)
                for _ in range(4):
                    name = f"r{rid}"
                    rid += 1
                    c.sequences[name] = seq
                    rec = SeqRecord(Seq(seq), id=name)
                    rec.letter_annotations["phred_quality"] = [30] * len(seq)
                    c.records[name] = rec
                    read_ids.add(name)
            else:
                # Noisy: 4 reads with different mutation profiles.
                for j in range(4):
                    seq = _mutate(base, seed=k * 100 + j * 11, n_subs=20)
                    name = f"r{rid}"
                    rid += 1
                    c.sequences[name] = seq
                    rec = SeqRecord(Seq(seq), id=name)
                    rec.letter_annotations["phred_quality"] = [30] * len(seq)
                    c.records[name] = rec
                    read_ids.add(name)
            subclusters.append({'read_ids': read_ids})
        c.total_input_reads = rid
        return subclusters

    c_seq = build_clusterer(threads=1)
    sub_seq = populate(c_seq)
    out_seq = c_seq._filter_noisy_clusters([{'read_ids': s['read_ids'].copy()} for s in sub_seq])
    disbanded_seq = c_seq.discarded_read_ids.copy()

    c_par = build_clusterer(threads=4)
    sub_par = populate(c_par)
    out_par = c_par._filter_noisy_clusters([{'read_ids': s['read_ids'].copy()} for s in sub_par])
    disbanded_par = c_par.discarded_read_ids.copy()

    # Same number of survivors.
    assert len(out_seq) == len(out_par), \
        f"survivor count differs: seq={len(out_seq)} par={len(out_par)}"

    # Same disbanded read set (read IDs are deterministic).
    assert disbanded_seq == disbanded_par, "disbanded read sets differ"


# ---------- _form_identity_groups parity ----------

def test_form_identity_groups_dense_vs_sparse_match():
    """The sparse vsearch-backed path must produce identical groupings to the
    dense O(n²) path on a 60-cluster input."""
    consensuses = _build_consensuses(n_per_group=20)  # 60 total
    assert len(consensuses) == 60

    # Dense: scale_threshold=0 → scalability disabled.
    dense_groups = _run_form_identity_groups(consensuses, scale_threshold=0)

    # Sparse: scale_threshold=10 → activation gate (n>50 AND n>=10) trips.
    sparse_groups = _run_form_identity_groups(consensuses, scale_threshold=10)

    assert dense_groups == sparse_groups, \
        f"identity-group memberships differ:\n  dense:  {dense_groups}\n  sparse: {sparse_groups}"


def test_form_identity_groups_sparse_recovers_three_groups():
    """The sparse path must recover the three obvious prototypes."""
    consensuses = _build_consensuses(n_per_group=20)
    sparse_groups = _run_form_identity_groups(consensuses, scale_threshold=10)
    sizes = sorted(len(s) for s in sparse_groups.values())
    assert sizes == [20, 20, 20], f"expected three groups of 20, got sizes {sizes}"


def test_form_identity_groups_skips_below_threshold():
    """Below the n>50 threshold, the dense path runs even with scalability
    enabled — sanity check that small-input behavior is unchanged."""
    consensuses = _build_consensuses(n_per_group=10)  # 30 total, below n>50
    dense_groups = _run_form_identity_groups(consensuses, scale_threshold=0)
    sparse_attempt_groups = _run_form_identity_groups(consensuses, scale_threshold=10)
    assert dense_groups == sparse_attempt_groups
