"""Tests for speconsense/summarize/tree.py — per-specimen variant tree.

Covers:
- Parent selection prefers the larger-size peer with highest pairwise identity.
- Anchor (largest) is the root.
- Status mix (passed + .ns + .lq) renders together, grouped by core gid.
- Singleton groups produce just an anchor line.
- File contents include the expected sections, tree glyphs, and diff lines.
"""

import os
import tempfile

from speconsense.types import ConsensusInfo
from speconsense.summarize.tree import (
    _Node,
    _build_tree,
    write_specimen_variant_tree,
)


def _ci(name: str, sequence: str, size: int, ric: int = None,
        gid: int = 1, vid: int = 1, cer=None, err=None) -> ConsensusInfo:
    return ConsensusInfo(
        sample_name=name,
        cluster_id=name,
        sequence=sequence,
        ric=ric if ric is not None else size,
        size=size,
        file_path="/tmp/fake.fastq",
        cer_factor=cer,
        err_factor=err,
        group_rank=gid,
        variant_rank=vid,
    )


def test_build_tree_picks_higher_identity_peer():
    """When two larger peers exist, the candidate hangs off the closer one."""
    # Build sequences with non-HP differences (so HP normalization doesn't
    # collapse them) at distinct positions to make distances clearly ordered.
    anchor_seq = "ACGTACGTACGT" * 20  # 240 bp
    # near_seq differs from anchor at 4 substitution positions
    near_seq = list(anchor_seq)
    for pos in (10, 50, 100, 150):
        near_seq[pos] = 'T' if anchor_seq[pos] != 'T' else 'A'
    near_seq = "".join(near_seq)
    # cand_seq starts from near_seq and adds one more substitution.
    # → cand vs near = 1 edit ; cand vs anchor = 5 edits.
    cand_seq = list(near_seq)
    cand_seq[200] = 'A' if near_seq[200] != 'A' else 'G'
    cand_seq = "".join(cand_seq)

    nodes = [
        _Node(info=_ci("S-1.v1", anchor_seq, size=100), status='passed', short_name="-1.v1"),
        _Node(info=_ci("S-1.v2", near_seq, size=50), status='passed', short_name="-1.v2"),
        _Node(info=_ci("S-1.v3", cand_seq, size=10), status='passed', short_name="-1.v3"),
    ]
    _build_tree(nodes, hp_normalization_length=6)

    # cand_seq differs from anchor by ~5 edits, from near_seq by ~1 edit.
    # Higher-identity peer for cand_seq is near_seq (index 1), not anchor (index 0).
    assert nodes[2].parent_idx == 1, \
        f"Expected parent=1 (near_seq), got {nodes[2].parent_idx}"


def test_build_tree_singleton_no_parent():
    nodes = [
        _Node(info=_ci("S-1.v1", "ACGT" * 50, size=100), status='passed', short_name="-1.v1"),
    ]
    _build_tree(nodes, hp_normalization_length=6)
    assert nodes[0].parent_idx is None
    assert nodes[0].children == []


def test_write_specimen_variant_tree_renders_passed_ns_lq_together():
    base = "ACGT" * 80
    passed = [
        _ci("SP-1.v1", base, size=100, gid=1),
        _ci("SP-1.v2", base[:-1] + "T", size=50, gid=1, cer=4.2, err=0.7),
    ]
    ns = [
        _ci("SP-1.v3", base[:-2] + "TT", size=10, gid=1, cer=0.85, err=0.9),
    ]
    lq = [
        _ci("SP-1.v4", "G" * len(base), size=5, gid=1, cer=0.5, err=1.7),
    ]

    with tempfile.TemporaryDirectory() as tmp:
        write_specimen_variant_tree(
            specimen_id="SP",
            passed=passed,
            ns=ns,
            lq=lq,
            output_dir=tmp,
            hp_normalization_length=6,
        )
        out_path = os.path.join(tmp, "SP.txt")
        assert os.path.exists(out_path), "Tree file not written"
        contents = open(out_path).read()

    assert "VARIANT TREE — SP" in contents
    assert "4 variants total" in contents
    assert "Group 1" in contents
    assert "[anchor]" in contents
    # Status indicators present for each non-anchor variant
    assert " passed " in contents
    assert " ns " in contents
    assert " lq " in contents
    # Tree glyphs
    assert "├─" in contents or "└─" in contents
    # Diff line format
    assert "vs " in contents and "id) —" in contents


def test_groups_separated_by_gid():
    base1 = "ACGT" * 80
    base2 = "GGGG" * 80
    passed = [
        _ci("SP-1.v1", base1, size=100, gid=1),
        _ci("SP-2.v1", base2, size=80, gid=2),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        write_specimen_variant_tree(
            specimen_id="SP",
            passed=passed,
            ns=[], lq=[],
            output_dir=tmp,
        )
        contents = open(os.path.join(tmp, "SP.txt")).read()

    assert "Group 1" in contents
    assert "Group 2" in contents
    assert "2 identity groups" in contents


def test_no_groups_no_file():
    """Variants without group_rank produce no file (legacy / malformed input)."""
    base = "ACGT" * 80
    passed = [_ci("SP-1.v1", base, size=100, gid=None)]
    with tempfile.TemporaryDirectory() as tmp:
        write_specimen_variant_tree(
            specimen_id="SP",
            passed=passed,
            ns=[], lq=[],
            output_dir=tmp,
        )
        # No tree file should be written.
        assert not os.path.exists(os.path.join(tmp, "SP.txt"))
