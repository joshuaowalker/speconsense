"""Tests for two-parent recombinant (PCR chimera) detection and routing.

Covers the pure scan functions in speconsense.chimera, the SPOA-backed
detect_chimera entry point, core's Phase 11 group detection + header
emission, and summarize's chimera= parsing and .chimera track routing.
"""

import io
import os
import random

import pytest

from speconsense.chimera import (
    MIN_DIAGNOSTIC_SITES,
    ChimeraHit,
    detect_chimera,
    diagnostic_columns,
    evaluate_triple,
    scan_breakpoint,
)


def _random_seq(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def _substitute(seq, positions):
    """Rotate the base at each position so it differs from the original."""
    rotation = {"A": "C", "C": "G", "G": "T", "T": "A"}
    out = list(seq)
    for pos in positions:
        out[pos] = rotation[out[pos]]
    return "".join(out)


def _make_parents_and_chimera(seed=42, length=400, n_diffs=12, breakpoint=200):
    """Two parents differing at n_diffs spread positions; chimera crosses over
    at `breakpoint` (prefix from A, suffix from B)."""
    rng = random.Random(seed)
    parent_a = _random_seq(rng, length)
    step = length // n_diffs
    diff_positions = [step // 2 + i * step for i in range(n_diffs)]
    parent_b = _substitute(parent_a, diff_positions)
    chimera = parent_a[:breakpoint] + parent_b[breakpoint:]
    return parent_a, parent_b, chimera


class TestScanFunctions:
    def test_diagnostic_columns_excludes_gaps_and_iupac(self):
        a = "ACGT-AR"
        b = "AGGTTAG"
        c = "AGGTTAG"
        # col 1: A/G diag; col 4: gap in A -> excluded; col 6: R in A -> excluded
        diag = diagnostic_columns(a, b, c)
        assert diag == [("C", "G", "G")]

    def test_scan_finds_breakpoint_prefix_a(self):
        # 8 diagnostic sites, C matches A on first 4, B on last 4
        diag = ([("A", "G", "A")] * 4) + ([("A", "G", "G")] * 4)
        m_single, m_chimera, left, right, orient = scan_breakpoint(diag)
        assert m_single == 4
        assert m_chimera == 0
        assert (left, right) == (4, 4)
        assert orient == "A|B"

    def test_scan_finds_breakpoint_prefix_b(self):
        diag = ([("A", "G", "G")] * 4) + ([("A", "G", "A")] * 4)
        m_single, m_chimera, left, right, orient = scan_breakpoint(diag)
        assert m_chimera == 0
        assert orient == "B|A"

    def test_scan_single_parent_candidate(self):
        # C matches A everywhere: chimera model can't improve on single parent
        diag = [("A", "G", "A")] * 8
        m_single, m_chimera, left, right, orient = scan_breakpoint(diag)
        assert m_single == 0
        assert m_chimera == 0

    def test_evaluate_rejects_below_min_sites(self):
        diag_rows = ("A" * (MIN_DIAGNOSTIC_SITES - 1),
                     "G" * (MIN_DIAGNOSTIC_SITES - 1),
                     "A" * (MIN_DIAGNOSTIC_SITES - 1))
        assert evaluate_triple(*diag_rows) is None

    def test_evaluate_rejects_weak_side_support(self):
        # Perfect split but only 2 sites on the left: below MIN_SIDE_SUPPORT
        a = "A" * 12
        b = "G" * 12
        c = "G" * 2 + "A" * 10
        assert evaluate_triple(a, b, c) is None

    def test_evaluate_accepts_clean_chimera(self):
        a = "A" * 12
        b = "G" * 12
        c = "A" * 6 + "G" * 6
        verdict = evaluate_triple(a, b, c)
        assert verdict is not None
        score, k, m_single, m_chimera, left, right, orient = verdict
        assert k == 12
        assert m_single == 6
        assert m_chimera == 0
        assert score == 6
        assert (left, right) == (6, 6)
        assert orient == "A|B"

    def test_evaluate_rejects_noisy_partition(self):
        # Interleaved pattern: no single breakpoint explains it
        a = "A" * 12
        b = "G" * 12
        c = "AG" * 6
        assert evaluate_triple(a, b, c) is None


class TestDetectChimera:
    """SPOA-backed end-to-end detection on realistic sequences."""

    def test_detects_synthetic_chimera(self):
        parent_a, parent_b, chimera = _make_parents_and_chimera()
        hit = detect_chimera(chimera, [(0, parent_a), (1, parent_b)])
        assert hit is not None
        assert {hit.parent_a, hit.parent_b} == {0, 1}
        assert hit.diagnostic_sites == 12
        assert hit.single_mismatches == 6
        assert hit.chimera_mismatches == 0
        assert hit.left_support == 6
        assert hit.right_support == 6
        assert hit.orientation == "A|B"

    def test_no_hit_for_noise_variant(self):
        # Candidate is parent A with a couple of novel substitutions —
        # single-parent model already explains all diagnostic sites.
        parent_a, parent_b, _ = _make_parents_and_chimera()
        noisy = _substitute(parent_a, [3, 7])  # not at diagnostic positions
        hit = detect_chimera(noisy, [(0, parent_a), (1, parent_b)])
        assert hit is None

    def test_no_hit_for_parent_itself(self):
        parent_a, parent_b, _ = _make_parents_and_chimera()
        hit = detect_chimera(parent_a, [(0, parent_a), (1, parent_b)])
        assert hit is None

    def test_single_parent_pool_returns_none(self):
        parent_a, _, chimera = _make_parents_and_chimera()
        assert detect_chimera(chimera, [(0, parent_a)]) is None


class TestCoreGroupDetection:
    def _make_clusterer(self, **kwargs):
        from speconsense.core.clusterer import SpecimenClusterer
        return SpecimenClusterer(**kwargs)

    def test_detect_group_chimeras_stamps_candidate(self):
        parent_a, parent_b, chimera = _make_parents_and_chimera()
        clusterer = self._make_clusterer()
        subclusters = [
            {'read_ids': {f"a{i}" for i in range(40)}},
            {'read_ids': {f"b{i}" for i in range(30)}},
            {'read_ids': {f"c{i}" for i in range(5)}},
        ]
        consensuses = {0: parent_a, 1: parent_b, 2: chimera}
        clusterer._detect_group_chimeras(subclusters, consensuses, [0, 1, 2])

        assert 'chimera' not in subclusters[0]
        assert 'chimera' not in subclusters[1]
        stamped = subclusters[2].get('chimera')
        assert stamped is not None
        assert stamped['diagnostic_sites'] == 12
        assert stamped['chimera_mismatches'] == 0
        assert stamped['score'] == 6
        # prefix parent is A (largest), suffix parent B
        parents = subclusters[2]['_chimera_parent_dicts']
        assert parents[0] is subclusters[0]
        assert parents[1] is subclusters[1]

    def test_abskew_excludes_equal_size_parents(self):
        # Candidate size 30: parent B (also size 30) fails the strictly-
        # larger abundance skew, leaving only one eligible parent — no test
        # possible. (A chimera cannot outnumber the templates it recombined.)
        parent_a, parent_b, chimera = _make_parents_and_chimera()
        clusterer = self._make_clusterer()
        subclusters = [
            {'read_ids': {f"a{i}" for i in range(40)}},
            {'read_ids': {f"b{i}" for i in range(30)}},
            {'read_ids': {f"c{i}" for i in range(30)}},
        ]
        consensuses = {0: parent_a, 1: parent_b, 2: chimera}
        clusterer._detect_group_chimeras(subclusters, consensuses, [0, 1, 2])
        assert 'chimera' not in subclusters[2]

    def test_barely_smaller_candidate_still_tested(self):
        # Candidate size 5 with a size-9 minor-haplotype parent: the real
        # ONT01.25 case that motivated abskew=1.0 over uchime's 2.0.
        parent_a, parent_b, chimera = _make_parents_and_chimera()
        clusterer = self._make_clusterer()
        subclusters = [
            {'read_ids': {f"a{i}" for i in range(777)}},
            {'read_ids': {f"b{i}" for i in range(9)}},
            {'read_ids': {f"c{i}" for i in range(5)}},
        ]
        consensuses = {0: parent_a, 1: parent_b, 2: chimera}
        clusterer._detect_group_chimeras(subclusters, consensuses, [0, 1, 2])
        assert subclusters[2].get('chimera') is not None

    def test_header_emission_includes_chimera_field(self, tmp_path):
        clusterer = self._make_clusterer(sample_name="t",
                                         output_dir=str(tmp_path))
        clusterer.debug_dir = str(tmp_path)
        buf = io.StringIO()
        clusterer.write_cluster_files(
            cluster_name="1.v9", cluster=set(), consensus="ACGT",
            consensus_fasta_handle=buf,
            identity_group_rank=1, identity_variant_rank=9,
            chimera="v1+v3",
        )
        header = buf.getvalue().splitlines()[0]
        assert "chimera=v1+v3" in header
        # chimera= sits before gid=/vid= per house field ordering
        assert header.index("chimera=") < header.index("gid=")


class TestSummarizeRouting:
    def _write_all_fasta(self, tmp_path):
        recs = [
            (">spec-1.v1 size=40 ric=40 err_factor=0.500 gid=1 vid=1", "A" * 60),
            (">spec-1.v2 size=30 ric=30 cer_factor=8.000 err_factor=0.600 gid=1 vid=2", "C" * 60),
            (">spec-1.v9 size=5 ric=5 cer_factor=15.000 err_factor=0.700 chimera=v1+v2 gid=1 vid=9", "G" * 60),
            (">spec-1.v10 size=4 ric=4 cer_factor=12.000 err_factor=9.000 chimera=v1+v2 gid=1 vid=10", "T" * 60),
        ]
        path = tmp_path / "spec-all.fasta"
        with open(path, "w") as f:
            for header, seq in recs:
                f.write(f"{header}\n{seq}\n")
        return str(tmp_path)

    def test_parse_header_extracts_chimera(self):
        from speconsense.summarize.io import parse_consensus_header
        result = parse_consensus_header(
            ">spec-1.v9 size=5 ric=5 chimera=v1+v2 gid=1 vid=9")
        assert result[13] == "v1+v2"
        result = parse_consensus_header(">spec-1.v1 size=40 ric=40 gid=1 vid=1")
        assert result[13] is None

    def test_report_only_default_keeps_pass_track(self, tmp_path):
        from speconsense.summarize.io import load_consensus_sequences
        source = self._write_all_fasta(tmp_path)
        passing, ns, lq, chim = load_consensus_sequences(source, min_ric=3)
        assert chim == []
        flagged = [c for c in passing if c.chimera]
        assert {c.sample_name for c in flagged} == {"spec-1.v9"}
        assert flagged[0].chimera == "v1+v2"
        # the high-err chimera record still routes to lq (lq precedence)
        assert {c.sample_name for c in lq} == {"spec-1.v10"}

    def test_filter_chimeras_routes_to_chimera_track(self, tmp_path):
        from speconsense.summarize.io import load_consensus_sequences
        source = self._write_all_fasta(tmp_path)
        passing, ns, lq, chim = load_consensus_sequences(
            source, min_ric=3, filter_chimeras=True)
        assert {c.sample_name for c in chim} == {"spec-1.v9"}
        # lq still outranks chimera
        assert {c.sample_name for c in lq} == {"spec-1.v10"}
        assert {c.sample_name for c in passing} == {"spec-1.v1", "spec-1.v2"}

    def test_chimera_track_files_written(self, tmp_path):
        from speconsense.summarize.io import write_chimera_variant_files
        from speconsense.summarize.fields import parse_fasta_fields
        from speconsense.types import ConsensusInfo
        rec = ConsensusInfo(
            sample_name="spec-1.v9", cluster_id="1.v9", sequence="G" * 60,
            ric=5, size=5, file_path="x", chimera="v1+v2",
            group_rank=1, variant_rank=9)
        summary = tmp_path / "__Summary__"
        os.makedirs(summary / "variants" / "FASTQ Files")
        write_chimera_variant_files(
            [rec], str(summary), {}, parse_fasta_fields("full"))
        out = summary / "variants" / "spec-1.v9.chimera-RiC5.fasta"
        assert out.exists()
        content = out.read_text()
        assert content.startswith(">spec-1.v9.chimera ")
        assert "chimera=v1+v2" in content


class TestChimeraField:
    def test_format_and_suffix_parsing(self):
        from speconsense.summarize.fields import FASTA_FIELDS
        from speconsense.types import ConsensusInfo
        rec = ConsensusInfo(
            sample_name="spec-1.v9.chimera", cluster_id="1.v9", sequence="G",
            ric=5, size=5, file_path="x", chimera="v1+v2")
        assert FASTA_FIELDS['chimera'].format_value(rec) == "chimera=v1+v2"
        assert FASTA_FIELDS['group'].format_value(rec) == "group=1"
        assert FASTA_FIELDS['variant'].format_value(rec) == "variant=v9"
        unflagged = rec._replace(chimera=None)
        assert FASTA_FIELDS['chimera'].format_value(unflagged) is None

    def test_full_preset_includes_chimera(self):
        from speconsense.summarize.fields import FASTA_FIELD_PRESETS
        assert 'chimera' in FASTA_FIELD_PRESETS['full']
