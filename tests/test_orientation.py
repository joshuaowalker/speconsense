#!/usr/bin/env python3
"""
Tests for sequence orientation functionality.

Tests cover:
- Primer-based orientation (--orient-mode=primer)
- pyitsx HMM-based orientation (--orient-mode=pyitsx, if available)
- Primer loading with and without position annotations
"""

import tempfile
import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from speconsense.core import SpecimenClusterer

try:
    import pyitsx
    HAS_PYITSX = True
except ImportError:
    HAS_PYITSX = False


def create_test_primers():
    """Create a simple primer file with forward and reverse primers."""
    primers_content = """>Forward1   position=forward
AAAAAAAAAA
>Reverse1   position=reverse
TTTTTTTTTT
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(primers_content)
        return f.name


def test_primer_orientation_basic():
    """Test that primer-based orientation correctly identifies and flips reads."""
    forward_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'
    reverse_seq = 'TTTTTTTTTT' + 'C' * 100 + 'TTTTTTTTTT'

    records = [
        SeqRecord(Seq(forward_seq), id='seq1', letter_annotations={'phred_quality': [30] * len(forward_seq)}),
        SeqRecord(Seq(reverse_seq), id='seq2', letter_annotations={'phred_quality': [30] * len(reverse_seq)})
    ]

    clusterer = SpecimenClusterer()
    clusterer.add_sequences(records)

    primer_file = create_test_primers()
    try:
        clusterer.load_primers(primer_file)

        original_seq1 = clusterer.sequences['seq1']
        original_seq2 = clusterer.sequences['seq2']

        failed = clusterer.orient_sequences()

        assert clusterer.sequences['seq1'] == original_seq1
        assert clusterer.sequences['seq2'] == str(reverse_complement(original_seq2))
        assert len(failed) == 0
    finally:
        os.unlink(primer_file)


def test_primer_orientation_with_failures():
    """Test that primer-based orientation returns failed read IDs."""
    forward_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'
    reverse_seq = 'TTTTTTTTTT' + 'C' * 100 + 'TTTTTTTTTT'
    no_primer_seq = 'G' * 120

    records = [
        SeqRecord(Seq(forward_seq), id='forward', letter_annotations={'phred_quality': [30] * len(forward_seq)}),
        SeqRecord(Seq(reverse_seq), id='reverse', letter_annotations={'phred_quality': [30] * len(reverse_seq)}),
        SeqRecord(Seq(no_primer_seq), id='none', letter_annotations={'phred_quality': [30] * len(no_primer_seq)})
    ]

    clusterer = SpecimenClusterer()
    clusterer.add_sequences(records)

    primer_file = create_test_primers()
    try:
        clusterer.load_primers(primer_file)

        original_forward = clusterer.sequences['forward']
        original_reverse = clusterer.sequences['reverse']

        failed = clusterer.orient_sequences()

        assert clusterer.sequences['forward'] == original_forward
        assert clusterer.sequences['reverse'] == str(reverse_complement(original_reverse))
        assert 'none' in failed
        assert len(clusterer.sequences) == 3
    finally:
        os.unlink(primer_file)


def test_primer_orientation_filter_discards():
    """Test that failed reads can be discarded after primer orientation."""
    good_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'
    no_primer_seq = 'G' * 120

    records = [
        SeqRecord(Seq(good_seq), id='good', letter_annotations={'phred_quality': [30] * len(good_seq)}),
        SeqRecord(Seq(no_primer_seq), id='no_primer', letter_annotations={'phred_quality': [30] * len(no_primer_seq)})
    ]

    clusterer = SpecimenClusterer()
    clusterer.add_sequences(records)

    primer_file = create_test_primers()
    try:
        clusterer.load_primers(primer_file)

        failed = clusterer.orient_sequences()

        assert 'no_primer' in failed
        assert 'good' not in failed

        for seq_id in failed:
            del clusterer.sequences[seq_id]
            del clusterer.records[seq_id]

        assert len(clusterer.sequences) == 1
        assert 'good' in clusterer.sequences
    finally:
        os.unlink(primer_file)


def test_primers_without_position():
    """Test that primers without position annotation are handled gracefully."""
    primers_content = """>Primer1
AAAAAAAAAA
>Primer2
TTTTTTTTTT
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(primers_content)
        primer_file = f.name

    try:
        clusterer = SpecimenClusterer()
        clusterer.load_primers(primer_file)

        assert len(clusterer.primers) > 0
        assert len(clusterer.forward_primers) == 2
        assert len(clusterer.reverse_primers) == 2
    finally:
        os.unlink(primer_file)


def test_primer_orientation_subset():
    """Test that orient_sequences can operate on a subset of sequence IDs."""
    forward_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'
    reverse_seq = 'TTTTTTTTTT' + 'C' * 100 + 'TTTTTTTTTT'
    other_seq = 'CCCC' * 30

    records = [
        SeqRecord(Seq(forward_seq), id='seq1', letter_annotations={'phred_quality': [30] * len(forward_seq)}),
        SeqRecord(Seq(reverse_seq), id='seq2', letter_annotations={'phred_quality': [30] * len(reverse_seq)}),
        SeqRecord(Seq(other_seq), id='seq3', letter_annotations={'phred_quality': [30] * len(other_seq)}),
    ]

    clusterer = SpecimenClusterer()
    clusterer.add_sequences(records)

    primer_file = create_test_primers()
    try:
        clusterer.load_primers(primer_file)

        original_seq1 = clusterer.sequences['seq1']
        original_seq2 = clusterer.sequences['seq2']
        original_seq3 = clusterer.sequences['seq3']

        # Only orient seq2 and seq3 — seq1 should be untouched
        failed = clusterer.orient_sequences(seq_ids={'seq2', 'seq3'})

        assert clusterer.sequences['seq1'] == original_seq1
        assert clusterer.sequences['seq2'] == str(reverse_complement(original_seq2))
        assert clusterer.sequences['seq3'] == original_seq3
        assert 'seq3' in failed
        assert 'seq2' not in failed
        assert 'seq1' not in failed
    finally:
        os.unlink(primer_file)


# --- pyitsx tests ---

@pytest.mark.skipif(not HAS_PYITSX, reason="pyitsx not installed")
def test_pyitsx_orientation_import_error():
    """Test clear error when pyitsx is not available."""
    from unittest.mock import patch

    clusterer = SpecimenClusterer()
    clusterer.sequences = {'seq1': 'ATCG'}
    clusterer.records = {}

    with patch.dict('sys.modules', {'pyitsx': None}):
        with pytest.raises(ImportError, match="pyitsx is required"):
            clusterer.orient_sequences_pyitsx()


@pytest.mark.skipif(not HAS_PYITSX, reason="pyitsx not installed")
def test_pyitsx_orientation_returns_failed_and_chimeric():
    """Test that pyitsx orientation returns separate failed and chimeric sets."""
    short_seq = 'ATCG' * 3

    records = [
        SeqRecord(Seq(short_seq), id='too_short', letter_annotations={'phred_quality': [30] * len(short_seq)}),
    ]

    clusterer = SpecimenClusterer()
    clusterer.add_sequences(records)

    failed_ids, chimeric_ids = clusterer.orient_sequences_pyitsx(organism="F")

    assert isinstance(failed_ids, set)
    assert isinstance(chimeric_ids, set)
    assert 'too_short' in failed_ids
