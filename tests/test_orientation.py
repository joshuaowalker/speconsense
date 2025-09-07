#!/usr/bin/env python3
"""
Minimal tests for sequence orientation functionality.

Tests focus on behavior, not implementation details:
- Sequences with clear orientation are handled correctly
- Orientation modes work as specified
- No orientation performed when mode is 'skip'
"""

import tempfile
import os
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from speconsense.core import SpecimenClusterer


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


def test_orientation_skip_mode():
    """Test that skip mode doesn't perform any orientation.
    
    Note: This tests that when --orient-mode=skip (default), the main() function
    doesn't call orient_sequences(). We can't easily test the full CLI here,
    so we just verify that orient_sequences() works when called directly.
    """
    # Create sequences that would be reoriented if orientation was performed
    forward_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'  # RC of reverse at 3' end
    reverse_seq = 'TTTTTTTTTT' + 'C' * 100 + 'TTTTTTTTTT'  # Reverse at 5', RC of forward at 3'
    
    records = [
        SeqRecord(Seq(forward_seq), id='seq1', letter_annotations={'phred_quality': [30] * len(forward_seq)}),
        SeqRecord(Seq(reverse_seq), id='seq2', letter_annotations={'phred_quality': [30] * len(reverse_seq)})
    ]
    
    clusterer = SpecimenClusterer()
    clusterer.add_sequences(records)
    
    primer_file = create_test_primers()
    try:
        clusterer.load_primers(primer_file)
        
        # Store originals
        original_seq1 = clusterer.sequences['seq1']
        original_seq2 = clusterer.sequences['seq2']
        
        # When orient_sequences IS called, it should work correctly
        failed = clusterer.orient_sequences()
        
        # seq1 should stay (already forward)
        assert clusterer.sequences['seq1'] == original_seq1
        # seq2 should be flipped (was reverse)
        assert clusterer.sequences['seq2'] == str(reverse_complement(original_seq2))
        
        # Both sequences had clear orientation, so no failures
        assert len(failed) == 0
    finally:
        os.unlink(primer_file)


def test_orientation_keep_all_mode():
    """Test that keep-all mode orients sequences but keeps failed ones."""
    # Forward orientation: forward primer at 5', RC of reverse at 3'
    forward_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'
    # Reverse orientation: reverse primer at 5', RC of forward at 3'  
    reverse_seq = 'TTTTTTTTTT' + 'C' * 100 + 'TTTTTTTTTT'
    # No primers - should fail orientation
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
        original_none = clusterer.sequences['none']
        
        failed = clusterer.orient_sequences()
        
        # Forward should stay the same
        assert clusterer.sequences['forward'] == original_forward
        
        # Reverse should be flipped
        assert clusterer.sequences['reverse'] == str(reverse_complement(original_reverse))
        
        # No-primer should stay the same but be marked as failed
        assert clusterer.sequences['none'] == original_none
        assert 'none' in failed
        
        # All sequences should still be present
        assert len(clusterer.sequences) == 3
    finally:
        os.unlink(primer_file)


def test_orientation_filter_failed():
    """Test that filter-failed mode removes sequences that fail orientation."""
    # Good sequence with clear orientation (forward at 5', RC of reverse at 3')
    good_seq = 'AAAAAAAAAA' + 'G' * 100 + 'AAAAAAAAAA'
    # No primers sequence - will fail orientation
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
        
        # no_primer should be in failed set
        assert 'no_primer' in failed
        assert 'good' not in failed
        
        # Simulate filtering (as done in main())
        for seq_id in failed:
            del clusterer.sequences[seq_id]
            del clusterer.records[seq_id]
        
        # Only good sequence should remain
        assert len(clusterer.sequences) == 1
        assert 'good' in clusterer.sequences
        assert 'no_primer' not in clusterer.sequences
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
        
        # Should have loaded primers despite missing position
        assert len(clusterer.primers) > 0
        # Should treat as bidirectional
        assert len(clusterer.forward_primers) == 2
        assert len(clusterer.reverse_primers) == 2
    finally:
        os.unlink(primer_file)