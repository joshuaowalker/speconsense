#!/usr/bin/env python3
"""
Tests for speconsense-synth synthetic read generator.
"""

import os
import tempfile
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from speconsense.synth import (
    parse_ratios,
    error_rate_to_phred,
    introduce_errors,
    generate_reads
)


def test_parse_ratios():
    """Test ratio parsing and normalization."""
    # Equal distribution by default
    assert parse_ratios("", 3) == [1/3, 1/3, 1/3]
    
    # Normalize to sum to 1
    assert parse_ratios("50,30,20", 3) == [0.5, 0.3, 0.2]
    assert parse_ratios("2,2,1", 3) == [0.4, 0.4, 0.2]
    
    # Single sequence
    assert parse_ratios("100", 1) == [1.0]
    
    # Test error cases
    try:
        parse_ratios("50,50", 3)  # Wrong number
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_error_rate_to_phred():
    """Test Phred score calculation."""
    # Common error rates
    assert error_rate_to_phred(0.1) == 10    # Q10 = 90% accuracy
    assert error_rate_to_phred(0.01) == 20   # Q20 = 99% accuracy
    assert error_rate_to_phred(0.001) == 30  # Q30 = 99.9% accuracy
    
    # Edge cases
    assert error_rate_to_phred(0) == 40      # Cap at Q40
    assert error_rate_to_phred(1) == 0       # Q0 = 0% accuracy
    assert error_rate_to_phred(0.5) == 3     # Q3 = 50% accuracy


def test_introduce_errors_no_errors():
    """Test that no errors are introduced at 0% rate."""
    rng = random.Random(42)
    sequence = "ACGTACGTACGT"
    
    result = introduce_errors(sequence, 0.0, rng)
    assert result == sequence


def test_introduce_errors_deterministic():
    """Test that errors are reproducible with same seed."""
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    
    rng1 = random.Random(42)
    result1 = introduce_errors(sequence, 0.1, rng1)
    
    rng2 = random.Random(42)
    result2 = introduce_errors(sequence, 0.1, rng2)
    
    assert result1 == result2


def test_introduce_errors_types():
    """Test that all error types occur."""
    sequence = "A" * 1000  # Long sequence of A's
    rng = random.Random(42)
    
    result = introduce_errors(sequence, 0.1, rng)
    
    # Should have substitutions (non-A bases)
    assert any(b != 'A' for b in result)
    
    # Length should differ due to indels
    assert len(result) != len(sequence)


def test_generate_reads_basic():
    """Test basic read generation."""
    sequences = [
        SeqRecord(Seq("ACGTACGTACGT"), id="seq1"),
        SeqRecord(Seq("TTTTGGGGCCCC"), id="seq2")
    ]
    
    reads = generate_reads(
        sequences=sequences,
        num_reads=100,
        error_rate=0.1,
        ratios=[0.5, 0.5],
        seed=42
    )
    
    assert len(reads) == 100
    
    # Check that reads have quality scores
    for read in reads:
        assert 'phred_quality' in read.letter_annotations
        assert len(read.letter_annotations['phred_quality']) == len(read.seq)
    
    # Check provenance in headers
    seq1_reads = [r for r in reads if 'seq1' in r.id]
    seq2_reads = [r for r in reads if 'seq2' in r.id]
    
    # Should be roughly 50/50 split
    assert 40 <= len(seq1_reads) <= 60
    assert 40 <= len(seq2_reads) <= 60


def test_generate_reads_single_sequence():
    """Test generation from single sequence."""
    sequences = [
        SeqRecord(Seq("ACGTACGTACGTACGT"), id="only_seq")
    ]
    
    reads = generate_reads(
        sequences=sequences,
        num_reads=10,
        error_rate=0.05,
        ratios=[1.0],
        seed=42
    )
    
    assert len(reads) == 10
    assert all('only_seq' in r.id for r in reads)


def test_generate_reads_uneven_ratios():
    """Test generation with uneven ratios."""
    sequences = [
        SeqRecord(Seq("AAAAAAAAAA"), id="major"),
        SeqRecord(Seq("TTTTTTTTTT"), id="minor")
    ]
    
    reads = generate_reads(
        sequences=sequences,
        num_reads=100,
        error_rate=0.1,
        ratios=[0.9, 0.1],  # 90/10 split
        seed=42
    )
    
    major_reads = [r for r in reads if 'major' in r.id]
    minor_reads = [r for r in reads if 'minor' in r.id]
    
    # Check approximate ratios
    assert 85 <= len(major_reads) <= 95
    assert 5 <= len(minor_reads) <= 15


def test_output_format():
    """Test that output can be written as valid FASTQ."""
    sequences = [SeqRecord(Seq("ACGTACGTACGT"), id="test")]
    
    reads = generate_reads(
        sequences=sequences,
        num_reads=5,
        error_rate=0.1,
        ratios=[1.0],
        seed=42
    )
    
    # Write and read back
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        SeqIO.write(reads, f, 'fastq')
        temp_file = f.name
    
    try:
        # Should be readable as FASTQ
        read_back = list(SeqIO.parse(temp_file, 'fastq'))
        assert len(read_back) == 5
        
        # Check quality scores preserved
        for original, read in zip(reads, read_back):
            assert original.letter_annotations['phred_quality'] == read.letter_annotations['phred_quality']
    finally:
        os.unlink(temp_file)


if __name__ == "__main__":
    # Run tests
    test_parse_ratios()
    test_error_rate_to_phred()
    test_introduce_errors_no_errors()
    test_introduce_errors_deterministic()
    test_introduce_errors_types()
    test_generate_reads_basic()
    test_generate_reads_single_sequence()
    test_generate_reads_uneven_ratios()
    test_output_format()
    
    print("All tests passed!")