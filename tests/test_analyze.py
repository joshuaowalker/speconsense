#!/usr/bin/env python3
"""
Tests for speconsense-analyze empirical error analysis functionality.

Tests focus on:
- Parsing consensus headers correctly
- Finding and matching cluster files
- Aligning reads to consensus
- Calculating error statistics
"""

import tempfile
import os
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from speconsense.analyze import (
    parse_consensus_header,
    find_cluster_files,
    load_consensus_sequences,
    align_read_to_consensus,
    analyze_cluster,
    ClusterInfo
)


def test_parse_consensus_header_basic():
    """Test parsing of basic consensus header."""
    header = ">sample-c1 size=200 ric=100"

    sample_name, cluster_num, size, ric, p50_diff, p95_diff = parse_consensus_header(header)

    assert sample_name == "sample-c1"
    assert cluster_num == 1
    assert size == 200
    assert ric == 100
    assert p50_diff is None
    assert p95_diff is None


def test_parse_consensus_header_with_stability():
    """Test parsing of header with stability metrics."""
    header = ">sample-c2 size=150 ric=75 p50diff=0.5 p95diff=1.2"

    sample_name, cluster_num, size, ric, p50_diff, p95_diff = parse_consensus_header(header)

    assert sample_name == "sample-c2"
    assert cluster_num == 2
    assert size == 150
    assert ric == 75
    assert p50_diff == 0.5
    assert p95_diff == 1.2


def test_parse_consensus_header_no_cluster_number():
    """Test that parsing fails for headers without cluster number."""
    header = ">sample size=100 ric=50"

    try:
        parse_consensus_header(header)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "cluster number" in str(e).lower()


def test_align_read_to_consensus_perfect_match():
    """Test alignment of identical read and consensus."""
    consensus = "ACGTACGTACGT"
    read = "ACGTACGTACGT"

    alignment = align_read_to_consensus(read, consensus, "read1")

    assert alignment.read_id == "read1"
    assert alignment.edit_distance == 0
    assert alignment.read_length == 12
    assert alignment.num_insertions == 0
    assert alignment.num_deletions == 0
    assert alignment.num_substitutions == 0


def test_align_read_to_consensus_with_errors():
    """Test alignment with various error types."""
    consensus = "ACGTACGTACGT"
    # Introduce errors: substitution at pos 2, deletion at pos 6, insertion after pos 8
    read = "ACTTACGTAACGT"  # C->T substitution, missing G, extra A

    alignment = align_read_to_consensus(read, consensus, "read2")

    assert alignment.read_id == "read2"
    assert alignment.edit_distance > 0
    assert alignment.read_length == len(read)


def test_find_cluster_files():
    """Test finding cluster FASTQ files in output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create cluster_debug directory structure
        cluster_debug = os.path.join(tmpdir, 'cluster_debug')
        os.makedirs(cluster_debug)

        # Create sample cluster files
        files = [
            'sample1-c1-RiC100-reads.fastq',
            'sample1-c2-RiC50-reads.fastq',
            'sample2-c1-RiC75-reads.fastq',
            'sample1-c1-RiC100-sampled.fastq',  # Sampled version
        ]

        for filename in files:
            Path(os.path.join(cluster_debug, filename)).touch()

        # Test finding reads files
        cluster_files = find_cluster_files(tmpdir, use_sampled=False)

        assert ('sample1', 1) in cluster_files
        assert ('sample1', 2) in cluster_files
        assert ('sample2', 1) in cluster_files
        assert cluster_files[('sample1', 1)].endswith('reads.fastq')

        # Test finding sampled files
        cluster_files_sampled = find_cluster_files(tmpdir, use_sampled=True)

        assert ('sample1', 1) in cluster_files_sampled
        assert cluster_files_sampled[('sample1', 1)].endswith('sampled.fastq')


def test_load_consensus_sequences():
    """Test loading consensus sequences from FASTA files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create cluster_debug directory
        cluster_debug_dir = os.path.join(tmpdir, 'cluster_debug')
        os.makedirs(cluster_debug_dir)

        # Create untrimmed consensus FASTA files
        fasta_content = """>sample1-c1 size=200 ric=100 p50diff=0.0 p95diff=1.0
ACGTACGTACGTACGTACGT
>sample1-c2 size=150 ric=75
TGCATGCATGCATGCA
"""
        fasta_file = os.path.join(cluster_debug_dir, 'sample1-c1-RiC100-untrimmed.fasta')
        with open(fasta_file, 'w') as f:
            f.write(fasta_content.split('\n>')[0] + '\n')

        fasta_file2 = os.path.join(cluster_debug_dir, 'sample1-c2-RiC75-untrimmed.fasta')
        with open(fasta_file2, 'w') as f:
            f.write('>' + fasta_content.split('\n>')[1])

        # Load consensus sequences
        consensus_map = load_consensus_sequences(tmpdir)

        assert ('sample1', 1) in consensus_map
        assert ('sample1', 2) in consensus_map

        c1 = consensus_map[('sample1', 1)]
        assert c1.sample_name == 'sample1-c1'
        assert c1.cluster_num == 1
        assert c1.consensus_seq == 'ACGTACGTACGTACGTACGT'
        assert c1.ric == 100
        assert c1.size == 200
        assert c1.p50_diff == 0.0
        assert c1.p95_diff == 1.0

        c2 = consensus_map[('sample1', 2)]
        assert c2.sample_name == 'sample1-c2'
        assert c2.cluster_num == 2
        assert c2.ric == 75


def test_analyze_cluster_integration():
    """Integration test of full cluster analysis."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create consensus sequence
        consensus_seq = "ACGTACGTACGTACGTACGT"

        # Create FASTQ file with reads (some perfect, some with errors)
        reads_file = os.path.join(tmpdir, 'test-reads.fastq')
        reads = [
            SeqRecord(Seq(consensus_seq), id='read1',
                     letter_annotations={'phred_quality': [30] * len(consensus_seq)}),
            SeqRecord(Seq(consensus_seq), id='read2',
                     letter_annotations={'phred_quality': [30] * len(consensus_seq)}),
            # Read with one substitution
            SeqRecord(Seq("ACTTACGTACGTACGTACGT"), id='read3',
                     letter_annotations={'phred_quality': [30] * 20}),
            # Read with deletion
            SeqRecord(Seq("ACGTACGTACGTACGACGT"), id='read4',
                     letter_annotations={'phred_quality': [30] * 19}),
        ]

        with open(reads_file, 'w') as f:
            SeqIO.write(reads, f, 'fastq')

        # Create ClusterInfo
        cluster_info = ClusterInfo(
            sample_name='test-c1',
            cluster_num=1,
            consensus_seq=consensus_seq,
            reads_file=reads_file,
            ric=4,
            size=4,
            p50_diff=0.0,
            p95_diff=1.0
        )

        # Analyze cluster
        stats, alignments = analyze_cluster(cluster_info)

        # Verify statistics
        assert stats is not None
        assert stats.sample_name == 'test-c1'
        assert stats.cluster_num == 1
        assert stats.num_reads_analyzed == 4
        assert stats.consensus_length == 20
        assert stats.mean_edit_distance > 0  # Should have some errors
        assert stats.mean_error_rate > 0

        # Verify alignments
        assert len(alignments) == 4

        # First two reads should be perfect
        assert alignments[0].edit_distance == 0
        assert alignments[1].edit_distance == 0

        # Third and fourth should have errors
        assert alignments[2].edit_distance > 0
        assert alignments[3].edit_distance > 0


def test_analyze_cluster_missing_file():
    """Test that analyze_cluster handles missing reads file gracefully."""
    cluster_info = ClusterInfo(
        sample_name='test-c1',
        cluster_num=1,
        consensus_seq="ACGT",
        reads_file="/nonexistent/file.fastq",
        ric=10,
        size=10
    )

    stats, alignments = analyze_cluster(cluster_info)

    # Should return None for stats and empty list for alignments
    assert stats is None
    assert alignments == []


def test_cluster_error_rate_calculation():
    """Test that error rates are calculated correctly."""
    with tempfile.TemporaryDirectory() as tmpdir:
        consensus_seq = "A" * 100  # 100bp consensus

        # Create reads with known error counts
        reads_file = os.path.join(tmpdir, 'test-reads.fastq')
        reads = [
            # Perfect read (0% error)
            SeqRecord(Seq("A" * 100), id='read1',
                     letter_annotations={'phred_quality': [30] * 100}),
            # Read with 2 errors (2% error rate)
            SeqRecord(Seq("C" + "A" * 98 + "G"), id='read2',
                     letter_annotations={'phred_quality': [30] * 100}),
            # Read with 5 errors (5% error rate)
            SeqRecord(Seq("CCCCC" + "A" * 95), id='read3',
                     letter_annotations={'phred_quality': [30] * 100}),
        ]

        with open(reads_file, 'w') as f:
            SeqIO.write(reads, f, 'fastq')

        cluster_info = ClusterInfo(
            sample_name='test-c1',
            cluster_num=1,
            consensus_seq=consensus_seq,
            reads_file=reads_file,
            ric=3,
            size=3
        )

        stats, alignments = analyze_cluster(cluster_info)

        # Mean error rate should be (0 + 2 + 5) / 3 / 100 = 0.0233 (2.33%)
        # Due to alignment complexities, check it's in reasonable range
        assert stats.mean_error_rate > 0.01  # At least 1%
        assert stats.mean_error_rate < 0.05  # Less than 5%

        # Verify individual alignments
        assert alignments[0].edit_distance == 0  # Perfect
        assert alignments[1].edit_distance >= 1  # Has errors
        assert alignments[2].edit_distance >= 1  # Has errors


if __name__ == '__main__':
    # Run basic tests
    test_parse_consensus_header_basic()
    test_parse_consensus_header_with_stability()
    test_parse_consensus_header_no_cluster_number()
    test_align_read_to_consensus_perfect_match()
    test_align_read_to_consensus_with_errors()
    test_find_cluster_files()
    test_load_consensus_sequences()
    test_analyze_cluster_integration()
    test_analyze_cluster_missing_file()
    test_cluster_error_rate_calculation()

    print("All tests passed!")
