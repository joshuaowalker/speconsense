#!/usr/bin/env python3
"""
Tests for speconsense-summarize functionality.

Tests focus on merge behavior with different sequence variants.
"""

import tempfile
import os
import shutil
import subprocess
import re
from Bio import SeqIO


def test_merge_behavior_with_full_hac_context():
    """Test merge behavior using real specimen file with MSA-based merging.

    This uses the ONT01.06-F01--iNat233404862-all.fasta file which contains
    2 clusters (c1 and c2). This test demonstrates that sequences with only
    SNP differences (no structural indels) merge with MSA-based consensus.

    The file contains:
    - c1: main cluster (ric=500, size=671)
    - c2: second cluster (ric=250, size=250)

    Both sequences are nearly identical with only a few SNP differences,
    so they should merge into a single IUPAC consensus sequence.

    Expected behavior: c1 and c2 merge into a single variant with rawric=500+250
    """
    # Use the real file from the data directory
    real_file = "/Users/josh/mm/data/ont98/demux20251106/clusters/ONT01.06-F01--iNat233404862-all.fasta"

    # Skip test if file doesn't exist (e.g., on CI)
    if not os.path.exists(real_file):
        import pytest
        pytest.skip(f"Test file not found: {real_file}")

    # Create temporary directory for output
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # Copy the real file to our temp directory
        import shutil as shutil_module
        dest_file = os.path.join(source_dir, "ONT01.06-F01--iNat233404862-all.fasta")
        shutil_module.copy(real_file, dest_file)

        # Run speconsense-summarize with default parameters
        result = subprocess.run(
            [
                "speconsense-summarize",
                "--source", source_dir,
                "--summary-dir", summary_dir,
                "--min-ric", "3"  # Include c4, c7, c8, c9 (all have ric >= 3)
            ],
            capture_output=True,
            text=True
        )

        # Check that the command succeeded
        assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

        # Read the main output FASTA file
        output_fasta = os.path.join(summary_dir, "summary.fasta")
        assert os.path.exists(output_fasta), \
            f"Expected output file not found: {output_fasta}"

        # Read all sequences from the output
        output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

        # Print diagnostic information
        print(f"\nOutput sequences: {len(output_sequences)}")
        for seq in output_sequences:
            print(f"  {seq.id}: {seq.description}")

        # c1 and c2 should merge into a single sequence since they only differ by SNPs
        assert len(output_sequences) == 1, \
            f"Expected 1 merged sequence in output, got {len(output_sequences)}"

        # Verify the merged sequence has the expected rawric (500+250 or 250+500)
        merged_seq = output_sequences[0]
        assert 'rawric=' in merged_seq.description, \
            f"Expected merged sequence to have rawric field"

        # Extract and verify rawric values
        rawric_match = re.search(r'rawric=([\d+]+)', merged_seq.description)
        assert rawric_match, "Could not find rawric in merged sequence"

        rawric_str = rawric_match.group(1)
        ric_values = set(int(x) for x in rawric_str.split('+'))

        # Should contain both original RiC values (500 and 250)
        assert ric_values == {500, 250}, \
            f"Expected rawric to contain 500 and 250, got: {ric_values}"

        # Verify SNP count is reported (sequences differ by a few SNPs)
        assert 'snp=' in merged_seq.description, \
            f"Expected merged sequence to have snp field indicating IUPAC consensus"

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)


def test_merge_with_homopolymer_only_differences():
    """Test that sequences differing only in homopolymer lengths DO merge.

    This test verifies that sequences with identical structure but different
    homopolymer lengths will merge with the homopolymer-aware algorithm.
    """
    # Create temporary directory structure
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # Create two sequences that differ only in homopolymer length
        # Base sequence with A homopolymer of length 5
        seq1 = "ATCGAAAAATCGATCGATCGATCG"
        # Same sequence with A homopolymer of length 8
        seq2 = "ATCGAAAAAAATCGATCGATCGATCG"

        fasta_content = f""">test-seq1 size=10 ric=10 primers=test
{seq1}
>test-seq2 size=8 ric=8 primers=test
{seq2}
"""

        fasta_file = os.path.join(source_dir, "test-homopoly-all.fasta")
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)

        # Run speconsense-summarize with default parameters
        result = subprocess.run(
            [
                "speconsense-summarize",
                "--source", source_dir,
                "--summary-dir", summary_dir,
                "--min-ric", "3"
            ],
            capture_output=True,
            text=True
        )

        # Check that the command succeeded
        assert result.returncode == 0, f"speconsense-summarize failed: {result.stderr}"

        # Read the main output FASTA file (summary.fasta combines all specimens)
        output_fasta = os.path.join(summary_dir, "summary.fasta")
        assert os.path.exists(output_fasta), \
            f"Expected output file not found: {output_fasta}"

        # Count sequences in output
        output_sequences = list(SeqIO.parse(output_fasta, "fasta"))

        # Should have 1 sequence (merged due to homopolymer equivalence)
        assert len(output_sequences) == 1, \
            f"Expected 1 merged sequence, but got {len(output_sequences)}"

        # The merged sequence should have the combined size
        # Check the header for size information
        header = output_sequences[0].description
        assert "size=" in header, "Expected size field in output header"

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)
