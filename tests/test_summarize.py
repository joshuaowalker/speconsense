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
    """Test merge behavior using complete real specimen file with multiple clusters.

    This uses the ONT01.06-F01--iNat233404862-all.fasta file which contains
    9 clusters. This test demonstrates that merge decisions depend on the full HAC
    group context, not just pairwise comparisons.

    The file contains:
    - c1: main cluster (ric=500) - majority pattern
    - c2: second major cluster (ric=250) - majority pattern
    - c3: contamination (ric=9) - separate group
    - c4: variant (ric=6) - ends with TAG, structural variant
    - c5: variant (ric=6) - structural variant
    - c6: contamination (ric=4) - separate group
    - c7: variant (ric=3) - ends with TAA
    - c8: variant (ric=3) - majority pattern
    - c9: variant (ric=3) - ends with TAA, homopolymer variation from c7

    Expected behavior:
    - c1 + c2 merge into 1.v1 with rawric=500+250
    - c7 + c9 merge into 1.v4 with rawric=3+3 (both end with TAA, differ only by homopolymers)
    - c4 stays separate as 1.v3 (ends with TAG, not compatible with TAA sequences)
    - Total: 7 sequences across 3 HAC groups
    """
    # Use the test data file
    test_file = os.path.join(os.path.dirname(__file__), "data", "ONT01.06-F01--iNat233404862-all.fasta")

    # Skip test if file doesn't exist (e.g., on CI)
    if not os.path.exists(test_file):
        import pytest
        pytest.skip(f"Test file not found: {test_file}")

    # Create temporary directory for output
    temp_dir = tempfile.mkdtemp()
    source_dir = os.path.join(temp_dir, "clusters")
    summary_dir = os.path.join(temp_dir, "__Summary__")
    os.makedirs(source_dir)

    try:
        # Copy the test file to our temp directory
        import shutil as shutil_module
        dest_file = os.path.join(source_dir, "ONT01.06-F01--iNat233404862-all.fasta")
        shutil_module.copy(test_file, dest_file)

        # Run speconsense-summarize with default parameters
        # Disable overlap merge (--min-merge-overlap 0) to test original behavior
        result = subprocess.run(
            [
                "speconsense-summarize",
                "--source", source_dir,
                "--summary-dir", summary_dir,
                "--min-ric", "3",  # Include c4, c7, c8, c9 (all have ric >= 3)
                "--min-merge-overlap", "0"  # Disable overlap merge for this test
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

        # Verify the expected number of output sequences (7 total for this specimen)
        assert len(output_sequences) == 7, \
            f"Expected 7 output sequences, got {len(output_sequences)}"

        # Check if c7 and c9 were merged by examining the output sequences
        # Look for sequences with rawric field indicating a merge
        c7_c9_merged = False
        merged_into = None

        for seq in output_sequences:
            # Check if this sequence has rawric=3+3 (c7 and c9)
            if 'rawric=' in seq.description:
                # Extract rawric values
                rawric_match = re.search(r'rawric=([\d+]+)', seq.description)
                if rawric_match:
                    rawric_str = rawric_match.group(1)
                    ric_values = [int(x) for x in rawric_str.split('+')]
                    # Check if both values are 3 (c7 and c9)
                    if ric_values == [3, 3]:
                        c7_c9_merged = True
                        merged_into = seq.id
                        break

        # Alternative check: Look at .raw files in variants directory
        # If c7 and c9 are in the same variant group, they were merged
        variants_dir = os.path.join(summary_dir, "variants")
        if os.path.exists(variants_dir) and not c7_c9_merged:
            specimen_raw_files = sorted([f for f in os.listdir(variants_dir)
                                         if f.startswith('ONT01.06-F01--iNat233404862') and '.raw' in f])

            # Group raw files by their variant (e.g., "1.v4")
            variant_groups = {}
            for raw_file in specimen_raw_files:
                # Extract variant identifier (e.g., "1.v4" from "...1.v4.raw1...")
                match = re.search(r'-(\d+\.v\d+)\.raw', raw_file)
                if match:
                    variant_id = match.group(1)
                    if variant_id not in variant_groups:
                        variant_groups[variant_id] = []
                    variant_groups[variant_id].append(raw_file)

            # Check each variant group for both c7 and c9 (both end with TAA)
            for variant_id, raw_files in variant_groups.items():
                has_c7_or_c9 = 0

                for raw_file in raw_files:
                    raw_path = os.path.join(variants_dir, raw_file)
                    raw_seqs = list(SeqIO.parse(raw_path, "fasta"))

                    for seq in raw_seqs:
                        seq_str = str(seq.seq)
                        # Both c7 and c9 end with TAA
                        if seq_str.endswith('GACCTCAAATCAGGTAGGACTACCCGCTGAACTTAA'):
                            has_c7_or_c9 += 1

                if has_c7_or_c9 >= 2:  # Found at least 2 sequences ending with TAA
                    c7_c9_merged = True
                    merged_into = variant_id
                    break

        # Key assertion: c7 and c9 SHOULD be merged when in HAC group context
        # This is because the multi-sequence alignment reveals their differences
        # are homopolymer variations (both end with TAA)
        assert c7_c9_merged, \
            f"c7 (ric=3) and c9 (ric=3) should be merged in HAC group context, " \
            f"but they were not merged"

        # Verify they merged into the expected variant group (1.v4)
        assert merged_into == '1.v4' or merged_into == 'ONT01.06-F01--iNat233404862-1.v4', \
            f"Expected c7 and c9 to merge into variant 1.v4, but merged into: {merged_into}"

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
