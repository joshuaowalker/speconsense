#!/usr/bin/env python3
"""
Validate consistency of speconsense-summarize output.

Checks:
1. RiC values in FASTA headers match actual read counts in FASTQ files
2. For merged consensuses, merged_ric values sum to total ric
3. For .raw files, parent consensus has matching merged_ric components
4. All expected .raw files exist for merged consensuses
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_fasta_header(header_line):
    """Parse FASTA header to extract metadata."""
    # Format: >sample_name size=N ric=N [merged_ric=N+N+...] [snp=N] [primers=...]
    match = re.match(r'>(\S+)\s+(.+)', header_line)
    if not match:
        return None

    sample_name = match.group(1)
    metadata_str = match.group(2)

    metadata = {'sample_name': sample_name}

    # Extract size
    size_match = re.search(r'size=(\d+)', metadata_str)
    if size_match:
        metadata['size'] = int(size_match.group(1))

    # Extract ric
    ric_match = re.search(r'ric=(\d+)', metadata_str)
    if ric_match:
        metadata['ric'] = int(ric_match.group(1))

    # Extract merged_ric if present
    merged_ric_match = re.search(r'merged_ric=([\d+]+)', metadata_str)
    if merged_ric_match:
        ric_values = [int(x) for x in merged_ric_match.group(1).split('+')]
        metadata['merged_ric'] = ric_values

    # Extract snp if present
    snp_match = re.search(r'snp=(\d+)', metadata_str)
    if snp_match:
        metadata['snp'] = int(snp_match.group(1))

    return metadata


def count_fastq_reads(fastq_path):
    """Count reads in FASTQ file (lines / 4)."""
    if not os.path.exists(fastq_path):
        return None

    try:
        with open(fastq_path, 'r') as f:
            line_count = sum(1 for _ in f)
        return line_count // 4
    except Exception as e:
        print(f"Error reading {fastq_path}: {e}")
        return None


def validate_dataset(summary_dir):
    """Run all validation checks on the dataset."""

    summary_fasta = os.path.join(summary_dir, 'summary.fasta')
    fastq_dir = os.path.join(summary_dir, 'FASTQ Files')

    if not os.path.exists(summary_fasta):
        print(f"Error: {summary_fasta} not found")
        return False

    # Parse all sequences from summary.fasta
    sequences = []
    with open(summary_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                metadata = parse_fasta_header(line.strip())
                if metadata:
                    sequences.append(metadata)

    print(f"Loaded {len(sequences)} sequences from summary.fasta")
    print()

    # Categorize sequences
    final_consensus = [s for s in sequences if '.raw' not in s['sample_name']]
    raw_files = [s for s in sequences if '.raw' in s['sample_name']]
    merged_consensus = [s for s in final_consensus if 'merged_ric' in s]

    print(f"  {len(final_consensus)} final consensus sequences")
    print(f"  {len(merged_consensus)} merged consensus sequences")
    print(f"  {len(raw_files)} .raw file sequences")
    print()

    # Group .raw files by parent
    raw_by_parent = defaultdict(list)
    for raw in raw_files:
        # Extract parent name (remove .rawN suffix)
        parent_name = re.sub(r'\.raw\d+$', '', raw['sample_name'])
        raw_by_parent[parent_name].append(raw)

    # Validation checks
    errors = []
    warnings = []

    print("=" * 70)
    print("CHECK 1: RiC values in FASTA headers match FASTQ read counts")
    print("=" * 70)

    checked = 0
    mismatches = []
    missing_fastq = []

    for seq in sequences:
        sample_name = seq['sample_name']
        ric = seq['ric']

        fastq_path = os.path.join(fastq_dir, f"{sample_name}-RiC{ric}.fastq")
        actual_reads = count_fastq_reads(fastq_path)

        if actual_reads is None:
            missing_fastq.append(sample_name)
        elif actual_reads != ric:
            mismatches.append((sample_name, ric, actual_reads))

        checked += 1

    if missing_fastq:
        print(f"⚠️  WARNING: {len(missing_fastq)} FASTQ files not found")
        for name in missing_fastq[:10]:
            print(f"    Missing: {name}")
        if len(missing_fastq) > 10:
            print(f"    ... and {len(missing_fastq) - 10} more")
        warnings.append(f"{len(missing_fastq)} missing FASTQ files")

    if mismatches:
        print(f"❌ ERROR: {len(mismatches)} RiC/read count mismatches found!")
        for name, ric, actual in mismatches[:10]:
            print(f"    {name}: header says ric={ric}, FASTQ has {actual} reads")
        if len(mismatches) > 10:
            print(f"    ... and {len(mismatches) - 10} more")
        errors.append(f"{len(mismatches)} RiC mismatches")
    else:
        print(f"✅ PASS: All {checked} sequences have matching RiC values and read counts")

    print()

    print("=" * 70)
    print("CHECK 2: merged_ric values sum to total ric")
    print("=" * 70)

    merged_sum_errors = []
    for seq in merged_consensus:
        ric = seq['ric']
        merged_ric = seq['merged_ric']
        merged_sum = sum(merged_ric)

        if merged_sum != ric:
            merged_sum_errors.append((seq['sample_name'], ric, merged_ric, merged_sum))

    if merged_sum_errors:
        print(f"❌ ERROR: {len(merged_sum_errors)} merged_ric sum mismatches!")
        for name, ric, merged_ric, merged_sum in merged_sum_errors[:10]:
            print(f"    {name}: ric={ric}, merged_ric={'+'.join(map(str, merged_ric))}={merged_sum}")
        if len(merged_sum_errors) > 10:
            print(f"    ... and {len(merged_sum_errors) - 10} more")
        errors.append(f"{len(merged_sum_errors)} merged_ric sum mismatches")
    else:
        print(f"✅ PASS: All {len(merged_consensus)} merged consensuses have correct merged_ric sums")

    print()

    print("=" * 70)
    print("CHECK 3: .raw files match merged_ric components")
    print("=" * 70)

    raw_component_errors = []
    for seq in merged_consensus:
        parent_name = seq['sample_name']
        merged_ric = sorted(seq['merged_ric'], reverse=True)

        # Get .raw files for this parent
        raw_files_for_parent = sorted(raw_by_parent.get(parent_name, []),
                                      key=lambda x: x['ric'], reverse=True)

        if len(raw_files_for_parent) != len(merged_ric):
            raw_component_errors.append(
                f"{parent_name}: has {len(merged_ric)} merged_ric components but {len(raw_files_for_parent)} .raw files"
            )
            continue

        # Check each .raw file's ric matches merged_ric component
        for raw_file, expected_ric in zip(raw_files_for_parent, merged_ric):
            if raw_file['ric'] != expected_ric:
                raw_component_errors.append(
                    f"{raw_file['sample_name']}: ric={raw_file['ric']}, expected {expected_ric} from parent merged_ric"
                )

    if raw_component_errors:
        print(f"❌ ERROR: {len(raw_component_errors)} .raw file component mismatches!")
        for error in raw_component_errors[:10]:
            print(f"    {error}")
        if len(raw_component_errors) > 10:
            print(f"    ... and {len(raw_component_errors) - 10} more")
        errors.append(f"{len(raw_component_errors)} .raw component mismatches")
    else:
        print(f"✅ PASS: All {len(merged_consensus)} merged consensuses have correct .raw components")

    print()

    print("=" * 70)
    print("CHECK 4: All merged consensuses have .raw files")
    print("=" * 70)

    missing_raw_files = []
    for seq in merged_consensus:
        parent_name = seq['sample_name']
        merged_ric = seq['merged_ric']

        if len(merged_ric) <= 1:
            # Should not have merged_ric if only 1 component
            missing_raw_files.append(f"{parent_name}: merged_ric has only {len(merged_ric)} component(s)")
            continue

        raw_files_for_parent = raw_by_parent.get(parent_name, [])
        if len(raw_files_for_parent) == 0:
            missing_raw_files.append(f"{parent_name}: has merged_ric={merged_ric} but no .raw files")

    if missing_raw_files:
        print(f"❌ ERROR: {len(missing_raw_files)} merged consensuses missing .raw files!")
        for error in missing_raw_files[:10]:
            print(f"    {error}")
        if len(missing_raw_files) > 10:
            print(f"    ... and {len(missing_raw_files) - 10} more")
        errors.append(f"{len(missing_raw_files)} missing .raw files")
    else:
        print(f"✅ PASS: All {len(merged_consensus)} merged consensuses have .raw files")

    print()

    print("=" * 70)
    print("CHECK 5: No .raw files for non-merged consensuses")
    print("=" * 70)

    unexpected_raw_files = []
    non_merged = [s for s in final_consensus if 'merged_ric' not in s]
    for seq in non_merged:
        name = seq['sample_name']
        if name in raw_by_parent:
            unexpected_raw_files.append(f"{name}: not merged but has {len(raw_by_parent[name])} .raw files")

    if unexpected_raw_files:
        print(f"❌ ERROR: {len(unexpected_raw_files)} non-merged consensuses have .raw files!")
        for error in unexpected_raw_files[:10]:
            print(f"    {error}")
        if len(unexpected_raw_files) > 10:
            print(f"    ... and {len(unexpected_raw_files) - 10} more")
        errors.append(f"{len(unexpected_raw_files)} unexpected .raw files")
    else:
        print(f"✅ PASS: No .raw files for {len(non_merged)} non-merged consensuses")

    print()

    # Summary
    print("=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    if errors:
        print(f"❌ FAILED: {len(errors)} error(s) found:")
        for error in errors:
            print(f"    - {error}")
        return False
    elif warnings:
        print(f"⚠️  WARNINGS: {len(warnings)} warning(s):")
        for warning in warnings:
            print(f"    - {warning}")
        print()
        print("✅ All validation checks passed (with warnings)")
        return True
    else:
        print("✅ All validation checks passed!")
        return True


if __name__ == '__main__':
    if len(sys.argv) > 1:
        summary_dir = sys.argv[1]
    else:
        summary_dir = '__Summary__'

    success = validate_dataset(summary_dir)
    sys.exit(0 if success else 1)
