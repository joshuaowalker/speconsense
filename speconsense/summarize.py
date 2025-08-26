#! /usr/bin/env python3

import os
import re
import glob
import csv
import shutil
import argparse
from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Speconsense output and generate summary.")
    parser.add_argument("--min-ric", type=int, default=3,
                        help="Minimum Reads in Consensus (RiC) threshold (default: 3)")
    parser.add_argument("--source", type=str, default=".",
                        help="Source directory containing Speconsense output (default: current directory)")
    parser.add_argument("--summary-dir", type=str, default="__Summary__",
                        help="Output directory for summary files (default: __Summary__)")
    return parser.parse_args()


def parse_consensus_header(header):
    """Extract information from Speconsense consensus FASTA header."""
    sample_match = re.match(r'>([^;]+);(.+)', header)
    if not sample_match:
        return None, None, None

    sample_name = sample_match.group(1)
    info_string = sample_match.group(2)

    # Extract RiC value
    ric_match = re.search(r'ric=(\d+)', info_string)
    ric = int(ric_match.group(1)) if ric_match else 0

    # Extract size value
    size_match = re.search(r'size=(\d+)', info_string)
    size = int(size_match.group(1)) if size_match else 0

    return sample_name, ric, size


def process_consensus_file(fasta_file, min_ric, summary_folder):
    """Process a single Speconsense consensus FASTA file."""
    with open(fasta_file, 'r') as f:
        lines = f.read().strip().split('\n')
        if len(lines) < 2:
            return None

        header = lines[0]
        sequence = lines[1]

        sample_name, ric, size = parse_consensus_header(header)
        if not sample_name or ric < min_ric:
            return None

        # Copy file to summary folder with the same name
        target_path = os.path.join(summary_folder, os.path.basename(fasta_file))
        shutil.copy(fasta_file, target_path)

        # Check if there's a corresponding FASTQ file in cluster_debug
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        debug_dir = os.path.dirname(os.path.dirname(fasta_file)) if os.path.dirname(fasta_file) != "." else "."
        fastq_file = os.path.join(debug_dir, "cluster_debug", f"{base_name}-reads.fastq")

        if os.path.exists(fastq_file):
            fastq_target = os.path.join(summary_folder, "FASTQ Files", f"{base_name}.fastq")
            shutil.copy(fastq_file, fastq_target)

        return [sample_name, len(sequence), ric, 1]  # 1 is for "multiple" value


def collect_stats(source_folder, min_ric, summary_folder):
    """Collect stats from all Speconsense consensus files."""
    stats = []

    # Find all FASTA files in the immediate source folder matching Speconsense output pattern
    # Note: Using non-recursive glob to avoid subdirectories like cluster_debug
    fasta_pattern = os.path.join(source_folder, "*-RiC*.fasta")
    for fasta_file in sorted(glob.glob(fasta_pattern)):
        # Skip files already in summary folder
        if summary_folder in fasta_file:
            continue

        result = process_consensus_file(fasta_file, min_ric, summary_folder)
        if result:
            stats.append(result)

    return stats


def generate_summary(stats, summary_file):
    """Generate summary statistics in a TSV file."""
    summary_totals = {
        'Total Unique Samples': 0,
        'Total Consensus Sequences': 0,
        'Total Reads in Consensus Sequences': 0,
    }

    # Track unique samples by their base name (before first dash)
    unique_samples = set()

    with open(summary_file, 'w', encoding='utf-8') as f_out:
        writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
        writer.writerow(['Filename', 'Length', 'Reads in Consensus', 'Multiple'])

        for row in stats:
            writer.writerow(row)

            # Extract base sample name (part before first dash)
            base_name = row[0].split('-')[0]
            unique_samples.add(base_name)

            # Update summary totals
            summary_totals['Total Consensus Sequences'] += 1
            summary_totals['Total Reads in Consensus Sequences'] += row[2]  # RiC value

        summary_totals['Total Unique Samples'] = len(unique_samples)

        # Write summary totals
        writer.writerow([])
        writer.writerow(['Total Unique Samples', summary_totals['Total Unique Samples']])
        writer.writerow(['Total Consensus Sequences', summary_totals['Total Consensus Sequences']])
        writer.writerow(['Total Reads in Consensus Sequences', summary_totals['Total Reads in Consensus Sequences']])


def merge_fasta_files(summary_folder, output_file):
    """Concatenate all FASTA files in summary folder into one file."""
    with open(os.path.join(summary_folder, output_file), 'w') as f_out:
        for fasta_file in sorted(glob.glob(f'{summary_folder}/*.fasta')):
            if fasta_file.endswith(output_file):
                continue
            with open(fasta_file, 'r') as f:
                f_out.write(f.read() + "\n")


def main():
    """Main function to process command line arguments and run the summarization."""
    args = parse_arguments()

    source_folder = args.source
    summary_folder = args.summary_dir
    min_ric = args.min_ric

    # Create summary folder and its subdirectories
    if not os.path.exists(summary_folder):
        os.mkdir(summary_folder)

    fastq_folder = os.path.join(summary_folder, 'FASTQ Files')
    if not os.path.exists(fastq_folder):
        os.mkdir(fastq_folder)

    # Process files and collect stats
    print(f"Processing files with minimum RiC threshold of {min_ric}...")
    stats = collect_stats(source_folder, min_ric, summary_folder)

    # Generate summary file
    summary_file = os.path.join(summary_folder, 'summary.txt')
    generate_summary(stats, summary_file)

    # Create merged FASTA file
    merge_fasta_files(summary_folder, 'summary.fasta')

    print(f"Done! Processed {len(stats)} consensus sequences.")
    print(f"Summary files written to {summary_folder}")


if __name__ == "__main__":
    main()