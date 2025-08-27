# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Speconsense is a Python tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads. It's designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

## Key Commands

### Running the main tool
```bash
# Basic usage
speconsense input.fastq

# With common options
speconsense input.fastq --algorithm greedy --min-size 10 --primers primers.fasta --variant-merge-threshold 2

# Test with sample data
speconsense input.fastq --min-identity 0.85 --max-sample-size 500 --presample 1000
```

### Post-processing
```bash
# Generate summary from output files
speconsense-summarize

# With custom parameters
speconsense-summarize --min-ric 5 --source /path/to/output --summary-dir MyResults
```

### Installation
```bash
# Install from GitHub (recommended)
pip install git+https://github.com/joshuaowalker/speconsense.git

# Or for development
git clone https://github.com/joshuaowalker/speconsense.git
cd speconsense
pip install -e .

# External tools (must be in PATH):
# - SPOA: conda install bioconda::spoa
# - MCL: conda install bioconda::mcl (optional but recommended)
```

## Code Architecture

### Main Components

**speconsense/core.py** - Core clustering and consensus tool (formerly speconsense.py):
- `SpecimenClusterer` class: Main orchestrator for the clustering pipeline
- `AmbiguityMode` class: Controls handling of ambiguous nucleotides using IUPAC codes
- Two clustering algorithms: graph-based MCL (default) and greedy clustering
- Consensus generation using external SPOA tool
- Cluster merging based on consensus sequence similarity

**speconsense/summarize.py** - Post-processing utility (formerly summarize_speconsense.py):
- Parses Speconsense FASTA output files
- Generates summary reports and filtered datasets
- Handles RiC (Reads in Consensus) filtering

**speconsense/cli.py** - Command-line interface entry point for main speconsense command

### Key Processing Steps

1. **Input processing**: Read FASTQ files, optional presampling
2. **Clustering**: Either MCL graph-based or greedy algorithm
3. **Consensus generation**: Uses SPOA for multiple sequence alignment
4. **Cluster merging**: Combines clusters with similar consensus sequences
5. **Primer trimming**: Optional primer removal using provided FASTA
6. **Stability assessment**: Evaluates consensus reliability through subsampling
7. **Output generation**: FASTA files with detailed headers including stability metrics

### External Dependencies

- **SPOA**: Required for consensus generation, must be available in PATH
- **MCL**: Optional but recommended for graph clustering
- **edlib**: Python library for edit distance calculations
- **BioPython**: Sequence handling and file I/O

### Output Structure

- Main consensus files: `{sample_name}-{cluster_num}-RiC{size}.fasta`
- Debug directory: `cluster_debug/` containing reads and untrimmed consensus
- Stability metrics included in FASTA headers (median/95th percentile distances)

## Development Notes

### No Testing Framework
This codebase does not include automated tests. Manual testing should use sample FASTQ files with known clustering expectations.

### Integration Context
Designed to replace NGSpeciesID in the ONT fungal barcoding pipeline. The tool processes demultiplexed FASTQ files and generates consensus sequences suitable for taxonomic identification.

### Configuration
All parameters are controlled via command-line arguments. No configuration files are used. Key parameters include identity thresholds, clustering algorithm choice, and sample size limits.