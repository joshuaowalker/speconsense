# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Speconsense is a Python tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads. It's designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

## Key Commands

### Running the main tool
```bash
# Basic usage
python speconsense.py input.fastq

# With common options
python speconsense.py input.fastq --algorithm greedy --min-size 10 --primers primers.fasta --medaka

# Test with sample data
python speconsense.py --min-identity 0.85 --max-sample-size 500 --presample 1000
```

### Post-processing
```bash
# Generate summary from output files
python summarize_speconsense.py

# With custom parameters
python summarize_speconsense.py --min-ric 5 --source /path/to/output --summary-dir MyResults
```

### Dependencies installation
```bash
# Install Python dependencies
pip install -r requirements.txt

# External tools (must be in PATH):
# - SPOA (SIMD POA): https://github.com/rvaser/spoa
# - MCL: https://micans.org/mcl/ (optional but recommended)
# - Medaka: pip install medaka (optional)
```

## Code Architecture

### Main Components

**speconsense.py** - Core clustering and consensus tool:
- `SpecimenClusterer` class: Main orchestrator for the clustering pipeline
- `AmbiguityMode` class: Controls handling of ambiguous nucleotides using IUPAC codes
- Two clustering algorithms: graph-based MCL (default) and greedy clustering
- Consensus generation using external SPOA tool
- Optional Medaka polishing integration
- Cluster merging based on consensus sequence similarity

**summarize_speconsense.py** - Post-processing utility:
- Parses Speconsense FASTA output files
- Generates summary reports and filtered datasets
- Handles RiC (Reads in Consensus) filtering

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
- **Medaka**: Optional for consensus polishing
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