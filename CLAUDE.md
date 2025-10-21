# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Speconsense is a Python tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads. It's designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

## Key Commands

### Running the tools
```bash
# Main clustering tool (basic usage)
speconsense input.fastq

# With common options
speconsense input.fastq --algorithm greedy --min-size 10 --primers primers.fasta

# Post-processing tool
speconsense-summarize --min-ric 5 --source /path/to/output --summary-dir MyResults

# Synthetic data generator for testing
speconsense-synth reference.fasta --num-reads 1000 --error-rate 0.05 --output synthetic.fastq
```

### Development setup
```bash
# Install in development mode
pip install -e .

# External dependencies (must be in PATH)
conda install bioconda::spoa  # Required
conda install bioconda::mcl   # Optional but recommended
```

### Testing
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_orientation.py

# Run tests with specific markers
pytest -m "not slow"
```

## Code Architecture

### Main Components

**speconsense/core.py** - Core clustering and consensus generation:
- `SpecimenClusterer` class: Main orchestrator for the clustering pipeline
- Two clustering algorithms: graph-based MCL (default) and greedy clustering
- Consensus generation via external SPOA tool
- Cluster merging based on consensus sequence similarity
- Stability assessment through subsampling

**speconsense/summarize.py** - Post-processing utility:
- Parses Speconsense FASTA output files
- SNP-based variant merging using IUPAC ambiguity codes
- Hierarchical Agglomerative Clustering (HAC) for variant grouping
- Variant selection strategies (size-based or diversity-based)
- Generates summary reports and filtered datasets

**speconsense/cli.py** - Command-line interface entry point for main tool

**speconsense/synth.py** - Synthetic read generator for testing consensus algorithms

### Key Processing Pipeline

1. **Input processing**: Read FASTQ files, optional presampling and augmentation
2. **Orientation (optional)**: Detect and correct sequence orientation based on primers
3. **Clustering**: Either MCL graph-based or greedy algorithm
4. **Consensus generation**: Uses SPOA for multiple sequence alignment
5. **Cluster merging**: Combines clusters with identical/homopolymer-equivalent consensus sequences
6. **Primer trimming**: Optional primer removal using provided FASTA
7. **Stability assessment**: Evaluates consensus reliability through subsampling
8. **Output generation**: FASTA files with detailed headers including stability metrics

### Post-Processing Pipeline (speconsense-summarize)

1. Load sequences with RiC filtering
2. SNP-based merging within each specimen (creates IUPAC consensus sequences)
3. HAC variant grouping by sequence identity
4. Variant selection within each group (size-based or diversity-based)
5. Output generation with full traceability

### External Dependencies

- **SPOA**: Required for consensus generation, must be in PATH
- **MCL**: Optional but recommended for graph clustering
- **edlib**: Python library for edit distance calculations
- **adjusted-identity**: Custom library for IUPAC-aware sequence alignment (from GitHub)
- **BioPython**: Sequence handling and file I/O
- **NumPy**: Numerical operations
- **tqdm**: Progress bars

### Output Structure

**Main consensus files** (from speconsense):
- `{sample_name}-all.fasta` - All consensus sequences for a specimen
- `cluster_debug/` directory:
  - `{sample_name}-c{num}-RiC{size}-reads.fastq` - Original reads
  - `{sample_name}-c{num}-RiC{size}-sampled.fastq` - Sampled reads for consensus
  - `{sample_name}-c{num}-RiC{size}-untrimmed.fasta` - Untrimmed consensus

**Summary files** (from speconsense-summarize):
- `__Summary__/` directory with processed outputs
- Dual namespace system:
  - Original clustering: `-c1`, `-c2`, `-c3` format
  - Summarization: `-1`, `-1.v1`, `-2` format for groups and variants

### IUPAC Ambiguity Code Handling

The codebase uses IUPAC nucleotide ambiguity codes throughout:
- `IUPAC_CODES` dict maps nucleotide sets to codes (R=A/G, Y=C/T, etc.)
- `IUPAC_EQUIV` list enables edlib alignment to treat ambiguity codes as matching their constituent bases
- `STANDARD_ADJUSTMENT_PARAMS` defines consistent sequence comparison parameters:
  - Homopolymer normalization enabled (treats "AAA" = "AAAAA")
  - IUPAC overlap disabled (uses standard IUPAC semantics: Y≠M)
  - No end trimming (`end_skip_distance=0`)
  - Single-base repeats for homopolymer normalization

### Algorithm Selection

**Graph-based MCL (default)**:
- Better at detecting sequence variants within specimens
- More clusters, higher granularity
- Requires MCL tool in PATH (falls back to greedy if missing)

**Greedy clustering (`--algorithm greedy`)**:
- Faster, simpler output
- Fewer clusters, focuses on well-separated sequences
- Good for detecting distinct targets vs. contaminants

### Integration Context

Designed to replace NGSpeciesID in the ONT fungal barcoding pipeline from protocols.io. Processes demultiplexed FASTQ files and generates consensus sequences suitable for taxonomic identification. Typically used downstream of specimux for demultiplexing.

### Configuration

All parameters controlled via command-line arguments. No configuration files. Key parameters:
- Identity thresholds (`--min-identity`)
- Clustering algorithm choice (`--algorithm`)
- Sample size limits (`--max-sample-size`, `--presample`)
- Cluster size filtering (`--min-size`, `--min-cluster-ratio`)
- Primer handling (`--primers`, `--orient-mode`)
