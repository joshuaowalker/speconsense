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

**speconsense/core/** - Core clustering and consensus generation subpackage:
- `clusterer.py`: `SpecimenClusterer` class - main orchestrator for the clustering pipeline
- `workers.py`: Worker functions for parallel processing (SPOA, cluster processing)
- `cli.py`: Command-line interface and argument parsing
- Two clustering algorithms: graph-based MCL (default) and greedy clustering
- Consensus generation via external SPOA tool
- Cluster merging based on consensus sequence similarity
- Stability assessment through subsampling

**speconsense/summarize/** - Post-processing utility subpackage:
- `cli.py`: Command-line interface and main entry point
- `iupac.py`: IUPAC ambiguity code utilities and distance calculations
- `fields.py`: FASTA header field classes and formatting
- `analysis.py`: MSA analysis and quality assessment
- `merging.py`: MSA-based variant merging with IUPAC consensus
- `clustering.py`: HAC clustering and variant selection
- `io.py`: File I/O operations (loading sequences, writing outputs)

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
2. HAC variant grouping by sequence identity (separates dissimilar sequences)
3. Homopolymer-aware MSA-based merging within each group (creates IUPAC consensus sequences)
4. Selection size ratio filtering (`--select-min-size-ratio`) — removes tiny post-merge variants
5. Variant selection within each group (size-based or diversity-based)
6. Optional full consensus generation per group (`--enable-full-consensus`) — IUPAC consensus from all pre-merge variants, gaps never win
7. Output generation with full traceability

Note: HAC grouping occurs BEFORE merging (since 0.4.0) to prevent inappropriate merging of dissimilar sequences (e.g., contaminants with primary targets). Homopolymer-aware merging (since 0.5.0) distinguishes structural indels from homopolymer length differences.

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
  - Summarization: `-1.v1`, `-1.v2`, `-2.v1` format for groups and variants
  - Full consensus: `-1.full` format (one per group, when `--enable-full-consensus` is used)

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

### MCL Graph Construction (Design Decision)

The K-NN similarity graph for MCL uses an **asymmetric edge storage pattern** where `similarities[id1]` only contains entries for neighbors `id2 > id1` (lexicographically). This is a weakly-held design decision that produces good clustering results despite the MCL documentation recommending symmetric graphs.

**Why asymmetric?** The pattern emerged from the original implementation and affects tie-breaking when multiple neighbors have identical similarity scores. Changing to symmetric storage would alter which neighbors are selected during K-NN construction, potentially changing clustering results.

**Key implementation details** (in `scalability/base.py`):
```python
# Asymmetric pattern - matches main branch behavior
similarities[id1][id2] = score
similarities.setdefault(id2, {})[id1] = score  # Gets overwritten when id2 is processed
```

**If considering symmetric graphs in the future:**
- Would require careful validation against existing test cases
- May improve MCL convergence (per MCL documentation)
- Would change clustering results - requires re-tuning or acceptance of different outputs

### Integration Context

Designed to replace NGSpeciesID in the ONT fungal barcoding pipeline from protocols.io. Processes demultiplexed FASTQ files and generates consensus sequences suitable for taxonomic identification. Typically used downstream of specimux for demultiplexing.

### Configuration

All parameters controlled via command-line arguments. No configuration files. Key parameters:
- Identity thresholds (`--min-identity`)
- Clustering algorithm choice (`--algorithm`)
- Sample size limits (`--max-sample-size`, `--presample`)
- Cluster size filtering (`--min-size`, `--min-cluster-ratio`)
- Primer handling (`--primers`, `--orient-mode`)

### Local Development Scripts

The `scripts/` directory (git-ignored) contains local development utilities:
- `verify_refactor.py`: AST-based verification of refactoring equivalence against current code
- `verify_refactor_at_commit.py`: Same verification against specific git commits

These scripts compare original monolithic files against refactored subpackages to detect transcription errors. Useful for future refactoring work.
