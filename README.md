# Speconsense

A tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads.

## Overview

Speconsense is a specialized tool for generating high-quality consensus sequences from Oxford Nanopore Technology (ONT) amplicon data. It is specifically designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

The key features of Speconsense include:
- Robust clustering of amplicon reads using either MCL graph-based or greedy algorithms
- Automatic merging of clusters with identical or similar consensus sequences
- High-quality consensus generation using SPOA
- Primer trimming for clean consensus sequences
- Optional polishing with Medaka
- Stability assessment through subsampling
- Optimized for fungal amplicon datasets but suitable for any amplicon sequencing application

## Installation

### Requirements

- Python 3.8 or higher
- External dependencies:
  - [SPOA (SIMD POA)](https://github.com/rvaser/spoa) - Must be in PATH
  - [MCL](https://micans.org/mcl/) - Optional but recommended for graph-based clustering
  - [Medaka](https://github.com/nanoporetech/medaka) - Optional for consensus polishing

### Install from GitHub (Recommended)

The easiest way to install speconsense is directly from GitHub using pip. We recommend using a virtual environment to avoid dependency conflicts:

```bash
# Create and activate a virtual environment (recommended)
python -m venv speconsense-env
source speconsense-env/bin/activate  # On Windows: speconsense-env\Scripts\activate

# Install directly from GitHub
pip install git+https://github.com/joshuaowalker/speconsense.git

# External dependencies need to be installed separately
# SPOA can be installed from https://github.com/rvaser/spoa
# MCL can be installed from https://micans.org/mcl/
# Medaka can be installed with: pip install medaka
```

After installation, the tools will be available as command-line programs:
- `speconsense` - Main clustering and consensus tool
- `speconsense-summarize` - Post-processing and summary tool

To deactivate the virtual environment when done:
```bash
deactivate
```

### Development Installation

For development or if you want to modify the code:

```bash
# Create and activate a virtual environment (recommended)
python -m venv speconsense-dev
source speconsense-dev/bin/activate  # On Windows: speconsense-dev\Scripts\activate

# Clone the repository
git clone https://github.com/joshuaowalker/speconsense.git
cd speconsense

# Install in editable mode
pip install -e .
```

## Usage

### Basic Usage

```bash
speconsense input.fastq
```

By default, this will:
1. Cluster reads using graph-based MCL algorithm
2. Merge clusters with identical consensus sequences
3. Generate a consensus sequence for each cluster
4. Output FASTA files containing consensus sequences

### Common Options

```bash
# Use greedy clustering algorithm instead of MCL
speconsense input.fastq --algorithm greedy

# Set minimum cluster size
speconsense input.fastq --min-size 10

# Set minimum identity threshold for clustering
speconsense input.fastq --min-identity 0.85

# Enable primer trimming
speconsense input.fastq --primers primers.fasta

# Enable Medaka polishing
speconsense input.fastq --medaka

# Control the maximum sample size for consensus generation
speconsense input.fastq --max-sample-size 500

# Filter out reads before clustering
speconsense input.fastq --presample 1000

# Merge clusters with consensus sequences within 3 edits of each other
speconsense input.fastq --max-consensus-distance 3

# Disable cluster merging
speconsense input.fastq --max-consensus-distance -1
```

### Post-processing and Summary

After running speconsense, use the summarize tool to process outputs:

```bash
# Generate summary with default settings
speconsense-summarize

# Custom minimum RiC threshold and output directory
speconsense-summarize --min-ric 5 --summary-dir MyResults

# Process specific source directory
speconsense-summarize --source /path/to/speconsense/output
```

### Full Command Line Options

```
usage: speconsense.py [-h] [--augment-input AUGMENT_INPUT] [--algorithm {graph,greedy}] [--min-identity MIN_IDENTITY] 
                     [--inflation INFLATION] [--min-size MIN_SIZE] [--min-cluster-ratio MIN_CLUSTER_RATIO]
                     [--max-sample-size MAX_SAMPLE_SIZE] [--presample PRESAMPLE] 
                     [--k-nearest-neighbors K_NEAREST_NEIGHBORS] [--primers PRIMERS]
                     [--stability-trials STABILITY_TRIALS] [--stability-sample STABILITY_SAMPLE]
                     [--disable-stability] [--medaka] [--max-consensus-distance MAX_CONSENSUS_DISTANCE]
                     [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--version]
                     input_file

MCL-based clustering of nanopore amplicon reads

positional arguments:
  input_file            Input FASTQ file

optional arguments:
  -h, --help            show this help message and exit
  --augment-input AUGMENT_INPUT
                        Additional input FASTQ file with mined sequences
  --algorithm {graph,greedy}
                        Clustering algorithm to use (default: graph)
  --min-identity MIN_IDENTITY
                        Minimum sequence identity threshold (default: 0.85)
  --inflation INFLATION
                        MCL inflation parameter (default: 4.0)
  --min-size MIN_SIZE   Minimum cluster size (default: 5)
  --min-cluster-ratio MIN_CLUSTER_RATIO
                        Minimum size ratio between a cluster and the largest cluster (default: 0.2)
  --max-sample-size MAX_SAMPLE_SIZE
                        Maximum cluster size for consensus (default: 500)
  --presample PRESAMPLE
                        Presample size for initial reads (default: 1000, 0 to disable)
  --k-nearest-neighbors K_NEAREST_NEIGHBORS
                        Number of nearest neighbors for graph construction (default: 5)
  --primers PRIMERS     FASTA file containing primer sequences
  --stability-trials STABILITY_TRIALS
                        Number of sampling trials to assess stability (default: 100)
  --stability-sample STABILITY_SAMPLE
                        Size of stability samples (default: 20)
  --disable-stability   Disable stability assessment
  --medaka              Enable consensus polishing with medaka
  --max-consensus-distance MAX_CONSENSUS_DISTANCE
                        Maximum edit distance between consensus sequences to merge clusters (default: 0, -1 to disable)
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set logging level (default: INFO)
  --version             Show program's version number and exit
```

## Integration with ONT Fungal Barcoding Pipeline

SpecConsense is designed to replace the NGSpeciesID step in the [ONT DNA Barcoding Fungal Amplicons protocol](https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v4). To use SpecConsense in this pipeline:

1. Follow the protocol through the "Demultiplex the Reads with Specimux" step
2. Instead of running NGSpeciesID, run SpecConsense on each demultiplexed FASTQ file:

```bash
# Replace this command from the protocol:
ls *.fastq | parallel NGSpeciesID --ont --consensus --t 1 --abundance_ratio 0.2 --top_reads --sample_size 500 --symmetric_map_align_thresholds --aligned_threshold 0.75 --mapped_threshold 1.0 --medaka --fastq {} --outfolder {.}

# With this command using Speconsense:
ls *.fastq | parallel speconsense {} --medaka --primers primers.fasta
```
3. Process the output FASTA files with the speconsense-summarize tool to prepare them for downstream analysis:
```bash
speconsense-summarize
```

## Algorithm Details

### Clustering Methods

SpecConsense offers two clustering approaches:

1. **Graph-based clustering (MCL)**: Constructs a similarity graph between reads and applies the Markov Cluster algorithm to identify clusters. This is the default and recommended method for most datasets.

2. **Greedy clustering**: A simpler approach that iteratively finds seeds with the most connections above the similarity threshold to form clusters. This is provided as a fallback when MCL is not available.

### Cluster Merging

After initial clustering, Speconsense can merge clusters with identical or similar consensus sequences:

1. **Identical consensus merging**: By default (`--max-consensus-distance 0`), only clusters with identical consensus sequences are merged.

2. **Similar consensus merging**: When specifying a positive distance threshold (`--max-consensus-distance N`), clusters whose consensus sequences differ by N or fewer edits will be merged.

3. **Disable merging**: Set `--max-consensus-distance -1` to disable the merging step entirely.

This step helps eliminate redundant clusters that represent the same biological sequence but were separated during the initial clustering.

### Consensus Generation

Consensus sequences are generated using SPOA (SIMD Partial Order Alignment) which efficiently handles the error profile of nanopore reads. For larger clusters, a random subset of reads (controlled by `--max-sample-size`) is used to generate the consensus.

### Stability Assessment

To evaluate the reliability of consensus sequences, Speconsense performs stability assessment by:
1. Generating multiple consensus sequences from random subsets of the cluster
2. Measuring the edit distance between each subsample consensus and the full consensus
3. Reporting the median and 95th percentile distances as stability metrics

Sequences with higher than expected distance should be checked for bioinformatic contamination and/or underlying biological variation, as these may indicate mixed clusters or other issues with the data.

### Primer Trimming

When a primers file is provided, Speconsense will identify and trim primer sequences from the 5' and 3' ends of consensus sequences, producing clean amplicon sequences for downstream analysis.

## Output Files

For each cluster, Speconsense generates:

1. **Main consensus FASTA**: `{sample_name}-{cluster_num}-RiC{size}.fasta`
   - Contains the final (trimmed, if primers provided) consensus sequence
   - Header includes cluster size, stability metrics, and found primers

2. **Debug directory files**:
   - `{sample_name}-{cluster_num}-RiC{size}-reads.fastq`: Original reads in the cluster
   - `{sample_name}-{cluster_num}-RiC{size}-untrimmed.fasta`: Untrimmed consensus sequence

## Citations

This project uses and builds upon:

- **MCL clustering algorithm**: van Dongen, Stijn, *Graph clustering via a discrete uncoupling process*, Siam Journal on Matrix Analysis and Applications 30-1, p121-141, 2008. (https://doi.org/10.1137/040608635)

- **ONT fungal barcoding protocol**: Russell, S.D., Geurin, Z., Walker, J. (2024). *Primary Data Analysis - Basecalling, Demultiplexing, and Consensus Building for ONT Fungal Barcodes*. protocols.io. https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v4

## License

[BSD 3-Clause License](LICENSE)

## Name

Speconsense is a playful portmanteau of "Specimen" and "Consensus", reflecting the tool's focus on generating high-quality consensus sequences from specimen amplicon reads.