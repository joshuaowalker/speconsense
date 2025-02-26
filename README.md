# Speconsense

A tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads.

## Overview

Speconsense is a specialized tool for generating high-quality consensus sequences from Oxford Nanopore Technology (ONT) amplicon data. It is specifically designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

The key features of Speconsense include:
- Robust clustering of amplicon reads using either MCL graph-based or greedy algorithms
- High-quality consensus generation using SPOA
- Primer trimming for clean consensus sequences
- Optional polishing with Medaka
- Stability assessment through subsampling
- Optimized for fungal amplicon datasets but suitable for any amplicon sequencing application

## Installation

### Requirements

- Python 3.8 or higher
- Required Python packages:
  - edlib
  - numpy
  - biopython
  - tqdm
- External dependencies:
  - [SPOA (SIMD POA)](https://github.com/rvaser/spoa) - Must be in PATH
  - [MCL](https://micans.org/mcl/) - Optional but recommended for graph-based clustering
  - [Medaka](https://github.com/nanoporetech/medaka) - Optional for consensus polishing

```bash
# Clone the repository
git clone https://github.com/joshuaowalker/speconsense.git
cd speconsense

# Install required Python packages
pip install -r requirements.txt

# External dependencies need to be installed separately
# SPOA can be installed from https://github.com/rvaser/spoa
# MCL can be installed from https://micans.org/mcl/
# Medaka can be installed with: pip install medaka
```

## Usage

### Basic Usage

```bash
python speconsense.py input.fastq
```

By default, this will:
1. Cluster reads using graph-based MCL algorithm
2. Generate a consensus sequence for each cluster
3. Output FASTA files containing consensus sequences

### Common Options

```bash
# Use greedy clustering algorithm instead of MCL
python speconsense.py input.fastq --algorithm greedy

# Set minimum cluster size
python speconsense.py input.fastq --min-size 10

# Set minimum identity threshold for clustering
python speconsense.py input.fastq --min-identity 0.85

# Enable primer trimming
python speconsense.py input.fastq --primers primers.fasta

# Enable Medaka polishing
python speconsense.py input.fastq --medaka

# Control the maximum sample size for consensus generation
python speconsense.py input.fastq --max-sample-size 500

# Filter out reads before clustering
python speconsense.py input.fastq --presample 1000
```

### Full Command Line Options

```
usage: speconsense.py [-h] [--algorithm {graph,greedy}] [--min-identity MIN_IDENTITY] [--inflation INFLATION]
                     [--min-size MIN_SIZE] [--min-yield MIN_YIELD] [--max-sample-size MAX_SAMPLE_SIZE]
                     [--presample PRESAMPLE] [--k-nearest-neighbors K_NEAREST_NEIGHBORS] [--primers PRIMERS]
                     [--stability-trials STABILITY_TRIALS] [--stability-sample STABILITY_SAMPLE]
                     [--disable-stability] [--medaka] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                     input_file

MCL-based clustering of nanopore amplicon reads

positional arguments:
  input_file            Input FASTQ file

optional arguments:
  -h, --help            show this help message and exit
  --algorithm {graph,greedy}
                        Clustering algorithm to use (default: graph)
  --min-identity MIN_IDENTITY
                        Minimum sequence identity threshold (default: 0.85)
  --inflation INFLATION
                        MCL inflation parameter (default: 4.0)
  --min-size MIN_SIZE   Minimum cluster size (default: 5)
  --min-yield MIN_YIELD
                        Minimum fraction of sequences to be represented in clusters (default: 0.8)
  --max-sample-size MAX_SAMPLE_SIZE
                        Maximum cluster size for consensus (default: 500)
  --presample PRESAMPLE
                        Presample size for initial reads (default: 1000, 0 to disable)
  --k-nearest-neighbors K_NEAREST_NEIGHBORS
                        Number of nearest neighbors for graph construction (default: 20)
  --primers PRIMERS     FASTA file containing primer sequences
  --stability-trials STABILITY_TRIALS
                        Number of sampling trials to assess stability (default: 100)
  --stability-sample STABILITY_SAMPLE
                        Size of stability samples (default: 20)
  --disable-stability   Disable stability assessment
  --medaka              Enable consensus polishing with medaka
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set logging level (default: INFO)
```

## Integration with ONT Fungal Barcoding Pipeline

SpecConsense is designed to replace the NGSpeciesID step in the [ONT DNA Barcoding Fungal Amplicons protocol](https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v4). To use SpecConsense in this pipeline:

1. Follow the protocol through the "Demultiplex the Reads with Specimux" step
2. Instead of running NGSpeciesID, run SpecConsense on each demultiplexed FASTQ file:

```bash
# Replace this command from the protocol:
ls *.fastq | parallel NGSpeciesID --ont --consensus --t 1 --abundance_ratio 0.2 --top_reads --sample_size 500 --symmetric_map_align_thresholds --aligned_threshold 0.75 --mapped_threshold 1.0 --medaka --fastq {} --outfolder {.}

# With this command using Speconsense:
ls *.fastq | parallel python speconsense.py {} --medaka --primers primers.fasta
```

3. Process the output files as needed for your downstream analysis

## Algorithm Details

### Clustering Methods

SpecConsense offers two clustering approaches:

1. **Graph-based clustering (MCL)**: Constructs a similarity graph between reads and applies the Markov Cluster algorithm to identify clusters. This is the default and recommended method for most datasets.

2. **Greedy clustering**: A simpler approach that iteratively finds seeds with the most connections above the similarity threshold to form clusters. This is provided as a fallback when MCL is not available.

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