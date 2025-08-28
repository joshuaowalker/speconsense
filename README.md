# Speconsense

A tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads.

## Overview

Speconsense is a specialized tool for generating high-quality consensus sequences from Oxford Nanopore Technology (ONT) amplicon data. It is specifically designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

The key features of Speconsense include:
- Robust clustering of amplicon reads using either MCL graph-based or greedy algorithms
- Automatic merging of clusters with identical or similar consensus sequences
- High-quality consensus generation using SPOA
- Primer trimming for clean consensus sequences
- Stability assessment through subsampling
- Optimized for fungal amplicon datasets but suitable for any amplicon sequencing application

## Installation

### Requirements

- Python 3.8 or higher
- External dependencies:
  - [SPOA (SIMD POA)](https://github.com/rvaser/spoa) - Required (install via conda)
  - [MCL](https://micans.org/mcl/) - Optional but recommended for graph-based clustering (install via conda)

### Install from GitHub (Recommended)

The easiest way to install speconsense is directly from GitHub using pip. We recommend using a virtual environment to avoid dependency conflicts:

```bash
# Create and activate a virtual environment (recommended)
python -m venv speconsense-env
source speconsense-env/bin/activate  # On Windows: speconsense-env\Scripts\activate

# Install directly from GitHub
pip install git+https://github.com/joshuaowalker/speconsense.git

# External dependencies need to be installed separately
# SPOA: conda install bioconda::spoa
# MCL: conda install bioconda::mcl (optional but recommended)
```

After installation, the tools will be available as command-line programs:
- `speconsense` - Main clustering and consensus tool
- `speconsense-summarize` - Post-processing and summary tool

To deactivate the virtual environment when done:
```bash
deactivate
```

### Installing External Dependencies

**MCL (Markov Clustering) - Recommended:**
```bash
# Via conda (easiest)
conda install bioconda::mcl

# Or from source (more complex)
# See https://micans.org/mcl/ for source installation
```

**SPOA (SIMD POA) - Required:**
```bash
# Via conda (easiest)
conda install bioconda::spoa

# Or install from GitHub releases or build from source
# See https://github.com/rvaser/spoa for source installation instructions
```

**Note:** If MCL is not available, speconsense will automatically fall back to the greedy clustering algorithm.

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

### Running Tests

The project includes pytest-based integration tests. To run the tests:

```bash
# Set up test virtual environment
python -m venv test-venv
source test-venv/bin/activate  # On Windows: test-venv\Scripts\activate

# Install project and test dependencies
pip install -e .
pip install pytest pytest-cov

# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_augment_input.py -v

# Run with coverage
python -m pytest tests/ --cov=speconsense --cov-report=html
```

**Available Tests:**
- `tests/test_augment_input.py`: Integration tests for --augment-input functionality

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

# Control the maximum sample size for consensus generation
speconsense input.fastq --max-sample-size 500

# Filter out reads before clustering
speconsense input.fastq --presample 1000

# Merge sequence variants within 3 differences of each other
speconsense input.fastq --variant-merge-threshold 3

# Disable variant merging
speconsense input.fastq --variant-merge-threshold -1
```

### Augmenting Input with Additional Sequences

The `--augment-input` parameter allows you to include additional sequences that were not successfully matched during primary demultiplexing but were recovered through other means. This is particularly useful in the specimux workflow when using sequence recovery tools.

**Primary Use Case:**
After primary demultiplexing with specimux, some sequences may remain unmatched due to sequencing errors, primer degradation, or other issues. These sequences can be recovered using tools like `specimine` from the specimux suite and then included in speconsense clustering to maximize data utilization.

**Workflow Integration:**
```bash
# 1. Run primary demultiplexing with specimux
specimux primers.fasta specimens.txt input.fastq -O results/

# 2. Recover additional sequences from unmatched reads
specimine results/unknown/ specimen_name --output recovered.fastq

# 3. Cluster with both primary and recovered sequences
speconsense results/full/specimen_name.fastq --augment-input recovered.fastq
```

**Usage:**
```bash
# Basic usage with recovered sequences
speconsense primary_demux.fastq --augment-input recovered_sequences.fastq

# Can be combined with other options
speconsense primary_demux.fastq --augment-input recovered.fastq --primers primers.fasta
```

**Key Features:**
- Supports both FASTQ and FASTA formats (auto-detected by file extension)
- Augmented sequences participate equally in clustering with primary sequences
- During presampling, primary sequences are prioritized over augmented sequences
- All sequences (primary + augmented) appear in output files and contribute to final consensus
- Augmented sequences are fully traceable through the entire pipeline

**Important Notes:**
- Augmented sequences will only cluster with primary sequences if they meet the similarity threshold
- Recovered sequences may form separate clusters if they represent different taxa
- The final consensus sequence headers will show the total count including augmented sequences
- This approach maximizes data utilization by including sequences that would otherwise be discarded

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

#### Variant Grouping and Selection

When multiple variants exist per specimen, `speconsense-summarize` first groups similar variants together, then applies selection strategies within each group:

**Variant Grouping:**
```bash
speconsense-summarize --variant-group-identity 0.9
```
- Uses **Hierarchical Agglomerative Clustering (HAC)** to group variants with sequence identity ≥ threshold
- Default threshold is 0.9 (90% identity) 
- Variants within each group are considered similar enough to represent the same biological entity
- Each group will contribute one primary variant plus additional variants based on selection strategy

**Variant Selection (within each group):**

When multiple variants exist per specimen, `speconsense-summarize` offers two strategies for selecting which variants to output from each group:

**Size-based selection (default):**
```bash
speconsense-summarize --variant-selection size --max-variants 2
```
- Selects variants by cluster size (largest first)
- Primary variant is always the largest cluster
- Additional variants are selected in order of decreasing size
- Best for identifying the most abundant sequence variants
- Suitable when read count reflects biological abundance

**Diversity-based selection:**
```bash
speconsense-summarize --variant-selection diversity --max-variants 2
```
- Uses a **maximum distance algorithm** to select variants that are most genetically different from each other
- Primary variant is still the largest cluster
- Additional variants are selected to maximize sequence diversity in the output
- Iteratively selects the variant with the **maximum minimum distance** to all previously selected variants
- Best for capturing the full range of genetic variation in your sample
- Suitable when you want to detect distinct sequence types regardless of their abundance

**Algorithm Details for Diversity Selection (within each group):**
1. Primary variant = largest cluster within the group (by read count)
2. For each additional variant slot in the group:
   - Calculate the minimum sequence distance from each remaining candidate to all already-selected variants in this group
   - Select the candidate with the largest minimum distance (farthest from all selected in this group)
   - Repeat until max_variants reached for this group

**Overall Process:**
1. Group variants by sequence identity using HAC clustering
2. For each group independently:
   - Apply variant selection strategy (size or diversity)
   - Output up to max_variants per group
3. Final output contains representatives from all groups, ensuring both biological diversity (between groups) and appropriate sampling within each biological entity (within groups)

This two-stage process ensures that distinct biological sequences are preserved as separate groups, while providing control over variant complexity within each group.

#### Additional Summarize Options

**Quality Filtering:**
```bash
speconsense-summarize --min-ric 5
```
- Filters out consensus sequences with fewer than the specified number of Reads in Consensus (RiC)
- Default is 3 - only sequences supported by at least 3 reads are processed
- Higher values provide more stringent quality control but may exclude valid low-abundance variants

**SNP-based Variant Merging:**
```bash
speconsense-summarize --snp-merge-limit 2
```
- **Occurs before variant grouping** - merges variants within each specimen that differ by ≤ N SNP positions
- Uses a **greedy merging approach**: repeatedly finds the best pairwise merge (largest combined cluster size) until no valid merges remain
- Creates **IUPAC consensus sequences** with ambiguity codes at polymorphic positions
- Only merges variants with **identical primer sets** to maintain biological validity
- Prevents merging if the result would exceed the SNP limit or contain complex indels
- Helps consolidate very similar variants that likely represent the same biological sequence
- **Note**: This is distinct from the `--variant-merge-threshold` in the main clustering step

**Directory Control:**
```bash
speconsense-summarize --source /path/to/speconsense/output --summary-dir MyResults
```
- `--source`: Directory containing speconsense output files (default: current directory)
- `--summary-dir`: Output directory name (default: `__Summary__`)

#### Processing Workflow Summary

The complete speconsense-summarize workflow operates in this order:

1. **Load sequences** with RiC filtering (`--min-ric`)
2. **SNP-based merging** within each specimen (`--snp-merge-limit`)  
3. **HAC variant grouping** by sequence identity (`--variant-group-identity`)
4. **Variant selection** within each group (`--max-variants`, `--variant-selection`)
5. **Output generation** with full traceability and statistics

#### Enhanced Logging and Traceability

Speconsense-summarize provides comprehensive logging to help users understand processing decisions:

**Variant Analysis Logging:**
- **Complete variant summaries** for every variant in each group, including those that are skipped
- **Detailed difference categorization**: substitutions, single-nt indels, short (≤3nt) indels, and long indels
- **IUPAC-aware comparisons**: treats ambiguity codes as matches (e.g., R matches A or G)
- **Group context**: clearly shows which variants belong to each HAC clustering group
- **Selection rationale**: explains why variants were included or excluded

**Example log output:**
```
HAC clustering created 2 groups
Group 1: ['sample-c3']  
Group 2: ['sample-c1', 'sample-c2']
Processing Variants in Group 2
Primary: sample-c1 (size=403, ric=403)
Variant 1: (size=269, ric=269) - 1 short (<= 3nt) indel
Variant 2: (size=180, ric=180) - 3 substitutions, 1 single-nt indel - skipping
```

**Traceability Features:**
- **Merge history**: tracks which original clusters were combined during SNP merging
- **File lineage**: maintains connection between final outputs and original speconsense clusters  
- **Read aggregation**: `FASTQ Files/` directory contains all reads that contributed to each final consensus
- **Raw preservation**: `raw_clusters/` directory preserves original speconsense output with stability metrics

This comprehensive logging allows users to understand exactly how the pipeline processed their data and make informed decisions about parameter tuning.

### Full Command Line Options

```
usage: speconsense.py [-h] [--augment-input AUGMENT_INPUT] [--algorithm {graph,greedy}] [--min-identity MIN_IDENTITY] 
                     [--inflation INFLATION] [--min-size MIN_SIZE] [--min-cluster-ratio MIN_CLUSTER_RATIO]
                     [--max-sample-size MAX_SAMPLE_SIZE] [--presample PRESAMPLE] 
                     [--k-nearest-neighbors K_NEAREST_NEIGHBORS] [--primers PRIMERS]
                     [--stability-trials STABILITY_TRIALS] [--stability-sample STABILITY_SAMPLE]
                     [--disable-stability] [--variant-merge-threshold VARIANT_MERGE_THRESHOLD]
                     [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--version]
                     input_file

MCL-based clustering of nanopore amplicon reads

positional arguments:
  input_file            Input FASTQ file

optional arguments:
  -h, --help            show this help message and exit
  --augment-input AUGMENT_INPUT
                        Additional FASTQ/FASTA file with sequences recovered after primary demultiplexing
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
  --variant-merge-threshold VARIANT_MERGE_THRESHOLD
                        Maximum distance between consensus sequences to merge sequence variants (default: 0, -1 to disable)
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
ls *.fastq | parallel speconsense {} --primers primers.fasta
```
3. Process the output FASTA files with the speconsense-summarize tool to prepare them for downstream analysis:
```bash
speconsense-summarize
```

## Algorithm Details

### Clustering Methods

SpecConsense offers two clustering approaches with different characteristics:

#### **Graph-based clustering (MCL)** - Default and recommended
- Constructs a similarity graph between reads and applies the Markov Cluster algorithm to identify clusters
- **Tends to produce more clusters** by discriminating between sequence variants within a specimen
- Suitable when you want to identify multiple variants per specimen and are willing to interpret complex results
- Excellent for detecting subtle sequence differences and biological variation
- Relatively fast with high-quality results for most datasets

#### **Greedy clustering** - Fast and simple alternative  
- Uses greedy star clustering that iteratively selects the sequence with the most similar neighbors as cluster centers
- **Tends to produce fewer, simpler clusters** by focusing on well-separated sequences
- Excellent at discriminating between distinct targets (e.g., primary sequence vs. contaminant)
- Usually does not separate variants within the same biological sequence
- Suitable when you want a fast algorithm with minimal complexity in interpreting results
- Ideal for applications where one consensus per specimen is sufficient

#### **Choosing the Right Algorithm**

**Use `--algorithm greedy` when:**
- You want one clear consensus sequence per specimen
- Speed is prioritized over detailed variant detection  
- You prefer simpler output interpretation
- Your dataset has well-separated sequences (targets vs. contaminants)

**Use `--algorithm graph` (default) when:**
- You want to detect sequence variants within specimens
- You're willing to evaluate multiple consensus sequences per specimen
- You need fine-grained discrimination between similar sequences
- You want the most comprehensive analysis of sequence diversity
- **Note**: Use `speconsense-summarize` after clustering to manage multiple variants per specimen. Key options include `--snp-merge-limit` for merging variants differing by few SNPs, `--variant-group-identity` for grouping similar variants, and `--max-variants`/`--variant-selection` for controlling which variants to output

### Cluster Merging

After initial clustering, Speconsense can merge clusters with identical or similar consensus sequences:

1. **Identical sequence merging**: By default (`--variant-merge-threshold 0`), only clusters with identical consensus sequences are merged.

2. **Similar sequence merging**: When specifying a positive threshold (`--variant-merge-threshold N`), clusters whose consensus sequences differ by N or fewer **substitutions** will be merged. Important notes:
   - Only substitutions (SNPs) count toward the threshold
   - Homopolymer length differences are ignored (treated as sequencing artifacts)
   - Non-homopolymer indels prevent merging to ensure IUPAC consensus validity

3. **Disable merging**: Set `--variant-merge-threshold -1` to disable the merging step entirely.

This step helps eliminate redundant clusters that represent the same biological sequence but were separated during the initial clustering. The adjusted identity scoring with homopolymer normalization provides more accurate assessment of sequence similarity for merging decisions, especially for nanopore data where homopolymer length variation is common.

**Merged Cluster Output**: When clusters are merged, the consensus sequence includes IUPAC ambiguity codes at SNP positions. The FASTA header includes `snp=N` indicating the number of positions with ambiguous nucleotides. Note that due to transitive merging, the total SNP count can exceed the merge threshold.

### Consensus Generation

Consensus sequences are generated using SPOA (SIMD Partial Order Alignment) which efficiently handles the error profile of nanopore reads. For larger clusters, a random subset of reads (controlled by `--max-sample-size`) is used to generate the consensus.

### Stability Assessment

To evaluate the reliability of consensus sequences, Speconsense performs stability assessment by:
1. Generating multiple consensus sequences from random subsets of the cluster
2. Measuring the adjusted identity between each subsample consensus and the full consensus
3. Converting identity scores to distance metrics for stability reporting
4. Reporting the median and 95th percentile distances as stability metrics

Sequences with higher than expected distance should be checked for bioinformatic contamination and/or underlying biological variation, as these may indicate mixed clusters or other issues with the data.

**Adjusted Identity Scoring**: Speconsense uses the adjusted-identity algorithm with homopolymer normalization for more accurate sequence comparisons. This means that differences in homopolymer run lengths (e.g., AAA vs AAAAA) are treated as identical, which is particularly important for nanopore sequencing data where homopolymer length calling can be inconsistent.

### Primer Trimming

When a primers file is provided, Speconsense will identify and trim primer sequences from the 5' and 3' ends of consensus sequences, producing clean amplicon sequences for downstream analysis.

## Output Files

### Speconsense Core Output

For each specimen, Speconsense generates:

1. **Main consensus FASTA**: `{sample_name}-all.fasta`
   - Contains all consensus sequences for the specimen (one per cluster)
   - Uses original cluster numbering: `{sample_name}-c{cluster_num}-RiC{size}`

2. **Debug directory** (`cluster_debug/`):
   - `{sample_name}-c{cluster_num}-RiC{size}-reads.fastq`: Original reads in each cluster
   - `{sample_name}-c{cluster_num}-RiC{size}-untrimmed.fasta`: Untrimmed consensus sequences

### Speconsense-Summarize Output

When using `speconsense-summarize` for post-processing, creates `__Summary__/` directory with:

#### **Main Output Files** (summarization numbering: `-1`, `-1.v1`, `-2`):
- **Individual FASTA files**: `{sample_name}-{group}.fasta`, `{sample_name}-{group}.v{variant}.fasta`
- **Combined file**: `summary.fasta` - all final consensus sequences
- **Statistics**: `summary.txt` - sequence counts and metrics
- **Log file**: `summarize_log.txt` - complete processing log

#### **FASTQ Files/** (aggregated reads for final consensus):
- `{sample_name}-{group}.fastq` - all reads contributing to main consensus
- `{sample_name}-{group}.v{variant}.fastq` - all reads contributing to variant consensus

#### **raw_clusters/** (original speconsense numbering: `-c1`, `-c2`):
- `{sample_name}-all.fasta` - original consensus files with stability metrics
- `{sample_name}-c{cluster_num}-RiC{size}-reads.fastq` - individual cluster reads

### Naming Convention Summary

**Two distinct namespaces maintain traceability:**

| **Namespace** | **Used In** | **Format** | **Purpose** |
|---------------|-------------|------------|-------------|
| **Original** | `raw_clusters/`, `cluster_debug/` | `-c1`, `-c2`, `-c3` | Preserves speconsense clustering results |
| **Summarization** | Main output, `FASTQ Files/` | `-1`, `-1.v1`, `-2` | Post-processing groups and variants |

### Example Directory Structure
```
__Summary__/
├── sample-1.fasta                           # Main consensus (group 1)
├── sample-1.v1.fasta                        # Additional variant
├── sample-2.fasta                           # Second organism group
├── summary.fasta                            # All sequences combined
├── summary.txt                              # Statistics
├── summarize_log.txt                        # Processing log
├── FASTQ Files/                             # Reads for final consensus
│   ├── sample-1.fastq                       # All reads for main consensus
│   ├── sample-1.v1.fastq                    # All reads for variant
│   └── sample-2.fastq                       # All reads for second group
└── raw_clusters/                            # Original clustering results
    ├── sample-all.fasta                     # Original consensus with stability metrics
    ├── sample-c1-RiC500-reads.fastq         # Original cluster 1 reads
    └── sample-c2-RiC300-reads.fastq         # Original cluster 2 reads
```

### FASTA Header Metadata

Consensus sequence headers contain metadata fields separated by spaces:

**Core Fields (Always Present):**
- `size=N` - Total number of raw sequence reads in the cluster
- `ric=N` - **Reads in Consensus** - Number of reads actually used for consensus generation (may be less than size due to sampling limits)

**Optional Fields:**
- `snp=N` - Number of SNP positions with IUPAC ambiguity codes (only present in merged variants from speconsense-summarize)
- `primers=list` - Comma-separated list of detected primer names (e.g., `primers=ITS1F,ITS4`)
- `median_diff=X.X` - Median edit distance from stability assessment (debug files only)  
- `p95_diff=X.X` - 95th percentile edit distance from stability assessment (debug files only)

**Example Headers:**
```
# Simple consensus from speconsense:
>sample-c1-RiC45;size=50 ric=45 median_diff=2.1 p95_diff=4.8 primers=ITS1F,ITS4

# Merged variant from speconsense-summarize (stability metrics removed):
>sample-1;size=95 ric=78 snp=2 primers=ITS1F,ITS4
```

**Notes:**
- Length is not included (easily calculated from sequence)
- Stability metrics (`median_diff`, `p95_diff`) are preserved in debug files but removed from final speconsense-summarize output
- Variant merging only occurs between sequences with identical primer sets
- SNP counts reflect IUPAC ambiguity positions in consensus sequences

## Changelog

### Version 0.3.0 (2025-08-28)

**Major Performance Optimizations:**
- **Dramatically improved file I/O performance** - Replaced BioPython FASTQ parsing with direct file concatenation, achieving orders of magnitude speedup
- **Eliminated directory scanning bottleneck** - Single lookup table build replaces hundreds of glob operations
- **Optimized raw file copying** - Pre-built file lookup system scales efficiently with large datasets

**Enhanced Variant Processing:**
- **Complete variant selection framework** - Added size-based and diversity-based selection strategies with configurable limits
- **Hierarchical variant grouping** - Implemented HAC clustering to separate distinct biological sequences
- **SNP-based variant merging** - Greedy merging approach with IUPAC consensus generation
- **Comprehensive variant analysis** - Detailed logging with difference categorization (substitutions, indels by length)

**Algorithm Selection Guide:**
- **Expanded clustering documentation** - Clear guidance on when to use greedy vs. graph clustering
- **Technical algorithm descriptions** - Proper terminology and algorithmic details
- **User decision framework** - Specific use cases and parameter recommendations

**Documentation Improvements:**
- **Complete workflow documentation** - Step-by-step processing pipeline with all options explained
- **Enhanced logging features** - IUPAC-aware comparisons and full traceability through processing steps
- **Parameter reference** - Comprehensive guide to all summarize options with examples
- **Example outputs** - Real log examples and file structure illustrations

**Code Quality:**
- **Removed legacy code** - Eliminated all non-optimized function versions
- **Simplified interfaces** - Clean function signatures and reduced conditional logic
- **Maintained backward compatibility** - All existing functionality preserved

### Version 0.2.1 (2025-08-27)

**Core Functionality:**
- Stable clustering and consensus generation with MCL and greedy algorithms
- Comprehensive variant merging with adjusted identity scoring
- SPOA-based consensus generation with stability assessment
- Complete output structure with traceability features

## Citations

This project uses and builds upon:

- **MCL clustering algorithm**: van Dongen, Stijn, *Graph clustering via a discrete uncoupling process*, Siam Journal on Matrix Analysis and Applications 30-1, p121-141, 2008. (https://doi.org/10.1137/040608635)

- **ONT fungal barcoding protocol**: Russell, S.D., Geurin, Z., Walker, J. (2024). *Primary Data Analysis - Basecalling, Demultiplexing, and Consensus Building for ONT Fungal Barcodes*. protocols.io. https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v4

## License

[BSD 3-Clause License](LICENSE)

## Name

Speconsense is a playful portmanteau of "Specimen" and "Consensus", reflecting the tool's focus on generating high-quality consensus sequences from specimen amplicon reads.