# Speconsense

A tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads.

## Overview

Speconsense is a specialized tool for generating high-quality consensus sequences from Oxford Nanopore Technology (ONT) amplicon data. It is specifically designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

The key features of Speconsense include:
- Robust clustering of amplicon reads using either Markov Clustering (MCL) graph-based or greedy algorithms
- **Automatic variant phasing**: Detects variant positions within clusters and splits reads into separate haplotypes
- Automatic merging of clusters with identical or similar consensus sequences
- High-quality consensus generation using SPOA
- Primer trimming for clean consensus sequences
- Read identity metrics for quality assessment
- IUPAC ambiguity codes for unphased heterozygous positions
- Optimized for fungal amplicon datasets but suitable for any amplicon sequencing application

## Installation

### Requirements

- Python 3.8 or higher
- External dependencies:
  - [SPOA (SIMD POA)](https://github.com/rvaser/spoa) - Required (install via conda)
  - [MCL](https://micans.org/mcl/) - Optional but recommended for graph-based clustering (install via conda)
  - [vsearch](https://github.com/torognes/vsearch) - Optional, required for scalability mode with large datasets (install via conda)

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
- `speconsense-synth` - Synthetic read generator for testing

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

**vsearch - Optional (for scalability mode):**
```bash
# Via conda (easiest)
conda install bioconda::vsearch

# Or download from https://github.com/torognes/vsearch/releases
```

**Note:** If the mcl tool is not available, speconsense will automatically fall back to the greedy clustering algorithm.

## Usage

### Basic Usage

```bash
speconsense input.fastq
```

By default, this will:
1. Cluster reads using graph-based Markov Clustering (MCL) algorithm
2. Merge clusters with identical consensus sequences
3. Generate a consensus sequence for each cluster
4. Output FASTA files containing consensus sequences

### Common Options

```bash
# Use greedy clustering algorithm instead of Markov Clustering
speconsense input.fastq --algorithm greedy

# Set minimum cluster size
speconsense input.fastq --min-size 10

# Set minimum identity threshold for clustering (default: 0.9)
speconsense input.fastq --min-identity 0.85

# Control the maximum sample size for consensus generation (default: 100)
speconsense input.fastq --max-sample-size 200

# Specify output directory (default: clusters)
speconsense input.fastq --output-dir results/

# Disable automatic variant phasing
speconsense input.fastq --disable-position-phasing

# Using short form for output directory
speconsense input.fastq -O my_results/
```

### Using Profiles

Profiles save parameter configurations for different workflows:

```bash
# List available profiles
speconsense --list-profiles

# Use a profile (CLI arguments override profile values)
speconsense input.fastq -p herbarium
speconsense input.fastq -p herbarium --min-size 10
```

**Bundled profiles:**
- `herbarium` — High-recall for degraded DNA/type specimens
- `largedata` — Experimental settings for large input files
- `nostalgia` — Simulate older bioinformatics pipelines
- `strict` — High-precision for confident results

On first use, an `example.yaml` template is created in `~/.config/speconsense/profiles/` — copy and edit it to create custom profiles.

### Two-Phase Processing Philosophy

**speconsense** extracts distinct biological sequences from raw reads within a single specimen. It separates true biological variation (alleles, loci, contaminants) from sequencing noise by clustering reads and generating consensus sequences. Philosophy: *split rather than conflate* — better to produce separate clusters that can be merged downstream than to lose real biological distinctions.

**speconsense-summarize** curates and consolidates consensus sequences for downstream analysis. It handles residual ambiguity that consensus generation cannot resolve — primarily minor SNP differences that may represent the same biological entity — groups related variants, selects representatives, and can analyze patterns across a run. Philosophy: *conservative simplification* — merge only what is clearly equivalent, preserve distinctions when uncertain.

| Aspect | speconsense | speconsense-summarize |
|--------|-------------|----------------------|
| Input | Raw reads | Consensus sequences |
| Scope | Single specimen | Single or multiple specimens |
| Goal | Discovery/enumeration | Curation/simplification |
| Error model | Read-level noise | Residual ambiguity |
| Bias | Prefer over-splitting | Conservative merging |

### Post-processing with Speconsense-Summarize

After running speconsense, use the summarize tool to process and refine outputs:

```bash
# Generate summary with default settings
speconsense-summarize

# Custom minimum RiC (Reads in Consensus) threshold
speconsense-summarize --min-ric 5

# Process specific source directory and custom output
speconsense-summarize --source /path/to/speconsense/output --summary-dir MyResults
```

**Variant Detection and Haplotype Phasing:**

Speconsense can detect and isolate sequence variants within specimens (aka "haplotype phasing"). The graph-based clustering algorithm excels at discriminating between variants, and `speconsense-summarize` provides sophisticated tools for managing multiple variants per specimen, including:
- MSA-based variant merging with IUPAC ambiguity codes and size-weighted consensus
- Hierarchical variant grouping to separate contaminants from primary targets
- Size-based or diversity-based variant selection strategies

For detailed information on variant handling options, see the [Advanced Post-Processing](#advanced-post-processing) section below.

## Quick Start: Complete Workflow

Speconsense is designed to replace the NGSpeciesID step in the [ONT DNA Barcoding Fungal Amplicons protocol](https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v4). Here's a complete workflow:

**1. Demultiplex with Specimux**

Follow the protocol through the "Demultiplex the Reads with Specimux" step.

**2. Run Speconsense**

Replace the NGSpeciesID step with speconsense:

```bash
# Instead of NGSpeciesID:
# ls *.fastq | parallel NGSpeciesID --ont --consensus --t 1 --abundance_ratio 0.2 ...

# Use Speconsense:
ls *.fastq | parallel speconsense {}
```

**3. Post-process with Speconsense-Summarize**

Process the outputs to prepare for downstream analysis:

```bash
speconsense-summarize
```

This will create a `__Summary__/` directory with:
- `summary.fasta` - All final consensus sequences (merged variants only)
- Individual FASTA files per specimen
- `summary.txt` - Statistics and metrics
- `quality_report.txt` - Prioritized list of sequences with potential quality concerns
- `FASTQ Files/` - Reads contributing to each final consensus
- `variants/` - Pre-merge variant FASTA and FASTQ files (for merged sequences)

## Output Files

### Speconsense Core Output

For each specimen, Speconsense generates:

1. **Main consensus FASTA**: `{sample_name}-all.fasta`
   - Contains all consensus sequences for the specimen (one per cluster)
   - Uses original cluster numbering: `{sample_name}-c{cluster_num}-RiC{size}`

2. **Debug directory** (`cluster_debug/`):
   - `{sample_name}-c{cluster_num}-RiC{size}-reads.fastq`: Original reads in each cluster
   - `{sample_name}-c{cluster_num}-RiC{size}-sampled.fastq`: Sampled reads used for consensus generation
   - `{sample_name}-c{cluster_num}-RiC{size}-untrimmed.fasta`: Untrimmed consensus sequences

### Speconsense-Summarize Output

When using `speconsense-summarize` for post-processing, creates `__Summary__/` directory with:

#### **Main Output Files** (summarization numbering: `-1.v1`, `-1.v2`, `-2.v1`):
- **Individual FASTA files**: `{sample_name}-{group}.v{variant}-RiC{reads}.fasta` (all final consensus variants)
- **Combined file**: `summary.fasta` - all final consensus sequences (merged variants only, excludes .raw files)
- **Statistics**: `summary.txt` - sequence counts and metrics
- **Quality report**: `quality_report.txt` - highlights sequences with potential quality concerns
- **Log file**: `summarize_log.txt` - complete processing log

#### **FASTQ Files/** (aggregated reads for final consensus):
- `{sample_name}-{group}.v{variant}-RiC{reads}.fastq` - all reads contributing to final variant consensus

#### **variants/** (pre-merge variant files for traceability):
- `{sample_name}-{group}.v{variant}.raw{N}-RiC{reads}.fasta` - FASTA for pre-merge variant N (when variants were merged)
- **FASTQ Files/**:
  - `{sample_name}-{group}.v{variant}.raw{N}-RiC{reads}.fastq` - reads for pre-merge variant N

### Naming Convention Summary

**Two distinct namespaces maintain traceability:**

| **Namespace** | **Used In** | **Format** | **Purpose** |
|---------------|-------------|------------|-------------|
| **Original** | Source `cluster_debug/` | `-c1`, `-c2`, `-c3` | Preserves speconsense clustering results |
| **Summarization** | `__Summary__/`, `FASTQ Files/`, `variants/` | `-1.v1`, `-1.v2`, `-2.v1`, `.raw1` | Post-processing groups and variants |

### Example Directory Structure
```
__Summary__/
├── sample-1.v1-RiC45.fasta                  # Primary variant (group 1, merged)
├── sample-1.v2-RiC23.fasta                  # Additional variant (not merged)
├── sample-2.v1-RiC30.fasta                  # Second organism group, primary variant
├── summary.fasta                            # All final consensus sequences (excludes .raw)
├── summary.txt                              # Statistics
├── quality_report.txt                       # Quality assessment report
├── summarize_log.txt                        # Processing log
├── FASTQ Files/                             # Reads for final consensus
│   ├── sample-1.v1-RiC45.fastq              # All reads for merged consensus
│   ├── sample-1.v2-RiC23.fastq              # Reads for additional variant
│   └── sample-2.v1-RiC30.fastq              # All reads for second group
└── variants/                                # Pre-merge variant files (for traceability)
    ├── sample-1.v1.raw1-RiC30.fasta         # Pre-merge variant 1 FASTA
    ├── sample-1.v1.raw2-RiC15.fasta         # Pre-merge variant 2 FASTA
    └── FASTQ Files/                         # Pre-merge variant FASTQ files
        ├── sample-1.v1.raw1-RiC30.fastq     # Reads for pre-merge variant 1
        └── sample-1.v1.raw2-RiC15.fastq     # Reads for pre-merge variant 2
```

### FASTA Header Metadata

Consensus sequence headers contain metadata fields separated by spaces:

**Core Fields (Always Present):**
- `size=N` - Total number of raw sequence reads in the cluster
- `ric=N` - **Reads in Consensus** - Number of reads actually used for consensus generation (may be less than size due to sampling limits)

**Optional Fields:**
- `rawric=N+N+...` - RiC values of .raw source variants (pre-merge, largest-first, only present in merged variants from speconsense-summarize)
- `rawlen=N+N+...` - Original sequence lengths before overlap merging (largest-first, only present in overlap-merged variants)
- `snp=N` - Number of SNP positions from IUPAC merging (only present in merged variants from speconsense-summarize)
- `ambig=N` - Count of IUPAC ambiguity codes in consensus (Y, R, W, S, K, M, etc.)
- `length=N` - Sequence length in bases (available via --fasta-fields option)
- `primers=list` - Comma-separated list of detected primer names (e.g., `primers=ITS1F,ITS4`)
- `rid=X.X` - Mean read identity percentage (0-100%, available via --fasta-fields qc preset)
- `rid_min=X.X` - Minimum read identity percentage (worst-case read, available via --fasta-fields qc preset)
- `group=N` - Variant group number (available via --fasta-fields option)
- `variant=vN` - Variant identifier within group (available via --fasta-fields option, only for variants)

**Example Headers:**
```
# Simple consensus from speconsense (debug files):
>sample-c1 size=50 ric=45 rid=97.2 rid_min=94.1 primers=ITS1F,ITS4

# Merged variant from speconsense-summarize (default fields):
>sample-1.v1 size=250 ric=250 rawric=100+89+61 snp=2 ambig=2 primers=ITS1F,ITS4

# Overlap-merged variant (different-length sequences merged):
>sample-1.v1 size=471 ric=471 rawric=248+223 rawlen=630+361 primers=ITS1F,ITS4

# With QC preset (--fasta-fields qc):
>sample-1.v1 size=250 ric=250 length=589 rid=98.5 ambig=1
```

**Notes:**
- Use `--fasta-fields` option to customize which fields appear in output headers (see [Customizing FASTA Header Fields](#customizing-fasta-header-fields))
- Read identity metrics (`rid`, `rid_min`) reflect homopolymer-normalized sequence identity
- Field name changes from earlier versions: `merged_ric` → `rawric`, `p50diff`/`p95diff` → `rid`/`rid_min`
- Variant merging only occurs between sequences with identical primer sets
- `snp` counts positions where IUPAC codes were introduced during merging; `ambig` counts total ambiguity codes in the final sequence

## Scalability Features

For large datasets (thousands of sequences), speconsense uses external tools to accelerate pairwise comparison operations. **Scalability mode is enabled by default** for datasets with 1000+ sequences.

### Configuration

```bash
# Default: scalability enabled for datasets with 1000+ sequences
speconsense input.fastq

# Custom threshold: enable scalability only for datasets with 500+ sequences
speconsense input.fastq --scale-threshold 500

# Disable scalability entirely (use brute-force for all sizes)
speconsense input.fastq --scale-threshold 0

# Same options work for speconsense-summarize
speconsense-summarize --scale-threshold 500 --source clusters/
```

The threshold parameter allows you to run a single command across multiple specimens of varying sizes. Smaller specimens will use brute-force (which has less overhead for small inputs), while larger specimens will use vsearch acceleration.

### Requirements

Scalability mode requires **vsearch** to be installed:

```bash
conda install bioconda::vsearch
```

If vsearch is not available, speconsense will automatically fall back to brute-force computation.

### How It Works

Scalability mode uses a two-stage approach to reduce the O(n^2) pairwise comparison cost:

1. **Fast candidate finding**: vsearch quickly identifies approximate sequence matches using heuristic search
2. **Exact scoring**: Only candidate pairs are scored using the full alignment algorithm

This reduces the computational complexity from O(n^2) to approximately O(n x k), where k is the number of candidates per sequence.

### Performance

Based on benchmarking, the break-even point is approximately 1000 sequences:
- Below 1000 sequences: brute-force is faster (less overhead)
- Above 1000 sequences: vsearch acceleration provides significant speedup
- At 10,000 sequences: ~10x faster than brute-force

### Implementation Notes

- The scalability layer is designed to be backend-agnostic; future versions may support alternative tools (BLAST, minimap2, etc.)
- For summarize's HAC clustering, scalability is most effective when there are many consensus sequences to cluster
- Scalability mode produces identical or near-identical results to brute-force mode

### Concurrency Control

The `--threads` option controls internal parallelism (vsearch, SPOA consensus generation). Use `0` for auto-detect.

**Default behavior differs between tools:**
- `speconsense`: defaults to `--threads 1` (safe for GNU parallel workflows)
- `speconsense-summarize`: defaults to `--threads 0` (auto-detect, since it runs once across all specimens)

```bash
# Many speconsense jobs via parallel - single-threaded by default (safe)
ls *.fastq | parallel speconsense {}

# Single large speconsense job - enable internal parallelism
speconsense large_dataset.fastq --threads 0

# speconsense-summarize auto-detects threads by default
speconsense-summarize --source clusters/
```

## Early Filtering

Speconsense can optionally apply size filtering early in the pipeline (after pre-phasing merge) to skip expensive variant phasing on small clusters that will ultimately be filtered out. This can significantly improve performance for large datasets.

### Options

- `--enable-early-filter`: Apply min-size and min-cluster-ratio filtering after pre-phasing merge. Small clusters skip variant phasing, improving performance.
- `--collect-discards`: Write all discarded reads (outliers + filtered clusters) to `cluster_debug/{sample}-discards.fastq` for inspection.

### Performance Optimization

For large datasets where performance is critical:

```bash
speconsense input.fastq --enable-early-filter
```

This skips variant phasing for clusters that will be filtered out anyway, while still using scalability optimizations for O(n²) comparisons.

## Algorithm Details

### Clustering Methods

Speconsense offers two clustering approaches with different characteristics:

#### **Graph-based clustering with Markov Clustering (MCL)** - Default and recommended
- Constructs a similarity graph between reads and applies the Markov Clustering algorithm to identify clusters
- **Tends to produce more clusters** by discriminating between sequence variants within a specimen
- Suitable when you want to identify multiple variants per specimen and are willing to interpret complex results
- Excellent for detecting subtle sequence differences and biological variation
- Relatively fast with high-quality results for most datasets

**Why MCL?** Markov Clustering was chosen for its ability to detect natural communities in densely connected graphs without requiring a priori specification of cluster numbers. The algorithm simulates flow in the graph, causing flow to spread within natural clusters and evaporate between different clusters. By varying a single parameter (inflation), MCL can detect clusters at different scales of granularity. MCL has proven highly effective in bioinformatics applications, particularly for protein family detection in large-scale sequence databases, making it well-suited for identifying sequence variants in amplicon data.

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
- **Note**: Use `speconsense-summarize` after clustering to manage multiple variants per specimen. Key options include `--merge-position-count` for merging variants differing by few SNPs/indels, `--group-identity` for grouping similar variants, and `--select-max-variants`/`--select-strategy` for controlling which variants to output

### Cluster Merging

After initial clustering, Speconsense automatically merges clusters with identical or homopolymer-equivalent consensus sequences:

**Default behavior (homopolymer-aware merging):**
- Clusters with identical consensus sequences are automatically merged
- Clusters with homopolymer-equivalent sequences are also merged (e.g., "AAA" vs "AAAAA" treated as identical)
- Homopolymer length differences are ignored as they typically represent ONT sequencing artifacts
- This helps eliminate redundant clusters that represent the same biological sequence but differ only in homopolymer lengths

**Strict identity merging (`--disable-homopolymer-equivalence`):**
- Only clusters with perfectly identical consensus sequences are merged
- Use this flag when homopolymer differences may have biological significance
- Results in more clusters but ensures no information loss from homopolymer variation

This automatic merging step helps consolidate redundant clusters that were separated during initial clustering. The adjusted identity scoring with homopolymer normalization provides more accurate assessment of sequence similarity for merging decisions, especially for nanopore data where homopolymer length calling can be inconsistent.

**Variant merging with IUPAC codes:** For more aggressive merging based on SNP and indel thresholds (which creates IUPAC consensus sequences with ambiguity codes), use the `--merge-position-count` and `--merge-indel-length` options in `speconsense-summarize` during post-processing. See the [Advanced Post-Processing](#advanced-post-processing) section.

### Cluster Size Filtering

Speconsense provides two complementary filters to control which clusters are output:

**Absolute size filtering (`--min-size`, default: 5):**
- Filters clusters by absolute number of sequences
- Applied **after merging** identical/homopolymer-equivalent clusters
- Set to 0 to disable and output all clusters regardless of size

**Relative size filtering (`--min-cluster-ratio`, default: 0.01):**
- Filters clusters based on size relative to the largest cluster
- Also applied **after merging** identical/homopolymer-equivalent clusters (post-merge sizes)
- **Based on original cluster sizes**, not sampled sizes from `--max-sample-size`
- Set to 0 to disable and keep all clusters that pass `--min-size`

**Processing order:**
1. Initial clustering produces raw clusters
2. Merge identical/homopolymer-equivalent clusters
3. Filter by `--min-size` (absolute threshold, using post-merge sizes)
4. Filter by `--min-cluster-ratio` (relative threshold, using post-merge sizes)
5. Sample sequences for consensus generation if cluster > `--max-sample-size`

This order ensures that small clusters with identical/homopolymer-equivalent consensus sequences can merge before size filtering is applied. This allows, for example, two clusters of size 3 with identical consensus to merge into a size-6 cluster that passes the `--min-size=5` threshold, rather than being discarded prematurely.

**Deferred filtering strategy:**
For maximum flexibility in detecting rare variants and contaminants, disable filtering in speconsense (`--min-size 0 --min-cluster-ratio 0`) and apply final quality thresholds using `--min-ric` in speconsense-summarize. This allows you to run expensive clustering once and experiment with different quality thresholds during post-processing. However, be aware that permissive filtering may allow more bioinformatic contamination through the pipeline. When using this approach, consider stricter filtering during upstream demultiplexing or perform careful manual review of low-abundance clusters.

### Consensus Generation

Consensus sequences are generated using SPOA (SIMD Partial Order Alignment) which efficiently handles the error profile of nanopore reads. For larger clusters, a random subset of reads (controlled by `--max-sample-size`) is used to generate the consensus.

### Read Identity and Variance Metrics

To evaluate the reliability of consensus sequences, Speconsense calculates read identity metrics by:
1. Aligning all reads in a cluster to the consensus sequence using SPOA multiple sequence alignment
2. Computing per-read identity scores with homopolymer normalization
3. Reporting mean read identity (`rid`) and minimum read identity (`rid_min`)

**Key metrics:**
- `rid` - Mean read identity: The average identity of all reads to the consensus (0-100%). Higher values indicate more homogeneous clusters.
- `rid_min` - Minimum read identity: The worst-case read identity. Low values may indicate outliers or mixed clusters.

**Homopolymer-normalized identity**: Speconsense uses the adjusted-identity algorithm with homopolymer normalization for more accurate sequence comparisons. This means that differences in homopolymer run lengths (e.g., AAA vs AAAAA) are treated as identical, which is particularly important for nanopore sequencing data where homopolymer length calling can be inconsistent. The identity metrics exclude homopolymer length differences, so low `rid` values indicate true substitutions or structural indels rather than sequencing artifacts.

**Automatic variant detection**: By default, Speconsense analyzes positional variation within clusters to detect true biological variants vs sequencing errors. Positions where variant alleles exceed the minimum frequency threshold are flagged as variant positions, and reads are automatically separated into distinct haplotypes (see [Automatic Variant Phasing](#automatic-variant-phasing)).

### Primer Trimming

When a primers file is provided via `--primers`, Speconsense will identify and trim primer sequences from the 5' and 3' ends of consensus sequences, producing clean amplicon sequences for downstream analysis.

**Automatic primer detection**: If `--primers` is not specified, Speconsense will automatically look for `primers.fasta` in the same directory as the input FASTQ file. If found, primer trimming will be enabled automatically.

### Automatic Variant Phasing

By default, Speconsense automatically detects and separates biological variants within clusters. This feature is particularly useful for heterozygous samples or mixed-species amplicons.

**How variant phasing works:**

1. **Variant detection**: After initial clustering, Speconsense analyzes positional variation using multiple sequence alignment. Positions where the minor allele frequency exceeds the threshold (default 10%) and meets minimum read count requirements are identified as variant positions.

2. **Position selection**: When multiple variant positions are detected, Speconsense selects the single best position for splitting (minimizing within-cluster error), then recursively regenerates MSA for each subcluster to discover additional variant positions. This hierarchical approach prevents over-fragmentation while allowing deep phasing when supported by the data.

3. **Haplotype separation**: Reads are grouped by their allele combinations at selected variant positions. Each unique combination becomes a separate sub-cluster (haplotype).

4. **Haplotype filtering**: Small haplotypes that don't meet minimum thresholds are reassigned to the nearest qualifying haplotype, preventing read loss.

5. **IUPAC ambiguity calling**: When only one haplotype qualifies (insufficient support for phasing), variant positions are encoded using IUPAC ambiguity codes (e.g., `Y` for C/T, `R` for A/G) rather than forcing a potentially incorrect consensus call. Ambiguity calling uses a lower count threshold (3 reads vs 5 for splitting), capturing variation that doesn't have enough support to form a viable separate cluster.

**Key parameters for variant phasing:**
- `--min-variant-frequency` - Minimum minor allele frequency to trigger cluster splitting (default: 0.10 = 10%)
- `--min-variant-count` - Minimum read count for minor allele to trigger splitting (default: 5)
- `--disable-position-phasing` - Disable variant phasing entirely

**Key parameters for IUPAC ambiguity calling:**
- `--min-ambiguity-frequency` - Minimum minor allele frequency for IUPAC codes (default: 0.10 = 10%)
- `--min-ambiguity-count` - Minimum read count for minor allele for IUPAC codes (default: 3)
- `--disable-ambiguity-calling` - Disable IUPAC codes for unphased variants

The frequency thresholds are unified at 10%, but the count thresholds differ: splitting requires 5 reads (matching `--min-size`) to ensure viable clusters, while ambiguity calling uses 3 reads since it doesn't create new clusters.

**Example:**
```bash
# Default behavior: variant phasing enabled
speconsense input.fastq

# More permissive variant detection (lower frequency threshold)
speconsense input.fastq --min-variant-frequency 0.05

# More aggressive IUPAC ambiguity calling
speconsense input.fastq --min-ambiguity-frequency 0.05 --min-ambiguity-count 2

# Disable variant phasing but keep ambiguity calling
speconsense input.fastq --disable-position-phasing

# Disable both phasing and ambiguity calling
speconsense input.fastq --disable-position-phasing --disable-ambiguity-calling
```

**When to use variant phasing (default):**
- Heterozygous specimens where you want to separate alleles
- Samples with potential mixed-species amplification
- Quality control to identify clusters with internal heterogeneity

**When to disable variant phasing:**
- Homozygous samples where variation indicates sequencing errors only
- When you prefer merged IUPAC consensus over separate haplotypes
- For faster processing when variant separation is not needed

## Advanced Post-Processing

The `speconsense-summarize` tool provides sophisticated options for managing multiple variants per specimen. This section covers advanced variant handling - for basic usage, see the [Usage](#usage) section above.

### Variant Grouping and Selection

When multiple variants exist per specimen, `speconsense-summarize` first groups similar variants together, then applies selection strategies within each group:

**Variant Grouping:**
```bash
speconsense-summarize --group-identity 0.9
```
- Uses **Hierarchical Agglomerative Clustering (HAC)** to group variants with sequence identity ≥ threshold
- Default threshold is 0.9 (90% identity)
- Variants within each group are considered similar enough to represent the same biological entity
- **Occurs before merging** to separate dissimilar sequences (e.g., contaminants from primary targets)
- Each group will contribute one primary variant plus additional variants based on selection strategy

**Variant Selection (within each group):**

When multiple variants exist per specimen, `speconsense-summarize` offers two strategies for selecting which variants to output from each group:

**Size-based selection (default):**
```bash
speconsense-summarize --select-strategy size --select-max-variants 2
```
- Selects variants by cluster size (largest first)
- Primary variant is always the largest cluster
- Additional variants are selected in order of decreasing size
- Best for identifying the most abundant sequence variants
- Suitable when read count reflects biological abundance
- **Default**: `--select-max-variants=-1` (outputs all variants, no limit)

**Diversity-based selection:**
```bash
speconsense-summarize --select-strategy diversity --select-max-variants 2
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
   - Repeat until select_max_variants reached for this group

**Overall Process:**
1. Group variants by sequence identity using HAC clustering
2. For each group independently:
   - Apply MSA-based merging to find largest compatible subsets
   - Apply variant selection strategy (size or diversity)
   - Output up to select_max_variants per group
3. Final output contains representatives from all groups, ensuring both biological diversity (between groups) and appropriate sampling within each biological entity (within groups)

This two-stage process ensures that distinct biological sequences are preserved as separate groups, while providing control over variant complexity within each group.

### Customizing FASTA Header Fields

Control which metadata fields appear in FASTA headers using the `--fasta-fields` option:

**Presets:**
```bash
# Default preset (current behavior)
speconsense-summarize --fasta-fields default
# Output: >sample-1 size=638 ric=638 rawric=333+305 snp=1 ambig=1 primers=5'-ITS1F,3'-ITS4_RC

# Minimal headers - just the essentials
speconsense-summarize --fasta-fields minimal
# Output: >sample-1 size=638 ric=638

# QC preset - includes read identity and length
speconsense-summarize --fasta-fields qc
# Output: >sample-1 size=638 ric=638 length=589 rid=98.5 ambig=1

# Full metadata (all available fields)
speconsense-summarize --fasta-fields full
# Output: >sample-1 size=638 ric=638 length=589 rawric=333+305 snp=1 ambig=1 rid=98.5 primers=...

# ID only (no metadata)
speconsense-summarize --fasta-fields id-only
# Output: >sample-1
```

**Custom field selection:**
```bash
# Specify individual fields (comma-separated)
speconsense-summarize --fasta-fields size,ric,primers

# Combine presets and fields
speconsense-summarize --fasta-fields minimal,rid

# Combine multiple presets
speconsense-summarize --fasta-fields minimal,qc
```

**Available fields:**
- `size` - Total reads across merged variants
- `ric` - Reads in consensus
- `length` - Sequence length in bases
- `rawric` - RiC values of .raw source variants (only when merged)
- `rawlen` - Original sequence lengths before overlap merging (only when overlap-merged)
- `snp` - Number of IUPAC positions from merging (only when >0)
- `ambig` - Count of IUPAC ambiguity codes in consensus
- `rid` - Mean read identity percentage (when available)
- `rid_min` - Minimum read identity percentage (when available)
- `primers` - Detected primer names (when detected)
- `group` - Variant group number
- `variant` - Variant identifier within group (only for variants)

**Use cases:**
- **Downstream tool compatibility**: Use `minimal` or `id-only` for tools expecting simple headers
- **Quality control**: Use `qc` preset to include read identity metrics for assessing consensus quality
- **File size optimization**: Use `minimal` to reduce file size for large datasets
- **Custom workflows**: Combine presets and fields for workflow-specific needs

### Quality Assessment and Reporting

Speconsense-summarize automatically generates a `quality_report.txt` file to help prioritize manual review of sequences with potential quality concerns. This is particularly valuable for high-throughput workflows where human time for manual inspection is limited.

**Report Generation:**
- Created automatically in the summary output directory
- No configuration required - generated regardless of `--fasta-fields` settings
- Focuses exclusively on sequences that may need review (no "excellent" sequences listed)

**What the Report Includes:**

**1. Statistical Outliers (Low Read Identity):**
- Sequences with mean read identity (`rid`) below the statistical threshold
- Threshold is calculated as: mean - 2×standard deviation across all sequences
- These represent the ~2.5% of sequences with lowest internal consistency
- May indicate mixed clusters, contamination, or problematic consensus

**2. Sequences with IUPAC Ambiguity Codes:**
- Sequences containing IUPAC ambiguity codes (Y, R, W, S, K, M, etc.)
- May result from merging or from unphased heterozygous positions
- Count shown in `ambig` field

**3. Overlap Merge Analysis:**
- Lists specimens where partial overlap merging occurred (different-length sequences with shared overlap region)
- Shows merge iterations, overlap percentages, and prefix/suffix extensions
- Flags edge cases (overlap near threshold, large length ratios >3:1)
- Only includes partial overlaps (prefix or suffix extension >0); excludes full containment merges

**Understanding Read Identity Metrics:**

Read identity assessment works by:
1. Aligning all reads in a cluster to the consensus using SPOA multiple sequence alignment
2. Computing per-read identity with homopolymer normalization
3. Reporting mean identity (`rid`) and minimum identity (`rid_min`)

**Key points about identity calculation:**
- Uses homopolymer normalization - differences like AAA vs AAAAA don't reduce identity
- Therefore, low identity indicates **true substitutions or structural indels**
- `rid` typically ranges from 95-99% for good quality clusters
- `rid_min` shows the worst-case read - very low values may indicate outlier reads

**Example Report Entry:**
```
Type    Sequence                                              RiC    rid      Notes
------------------------------------------------------------------------------------------------
STAT    specimen-1                                            450    93.2     Statistical outlier (rid below threshold)
        specimen-2                                            320    97.8     Contains 2 ambiguity codes
```

**Recommended Actions:**

**For Statistical Outliers (low rid):**
- Review cluster using `cluster_debug/` FASTQ files in source directory
- Check for biological variation (multiple true variants) vs bioinformatic contamination
- Variant phasing is enabled by default; consider adjusting `--min-variant-frequency` if variants weren't separated
- May require manual curation or re-demultiplexing

**For Sequences with Ambiguity Codes:**
- Review whether ambiguity represents true heterozygosity or merging artifact
- Consider adjusting merge parameters if too many variants are being merged
- Check if variant phasing would produce cleaner separate haplotypes

**Workflow Integration:**

The quality report is designed for efficient triage:
1. Statistical outliers are flagged automatically based on global thresholds
2. Review flagged sequences starting with lowest `rid` values
3. Use `rid_min` to identify sequences with problematic individual reads
4. Consider re-processing with variant phasing for mixed clusters

For high-throughput workflows (e.g., 100K sequences/year), this prioritization ensures human review time focuses on the most actionable quality issues.

### Additional Summarize Options

**Quality Filtering:**
```bash
speconsense-summarize --min-ric 5
```
- Filters out consensus sequences with fewer than the specified number of Reads in Consensus (RiC)
- Default is 3 - only sequences supported by at least 3 reads are processed
- Higher values provide more stringent quality control but may exclude valid low-abundance variants

**Length Filtering:**
```bash
speconsense-summarize --min-len 400 --max-len 800
```
- `--min-len`: Minimum sequence length in bp (default: 0 = disabled)
- `--max-len`: Maximum sequence length in bp (default: 0 = disabled)
- Applied during initial sequence loading, before HAC grouping
- Useful for filtering incomplete amplicons or chimeric sequences

**Variant Merging:**
```bash
# Basic SNP-only merging (default)
speconsense-summarize --merge-position-count 2

# Enable indel merging (up to 3bp indels)
speconsense-summarize --merge-position-count 3 --merge-indel-length 3

# Disable SNP merging (only merge identical sequences)
speconsense-summarize --merge-snp false

# Disable all merging (fastest, skip MSA evaluation entirely)
speconsense-summarize --disable-merging

# Control merge evaluation effort (performance vs thoroughness)
speconsense-summarize --merge-effort fast      # Faster, may miss some merges in large groups
speconsense-summarize --merge-effort balanced  # Default
speconsense-summarize --merge-effort thorough  # More exhaustive for large variant groups

# Legacy parameter (still supported)
speconsense-summarize --snp-merge-limit 2  # Equivalent to --merge-position-count 2
```
- **Occurs within each HAC group** - merges variants that differ by small numbers of SNPs and/or indels
- Uses **MSA-based approach with SPOA**: evaluates all possible subsets to find the largest compatible group for merging
- Creates **IUPAC consensus sequences** with size-weighted majority voting at polymorphic positions
- **Order-independent**: produces identical results regardless of input order (unlike pairwise greedy approaches)
- Only merges variants with **identical primer sets** to maintain biological validity

**Merge Parameters:**
- `--disable-merging`: Skip MSA-based merge evaluation entirely (fastest option when merging is not needed)
- `--merge-effort LEVEL`: Control merge evaluation thoroughness. Presets: `fast` (8), `balanced` (10, default), `thorough` (12), or numeric 6-14. Higher values use larger batch sizes for exhaustive subset search, improving merge quality for large variant groups at the cost of runtime
- `--merge-position-count N`: Maximum total SNP + structural indel positions allowed (default: 2)
- `--merge-indel-length N`: Maximum length of individual structural indels allowed (default: 0 = disabled)
- `--merge-snp`: Enable/disable SNP merging (default: True)
- `--merge-min-size-ratio R`: Minimum size ratio (smaller/larger) for merging clusters (default: 0.1, 0 to disable)
- `--disable-homopolymer-equivalence`: Treat homopolymer length differences as structural indels (default: disabled, meaning homopolymer equivalence is enabled)
- `--snp-merge-limit N`: Legacy parameter, equivalent to `--merge-position-count` (deprecated)

**Homopolymer-Aware Merging (Default Behavior):**

By default, speconsense-summarize uses **homopolymer-aware merging** that distinguishes between structural indels (true insertions/deletions) and homopolymer length differences (e.g., AAA vs AAAA). This matches the semantics of adjusted-identity alignment used throughout the pipeline.

**How it works:**
- Analyzes SPOA multiple sequence alignment to classify each indel column
- **Homopolymer indel**: All non-gap bases are the same, and all sequences agree on that base in adjacent columns
  - Example: `ATA-AAGC` vs `ATAAAAGC` → homopolymer (all have A's flanking the gap)
- **Structural indel**: True insertion/deletion, or indel adjacent to SNP
  - Example: `CTAA-GC` vs `CTG-AGC` → structural (A vs G flanking the gap)
- **Homopolymer indels are ignored** when checking merge compatibility (treated as equivalent)
- **Structural indels count** against `--merge-position-count` and `--merge-indel-length` limits

**Default behavior examples:**
- Sequences differing only by homopolymer length (AAA vs AAAA): **Merge** ✓
- Sequences with 2 SNPs + 5 homopolymer indels: **Merge** ✓ (only SNPs count)
- Sequences with 1 structural indel: **Blocked** (default `--merge-indel-length=0`)

**Strict identity merging (`--disable-homopolymer-equivalence`):**
```bash
speconsense-summarize --disable-homopolymer-equivalence
```
- Treats homopolymer length differences as structural indels
- All indels (both homopolymer and structural) count against merge limits
- Only sequences identical (or differing by SNPs if `--merge-snp` enabled) can merge
- Use when you want to preserve homopolymer length variation as distinct variants

**Example comparison:**
```bash
# Default (homopolymer equivalence enabled)
# Variants: ATAAAGC (3 A's) vs ATAAAAGC (4 A's)
speconsense-summarize
# Result: Merge ✓ (homopolymer difference ignored)

# Strict mode (homopolymer equivalence disabled)
speconsense-summarize --disable-homopolymer-equivalence
# Result: Do not merge (treated as structural indel)
```

**Position counting examples:**

With default settings (`--merge-position-count 2 --merge-indel-length 0`):
- 2 SNPs + 0 structural indels + 10 homopolymer indels ✓ (only structural counted)
- 0 SNPs + 1 structural indel + 5 homopolymer indels ✗ (structural indel blocked)
- 1 SNP + 1 structural indel (≤2bp) ✗ (indel blocked by length limit)

With indels enabled (`--merge-position-count 3 --merge-indel-length 2`):
- 3 SNPs + 0 structural indels + any homopolymer indels ✓
- 2 SNPs + 1 structural indel (≤2bp) + any homopolymer indels ✓
- 0 SNPs + 3 structural indels (each ≤2bp) + any homopolymer indels ✓
- 2 SNPs + 2 structural indels (each ≤2bp) ✗ (total=4 > 3)
- 2 SNPs + 1 structural indel (3bp) ✗ (indel too long)

**Merged consensus tracking:**
- Merged sequences include `rawric` header field showing RiC values of merged .raw variants
- Example: `>sample-1 size=250 ric=250 rawric=100+89+61 snp=2`
- Helps trace which original clusters contributed to merged consensus

**Note**: Homopolymer-aware merging in speconsense-summarize complements the automatic homopolymer-equivalent cluster merging that occurs during the main clustering step in speconsense

**Understanding MSA-Based Merge Evaluation:**

A key architectural feature of speconsense-summarize is that **all sequences within a HAC group are aligned together by SPOA first**, creating a single multi-sequence alignment (MSA) that provides biological and evolutionary context for all subsequent merge decisions within that group.

**Why this matters:**

When evaluating whether two sequences can merge, the algorithm doesn't perform a simple pairwise alignment. Instead, it examines how those sequences align within the context of all other sequences in their HAC group. This shared alignment context can reveal that apparent structural differences are actually homopolymer length variations relative to the majority pattern.

**Example:**

Consider two sequences (c4 and c9) from the same HAC group. When aligned in isolation, they appear incompatible:

```
Pairwise alignment (c4 vs c9 only):
c4: ..504..GTAAATT..242..AG
c9: ..504..G-GTACT..242..AA
```

This pairwise alignment shows 4 SNPs + 1 structural indel, making them incompatible for merging with default parameters.

However, when aligned with all sequences in their HAC group, a different picture emerges:

```
Multi-sequence alignment (c1, c2, c4, c5, c7, c8, c9):
c1: ..153..AAA..57..AGC..151..CTCTT..131..A-GTAAATT.TCA..213..TCAG..21..AA
c2: ..153..ATA..57..AAC..151..CTCTT..131..A-GTAAATT.TCA..213..TCAG..21..AA
c4: ..153..AAA..57..AGC..151..CTCTT..131..A-GTAAATT.TCA..213..TCAG..21..AG
c5: ..153..AAA..57..AGC..151..C---T..131..A-GTAAATT.TCA..213..TCAG..21..AA
c7: ..153..AAA..57..AGC..151..CTCTT..131..A-GT-ACTT.T-A..213..TCAG..21..AA
c8: ..153..AAA..57..AGC..151..CTCTT..131..A-GTAAATT.TCA..213..TAGG..21..AA
c9: ..153..AAA..57..AGC..151..CTCTT..131..AGGT-AC-T.TCA..213..TCAG..21..AA
```

In this multi-sequence context, the majority pattern (seen in c1, c2, c4, c5, c8) is `A-GTAAATT` in the critical region around position 131. The c9 sequence (`AGGT-AC-T`) differs from c4, but the gaps occur at different positions in the homopolymer-rich region. The multi-sequence alignment reveals that c9's differences are better explained as homopolymer length variations (extra G's and different gap placements in poly-A and poly-T regions) rather than true structural changes. This results in only 2 SNPs + 3 homopolymer indels, making them compatible for merging.

This biological interpretation, informed by the evolutionary context of multiple related sequences, produces more accurate merge decisions than pairwise comparisons alone.

**Practical implication:**

Merge decisions depend on the complete HAC group composition, not just the sequences being evaluated. Two sequences that appear incompatible in isolation may merge when analyzed in their full biological context, and vice versa. This is by design—the multi-sequence context provides more accurate biological interpretation of sequence variation.

**Merge Size Ratio Filtering:**
```bash
speconsense-summarize --merge-min-size-ratio 0.1
```
- Prevents merging clusters with very different sizes (e.g., well-supported variant + poorly-supported variant)
- Ratio calculated as `smaller_size / larger_size` - must be ≥ threshold to merge
- Example: `--merge-min-size-ratio 0.1` means smaller cluster must be ≥10% size of larger
- Default is 0.1 - smaller cluster must be ≥10% size of larger to merge
- **Use cases:**
  - Prevent poorly-supported variants (low read count) from introducing ambiguities into well-supported sequences
  - Avoid merging weakly-supported bioinformatic artifacts into high-confidence sequences
- Applied during MSA-based merging step, within each HAC group

**Group Output Limiting:**
```bash
speconsense-summarize --select-max-groups 2
```
- Limits output to top N groups per specimen (by size of largest member in each group)
- Default is -1 (output all groups)
- Applied after HAC clustering but before MSA-based merging and variant selection
- Useful for focusing on primary specimens and ignoring small contaminant groups
- Example: `--select-max-groups=1` outputs only the largest variant group per specimen

**Directory Control:**
```bash
speconsense-summarize --source /path/to/speconsense/output --summary-dir MyResults
```
- `--source`: Directory containing speconsense output files (default: clusters)
- `--summary-dir`: Output directory name (default: `__Summary__`)

### Overlap Merging for Primer Pools

When using multiple primers targeting the same locus ("primer pools"), reads may have overlapping but different coverage. The overlap merge feature (enabled by default) allows merging such sequences when they share sufficient overlap.

**How it works:**
1. During HAC clustering, uses single-linkage to group overlapping sequences
2. Identifies sequences with sufficient overlap meeting identity threshold
3. Creates consensus from union of coverage (overlap region uses majority voting)
4. Supports iterative merging for 3+ overlapping sequences
5. **Primer-constrained**: Overlap-aware distance only applies when sequences have different primer pairs (legitimate primer pool scenarios). Same-primer sequences use global distance to prevent chimeras from incorrectly merging with shorter amplicons.

**Parameters:**
- `--min-merge-overlap N`: Minimum overlap in bp (default: 200, 0 to disable)
- `--group-identity`: Identity threshold for HAC grouping, also used for overlap region identity (default: 0.9)

**Example:**
```bash
# Default behavior (overlap merging enabled)
speconsense-summarize --source clusters

# Disable overlap merging (original behavior)
speconsense-summarize --source clusters --min-merge-overlap 0

# More permissive overlap (allow smaller overlaps)
speconsense-summarize --source clusters --min-merge-overlap 100
```

**Use cases:**
- ITS2 sequence merging with full ITS sequence
- Overlapping amplicons from primer pools
- Partial sequences merging with complete references

**Output indicators:**
- Log messages show `(overlap=Xbp, prefix=Ybp, suffix=Zbp)` for overlap merges
- FASTA headers include `rawlen=X+Y` showing original sequence lengths
- Quality report includes "OVERLAP MERGE ANALYSIS" section with details on partial overlaps

**Containment handling:**
When a shorter sequence is fully contained within a longer one (e.g., ITS2 within full ITS), the merge is allowed if `overlap >= min(threshold, shorter_length)`.

### Processing Workflow Summary

The complete speconsense-summarize workflow operates in this order:

1. **Load sequences** with RiC filtering (`--min-ric`) and length filtering (`--min-len`, `--max-len`)
2. **HAC variant grouping** by sequence identity to separate dissimilar sequences (`--group-identity`); uses single-linkage when overlap merging is enabled
3. **Group filtering** to limit output groups (`--select-max-groups`)
4. **Homopolymer-aware MSA-based variant merging** within each group, including **overlap merging** for different-length sequences (`--disable-merging`, `--merge-effort`, `--merge-position-count`, `--merge-indel-length`, `--min-merge-overlap`, `--merge-snp`, `--merge-min-size-ratio`, `--disable-homopolymer-equivalence`)
5. **Variant selection** within each group (`--select-max-variants`, `--select-strategy`)
6. **Output generation** with customizable header fields (`--fasta-fields`) and full traceability

**Key architectural features**:
- HAC grouping occurs BEFORE merging to prevent inappropriate merging of dissimilar sequences (e.g., contaminants with primary targets)
- Merging is applied independently within each group using MSA-based consensus generation
- Homopolymer-aware merging by default (AAA ≈ AAAA) to match pipeline-wide adjusted-identity semantics
- Overlap merging enabled by default (`--min-merge-overlap 200`) for primer pool support

### Enhanced Logging and Traceability

Speconsense-summarize provides comprehensive logging to help users understand processing decisions:

**Variant Analysis Logging:**
- **Complete variant summaries** for every variant in each group, including those that are skipped
- **Detailed difference categorization**: substitutions, single-nt indels, short (≤3nt) indels, and long indels
- **IUPAC-aware comparisons**: treats ambiguity codes as matches (e.g., R matches A or G)
- **Homopolymer-aware merge reporting**: distinguishes structural indels from homopolymer indels
- **Group context**: clearly shows which variants belong to each HAC clustering group
- **Selection rationale**: explains why variants were included or excluded

**Example log output:**
```
HAC clustering created 2 groups
Group 1: ['sample-c3']
Group 2: ['sample-c1', 'sample-c2']
Found mergeable subset of 2 variants: 2 SNPs, 1 homopolymer indels
Processing Variants in Group 2
Primary: sample-c1 (size=403, ric=403)
Variant 1: (size=269, ric=269) - 1 short (<= 3nt) indel
Variant 2: (size=180, ric=180) - 3 substitutions, 1 single-nt indel - skipping
```

**Traceability Features:**
- **Merge history**: tracks which original clusters were combined during variant merging
- **File lineage**: maintains connection between final outputs and original speconsense clusters
- **Read aggregation**: `FASTQ Files/` directory contains all reads that contributed to each final consensus
- **Pre-merge preservation**: `variants/` directory contains `.raw` files that preserve individual pre-merge variants with their original sequences and reads
- **Cluster boundaries in merged FASTQ**: When multiple clusters are merged, synthetic delimiter records mark boundaries between clusters for easy identification in sequence viewers (e.g., UGENE). Format: `@CLUSTER_BOUNDARY_{n}:{cluster}:RiC={ric}:reads={count}`

This comprehensive logging allows users to understand exactly how the pipeline processed their data and make informed decisions about parameter tuning.

## Full Command Line Options

```
usage: speconsense [-h] [-O OUTPUT_DIR] [--primers PRIMERS]
                   [--augment-input AUGMENT_INPUT]
                   [--algorithm {graph,greedy}] [--min-identity MIN_IDENTITY]
                   [--inflation INFLATION]
                   [--k-nearest-neighbors K_NEAREST_NEIGHBORS]
                   [--min-size MIN_SIZE]
                   [--min-cluster-ratio MIN_CLUSTER_RATIO]
                   [--max-sample-size MAX_SAMPLE_SIZE]
                   [--outlier-identity OUTLIER_IDENTITY]
                   [--disable-position-phasing]
                   [--min-variant-frequency MIN_VARIANT_FREQUENCY]
                   [--min-variant-count MIN_VARIANT_COUNT]
                   [--disable-ambiguity-calling]
                   [--min-ambiguity-frequency MIN_AMBIGUITY_FREQUENCY]
                   [--min-ambiguity-count MIN_AMBIGUITY_COUNT]
                   [--disable-cluster-merging]
                   [--disable-homopolymer-equivalence]
                   [--orient-mode {skip,keep-all,filter-failed}]
                   [--presample PRESAMPLE] [--scale-threshold SCALE_THRESHOLD]
                   [--threads N] [--enable-early-filter] [--collect-discards]
                   [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                   [--version] [-p NAME] [--list-profiles]
                   input_file

MCL-based clustering of nanopore amplicon reads

options:
  -h, --help            show this help message and exit
  --version             Show program's version number and exit
  -p NAME, --profile NAME
                        Load parameter profile (use --list-profiles to see
                        available)
  --list-profiles       List available profiles and exit

Input/Output:
  input_file            Input FASTQ file
  -O OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for all files (default: clusters)
  --primers PRIMERS     FASTA file containing primer sequences (default: looks
                        for primers.fasta in input file directory)
  --augment-input AUGMENT_INPUT
                        Additional FASTQ/FASTA file with sequences recovered
                        after primary demultiplexing (e.g., from specimine)

Clustering:
  --algorithm {graph,greedy}
                        Clustering algorithm to use (default: graph)
  --min-identity MIN_IDENTITY
                        Minimum sequence identity threshold for clustering
                        (default: 0.9)
  --inflation INFLATION
                        MCL inflation parameter (default: 4.0)
  --k-nearest-neighbors K_NEAREST_NEIGHBORS
                        Number of nearest neighbors for graph construction
                        (default: 5)

Filtering:
  --min-size MIN_SIZE   Minimum cluster size (default: 5, 0 to disable)
  --min-cluster-ratio MIN_CLUSTER_RATIO
                        Minimum size ratio between a cluster and the largest
                        cluster (default: 0.01, 0 to disable)
  --max-sample-size MAX_SAMPLE_SIZE
                        Maximum cluster size for consensus (default: 100)
  --outlier-identity OUTLIER_IDENTITY
                        Minimum read-to-consensus identity to keep a read
                        (default: auto). Reads below this threshold are
                        removed as outliers before final consensus generation.
                        Auto-calculated as (1 + min_identity) / 2. This
                        threshold is typically higher than --min-identity
                        because the consensus is error-corrected through
                        averaging.

Variant Phasing:
  --disable-position-phasing
                        Disable position-based variant phasing (enabled by
                        default). MCL graph clustering already separates most
                        variants; this second pass analyzes MSA positions to
                        phase remaining variants.
  --min-variant-frequency MIN_VARIANT_FREQUENCY
                        Minimum alternative allele frequency to call variant
                        (default: 0.10 for 10%)
  --min-variant-count MIN_VARIANT_COUNT
                        Minimum alternative allele read count to call variant
                        (default: 5)

Ambiguity Calling:
  --disable-ambiguity-calling
                        Disable IUPAC ambiguity code calling for unphased
                        variant positions
  --min-ambiguity-frequency MIN_AMBIGUITY_FREQUENCY
                        Minimum alternative allele frequency for IUPAC
                        ambiguity calling (default: 0.10 for 10%)
  --min-ambiguity-count MIN_AMBIGUITY_COUNT
                        Minimum alternative allele read count for IUPAC
                        ambiguity calling (default: 3)

Cluster Merging:
  --disable-cluster-merging
                        Disable merging of clusters with identical consensus
                        sequences
  --disable-homopolymer-equivalence
                        Disable homopolymer equivalence in cluster merging
                        (only merge identical sequences)

Orientation:
  --orient-mode {skip,keep-all,filter-failed}
                        Sequence orientation mode: skip (default, no
                        orientation), keep-all (orient but keep failed), or
                        filter-failed (orient and remove failed)

Performance:
  --presample PRESAMPLE
                        Presample size for initial reads (default: 1000, 0 to
                        disable)
  --scale-threshold SCALE_THRESHOLD
                        Sequence count threshold for scalable mode (requires
                        vsearch). Set to 0 to disable. Default: 1001
  --threads N           Max threads for internal parallelism (vsearch, SPOA).
                        0=auto-detect, default=1 (safe for parallel
                        workflows).
  --enable-early-filter
                        Enable early filtering to skip small clusters before
                        variant phasing (improves performance for large
                        datasets)

Debugging:
  --collect-discards    Write discarded reads (outliers and filtered clusters)
                        to cluster_debug/{sample}-discards.fastq
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
```

### speconsense-summarize Options

```
usage: speconsense-summarize [-h] [--source SOURCE]
                             [--summary-dir SUMMARY_DIR]
                             [--fasta-fields FASTA_FIELDS] [--min-ric MIN_RIC]
                             [--min-len MIN_LEN] [--max-len MAX_LEN]
                             [--group-identity GROUP_IDENTITY] [--merge-snp]
                             [--merge-indel-length MERGE_INDEL_LENGTH]
                             [--merge-position-count MERGE_POSITION_COUNT]
                             [--merge-min-size-ratio MERGE_MIN_SIZE_RATIO]
                             [--min-merge-overlap MIN_MERGE_OVERLAP]
                             [--disable-homopolymer-equivalence]
                             [--select-max-groups SELECT_MAX_GROUPS]
                             [--select-max-variants SELECT_MAX_VARIANTS]
                             [--select-strategy {size,diversity}]
                             [--scale-threshold SCALE_THRESHOLD] [--threads N]
                             [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                             [--version] [-p NAME] [--list-profiles]

Process Speconsense output with advanced variant handling.

options:
  -h, --help            show this help message and exit
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level
  --version             Show program's version number and exit
  -p NAME, --profile NAME
                        Load parameter profile (use --list-profiles to see
                        available)
  --list-profiles       List available profiles and exit

Input/Output:
  --source SOURCE       Source directory containing Speconsense output
                        (default: clusters)
  --summary-dir SUMMARY_DIR
                        Output directory for summary files (default:
                        __Summary__)
  --fasta-fields FASTA_FIELDS
                        FASTA header fields to output. Can be: (1) a preset
                        name (default, minimal, qc, full, id-only), (2) comma-
                        separated field names (size, ric, length, rawric, snp,
                        rid, rid_min, primers, group, variant), or (3) a
                        combination of presets and fields (e.g., minimal,qc or
                        minimal,rid). Duplicates removed, order preserved left
                        to right. Default: default

Filtering:
  --min-ric MIN_RIC     Minimum Reads in Consensus (RiC) threshold (default:
                        3)
  --min-len MIN_LEN     Minimum sequence length in bp (default: 0 = disabled)
  --max-len MAX_LEN     Maximum sequence length in bp (default: 0 = disabled)

Grouping:
  --group-identity GROUP_IDENTITY, --variant-group-identity GROUP_IDENTITY
                        Identity threshold for variant grouping using HAC
                        (default: 0.9)

Merging:
  --disable-merging     Disable all variant merging (skip MSA-based merge
                        evaluation entirely)
  --merge-effort LEVEL  Merging effort level: fast (8), balanced (10),
                        thorough (12), or numeric 6-14. Higher values allow
                        larger batch sizes for exhaustive subset search.
                        Default: balanced
  --merge-snp, --no-merge-snp
                        Enable SNP-based merging (default: True, use --no-
                        merge-snp to disable)
  --merge-indel-length MERGE_INDEL_LENGTH
                        Maximum length of individual indels allowed in merging
                        (default: 0 = disabled)
  --merge-position-count MERGE_POSITION_COUNT
                        Maximum total SNP+indel positions allowed in merging
                        (default: 2)
  --merge-min-size-ratio MERGE_MIN_SIZE_RATIO
                        Minimum size ratio (smaller/larger) for merging
                        clusters (default: 0.1, 0 to disable)
  --min-merge-overlap MIN_MERGE_OVERLAP
                        Minimum overlap in bp for merging sequences of
                        different lengths (default: 200, 0 to disable)
  --disable-homopolymer-equivalence
                        Disable homopolymer equivalence in merging (treat AAA
                        vs AAAA as different)

Selection:
  --select-max-groups SELECT_MAX_GROUPS, --max-groups SELECT_MAX_GROUPS
                        Maximum number of groups to output per specimen
                        (default: -1 = all groups)
  --select-max-variants SELECT_MAX_VARIANTS, --max-variants SELECT_MAX_VARIANTS
                        Maximum total variants to output per group (default:
                        -1 = no limit, 0 also means no limit)
  --select-strategy {size,diversity}, --variant-selection {size,diversity}
                        Variant selection strategy: size or diversity
                        (default: size)

Performance:
  --scale-threshold SCALE_THRESHOLD
                        Sequence count threshold for scalable mode in HAC
                        clustering (requires vsearch). Set to 0 to disable.
                        Default: 1001
  --threads N           Max threads for internal parallelism. 0=auto-detect
                        (default), N>0 for explicit count.
```

## Specialized Workflows

The following features support less common workflows and specialized use cases.

### Sequence Orientation Normalization

**Use Case:** Reprocessing output from older bioinformatics pipelines (such as minibar) that do not automatically normalize sequence orientation.

**Note:** When using speconsense downstream of specimux (recommended workflow), orientation normalization is unnecessary as specimux automatically orients all sequences during demultiplexing.

The `--orient-mode` parameter enables automatic detection and correction of sequence orientation based on primer positions:

```bash
# Keep all sequences, including those with ambiguous orientation
speconsense input.fastq --primers primers.fasta --orient-mode keep-all

# Filter out sequences that couldn't be reliably oriented
speconsense input.fastq --primers primers.fasta --orient-mode filter-failed

# Skip orientation (default, appropriate when using specimux output)
speconsense input.fastq --primers primers.fasta
```

**Available modes:**
- `skip` (default): No orientation performed, sequences processed as-is
- `keep-all`: Perform orientation but keep all sequences, including those that failed
- `filter-failed`: Perform orientation and remove sequences with failed/ambiguous orientation

**How it works:**
- Detects orientation by checking for forward and reverse primers at expected positions
- Uses binary scoring: +1 for forward primer match, +1 for reverse primer match
- Sequences with clear orientation (score >0 in one direction, 0 in the other) are reoriented if needed
- Sequences with ambiguous orientation (both 0 or both >0) are marked as failed

**Requirements:**
- Primers FASTA file must include position annotations: `position=forward` or `position=reverse`
- Example format:
  ```
  >ITS1F  pool=ITS    position=forward
  CTTGGTCATTTAGAGGAAGTAA
  >ITS4   pool=ITS    position=reverse
  TCCTCCGCTTATTGATATGC
  ```

**Technical details:**
- Orientation occurs before clustering, ensuring all sequences are in the same direction
- Failed orientations typically indicate: no primers found, chimeric sequences, or degraded primers
- Quality scores are reversed when sequences are reverse-complemented

### Augmenting Clusters with Recovered Sequences

**Use Case:** Increasing cluster abundance by including sequences that were partially demultiplexed or recovered through mining tools like `specimine`.

The `--augment-input` parameter allows you to supplement primary demultiplexing results with additional sequences recovered from unmatched reads:

```bash
# 1. Run primary demultiplexing with specimux
specimux primers.fasta specimens.txt input.fastq -O results/

# 2. Recover additional sequences from unmatched reads
specimine results/unknown/ specimen_name --output recovered.fastq

# 3. Cluster with both primary and recovered sequences
speconsense results/full/specimen_name.fastq --augment-input recovered.fastq
```

**How it works:**
- Augmented sequences participate equally in clustering with primary sequences
- During presampling, primary sequences are prioritized to ensure representative sampling
- All sequences (primary + augmented) contribute to final consensus generation
- Final output headers show total read counts including augmented sequences

**Key features:**
- Supports both FASTQ and FASTA formats (auto-detected by file extension)
- Augmented sequences fully traceable through the entire pipeline
- Sequences only cluster together if they meet the similarity threshold
- Recovered sequences may form separate clusters if they represent different taxa
- Maximizes data utilization by including sequences that would otherwise be discarded

**Typical workflow:**
After primary demultiplexing, some sequences may remain unmatched due to sequencing errors, primer degradation, or edge cases in barcode detection. Mining tools like `specimine` can recover these sequences based on sequence composition or other characteristics, allowing them to be included in consensus generation and increase cluster support.

### Testing with Synthetic Data

For empirical testing of consensus quality, variant detection, contamination scenarios, and understanding toolchain behavior with controlled datasets, see [Testing with Speconsense-Synth](docs/synthetic-testing.md).

## Future Enhancements

The following features are under consideration for future development:

- **Background contamination detection tool**: A utility to identify contamination patterns affecting multiple specimens within the same sequencing run, helping to distinguish systematic contamination from genuine biological sequences

- **Scalability improvements**: Algorithm enhancements to enable graph-based clustering to efficiently handle datasets with hundreds of thousands of sequences while maintaining accuracy

These features are being explored based on user feedback. Implementation timelines and feasibility are still being evaluated.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for detailed version history.

## Citations

This project uses and builds upon:

- **Markov Clustering (MCL) algorithm**:
  - van Dongen, S. (2000). *Graph Clustering by Flow Simulation*. PhD thesis, University of Utrecht. https://micans.org/mcl/
  - van Dongen, S. (2008). *Graph clustering via a discrete uncoupling process*. SIAM Journal on Matrix Analysis and Applications 30(1):121-141. https://doi.org/10.1137/040608635

- **MCL in bioinformatics**: Enright, A.J., Van Dongen, S., Ouzounis, C.A. (2002). *An efficient algorithm for large-scale detection of protein families*. Nucleic Acids Research 30(7):1575-1584. https://doi.org/10.1093/nar/30.7.1575 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC101833/)

- **ONT fungal barcoding protocol**: Russell, S.D., Geurin, Z., Walker, J. (2024). *Primary Data Analysis - Basecalling, Demultiplexing, and Consensus Building for ONT Fungal Barcodes*. protocols.io. https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v4

## Contributing

Contributions are welcome! For development setup, testing guidelines, and contribution workflow, please see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

[BSD 3-Clause License](LICENSE)

## Name

Speconsense is a playful portmanteau of "Specimen" and "Consensus", reflecting the tool's focus on generating high-quality consensus sequences from specimen amplicon reads.
