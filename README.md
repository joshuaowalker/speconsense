# Speconsense

A tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads.

## Overview

Speconsense is a specialized tool for generating high-quality consensus sequences from Oxford Nanopore Technology (ONT) amplicon data. It is specifically designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

The key features of Speconsense include:
- Robust clustering of amplicon reads using either Markov Clustering (MCL) graph-based or greedy algorithms
- **Automatic variant phasing**: Detects variant positions within clusters and splits reads into separate haplotypes
- **Post-phasing refinement**: Reads can move between clusters within an identity group, and previously-discarded reads can be re-admitted, when consensus concordance supports it
- **Statistical variant validation (CER)**: Each cluster is annotated with a context-aware Critical Error Rate factor measuring how implausible it is as basecaller-noise replicates of a larger peer
- **Cluster homogeneity (err_factor)**: A peer-independent metric for whether a cluster's reads are internally consistent with the assumed sequencing-noise profile
- Automatic merging of clusters with identical or homopolymer-equivalent consensus sequences
- High-quality consensus generation using SPOA, with MAD-based outlier removal at final consensus
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
- `speconsense-fit-error-model` - Fit a custom q_ctx error model from a finished speconsense run (see [Fitting a Custom Error Model](#fitting-a-custom-error-model))

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
- `compressed` — Compress variants into minimal IUPAC consensus sequences (aggressive merging with indels, 20% thresholds, 20% selection size ratio)
- `herbarium` — High-recall for degraded DNA/type specimens (CER and err_factor filters disabled to surface weak-signal variants for review)
- `largedata` — Experimental settings for large input files (tighter `group-identity` to keep identity groups small)
- `nostalgia` — Simulate older bioinformatics pipelines (disables read reassignment, discard recovery, CER and err_factor filters)
- `strict` — High-precision for confident results (per-group `select-min-size-ratio`)

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

Speconsense can detect and isolate sequence variants within specimens (aka "haplotype phasing"). The graph-based clustering algorithm excels at discriminating between variants, and a multi-phase post-processing pipeline (read reassignment, discard recovery, second phasing pass, CER validation) refines the result. Core groups every variant into a complete-linkage **identity group** (gid) and assigns each variant a within-group rank (vid); these labels travel through to summarize via FASTA-header fields. `speconsense-summarize` honors that grouping rather than re-clustering, and provides:
- MSA-based variant merging with IUPAC ambiguity codes and size-weighted consensus, applied within each identity group
- Cross-primer overlap merging across identity groups (for primer-pool workflows)
- CER and err_factor filtering to route low-confidence variants to a separate `.ns`/`.lq` track
- Size-based variant selection with configurable limits per group

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

Filenames use the identity group rank (`gid`) and within-group variant rank (`vid`) end-to-end: `{sample_name}-{gid}.v{vid}-RiC{size}`. There is no separate "original" namespace.

For each specimen, Speconsense generates:

1. **Main consensus FASTA**: `{sample_name}-all.fasta`
   - Contains all consensus sequences for the specimen (one per cluster), tagged with `gid=`/`vid=` and `cer_factor=`/`err_factor=` headers

2. **Debug directory** (`cluster_debug/`):
   - `{sample_name}-{gid}.v{vid}-RiC{size}-reads.fastq`: Original reads in each cluster
   - `{sample_name}-{gid}.v{vid}-RiC{size}-sampled.fastq`: Sampled reads used for consensus generation
   - `{sample_name}-{gid}.v{vid}-RiC{size}-untrimmed.fasta`: Untrimmed consensus sequences
   - `{sample_name}-discards.fastq`: Reads discarded across the pipeline (only when `--collect-discards` is set)

3. **Metadata directory** (`metadata/`):
   - `{sample_name}_metadata.json`: Full per-cluster diagnostics including pipeline phase counts, MAD outlier drops, per-cluster `cer_details` (peer id, K, p*, joint q, per-position context tags and q_ctx values, reference idx) for replaying CER decisions offline, and raw `obs_sum`/`exp_sum`/`cols` for `err_factor` reproduction

### Speconsense-Summarize Output

`speconsense-summarize` honors core's identity groups (parsed from `gid=`/`vid=` headers) and applies cross-primer overlap merging, MSA-based variant merging within each group, CER and err_factor filtering, and variant selection. It refuses to process FASTA inputs that lack `gid=`/`vid=` headers — re-run `speconsense` if upgrading from a pre-0.8 output directory.

Creates the `__Summary__/` directory with:

#### **Main Output Files** (passing variants):
- **Individual FASTA files**: `{sample_name}-{gid}.v{vid}-RiC{reads}.fasta` (one per passing variant)
- **Combined file**: `summary.fasta` — all passing final consensus sequences (excludes `.raw`, `.ns`, and `.lq` files)
- **Statistics**: `summary.txt` — sequence counts and metrics
- **Quality report**: `quality_report.txt` — inspection-oriented report (executive summary, passed variants worth inspecting, possible rescues, overlap merge details, run-wide parameter signals, pipeline activity)
- **Log file**: `summarize_log.txt` — complete processing log

#### **FASTQ Files/** (reads for passing consensus):
- `{sample_name}-{gid}.v{vid}-RiC{reads}.fastq`

#### **variants/** (filtered + pre-merge variants for traceability):
- `{sample_name}-{gid}.v{vid}.raw{N}-RiC{reads}.fasta` — FASTA for pre-merge variant N (when variants were merged)
- `{sample_name}-{gid}.v{vid}.ns-RiC{reads}.fasta` — variants routed off the primary track because their `cer_factor` fell below `--min-cer-factor` (default `1.0`)
- `{sample_name}-{gid}.v{vid}.lq-RiC{reads}.fasta` — variants routed off the primary track because their `err_factor` exceeded `--max-err-factor` (default `1.5`); `.lq` takes precedence over `.ns` when both fire
- **FASTQ Files/** (matching FASTQs for `.raw`, `.ns`, `.lq` records)

#### **trees/** (per-specimen variant tree):
- `{specimen}.txt` — ASCII hierarchy of every variant in a specimen (passing, `.ns`, `.lq` together), grouped by core identity group. Each non-anchor variant branches from the larger-size peer with the highest pairwise identity, with a one-line edit summary (substitutions, single-nt indels, short ≤3 nt indels, long indels) versus that parent. Filtered variants carry the on-disk status marker (`-1.v4.ns`, `-1.v8.lq`)

### Naming Convention

The same `-{gid}.v{vid}` schema is used everywhere. There is no separate "original" namespace:

| **Field** | **Meaning** | **Source** |
|-----------|-------------|------------|
| `gid` | Identity group rank within the specimen (1 = largest) | Core (`scipy.cluster.hierarchy.fcluster` complete linkage at `--group-identity`) |
| `vid` | Variant rank within the identity group (1 = anchor / largest) | Core (size-ordered) |
| `.raw{N}` | Pre-merge component of a merged variant | Summarize |
| `.ns` | Variant routed off the primary track by `--min-cer-factor` | Summarize |
| `.lq` | Variant routed off the primary track by `--max-err-factor` | Summarize |

### Example Directory Structure
```
clusters/
├── sample-all.fasta                         # All consensus sequences
├── cluster_debug/
│   ├── sample-1.v1-RiC120-reads.fastq
│   ├── sample-1.v1-RiC120-sampled.fastq
│   └── sample-1.v1-RiC120-untrimmed.fasta
└── metadata/
    └── sample_metadata.json                  # Pipeline diagnostics + cer_details

__Summary__/
├── sample-1.v1-RiC120.fasta                 # Anchor of group 1 (passing)
├── sample-1.v2-RiC45.fasta                  # Second variant in group 1
├── sample-2.v1-RiC30.fasta                  # Anchor of group 2
├── summary.fasta                            # All passing consensus sequences
├── summary.txt                              # Statistics
├── quality_report.txt                       # Quality assessment report
├── summarize_log.txt                        # Processing log
├── FASTQ Files/
│   ├── sample-1.v1-RiC120.fastq
│   ├── sample-1.v2-RiC45.fastq
│   └── sample-2.v1-RiC30.fastq
├── trees/
│   └── sample.txt                           # Per-specimen variant hierarchy
└── variants/
    ├── sample-1.v1.raw1-RiC75.fasta         # Pre-merge component of merged variant
    ├── sample-1.v1.raw2-RiC45.fasta
    ├── sample-1.v3.ns-RiC18.fasta           # Failed --min-cer-factor
    ├── sample-2.v2.lq-RiC22.fasta           # Failed --max-err-factor
    └── FASTQ Files/                         # Matching FASTQs for .raw / .ns / .lq
```

### FASTA Header Metadata

Consensus sequence headers contain metadata fields separated by spaces:

**Core Fields (Always Present):**
- `size=N` — Total number of raw sequence reads in the cluster
- `ric=N` — **Reads in Consensus** — Number of reads actually used for consensus generation (may be less than `size` due to sampling limits)
- `gid=N` — Identity group rank within the specimen (1 = largest)
- `vid=N` — Variant rank within the identity group (1 = anchor / largest)

**Statistical Validation Fields:**
- `cer_factor=X.X` — Per-position multiplicative error rate inflation needed for the cluster to be plausible as basecaller-noise replicates of a larger peer. Higher = more confident the cluster is a real variant. `cer_factor=None` means anchor or no comparison peer found (always passes filtering)
- `err_factor=X.X` — Cluster homogeneity ratio (observed disagreement / expected from q_ctx). Values near 1.0 indicate reads consistent with basecaller noise; values ≫ 1.0 indicate internal heterogeneity beyond what sequencing noise produces

**Optional Fields:**
- `rawric=N+N+...` — RiC values of `.raw` source variants (pre-merge, largest-first; only present in merged variants from `speconsense-summarize`)
- `rawlen=N+N+...` — Original sequence lengths before overlap merging (largest-first; only present in overlap-merged variants)
- `snp=N` — Number of SNP positions from IUPAC merging (only present in merged variants)
- `ambig=N` — Count of IUPAC ambiguity codes in consensus (Y, R, W, S, K, M, etc.)
- `length=N` — Sequence length in bases (via `--fasta-fields` option)
- `primers=list` — Comma-separated list of detected primer names (e.g., `primers=ITS1F,ITS4`)
- `rid=X.X` — Mean read identity percentage (0–100%; via `--fasta-fields qc` preset)
- `rid_min=X.X` — Minimum read identity percentage (worst-case read; via `--fasta-fields qc` preset)

**Example Headers:**
```
# Anchor of group 1, passing both CER and err_factor checks:
>sample-1.v1 size=120 ric=100 gid=1 vid=1 cer_factor=None err_factor=1.02 primers=ITS1F,ITS4

# Second variant in group 1, statistically validated against the anchor:
>sample-1.v2 size=45 ric=45 gid=1 vid=2 cer_factor=12.4 err_factor=0.97 primers=ITS1F,ITS4

# Merged variant from speconsense-summarize:
>sample-1.v1 size=250 ric=250 gid=1 vid=1 rawric=100+89+61 snp=2 ambig=2 cer_factor=None err_factor=1.05 primers=ITS1F,ITS4

# Overlap-merged variant (different-length sequences merged):
>sample-1.v1 size=471 ric=471 gid=1 vid=1 rawric=248+223 rawlen=630+361 primers=ITS1F,ITS4

# QC preset (--fasta-fields qc):
>sample-1.v1 size=250 ric=250 gid=1 vid=1 length=589 rid=98.5 ambig=1
```

**Notes:**
- Use `--fasta-fields` to customize which fields appear in output headers (see [Customizing FASTA Header Fields](#customizing-fasta-header-fields))
- Read identity metrics (`rid`, `rid_min`) reflect homopolymer-normalized sequence identity
- `cer_factor` is filtered by `--min-cer-factor` (default `1.0`) in summarize; `err_factor` by `--max-err-factor` (default `1.5`)
- Variant merging only occurs between sequences with identical primer sets, within the same identity group (or via cross-primer overlap merging across groups)
- `snp` counts positions where IUPAC codes were introduced during merging; `ambig` counts total ambiguity codes in the final sequence
- Per-position CER reproduction data (peer id, K, p*, joint q, context tags, q_ctx values) lives in `clusters/metadata/{specimen}_metadata.json`, not in the FASTA header

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

## Discard Inspection

`--collect-discards` writes all discarded reads (unclustered, MAD outliers, hard-floor or noise-filter dropouts, size-filtered clusters) to `cluster_debug/{sample}-discards.fastq` for inspection. Useful when investigating why a specimen lost reads.

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
- **Note**: Use `speconsense-summarize` after clustering to manage multiple variants per specimen. Key options include `--merge-position-count` for merging variants differing by few SNPs/indels, `--group-identity` for grouping similar variants, and `--select-max-variants` for controlling which variants to output

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

**Absolute size filtering (`--min-size`, default: 3):**
- Filters clusters by absolute number of reads
- Applied at the post-phasing size-filtering phase (after merging, phasing, reassignment, recovery, and CER validation)
- Set to 0 to disable and output all clusters regardless of size

**Processing order (high level):**
1. Initial clustering produces raw clusters
2. Pre-phasing merge of identical/homopolymer-equivalent clusters
3. Variant phasing splits clusters by haplotype
4. Post-phasing merge collapses any newly-equivalent subclusters
5. Noise filter drops clusters with no-majority MSA columns
6. Read reassignment moves reads between clusters within identity groups
7. Discard recovery re-admits previously-dropped reads to surviving clusters
8. Second phasing pass re-phases any clusters that gained reads
9. CER validation annotates each non-anchor cluster with a `cer_factor`
10. **Size filtering: apply `--min-size`**
11. Final consensus + MAD outlier removal + FASTA writing

**Deferred filtering strategy:**
For maximum flexibility in detecting rare variants and contaminants, set `--min-size 0` in `speconsense` and apply final quality thresholds in summarize via `--min-ric`, `--min-cer-factor`, `--max-err-factor`, and `--prune-group-ratio`/`--prune-group-abs`. This lets you run expensive clustering once and experiment with thresholds during post-processing. The CER and err_factor filters route low-confidence clusters to the `.ns`/`.lq` track rather than dropping them entirely, so you can revisit decisions without re-running core.

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

1. **Variant detection**: After initial clustering, Speconsense analyzes positional variation using multiple sequence alignment. A position qualifies as a variant either across the full cluster or within a quality-biased sample of top-quality reads (whichever fires) when the minor allele frequency exceeds the threshold (default 10%) and meets minimum read count requirements.

2. **Position selection**: When multiple variant positions are detected, Speconsense selects the single best position for splitting (minimizing within-cluster error), then recursively regenerates MSA for each subcluster to discover additional variant positions. This hierarchical approach prevents over-fragmentation while allowing deep phasing when supported by the data.

3. **Haplotype separation**: Reads are grouped by their allele combinations at selected variant positions. Each unique combination becomes a separate sub-cluster (haplotype). Reads that don't qualify for any haplotype are discarded rather than force-reassigned, and may be recovered later by phase 7.

4. **IUPAC ambiguity calling**: When only one haplotype qualifies (insufficient support for phasing), variant positions are encoded using IUPAC ambiguity codes (e.g., `Y` for C/T, `R` for A/G). The same machinery handles "no-majority" columns where phasing produced a single haplotype but multiple alleles still exceed thresholds. Ambiguity calling uses a lower count threshold (3 reads vs 5 for splitting), capturing variation that doesn't have enough support to form a viable separate cluster.

**Key parameters for variant phasing:**
- `--min-variant-frequency` — Minimum minor allele frequency to trigger cluster splitting (default: 0.10 = 10%)
- `--min-variant-count` — Minimum read count for minor allele to trigger splitting (default: 3)
- `--disable-position-phasing` — Disable variant phasing entirely (also skips the second phasing pass)
- `--disable-read-reassignment` — Disable cross-cluster reassignment within identity groups (phase 6)
- `--disable-discard-recovery` — Disable re-admission of previously discarded reads to surviving clusters (phase 7; requires `--enable-read-reassignment`)

**Key parameters for IUPAC ambiguity calling:**
- `--min-ambiguity-frequency` — Minimum minor allele frequency for IUPAC codes (default: 0.10 = 10%)
- `--min-ambiguity-count` — Minimum read count for minor allele for IUPAC codes (default: 3)
- `--disable-ambiguity-calling` — Disable IUPAC codes for unphased variants

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

### Pipeline Phases (Core)

After initial clustering, speconsense runs a 12-phase pipeline. Phases 5–9 are the post-phasing refinement work that distinguishes 0.8.x from earlier versions:

| Phase | Name | Purpose |
|------:|------|---------|
| 1 | Initial clustering | MCL graph-based or greedy |
| 2 | Pre-phasing merge | Combine HP-equivalent initial clusters |
| 3 | Variant phasing | Split clusters by haplotype via MSA position analysis |
| 4 | Post-phasing merge | Combine HP-equivalent subclusters |
| 5 | Noise filter | Drop small clusters with no-majority MSA columns |
| 6 | Read reassignment | Move reads between clusters within an identity group when consensus concordance supports it |
| 7 | Discard recovery | Re-admit previously-dropped reads to surviving clusters |
| 8 | Second phasing pass | Re-phase clusters whose membership changed via reassignment/recovery |
| 9 | CER validation | Annotate each non-anchor cluster with a `cer_factor` |
| 10 | Size filtering | Apply `--min-size` |
| 11 | Output generation | Final consensus, MAD outlier removal, FASTA writing |
| 12 | Discard writeout | Optional, with `--collect-discards` |

**Identity groups**: phases 6, 7, and 9 operate within identity groups computed via complete linkage at `--group-identity` (default 0.85). Every pair within a group must meet the threshold, preventing transitive collapse of close-but-distinct variants. Group ranks (`gid`) and within-group variant ranks (`vid`) are emitted into the FASTA header and travel through to summarize.

### Statistical Variant Validation (CER)

Each non-anchor cluster is annotated with a **CER factor** (`cer_factor=`) measuring how implausible it is as basecaller-noise replicates of a larger peer in its identity group. This is a per-position multiplicative inflation of the assumed error rate that would be required for the peer's read counts to plausibly produce the candidate's read counts at the observed disagreement positions, under a binomial test with combinatorial Bonferroni correction over `C(L, K)` site combinations (where `L` is the HP-compressed length and `K` is the number of disagreement sites).

A higher `cer_factor` means more confidence the candidate is a real variant. Anchors and clusters that fail to find a valid pairwise comparison carry `cer_factor=None` and always pass filtering.

**Context-aware error model**: Per-position error rates come from a shipped `dorado-v5.0` model (other models can be selected with `--error-model NAME`, listed via `--list-error-models`). Each variant event is classified as `non-hp-sub`, `non-hp-indel`, or `hp-l{N}` (HP run of length N) using the **reference** consensus's HP context (the artifact hypothesis under test is that the candidate's reads are miscalled copies of the reference). HP runs longer than the table's max length route to blanket HP normalization rather than CER evaluation; `--hp-normalization-length` (default 6) sets the cutoff.

**No filtering happens in core** — the CER decision is a summarize concern, applied via `--min-cer-factor` (default `1.0`, set to `0` to disable). Records below threshold route to `__Summary__/variants/{name}.ns-RiC{ric}.fasta`.

### Cluster Homogeneity (err_factor)

Complementary to CER, **`err_factor`** is a peer-independent metric: for each non-gap consensus column, the fraction of reads disagreeing with the consensus divided by the q_ctx rate predicted for that column's context (HP run length or non-HP). Values near 1.0 mean reads consistent with basecaller noise; values ≫ 1.0 indicate internal heterogeneity beyond what sequencing noise produces.

Unlike `cer_factor`, `err_factor` distinguishes novel-but-real sequences (low `err_factor`, possibly high `cer_factor`) from noise combinations (high `err_factor` regardless of CER).

Filtering happens in summarize via `--max-err-factor` (default `1.5`, set to `0` to disable). Records above threshold route to `__Summary__/variants/{name}.lq-RiC{ric}.fasta`. The `1.5` default is safe because MAD outlier removal at final consensus removes single-read outliers that would otherwise inflate `err_factor` on real clusters.

### Error Models (q_ctx tables)

Both `cer_factor` and `err_factor` rely on a per-context error model that says "what's the basecaller's error rate at a position with this HP context?" These models live as YAML files keyed by event type:

```yaml
name: dorado-v5.0
chemistry: R10.4.1
basecaller: Dorado v5.0
source: Re-estimated on ont98 dataset under speconsense 0.8.0
rates:
  non-hp-sub:  0.005    # Substitution at a non-HP position
  non-hp-indel: 0.010   # Insertion + deletion at a non-HP position
  hp-l1: 0.006          # Length change in an HP run of length 1
  hp-l2: 0.008
  hp-l3: 0.010
  hp-l4: 0.011
  hp-l5: 0.012
```

**Selecting a model:**
```bash
speconsense input.fastq --error-model dorado-v5.0    # Default
speconsense input.fastq --list-error-models          # Show installed models
speconsense input.fastq --error-model /path/to/custom.yaml   # Filesystem path
```

Resolution order: filesystem path → `~/.config/speconsense/error_models/` → bundled. Bundled models include `dorado-v5.0` (default) and `dorado-v3.5`.

**Creating a custom model:**

The easiest way is to **fit one from a finished run** using `speconsense-fit-error-model`:

```bash
speconsense-fit-error-model /path/to/clusters --name my-basecaller
```

This walks the run's `cluster_debug/` and writes `~/.config/speconsense/error_models/my-basecaller.yaml`, then prints a comparison against the model the run was clustered under. See [Fitting a Custom Error Model](#fitting-a-custom-error-model) for the full workflow and the underlying procedure (HP paper §8 re-estimation: approach-1 mode-as-ground-truth at qualifying primary anchors + approach-2 all-cluster non-HP pooling, with biological-variant filtering).

Alternatively, drop a YAML file into `~/.config/speconsense/error_models/my-basecaller.yaml` by hand with the format above, then select it via `--error-model my-basecaller`. HP run lengths absent from the model route to blanket HP normalization rather than context-aware CER evaluation, so models can be sparse.

**HP normalization length** (`--hp-normalization-length`, default 6):
- HP runs of length ≥ this value are blanket-normalized in distance/MSA comparisons (treated as noise)
- HP runs shorter than this are surfaced as candidates and evaluated by context-aware CER (using `hp-l1` ... `hp-l5` from the model)
- Set to 1 for legacy blanket-normalize-all behavior; matches the HP error rate paper's recommendation of CER-evaluating L≤5 HP variants

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

Identity groups are computed once, in core, via complete linkage at `--group-identity` (default `0.85`). `speconsense-summarize` honors those groups via the `gid=`/`vid=` headers and does not re-cluster. Within each group, summarize applies MSA-based variant merging and selection.

Output naming round-trips core's `gid.vid` for every variant unless cross-primer overlap conflation moved that variant into a different group. Filtered-out or merged-away variants leave their vids as gaps under the surviving gid. Moved variants adopt the survivor group's `gid` and receive a freshly-minted `vid` above the highest vid core wrote under that gid (across passed, `.ns`, and `.lq` records); the allocator handles 2+ group conflation without collisions.

**Group identity (cross-primer overlap merging only):**
```bash
speconsense-summarize --group-identity 0.85
```
- Used by summarize **only** for cross-primer overlap merging — bridging two core groups whose anchors overlap above the identity threshold (e.g., ITS and ITS2 amplicons of the same biological entity)
- Default `0.85` matches core's `--group-identity` so cross-primer pairs use the same anchor distance as core's complete-linkage groups
- For raising/lowering core's grouping aggressiveness, set `--group-identity` on `speconsense` itself, then re-run summarize

**Variant Selection (within each group):**

When multiple variants exist per identity group, `speconsense-summarize` selects variants by size (largest first). Use `--select-max-variants N` to cap the number of variants per group (default: no limit).

**Overall Process:**
1. Read FASTA outputs from `speconsense`, parsing `gid=`/`vid=` to recover core's identity groups
2. Cross-primer overlap merger across groups (matches different-primer anchors that overlap above `--group-identity`)
3. For each (possibly conflated) group independently:
   - Apply MSA-based merging to find largest compatible subsets
   - Apply selection size ratio filtering (`--select-min-size-ratio`)
   - Apply variant selection strategy (size or diversity)
   - Output up to `select_max_variants` per group
4. Apply CER filter (`--min-cer-factor`) — route failing records to `.ns`
5. Apply err_factor filter (`--max-err-factor`) — route failing records to `.lq` (`.lq` takes precedence over `.ns` if both fire)
6. Final output contains representatives from all groups; filtered records remain accessible under `__Summary__/variants/`

**Selection Size Ratio Filtering:**
```bash
speconsense-summarize --select-min-size-ratio 0.2
```
- Filters out post-merge variants whose size is too small relative to the group total
- Ratio calculated as `variant_size / group_total` — must be ≥ threshold to keep
- The largest variant in each group is always kept
- Example: `--select-min-size-ratio 0.2` means a variant must represent ≥20% of the group's total reads
- Default is 0 (disabled) — all post-merge variants pass through to selection
- Applied after merging but before variant selection
- Useful for suppressing noise variants that survived merging but are too small to be meaningful
- Set to 0.2 in the `compressed` profile to match the 20% calling threshold theme

This two-stage process ensures that distinct biological sequences are preserved as separate groups, while providing control over variant complexity within each group.

**Group "full" Consensus (query-optimized):**
```bash
speconsense-summarize --enable-full-consensus
```
- Per identity group with ≥2 selected variants on the pass track, emit an additional `-{gid}-full` consensus alongside the phased variants
- Built from a size-weighted, top-mean-Phred read sample (budget = core's `--max-sample-size`, default 100) drawn across the pre-merge core variants whose size clears the running-total gate at `--min-ambiguity-frequency` (also inherited from the core run, default 0.10)
- SPOA with linear gap scoring; local alignment (`-l 0`) for cross-primer-conflated groups, global (`-l 1`) otherwise
- Consensus uses one-vote-per-read, majority-wins gaps, and IUPAC ambiguity codes at columns where ≥2 bases each clear `--min-ambiguity-frequency`
- Suppressed when fewer than 2 pass-track variants exist or fewer than 2 contributors clear the gate; `.ns` / `.lq` records are not eligible to contribute
- Intended use: BLAST query against legacy unphased ITS references. The `-full` output is IUPAC-bearing by construction and should be scored with adjusted-identity (MycoBLAST) tools — under raw BLAST it will silently degrade because every IUPAC code counts as a mismatch
- Enabled by default in the `compressed` profile

**Consensus Column Retention (gap handling):**
```bash
speconsense-summarize --min-position-frequency 0.1 --min-position-count 3
```
- Controls how gap vs. base disagreements are resolved in merged and `-full` consensus sequences
- A column is retained when the fraction of non-gap content is ≥ `--min-position-frequency` AND the absolute non-gap support is ≥ `--min-position-count`
- Default `0.5` matches majority-wins behavior: positions where more than half the (size-weighted) votes are gaps are omitted
- Lower values (e.g. `0.1`) preserve positions where a minority of contributors carry content — useful when merging variants of different lengths or when indel events should not delete content from the merged consensus
- `--min-position-count` (default 3) prevents columns with negligible support from surviving at low frequency thresholds
- Set to `0.1` in the `compressed` profile alongside its aggressive merge settings
- Applies to within-group MSA merging, overlap merging (overlap region only — non-overlap regions always preserve content from contributing sequences), and `-full` group consensus

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
- `size` — Total reads across merged variants
- `ric` — Reads in consensus
- `length` — Sequence length in bases
- `rawric` — RiC values of `.raw` source variants (only when merged)
- `rawlen` — Original sequence lengths before overlap merging (only when overlap-merged)
- `snp` — Number of IUPAC positions from merging (only when >0)
- `ambig` — Count of IUPAC ambiguity codes in consensus
- `rid` — Mean read identity percentage (when available)
- `rid_min` — Minimum read identity percentage (when available)
- `cer_factor` — Statistical variant validation factor (None for anchors)
- `err_factor` — Cluster homogeneity ratio (observed/expected disagreement)
- `primers` — Detected primer names (when detected)
- `group` — Identity group number (extracted from `-{gid}.v{vid}` filename; superseded by `gid=` from core)
- `variant` — Variant identifier within group (extracted from filename; superseded by `vid=` from core)

**Preset definitions:**
- `default`: `size, ric, rawric, rawlen, snp, ambig, primers`
- `minimal`: `size, ric`
- `qc`: `size, ric, length, rid, ambig, cer_factor, err_factor`
- `full`: `size, ric, length, rawric, rawlen, snp, ambig, rid, cer_factor, err_factor, primers`
- `id-only`: (no fields)

The `default` preset does not include `cer_factor` / `err_factor` / `gid` / `vid`. If you want CER and homogeneity metrics in summarize-emitted FASTAs, use `--fasta-fields qc` or `--fasta-fields full`. Core's own FASTA always includes `gid`, `vid`, `cer_factor`, and `err_factor` regardless of the summarize preset.

**Use cases:**
- **Downstream tool compatibility**: Use `minimal` or `id-only` for tools expecting simple headers
- **Quality control**: Use `qc` preset to include read identity metrics for assessing consensus quality
- **File size optimization**: Use `minimal` to reduce file size for large datasets
- **Custom workflows**: Combine presets and fields for workflow-specific needs

### Quality Assessment and Reporting

Speconsense-summarize automatically generates a `quality_report.txt` to help prioritize manual review of sequences with potential quality concerns. The report is inspection-oriented — it highlights specimens and variants that warrant a closer look, rather than dumping descriptive statistics.

**Report Generation:**
- Created automatically in the summary output directory
- No configuration required — generated regardless of `--fasta-fields` settings
- Focuses on actionable issues, not exhaustive enumeration

**Report Sections:**

**1. Executive Summary:**
- Pass/`.ns`/`.lq` counts across the run
- Specimens that had input reads but produced no clusters (often a demultiplexing or input-quality signal)
- Filter decisions from `--min-cer-factor` and `--max-err-factor`

**2. Passed Variants Worth Inspecting:**
- Run-relative outliers on `err_factor` and ambiguity count (mean ± 2σ)
- Variants that passed filters but stand out from the rest of the run

**3. Possible Rescues:**
- Low-yield specimens whose filtered (`.ns`/`.lq`) variants are close to the threshold
- Suggests manual review of variants that might warrant a relaxed threshold for that specimen

**4. Overlap Merge Details** (when present):
- Lists cross-primer overlap merges, with overlap bp and prefix/suffix extensions
- Flags edge cases (overlap near threshold, large length ratios)

**5. Run-Wide Parameter Signals:**
- Suggests parameter tuning when a metric is systematically high or low across the run

**6. Pipeline Activity:**
- Phase-by-phase activity counts (clusters disbanded, reads reassigned, reads recovered, second-pass splits, CER-filtered, err_factor-filtered) aggregated across specimens
- Helps surface unexpected pipeline behavior on a new dataset

**Recommended Actions:**

**For high `err_factor` outliers:**
- The cluster's reads disagree with the consensus more than the basecaller's noise model predicts
- Review `cluster_debug/` FASTQ files in source directory
- Check for biological heterogeneity (multiple true variants in one cluster) — if so, consider lowering `--min-variant-frequency` or `--min-variant-count` and re-running core
- May indicate cross-specimen contamination if the disagreement pattern aligns with another specimen

**For high ambiguity count outliers:**
- Review whether ambiguity represents true heterozygosity, merging artifact, or unphased variation
- Check if variant phasing would produce cleaner separate haplotypes (lower `--min-variant-frequency`)
- Cross-reference with `__Summary__/trees/{specimen}.txt` to see the variant hierarchy

**For `.ns` records (failed CER):**
- The cluster is statistically explained as basecaller noise replicated from a larger peer
- Most are correctly filtered. To rescue: lower `--min-cer-factor` (default `1.0`; `0.5` is a common relaxation) and re-run summarize

**For `.lq` records (failed err_factor):**
- The cluster's internal heterogeneity exceeds the basecaller noise model
- Consider whether the cluster is real but degraded (herbarium samples often inflate `err_factor`) — the `herbarium` profile disables this filter

### Additional Summarize Options

**Quality Filtering:**
```bash
speconsense-summarize --min-ric 5
```
- Filters out consensus sequences with fewer than the specified number of Reads in Consensus (RiC)
- Default is 3 — only sequences supported by at least 3 reads are processed
- Higher values provide more stringent quality control but may exclude valid low-abundance variants

**CER Filtering (statistical variant validation):**
```bash
speconsense-summarize --min-cer-factor 1.0      # Default
speconsense-summarize --min-cer-factor 0.5      # More permissive
speconsense-summarize --min-cer-factor 0        # Disable
```
- `--min-cer-factor`: Minimum per-position CER factor for a variant to be kept on the primary track (default: `1.0`)
- Variants with `cer_factor` below threshold route to `__Summary__/variants/{name}.ns-RiC{ric}.fasta` (and matching FASTQ) — they are not deleted, just sidelined
- `cer_factor=None` (anchors and clusters without a valid pairwise comparison) always passes
- See [Statistical Variant Validation (CER)](#statistical-variant-validation-cer) for the underlying model

**Cluster Homogeneity Filtering (err_factor):**
```bash
speconsense-summarize --max-err-factor 1.5      # Default
speconsense-summarize --max-err-factor 2.0      # More permissive (e.g., for herbarium specimens)
speconsense-summarize --max-err-factor 0        # Disable
```
- `--max-err-factor`: Maximum cluster `err_factor` (observed/q_ctx-expected disagreement) for primary-track inclusion (default: `1.5`)
- Variants above threshold route to `__Summary__/variants/{name}.lq-RiC{ric}.fasta`
- `.lq` takes precedence over `.ns` when both fire — a cluster filtered by both labels appears only as `.lq`
- See [Cluster Homogeneity (err_factor)](#cluster-homogeneity-err_factor) for the underlying metric

**Length Filtering:**
```bash
speconsense-summarize --min-len 400 --max-len 800
```
- `--min-len`: Minimum sequence length in bp (default: 0 = disabled)
- `--max-len`: Maximum sequence length in bp (default: 0 = disabled)
- Applied during initial sequence loading
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
- `--merge-min-size-ratio R`: Minimum size ratio (contributor/merged total) for merging clusters (default: 0.1, 0 to disable)
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
- For each candidate merge subset, every contributor must be ≥ threshold fraction of the subset total
- Example: `--merge-min-size-ratio 0.1` means every contributor must represent ≥10% of the merged total
- Default is 0.1
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

1. **Load sequences** with RiC filtering (`--min-ric`) and length filtering (`--min-len`, `--max-len`); parse `gid=`/`vid=` from FASTA headers (hard-fails if missing)
2. **Bucket by core-assigned identity group** — no re-grouping
3. **Cross-primer overlap conflation** (`--min-merge-overlap`, `--group-identity`) merges different-primer core groups whose anchors overlap above the identity threshold
4. **Homopolymer-aware MSA-based variant merging** within each (possibly-conflated) group (`--disable-merging`, `--merge-effort`, `--merge-position-count`, `--merge-indel-length`, `--min-merge-overlap`, `--merge-snp`, `--merge-min-size-ratio`, `--disable-homopolymer-equivalence`, `--hp-normalization-length`)
5. **Selection size ratio filtering** to remove tiny post-merge variants (`--select-min-size-ratio`)
6. **Variant selection** within each group (`--select-max-variants`, `--select-max-groups`)
7. **CER and err_factor filtering** — route failing records to `.ns` / `.lq` (`--min-cer-factor`, `--max-err-factor`)
8. **Output generation** — passing FASTA + FASTQ, filtered records under `variants/`, per-specimen tree under `trees/`, customizable header fields (`--fasta-fields`) and full traceability

**Key architectural features**:
- Identity grouping happens once, in core, via complete-linkage HAC at `--group-identity` (default 0.85). Summarize honors that grouping via `gid`/`vid` headers — no independent HAC pass
- Cross-primer overlap conflation is the only operation that crosses core-group boundaries, gated by anchor identity AND sufficient sequence overlap
- Merging is applied independently within each (possibly-conflated) group using MSA-based consensus generation
- Homopolymer-aware merging by default (AAA ≈ AAAA, runs ≥ `--hp-normalization-length` blanket-normalized) to match pipeline-wide adjusted-identity semantics
- CER and err_factor filtering routes records rather than deleting them — `.ns` / `.lq` files preserve the data for review and threshold tuning

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
                   [--algorithm {graph,greedy}] [--min-identity MIN_IDENTITY]
                   [--inflation INFLATION]
                   [--k-nearest-neighbors K_NEAREST_NEIGHBORS]
                   [--min-size MIN_SIZE]
                   [--max-sample-size MAX_SAMPLE_SIZE]
                   [--disable-position-phasing] [--enable-position-phasing]
                   [--disable-read-reassignment] [--enable-read-reassignment]
                   [--disable-discard-recovery] [--enable-discard-recovery]
                   [--min-variant-frequency MIN_VARIANT_FREQUENCY]
                   [--min-variant-count MIN_VARIANT_COUNT]
                   [--significance-level SIGNIFICANCE_LEVEL]
                   [--group-identity GROUP_IDENTITY]
                   [--hp-normalization-length HP_NORMALIZATION_LENGTH]
                   [--error-model ERROR_MODEL]
                   [--disable-ambiguity-calling] [--enable-ambiguity-calling]
                   [--min-ambiguity-frequency MIN_AMBIGUITY_FREQUENCY]
                   [--min-ambiguity-count MIN_AMBIGUITY_COUNT]
                   [--disable-cluster-merging] [--enable-cluster-merging]
                   [--disable-homopolymer-equivalence]
                   [--enable-homopolymer-equivalence]
                   [--orient-mode {skip,keep-all,filter-failed}]
                   [--presample PRESAMPLE] [--scale-threshold SCALE_THRESHOLD]
                   [--threads N] [--collect-discards]
                   [--no-collect-discards]
                   [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                   [--version] [-p NAME] [--list-profiles]
                   [--list-error-models]
                   input_file

MCL-based clustering of nanopore amplicon reads

options:
  -h, --help            show this help message and exit
  --version             Show program's version number and exit
  -p NAME, --profile NAME
                        Load parameter profile (use --list-profiles to see
                        available)
  --list-profiles       List available profiles and exit
  --list-error-models   List available error models and exit

Input/Output:
  input_file            Input FASTQ file
  -O OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for all files (default: clusters)
  --primers PRIMERS     FASTA file containing primer sequences (default: looks
                        for primers.fasta in input file directory)

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
  --min-size MIN_SIZE   Minimum cluster size (default: 3, 0 to disable)
  --max-sample-size MAX_SAMPLE_SIZE
                        Maximum cluster size for consensus (default: 100)

Variant Phasing:
  --disable-position-phasing
                        Disable position-based variant phasing (enabled by
                        default). MCL graph clustering already separates most
                        variants; this second pass analyzes MSA positions to
                        phase remaining variants.
  --enable-position-phasing
                        Override --disable-position-phasing or profile setting
  --disable-read-reassignment
                        Disable post-phasing read reassignment within identity
                        groups. Reassignment moves reads between clusters
                        based on consensus concordance; disable to preserve
                        original phasing-time membership.
  --enable-read-reassignment
                        Override --disable-read-reassignment or profile
                        setting
  --disable-discard-recovery
                        Disable recovery of previously discarded reads into
                        existing clusters. Has no effect if
                        --disable-read-reassignment is also set (recovery
                        requires reassignment).
  --enable-discard-recovery
                        Override --disable-discard-recovery or profile
                        setting
  --min-variant-frequency MIN_VARIANT_FREQUENCY
                        Minimum alternative allele frequency to call variant
                        (default: 0.10 for 10%)
  --min-variant-count MIN_VARIANT_COUNT
                        Minimum alternative allele read count to call variant
                        (default: 3)
  --significance-level SIGNIFICANCE_LEVEL
                        Significance level (alpha) for variant significance
                        testing (default: 1e-5)
  --group-identity GROUP_IDENTITY
                        Minimum pairwise identity to group clusters for read
                        reassignment, discard recovery, and CER validation.
                        Grouping uses complete linkage: every pair within a
                        group must meet this threshold. (default: 0.85)
  --hp-normalization-length HP_NORMALIZATION_LENGTH
                        Minimum homopolymer run length at/above which HP
                        length variants are blanket-normalized (treated as
                        noise). Runs of length >= this value have their
                        length variants suppressed in MSA variant detection;
                        runs shorter than this surface as candidates and are
                        evaluated by context-aware CER. Default 6 matches
                        the HP paper recommendation of CER-evaluating
                        L <= 5 HP variants. (default: 6)
  --error-model ERROR_MODEL
                        Per-basecaller error model used for context-aware
                        variant validation. Either a shipped model name (use
                        --list-error-models to see available), a user model in
                        ~/.config/speconsense/error_models/, or a filesystem
                        path to a YAML file with a 'rates' mapping. (default:
                        dorado-v5.0)

Ambiguity Calling:
  --disable-ambiguity-calling
                        Disable IUPAC ambiguity code calling for unphased
                        variant positions
  --enable-ambiguity-calling
                        Override --disable-ambiguity-calling or profile
                        setting
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
  --enable-cluster-merging
                        Override --disable-cluster-merging or profile setting
  --disable-homopolymer-equivalence
                        Disable homopolymer equivalence in cluster merging
                        (only merge identical sequences)
  --enable-homopolymer-equivalence
                        Override --disable-homopolymer-equivalence or profile
                        setting

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

Debugging:
  --collect-discards    Write discarded reads (outliers and filtered clusters)
                        to cluster_debug/{sample}-discards.fastq
  --no-collect-discards
                        Override --collect-discards or profile setting
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
```

### speconsense-summarize Options

```
usage: speconsense-summarize [-h] [--source SOURCE]
                             [--summary-dir SUMMARY_DIR] [--specimen SPECIMEN]
                             [--aggregate-only] [--fasta-fields FASTA_FIELDS]
                             [--min-ric MIN_RIC] [--min-len MIN_LEN]
                             [--max-len MAX_LEN]
                             [--min-cer-factor MIN_CER_FACTOR]
                             [--max-err-factor MAX_ERR_FACTOR]
                             [--group-identity GROUP_IDENTITY]
                             [--disable-merging] [--enable-merging]
                             [--merge-snp | --no-merge-snp]
                             [--merge-indel-length MERGE_INDEL_LENGTH]
                             [--merge-position-count MERGE_POSITION_COUNT]
                             [--merge-min-size-ratio MERGE_MIN_SIZE_RATIO]
                             [--min-merge-overlap MIN_MERGE_OVERLAP]
                             [--disable-homopolymer-equivalence]
                             [--enable-homopolymer-equivalence]
                             [--merge-effort LEVEL]
                             [--hp-normalization-length HP_NORMALIZATION_LENGTH]
                             [--select-max-groups SELECT_MAX_GROUPS]
                             [--select-max-variants SELECT_MAX_VARIANTS]
                             [--select-min-size-ratio SELECT_MIN_SIZE_RATIO]
                             [--min-position-frequency MIN_POSITION_FREQUENCY]
                             [--min-position-count MIN_POSITION_COUNT]
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
  --specimen SPECIMEN   Process only this specimen. Loads only
                        <specimen>-all.fasta from --source.
  --aggregate-only      Skip processing. Generate aggregate summary from
                        existing per-specimen outputs.
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
  --min-cer-factor MIN_CER_FACTOR
                        Minimum per-position CER factor for a variant to be
                        kept as a primary output. Variants with cer_factor
                        below this are routed to __Summary__/variants/ as .ns
                        records. Variants with cer_factor=None (anchors,
                        clusters without a valid pairwise comparison) always
                        pass. Set to 0 to disable CER filtering. (default:
                        1.0)
  --max-err-factor MAX_ERR_FACTOR
                        Maximum cluster err_factor (observed/q_ctx-expected
                        disagreement ratio). Clusters above this threshold are
                        routed to __Summary__/variants/ as .lq records.
                        Variants with err_factor=None (legacy output) always
                        pass. Set to 0 to disable err_factor filtering.
                        (default: 1.5)

Grouping:
  --group-identity GROUP_IDENTITY, --variant-group-identity GROUP_IDENTITY
                        Anchor-to-anchor identity threshold for cross-primer
                        overlap merging between core-assigned groups. Matches
                        core's --group-identity default. (default: 0.85)

Merging:
  --disable-merging     Disable all variant merging (skip MSA-based merge
                        evaluation entirely)
  --enable-merging      Override --disable-merging or profile setting
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
                        Minimum size ratio (contributor/merged total) for
                        merging clusters. Subsets where any contributor is
                        below this fraction are skipped. (default: 0.1,
                        0 to disable)
  --min-merge-overlap MIN_MERGE_OVERLAP
                        Minimum overlap in bp for merging sequences of
                        different lengths (default: 200, 0 to disable)
  --disable-homopolymer-equivalence
                        Disable homopolymer equivalence in merging (treat AAA
                        vs AAAA as different)
  --enable-homopolymer-equivalence
                        Override --disable-homopolymer-equivalence or profile
                        setting
  --merge-effort LEVEL  Merging effort level: fast (8), balanced (10),
                        thorough (12), or numeric 6-14. Higher values allow
                        larger batch sizes for exhaustive subset search.
                        Default: balanced
  --hp-normalization-length HP_NORMALIZATION_LENGTH
                        Minimum homopolymer run length at/above which HP
                        length differences are blanket-normalized (treated as
                        noise). Runs shorter than this are surfaced as real
                        edits in both distance calculations and MSA merging.
                        Matches core's --hp-normalization-length default.
                        (default: 6)

Selection:
  --select-max-groups SELECT_MAX_GROUPS, --max-groups SELECT_MAX_GROUPS
                        Maximum number of groups to output per specimen
                        (default: -1 = all groups)
  --select-max-variants SELECT_MAX_VARIANTS, --max-variants SELECT_MAX_VARIANTS
                        Maximum total variants to output per group (default:
                        -1 = no limit, 0 also means no limit)
  --select-min-size-ratio SELECT_MIN_SIZE_RATIO
                        Minimum size ratio (variant/group total) to include in
                        output. The largest variant in each group is always
                        kept. (default: 0 = disabled, e.g. 0.2 for 20% cutoff)

Consensus output:
  --min-position-frequency MIN_POSITION_FREQUENCY
                        Minimum fraction of non-gap content at a column to
                        retain the position in merged and -full consensus
                        sequences. Default 0.5 matches majority-wins behavior.
                        Set lower (e.g. 0.1) to preserve positions where a
                        minority of contributors carry content. (default: 0.5)
  --min-position-count MIN_POSITION_COUNT
                        Minimum absolute non-gap support (size-weighted votes
                        for merging, read count for -full) to retain a column.
                        Both --min-position-frequency and --min-position-count
                        must be met. (default: 3)

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

### Fitting a Custom Error Model

**Use Case:** Calibrating CER and `err_factor` against your own basecaller version, chemistry, library prep, or amplicon system — without waiting for a shipped update. Typical scenarios: a new Dorado release the suite doesn't yet ship a model for, a non-fungal amplicon system whose error profile differs from the bundled ITS-tuned models, or confirming that a custom flowcell/library prep behaves like the shipped table predicts.

`speconsense-fit-error-model` productizes the offline q_ctx re-estimation procedure documented in the HP paper §8 / CER paper §4.2 ("Phase 1 deployment regime"). It walks a finished speconsense output tree, applies the paper's selection and filtering, and writes a YAML error model to `~/.config/speconsense/error_models/`.

**Basic usage:**

```bash
# Default thresholds (RiC>=200, err_factor<1.0) — matches the v5.0 retune
speconsense-fit-error-model /path/to/clusters --name my-model
```

The output model becomes pickable in subsequent runs via `--error-model my-model`. The tool also prints a comparison against the model the run was clustered under (or `--compare-against` if explicit), with loud-warns when any key drifts more than 2× from the source — useful for catching cross-basecaller comparisons or material library-prep effects.

**Common options:**

```bash
# Lower thresholds for shallower datasets (paper's ont37 setting)
speconsense-fit-error-model /path/to/clusters --name my-model --min-ric 50

# Tighten the approach-2 non-HP RiC floor (matches paper §8 ont37 methods)
speconsense-fit-error-model /path/to/clusters --name my-model --nonhp-min-ric 20

# Override the auto-detected comparison model
speconsense-fit-error-model /path/to/clusters --name my-model --compare-against dorado-v3.5

# Annotate the YAML frontmatter (chemistry/basecaller default to source-model values)
speconsense-fit-error-model /path/to/clusters --name my-model \
  --chemistry R10.4.1 --basecaller "Dorado SUP v5.1" --dataset run123 \
  --description "Custom calibration for run123 LSU amplicons"

# Preview without writing
speconsense-fit-error-model /path/to/clusters --name my-model --dry-run

# Write diagnostic TSVs alongside the model
speconsense-fit-error-model /path/to/clusters --name my-model \
  --debug-dir /path/to/debug-output
```

**Selection criteria** (paper defaults):
- `--min-ric` (default 200): primary anchor must have at least this many reads. Use 50 for lower-depth datasets — the paper's ont37 retune used 50.
- `--max-err-factor` (default 1.0): primary anchor's `err_factor` must be below this. Acts as a single-template plausibility gate; set 0 to disable.
- `--nonhp-min-ric` (default 5): cluster must have at least this many reads to contribute to approach-2 non-HP rate pooling.

**How it works** (HP paper §8):
1. **Indexes specimens** via `cluster_debug/*-metadata.json` (requires schema_version 2.0+ — re-cluster on a current speconsense version if you see schema warnings).
2. **Approach 1 (HP rates)**: For each qualifying specimen's primary-anchor MSA, identifies HP runs in the consensus, extracts each read's called HP length via flanking-anchor columns, takes the mode across reads as ground truth, and aggregates per `(base, hp_length)` context. Bonferroni-corrected binomial outlier and bimodal-distribution filters exclude positions that look like biological HP-only variants. The final `hp-l{N}` value is the **implied per-position error rate** `p` such that `(1-p)^N` equals the observed run-level fraction-correct (HP paper §3.5), pooled across bases.
3. **Approach 2 (non-HP rates)**: Walks every cluster MSA in the tree (regardless of qualification), computes per-position errors against each cluster's own consensus excluding HP runs of length > 1, and pools across clusters. The pooled rates match the operational distribution CER evaluates against in production.
4. **Source-model comparison**: Auto-detects the model the run was clustered under (the most common `parameters.error_model` across metadata JSONs), loads it via `qctx.load_table`, and prints a per-key diff.
5. **YAML output**: Writes the model to `~/.config/speconsense/error_models/{name}.yaml`, creating the directory if needed. Refuses to overwrite without `--force`.

**Reproducibility check:**

Run against the same dataset that produced a shipped model and compare. On `/Users/josh/mm/data/ont37/hp-analysis/clusters-2026-05-26-pass2` with paper-matching settings (`--min-ric 50 --nonhp-min-ric 20`), the tool reproduces `dorado-v3.5.yaml` essentially exactly (all eight keys within 3% of the shipped values).

**See also:** the HP error rate paper (`docs/hp_error_rate/hp_error_rate_report.pdf`) for the empirical foundation and the CER paper (`docs/cer_in_practice/cer_in_practice.pdf`) §4.2 for the framework's "Phase 1" deployment regime that this tool implements.

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

### Testing with Synthetic Data

For empirical testing of consensus quality, variant detection, contamination scenarios, and understanding toolchain behavior with controlled datasets, see [Testing with Speconsense-Synth](docs/synthetic-testing.md).

## Future Enhancements

The following features are under consideration for future development:

- **Background contamination detection tool**: A utility to identify contamination patterns affecting multiple specimens within the same sequencing run, helping to distinguish systematic contamination from genuine biological sequences

- **Further scalability improvements**: 0.8.x added vsearch-based sparse identity grouping, sparse discard screening, and within-group top-K CER, which collectively make 100K-read inputs tractable. Continued work targets graph clustering for inputs in the hundreds of thousands of reads

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
