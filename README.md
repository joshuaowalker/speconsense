# Speconsense

A tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads.

## Overview

Speconsense is a specialized tool for generating high-quality consensus sequences from Oxford Nanopore Technology (ONT) amplicon data. It is specifically designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

The key features of Speconsense include:
- Robust clustering of amplicon reads using either Markov Clustering (MCL) graph-based or greedy algorithms
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

# Set minimum identity threshold for clustering
speconsense input.fastq --min-identity 0.85

# Control the maximum sample size for consensus generation
speconsense input.fastq --max-sample-size 500

# Specify output directory (default: clusters)
speconsense input.fastq --output-dir results/

# Using short form for output directory
speconsense input.fastq -O my_results/
```

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
- `summary.fasta` - All final consensus sequences
- Individual FASTA files per specimen (including `.raw` files for merged variants)
- `summary.txt` - Statistics and metrics
- `quality_report.txt` - Prioritized list of sequences with potential quality concerns
- `FASTQ Files/` - Reads contributing to each consensus

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
- **Individual FASTA files**: `{sample_name}-{group}.v{variant}-RiC{reads}.fasta` (all variants including primary)
- **Combined file**: `summary.fasta` - all final consensus sequences
- **Statistics**: `summary.txt` - sequence counts and metrics
- **Quality report**: `quality_report.txt` - highlights sequences with potential quality concerns
- **Log file**: `summarize_log.txt` - complete processing log

#### **FASTQ Files/** (aggregated reads for final consensus):
- `{sample_name}-{group}.v{variant}-RiC{reads}.fastq` - all reads contributing to variant consensus
- `{sample_name}-{group}.v{variant}.raw{N}-RiC{reads}.fastq` - reads for pre-merge variant N (when variants were merged)

### Naming Convention Summary

**Two distinct namespaces maintain traceability:**

| **Namespace** | **Used In** | **Format** | **Purpose** |
|---------------|-------------|------------|-------------|
| **Original** | Source `cluster_debug/` | `-c1`, `-c2`, `-c3` | Preserves speconsense clustering results |
| **Summarization** | `__Summary__/`, `FASTQ Files/` | `-1.v1`, `-1.v2`, `-2.v1`, `.raw1` | Post-processing groups and variants |

### Example Directory Structure
```
__Summary__/
├── sample-1.v1-RiC45.fasta                  # Primary variant (group 1, merged)
├── sample-1.v1.raw1-RiC30.fasta             # Pre-merge variant 1
├── sample-1.v1.raw2-RiC15.fasta             # Pre-merge variant 2
├── sample-1.v2-RiC23.fasta                  # Additional variant (not merged)
├── sample-2.v1-RiC30.fasta                  # Second organism group, primary variant
├── summary.fasta                            # All sequences combined
├── summary.txt                              # Statistics
├── quality_report.txt                       # Quality assessment report
├── summarize_log.txt                        # Processing log
└── FASTQ Files/                             # Reads for final consensus
    ├── sample-1.v1-RiC45.fastq              # All reads for merged consensus
    ├── sample-1.v1.raw1-RiC30.fastq         # Reads for pre-merge variant 1
    ├── sample-1.v1.raw2-RiC15.fastq         # Reads for pre-merge variant 2
    ├── sample-1.v2-RiC23.fastq              # Reads for additional variant
    └── sample-2.v1-RiC30.fastq              # All reads for second group
```

### FASTA Header Metadata

Consensus sequence headers contain metadata fields separated by spaces:

**Core Fields (Always Present):**
- `size=N` - Total number of raw sequence reads in the cluster
- `ric=N` - **Reads in Consensus** - Number of reads actually used for consensus generation (may be less than size due to sampling limits)

**Optional Fields:**
- `rawric=N+N+...` - RiC values of .raw source variants (pre-merge, largest-first, only present in merged variants from speconsense-summarize)
- `snp=N` - Number of SNP positions with IUPAC ambiguity codes (only present in merged variants from speconsense-summarize)
- `length=N` - Sequence length in bases (available via --fasta-fields option)
- `primers=list` - Comma-separated list of detected primer names (e.g., `primers=ITS1F,ITS4`)
- `p50diff=X.X` - Median (p50) edit distance from stability assessment (available via --fasta-fields qc preset)
- `p95diff=X.X` - 95th percentile (p95) edit distance from stability assessment (available via --fasta-fields qc preset)
- `group=N` - Variant group number (available via --fasta-fields option)
- `variant=vN` - Variant identifier within group (available via --fasta-fields option, only for variants)

**Example Headers:**
```
# Simple consensus from speconsense (debug files):
>sample-c1 size=50 ric=45 p50diff=2.1 p95diff=4.8 primers=ITS1F,ITS4

# Merged variant from speconsense-summarize (default fields):
>sample-1 size=250 ric=250 rawric=100+89+61 snp=2 primers=ITS1F,ITS4

# With QC preset (--fasta-fields qc):
>sample-1 size=250 ric=250 length=589 p50diff=0.0 p95diff=2.0
```

**Notes:**
- Use `--fasta-fields` option to customize which fields appear in output headers (see [Customizing FASTA Header Fields](#customizing-fasta-header-fields))
- Stability metrics (`p50diff`, `p95diff`) are available in speconsense debug files and can be included in summarize output via `--fasta-fields qc`
- Field name changes from earlier versions: `merged_ric` → `rawric`, `median_diff` → `p50diff`, `p95_diff` → `p95diff`
- Variant merging only occurs between sequences with identical primer sets
- SNP counts reflect IUPAC ambiguity positions in consensus sequences

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

**Relative size filtering (`--min-cluster-ratio`, default: 0.2):**
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

### Stability Assessment

To evaluate the reliability of consensus sequences, Speconsense performs stability assessment by:
1. Generating multiple consensus sequences from random subsets of the cluster
2. Measuring the adjusted identity between each subsample consensus and the full consensus
3. Converting identity scores to distance metrics for stability reporting
4. Reporting the median and 95th percentile distances as stability metrics

Sequences with higher than expected distance should be checked for bioinformatic contamination and/or underlying biological variation, as these may indicate mixed clusters or other issues with the data.

**Adjusted Identity Scoring**: Speconsense uses the adjusted-identity algorithm with homopolymer normalization for more accurate sequence comparisons. This means that differences in homopolymer run lengths (e.g., AAA vs AAAAA) are treated as identical, which is particularly important for nanopore sequencing data where homopolymer length calling can be inconsistent.

### Primer Trimming

When a primers file is provided via `--primers`, Speconsense will identify and trim primer sequences from the 5' and 3' ends of consensus sequences, producing clean amplicon sequences for downstream analysis.

**Automatic primer detection**: If `--primers` is not specified, Speconsense will automatically look for `primers.fasta` in the same directory as the input FASTQ file. If found, primer trimming will be enabled automatically.

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
# Output: >sample-1 size=638 ric=638 rawric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC

# Minimal headers - just the essentials
speconsense-summarize --fasta-fields minimal
# Output: >sample-1 size=638 ric=638

# QC preset - includes stability metrics and length
speconsense-summarize --fasta-fields qc
# Output: >sample-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=2.0

# Full metadata (all available fields)
speconsense-summarize --fasta-fields full
# Output: >sample-1 size=638 ric=638 length=589 rawric=333+305 snp=1 p50diff=0.0 p95diff=2.0 primers=...

# ID only (no metadata)
speconsense-summarize --fasta-fields id-only
# Output: >sample-1
```

**Custom field selection:**
```bash
# Specify individual fields (comma-separated)
speconsense-summarize --fasta-fields size,ric,primers

# Combine presets and fields
speconsense-summarize --fasta-fields minimal,p50diff,p95diff

# Combine multiple presets
speconsense-summarize --fasta-fields minimal,qc
```

**Available fields:**
- `size` - Total reads across merged variants
- `ric` - Reads in consensus
- `length` - Sequence length in bases
- `rawric` - RiC values of .raw source variants (only when merged)
- `snp` - Number of IUPAC ambiguity positions (only when >0)
- `p50diff` - Median stability difference (when available)
- `p95diff` - 95th percentile stability difference (when available)
- `primers` - Detected primer names (when detected)
- `group` - Variant group number
- `variant` - Variant identifier within group (only for variants)

**Use cases:**
- **Downstream tool compatibility**: Use `minimal` or `id-only` for tools expecting simple headers
- **Quality control**: Use `qc` preset to include stability metrics for assessing consensus quality
- **File size optimization**: Use `minimal` to reduce file size for large datasets
- **Custom workflows**: Combine presets and fields for workflow-specific needs

### Quality Assessment and Reporting

Speconsense-summarize automatically generates a `quality_report.txt` file to help prioritize manual review of sequences with potential quality concerns. This is particularly valuable for high-throughput workflows where human time for manual inspection is limited.

**Report Generation:**
- Created automatically in the summary output directory
- No configuration required - generated regardless of `--fasta-fields` settings
- Focuses exclusively on sequences that may need review (no "excellent" sequences listed)

**What the Report Includes:**

**1. Elevated Variation (High Priority Section):**
- Sequences where p50diff > 0 or p95diff > 0
- Includes both unmerged sequences and merged sequences with problematic components
- Sorted by (p50diff, p95diff) descending for easy prioritization
- p50diff > 0 indicates **systematic heterogeneity** - more than half of stability trials produced variant consensuses
- p95diff > 0 (with p50diff = 0) indicates **outlier variation** - only worst 5% of trials showed variation

**2. Small Components (Lower Priority Section):**
- Merged sequences where at least one component had RiC < 21
- Small components lack stability metrics (too few reads for subsampling)
- May represent low-abundance variants or contamination

**Understanding Stability Metrics:**

Stability assessment works by:
1. Generating 100 consensus sequences from random subsets of the cluster (20 reads each)
2. Measuring edit distance between each subsample consensus and the full consensus
3. Reporting median (p50diff) and 95th percentile (p95diff) distances

**Key points about distance calculation:**
- Uses homopolymer normalization - differences like AAA vs AAAAA don't count as variation
- Therefore, elevated metrics indicate **substitutions or indels**, not homopolymer issues
- p50diff > 0 is much more concerning than p95diff > 0
- Both metrics = 0 indicates excellent consensus stability

**For Merged Sequences:**
- Merged sequences themselves don't have stability metrics (can't subsample IUPAC consensus)
- Report shows stability of the **worst component** (.raw file) that was merged
- Only included if: (A) a component has elevated variation, OR (B) a component is too small (RiC < 21)

**Example Report Entry:**
```
Type    Sequence                                              RiC    p50diff  p95diff  Notes
------------------------------------------------------------------------------------------------
        specimen-1                                            450        2.1      4.8  Consensus instability
MERGED  specimen-2                                            320        1.5      3.2  raw1: p50=1.5, p95=3.2 (RiC=180)
        specimen-3                                            215        0.0      2.0  Outlier variation
```

**Recommended Actions:**

**For p50diff > 0 (HIGH PRIORITY):**
- Review cluster using `cluster_debug/` FASTQ files in source directory
- Check for biological variation (multiple true variants) vs bioinformatic contamination
- Consider stricter clustering parameters (`--min-identity` or `--inflation` in speconsense)
- May require manual curation or re-demultiplexing

**For p95diff > 0 only (outlier variation):**
- Review if sequence is critical for your analysis
- Often acceptable - may just be rare sequencing errors
- Consider the biological context and downstream use case

**For Merged with Small Components:**
- Review whether small components should have been filtered earlier
- Consider adjusting `--min-size` or `--min-cluster-ratio` in speconsense
- Consider stricter `--min-ric` threshold in speconsense-summarize
- Check if small components represent real low-abundance variants (sensical by themselves)

**Workflow Integration:**

The quality report is designed for efficient triage:
1. Scan from top to bottom - most critical issues appear first
2. Focus on HIGH PRIORITY sequences with p50diff > 0
3. Use component information to decide whether to review .raw files
4. Defer small component issues to lower priority unless critical to analysis

For high-throughput workflows (e.g., 100K sequences/year), this prioritization ensures human review time focuses on the most actionable quality issues.

### Additional Summarize Options

**Quality Filtering:**
```bash
speconsense-summarize --min-ric 5
```
- Filters out consensus sequences with fewer than the specified number of Reads in Consensus (RiC)
- Default is 3 - only sequences supported by at least 3 reads are processed
- Higher values provide more stringent quality control but may exclude valid low-abundance variants

**Variant Merging:**
```bash
# Basic SNP-only merging (default)
speconsense-summarize --merge-position-count 2

# Enable indel merging (up to 3bp indels)
speconsense-summarize --merge-position-count 3 --merge-indel-length 3

# Disable SNP merging (only merge identical sequences)
speconsense-summarize --merge-snp false

# Legacy parameter (still supported)
speconsense-summarize --snp-merge-limit 2  # Equivalent to --merge-position-count 2
```
- **Occurs within each HAC group** - merges variants that differ by small numbers of SNPs and/or indels
- Uses **MSA-based approach with SPOA**: evaluates all possible subsets to find the largest compatible group for merging
- Creates **IUPAC consensus sequences** with size-weighted majority voting at polymorphic positions
- **Order-independent**: produces identical results regardless of input order (unlike pairwise greedy approaches)
- Only merges variants with **identical primer sets** to maintain biological validity

**Merge Parameters:**
- `--merge-position-count N`: Maximum total SNP + indel positions allowed (default: 2)
- `--merge-indel-length N`: Maximum length of individual indels allowed (default: 0 = disabled)
- `--merge-snp`: Enable/disable SNP merging (default: True)
- `--merge-min-size-ratio R`: Minimum size ratio (smaller/larger) for merging clusters (default: 0.0 = disabled)
- `--snp-merge-limit N`: Legacy parameter, equivalent to `--merge-position-count` (deprecated)

**How it works:**
- Counts SNPs and indels separately during alignment
- Both must be within limits for merge to proceed
- Example: `--merge-position-count 3 --merge-indel-length 2` allows merging sequences with:
  - 3 SNPs + 0 indels ✓
  - 2 SNPs + 1 indel (≤2bp) ✓
  - 0 SNPs + 3 indels (each ≤2bp) ✓
  - 2 SNPs + 2 indels (each ≤2bp) ✗ (total=4 > 3)
  - 2 SNPs + 1 indel (3bp) ✗ (indel too long)

**Merged consensus tracking:**
- Merged sequences include `rawric` header field showing RiC values of merged .raw variants
- Example: `>sample-1 size=250 ric=250 rawric=100+89+61 snp=2`
- Helps trace which original clusters contributed to merged consensus

**Note**: This is distinct from the automatic homopolymer-aware merging that occurs during the main clustering step in speconsense

**Merge Size Ratio Filtering:**
```bash
speconsense-summarize --merge-min-size-ratio 0.1
```
- Prevents merging clusters with very different sizes (e.g., well-supported variant + poorly-supported variant)
- Ratio calculated as `smaller_size / larger_size` - must be ≥ threshold to merge
- Example: `--merge-min-size-ratio 0.1` means smaller cluster must be ≥10% size of larger
- Default is 0.0 (disabled) - all size combinations allowed
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

### Processing Workflow Summary

The complete speconsense-summarize workflow operates in this order:

1. **Load sequences** with RiC filtering (`--min-ric`)
2. **HAC variant grouping** by sequence identity to separate dissimilar sequences (`--group-identity`)
3. **Group filtering** to limit output groups (`--select-max-groups`)
4. **MSA-based variant merging** within each group (`--merge-position-count`, `--merge-indel-length`, `--merge-snp`, `--merge-min-size-ratio`)
5. **Variant selection** within each group (`--select-max-variants`, `--select-strategy`)
6. **Output generation** with customizable header fields (`--fasta-fields`) and full traceability

**Key architectural change**: HAC grouping now occurs BEFORE merging to prevent inappropriate merging of dissimilar sequences (e.g., contaminants with primary targets). Merging is then applied independently within each group using MSA-based consensus generation.

### Enhanced Logging and Traceability

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
- **Merge history**: tracks which original clusters were combined during variant merging
- **File lineage**: maintains connection between final outputs and original speconsense clusters
- **Read aggregation**: `FASTQ Files/` directory contains all reads that contributed to each final consensus
- **Pre-merge preservation**: `.raw` files preserve individual pre-merge variants with their original sequences and reads

This comprehensive logging allows users to understand exactly how the pipeline processed their data and make informed decisions about parameter tuning.

## Full Command Line Options

```
usage: speconsense.py [-h] [--augment-input AUGMENT_INPUT] [--algorithm {graph,greedy}] [--min-identity MIN_IDENTITY]
                     [--inflation INFLATION] [--min-size MIN_SIZE] [--min-cluster-ratio MIN_CLUSTER_RATIO]
                     [--max-sample-size MAX_SAMPLE_SIZE] [--presample PRESAMPLE]
                     [--k-nearest-neighbors K_NEAREST_NEIGHBORS] [--primers PRIMERS]
                     [-O OUTPUT_DIR] [--stability-trials STABILITY_TRIALS] [--stability-sample STABILITY_SAMPLE]
                     [--disable-stability] [--disable-homopolymer-equivalence]
                     [--orient-mode {skip,keep-all,filter-failed}]
                     [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--version]
                     input_file

Markov Clustering-based clustering of nanopore amplicon reads

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
  --min-size MIN_SIZE   Minimum cluster size (default: 5, 0 to disable)
  --min-cluster-ratio MIN_CLUSTER_RATIO
                        Minimum size ratio between a cluster and the largest cluster (default: 0.2, 0 to disable)
  --max-sample-size MAX_SAMPLE_SIZE
                        Maximum cluster size for consensus (default: 500)
  --presample PRESAMPLE
                        Presample size for initial reads (default: 1000, 0 to disable)
  --k-nearest-neighbors K_NEAREST_NEIGHBORS
                        Number of nearest neighbors for graph construction (default: 5)
  --primers PRIMERS     FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)
  -O OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for all files (default: clusters)
  --stability-trials STABILITY_TRIALS
                        Number of sampling trials to assess stability (default: 100)
  --stability-sample STABILITY_SAMPLE
                        Size of stability samples (default: 20)
  --disable-stability   Disable stability assessment
  --disable-homopolymer-equivalence
                        Disable homopolymer equivalence in cluster merging (only merge identical sequences)
  --orient-mode {skip,keep-all,filter-failed}
                        Sequence orientation mode: skip (default, no orientation), keep-all (orient but keep failed),
                        or filter-failed (orient and remove failed). Requires primers file with position annotations.
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set logging level (default: INFO)
  --version             Show program's version number and exit
```

### speconsense-summarize Options

```
usage: speconsense-summarize [-h] [--min-ric MIN_RIC] [--source SOURCE]
                             [--summary-dir SUMMARY_DIR]
                             [--fasta-fields FASTA_FIELDS] [--merge-snp]
                             [--merge-indel-length MERGE_INDEL_LENGTH]
                             [--merge-position-count MERGE_POSITION_COUNT]
                             [--merge-min-size-ratio MERGE_MIN_SIZE_RATIO]
                             [--group-identity GROUP_IDENTITY]
                             [--select-max-variants SELECT_MAX_VARIANTS]
                             [--select-max-groups SELECT_MAX_GROUPS]
                             [--select-strategy {size,diversity}]
                             [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Process Speconsense output with advanced variant handling.

options:
  -h, --help            show this help message and exit
  --min-ric MIN_RIC     Minimum Reads in Consensus (RiC) threshold (default: 3)
  --source SOURCE       Source directory containing Speconsense output (default: clusters)
  --summary-dir SUMMARY_DIR
                        Output directory for summary files (default: __Summary__)
  --fasta-fields FASTA_FIELDS
                        FASTA header fields to output. Can be: (1) a preset name (default,
                        minimal, qc, full, id-only), (2) comma-separated field names (size,
                        ric, length, rawric, snp, p50diff, p95diff, primers, group, variant),
                        or (3) a combination of presets and fields (e.g., minimal,qc or
                        minimal,p50diff,p95diff). Duplicates removed, order preserved left to
                        right. Default: default
  --merge-snp           Enable SNP-based merging (default: True)
  --merge-indel-length MERGE_INDEL_LENGTH
                        Maximum length of individual indels allowed in merging (default: 0 =
                        disabled)
  --merge-position-count MERGE_POSITION_COUNT
                        Maximum total SNP+indel positions allowed in merging (default: 2)
  --merge-min-size-ratio MERGE_MIN_SIZE_RATIO
                        Minimum size ratio (smaller/larger) for merging clusters (default:
                        0.0 = disabled)
  --group-identity GROUP_IDENTITY, --variant-group-identity GROUP_IDENTITY
                        Identity threshold for variant grouping using HAC (default: 0.9)
  --select-max-variants SELECT_MAX_VARIANTS, --max-variants SELECT_MAX_VARIANTS
                        Maximum number of additional variants to output per group (default:
                        -1 = no limit)
  --select-max-groups SELECT_MAX_GROUPS, --max-groups SELECT_MAX_GROUPS
                        Maximum number of groups to output per specimen (default: -1 = all
                        groups)
  --select-strategy {size,diversity}, --variant-selection {size,diversity}
                        Variant selection strategy: size or diversity (default: size)
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level
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
