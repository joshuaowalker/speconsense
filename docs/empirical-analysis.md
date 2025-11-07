# Empirical Error Analysis with Speconsense-Analyze

This guide explains how to use `speconsense-analyze` to estimate real-world Nanopore error rates and characteristics by comparing raw reads against their consensus sequences.

## Overview

While [synthetic testing](synthetic-testing.md) uses `speconsense-synth` to generate data with **known** error rates, empirical analysis flips this approach: it uses real-world data where the **consensus sequences serve as proxies for biological truth**, then measures how raw reads deviate from their consensus.

This approach allows you to:
- Estimate actual Nanopore error rates in your data
- Validate synthetic testing assumptions against real data
- Understand error type distributions (insertions vs deletions vs substitutions)
- Detect systematic error patterns (read-end degradation, homopolymer issues)
- Correlate read depth (RiC) with empirical consensus quality

## Quick Start

```bash
# Basic usage - analyze a speconsense output directory
speconsense-analyze /path/to/speconsense/output

# With options
speconsense-analyze /path/to/speconsense/output \
  --output-dir my_analysis \
  --min-ric 5 \
  --use-sampled \
  --per-read
```

## How It Works

### Core Concept

For each cluster in your speconsense output:
1. Load the consensus sequence from `{sample}-all.fasta`
2. Load the raw reads from `cluster_debug/{sample}-c{N}-RiC{M}-reads.fastq`
3. Align each raw read to its consensus using edlib
4. Calculate edit distances and classify error types
5. Aggregate statistics per-cluster and dataset-wide

### Assumptions and Limitations

**Key Assumption:** The consensus sequence represents biological truth. This is reasonable when:
- Clusters have sufficient read depth (RiC ≥ 10)
- Error correction via voting has converged (modern chemistry: N≥3-4 typically sufficient)
- No contamination or heterogeneity (single biological source)

**Limitations:**
- **Errors in consensus are invisible**: If consensus itself has errors, they won't be detected
- **Low RiC clusters**: With RiC < 5, consensus may not be fully error-corrected
- **Heterogeneous clusters**: Mixed populations or contamination confound results
- **Systematic errors**: If all reads share the same systematic error, consensus will too

**Best practices:**
- Filter to RiC ≥ 10 for most reliable results (use `--min-ric 10`)
- Cross-validate with synthetic testing
- Check stability metrics (p50diff, p95diff) - elevated values suggest heterogeneity
- Compare results across multiple specimens to identify outliers

## Command Line Options

### Required Arguments

- `output_dir`: Path to speconsense output directory (containing `cluster_debug/` and `*-all.fasta` files)

### Optional Arguments

- `--output-dir, -o DIR`: Output directory for analysis results (default: `analysis_results`)
- `--use-sampled`: Use sampled reads (`*-sampled.fastq`) instead of full cluster (`*-reads.fastq`)
  - Use this if you want to match the exact reads used for consensus generation
  - Otherwise, analyzes all reads in the cluster (more comprehensive)
- `--min-ric N`: Minimum RiC threshold for analysis (default: 0)
  - Recommended: 10 for reliable error estimates
- `--per-read`: Output detailed per-read alignment data (can create large files)
- `--log-level LEVEL`: Logging verbosity (DEBUG, INFO, WARNING, ERROR, CRITICAL)

## Output Files

All outputs are written to the specified `--output-dir` (default: `analysis_results/`):

### 1. `cluster_statistics.csv`

Per-cluster statistics, one row per cluster. Columns:

| Column | Description |
|--------|-------------|
| `sample_name` | Sample name with cluster number (e.g., `sample-c1`) |
| `cluster_num` | Cluster number |
| `ric` | Reads in Consensus (from speconsense) |
| `size` | Total cluster size (from speconsense) |
| `num_reads_analyzed` | Number of reads successfully aligned |
| `mean_edit_distance` | Mean edit distance from reads to consensus |
| `median_edit_distance` | Median edit distance |
| `p95_edit_distance` | 95th percentile edit distance |
| `mean_error_rate` | Mean error rate (edit_distance / consensus_length) |
| `median_error_rate` | Median error rate |
| `insertion_rate` | Proportion of errors that are insertions |
| `deletion_rate` | Proportion of errors that are deletions |
| `substitution_rate` | Proportion of errors that are substitutions |
| `consensus_length` | Length of consensus sequence |
| `mean_read_length` | Mean length of raw reads |
| `has_contiguous_errors` | Boolean: detected contiguous high-error regions |
| `p50_diff` | Stability metric from speconsense (if available) |
| `p95_diff` | Stability metric from speconsense (if available) |

**Use this file to:**
- Identify problematic clusters (high error rates, contiguous errors)
- Analyze RiC vs quality relationships
- Filter clusters for downstream analysis

### 2. `analysis_report.txt`

Human-readable summary of dataset-wide statistics:

- **Dataset Overview**: Number of clusters and reads analyzed
- **Error Rate Estimates**: Mean, median, and percentile error rates
- **Error Type Distribution**: Proportion of insertions, deletions, substitutions
- **Positional Error Patterns**: Frequency of contiguous high-error regions
- **Read Depth vs Quality**: Error rates binned by RiC
- **Comparison to Expectations**: How observed rates compare to synthetic testing predictions

**Use this file to:**
- Quickly assess overall data quality
- Compare to synthetic testing expectations (1-2% for modern chemistry)
- Identify systematic issues

### 3. `read_distances.csv` (Optional, `--per-read`)

Detailed per-read alignment data. Columns:

- `sample_name`, `cluster_num`: Cluster identification
- `read_id`: Read identifier from FASTQ
- `read_length`: Length of read
- `edit_distance`: Edit distance to consensus
- `error_rate`: Per-read error rate
- `insertions`, `deletions`, `substitutions`: Error type counts

**Use this file to:**
- Investigate specific outlier reads
- Build custom visualizations
- Detailed error modeling

**Warning:** Can be very large for high-coverage datasets. Only use when necessary.

## Interpreting Results

### Error Rate Estimates

Compare observed error rates to expectations from synthetic testing:

| Observed Mean Error Rate | Interpretation |
|-------------------------|----------------|
| < 2% | Excellent - consistent with modern R10.4.1 + V14 chemistry |
| 2-5% | Good - slightly above modern expectations, but acceptable |
| 5-10% | Fair - may indicate older chemistry (R9.4) or quality issues |
| > 10% | Poor - investigate data quality, may have systematic issues |

**Note:** These thresholds assume consensus sequences are accurate. With low RiC or heterogeneous clusters, apparent error rates will be inflated.

### Error Type Distribution

Expected from literature (see [Nanopore error characteristics](synthetic-testing.md#understanding-nanopore-error-characteristics)):

- **Deletions**: Most common (~40-50% of errors)
- **Substitutions**: Moderate (~25-35% of errors)
- **Insertions**: Least common (~20-30% of errors)

**Deviations from expected:**
- **Equal distribution** (~33% each): Suggests `speconsense-synth`-like uniform error model
- **High insertions** (>40%): Unusual, may indicate systematic issue
- **Very high deletions** (>60%): Consistent with homopolymer-rich sequences

### Contiguous High-Error Regions

The `has_contiguous_errors` flag indicates whether reads show localized regions of high error density (>20% errors in a 50bp window).

**Common causes:**
- **Read-end degradation**: Nanopore reads often have lower quality at ends
- **Homopolymer runs**: Long stretches of single nucleotide prone to length errors
- **Low-quality reads**: Individual reads with poor sequencing quality

**Interpretation:**
- **<10% of clusters**: Normal - occasional low-quality reads
- **10-30% of clusters**: Moderate - some systematic quality issues
- **>30% of clusters**: High - investigate sequencing quality, consider stricter filtering

### Read Depth vs Quality

The analysis bins clusters by RiC and reports mean error rates per bin.

**Expected pattern (from synthetic testing):**
- **RiC < 5**: Higher apparent error rate (consensus not fully converged)
- **RiC 5-20**: Error rate stabilizes
- **RiC > 20**: Error rate plateaus at true Nanopore error rate

**If you see:**
- **No trend**: Good - consensus quality not dependent on depth (sufficient reads throughout)
- **Improving with RiC**: Expected - validates read depth recommendations
- **Degrading with RiC**: Unexpected - investigate, may indicate technical issue

### Comparing to Stability Metrics

The output includes `p50_diff` and `p95_diff` from speconsense stability assessment.

**Expected relationship:**
- **Low error rate + low p50diff**: Excellent - homogeneous, high-quality cluster
- **Low error rate + high p50diff**: Warning - may indicate heterogeneity or contamination
- **High error rate + high p50diff**: Poor - likely mixed population or low quality
- **High error rate + low p50diff**: Unusual - investigate (may be real variant vs reference)

## Worked Examples

### Example 1: Assessing Data Quality

**Scenario:** You've sequenced 50 specimens with modern ONT chemistry and want to confirm data quality matches expectations.

```bash
# Run analysis with RiC filtering
speconsense-analyze /path/to/speconsense/output \
  --min-ric 10 \
  --output-dir quality_check

# Review summary
cat quality_check/analysis_report.txt
```

**Interpretation:**
```
Error Rate Estimates:
  Mean error rate: 1.87%
  Median error rate: 1.65%
  95th percentile: 3.12%
  5th percentile: 0.95%

Error Type Distribution:
  Mean insertion rate: 28.3%
  Mean deletion rate: 42.7%
  Mean substitution rate: 29.0%
```

**Conclusion:** ✓ Consistent with modern chemistry (< 2% mean error, deletions most common)

### Example 2: Validating Synthetic Testing Assumptions

**Scenario:** Synthetic tests used 2% error rate. Does real data match?

```bash
# Analyze both synthetic and real data
speconsense-analyze synthetic_output --output-dir synthetic_analysis
speconsense-analyze real_output --output-dir real_analysis

# Compare analysis_report.txt files
```

**Expected:** Real data error rates should be ≤ synthetic testing assumptions (if synthetic testing used conservative parameters).

### Example 3: Investigating RiC Thresholds

**Scenario:** Synthetic testing suggests N≥3 for reliable consensus. Does real data validate this?

```bash
# Analyze with no RiC filter to see all clusters
speconsense-analyze /path/to/output --min-ric 0 --output-dir ric_analysis

# Examine cluster_statistics.csv
# Plot mean_error_rate vs ric in spreadsheet or R/Python
```

**Look for:**
- Do clusters with RiC < 5 have higher error rates?
- At what RiC does error rate stabilize?
- Are there outliers (high RiC but high error rate)?

### Example 4: Detecting Problematic Clusters

**Scenario:** Identify clusters with unusual error patterns for manual review.

```bash
speconsense-analyze /path/to/output --output-dir problem_check

# Filter cluster_statistics.csv for issues:
# - mean_error_rate > 0.05 (5% error)
# - has_contiguous_errors == True
# - high p50_diff (> 2.0)

# Investigate flagged clusters in speconsense output
```

## Positional Error Analysis (Future Enhancement)

A future enhancement could add `error_positions.csv` with per-position error density across the consensus, enabling:
- Identification of problematic sequence regions (homopolymers, GC extremes)
- Read-end quality visualization
- Systematic error pattern detection

Currently, positional analysis is limited to contiguous error region detection.

## Comparison to Synthetic Testing

| Aspect | Synthetic Testing | Empirical Analysis |
|--------|------------------|-------------------|
| **Ground truth** | Known (input reference) | Assumed (consensus) |
| **Error rate** | Controlled | Measured |
| **Advantages** | Precise, reproducible, controlled | Real-world, validates assumptions |
| **Limitations** | Simplified error model | Consensus errors invisible |
| **Best for** | Algorithm development, parameter tuning | Quality assessment, validation |

**Recommendation:** Use both approaches:
1. **Synthetic testing** to establish expectations and parameter ranges
2. **Empirical analysis** to validate real-world performance

## Integration with Pipeline

### Typical Workflow

```bash
# 1. Run speconsense
speconsense input.fastq -o output --min-size 5

# 2. Run empirical analysis
speconsense-analyze output --output-dir analysis --min-ric 10

# 3. Review quality
cat analysis/analysis_report.txt

# 4. If quality is good, proceed with summarization
speconsense-summarize --source output --summary-dir results

# 5. Optionally analyze summarized data
speconsense-analyze output --output-dir final_analysis --min-ric 3
```

### Quality Filtering Based on Analysis

Use cluster statistics to filter:

```python
import pandas as pd

# Load cluster statistics
stats = pd.read_csv('analysis/cluster_statistics.csv')

# Define quality filters
good_clusters = stats[
    (stats['mean_error_rate'] < 0.03) &  # <3% error
    (stats['ric'] >= 10) &                # Sufficient depth
    (~stats['has_contiguous_errors']) &   # No error clustering
    (stats['p50_diff'] < 2.0)             # Stable consensus
]

# Use good_clusters to filter downstream
```

## Troubleshooting

### "No consensus sequences found"

**Cause:** Output directory doesn't contain `*-all.fasta` files

**Solution:**
- Check you're pointing to speconsense output directory
- Ensure speconsense completed successfully
- Look for files matching pattern `*-all.fasta`

### "No cluster files found"

**Cause:** Missing `cluster_debug/` directory or no FASTQ files

**Solution:**
- Ensure speconsense was run (not just speconsense-summarize)
- Check that `cluster_debug/` subdirectory exists
- If using `--use-sampled`, ensure `*-sampled.fastq` files exist

### "No matching clusters found"

**Cause:** Cluster numbers in FASTQ don't match consensus FASTA, or all filtered by `--min-ric`

**Solution:**
- Check file naming consistency
- Try `--min-ric 0` to include all clusters
- Verify speconsense output is complete

### Very high error rates (>10%)

**Possible causes:**
1. Low RiC - consensus not error-corrected (try `--min-ric 10`)
2. Heterogeneous clusters - mixing multiple variants
3. Actual poor sequencing quality
4. Wrong consensus matched to reads

**Diagnosis:**
- Check RiC distribution in `cluster_statistics.csv`
- Review `p50_diff` / `p95_diff` - elevated values suggest heterogeneity
- Run synthetic test to confirm pipeline working correctly
- Inspect problem clusters manually

### Error type distribution doesn't match literature

**If equal distribution (~33% each):**
- May indicate highly GC-balanced sequences
- Unusual but not necessarily wrong

**If very skewed:**
- Check for sequence composition biases (e.g., homopolymer-rich)
- Verify alignment quality (CIGAR string parsing)
- Consider biological factors (e.g., amplicon design)

## Next Steps

After running empirical analysis:

1. **Compare to synthetic testing**: Do real error rates match synthetic testing assumptions?
2. **Adjust parameters**: Use RiC vs quality curves to set `--min-ric` thresholds
3. **Filter data**: Exclude problematic clusters from downstream analysis
4. **Iterate**: Re-run speconsense with adjusted parameters if needed
5. **Document**: Record observed error rates for your sequencing platform/chemistry

## References

See [Synthetic Testing Guide](synthetic-testing.md) for:
- Expected error rates by chemistry
- Synthetic testing methodology
- Interpretation guidelines
- Literature references on Nanopore error characteristics
