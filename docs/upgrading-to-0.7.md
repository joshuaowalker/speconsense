# Upgrading to Speconsense 0.7

This document describes changes in Speconsense 0.7 that may affect your workflows.

## Philosophy Clarification

Version 0.7 clarifies the separation of concerns between the two tools:

**speconsense** extracts distinct biological sequences from raw reads within a single specimen. It separates true biological variation (alleles, loci, contaminants) from sequencing noise by clustering reads and generating consensus sequences. Philosophy: *split rather than conflate* — better to produce separate clusters that can be merged downstream than to lose real biological distinctions.

**speconsense-summarize** curates and consolidates consensus sequences for downstream analysis. It handles residual ambiguity that consensus generation cannot resolve — primarily minor SNP differences that may represent the same biological entity — groups related variants, selects representatives, and can analyze patterns across a run. Philosophy: *conservative simplification* — merge only what is clearly equivalent, preserve distinctions when uncertain.

## Default Changes in speconsense

### `--min-cluster-ratio`: 0.20 → 0.01

**What changed:** The minimum cluster size ratio (relative to the largest cluster) decreased from 20% to 1%.

**Why:** Statistical analysis shows that random sequencing errors (1-2%) cannot produce artifact clusters — the 90% identity clustering threshold filters them out. Small clusters represent real biological signal: minor variants, contaminants, or chimeric sequences. The previous 0.20 threshold was discarding this signal before speconsense-summarize could curate it.

**Impact:** You will see more small clusters in your output. These can be filtered during post-processing with speconsense-summarize using `--min-ric` or `--select-max-variants`.

**To restore previous behavior:** Use `--min-cluster-ratio 0.2`

### `--min-variant-frequency`: 0.20 → 0.10

**What changed:** The minimum allele frequency for variant phasing (cluster splitting) decreased from 20% to 10%, matching the ambiguity calling threshold.

**Why:** With the more permissive `--min-cluster-ratio`, splitting at lower frequencies is safe — small split clusters won't be discarded. The unified 10% threshold is easier to remember and aligns with intuition about when variation is "significant."

**Impact:** You may see more variant splitting in heterogeneous samples. Variants between 10-20% frequency that previously received IUPAC codes will now be split into separate clusters.

**To restore previous behavior:** Use `--min-variant-frequency 0.2`

### Threshold Summary

| Option | Old Default | New Default | Rationale |
|--------|-------------|-------------|-----------|
| `--min-cluster-ratio` | 0.20 | 0.01 | Preserve small clusters for downstream curation |
| `--min-variant-frequency` | 0.20 | 0.10 | Unify with ambiguity threshold |
| `--min-variant-count` | 5 | 5 | (unchanged) Must match `--min-size` to ensure viable clusters |
| `--min-ambiguity-frequency` | 0.10 | 0.10 | (unchanged) |
| `--min-ambiguity-count` | 3 | 3 | (unchanged) Can be lower since no cluster creation |

## Default Changes in speconsense-summarize

### `--merge-min-size-ratio`: 0.0 → 0.1

**What changed:** The minimum size ratio for merging clusters increased from 0.0 (disabled) to 0.1 (10%).

**Why:** Aligns with the 10% threshold used throughout the pipeline (variant phasing, ambiguity calling). Prevents poorly-supported clusters from introducing IUPAC ambiguities into well-supported consensus sequences. When manually inspecting merged FASTQ files, a 10% threshold is intuitive.

**Impact:** Very small clusters (less than 10% the size of the larger cluster) will no longer merge. They will remain as separate variants in the output.

**To restore previous behavior:** Use `--merge-min-size-ratio 0`

### Threshold Summary

| Option | Old Default | New Default | Rationale |
|--------|-------------|-------------|-----------|
| `--merge-min-size-ratio` | 0.0 | 0.1 | Prevent poorly-supported variants from affecting high-confidence sequences |

## New Features

### Scalability Mode (speconsense and speconsense-summarize)

Large datasets (1000+ sequences) now automatically use vsearch-based acceleration for pairwise comparisons, reducing complexity from O(n²) to approximately O(n log n).

**New options:**
- `--scale-threshold N` — Sequence count threshold for scalable mode (default: 1000, 0 to disable)
- `--threads N` — Thread count for internal parallelism (default: 1 for speconsense, 0/auto for speconsense-summarize)

**Requirements:** Install vsearch via `conda install bioconda::vsearch`. Falls back to brute-force if unavailable.

### Early Filtering (speconsense)

New `--enable-early-filter` option applies size filtering before variant phasing, improving performance for large datasets by skipping expensive phasing on clusters that will be filtered anyway.

### Discard Collection (speconsense)

New `--collect-discards` option writes all discarded reads (outliers and filtered clusters) to `cluster_debug/{sample}-discards.fastq` for inspection.

### Profile System (speconsense and speconsense-summarize)

New profile system allows saving and reusing parameter configurations for different workflows.

**Usage:**
```bash
# List available profiles
speconsense --list-profiles
speconsense-summarize --list-profiles

# Use a profile
speconsense input.fastq -p herbarium
speconsense-summarize -p strict

# CLI arguments override profile values
speconsense input.fastq -p herbarium --min-size 10
```

**Bundled profiles:**
- `herbarium` — High-recall settings for degraded DNA and type specimens
- `specimens` — Balanced settings for fresh tissue (matches defaults)
- `strict` — High-precision settings for confident results

**Custom profiles:**
On first use, an `example.yaml` template is created in `~/.config/speconsense/profiles/`. Copy and modify it to create your own profiles:

```bash
cd ~/.config/speconsense/profiles
cp example.yaml my-workflow.yaml
# Edit my-workflow.yaml with your settings
```

Profiles include version compatibility checking — profiles created for older versions will error with a helpful message when used with newer speconsense versions.

## Migration Checklist

1. **Review your `--min-cluster-ratio` usage:**
   - If you were using `--min-cluster-ratio 0` to see all clusters, the new default (0.01) should work well
   - If you were relying on 0.20 to filter contaminants, consider using speconsense-summarize's `--min-ric` instead

2. **Review your variant thresholds:**
   - If you have workflows that depend on the 10-20% "IUPAC zone," you may need to adjust
   - The unified 10% threshold means more splitting, less IUPAC marking

3. **Consider enabling scalability for large datasets:**
   - Install vsearch if processing specimens with 1000+ reads
   - The speedup is significant for large datasets

4. **Update any scripts that parse output:**
   - More clusters may be produced with the new defaults
   - Cluster numbering may differ due to increased splitting
