# Customizing FASTA Header Fields

This guide explains how to customize which metadata fields appear in FASTA headers when using `speconsense-summarize`, allowing you to tailor output for different downstream tools and workflows.

## Overview

By default, `speconsense-summarize` includes comprehensive metadata in FASTA headers:

```
>ONT01.01-A01-...-1 size=638 ric=638 rawric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC
```

The `--fasta-fields` option lets you customize which fields appear, from minimal ID-only headers to comprehensive metadata.

## Why Customize Headers?

**Common use cases:**

1. **Downstream tool compatibility** - Some tools expect minimal headers or specific formats
2. **Readability** - Include only essential fields for cleaner output
3. **File size optimization** - Minimal headers reduce file size for large datasets
4. **Quality control** - Include stability metrics for assessing consensus quality
5. **Workflow-specific needs** - Different analyses need different metadata

## Design Principles

The field customization system follows these principles:

1. **Field codes = field names** - No translation needed (e.g., `rawric` in command = `rawric=` in output)
2. **Preset composition** - Combine presets with union semantics: `--fasta-fields minimal,qc` expands to all fields from both
3. **Conditional output** - Fields with null/zero/empty values are automatically omitted for clean headers
4. **Order preservation** - Fields appear in the order specified (left to right)

## Quick Start

### Using Presets

```bash
# Default preset (current behavior, backward compatible)
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

# ID only (no metadata fields)
speconsense-summarize --fasta-fields id-only
# Output: >sample-1
```

### Custom Field Selection

```bash
# Specify individual fields (comma-separated, field codes = field names)
speconsense-summarize --fasta-fields size,ric,primers
# Output: >sample-1 size=638 ric=638 primers=5'-ITS1F,3'-ITS4_RC

# Just stability metrics
speconsense-summarize --fasta-fields p50diff,p95diff
# Output: >sample-1 p50diff=0.0 p95diff=2.0
```

### Combining Presets and Fields

```bash
# Combine presets (union semantics)
speconsense-summarize --fasta-fields minimal,qc
# Output: >sample-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=2.0

# Preset + individual fields
speconsense-summarize --fasta-fields minimal,p50diff,p95diff
# Output: >sample-1 size=638 ric=638 p50diff=0.0 p95diff=2.0

# Add group/variant info to default
speconsense-summarize --fasta-fields default,group,variant
# Output: >sample-1 size=638 ric=638 rawric=333+305 snp=1 primers=... group=1
# (variant field only present for within-group variants like sample-1.v1)
```

## Available Fields

All available metadata fields that can be included in headers:

| Field | Description | Example | Always Present? |
|-------|-------------|---------|-----------------|
| `size` | Total reads across merged variants | `size=638` | Yes |
| `ric` | Reads in Consensus (total) | `ric=638` | Yes |
| `length` | Sequence length in bases | `length=589` | Yes |
| `rawric` | RiC values of .raw source variants (pre-merge) | `rawric=333+305` | Only when merged |
| `snp` | Number of IUPAC ambiguity positions | `snp=1` | Only when >0 |
| `p50diff` | Median (p50) edit distance from stability assessment | `p50diff=0.0` | Only when available |
| `p95diff` | 95th percentile edit distance from stability | `p95diff=2.0` | Only when available |
| `primers` | Detected primer names | `primers=5'-ITS1F,3'-ITS4_RC` | Only when detected |
| `group` | Variant group number | `group=1` | Yes |
| `variant` | Variant identifier within group | `variant=v1` | Only for variants |

**Notes:**
- Conditional fields (marked "Only when...") are automatically omitted if not applicable
- Stability fields (`p50diff`, `p95diff`) come from speconsense's stability assessment and are available if stability was not disabled
- `group` and `variant` fields extract info from the sample name (e.g., `sample-1` → `group=1`, `sample-1.v1` → `variant=v1`)
- The `rawric` field shows the quality distribution of merged variants - see [Understanding RiC and Merging](understanding-ric-and-merging.md)

## Preset Reference

| Preset | Fields Included | Use Case |
|--------|-----------------|----------|
| `default` | `size,ric,rawric,snp,primers` | General use, current behavior (backward compatible) |
| `minimal` | `size,ric` | Quick inspection, essential info only |
| `qc` | `size,ric,length,p50diff,p95diff` | Quality control with stability metrics |
| `full` | `size,ric,length,rawric,snp,p50diff,p95diff,primers` | Complete metadata (excludes group/variant) |
| `id-only` | *(none)* | Clean IDs only, no metadata fields |

**Preset composition:**
- Presets can be combined: `minimal,qc` → `size,ric,length,p50diff,p95diff` (union of both)
- Mix presets and individual fields: `minimal,primers` → `size,ric,primers`
- Duplicates are automatically removed, order preserved left to right

## Common Workflows

### For Downstream Tool Compatibility

Many sequence analysis tools expect simple headers:

```bash
# BLAST-like minimal headers
speconsense-summarize --fasta-fields minimal

# ID only for tools that parse metadata separately
speconsense-summarize --fasta-fields id-only
```

### For Quality Control

Include stability metrics to assess consensus reliability:

```bash
# QC-focused headers
speconsense-summarize --fasta-fields qc
# Output: >sample-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=2.0

# Add stability to default fields
speconsense-summarize --fasta-fields default,p50diff,p95diff
```

Sequences with higher-than-expected stability differences should be checked for bioinformatic contamination or biological variation. See speconsense documentation for interpreting stability metrics.

### For File Size Optimization

Large datasets with thousands of sequences benefit from minimal headers:

```bash
# Minimal headers reduce file size
speconsense-summarize --fasta-fields minimal

# Or just essential quality info
speconsense-summarize --fasta-fields size,ric,p50diff
```

### For Publication/GenBank Submission

Include comprehensive metadata for transparency:

```bash
# Full metadata for archival
speconsense-summarize --fasta-fields full

# Or default with added quality metrics
speconsense-summarize --fasta-fields default,p50diff,p95diff
```

### For Manual Review

Include grouping information to understand variant structure:

```bash
# Show variant grouping structure
speconsense-summarize --fasta-fields default,group,variant
# Output for merged: >sample-1 size=638 ric=638 rawric=333+305 snp=1 group=1
# Output for variant: >sample-1.v1 size=100 ric=100 group=1 variant=v1
```

## Field Details

### Size vs RiC

- **`size`** - Total reads in the cluster(s) before any sampling
- **`ric`** - Reads actually used for consensus generation (may be sampled down)

For single (non-merged) variants: `size` ≥ `ric` (equality when cluster smaller than `--max-sample-size`)

For merged variants: `ric` can exceed `--max-sample-size` because it's the sum of RiC values from multiple pre-merge variants. See [Understanding RiC and Merging](understanding-ric-and-merging.md) for details.

### The rawric Field

Present only when variants were merged:

```
rawric=333+305
```

This shows:
- Two variants were merged to create this consensus
- The first variant had ric=333, the second had ric=305
- Values are ordered largest-first
- Total ric = 333+305 = 638

Use `rawric` to assess merge quality:
- **Balanced** (good): `rawric=60+55+50` - all variants have similar support
- **Dominated** (review): `rawric=100+8+5` - one strong variant, others weak

See [Understanding RiC and Merging](understanding-ric-and-merging.md) for interpreting merged consensuses.

### Stability Fields

Fields `p50diff` and `p95diff` come from speconsense's stability assessment:

- **`p50diff`** - Median edit distance between subsample consensuses and full consensus
- **`p95diff`** - 95th percentile edit distance

Lower values indicate more stable consensuses. Values are edit distances, so:
- `p50diff=0.0` - very stable (median subsample matches full consensus exactly)
- `p50diff=2.0` - typical for good quality (median differs by ~2 bases)
- `p50diff=10.0+` - less stable, review for contamination or biological variation

These fields are available if:
1. Speconsense was run with stability assessment (default behavior)
2. You use presets/fields that include them (`qc`, `full`, or manual specification)

### Length Field

Simple sequence length in bases:
```
length=589
```

Useful for:
- Filtering by amplicon length
- Quality control (expected length ranges)
- Size distribution analysis

### Primers Field

Shows detected primers:
```
primers=5'-ITS1F,3'-ITS4_RC
```

Format: `<position>-<primer_name>` where position is `5'` (forward) or `3'` (reverse). The `_RC` suffix indicates the primer was detected as reverse complement.

### Group and Variant Fields

Extract grouping information from sample names:

```
# Main merged consensus: sample-1
group=1

# Within-group variant: sample-1.v1
group=1 variant=v1
```

These fields make the variant structure explicit without parsing sample names. Useful when you need to programmatically identify variant relationships.

## Examples

### Example 1: Minimal for BLAST

```bash
speconsense-summarize --fasta-fields minimal --summary-dir blast_ready
```

Output:
```
>specimen-A-1 size=638 ric=638
ATCGATCG...
>specimen-B-1 size=203 ric=203
GCTAGCTA...
```

Clean, simple headers compatible with most sequence analysis tools.

### Example 2: QC Report

```bash
speconsense-summarize --fasta-fields qc --summary-dir qc_report
```

Output:
```
>specimen-A-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=2.0
ATCGATCG...
```

Includes all quality metrics for assessing consensus reliability.

### Example 3: Custom Combination

```bash
speconsense-summarize --fasta-fields size,ric,rawric,p50diff,primers
```

Output for merged variant:
```
>specimen-A-1 size=638 ric=638 rawric=333+305 p50diff=0.0 primers=5'-ITS1F,3'-ITS4_RC
ATCGATCG...
```

Output for single variant (no rawric since not merged):
```
>specimen-B-1 size=203 ric=203 p50diff=1.2 primers=5'-ITS1F,3'-ITS4_RC
GCTAGCTA...
```

### Example 4: Preset Composition

```bash
speconsense-summarize --fasta-fields minimal,qc
```

Output (union of minimal + qc):
```
>specimen-A-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=2.0
ATCGATCG...
```

This combines essential metrics (size, ric) with quality control fields (length, stability).

## Technical Notes

### Field Order

Fields appear in the order specified:

```bash
# Order matters
speconsense-summarize --fasta-fields primers,ric,size
# Output: >sample-1 primers=... ric=638 size=638

# Different order
speconsense-summarize --fasta-fields size,ric,primers
# Output: >sample-1 size=638 ric=638 primers=...
```

### Duplicate Removal

When combining presets, duplicates are automatically removed:

```bash
# minimal = size,ric
# qc = size,ric,length,p50diff,p95diff
speconsense-summarize --fasta-fields minimal,qc
# Result: size,ric,length,p50diff,p95diff (size,ric not duplicated)
```

### Conditional Fields

Some fields are only output when applicable:

- `rawric` - only for merged variants
- `snp` - only when >0
- `p50diff`, `p95diff` - only when stability data available
- `primers` - only when primers detected
- `variant` - only for within-group variants (e.g., `.v1`, `.v2`)

This keeps headers clean and avoids clutter.

### Error Handling

Invalid field names produce clear error messages:

```bash
speconsense-summarize --fasta-fields size,invalid_field,ric
# Error: Unknown preset or field name: 'invalid_field'
# Available presets: default, full, id-only, minimal, qc
# Available fields: group, length, p50diff, p95diff, primers, rawric, ric, size, snp, variant
```

## Backward Compatibility

The default behavior is unchanged:

```bash
# These are equivalent
speconsense-summarize
speconsense-summarize --fasta-fields default
```

Both produce the same output format as previous versions. Existing scripts and workflows continue to work without modification.

## Related Documentation

- **[Understanding RiC and Merging](understanding-ric-and-merging.md)** - Explains RiC accumulation, rawric interpretation, and variant merging philosophy
- **README.md** - See "Customizing FASTA Header Fields" section for basic usage
- **README.md** - See "FASTA Header Metadata" section for field definitions

## Summary

The `--fasta-fields` option provides flexible control over FASTA header metadata:

- **Presets** for common use cases (minimal, qc, full, etc.)
- **Custom field selection** with field codes matching field names
- **Preset composition** to combine multiple presets
- **Conditional output** to keep headers clean
- **Backward compatible** with default preset

Choose the right combination for your workflow - from minimal ID-only headers to comprehensive metadata with quality metrics.
