# Customizing FASTA Header Fields

This guide explains how to customize which metadata fields appear in FASTA headers when using `speconsense-summarize`, allowing you to tailor output for different downstream tools and workflows.

## Overview

By default, `speconsense-summarize` includes a comprehensive but readable metadata payload in FASTA headers:

```
>ONT01.01-A01-...-1.v1 size=638 ric=638 rawric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC
```

The `--fasta-fields` option lets you customize which fields appear, from minimal ID-only headers to full metadata including statistical validation (`cer_factor`) and cluster homogeneity (`err_factor`).

> **Note**: Speconsense's *core* output (`clusters/{sample}-all.fasta`) always includes `gid=`, `vid=`, `cer_factor=`, and `err_factor=` regardless of any summarize-side preset. The `--fasta-fields` option only controls headers written by `speconsense-summarize`.

## Why Customize Headers?

**Common use cases:**

1. **Downstream tool compatibility** — Some tools expect minimal headers or specific formats
2. **Readability** — Include only essential fields for cleaner output
3. **File size optimization** — Minimal headers reduce file size for large datasets
4. **Quality control** — Include `cer_factor`, `err_factor`, and read identity metrics for assessing consensus quality
5. **Workflow-specific needs** — Different analyses need different metadata

## Design Principles

The field customization system follows these principles:

1. **Field codes = field names** — No translation needed (e.g., `rawric` in command = `rawric=` in output)
2. **Preset composition** — Combine presets with union semantics: `--fasta-fields minimal,qc` expands to all fields from both
3. **Conditional output** — Fields with null/zero/empty values are automatically omitted for clean headers
4. **Order preservation** — Fields appear in the order specified (left to right)

## Quick Start

### Using Presets

```bash
# Default preset
speconsense-summarize --fasta-fields default
# Output: >sample-1.v1 size=638 ric=638 rawric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC

# Minimal headers - just the essentials
speconsense-summarize --fasta-fields minimal
# Output: >sample-1.v1 size=638 ric=638

# QC preset - includes statistical validation, homogeneity, length, and read identity
speconsense-summarize --fasta-fields qc
# Output: >sample-1.v1 size=638 ric=638 length=589 rid=98.5 cer_factor=12.4 err_factor=1.05

# Full metadata (most fields)
speconsense-summarize --fasta-fields full
# Output: >sample-1.v1 size=638 ric=638 length=589 rawric=333+305 snp=1 rid=98.5 cer_factor=12.4 err_factor=1.05 primers=...

# ID only (no metadata fields)
speconsense-summarize --fasta-fields id-only
# Output: >sample-1.v1
```

### Custom Field Selection

```bash
# Specify individual fields (comma-separated)
speconsense-summarize --fasta-fields size,ric,primers
# Output: >sample-1.v1 size=638 ric=638 primers=5'-ITS1F,3'-ITS4_RC

# Just the statistical validation fields
speconsense-summarize --fasta-fields cer_factor,err_factor
# Output: >sample-1.v1 cer_factor=12.4 err_factor=1.05
```

### Combining Presets and Fields

```bash
# Combine presets (union semantics)
speconsense-summarize --fasta-fields minimal,qc
# Output: >sample-1.v1 size=638 ric=638 length=589 rid=98.5 cer_factor=12.4 err_factor=1.05

# Preset + individual fields
speconsense-summarize --fasta-fields minimal,cer_factor,err_factor
# Output: >sample-1.v1 size=638 ric=638 cer_factor=12.4 err_factor=1.05
```

## Available Fields

| Field | Description | Example | Always Present? |
|-------|-------------|---------|-----------------|
| `size` | Total reads across merged variants | `size=638` | Yes |
| `ric` | Reads in Consensus | `ric=638` | Yes |
| `length` | Sequence length in bases | `length=589` | Yes |
| `rawric` | RiC values of `.raw` source variants (pre-merge) | `rawric=333+305` | Only when merged |
| `rawlen` | Original sequence lengths before overlap merging | `rawlen=630+361` | Only when overlap-merged |
| `snp` | Number of IUPAC ambiguity positions from merging | `snp=1` | Only when >0 |
| `ambig` | Total IUPAC ambiguity codes in consensus | `ambig=2` | Only when >0 |
| `rid` | Mean read identity (0–100%) | `rid=98.5` | When available |
| `rid_min` | Minimum read identity (0–100%) | `rid_min=94.1` | When available |
| `cer_factor` | Per-position CER factor (statistical validation) | `cer_factor=12.4` | When applicable (None for anchors) |
| `err_factor` | Cluster homogeneity ratio | `err_factor=1.05` | When applicable |
| `primers` | Detected primer names | `primers=ITS1F,ITS4` | When detected |
| `group` | Identity group number (extracted from filename) | `group=1` | Yes (superseded by `gid=` from core) |
| `variant` | Variant identifier within group | `variant=v1` | Yes (superseded by `vid=` from core) |

**Notes:**
- Conditional fields (marked "Only when…") are automatically omitted if not applicable
- `cer_factor` is `None` for anchors (largest cluster in a group) and clusters that fail to find a comparison peer; the field is dropped from the header in that case
- `group` and `variant` are vestigial — they parse the `-{gid}.v{vid}` pattern from the filename and emit redundant labels. Core's own `clusters/{sample}-all.fasta` already carries `gid=N` and `vid=N` directly, and summarize preserves `gid`/`vid` from input headers when re-writing
- See [Understanding RiC and Merging](understanding-ric-and-merging.md) for interpreting `rawric`

## Preset Reference

| Preset | Fields Included | Use Case |
|--------|-----------------|----------|
| `default` | `size, ric, rawric, rawlen, snp, ambig, primers` | General use; balanced metadata |
| `minimal` | `size, ric` | Quick inspection, essential info only |
| `qc` | `size, ric, length, rid, ambig, cer_factor, err_factor` | Quality control with statistical validation |
| `full` | `size, ric, length, rawric, rawlen, snp, ambig, rid, cer_factor, err_factor, primers` | Comprehensive metadata for archival or review |
| `id-only` | *(none)* | Clean IDs only, no metadata fields |

**Preset composition:**
- Presets can be combined: `minimal,qc` → union of both
- Mix presets and individual fields: `minimal,primers` → `size, ric, primers`
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

Include the statistical validation and homogeneity fields to assess consensus reliability:

```bash
# QC-focused headers
speconsense-summarize --fasta-fields qc
# Output: >sample-1.v1 size=638 ric=638 length=589 rid=98.5 cer_factor=12.4 err_factor=1.05

# Add CER + err_factor to default fields
speconsense-summarize --fasta-fields default,cer_factor,err_factor,rid
```

Sequences with high `err_factor` (≫ 1.0) or low `cer_factor` (close to 1.0) warrant review. The default summarize filters route extreme outliers to `__Summary__/variants/{name}.lq-` and `.ns-` respectively, but you can include the metrics in passing-track headers for inspection.

### For File Size Optimization

Large datasets with thousands of sequences benefit from minimal headers:

```bash
# Minimal headers reduce file size
speconsense-summarize --fasta-fields minimal

# Or just essential quality info
speconsense-summarize --fasta-fields size,ric,rid
```

### For Publication / GenBank Submission

Include comprehensive metadata for transparency:

```bash
# Full metadata for archival
speconsense-summarize --fasta-fields full

# Or default with added quality metrics
speconsense-summarize --fasta-fields default,rid,rid_min,cer_factor,err_factor
```

## Field Details

### Size vs RiC

- **`size`** — Total reads in the cluster(s) before any sampling
- **`ric`** — Reads actually used for consensus generation (may be sampled down via `--max-sample-size`)

For single (non-merged) variants: `size` ≥ `ric` (equality when cluster smaller than `--max-sample-size`).

For merged variants: `ric` can exceed `--max-sample-size` because it's the sum of RiC values from multiple pre-merge variants. See [Understanding RiC and Merging](understanding-ric-and-merging.md) for details.

### The rawric Field

Present only when variants were merged in summarize:

```
rawric=333+305
```

This shows:
- Two variants were merged to create this consensus
- The first variant had ric=333, the second had ric=305
- Values are ordered largest-first
- Total ric = 333+305 = 638

Use `rawric` to assess merge quality:
- **Balanced** (good): `rawric=60+55+50` — all variants have similar support
- **Dominated** (review): `rawric=100+8+5` — one strong variant, others weak

### The cer_factor Field

```
cer_factor=12.4
```

The per-position multiplicative inflation of the basecaller's error rate that would be required for the cluster's read counts at disagreement positions to be plausibly explained as noise replicates of a larger peer in the same identity group.

- **High `cer_factor` (≫ 1)**: the cluster's variants are statistically inconsistent with noise replication of any larger peer — likely a real biological variant
- **Low `cer_factor` (close to 1)**: the cluster's variants could be plausibly explained as basecaller noise; summarize routes these to `.ns-` files at `--min-cer-factor` (default 1.0)
- **Field omitted** (`cer_factor=None` internally): anchor (largest cluster in group) or no comparison peer found — always passes filtering

### The err_factor Field

```
err_factor=1.05
```

Peer-independent cluster homogeneity: observed disagreement at consensus columns divided by the expected disagreement under the basecaller's error model.

- **`err_factor ≈ 1.0`**: reads are consistent with basecaller noise around the consensus
- **`err_factor ≫ 1.0`**: cluster contains internal heterogeneity beyond what sequencing noise would produce — possibly a chimera, contamination, or two true haplotypes that escaped phasing

Routed to `.lq-` files when above `--max-err-factor` (default 1.5).

### Read Identity (rid, rid_min)

Computed from the cluster's MSA against the consensus, with homopolymer normalization:

- **`rid`** — Mean read identity (typical: 95–99% for good clusters)
- **`rid_min`** — Worst-case read identity

Low values can indicate outlier reads. Note that `err_factor` largely subsumes the role of these fields for outlier detection — it normalizes by the basecaller's expected error rate, whereas `rid` is a raw identity number that varies by basecaller and chemistry.

### Length Field

Simple sequence length in bases:
```
length=589
```

Useful for filtering by amplicon length, quality control, or size distribution analysis.

### Primers Field

Shows detected primers:
```
primers=ITS1F,ITS4
```

Primers are detected by Speconsense during consensus generation when a primers FASTA file is provided.

## Examples

### Example 1: Minimal for BLAST

```bash
speconsense-summarize --fasta-fields minimal --summary-dir blast_ready
```

Output:
```
>specimen-A-1.v1 size=638 ric=638
ATCGATCG...
>specimen-B-1.v1 size=203 ric=203
GCTAGCTA...
```

Clean, simple headers compatible with most sequence analysis tools.

### Example 2: QC Report

```bash
speconsense-summarize --fasta-fields qc --summary-dir qc_report
```

Output:
```
>specimen-A-1.v1 size=638 ric=638 length=589 rid=98.5 cer_factor=12.4 err_factor=1.05
ATCGATCG...
```

Includes all quality metrics for assessing consensus reliability.

### Example 3: Custom Combination

```bash
speconsense-summarize --fasta-fields size,ric,rawric,cer_factor,primers
```

Output for merged variant:
```
>specimen-A-1.v1 size=638 ric=638 rawric=333+305 cer_factor=15.2 primers=ITS1F,ITS4
ATCGATCG...
```

Output for anchor (no `cer_factor`):
```
>specimen-B-1.v1 size=203 ric=203 primers=ITS1F,ITS4
GCTAGCTA...
```

### Example 4: Preset Composition

```bash
speconsense-summarize --fasta-fields minimal,qc
```

Output (union of minimal + qc):
```
>specimen-A-1.v1 size=638 ric=638 length=589 rid=98.5 cer_factor=12.4 err_factor=1.05
ATCGATCG...
```

This combines essential metrics (size, ric) with quality control fields.

## Technical Notes

### Field Order

Fields appear in the order specified:

```bash
# Order matters
speconsense-summarize --fasta-fields primers,ric,size
# Output: >sample-1.v1 primers=... ric=638 size=638

# Different order
speconsense-summarize --fasta-fields size,ric,primers
# Output: >sample-1.v1 size=638 ric=638 primers=...
```

### Duplicate Removal

When combining presets, duplicates are automatically removed:

```bash
# minimal = size,ric
# qc = size, ric, length, rid, ambig, cer_factor, err_factor
speconsense-summarize --fasta-fields minimal,qc
# Result: size,ric,length,rid,ambig,cer_factor,err_factor (size,ric not duplicated)
```

### Conditional Fields

Some fields are only output when applicable:

- `rawric` — only for merged variants
- `rawlen` — only for overlap-merged variants
- `snp` — only when >0
- `ambig` — only when >0
- `rid`, `rid_min` — only when read identity data is available
- `cer_factor` — only when not None (i.e., not an anchor and a comparison peer was found)
- `err_factor` — only when computed
- `primers` — only when primers were detected
- `variant` — only for within-group variants (e.g., `.v2`, `.v3`); omitted for `.v1` anchors

This keeps headers clean and avoids clutter.

### Error Handling

Invalid field names produce clear error messages:

```bash
speconsense-summarize --fasta-fields size,invalid_field,ric
# Error: Unknown preset or field name: 'invalid_field'
```

## Backward Compatibility

The default behavior is unchanged in 0.8.x:

```bash
# These are equivalent
speconsense-summarize
speconsense-summarize --fasta-fields default
```

Both produce the same output format. **However**, the `default` preset does not include `cer_factor` / `err_factor` / `rid`. If you want those fields in summarize-emitted FASTAs, switch to `qc` or `full`, or compose explicitly: `--fasta-fields default,cer_factor,err_factor,rid`.

## Related Documentation

- **[Understanding RiC and Merging](understanding-ric-and-merging.md)** — Explains RiC accumulation, rawric interpretation, and variant merging philosophy
- **README.md** — See "Customizing FASTA Header Fields" section for basic usage
- **README.md** — See "FASTA Header Metadata" and "Statistical Variant Validation (CER)" sections for field semantics

## Summary

The `--fasta-fields` option provides flexible control over FASTA header metadata:

- **Presets** for common use cases (minimal, qc, full, etc.)
- **Custom field selection** with field codes matching field names
- **Preset composition** to combine multiple presets
- **Conditional output** to keep headers clean
- **0.8.x adds `cer_factor` and `err_factor`** for statistical validation and homogeneity inspection
