# Speconsense Profile Parameters Reference

This document provides a complete reference of all configurable parameters across the bundled profiles. Use this to understand profile differences and customize settings for your workflow.

## Profile Descriptions

| Profile | Purpose |
|---------|---------|
| **default** | Balanced settings suitable for most fungal amplicon workflows |
| **compressed** | Compress variants into minimal IUPAC consensus sequences (aggressive merging with indels, 20% thresholds) |
| **herbarium** | High-recall for degraded DNA/type specimens where false negatives are costly |
| **largedata** | Experimental settings optimized for large input files |
| **nostalgia** | Simulate older bioinformatics pipelines (for comparison only) |
| **strict** | High-precision for confident results where false positives are costly |

---

## speconsense Parameters

Parameters for the main clustering and consensus tool.

| Parameter | Description | Default | compressed | herbarium | largedata | nostalgia | strict |
|-----------|-------------|---------|------------|-----------|-----------|-----------|--------|
| `algorithm` | Clustering algorithm: `graph` (MCL) or `greedy` | graph | — | — | — | greedy | — |
| `min-identity` | Minimum sequence identity threshold for clustering | 0.9 | — | 0.85 | — | 0.85 | 0.95 |
| `inflation` | MCL inflation parameter (higher = more clusters) | 4.0 | — | — | — | — | — |
| `k-nearest-neighbors` | Number of nearest neighbors for graph construction | 5 | — | — | — | — | — |
| `min-size` | Minimum cluster size to output (0 = disabled) | 5 | — | 3 | — | 5 | 10 |
| `min-cluster-ratio` | Minimum size relative to largest cluster (0 = disabled) | 0.01 | — | 0 | 0 | 0.2 | 0.05 |
| `max-sample-size` | Maximum reads sampled for consensus generation | 100 | — | 100 | — | 500 | 100 |
| `outlier-identity` | Min read-to-consensus identity (`auto` = (1+min-identity)/2) | auto | — | auto | — | 0.85 | 0.98 |
| `presample` | Presample size for initial reads (0 = use all) | 1000 | — | 0 | 0 | 500 | 0 |
| `disable-position-phasing` | Disable variant phasing within clusters | false | — | — | — | true | — |
| `min-variant-frequency` | Min minor allele frequency for variant phasing | 0.10 | 0.20 | 0.05 | — | — | 0.25 |
| `min-variant-count` | Min read count for minor allele to trigger phasing | 5 | — | — | — | — | — |
| `disable-ambiguity-calling` | Disable IUPAC codes for unphased variants | false | — | — | — | true | true |
| `min-ambiguity-frequency` | Min minor allele frequency for IUPAC calling | 0.10 | 0.20 | — | — | — | — |
| `min-ambiguity-count` | Min read count for IUPAC ambiguity calling | 3 | — | — | — | — | — |
| `disable-cluster-merging` | Disable merging identical consensus sequences | false | — | — | — | — | — |
| `disable-homopolymer-equivalence` | Require exact match for cluster merging | false | — | — | — | — | — |
| `orient-mode` | Sequence orientation: `skip`, `keep-all`, `filter-failed` | skip | — | — | — | — | — |
| `scale-threshold` | Sequence count to enable vsearch acceleration (0 = disabled) | 1001 | — | — | — | — | — |
| `threads` | Max threads for internal parallelism (0 = auto) | 1 | — | — | 0 | — | — |

**Notes:**
- "—" indicates the parameter uses the default value
- `outlier-identity: auto` calculates as (1 + min-identity) / 2, e.g., 0.95 for default min-identity of 0.9

---

## speconsense-summarize Parameters

Parameters for the post-processing and summarization tool.

| Parameter | Description | Default | compressed | herbarium | largedata | nostalgia | strict |
|-----------|-------------|---------|------------|-----------|-----------|-----------|--------|
| `min-ric` | Minimum Reads in Consensus threshold | 3 | — | 3 | — | 5 | 5 |
| `min-len` | Minimum sequence length in bp (0 = disabled) | 0 | — | — | — | — | — |
| `max-len` | Maximum sequence length in bp (0 = disabled) | 0 | — | — | — | — | — |
| `group-identity` | Identity threshold for HAC variant grouping | 0.9 | — | — | 0.95 | — | — |
| `disable-merging` | Skip MSA-based merge evaluation entirely | false | — | — | — | true | true |
| `merge-effort` | Merge thoroughness: `fast`, `balanced`, `thorough`, or 6-14 | balanced | — | — | fast | — | — |
| `merge-snp` | Enable SNP-based variant merging | true | — | — | — | — | — |
| `merge-indel-length` | Max individual indel length for merging (0 = disabled) | 0 | 5 | — | — | — | — |
| `merge-position-count` | Max total SNP+indel positions for merging | 2 | 10 | — | — | — | — |
| `merge-min-size-ratio` | Min size ratio (smaller/larger) for merging (0 = disabled) | 0.1 | 0.2 | — | — | — | — |
| `min-merge-overlap` | Min overlap in bp for different-length sequence merging | 200 | 0 | — | — | — | — |
| `disable-homopolymer-equivalence` | Treat homopolymer length differences as structural | false | — | — | — | — | — |
| `select-max-groups` | Max groups to output per specimen (-1 = all) | -1 | — | — | — | — | — |
| `select-max-variants` | Max variants per group (-1 = no limit) | -1 | — | — | — | — | — |
| `select-strategy` | Variant selection: `size` or `diversity` | size | — | — | — | — | — |
| `enable-full-consensus` | Generate full IUPAC consensus per variant group | false | true | — | — | — | — |
| `fasta-fields` | Header fields: preset or comma-separated list | default | — | — | — | — | — |
| `scale-threshold` | Sequence count to enable vsearch acceleration (0 = disabled) | 1001 | — | — | — | — | — |
| `threads` | Max threads for internal parallelism (0 = auto) | 0 | — | — | — | — | — |

**Notes:**
- "—" indicates the parameter uses the default value
- `fasta-fields` presets: `default`, `minimal`, `qc`, `full`, `id-only`

---

## Profile Selection Guide

| Scenario | Recommended Profile |
|----------|---------------------|
| Standard fungal amplicon workflow | **default** |
| Fewer sequences to review, variation captured as IUPAC codes | **compressed** |
| Herbarium/type specimens with degraded DNA | **herbarium** |
| Processing a single large FASTQ file | **largedata** |
| Comparing results to older pipeline outputs | **nostalgia** |
| Publication-quality sequences requiring high confidence | **strict** |

---

## Creating Custom Profiles

Profiles are YAML files stored in `~/.config/speconsense/profiles/`. Copy the example template created on first use:

```bash
cp ~/.config/speconsense/profiles/example.yaml ~/.config/speconsense/profiles/myprofile.yaml
```

Edit to customize parameters, then use with:

```bash
speconsense input.fastq -p myprofile
speconsense-summarize -p myprofile
```

CLI arguments always override profile values:

```bash
speconsense input.fastq -p herbarium --min-size 10  # Uses herbarium but overrides min-size
```

---

## Version Compatibility

Profiles specify a `speconsense-version` field (e.g., `"0.7.*"`) indicating compatible versions. Profiles will warn if used with incompatible versions.

---

*Document generated for speconsense 0.7.3*
