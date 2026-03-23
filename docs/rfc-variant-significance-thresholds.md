# RFC: Statistical Variant Significance Thresholds

**Status:** Draft — seeking feedback from sequence validators
**Author:** Josh Walker
**Date:** 2026-03-22
**Affects:** speconsense (core), speconsense-summarize

---

## 1. Problem

Speconsense performs recursive MSA-based variant phasing: it examines the partial order alignment within each cluster, identifies positions where reads disagree, and splits the cluster into haplotypes. The split decision currently uses two fixed thresholds applied locally within each subcluster:

- `--min-variant-frequency 0.10` (at least 10% of reads at that position)
- `--min-variant-count 5` (at least 5 reads supporting the alternative allele)

These thresholds are applied to the *current subcluster size*, not the original specimen read count. During recursive phasing, subclusters shrink at each level. A subcluster of 50 reads (split from an original 1000) can itself split when 5 reads (10%) show an alternative base — even though 5 out of 1000 total reads is well within the expected range of sequencing error.

This creates a risk of **error-driven splits**: the tool reports sequence variants that are statistically indistinguishable from nanopore substitution errors. These false variants consume analyst time during validation and, in automated pipelines, can inflate diversity estimates.

## 2. Solution

We add a statistical significance floor based on the **critical error rate** (p\*) framework described in the companion paper (`docs/variant_significance/`). Before accepting a variant split, speconsense now asks:

> "What per-position error rate would make this observation plausible as artifact?"

If the answer is below the assumed platform error rate, the variant is suppressed.

### How it works

At each candidate split position, the tool computes p\* by solving:

```
L * P(Binom(N, p*/3) >= M) = alpha
```

where:

| Symbol | Meaning | Source |
|--------|---------|--------|
| **N** | Total specimen input reads | Fixed at specimen level (not subcluster) |
| **M** | Reads supporting the minority allele at this position | From MSA at current recursion level |
| **L** | Consensus sequence length | From SPOA result at current recursion level |
| **alpha** | Significance level | User parameter (default 1e-5) |
| **p\*/3** | Per-base error probability under uniform model | Assumes errors distribute equally across 3 alternative bases |

The variant is **significant** (allowed to split) when p\* >= the assumed error rate. In other words: explaining this many reads as error would require an error rate at least as high as what we assume the platform actually produces. If p\* is below the assumed rate, the pattern is consistent with normal sequencing error and the split is suppressed.

### Key design choices

**N is total specimen reads, not subcluster size.** This is the conservative choice. A variant that appears in 5 out of 50 reads in a subcluster is really 5 out of the original 1000 specimen reads — the other 950 reads were assigned to other haplotypes at earlier recursion levels. Using the full N prevents deep recursion from manufacturing significance through denominator shrinkage.

**Uniform error model (q = p/3).** We assume substitution errors distribute equally across the three alternative bases. This is less conservative than the worst-case model (q = p, all errors produce the same base) but more realistic for nanopore chemistry. The uniform model yields p\* values roughly 3x higher than the conservative model, providing comfortable margins for well-supported variants while still catching weak ones.

**Per-variant significance, no cross-variant correction.** Each candidate split position is tested independently. We do not apply additional Bonferroni correction across the set of candidate positions within a cluster, because the positions are not independent (they share the same reads) and the existing frequency/count filters already reduce the candidate set.

**Bonferroni over positions (L).** The correction factor L accounts for multiple testing across all amplicon positions. This is conservative — most positions will not have variant signal — but avoids assumptions about which positions are "testable."

## 3. User interface

### speconsense (core tool)

Two new CLI parameters in the Variant Phasing group:

```
--assumed-error-rate FLOAT   Assumed per-position error rate for variant
                             significance (uniform model). Set to 0 to disable.
                             (default: 0.02)

--significance-level FLOAT   Significance level (alpha) for variant significance
                             testing. (default: 1e-5)
```

Both are also available as profile keys (`assumed-error-rate`, `significance-level`).

### speconsense-summarize

Two new CLI parameters in the Filtering group:

```
--assumed-error-rate FLOAT   Assumed per-position error rate. Variants with
                             critical error rate (cer) below this value are
                             filtered as potentially artifactual. Set to 0 to
                             disable. (default: 0 = disabled)

--no-cer-filter              Disable all CER-based filtering, even when
                             --assumed-error-rate is set.
```

The summarize tool's `--assumed-error-rate` defaults to 0 (disabled) because it may process output from runs with different parameters. When enabled, it filters input variants whose `cer` header value falls below the specified rate. Variants without a `cer` header (e.g., from older speconsense output) are never filtered.

Both are also available as profile keys (`assumed-error-rate`, `no-cer-filter`). The `assumed-error-rate` key can be set once in a profile to apply to both tools.

### FASTA header annotation

When variant phasing produces a split, the output FASTA headers include two new fields:

```
>sample-c1 size=200 ric=100 rid=95.2 rid_min=90.1 cer=0.17 cer.a=1e-05
```

| Field | Meaning |
|-------|---------|
| `cer=0.17` | Critical error rate (p\*) — the per-position error rate that would make this cluster's read count plausible as artifact. Higher is better. |
| `cer.a=1e-05` | The alpha (significance level) used to compute this p\*. |

The `cer` value is computed at output time from the final cluster's read count (M), total specimen reads (N), and consensus length (L). Every cluster gets a `cer` annotation when `--assumed-error-rate` is set, regardless of how it was formed (MCL clustering, phasing, merging). The annotation answers a simple question: "could this many reads be sequencing error?"

When `--assumed-error-rate 0` is set (CER disabled), no `cer` fields appear.

### FASTA field presets (speconsense-summarize)

The `cer` and `cer_alpha` fields are available in the `--fasta-fields` system:

- Added to the **qc** preset: `size, ric, length, rid, ambig, cer, cer_alpha`
- Added to the **full** preset: `size, ric, length, rawric, rawlen, snp, ambig, rid, cer, cer_alpha, primers`
- Available individually: `--fasta-fields default,cer,cer_alpha`

## 4. Defaults and behavior

### Default: 2% assumed error rate, alpha = 1e-5

At the standard operating point (N=1000 presample, L~700 for ITS), this requires roughly **M >= 19 reads** for a variant to pass the significance test at 1% assumed error, or **M >= 25** at 2%. (See the paper's Table 2 for the full matrix.)

For comparison, the existing `--min-variant-count 5` threshold allows splits with as few as 5 reads. The CER gate will suppress most of these at high N, while still allowing them at low N where 5 reads represents stronger evidence.

### Interaction with existing thresholds

The CER test is applied **after** the existing frequency/count filters (cheap) and **before** the within-cluster error calculation (expensive). A variant must pass all three gates:

1. **Frequency gate:** Alternative allele frequency >= `--min-variant-frequency`
2. **Count gate:** Alternative allele count >= `--min-variant-count`
3. **CER gate:** p\* >= `--assumed-error-rate` (new)

The existing thresholds remain useful as fast pre-filters. The CER gate adds the statistical dimension that accounts for total specimen size.

### Disabling the feature

Setting `--assumed-error-rate 0` disables the CER gate entirely. All variant splits that pass the existing frequency/count thresholds will proceed as before. This is the recommended setting for users who want to preserve exact backward compatibility.

## 5. Practical implications for validators

### What changes in output?

With the default 2% assumed error rate:

- **Well-supported variants (M >> 20):** No change. These variants have p\* well above any realistic error rate and will continue to be reported.
- **Marginal variants (M ~ 5-20 from large specimens):** Some will be suppressed. These are the variants most likely to be error-driven splits. Validators currently spending time on these can expect fewer of them.
- **Small specimens (N < 100):** Minimal change. When N is small, even M=5 represents a substantial fraction, and p\* will typically exceed the assumed error rate.

### New header fields for triage

The `cer` field provides a direct quality signal for variant triage:

- **cer > 0.10:** Very robust. Would require >10% per-position error to explain as artifact.
- **cer = 0.02-0.10:** Moderate. Significant at the default threshold, but closer to the boundary.
- **cer < 0.02:** Below default threshold (would be suppressed by CER gate). If seen, it means the variant was produced with CER disabled or a lower assumed error rate.

The `cer` value is interpretable independently of the specific alpha used — it represents the error rate needed, not a p-value. As platform chemistry improves and error rates drop, the same `cer` values become more convincing without changing any parameters.

### Comparison with existing quality signals

| Signal | What it measures | Limitation |
|--------|-----------------|------------|
| `ric` | Number of reads used for consensus | Doesn't account for total specimen size |
| `rid` | Mean read-to-consensus identity | Measures cluster homogeneity, not variant evidence |
| `cer` | Statistical strength of variant evidence | Requires `--assumed-error-rate > 0` |

## 6. Examples

### Example 1: Strong variant passes

Specimen with N=1000 reads. Phasing detects 80 reads supporting an alternative allele at a position in a 700bp amplicon.

```
p* = compute_critical_error_rate(N=1000, M=80, L=700, alpha=1e-5)
# p* ≈ 0.13 (13%)
```

Since 13% >> 2% (assumed error rate), this variant passes. Header: `cer=0.13 cer.a=1e-05`.

### Example 2: Weak variant suppressed

Same specimen, N=1000. Deep recursion produces a subcluster where 6 reads show an alternative allele.

```
p* = compute_critical_error_rate(N=1000, M=6, L=700, alpha=1e-5)
# p* ≈ 0.003 (0.3%)
```

Since 0.3% < 2% (assumed error rate), this variant is suppressed. The split does not occur, and those 6 reads remain in their parent cluster.

### Example 3: Small specimen, same count passes

Specimen with N=20 reads. 6 reads support an alternative allele.

```
p* = compute_critical_error_rate(N=20, M=6, L=700, alpha=1e-5)
# p* ≈ 0.07 (7%)
```

Since 7% > 2%, this variant passes despite having the same M=6 as Example 2. At N=20, 6 reads represent 30% of the specimen, which is much stronger evidence than 6/1000.

## 7. Migration and backward compatibility

- **Output format:** The `cer=` and `cer.a=` header fields are additive. Downstream tools that don't recognize them will ignore them (they appear after existing fields in the space-delimited header).
- **Default behavior change:** The default `--assumed-error-rate 0.02` means the CER gate is active by default. Some variants that were previously reported will now be suppressed. To preserve exact previous behavior, set `--assumed-error-rate 0`.
- **speconsense-summarize:** CER filtering is disabled by default (`--assumed-error-rate 0`). Existing summarize workflows are unaffected unless explicitly opted in.
- **Profiles:** The `assumed-error-rate` key works in both `speconsense` and `speconsense-summarize` sections of a profile, allowing a single profile to configure both tools consistently.
- **ConsensusInfo type:** Two new optional fields (`cer`, `cer_alpha`) with `None` defaults. Existing code using ConsensusInfo will not break.

## 8. Open questions for reviewers

1. **Is 2% the right default assumed error rate?** This is based on typical ONT per-position substitution rates for R10 chemistry. Should we default higher (more permissive, fewer suppressed variants) or lower (more aggressive filtering)?

2. **Should the default alpha be 1e-5 or something else?** The paper suggests 1e-5 as a "run-corrected" level (accounting for ~100 specimens per run). A per-specimen alpha of 0.05 would be more permissive.

3. **Should CER filtering in speconsense-summarize default to enabled?** Currently it defaults to 0 (disabled) to avoid surprising behavior when summarizing output from older runs. If most users will want it enabled, we could match the core tool's default of 0.02.

4. **Is the `cer` annotation useful for your validation workflows?** Would you use it for triage? Should it appear in additional presets (e.g., `default`)?

5. ~~**Edge case: MCL-only clusters.**~~ Resolved: all clusters now receive `cer` annotation computed from final cluster size, so MCL-only and phasing-split clusters are treated uniformly.

---

## Appendix: Reference values

### Minimum M for significance (uniform model, L=700)

From the companion paper, Table 2:

| N | p=1%, alpha=0.05 | p=1%, alpha=1e-5 | p=2%, alpha=1e-5 | p=3%, alpha=1e-5 |
|---|---|---|---|---|
| 100 | 5 | 8 | 11 | 11 |
| 500 | 9 | 14 | — | 22 |
| 1000 | 13 | 19 | — | 33 |

### Critical error rate p\* at common operating points (uniform model, L=700)

| N | M | alpha | p\* |
|---|---|-------|-----|
| 1000 | 100 | 0.05 | 20.25% |
| 1000 | 100 | 1e-5 | 16.59% |
| 100 | 10 | 0.05 | 6.55% |
| 100 | 10 | 1e-5 | 2.50% |
