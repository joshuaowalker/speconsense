# RFC: Statistical Variant Significance Thresholds

**Status:** Draft — seeking feedback from sequence validators
**Author:** Josh Walker
**Date:** 2026-03-25
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
                             (default: 0.015)

--significance-level FLOAT   Significance level (alpha) for variant significance
                             testing. (default: 1e-5)
```

Both are also available as profile keys (`assumed-error-rate`, `significance-level`).

### speconsense-summarize

Two new CLI parameters in the Filtering group:

```
--assumed-error-rate FLOAT   Assumed per-position error rate. Secondary variants
                             with critical error rate (cer) below this value are
                             filtered as potentially artifactual. Primary variants
                             are never filtered. Set to 0 to disable.
                             (default: 0.015)

--no-cer-filter              Disable all CER-based filtering, even when
                             --assumed-error-rate is set.
```

CER filtering in speconsense-summarize is applied **after HAC grouping**, only to **secondary variants** within each group. The primary (largest) variant in each group is never filtered by CER — the artifact hypothesis requires a stronger competing signal, and the primary has none. This prevents eliminating single-cluster specimens or dominant variants where the CER question ("artifact of what?") is not coherent. Variants without a `cer` header (e.g., from older speconsense output) are never filtered.

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

The `cer` value is computed at output time from the final cluster's read count (M), total specimen reads (N), and consensus length (L). The annotation answers a simple question: "could this many reads be sequencing error?"

CER is only reported when the value is physically meaningful — specifically, when p\* < 0.75. Under the uniform error model, p\* = 0.75 is the **signal destruction threshold**: each of the 4 bases becomes equally likely (25% each), and the implied reference population (reads carrying the true sequence) equals the variant count. Above this, the "true sequence" has fewer supporting reads than the variant, making the artifact hypothesis incoherent. In practice, this means dominant clusters (typically >33% of reads at N=1000) receive no `cer` annotation — they are the presumed true sequence, not artifact candidates.

When `--assumed-error-rate 0` is set (CER disabled), no `cer` fields appear.

### FASTA field presets (speconsense-summarize)

The `cer` and `cer_alpha` fields are available in the `--fasta-fields` system:

- Added to the **qc** preset: `size, ric, length, rid, ambig, cer, cer_alpha`
- Added to the **full** preset: `size, ric, length, rawric, rawlen, snp, ambig, rid, cer, cer_alpha, primers`
- Available individually: `--fasta-fields default,cer,cer_alpha`

## 4. Defaults and behavior

### Default: 1.5% assumed error rate, alpha = 1e-5

At the standard operating point (N=1000 presample, L~700 for ITS), this requires roughly **M >= 19 reads** for a variant to pass the significance test at 1% assumed error, or **M >= 17** at 1.5%. (See the paper's Table 2 for the full matrix.)

For comparison, the existing `--min-variant-count 5` threshold allows splits with as few as 5 reads. The CER gate will suppress most of these at high N, while still allowing them at low N where 5 reads represents stronger evidence.

The default of 1.5% is slightly below the typical ONT R10 per-position substitution rate of ~2%, providing a margin that avoids filtering variants near the boundary while still catching clear artifacts.

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

With the default 1.5% assumed error rate:

- **Well-supported variants (M >> 20):** No change. These variants have p\* well above any realistic error rate and will continue to be reported.
- **Marginal variants (M ~ 5-20 from large specimens):** Some will be suppressed. These are the variants most likely to be error-driven splits. Validators currently spending time on these can expect fewer of them.
- **Small specimens (N < 100):** Minimal change. When N is small, even M=5 represents a substantial fraction, and p\* will typically exceed the assumed error rate.

### New header fields for triage

The `cer` field provides a direct quality signal for variant triage:

- **cer > 0.10:** Very robust. Would require >10% per-position error to explain as artifact.
- **cer = 0.02-0.10:** Moderate. Significant at the default threshold, but closer to the boundary.
- **cer < 0.015:** Below default threshold (would be suppressed by CER gate). If seen, it means the variant was produced with CER disabled or a lower assumed error rate.
- **No `cer` field:** The cluster is too large relative to total specimen reads for the artifact hypothesis to be meaningful (p\* >= 0.75). This is expected for dominant clusters and is not a quality concern.

The `cer` value depends on the alpha used (recorded in the `cer.a` field), so values computed at different significance levels are not directly comparable. However, unlike a p-value, `cer` has a concrete physical interpretation: it is the per-position error rate needed to explain the observation as artifact. As platform chemistry improves and error rates drop, the same `cer` values become more convincing without changing any parameters.

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

Since 13% >> 1.5% (assumed error rate), this variant passes. Header: `cer=0.13 cer.a=1e-05`.

### Example 2: Weak variant suppressed

Same specimen, N=1000. Deep recursion produces a subcluster where 6 reads show an alternative allele.

```
p* = compute_critical_error_rate(N=1000, M=6, L=700, alpha=1e-5)
# p* ≈ 0.003 (0.3%)
```

Since 0.3% < 1.5% (assumed error rate), this variant is suppressed. The split does not occur, and those 6 reads remain in their parent cluster.

### Example 3: Small specimen, same count passes

Specimen with N=20 reads. 6 reads support an alternative allele.

```
p* = compute_critical_error_rate(N=20, M=6, L=700, alpha=1e-5)
# p* ≈ 0.07 (7%)
```

Since 7% > 1.5%, this variant passes despite having the same M=6 as Example 2. At N=20, 6 reads represent 30% of the specimen, which is much stronger evidence than 6/1000.

## 7. Migration and backward compatibility

- **Output format:** The `cer=` and `cer.a=` header fields are additive. Downstream tools that don't recognize them will ignore them (they appear after existing fields in the space-delimited header).
- **Default behavior change:** The default `--assumed-error-rate 0.015` means the CER gate is active by default in both tools. Some variants that were previously reported will now be suppressed. To preserve exact previous behavior, set `--assumed-error-rate 0`.
- **speconsense-summarize:** CER filtering is now active by default (`--assumed-error-rate 0.015`), applied after HAC grouping to secondary variants only — primary variants are always retained. Use `--assumed-error-rate 0` or `--no-cer-filter` to disable.
- **Profiles:** The `assumed-error-rate` key works in both `speconsense` and `speconsense-summarize` sections of a profile, allowing a single profile to configure both tools consistently.
- **ConsensusInfo type:** Two new optional fields (`cer`, `cer_alpha`) with `None` defaults. Existing code using ConsensusInfo will not break.

## 8. Real-world validation (ONT98 dataset)

The CER framework was validated on the ONT98 dataset (955 specimens across 10 ONT runs, ITS amplicons). Both speconsense and speconsense-summarize were run with `--assumed-error-rate 0.02`.

### Impact on variant counts

**Speconsense (core):**

| Metric | Baseline | CER-aware |
|--------|----------|-----------|
| Total clusters | 3,792 | 3,424 (-9.7%) |
| Specimens affected | — | 200 reduced, 3 gained |

The modest reduction at the cluster level reflects CER gating during recursive phasing — weak splits (M ~ 5-15 from N ~ 1000) are suppressed while well-supported variants pass unchanged.

**Speconsense-summarize (with `--assumed-error-rate 0.02`):**

| Metric | Baseline | CER-aware |
|--------|----------|-----------|
| Output files | 3,186 | 1,602 (-50%) |
| Single-variant specimens | 226 (27%) | 452 (54%) |
| 1-2 variant specimens | 373 (45%) | 658 (79%) |
| 5+ variant specimens | 269 (32%) | 41 (5%) |
| Specimens lost | — | 0 |

The combined effect of CER gating (phasing) and CER filtering (summarize) halves the output while preserving all specimens. The shift toward single-variant output is consistent with the hypothesis that many baseline variants were error-driven splits.

### Precision and recall

Evaluated against a verified reference set of 1,159 organisms:

- **Precision: 100%** at all RIC thresholds — no false positives introduced
- **Recall: 100%** at RIC >= 1 — all organisms recovered
- **All 1,600 variant matches at >= 99.5% identity**
- **PRAUC: 0.9991**

Zero organisms were lost. The 1,584 removed variants were all secondary sequences that did not contribute unique organism detections.

### Quality metrics

**CER value distribution** (cluster-level, 2,632 annotated variants):

| CER range | Count | Interpretation |
|-----------|-------|----------------|
| < 0.01 | 1,603 (61%) | Low support, many below 2% threshold |
| 0.01 - 0.02 | 128 (5%) | Near threshold boundary |
| 0.02 - 0.10 | 355 (13%) | Moderate to good support |
| 0.10 - 0.75 | 546 (21%) | Strong support |
| No CER (dominant) | 792 | p\* >= 0.75, not applicable |

**Read identity:** Mean RID improved slightly (99.37% → 99.40%) as noisier low-support clusters were eliminated.

**IUPAC ambiguities:** Total ambiguous bases in cluster output increased from 936 to 1,561 (+67%). When CER prevents a split, the variation that would have produced a separate variant is instead captured as IUPAC ambiguity codes in the parent cluster. This is the expected tradeoff — positional ambiguity in a strong consensus vs. a separate weak variant — and is generally preferable for taxonomic identification.

### Edge cases observed

**Small specimens (N < 25):** The CER gate is effectively inert because few variants can reach the significance threshold. The existing frequency/count filters remain the gatekeepers. This is by design — at very low N, even the frequency-based filters provide adequate protection.

**Correlated multi-position variation:** One cluster (ONT01.19, size=11) accumulated 35 IUPAC ambiguities after CER prevented splitting reads that differed at multiple correlated positions. The current framework evaluates positions independently (K=1), so correlated variation is not modeled. This is acknowledged in the companion paper as out of scope for this version — multi-position variants are at least as significant as single-position variants, so the K=1 analysis is conservative.

## 9. Open questions for reviewers

1. ~~**Is 2% the right default assumed error rate?**~~ Resolved: default set to 1.5%, slightly below the typical ONT R10 per-position substitution rate of ~2%. This provides a margin that avoids filtering variants near the boundary while still catching clear artifacts.

2. **Should the default alpha be 1e-5 or something else?** The paper suggests 1e-5 as a "run-corrected" level (accounting for ~100 specimens per run). A per-specimen alpha of 0.05 would be more permissive.

3. ~~**Should CER filtering in speconsense-summarize default to enabled?**~~ Resolved: yes, defaults to 0.015 to match the core tool. With primary-variant protection, the risk of surprising behavior is minimal — only secondary variants with weak statistical support are affected. Use `--no-cer-filter` to disable.

4. **Is the `cer` annotation useful for your validation workflows?** Would you use it for triage? Should it appear in additional presets (e.g., `default`)?

5. ~~**Edge case: MCL-only clusters.**~~ Resolved: all clusters now receive `cer` annotation computed from final cluster size, so MCL-only and phasing-split clusters are treated uniformly.

6. ~~**Should CER filtering protect primary variants?**~~ Resolved: CER filtering in speconsense-summarize now operates after HAC grouping and protects the primary (largest) variant in each group. Single-cluster specimens are never eliminated by CER filtering. This was validated on the ONT98 dataset — zero specimens lost with the updated filtering.

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
