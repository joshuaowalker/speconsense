# RFC: Statistical Variant Significance Thresholds

**Status:** Draft — seeking feedback from sequence validators
**Author:** Josh Walker
**Date:** 2026-04-03
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
C(L, K) * P(Binom(N, (p*/3)^K) >= M) = alpha
```

where:

| Symbol | Meaning | Source |
|--------|---------|--------|
| **N** | Total specimen input reads | Fixed at specimen level (not subcluster) |
| **M** | Reads supporting the minority allele(s) | From MSA (core) or variant read count (summarize) |
| **L** | Consensus sequence length | From SPOA result (core) or sequence length (summarize) |
| **K** | Number of variant positions evaluated jointly | 1 for single-position, >1 for correlated multi-position variants |
| **alpha** | Significance level | User parameter (default 1e-5) |
| **p\*/3** | Per-base error probability under uniform model | Assumes errors distribute equally across 3 alternative bases |

At K=1 this reduces to `L * P(Binom(N, p*/3) >= M) = alpha`. At K>1, the joint probability `(p/3)^K` models K independent errors in the same reads, and the Bonferroni correction becomes C(L, K) — the number of ways to choose K positions from L.

The variant is **significant** (allowed to split) when p\* >= the assumed error rate. In other words: explaining this many reads as error would require an error rate at least as high as what we assume the platform actually produces. If p\* is below the assumed rate, the pattern is consistent with normal sequencing error and the split is suppressed.

### Key design choices

**N is total specimen reads, not subcluster size.** This is the conservative choice. A variant that appears in 5 out of 50 reads in a subcluster is really 5 out of the original 1000 specimen reads — the other 950 reads were assigned to other haplotypes at earlier recursion levels. Using the full N prevents deep recursion from manufacturing significance through denominator shrinkage.

**Uniform error model (q = p/3).** We assume substitution errors distribute equally across the three alternative bases. This is less conservative than the worst-case model (q = p, all errors produce the same base) but more realistic for nanopore chemistry. The uniform model yields p\* values roughly 3x higher than the conservative model, providing comfortable margins for well-supported variants while still catching weak ones.

**Per-variant significance, no cross-variant correction.** Each candidate split position is tested independently at K=1. We do not apply additional Bonferroni correction across the set of candidate positions within a cluster, because the positions are not independent (they share the same reads) and the existing frequency/count filters already reduce the candidate set.

**Bonferroni over positions.** The correction factor C(L, K) accounts for multiple testing across all ways to choose K positions from L amplicon positions. At K=1, this is simply L. At K=2, it is L(L−1)/2 — still conservative, since most position pairs will not have correlated variant signal.

**K>1 multi-position significance.** The K=1 test evaluates each position independently, which is conservative: a variant that differs at multiple positions is dramatically more significant under the joint error model because independent errors must co-occur in the same reads. The K>1 extension exploits this by jointly evaluating correlated positions. At K=2, the minimum M for significance drops from 23 to 6 at production settings (N=1000, L=700, p=1.5%). This is handled differently in the two stages:

- **Core (phasing gate):** K>1 is a fallback when no single position passes K=1 CER. The tool discovers correlated positions by finding those whose minority alleles co-vary with a seed position. A proximity filter (configurable, default 10 consensus positions) excludes nearby positions that may share correlated error profiles within the same k-mer context. Only positions with high read overlap (configurable, default 90%) are considered correlated.

- **Summarize (pairwise CER):** K is computed directly as the number of differences (substitutions + indel events) between two consensus sequences. Each candidate variant is evaluated pairwise against all validated variants in its HAC group, and its CER is the minimum p\* across all pairwise comparisons.

## 3. User interface

### speconsense (core tool)

New CLI parameters in the Variant Phasing group:

```
--assumed-error-rate FLOAT   Assumed per-position error rate for variant
                             significance (uniform model). Set to 0 to disable.
                             (default: 0.015)

--significance-level FLOAT   Significance level (alpha) for variant significance
                             testing. (default: 1e-5)

--min-k-position-gap INT     Minimum gap in consensus positions between correlated
                             variant positions for K>1 CER evaluation. (default: 10)

--k-correlation-threshold FLOAT
                             Minimum fraction of minority reads that must overlap
                             for positions to be considered correlated in K>1 CER.
                             (default: 0.9)
```

All are also available as profile keys (`assumed-error-rate`, `significance-level`, `min-k-position-gap`, `k-correlation-threshold`).

### speconsense-summarize

New CLI parameters in the Filtering group:

```
--assumed-error-rate FLOAT   Assumed per-position error rate. Secondary variants
                             with critical error rate (cer) below this value are
                             filtered as potentially artifactual. Primary variants
                             are never filtered. Set to 0 to disable.
                             (default: 0.015)

--no-cer-filter              Disable all CER-based filtering, even when
                             --assumed-error-rate is set.
```

CER filtering in speconsense-summarize uses **pairwise evaluation**. After HAC grouping, each group's variants are processed in descending size order. The primary (largest) variant is always validated. Each subsequent candidate is compared pairwise against all validated variants: K is computed as the number of sequence differences (substitutions + indel events, where each contiguous indel counts as 1), and p\* is computed at that K. The candidate's CER is the minimum p\* across all pairwise comparisons — i.e., the error rate needed to explain the candidate as an artifact of its most similar validated neighbor.

This pairwise approach naturally incorporates K>1: two consensus sequences differing at 3 positions have K=3, dramatically reducing the p\* threshold needed for validation. Single-base variants (K=1) are held to the standard, and multi-position variants receive appropriate credit for the improbability of correlated errors.

The primary variant in each group is never filtered — the artifact hypothesis requires a stronger competing signal, and the primary has none. This prevents eliminating single-cluster specimens or dominant variants where the CER question ("artifact of what?") is not coherent. Variants without a `cer` header (e.g., from older speconsense output) are never filtered.

N for the pairwise computation is loaded from the specimen's metadata JSON (`total_input_reads`), falling back to the sum of variant read counts if metadata is unavailable.

All parameters are also available as profile keys (`assumed-error-rate`, `no-cer-filter`). The `assumed-error-rate` key can be set once in a profile to apply to both tools.

### FASTA header annotation

When variant phasing produces a split, the output FASTA headers include two new fields:

```
>sample-c1 size=200 ric=100 rid=95.2 rid_min=90.1 cer=0.17 cer.a=1e-05
```

| Field | Meaning |
|-------|---------|
| `cer=0.17` | Critical error rate (p\*) — the per-position error rate that would make this cluster's read count plausible as artifact. Higher is better. |
| `cer.a=1e-05` | The alpha (significance level) used to compute this p\*. |

In the core tool, the `cer` value is computed at output time from the final cluster's read count (M), total specimen reads (N), and consensus length (L) at K=1. In speconsense-summarize, the `cer` value is updated to reflect the pairwise minimum p\* (which may use K>1). In both cases, the annotation answers a simple question: "could this many reads be sequencing error?"

CER is only reported when the value is physically meaningful — specifically, when p\* < 0.75. Under the uniform error model, p\* = 0.75 is the **signal destruction threshold**: each of the 4 bases becomes equally likely (25% each), and the implied reference population (reads carrying the true sequence) equals the variant count. Above this, the "true sequence" has fewer supporting reads than the variant, making the artifact hypothesis incoherent. In practice, this means dominant clusters (typically >33% of reads at N=1000) receive no `cer` annotation — they are the presumed true sequence, not artifact candidates.

When `--assumed-error-rate 0` is set (CER disabled), no `cer` fields appear.

### FASTA field presets (speconsense-summarize)

The `cer` and `cer_alpha` fields are available in the `--fasta-fields` system:

- Added to the **qc** preset: `size, ric, length, rid, ambig, cer, cer_alpha`
- Added to the **full** preset: `size, ric, length, rawric, rawlen, snp, ambig, rid, cer, cer_alpha, primers`
- Available individually: `--fasta-fields default,cer,cer_alpha`

## 4. Defaults and behavior

### Default: 1.5% assumed error rate, alpha = 1e-5

At the standard operating point (N=1000 presample, L~700 for ITS), this requires roughly **M >= 23 reads** at K=1 for a single-position variant to pass the significance test at the default 1.5% assumed error rate. Multi-position variants require dramatically fewer reads: M >= 6 at K=2, M >= 4 at K=3 (see Appendix). At lower N, the threshold is a larger fraction of reads (M=9 at N=100, M=6 at N=20 for K=1), ensuring CER is conservative for small specimens.

For comparison, the existing `--min-variant-count 5` threshold allows splits with as few as 5 reads. The CER gate will suppress most single-position variants at high N, while still allowing them at low N where 5 reads represents stronger evidence. Multi-position variants that fail K=1 significance can still pass at K>1, preventing the CER gate from blocking well-supported multi-site variants.

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
- **Marginal single-position variants (M ~ 5-20 from large specimens):** Some will be suppressed at K=1. These are the variants most likely to be error-driven splits. Validators currently spending time on these can expect fewer of them.
- **Multi-position variants:** Variants differing at multiple positions receive K>1 evaluation, which dramatically lowers the minimum M. A variant with 3 substitutions relative to its nearest neighbor needs only M=4 reads at N=1000 — making it nearly impossible for a real multi-site variant to be filtered.
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
# p* ≈ 0.12 (12%)
```

Since 12% >> 1.5% (assumed error rate), this variant passes. Header: `cer=0.12 cer.a=1e-05`.

### Example 2: Weak variant suppressed

Same specimen, N=1000. Deep recursion produces a subcluster where 6 reads show an alternative allele.

```
p* = compute_critical_error_rate(N=1000, M=6, L=700, alpha=1e-5)
# p* ≈ 0.0005 (0.05%)
```

Since 0.05% < 1.5% (assumed error rate), this variant is suppressed. The split does not occur, and those 6 reads remain in their parent cluster.

### Example 3: Multi-position variant passes at K>1

Same specimen, N=1000. During summarize, a candidate variant's consensus sequence differs from the primary at 3 positions (2 substitutions + 1 indel event). The candidate has 8 supporting reads.

```
p* = compute_critical_error_rate(N=1000, M=8, L=700, alpha=1e-5, K=3)
# p* ≈ 0.14 (14%)
```

At K=1, this variant would have p\* ≈ 0.001 (0.1%) — well below the 1.5% threshold. But at K=3, requiring 3 independent errors in the same reads, the implied error rate jumps to 14%. The variant passes comfortably.

### Example 4: Small specimen, same count passes

Specimen with N=20 reads. 6 reads support an alternative allele.

```
p* = compute_critical_error_rate(N=20, M=6, L=700, alpha=1e-5)
# p* ≈ 0.026 (2.6%)
```

Since 2.6% > 1.5%, this variant passes despite having the same M=6 as Example 2. At N=20, 6 reads represent 30% of the specimen, which is much stronger evidence than 6/1000.

## 7. Migration and backward compatibility

- **Output format:** The `cer=` and `cer.a=` header fields are additive. Downstream tools that don't recognize them will ignore them (they appear after existing fields in the space-delimited header).
- **Default behavior change:** The default `--assumed-error-rate 0.015` means the CER gate is active by default in both tools. Some variants that were previously reported will now be suppressed. To preserve exact previous behavior, set `--assumed-error-rate 0`.
- **speconsense-summarize:** CER filtering is now active by default (`--assumed-error-rate 0.015`), applied after HAC grouping to secondary variants only — primary variants are always retained. Use `--assumed-error-rate 0` or `--no-cer-filter` to disable.
- **Profiles:** The `assumed-error-rate` key works in both `speconsense` and `speconsense-summarize` sections of a profile, allowing a single profile to configure both tools consistently.
- **ConsensusInfo type:** Two new optional fields (`cer`, `cer_alpha`) with `None` defaults. Existing code using ConsensusInfo will not break.

## 8. Real-world validation (ONT98 dataset)

The CER framework was validated on the ONT98 dataset (960 wells across 10 ONT runs, ITS amplicons). Of these, 955 produced clustering output and 837 survived summarize filtering (min-RiC). Both speconsense and speconsense-summarize were run with the default settings (`--assumed-error-rate 0.015`, default `--min-cluster-ratio 0.01`).

### Two-stage impact

CER operates at two stages:

1. **Split prevention** (speconsense core): The CER gate suppresses weak splits during recursive phasing. At K=1, each candidate position is evaluated independently. When no position passes K=1, the K>1 fallback discovers correlated positions and evaluates their joint significance. Reads that would have formed a separate cluster instead remain in the parent cluster.
2. **Pairwise filtering** (speconsense-summarize): Post-HAC pairwise CER filtering removes secondary variants whose minimum pairwise p\* falls below the assumed error rate. K is computed as the number of sequence differences between each pair of consensus sequences. Primary variants are never filtered.

### Overall impact

| Metric | Baseline | CER K>=1 | Delta |
|--------|----------|----------|-------|
| Raw clusters (pre-summarize) | 3,369 | 3,033 | -336 (-10.0%) |
| Summary variants | 2,831 | 2,121 | -710 (-25.1%) |
| Merged (multi-component) variants | 497 | 321 | -176 |
| Single-variant specimens | 246 (29%) | 319 (38%) | +73 |
| HAC groups | 1,146 | 1,146 | 0 |
| Specimens | 837 | 837 | 0 |

Zero specimens and zero HAC groups were lost. The removed variants are overwhelmingly small secondary clusters from high-read-count specimens.

### Impact by specimen read count

The CER framework's impact scales with specimen size, consistent with its statistical design:

| Specimen Reads | Specimens | Baseline Variants | CER Variants | Delta | Reduction |
|----------------|-----------|-------------------|--------------|-------|-----------|
| < 50 | 84 | 95 | 94 | -1 | 1.1% |
| 50–99 | 60 | 105 | 99 | -6 | 5.7% |
| 100–199 | 102 | 217 | 196 | -21 | 9.7% |
| 200–499 | 181 | 548 | 450 | -98 | 17.9% |
| 500–999 | 196 | 808 | 578 | -230 | 28.5% |
| 1000+ | 214 | 1,058 | 704 | -354 | 33.5% |

**Low-read specimens (< 100 reads):** CER has minimal effect (-4.2% combined). The minimum M for significance at 1.5% error requires 30–50% of reads at small N (e.g., M=5 at N=10, M=6 at N=20), so most variants that survive the existing frequency/count filters are already significant.

**Mid-range specimens (100–499 reads):** Moderate impact. The minimum M threshold drops to 5–14% of reads, catching more marginal single-position variants.

**High-read specimens (500+ reads):** Largest reduction (-33% for 1000+ reads). Pairwise CER filtering at summarize removes secondary variants with K=1 (identical or nearly identical to a validated neighbor) that have insufficient read support. Multi-position variants (K>1) are retained at much lower M thresholds.

### Sequence quality

Cross-comparison with the structured verified dataset (1,145 organisms across 20 ONT runs):

| Metric | Value |
|--------|-------|
| Precision (all RiC thresholds) | 100% |
| Recall (RiC >= 1) | 100% |
| PRAUC | 0.9991 |
| Mean identity (count-weighted) | 100.00% |
| Variants >= 99.5% identity | 100% (2,119/2,119) |
| RiC drift | -3.6% |

No organisms were lost or misidentified. All 710 filtered variants were small secondary clusters; no false positives were introduced. The -3.6% RiC drift reflects reads from filtered variants being absorbed into parent clusters.

### IUPAC ambiguities

| Level | Baseline | CER K>=1 | Delta |
|-------|----------|----------|-------|
| Raw clusters (pre-summarize) | 928 | 1,520 | +592 (+64%) |
| Summary output | 1,615 | 1,794 | +179 (+11%) |

Raw cluster ambiguities increase because split prevention absorbs variant reads into parent clusters, encoding the variation as IUPAC codes rather than separate clusters. At the summary level, pairwise CER filtering removes some secondary variants that would have introduced ambiguities during merging, partially offsetting the increase from larger parent clusters.

### Primary variant stability

95% of primary (v1) sequences are identical between baseline and CER (795/837). Of the 42 that changed, most gained or lost a small number of IUPAC positions due to reads shifting between clusters.

### Edge cases observed

**Small specimens (N < 50):** Of 84 specimens in this range, CER affected only 1 variant. The framework is appropriately conservative at low N, deferring to the existing frequency/count filters.

**K>1 fallback at core level:** The K>1 phasing gate fallback fires rarely — it requires all candidate positions to fail K=1 significance AND correlate with each other. The primary K>1 impact is at the summarize level, where pairwise CER naturally incorporates multi-position differences between consensus sequences.

## 9. Open questions for reviewers

1. ~~**Is 2% the right default assumed error rate?**~~ Resolved: default set to 1.5%, slightly below the typical ONT R10 per-position substitution rate of ~2%. This provides a margin that avoids filtering variants near the boundary while still catching clear artifacts. Validated on ONT98: 25% variant reduction with zero specimen loss at 1.5%.

2. **Should the default alpha be 1e-5 or something else?** The paper suggests 1e-5 as a "run-corrected" level (accounting for ~100 specimens per run). A per-specimen alpha of 0.05 would be more permissive.

3. ~~**Should CER filtering in speconsense-summarize default to enabled?**~~ Resolved: yes, defaults to 0.015 to match the core tool. With primary-variant protection, the risk of surprising behavior is minimal — only secondary variants with weak statistical support are affected. Use `--no-cer-filter` to disable.

4. **Is the `cer` annotation useful for your validation workflows?** Would you use it for triage? Should it appear in additional presets (e.g., `default`)?

5. ~~**Edge case: MCL-only clusters.**~~ Resolved: all clusters now receive `cer` annotation computed from final cluster size, so MCL-only and phasing-split clusters are treated uniformly.

6. ~~**Should CER filtering protect primary variants?**~~ Resolved: CER filtering in speconsense-summarize now operates after HAC grouping and protects the primary (largest) variant in each group. Single-cluster specimens are never eliminated by CER filtering. This was validated on the ONT98 dataset — zero specimens lost with the updated filtering.

7. ~~**Multi-position variants (K>1).**~~ Resolved: K>1 support is now implemented at both stages. The core phasing gate uses correlation-based K discovery as a fallback. The summarize pipeline computes K as the number of sequence differences between consensus pairs. Validated on ONT98 with no organisms lost and 100% precision.

---

## Appendix: Reference values

### Minimum M for significance at default settings (uniform model, L=700, alpha=1e-5, p=1.5%)

| N | K=1 | K=2 | K=3 |
|---|-----|-----|-----|
| 10 | 5 (50%) | 3 (30%) | 3 (30%) |
| 20 | 6 (30%) | 3 (15%) | 3 (15%) |
| 50 | 7 (14%) | 4 (8%) | 3 (6%) |
| 100 | 9 (9%) | 4 (4%) | 3 (3%) |
| 200 | 11 (5.5%) | 4 (2%) | 3 (1.5%) |
| 500 | 16 (3.2%) | 5 (1%) | 3 (0.6%) |
| 1000 | 23 (2.3%) | 6 (0.6%) | 4 (0.4%) |

At K=1, the minimum M is a large fraction of reads at low N — CER is conservative and defers to frequency/count filters. At high N, the threshold drops to a few percent, catching weak single-position variants. At K>=2, the threshold drops dramatically because K independent errors must co-occur in the same reads, making multi-position variants far easier to validate.

### Minimum M at other operating points (uniform model, L=700, K=1)

From the companion paper, Table 2:

| N | p=1%, alpha=0.05 | p=1%, alpha=1e-5 | p=1.5%, alpha=1e-5 | p=2%, alpha=1e-5 | p=3%, alpha=1e-5 |
|---|---|---|---|---|---|
| 100 | 5 | 8 | 9 | 10 | 11 |
| 500 | 9 | 14 | 16 | 18 | 22 |
| 1000 | 13 | 19 | 23 | 26 | 33 |

### Critical error rate p\* at common operating points (uniform model, L=700)

| N | M | K | alpha | p\* |
|---|---|---|-------|-----|
| 1000 | 100 | 1 | 1e-5 | 16.59% |
| 1000 | 100 | 2 | 1e-5 | 66.63% |
| 1000 | 8 | 1 | 1e-5 | 0.12% |
| 1000 | 8 | 3 | 1e-5 | 13.79% |
| 100 | 10 | 1 | 1e-5 | 2.50% |
| 100 | 10 | 2 | 1e-5 | 20.13% |
