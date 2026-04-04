# Exploration Notes: IUPAC Ambiguity Codes and CER-Gated Phasing

**Date:** 2026-04-04
**Context:** Investigating when speconsense core produces IUPAC ambiguity codes, and whether the K>1 CER phasing gate can reduce unnecessary ambiguities.

## 1. When does speconsense core produce IUPAC codes?

### Pipeline flow

1. MCL/greedy clustering → initial clusters
2. Pre-phasing merge (HP-equivalent consensus)
3. Outlier removal
4. Variant detection (`min_variant_frequency`, `min_variant_count`)
5. Recursive phasing — CER gate decides whether to split
6. Post-phasing merge (HP-equivalent consensus)
7. Size filtering (`min_size`)
8. Final consensus generation (re-runs SPOA on each cluster's reads)
9. **Ambiguity calling** — scans final MSA for variant positions using `min_ambiguity_frequency` (default 0.10) and `min_ambiguity_count` (default 3)

IUPAC codes appear in step 9 when variant positions survive in the final cluster MSA.

### Three mechanisms for IUPAC codes

**Mechanism 1: CER suppression (most common at default thresholds)**

A variant is detected during phasing but CER determines the minor allele count (M) is not statistically distinguishable from sequencing error at the assumed error rate. The variant stays in one cluster and ambiguity calling annotates it. This is the CER framework working as designed — borderline variants that can't be confidently phased get annotated instead.

Empirical signature: minor allele at 11-15% frequency, p* just below `assumed_error_rate`.

Example: ONT08.75 cluster 2 (209 reads), two positions with 11-12% minor alleles, p*=0.0135 and 0.0148 (both below 0.015 threshold). The two minor alleles are **uncorrelated** — independent error/variant positions, not a two-population mixture.

**Mechanism 2: Post-phasing merge combining heterogeneous populations**

Two clusters with HP-equivalent consensus sequences get merged after phasing, but they differ at non-HP positions. These differences become IUPAC codes in the merged consensus.

Empirical signature: the IUPAC-bearing cluster has reads from two distinct origins, with a small correlated minority (6-10 reads carrying the same set of alleles at multiple positions).

Example: ONT06.61 cluster 3 (76 reads) = post-phasing merge of:
- 61 reads from a phasing split of initial cluster 1
- 15 reads from initial cluster 3

The 61-read portion contained 6 reads with a correlated 4-position minority haplotype (T,G,G,T). These 6 reads were at 1.6% in the original 372-read cluster (invisible to variant detection at 10% threshold), 9.8% in the 61-read split half (just below 10%), and 7.9% in the merged 76-read cluster. But uncorrelated sequencing errors at each position add 3-6 additional minor allele counts, pushing per-position totals to 9-12/76 (12-16%) — above ambiguity calling thresholds.

Key insight: a genuine 6-read sub-population was never large enough relative to its container to cross the 10% frequency threshold at any pipeline stage. It only becomes visible as IUPAC because noise inflates per-position counts.

**Mechanism 3: Threshold gap between phasing and ambiguity calling**

Default thresholds differ:
- Phasing: `min_variant_count=5`, `min_variant_frequency=0.10`
- Ambiguity: `min_ambiguity_count=3`, `min_ambiguity_frequency=0.10`

A variant with 3-4 reads won't trigger phasing but will trigger ambiguity calling. This is by design (annotate what can't be confidently split), but the gap means some positions are never even *attempted* by phasing.

### Denominator-shrinking amplification effect

When phasing splits a large cluster into smaller sub-clusters, the same absolute minor allele count represents a larger fraction. A 3-read minority that was 3/755 (0.4%) in the parent becomes 3/25 (12%) in a sub-cluster, passing the ambiguity frequency threshold. More aggressive phasing creates more small clusters with higher IUPAC rates.

### Quantitative summary (ONT98 dataset, 955 specimens)

| Run | Clusters | IUPAC positions | IUPAC in clusters <50 |
|-----|----------|-----------------|----------------------|
| baseline-mcr0 (no CER) | 3,792 | 936 | 98% |
| cer0.015-kg1-mcr0 (default CER) | 3,491 | 1,367 | 91% |
| cer0.015-kg1-mcr0-mvc3 (aggressive) | 4,331 | 1,747 | 93% |

CER reduces clusters (blocks some splits) but increases IUPAC (unsplit variants get annotated). Aggressive phasing thresholds increase both clusters AND IUPAC via the denominator-shrinking effect.

## 2. K>1 correlation heuristic — why it failed

The original K>1 implementation used a correlation threshold (default 0.9 = 90% overlap) to discover co-varying positions. This was too restrictive.

### ONT08.31 case study

Cluster with 83 reads, 3 variant positions (pos 90, 698, 727). All failed K=1 CER individually. The correlation heuristic checked pairwise overlap of minority read sets:

| Seed → Other | Overlap | Ratio | 0.9 threshold |
|---|---|---|---|
| pos 90 → 698 | 3/8 | 0.38 | FAIL |
| pos 90 → 727 | 3/8 | 0.38 | FAIL |
| pos 698 → 90 | 3/5 | 0.60 | FAIL |
| pos 698 → 727 | 0/5 | 0.00 | FAIL |
| pos 727 → 90 | 3/8 | 0.38 | FAIL |
| pos 727 → 698 | 0/8 | 0.00 | FAIL |

All pairs rejected. The minority alleles are partially but not tightly correlated — multiple sub-haplotypes carry different subsets of variant alleles. This is genuine biological heterogeneity, not the clean two-population structure the correlation heuristic assumed.

## 3. Unified CER phasing — what we built

Replaced the two-phase algorithm (K=1 evaluation + K>1 correlation fallback) with a unified beam search that evaluates all (K, M) combinations.

### Algorithm

1. **K=1 level**: Evaluate all F candidate positions individually. Keep top B=20 by within-cluster error.
2. **K=2..5 levels**: For each beam entry, expand by trying remaining positions. Evaluate CER at K=|positions|, M=min qualifying haplotype size. Keep top B by M (highest M = best split potential).
3. **Selection**: Among all CER-significant subsets found at any level, return the one with lowest within-cluster error.

### Key design decisions

- **CER is the gate, within-cluster error is the selector.** We initially tried selecting by highest p* (CER value), but this drove pathological behavior: the search kept adding positions to maximize p* even when the split was already significant, producing K=200+ subsets with M=3.

- **Beam width B=20, max K=5.** These are sufficient for the empirical data. K=5 with M=5 gives p*≈0.05 at N=400, L=750 — already very strong evidence.

- **Candidate position cap at 50 for K>1 expansion.** Positions are sorted by minority set size (descending), so the top 50 include the most promising ones.

### Results on ONT08.31

| Metric | Old (K=1 + correlation) | New (unified beam) |
|--------|------------------------|-------------------|
| Clusters | 11 | 13 |
| c1 size / IUPAC | 89 / 3 | 77 / 1 |
| c4 K>1 split | None (0.55-0.67 overlap rejected) | K=2, M=5, p*=0.023 |
| Total IUPAC | 18 | 16 |

The K=3 split of cluster 1 (M=26, p*=0.66) was particularly notable — the old algorithm only found a K=1 split at a single position while missing the much stronger 3-position pattern.

### Parameters removed

- `--k-correlation-threshold` — no longer needed; direct subset evaluation replaces the heuristic

### Parameters retained

- `--min-k-position-gap` — proximity filter for nearby positions (correlated error profiles)

## 4. Pathological parameter case: F=755

With `--min-variant-count 3 --min-variant-frequency 0.0001`, nearly every MSA position with 3+ non-consensus reads becomes a variant candidate. A 1000-read cluster produced F=755 candidate positions.

This is inherently expensive regardless of algorithm — evaluating 755 positions at K=1 requires 755 group-by-position + filter + error calculations, each on 1000 reads. The beam search handles K>1 gracefully (O(K × B × min(F, 50) × R) per cluster), but the K=1 level is unavoidably O(F × R).

The F=755 case takes minutes (dominated by the K=1 evaluation of all 755 positions), compared to milliseconds for F=1-5 with default thresholds. This is not a bug — with a 0.01% frequency threshold, the user is asking the algorithm to consider essentially every noisy position as a potential variant.

## 5. Open questions for future work

1. **Should ambiguity calling use the CER framework?** Currently ambiguity calling has no statistical gate — it calls any position meeting count/frequency thresholds. Adding a CER test (is this minor allele distinguishable from error?) could reduce false IUPAC codes. The challenge: ambiguity calling runs on the *final* consensus MSA, which may have different read composition than the phasing MSA.

2. **Post-phasing merge and IUPAC.** Mechanism 2 (merge combining heterogeneous populations) is a pipeline-level issue. The merge correctly identifies HP-equivalent consensus sequences, but the merged cluster can have variant positions that didn't exist in either source. Options: run a variant check before merging, or re-phase after merge.

3. **Denominator-shrinking amplification.** The ambiguity calling frequency threshold is applied to the *cluster* read count, not the *specimen* read count. A position-level test using N (total specimen reads) rather than cluster size would be more robust against fragmentation artifacts.

4. **Selection criterion.** We settled on within-cluster error as the selector (CER as gate). This could be revisited — p* captures statistical confidence while error captures split cleanliness. A combined metric might be better.

5. **Beam width and depth.** B=20, K_max=5 were chosen pragmatically. The beam width affects whether the search finds non-greedy solutions (where the best K=2 pair doesn't contain the best K=1 position). The depth affects how many correlated positions can be jointly evaluated. Both could be tuned or made adaptive.
