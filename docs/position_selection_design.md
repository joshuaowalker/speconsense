# Variant Position Selection for Phasing - Design Document

**Date:** 2025-11-28
**Status:** Implemented

## Problem Statement

When variant phasing detects multiple variant positions, using all positions may not produce optimal results:
1. Combining all positions can fragment reads into too many haplotypes where none qualify
2. Even when all positions yield qualifying haplotypes, a subset may produce cleaner clusters
3. Some variant positions may be noise rather than true biological variation

**Goal:** Select the subset of variant positions that minimizes within-cluster error (produces the most homogeneous haplotype clusters).

**When it runs:** Position selection always runs when there are ≥2 variant positions. This ensures we find the optimal position subset regardless of whether all positions would technically work.

## Design Principles

1. **Minimize within-cluster error** - Choose positions such that resulting haplotype clusters have minimal internal variation (unphased variation at other positions)
2. **Use established metrics** - Within-cluster variance minimization is foundational to clustering
3. **Use NMI as heuristic** - Normalized mutual information efficiently identifies correlated positions to narrow search space

## Algorithm Overview

### Two-Stage Approach

**Stage 1: Candidate Generation (NMI-based, fast)**
- Compute pairwise NMI for all variant positions
- Use greedy forward selection to identify promising position combinations
- Generate a small set of candidate position sets to evaluate

**Stage 2: Candidate Evaluation (within-cluster error, accurate)**
- For each candidate position set:
  - Group reads into haplotypes by allele combination at those positions
  - Calculate within-cluster error (variation at non-phased positions)
  - Only consider position sets that produce ≥2 qualifying haplotypes
- Select the position set with minimum within-cluster error

### Within-Cluster Error Calculation

The error calculation simulates the full haplotype assignment logic, including reassignment
of non-qualifying haplotypes to the nearest qualifying haplotype. This gives an accurate
picture of what the final error rate would be after phasing completes.

```
For each candidate position set S:
    Group reads by allele combination at positions in S → haplotypes H₁, H₂, ...

    If fewer than 2 qualifying haplotypes: skip this candidate

    # Simulate reassignment (mirrors phase_reads_by_variants logic)
    For each non-qualifying haplotype:
        Find nearest qualifying haplotype by edit distance on combo string
        Reassign reads to that qualifying haplotype

    # Calculate error on reassigned groups
    For each (reassigned) haplotype Hᵢ:
        For each MSA position p NOT in S:
            error_p = fraction of reads in Hᵢ with non-consensus allele at p
        mean_error_Hᵢ = average(error_p for all p not in S)

    overall_error = weighted_average(mean_error_Hᵢ, weights=|Hᵢ|)

Select position set with lowest overall_error
```

### Candidate Generation Strategy

```
If num_variant_positions <= 3:
    Evaluate all non-empty subsets directly
Else:
    candidates = []

    # 1. Greedy NMI selection (current algorithm)
    candidates.append(greedy_nmi_forward_selection())

    # 2. Each single position (fallback for 2-position case)
    for pos in variant_positions:
        candidates.append({pos})

    # 3. Best pairwise combinations by NMI
    for top_pair in top_k_nmi_pairs(k=3):
        candidates.append(top_pair)
```

## Key Metrics

### Normalized Mutual Information (NMI)

Measures correlation between allele distributions at two positions:
- NMI = 1.0: Perfectly correlated (knowing one predicts the other)
- NMI = 0.0: Independent (no predictive relationship)

```python
def calculate_normalized_mi(alleles1, alleles2):
    mi = mutual_information(alleles1, alleles2)
    return mi / min(entropy(alleles1), entropy(alleles2))
```

### Within-Cluster Error

For a haplotype cluster, the average error rate at positions not used for phasing:
- Lower error = more homogeneous cluster
- Directly measures what we care about: cluster quality

## Implementation Location

- `core.py:561-648` - Entropy and MI calculation functions
- `core.py:651-714` - `calculate_within_cluster_error()` function
- `core.py:717-936` - `select_correlated_positions()` with reassignment simulation
- `core.py:2357-2400` - Integration in `phase_reads_by_variants()`

## Implementation Notes

All changes have been implemented:

1. ✓ Added `calculate_within_cluster_error()` function
2. ✓ Modified `select_correlated_positions()` to:
   - Generate multiple candidates (not just one greedy result)
   - Simulate haplotype reassignment before calculating error
   - Evaluate each candidate by within-cluster error
   - Return the best candidate
3. ✓ Handle the 2-position case (select single best position if needed)
4. ✓ Error calculation includes reassignment simulation (mirrors `phase_reads_by_variants` logic)

## Example

```
Cluster with 4 variant positions: [100, 200, 300, 400]
Positions 100, 200 are correlated; 300, 400 are noise

Candidates evaluated:
  {100, 200}     → 2 haplotypes, error=0.03  ← BEST
  {100}          → 2 haplotypes, error=0.05
  {200}          → 2 haplotypes, error=0.06
  {100, 200, 300} → 1 qualifying haplotype, skip
  {300, 400}     → 0 qualifying haplotypes, skip

Selected: positions {100, 200}
```

## References

- Inverse Simpson Index: https://en.wikipedia.org/wiki/Diversity_index
- Within-cluster variance: standard clustering quality metric
- NMI for categorical variables: sklearn.metrics.normalized_mutual_info_score
