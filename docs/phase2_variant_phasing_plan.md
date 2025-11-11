# Phase 2: Automatic Variant Phasing - Implementation Plan

**Date:** 2025-11-09
**Status:** Planning
**Goal:** Automatically split clusters into separate haplotypes based on detected variant positions

---

## Design Decisions

### 1. Phasing Strategy
- **All unique allele combinations get separate clusters**
  - Once a position is identified as variant (≥20% alt allele frequency), ALL alleles at that position are phased separately
  - A cluster with 2 variant positions having alleles (C/T) and (A/G) creates up to 4 combinations: CA, CG, TA, TG
  - Each observed combination becomes a separate sub-cluster
  - No minimum threshold for haplotypes - rely on size filters to drop rare combinations

### 2. Processing Order: "Splitting First, Then Merging"
```
Initial Clustering → Outlier Removal → Variant Phasing (NEW)
                         ↓
            Homopolymer-Equivalent Merging
                         ↓
                  Size Filtering
```

### 3. Size Filtering
- **Applied AFTER phasing and merging** (not before)
- Ensures tiny haplotypes are correctly filtered
- `--min-size` and `--min-cluster-ratio` both apply to final phased results

### 4. Naming Scheme
- **Final output**: Renumber as c1, c2, c3... (phasing invisible in cluster IDs)
- **Provenance tracking**: Phasing history stored in metadata and debug files
- Clean, simple naming for users

### 5. Outlier Removal Threshold
- **Current behavior (INCORRECT)**: Only attempt removal if P95 > threshold
- **Corrected behavior**: Always check individual reads when `--mean-error-rate` is set
- Remove the `p95_var > self.outlier_threshold` gate condition

---

## Proposed Workflow

### Current Workflow (for reference)
```
1. Initial clustering (MCL/greedy)
2. Merge clusters with identical/homopolymer-equivalent consensus
3. Filter by size
4. For each cluster:
   - Generate consensus, calculate variance
   - (Optional) Remove outliers, regenerate
   - (Optional) Detect variants
   - Write output
```

### New Workflow (Phase 2)
```
1. Initial clustering (MCL/greedy)

2. For each initial cluster (VARIANT DETECTION PHASE):
   a. Sample reads (if cluster > max_sample_size)
   b. Generate consensus + MSA via SPOA
   c. Calculate variance metrics
   d. Remove outliers (if --mean-error-rate set and reads exceed threshold)
      - Regenerate consensus + MSA with filtered reads
      - Recalculate variance
   e. Detect variant positions (if --mean-error-rate set)
   f. **Phase reads into sub-clusters (NEW):**
      - If no variants detected: create single sub-cluster with all reads
      - If variants detected:
        * For each read, extract alleles at variant positions
        * Group reads by unique allele combinations
        * Create one sub-cluster per unique combination
        * Track provenance (original cluster ID + allele combo)

3. Collect all sub-clusters from all initial clusters

4. Merge sub-clusters (MERGE PHASE):
   a. For each sub-cluster, generate temporary consensus
   b. Group by identical/homopolymer-equivalent consensus
   c. Merge groups, combine read sets
   d. Track merge provenance

5. Filter merged results (FILTER PHASE):
   a. Filter by absolute size (--min-size)
   b. Filter by relative size (--min-cluster-ratio)
   c. Renumber as c1, c2, c3...
   d. Store provenance in metadata

6. For each final cluster (OUTPUT PHASE):
   a. Sample reads (if needed)
   b. Generate final consensus + MSA
   c. Calculate final variance metrics
   d. Trim primers
   e. Write all output files
```

---

## Implementation Tasks

### Task 1: Fix Outlier Removal Logic
**File:** `speconsense/core.py`
**Current Code (lines 711-735):**
```python
if (self.outlier_threshold is not None and consensus and msa and
    p95_var is not None and p95_var > self.outlier_threshold):
```

**Change:** Remove `p95_var > self.outlier_threshold` check
```python
if self.outlier_threshold is not None and consensus and msa:
```

**Rationale:** Always check individual reads against threshold when enabled, don't gate on P95.

---

### Task 2: MSA Optimization - Accept String MSAs
**File:** `speconsense/analyze.py`
**Goal:** Modify `extract_alignments_from_msa()` to accept either file path or string

**Current signature:**
```python
def extract_alignments_from_msa(msa_file: str) -> Tuple[List[ReadAlignment], str, Dict[int, int]]:
```

**New signature:**
```python
def extract_alignments_from_msa(msa_input: Union[str, Path]) -> Tuple[List[ReadAlignment], str, Dict[int, int]]:
```

**Implementation:**
- Check if `msa_input` is a string containing newlines (MSA content) or a path (file)
- If string: parse directly from string
- If path: read file and parse
- Update all callers in `core.py` to pass MSA strings instead of writing temp files

**Affected functions in core.py:**
- `identify_outlier_reads()` - line 887
- `calculate_read_variance()` - line 951
- `detect_variant_positions()` - line 1012

**Benefits:**
- Eliminates 3+ temporary file writes per cluster
- Faster processing
- Cleaner code

---

### Task 3: Implement Read Phasing Function
**File:** `speconsense/core.py`
**New method:** `phase_reads_by_variants()`

```python
def phase_reads_by_variants(
    self,
    msa_string: str,
    consensus_seq: str,
    cluster_read_ids: Set[str],
    variant_positions: List[Dict]
) -> List[Tuple[str, Set[str]]]:
    """Phase reads into haplotypes based on variant positions.

    Args:
        msa_string: MSA in FASTA format from SPOA
        consensus_seq: Ungapped consensus sequence
        cluster_read_ids: Set of read IDs in this cluster
        variant_positions: List of variant position dicts from detect_variant_positions()

    Returns:
        List of (allele_combo_string, read_id_set) tuples
        e.g., [("C-T-A", {id1, id2}), ("T-C-A", {id3, id4})]
    """
```

**Implementation details:**
1. Extract alignments from MSA
2. Map MSA positions to variant positions
3. For each read:
   - Extract alleles at each variant position
   - Create allele combination string (e.g., "C-T-A")
   - Handle gaps/missing data (treat as distinct allele)
4. Group reads by allele combination
5. Return list of (combo, read_set) tuples

**Edge cases:**
- Reads with gaps at variant positions: treat gap as a distinct allele
- Reads that don't cover all variant positions: use partial combination
- All combinations observed in the data are kept (no filtering here)

---

### Task 4: Restructure Main Clustering Pipeline
**File:** `speconsense/core.py`
**Method:** `cluster()`

**Major changes:**
1. **Remove early merging** (currently line 626)
   - Don't merge before per-cluster processing
   - Merging happens AFTER phasing

2. **Remove early size filtering** (currently line 632)
   - Size filtering happens AFTER phasing and merging

3. **Add phasing step** in per-cluster loop:
   ```python
   # After variant detection (currently line 747)
   if variant_positions:
       # Phase reads into sub-clusters
       phased_subclusters = self.phase_reads_by_variants(
           msa, consensus, cluster, variant_positions
       )
   else:
       # No variants - keep as single sub-cluster
       phased_subclusters = [(None, cluster)]

   # Store sub-clusters with provenance
   for combo, read_ids in phased_subclusters:
       all_subclusters.append({
           'read_ids': read_ids,
           'original_cluster_num': i,
           'allele_combo': combo,
           'variant_positions': variant_positions if combo else None
       })
   ```

4. **Add post-phasing merge phase:**
   ```python
   # After all clusters are phased
   merged_subclusters = self.merge_similar_clusters(all_subclusters)
   ```

5. **Add size filtering:**
   ```python
   # Filter by absolute size
   large_clusters = [c for c in merged_subclusters if len(c['read_ids']) >= self.min_size]

   # Filter by relative size
   if self.min_cluster_ratio > 0:
       largest_size = max(len(c['read_ids']) for c in large_clusters)
       large_clusters = [c for c in large_clusters
                         if len(c['read_ids']) / largest_size >= self.min_cluster_ratio]

   # Renumber as c1, c2, c3...
   ```

6. **Update output loop** to use final clusters with provenance tracking

---

### Task 5: Update merge_similar_clusters()
**File:** `speconsense/core.py`
**Method:** `merge_similar_clusters()`

**Current signature:**
```python
def merge_similar_clusters(self, clusters: List[Set[str]]) -> List[Set[str]]:
```

**New signature:**
```python
def merge_similar_clusters(self, clusters: List[Dict]) -> List[Dict]:
```

**Changes:**
- Accept list of cluster dictionaries (with provenance) instead of just read ID sets
- Preserve provenance metadata through merging
- Track which clusters were merged together
- Return merged clusters with combined provenance

**Cluster dict structure:**
```python
{
    'read_ids': Set[str],                    # Read IDs in this cluster
    'original_cluster_num': int,             # Initial cluster number
    'allele_combo': Optional[str],           # Allele combination (e.g., "C-T-A")
    'variant_positions': Optional[List[Dict]], # Variant metadata
    'merged_from': Optional[List[int]]       # Track merge history
}
```

---

### Task 6: Update write_cluster_files() and Metadata
**File:** `speconsense/core.py`

**Add provenance to cluster metadata:**
- Original cluster ID before phasing
- Allele combination used for phasing
- Which clusters were merged together
- Total number of haplotypes created from this initial cluster

**Debug file naming:**
- Keep current naming: `{sample}-c{N}-RiC{M}-*.fasta`
- Store provenance in metadata JSON and variant files

**Add to metadata JSON:**
```python
"phasing": {
    "enabled": bool,
    "total_initial_clusters": int,
    "total_phased_subclusters": int,
    "total_after_merging": int,
    "total_after_filtering": int
}
```

---

### Task 7: Update write_metadata()
**File:** `speconsense/core.py`

**Add phasing statistics:**
- How many clusters had variants detected
- How many sub-clusters created by phasing
- How many sub-clusters merged
- Final cluster count after filtering

---

### Task 8: Testing
**File:** `tests/test_phasing.py` (NEW)

**Test cases:**
1. **No variants detected** - cluster stays as single cluster
2. **Single variant position** - cluster splits into 2 haplotypes
3. **Two variant positions** - cluster splits into 4 haplotypes
4. **Rare allele combinations** - all combinations created, size filter drops small ones
5. **Merging after phasing** - two phased haplotypes from different clusters merge
6. **Outlier removal before phasing** - outliers removed, then phasing happens
7. **MSA optimization** - verify no temp files created during processing

**Update existing tests:**
- `test_summarize.py` - may need updates for new phasing behavior
- `test_augment_input.py` - verify augmented sequences work with phasing

---

## Output Changes

### Current Output
```
clusters/
  sample-all.fasta              # All consensus sequences
  cluster_debug/
    sample-c1-RiC500-*.fastq/fasta
    sample-c2-RiC250-*.fastq/fasta
```

### New Output (Phase 2)
```
clusters/
  sample-all.fasta              # All consensus sequences (renumbered c1, c2, c3...)
  cluster_debug/
    sample-c1-RiC300-*.fastq/fasta   # May come from initial cluster c1 haplotype 1
    sample-c2-RiC200-*.fastq/fasta   # May come from initial cluster c1 haplotype 2
    sample-c3-RiC150-*.fastq/fasta   # May come from merged haplotypes
    sample-metadata.json             # Includes phasing provenance
```

**Metadata example:**
```json
{
  "version": "0.6.0",
  "phasing": {
    "enabled": true,
    "initial_clusters": 5,
    "phased_subclusters": 12,
    "after_merging": 8,
    "after_filtering": 7
  },
  "clusters": {
    "c1": {
      "final_size": 300,
      "final_ric": 300,
      "provenance": {
        "source": "initial_c1_phased",
        "allele_combo": "C-T",
        "variant_positions": [45, 123]
      }
    },
    "c2": {
      "final_size": 200,
      "final_ric": 200,
      "provenance": {
        "source": "initial_c1_phased",
        "allele_combo": "T-C",
        "variant_positions": [45, 123]
      }
    },
    "c3": {
      "final_size": 150,
      "final_ric": 150,
      "provenance": {
        "source": "merged",
        "merged_from": ["initial_c2_phased_A-G", "initial_c3_phased_A-G"]
      }
    }
  }
}
```

---

## Future Optimizations (Not in Phase 2)

1. **Reduce variance metrics** to just mean and max (defer to later)
2. **Optimize SPOA calls** - currently generate consensus multiple times per cluster
3. **Parallel processing** - phase multiple clusters concurrently
4. **Smarter haplotype assignment** - use probabilistic models for ambiguous reads
5. **Reference-guided phasing** - use known haplotypes if available

---

## Risk Assessment

### Low Risk
- MSA optimization (Task 2) - straightforward refactoring
- Outlier removal fix (Task 1) - simple logic change
- Metadata updates (Tasks 6-7) - additive changes

### Medium Risk
- Phasing function (Task 3) - new algorithm, edge cases to handle
- Cluster restructuring (Task 4) - changes main pipeline flow

### High Risk
- Merge function changes (Task 5) - complex logic, easy to introduce bugs
- Complete workflow reorganization - could break existing behavior

### Mitigation
- Extensive testing with real data
- Add `--disable-phasing` flag for backward compatibility
- Keep existing merge behavior as fallback
- Comprehensive unit tests for phasing logic

---

## Open Questions

1. **Should we add a `--disable-phasing` flag** for users who want old behavior?
   - Pro: Easy rollback if phasing causes issues
   - Con: More code to maintain

2. **How to handle very complex haplotypes** (e.g., 3+ variant positions creating 8+ combinations)?
   - Currently: Create all combinations, let size filters handle it
   - Alternative: Add max haplotypes per cluster limit?

3. **Should variant detection be enabled by default** now that it drives phasing?
   - Currently: Only when `--mean-error-rate` is set
   - Alternative: Always detect variants, use default threshold?

4. **Performance impact** of generating consensus multiple times?
   - Detection phase: once per initial cluster
   - Merge phase: once per sub-cluster
   - Output phase: once per final cluster
   - Could be 2-3x more SPOA calls

---

## Success Criteria

1. ✅ Clusters with detected variants are split into separate haplotypes
2. ✅ Each unique allele combination creates a separate sub-cluster
3. ✅ Homopolymer-equivalent haplotypes merge correctly
4. ✅ Size filters correctly drop small haplotypes
5. ✅ Final output is cleanly numbered c1, c2, c3...
6. ✅ Provenance tracking shows phasing history
7. ✅ All existing tests still pass
8. ✅ New tests verify phasing behavior
9. ✅ No temporary files created during MSA processing
10. ✅ Outlier removal works correctly (no P95 gate)

---

## Timeline Estimate

- **Task 1** (Outlier fix): 15 min
- **Task 2** (MSA optimization): 1-2 hours
- **Task 3** (Phasing function): 3-4 hours
- **Task 4** (Pipeline restructure): 4-6 hours
- **Task 5** (Merge update): 2-3 hours
- **Task 6** (Metadata): 1 hour
- **Task 7** (write_metadata): 30 min
- **Task 8** (Testing): 3-4 hours

**Total:** ~15-20 hours

---

## Notes from User

- Order clustering as "splitting first, then merging"
- Remove p95 gate on outlier removal - always check reads when threshold set
- Variance metrics may eventually reduce to mean and max (future work)
- MSA optimization needed - avoid temp file writes
- For multiple variants: split ALL combinations, let filters handle dropping minor ones
- Once position is polymorphic (≥20% alt allele), ALL alleles at that position get phased
