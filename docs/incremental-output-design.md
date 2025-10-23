# Incremental Output Writing Design

**Date:** 2025-10-22
**Status:** Approved for implementation

## Overview

Modify `speconsense-summarize` to write individual FASTA and FASTQ files incrementally as each specimen is processed, while keeping summary/report files at the end.

## Motivation

**User Benefits:**
- Can start analyzing individual specimen files before all processing completes
- Partial results preserved if process crashes or is terminated
- Better visibility into progress by watching output directory

**Implementation Simplicity:**
- Code is already structured per-specimen via `process_single_specimen()`
- Naming is specimen-local (groups numbered 1, 2, 3 within each specimen)
- Dependencies can be satisfied by moving setup before processing loop

## File Writing Strategy

### Write Incrementally (per specimen, during processing loop)
- Individual FASTA files: `{sample}-{group}-RiC{ric}.fasta`
- Individual FASTQ files: `FASTQ Files/{sample}-{group}-RiC{ric}.fastq`
- Individual .raw FASTA files: `{sample}-{group}.raw{N}-RiC{ric}.fasta`
- Individual .raw FASTQ files: `FASTQ Files/{sample}-{group}.raw{N}-RiC{ric}.fastq`

### Write at End (after all processing completes)
- `summary.fasta` - combined index of all sequences
- `summary.txt` - statistics and totals
- `summarize_log.txt` - log file copy

**Rationale:** Summary files require complete information from all specimens to compute totals and provide a complete index. Individual data files are self-contained per specimen.

## Implementation Plan

### Step 1: Create per-specimen file writer function

**New function:** `write_specimen_data_files()`

**Location:** Insert after `write_consensus_fastq()` (around line 1110)

**Extracted from:** `write_output_files()` lines 1184-1288

**Signature:**
```python
def write_specimen_data_files(specimen_consensus: List[ConsensusInfo],
                               merge_traceability: Dict[str, List[str]],
                               naming_info: Dict,
                               summary_folder: str,
                               fastq_dir: str,
                               fastq_lookup: Dict[str, List[str]],
                               original_consensus_lookup: Dict[str, ConsensusInfo]
                               ) -> List[Tuple[ConsensusInfo, str]]:
    """
    Write individual FASTA and FASTQ files for a single specimen.
    Does NOT write summary files (summary.fasta, summary.txt).

    Returns:
        List of (raw_consensus, original_cluster_name) tuples for summary.fasta
    """
```

**Contents:**
1. Generate .raw file consensuses (current lines 1184-1228)
2. Write individual FASTA files for final consensus (current lines 1233-1248)
3. Write FASTQ files for final consensus (current lines 1250-1253)
4. Write .raw FASTA files (current lines 1255-1264)
5. Write .raw FASTQ files (current lines 1266-1288)
6. Return raw_file_consensuses for later use in summary.fasta

### Step 2: Refactor write_output_files()

**New signature:**
```python
def write_output_files(final_consensus: List[ConsensusInfo],
                      merge_traceability: Dict[str, List[str]],
                      naming_info: Dict,
                      all_raw_consensuses: List[Tuple[ConsensusInfo, str]],
                      summary_folder: str,
                      temp_log_file: str = None):
    """Write summary files only. Individual files already written per-specimen."""
```

**Changes:**
- Remove individual file writing code (extracted to write_specimen_data_files)
- Add `all_raw_consensuses` parameter (collected from all specimens)
- Remove `fastq_lookup` and `consensus_list` parameters (no longer needed)
- Keep summary.fasta writing (lines 1290-1315)
- Keep summary.txt writing (lines 1317-1344)
- Keep log file copying (lines 1346-1356)

### Step 3: Modify main() function

**Move before processing loop:**
```python
# Create output directories
os.makedirs(args.summary_dir, exist_ok=True)
os.makedirs(os.path.join(args.summary_dir, 'FASTQ Files'), exist_ok=True)

# Build lookup tables once
fastq_lookup = build_fastq_lookup_table(args.source)
original_consensus_lookup = {cons.sample_name: cons for cons in consensus_list}
```

**Add to accumulation in processing loop:**
```python
all_raw_consensuses = []  # New: collect .raw files from all specimens

for file_path in tqdm(...):
    # Process specimen (unchanged)
    final_consensus, merge_traceability, naming_info, limited_count = process_single_specimen(...)

    # NEW: Write individual files immediately
    specimen_raw_consensuses = write_specimen_data_files(
        final_consensus, merge_traceability, naming_info,
        args.summary_dir,
        os.path.join(args.summary_dir, 'FASTQ Files'),
        fastq_lookup,
        original_consensus_lookup
    )

    # Accumulate for summary files
    all_final_consensus.extend(final_consensus)
    all_merge_traceability.update(merge_traceability)
    all_raw_consensuses.extend(specimen_raw_consensuses)  # NEW
    total_limited_merges += limited_count

    # Update naming info (unchanged)
    ...
```

**Update final write call:**
```python
# Write summary files at end (after all processing)
write_output_files(
    all_final_consensus,
    all_merge_traceability,
    all_naming_info,
    all_raw_consensuses,  # NEW parameter
    args.summary_dir,
    temp_log_file.name
)
```

## Dependencies

All dependencies can be satisfied before the processing loop:

| Dependency | Source | Timing |
|------------|--------|--------|
| `summary_folder` | args.summary_dir | Available immediately |
| `fastq_dir` | Derived from summary_folder | Created before loop |
| `fastq_lookup` | Built from args.source | Build before loop |
| `original_consensus_lookup` | Built from consensus_list | Build before loop |
| `specimen_consensus` | From process_single_specimen() | Per-specimen in loop |
| `merge_traceability` | From process_single_specimen() | Per-specimen in loop |
| `naming_info` | From process_single_specimen() | Per-specimen in loop |

## Testing Strategy

1. **Functional correctness:** Run on test dataset and compare summary files (summary.fasta, summary.txt) before and after change - should be identical
2. **Incremental behavior:** Monitor output directory during processing to verify files appear progressively
3. **Interruption handling:** Kill process mid-run and verify partial results are present and valid
4. **Edge cases:** Test with single specimen, empty specimens, merged vs non-merged variants

## Estimated Effort

| Task | Time |
|------|------|
| Extract write_specimen_data_files() | 45 min |
| Refactor write_output_files() signature | 15 min |
| Modify main() loop | 30 min |
| Testing | 1 hour |
| **Total** | **2.5 hours** |

## Risk Assessment

**Low Risk:**
- Code already structured for per-specimen processing
- Extraction refactoring, not algorithmic changes
- Dependencies straightforward
- Summary files unchanged (easy to verify correctness)

**Considerations:**
- If crash occurs mid-specimen write, that specimen's files may be incomplete
- Could add atomic writes (write to temp, rename on success) if needed
- No performance impact expected (same total I/O, just different timing)

## Rollback Plan

If issues arise, changes are isolated to:
1. New function `write_specimen_data_files()` - can be removed
2. Modified `write_output_files()` - revert to previous version
3. Modified `main()` - revert loop structure

Git history provides clean rollback path.
