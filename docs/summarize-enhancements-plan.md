# speconsense-summarize Enhancement Implementation Plan

**Developer Guide for Merging and Variant Management Improvements**

*Last Updated: 2025-10-21*

## Executive Summary

This document outlines the implementation plan for enhancing `speconsense-summarize` with improved variant merging capabilities, better parameter organization, and MSA-based consensus generation.

**Key Changes:**
- Add indel merging support (currently SNP-only)
- Implement MSA-based merging using SPOA with subset evaluation
- Reorganize parameters by processing phase
- Add group output limiting and size-ratio filtering
- Support outputting raw pre-merge variants

**Architecture Change:**
- **Current:** Filter → Merge → HAC Group → Select → Output
- **Proposed:** Filter → HAC Group → MSA Merge (per group) → Select → Output

---

## Current Architecture

### Processing Pipeline (`process_single_specimen`, line 985)

```python
def process_single_specimen(file_consensuses, args):
    # Phase 1: Merge variants within specimen
    merged_consensus, merge_traceability = merge_variants_within_specimen(
        file_consensuses, args.snp_merge_limit
    )

    # Phase 2: HAC clustering to separate variant groups
    variant_groups = perform_hac_clustering(
        merged_consensus, args.variant_group_identity
    )

    # Phase 3: Select representative variants per group
    final_consensus = []
    for group in variant_groups:
        selected = select_variants(
            group, args.max_variants, args.variant_selection
        )
        final_consensus.extend(selected)

    return final_consensus, merge_traceability, naming_info
```

### Current Merging (`merge_variants_within_specimen`, line 1057)

**Algorithm:** Iterative pairwise greedy merging
- Build list of all valid pairwise merges
- Select merge with largest combined size
- Create IUPAC consensus from 2 sequences
- Repeat until no valid merges remain

**Limitations:**
- Only handles SNPs (no indels)
- Pairwise only (order-dependent for >2 sequences)
- Merges across entire specimen (doesn't respect HAC groups)

---

## Phase 1: Simple Improvements

**Goal:** Quick wins with minimal code changes
**Effort:** 2-4 hours
**Lines Changed:** ~50

### 1.1 Change `--max-variants` Default

**File:** `speconsense/summarize.py`
**Location:** `parse_arguments()` (line 77)

```python
# Current
parser.add_argument("--max-variants", type=int, default=2,
                    help="Maximum number of additional variants...")

# Change to
parser.add_argument("--max-variants", type=int, default=-1,
                    help="Maximum number of additional variants to output per group (default: -1 = no limit)")
```

**Rationale:** Current default of 2 is arbitrary and surprising. -1 (unlimited) is more transparent.

### 1.2 Add `--max-groups` Parameter

**File:** `speconsense/summarize.py`
**Location:** `parse_arguments()` (line 77)

```python
parser.add_argument("--max-groups", type=int, default=-1,
                    help="Maximum number of groups to output per specimen (default: -1 = all groups)")
```

**Implementation:** In `process_single_specimen()` after HAC clustering:

```python
# After perform_hac_clustering()
if args.max_groups > 0 and len(variant_groups) > args.max_groups:
    # Sort groups by size of largest member
    sorted_groups = sorted(
        variant_groups.items(),
        key=lambda x: max(m.size for m in x[1]),
        reverse=True
    )
    # Keep only top N groups
    variant_groups = dict(sorted_groups[:args.max_groups])
    logging.info(f"Filtered to top {args.max_groups} groups by size")
```

### 1.3 Add `--merge-min-size-ratio` Parameter

**File:** `speconsense/summarize.py`
**Location:** `parse_arguments()` (line 77)

```python
parser.add_argument("--merge-min-size-ratio", type=float, default=0.0,
                    help="Minimum size ratio (smaller/larger) for merging clusters (default: 0.0 = disabled)")
```

**Implementation:** In merge evaluation logic:

```python
# During pairwise merge evaluation
size_ratio = min(consensus_i.size, consensus_j.size) / max(consensus_i.size, consensus_j.size)
if size_ratio < args.merge_min_size_ratio:
    logging.debug(f"Skipping merge: size ratio {size_ratio:.3f} < {args.merge_min_size_ratio}")
    continue
```

**Use case:** Prevents merges like `size=100 + size=5` (ratio=0.05) with `--merge-min-size-ratio=0.1`

---

## Phase 2: Enhanced Merging with Indels

**Goal:** Support indel merging with reorganized parameters
**Effort:** 1 week
**Lines Changed:** ~250

### 2.1 New Parameter Structure

**Merge Phase Parameters:**

```python
parser.add_argument("--merge-snp", type=bool, default=True,
                    help="Enable SNP-based merging (default: True)")

parser.add_argument("--merge-indel-length", type=int, default=0,
                    help="Maximum length of individual indels allowed in merging (default: 0 = disabled)")

parser.add_argument("--merge-position-count", type=int, default=2,
                    help="Maximum total SNP+indel positions allowed in merging (default: 2)")

# Already added in Phase 1
# parser.add_argument("--merge-min-size-ratio", ...)
```

**Backward Compatibility:**

```python
# Support old parameter name
parser.add_argument("--snp-merge-limit", dest="merge_position_count", type=int,
                    help=argparse.SUPPRESS)  # Hidden but functional

# Check if old parameter was used
if '--snp-merge-limit' in sys.argv:
    logging.warning("--snp-merge-limit is deprecated, use --merge-position-count instead")
```

### 2.2 Updated Distance Calculation

**Current:** `calculate_substitution_distance()` (line 194)
Returns: `int` (SNP count or -1 if indels present)

**Proposed:** `calculate_variant_distance()`
Returns: `dict` with detailed breakdown

```python
def calculate_variant_distance(seq1: str, seq2: str) -> dict:
    """
    Calculate distance between sequences counting SNPs and indels separately.

    Returns dict with:
        'snp_count': int - number of substitution positions
        'indel_count': int - number of indel positions
        'max_indel_length': int - length of longest consecutive indel
        'compatible': bool - whether sequences can be aligned
    """
    if not seq1 or not seq2:
        return {'compatible': False}

    # Get alignment from edlib with IUPAC awareness
    result = edlib.align(seq1, seq2, task="path", additionalEqualities=IUPAC_EQUIV)
    if result["editDistance"] == -1:
        return {'compatible': False}

    alignment = edlib.getNiceAlignment(result, seq1, seq2)
    if not alignment:
        return {'compatible': False}

    # Analyze alignment using adjusted_identity for scoring
    score_result = score_alignment(
        alignment['query_aligned'],
        alignment['target_aligned'],
        adjustment_params=STANDARD_ADJUSTMENT_PARAMS,
        scoring_format=ScoringFormat(
            match='|',
            substitution='X',
            indel_start='I',
            indel_extension='-',
            homopolymer_extension='=',
            end_trimmed='.'
        )
    )

    # Count variant types
    snp_count = score_result.score_aligned.count('X')
    indel_starts = score_result.score_aligned.count('I')
    indel_extensions = score_result.score_aligned.count('-')

    # Calculate indel positions and max length
    indel_positions = []
    if indel_starts > 0 or indel_extensions > 0:
        # Parse score_aligned to find indel runs
        i = 0
        while i < len(score_result.score_aligned):
            if score_result.score_aligned[i] in ('I', '-'):
                # Start of indel
                start = i
                while i < len(score_result.score_aligned) and score_result.score_aligned[i] in ('I', '-'):
                    i += 1
                indel_positions.append((start, i))
            else:
                i += 1

    max_indel_length = max((end - start for start, end in indel_positions), default=0)

    return {
        'compatible': True,
        'snp_count': snp_count,
        'indel_count': len(indel_positions),
        'max_indel_length': max_indel_length
    }
```

### 2.3 Updated Merge Evaluation

**Current:** `merge_variants_within_specimen()` uses `calculate_substitution_distance()`

**Updated logic:**

```python
# Check compatibility
dist = calculate_variant_distance(consensus_i.sequence, consensus_j.sequence)

if not dist['compatible']:
    continue

# Check SNP limit
if not args.merge_snp or dist['snp_count'] > 0:
    if not args.merge_snp:
        logging.debug(f"Skipping merge: SNP merging disabled")
        continue

# Check indel limits
if dist['indel_count'] > 0:
    if args.merge_indel_length == 0:
        logging.debug(f"Skipping merge: indel merging disabled")
        continue
    if dist['max_indel_length'] > args.merge_indel_length:
        logging.debug(f"Skipping merge: max indel length {dist['max_indel_length']} > {args.merge_indel_length}")
        continue

# Check total position count
total_positions = dist['snp_count'] + dist['indel_count']
if total_positions > args.merge_position_count:
    logging.debug(f"Skipping merge: {total_positions} positions > {args.merge_position_count}")
    continue

# Check size ratio (from Phase 1)
size_ratio = min(consensus_i.size, consensus_j.size) / max(consensus_i.size, consensus_j.size)
if size_ratio < args.merge_min_size_ratio:
    continue

# Valid merge!
valid_merges.append({...})
```

### 2.4 Add `merged_ric` Field

**Implementation:** Track merge provenance during merging

```python
# In merge_variants_within_specimen()
merged_info = ConsensusInfo(
    sample_name=consensus_i.sample_name,
    cluster_id=consensus_i.cluster_id,
    sequence=merged_sequence,
    ric=consensus_i.ric + consensus_j.ric,
    size=consensus_i.size + consensus_j.size,
    file_path=consensus_i.file_path,
    snp_count=snp_count,
    primers=consensus_i.primers,
    merged_ric=[consensus_i.ric, consensus_j.ric]  # NEW: track RiC values
)

# When writing headers (in write_output_files(), line 904)
header_parts = [f"size={consensus.size}", f"ric={consensus.ric}"]

if hasattr(consensus, 'merged_ric') and consensus.merged_ric:
    # Sort largest-first
    ric_values = sorted(consensus.merged_ric, reverse=True)
    header_parts.append(f"merged_ric={'+'.join(str(r) for r in ric_values)}")

if consensus.snp_count is not None and consensus.snp_count > 0:
    header_parts.append(f"snp={consensus.snp_count}")
```

**Update `ConsensusInfo` NamedTuple:**

```python
class ConsensusInfo(NamedTuple):
    sample_name: str
    cluster_id: str
    sequence: str
    ric: int
    size: int
    file_path: str
    snp_count: Optional[int] = None
    primers: Optional[List[str]] = None
    merged_ric: Optional[List[int]] = None  # NEW: RiC values of merged clusters
```

---

## Phase 3: MSA-Based Merging

**Goal:** Implement order-independent merging using SPOA MSA
**Effort:** 2-3 weeks
**Lines Changed:** ~400

### 3.1 Architecture Change

**Key insight:** HAC grouping must happen BEFORE merging to separate dissimilar sequences.

**New pipeline:**

```python
def process_single_specimen(file_consensuses, args):
    # Phase 1: HAC clustering to separate variant groups (primary vs contaminants)
    variant_groups = perform_hac_clustering(
        file_consensuses, args.group_identity
    )

    # Phase 2: MSA-based merging within each group
    merged_groups = {}
    all_merge_traceability = {}

    for group_id, group_members in variant_groups.items():
        merged, traceability = merge_group_with_msa(
            group_members, args
        )
        merged_groups[group_id] = merged
        all_merge_traceability.update(traceability)

    # Phase 3: Select representative variants per group
    final_consensus = []
    naming_info = {}

    for group_id, merged_variants in merged_groups.items():
        selected = select_variants(
            merged_variants, args.max_variants, args.variant_selection
        )
        # ... naming and output

    return final_consensus, all_merge_traceability, naming_info
```

### 3.2 SPOA MSA Integration

**Reference:** `msalign.py` (line 76)

**Important difference from msalign.py:** The msalign.py tool is designed for visual comparison of sequences from different specimens with potentially different primers, so it treats terminal gaps as missing data (converts them to '.' and ignores them). However, for our merging implementation, all variants within a group share the same primers. Terminal gaps represent real biological variation (deletions at amplicon ends), so we treat them as variant positions just like internal gaps. We do NOT distinguish terminal gaps from internal gaps.

**Implementation:**

```python
def run_spoa_msa(sequences: List[str]) -> List[SeqRecord]:
    """
    Run SPOA to create multiple sequence alignment.

    Args:
        sequences: List of DNA sequence strings

    Returns:
        List of SeqRecord objects with aligned sequences (including gaps)
    """
    import tempfile
    import subprocess
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
        try:
            # Write sequences to temporary file
            records = [
                SeqRecord(Seq(seq), id=f"seq{i}", description="")
                for i, seq in enumerate(sequences)
            ]
            SeqIO.write(records, temp_input, "fasta")
            temp_input.flush()

            # Run SPOA with alignment output (-r 2)
            result = subprocess.run(
                ['spoa', temp_input.name, '-r', '2'],
                capture_output=True,
                text=True,
                check=True
            )

            # Parse aligned sequences from SPOA output
            aligned_sequences = []
            lines = result.stdout.strip().split('\n')
            current_id = None
            current_seq = []

            for line in lines:
                if line.startswith('>'):
                    if current_id is not None:
                        # Skip consensus sequence (usually last)
                        if not current_id.startswith('Consensus'):
                            aligned_sequences.append(SeqRecord(
                                Seq(''.join(current_seq)),
                                id=current_id,
                                description=""
                            ))
                    current_id = line[1:]
                    current_seq = []
                elif line.strip():
                    current_seq.append(line.strip())

            # Add last sequence (if not consensus)
            if current_id is not None and not current_id.startswith('Consensus'):
                aligned_sequences.append(SeqRecord(
                    Seq(''.join(current_seq)),
                    id=current_id,
                    description=""
                ))

            return aligned_sequences

        finally:
            if os.path.exists(temp_input.name):
                os.unlink(temp_input.name)
```

**Note:** SPOA returns aligned sequences with gaps represented as `-`. We use all gaps directly in our analysis - no need to distinguish terminal vs internal gaps since all variants share the same primers.

### 3.3 Subset Evaluation Algorithm

**Core algorithm:**

```python
def merge_group_with_msa(variants: List[ConsensusInfo], args) -> Tuple[List[ConsensusInfo], Dict]:
    """
    Find largest mergeable subset of variants using MSA-based evaluation.

    Algorithm:
    1. Run SPOA MSA on all variants
    2. Generate subsets in descending size order
    3. Evaluate each subset for compatibility (SNP/indel limits)
    4. Return first (largest) compatible subset as merged consensus

    Args:
        variants: List of ConsensusInfo from HAC group
        args: Command-line arguments with merge parameters

    Returns:
        (merged_variants, merge_traceability) where merged_variants is list
        of merged ConsensusInfo objects, and traceability maps merged names
        to original cluster names
    """
    if len(variants) == 1:
        return variants, {}

    logging.info(f"MSA-based merging of {len(variants)} variants")

    # Step 1: Run SPOA MSA on all variants
    sequences = [v.sequence for v in variants]
    aligned_seqs = run_spoa_msa(sequences)

    logging.debug(f"Generated MSA with length {len(aligned_seqs[0].seq)}")

    # Step 2: Generate subsets in descending size order
    subsets_by_size = generate_subsets_by_total_size(variants, args)

    logging.debug(f"Evaluating {len(subsets_by_size)} candidate subsets")

    # Step 3: Find first (largest) compatible subset
    for subset_indices in subsets_by_size:
        subset_variants = [variants[i] for i in subset_indices]
        subset_aligned = [aligned_seqs[i] for i in subset_indices]

        # Analyze MSA for this subset
        variant_stats = analyze_msa_columns(subset_aligned)

        # Check compatibility against merge limits
        if is_compatible_subset(variant_stats, args):
            logging.info(f"Found mergeable subset of {len(subset_indices)} variants: "
                        f"{variant_stats['snp_count']} SNPs, "
                        f"{variant_stats['indel_count']} indels")

            # Create merged consensus
            merged_consensus = create_consensus_from_msa(
                subset_aligned, subset_variants
            )

            # Track merge provenance
            traceability = {
                merged_consensus.sample_name: [v.sample_name for v in subset_variants]
            }

            # If --output-raw-variants, include originals as additional variants
            if args.output_raw_variants:
                # Sort originals by size
                raw_variants = sorted(subset_variants, key=lambda v: v.size, reverse=True)
                return [merged_consensus] + raw_variants, traceability

            # Otherwise, continue merging remaining variants
            remaining_indices = [i for i in range(len(variants)) if i not in subset_indices]
            if remaining_indices:
                remaining_variants = [variants[i] for i in remaining_indices]
                remaining_merged, remaining_trace = merge_group_with_msa(remaining_variants, args)
                traceability.update(remaining_trace)
                return [merged_consensus] + remaining_merged, traceability
            else:
                return [merged_consensus], traceability

    # No compatible subsets found - return largest variant alone
    logging.info(f"No mergeable subsets found, outputting {len(variants)} separate variants")
    return variants, {}
```

**Subset generation with optimizations:**

```python
def generate_subsets_by_total_size(variants: List[ConsensusInfo], args) -> List[Tuple[int, ...]]:
    """
    Generate subsets of variant indices in descending order by total cluster size.

    Optimizations:
    - Pre-filter by size ratio
    - Start with largest subsets (r=N, then r=N-1, etc.)
    - Early termination: first valid subset is optimal
    """
    import itertools

    n = len(variants)
    sizes = [v.size for v in variants]

    # Build list of (total_size, subset_indices) tuples
    candidates = []

    # Generate from largest subsets to smallest
    for r in range(n, 0, -1):
        for subset_indices in itertools.combinations(range(n), r):
            # Pre-filter by size ratio if enabled
            if args.merge_min_size_ratio > 0 and len(subset_indices) > 1:
                subset_sizes = [sizes[i] for i in subset_indices]
                min_size = min(subset_sizes)
                max_size = max(subset_sizes)
                size_ratio = min_size / max_size

                if size_ratio < args.merge_min_size_ratio:
                    continue  # Skip this subset

            # Calculate total size
            total_size = sum(sizes[i] for i in subset_indices)
            candidates.append((total_size, subset_indices))

    # Sort by total size descending
    candidates.sort(reverse=True, key=lambda x: x[0])

    # Return just the subset indices
    return [subset for _, subset in candidates]
```

**MSA column analysis:**

```python
def analyze_msa_columns(aligned_seqs: List[SeqRecord]) -> dict:
    """
    Analyze aligned sequences to count SNPs and indels.

    Important: All gaps (including terminal gaps) count as variant positions
    since variants within a group share the same primers.

    Returns dict with:
        'snp_count': number of positions with >1 non-gap base
        'indel_count': number of positions with gaps mixed with bases
        'max_indel_length': length of longest consecutive indel run
    """
    alignment_length = len(aligned_seqs[0].seq)

    snp_positions = []
    indel_positions = []

    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]

        # Separate bases from gaps (all gaps count, including terminal)
        unique_bases = set(c for c in column if c != '-')
        has_gap = '-' in column

        # Skip all-gap columns (shouldn't happen with SPOA output)
        if not unique_bases and has_gap:
            continue

        # SNP position: multiple different bases (ignoring gaps)
        if len(unique_bases) > 1:
            snp_positions.append(col_idx)

        # Indel position: mix of gaps and bases
        if has_gap and unique_bases:
            indel_positions.append(col_idx)

    # Calculate max indel length (consecutive indel positions)
    max_indel_length = 0
    if indel_positions:
        current_run = 1
        for i in range(1, len(indel_positions)):
            if indel_positions[i] == indel_positions[i-1] + 1:
                current_run += 1
                max_indel_length = max(max_indel_length, current_run)
            else:
                current_run = 1
        max_indel_length = max(max_indel_length, current_run)

    return {
        'snp_count': len(snp_positions),
        'indel_count': len(indel_positions),
        'max_indel_length': max_indel_length
    }
```

**Compatibility checking:**

```python
def is_compatible_subset(variant_stats: dict, args) -> bool:
    """Check if variant statistics are within merge limits."""

    # Check SNP limit
    if variant_stats['snp_count'] > 0 and not args.merge_snp:
        return False

    # Check indel limits
    if variant_stats['indel_count'] > 0:
        if args.merge_indel_length == 0:
            return False
        if variant_stats['max_indel_length'] > args.merge_indel_length:
            return False

    # Check total position count
    total_positions = variant_stats['snp_count'] + variant_stats['indel_count']
    if total_positions > args.merge_position_count:
        return False

    return True
```

### 3.4 Size-Weighted Consensus Generation

```python
def create_consensus_from_msa(aligned_seqs: List[SeqRecord],
                             variants: List[ConsensusInfo]) -> ConsensusInfo:
    """
    Generate consensus from MSA using size-weighted majority voting.

    At each position:
    - Weight each variant by cluster size
    - Choose majority representation (base vs gap)
    - For multiple bases, generate IUPAC code representing all variants

    Important: All gaps (including terminal) count as variant positions
    since variants share the same primers.

    Args:
        aligned_seqs: MSA sequences with gaps as '-'
        variants: Original ConsensusInfo objects (for size weighting)

    Returns:
        ConsensusInfo with merged consensus sequence
    """
    from collections import defaultdict

    consensus_seq = []
    snp_count = 0
    alignment_length = len(aligned_seqs[0].seq)

    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]

        # Weight each base/gap by cluster size
        votes_with_size = [(base, variants[i].size) for i, base in enumerate(column)]

        # Count size-weighted votes (EXACT match only, no IUPAC expansion)
        votes = defaultdict(int)
        for base, size in votes_with_size:
            votes[base.upper()] += size

        # Separate gap votes from base votes
        gap_votes = votes.get('-', 0)
        base_votes = {b: v for b, v in votes.items() if b != '-'}

        # Determine if position should be included
        total_base_votes = sum(base_votes.values())

        if total_base_votes > gap_votes:
            # Majority wants a base - include position
            if len(base_votes) == 1:
                # Single base - no ambiguity
                consensus_seq.append(list(base_votes.keys())[0])
            else:
                # Multiple bases - generate IUPAC code
                represented_bases = set(base_votes.keys())
                iupac_code = IUPAC_CODES.get(frozenset(represented_bases), 'N')
                consensus_seq.append(iupac_code)
                snp_count += 1
        # else: majority wants gap, omit position

    # Create merged ConsensusInfo
    consensus_sequence = ''.join(consensus_seq)
    total_size = sum(v.size for v in variants)
    total_ric = sum(v.ric for v in variants)
    merged_ric_values = sorted([v.ric for v in variants], reverse=True)

    # Use name from largest variant
    largest_variant = max(variants, key=lambda v: v.size)

    return ConsensusInfo(
        sample_name=largest_variant.sample_name,
        cluster_id=largest_variant.cluster_id,
        sequence=consensus_sequence,
        ric=total_ric,
        size=total_size,
        file_path=largest_variant.file_path,
        snp_count=snp_count if snp_count > 0 else None,
        primers=largest_variant.primers,
        merged_ric=merged_ric_values
    )
```

### 3.5 Output Raw Variants

**Parameter:**

```python
parser.add_argument("--output-raw-variants", action="store_true",
                    help="Output raw pre-merge sequences as additional .v* variants")
```

**Implementation:** In `merge_group_with_msa()` (see above), when flag is set:

```python
if args.output_raw_variants:
    # Return merged consensus + all raw variants
    raw_variants = sorted(subset_variants, key=lambda v: v.size, reverse=True)
    return [merged_consensus] + raw_variants, traceability
```

**Naming:** Use existing `.v*` convention, add metadata to headers

```python
# Merged consensus
>specimen-1 size=250 ric=250 snp=2 merged_ric=100+89+61

# Raw variants (if --output-raw-variants)
>specimen-1.v1 size=100 ric=100 merged_from=c1
>specimen-1.v2 size=89 ric=89 merged_from=c2
>specimen-1.v3 size=61 ric=61 merged_from=c3
```

---

## Parameter Migration Plan

### New Parameter Names (Phase 3)

**Organized by processing phase:**

```python
# Filter phase
--min-ric                    # (unchanged)

# Merge phase
--merge-snp                  # Enable SNP merging (new)
--merge-indel-length         # Max indel length (new)
--merge-position-count       # Max positions (replaces --snp-merge-limit)
--merge-min-size-ratio       # Min size ratio (new)
--merge-strategy             # iterative|msa (future)

# Group phase
--group-identity             # (replaces --variant-group-identity)

# Select phase
--select-max-variants        # (replaces --max-variants)
--select-max-groups          # Max groups (new)
--select-strategy            # (replaces --variant-selection)

# Output options
--output-raw-variants        # Output raw pre-merge variants (new)
```

### Deprecation Strategy

**Phase 2:** Add aliases, keep old names working

```python
parser.add_argument("--merge-position-count", "--snp-merge-limit",
                    dest="merge_position_count", type=int, default=2)

# Detect which was used
if '--snp-merge-limit' in sys.argv:
    logging.warning("--snp-merge-limit is deprecated, use --merge-position-count")
```

**Phase 3:** Full migration with comprehensive aliases

```python
# Group phase
parser.add_argument("--group-identity", "--variant-group-identity",
                    dest="group_identity", type=float, default=0.9)

# Select phase
parser.add_argument("--select-max-variants", "--max-variants",
                    dest="select_max_variants", type=int, default=-1)

parser.add_argument("--select-strategy", "--variant-selection",
                    dest="select_strategy", choices=["size", "diversity"], default="size")
```

**Future (v0.5.0+):** Remove old names, keep only new prefixed names

---

## Testing Strategy

### Phase 1 Tests

**Test cases:**
1. `--max-variants=-1` outputs all variants (new default)
2. `--max-groups=2` outputs only top 2 groups by size
3. `--merge-min-size-ratio=0.2` prevents dominated merges

**Test data:** Use existing test specimens with multiple variants

### Phase 2 Tests

**Test cases:**
1. Merge with 1-2bp indels (`--merge-indel-length=2`)
2. Combined SNP+indel merging (`--merge-position-count=3`)
3. `merged_ric` field appears in headers
4. Backward compatibility: `--snp-merge-limit` still works

**Test data:**
- Create synthetic variants with known indels
- Verify indel resolution uses size-weighted majority

### Phase 3 Tests

**Test cases:**
1. MSA-based merging produces same results as iterative (for SNP-only)
2. Multi-sequence merging (3+ variants) produces correct IUPAC consensus
3. Subset evaluation finds largest compatible group
4. `--output-raw-variants` outputs merged + raw sequences
5. Order-independent: different input orders produce same output

**Test data:**
- 5 variants with known SNP patterns
- Verify largest mergeable subset is selected
- Check consensus matches expected IUPAC codes

**Performance tests:**
- 10 variants: should complete in <100ms
- 15 variants: should complete in <1s
- 20 variants: acceptable if <30s

### Integration Tests

**End-to-end tests:**
1. Real ONT data with multiple specimens
2. Compare Phase 3 output to Phase 2 output (should be compatible)
3. Verify FASTQ file concatenation still works
4. Check all header fields are present and correct

---

## Implementation Checklist

### Phase 1 (2-4 hours)
- [ ] Change `--max-variants` default to -1
- [ ] Add `--max-groups` parameter
- [ ] Implement group filtering logic
- [ ] Add `--merge-min-size-ratio` parameter
- [ ] Implement size ratio checking
- [ ] Write unit tests
- [ ] Update help text
- [ ] Test on real data

### Phase 2 (1 week)
- [ ] Add new merge parameters (--merge-snp, --merge-indel-length, --merge-position-count)
- [ ] Implement `calculate_variant_distance()` with indel support
- [ ] Update merge evaluation logic
- [ ] Add `merged_ric` field to ConsensusInfo
- [ ] Update header writing to include merged_ric
- [ ] Add backward compatibility aliases
- [ ] Write comprehensive tests
- [ ] Update documentation
- [ ] Test on real data with indels

### Phase 3 (2-3 weeks)
- [ ] Implement `run_spoa_msa()` function
- [ ] Implement `replace_terminal_gaps()`
- [ ] Implement `analyze_msa_columns()`
- [ ] Implement `generate_subsets_by_total_size()`
- [ ] Implement `is_compatible_subset()`
- [ ] Implement `create_consensus_from_msa()` with size-weighted voting
- [ ] Implement `merge_group_with_msa()` main algorithm
- [ ] Refactor `process_single_specimen()` to do HAC before merging
- [ ] Add `--output-raw-variants` support
- [ ] Add `--merge-strategy` parameter (iterative vs msa)
- [ ] Implement parameter renaming with aliases
- [ ] Write comprehensive unit tests
- [ ] Write integration tests
- [ ] Performance testing and optimization
- [ ] Update all documentation
- [ ] Test on diverse real data

---

## Open Questions / Decisions Needed

1. **Subset generation optimization:** For N>18, use exhaustive or greedy approximation?
   - Recommendation: Exhaustive unless performance issues observed

2. **IUPAC frequency threshold:** Always include all variants, or add `--merge-min-allele-freq`?
   - Recommendation: Start without threshold, add if users request it

3. **Indel ambiguity edge cases:** What if size-weighted vote is tied (50/50 gap vs base)?
   - Recommendation: Include base (conservative - preserve information)

4. **Phase 3 rollout:** Implement as `--merge-strategy=msa` experimental first?
   - Recommendation: Yes, test thoroughly before making default

5. **Version numbering:**
   - Phase 1: 0.3.6 (patch)
   - Phase 2: 0.4.0 (minor - new features)
   - Phase 3: 0.5.0 (minor - architecture change)

---

## References

- Current implementation: `speconsense/summarize.py`
- MSA reference: `msalign.py` (SPOA integration example)
- IUPAC codes: `IUPAC_CODES` dict (line 25)
- Adjusted identity: `adjusted_identity` package
- SPOA documentation: https://github.com/rvaser/spoa

---

## Notes

- All line numbers reference current `speconsense/summarize.py` (as of 2025-01-21)
- Estimated efforts assume single developer, include testing and documentation
- Real-world performance may vary; monitor on actual ONT datasets
- Consider adding `--debug-msa` flag to output aligned sequences for troubleshooting
