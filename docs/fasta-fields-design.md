# Configurable FASTA Field Output Design

**Date:** 2025-10-22
**Status:** Proposed for review

## Overview

Add a `--fasta-fields` option to `speconsense-summarize` allowing users to customize which metadata fields appear in FASTA headers, similar to BLAST's `-outfmt` or vsearch's `--userfields`.

## Motivation

**Current behavior:** All available fields are written to FASTA headers
```
>ONT01.01-A01-...-1 size=638 ric=638 merged_ric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC
```

**Use cases for customization:**
1. **Downstream tool compatibility** - Some tools expect minimal headers or specific formats
2. **Readability** - Users may want only essential fields for cleaner output
3. **File size** - Minimal headers reduce file size for large datasets
4. **Workflow-specific needs** - Different analyses need different metadata

## Design

### Design Principles

1. **Field codes = field names:** For clarity and consistency, the command-line field codes are identical to the field names that appear in FASTA output. No translation or mapping needed.

2. **Concise field names:** Field names are shortened for readability while maintaining clarity:
   - `rawric` (not `merged_ric`) - emphasizes these are from .raw files
   - `p50diff`, `p95diff` (not `median_diff`, `p95_diff`) - consistent percentile naming
   - `ric` stays as-is (established convention)

3. **Upstream consistency:** Stability field names changed in speconsense (core.py) to match:
   - Write `p50diff` and `p95diff` instead of `median_diff` and `p95_diff`
   - No backward compatibility issues (debug files only)

4. **Preset composition:** Presets can be combined with union semantics:
   - `--fasta-fields minimal,qc` expands to all fields from both presets
   - Duplicates removed, order preserved from left to right

5. **Backward compatibility:** Default behavior unchanged - existing scripts continue to work

6. **Conditional fields:** Don't output fields with null/zero/empty values (keeps headers clean)

### Command-line Interface

```bash
--fasta-fields SPEC
```

**SPEC can be:**
- A **preset name** (keyword): `default`, `minimal`, `full`, `qc`, `id-only`
- A **custom field list**: comma-separated field names (e.g., `size,ric,primers` or `size,ric,p50diff,p95diff`)
- A **combination**: presets and fields mixed (e.g., `minimal,p50diff,p95diff` or `minimal,qc`)
- Field order: left to right, duplicates removed

### Field Codes

**Design principle:** Field codes are identical to field names for clarity and consistency.

| Field Code | Description | Example Output | Always Present? | Included in Presets? |
|------------|-------------|----------------|-----------------|---------------------|
| `size` | Total reads across merged variants | `size=638` | Yes | default, minimal, full, qc |
| `ric` | Reads in Consensus (total) | `ric=638` | Yes | default, minimal, full, qc |
| `length` | Sequence length in bases | `length=589` | Yes | qc, full |
| `rawric` | RiC values of .raw source variants | `rawric=333+305` | Only when merged | default, full |
| `snp` | Number of IUPAC ambiguity positions | `snp=1` | Only when >0 | default, full |
| `p50diff` | Median (p50) edit distance from stability | `p50diff=0.0` | Only when available | qc, full |
| `p95diff` | 95th percentile edit distance | `p95diff=2.0` | Only when available | qc, full |
| `primers` | Detected primer names | `primers=5'-ITS1F,3'-ITS4_RC` | Only when detected | default, full |
| `group` | Variant group number | `group=1` | Yes | _(none - manual only)_ |
| `variant` | Variant identifier within group | `variant=v1` | Only for variants | _(none - manual only)_ |

**Notes:**
- Conditional fields (marked "Only when...") are omitted if not applicable - keeps headers clean
- Stability fields (`p50diff`, `p95diff`) are written by speconsense but currently dropped by speconsense-summarize - this makes them available again
- `group` and `variant` fields are available for selection but not included in any preset (info is in sample_name)

### Presets

| Preset | Fields Included | Use Case |
|--------|-----------------|----------|
| `default` | `size,ric,rawric,snp,primers` | General use, current behavior (backward compatible) |
| `minimal` | `size,ric` | Quick inspection, essential info only |
| `qc` | `size,ric,length,p50diff,p95diff` | Quality control with stability metrics |
| `full` | `size,ric,length,rawric,snp,p50diff,p95diff,primers` | Complete metadata (excludes group/variant) |
| `id-only` | _(none)_ | Clean IDs only, no metadata fields |

**Preset composition examples:**
- `minimal,qc` → `size,ric,length,p50diff,p95diff` (union of both)
- `minimal,primers` → `size,ric,primers` (preset + individual field)
- `default,group,variant` → default fields + group and variant info

**Note:** `default` maintains current output format (excludes stability/length/group/variant)

### Examples

```bash
# Current behavior (default preset)
speconsense-summarize --min-ric 3
speconsense-summarize --min-ric 3 --fasta-fields default
# Output: >ONT01.01-A01-...-1 size=638 ric=638 rawric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC

# Minimal headers - just the essentials
speconsense-summarize --min-ric 3 --fasta-fields minimal
# Output: >ONT01.01-A01-...-1 size=638 ric=638

# QC preset - includes stability metrics and length
speconsense-summarize --min-ric 3 --fasta-fields qc
# Output: >ONT01.01-A01-...-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=0.0

# Full metadata
speconsense-summarize --min-ric 3 --fasta-fields full
# Output: >ONT01.01-...-1 size=638 ric=638 length=589 rawric=333+305 snp=1 p50diff=0.0 p95diff=0.0 primers=...

# ID only (no metadata fields)
speconsense-summarize --min-ric 3 --fasta-fields id-only
# Output: >ONT01.01-A01-...-1

# Custom field selection (codes = field names)
speconsense-summarize --min-ric 3 --fasta-fields size,ric,primers
# Output: >ONT01.01-A01-...-1 size=638 ric=638 primers=5'-ITS1F,3'-ITS4_RC

# Combining presets (union semantics)
speconsense-summarize --min-ric 3 --fasta-fields minimal,qc
# Output: >ONT01.01-A01-...-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=0.0

# Preset + individual fields
speconsense-summarize --min-ric 3 --fasta-fields minimal,p50diff,p95diff
# Output: >ONT01.01-A01-...-1 size=638 ric=638 p50diff=0.0 p95diff=0.0

# Add group/variant info to default
speconsense-summarize --min-ric 3 --fasta-fields default,group,variant
# Output: >ONT01.01-...-1 size=638 ric=638 rawric=333+305 snp=1 primers=... group=1
# (variant field only present for within-group variants)

# Just stability metrics
speconsense-summarize --min-ric 3 --fasta-fields p50diff,p95diff
# Output: >ONT01.01-A01-...-1 p50diff=0.0 p95diff=0.0
```

## Implementation Plan

### 1. Update speconsense (core.py) to Write New Field Names

**Location:** speconsense/core.py, line ~419-422

```python
# OLD:
if median_diff is not None:
    info_parts.append(f"median_diff={median_diff:.1f}")
if p95_diff is not None:
    info_parts.append(f"p95_diff={p95_diff:.1f}")

# NEW:
if median_diff is not None:
    info_parts.append(f"p50diff={median_diff:.1f}")
if p95_diff is not None:
    info_parts.append(f"p95diff={p95_diff:.1f}")
```

**Note:** Variable names stay the same internally (`median_diff`, `p95_diff`), only output changes to `p50diff`, `p95diff`

### 2. Update ConsensusInfo to Include Stability Fields

**Location:** speconsense/summarize.py, line ~68

```python
class ConsensusInfo(NamedTuple):
    """Information about a consensus sequence from speconsense output."""
    sample_name: str
    cluster_id: str
    sequence: str
    ric: int
    size: int
    file_path: str
    snp_count: Optional[int] = None
    primers: Optional[List[str]] = None
    merged_ric: Optional[List[int]] = None
    p50_diff: Optional[float] = None  # NEW - median from stability assessment
    p95_diff: Optional[float] = None  # NEW - 95th percentile from stability
```

### 3. Update parse_consensus_header to Extract Stability Fields

**Location:** speconsense/summarize.py, line ~174

```python
def parse_consensus_header(header: str) -> Tuple[Optional[str], Optional[int], Optional[int],
                                                   Optional[List[str]], Optional[float], Optional[float]]:
    """
    Extract information from Speconsense consensus FASTA header.
    Now includes optional stability metrics (p50diff, p95diff).
    """
    # ... existing parsing code ...

    # Extract stability metrics with NEW field names
    p50_diff_match = re.search(r'p50diff=([\d.]+)', info_string)
    p50_diff = float(p50_diff_match.group(1)) if p50_diff_match else None

    p95_diff_match = re.search(r'p95diff=([\d.]+)', info_string)
    p95_diff = float(p95_diff_match.group(1)) if p95_diff_match else None

    return sample_name, ric, size, primers, p50_diff, p95_diff
```

**Note:** Must also update all call sites of `parse_consensus_header()` to handle the two new return values

### 4. Define Field Registry

**Location:** After `ConsensusInfo` class (around line 80 in summarize.py)

```python
class FastaField:
    """Base class for FASTA header field definitions."""

    def __init__(self, name: str, description: str):
        self.name = name  # Field name (same as code for clarity)
        self.description = description

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        """Format field value for this consensus. Returns None if not applicable."""
        raise NotImplementedError


class SizeField(FastaField):
    def __init__(self):
        super().__init__('size', 'Total reads')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"size={consensus.size}"


class RicField(FastaField):
    def __init__(self):
        super().__init__('ric', 'Reads in consensus')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"ric={consensus.ric}"


class LengthField(FastaField):
    def __init__(self):
        super().__init__('length', 'Sequence length')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"length={len(consensus.sequence)}"


class RawRicField(FastaField):
    def __init__(self):
        super().__init__('rawric', 'Raw variant RiC values')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.merged_ric and len(consensus.merged_ric) > 0:
            ric_values = sorted(consensus.merged_ric, reverse=True)
            return f"rawric={'+'.join(str(r) for r in ric_values)}"
        return None


class SnpField(FastaField):
    def __init__(self):
        super().__init__('snp', 'SNP count')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.snp_count is not None and consensus.snp_count > 0:
            return f"snp={consensus.snp_count}"
        return None


class P50DiffField(FastaField):
    def __init__(self):
        super().__init__('p50diff', 'Median stability difference')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.p50_diff is not None:
            return f"p50diff={consensus.p50_diff:.1f}"
        return None


class P95DiffField(FastaField):
    def __init__(self):
        super().__init__('p95diff', '95th percentile stability difference')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.p95_diff is not None:
            return f"p95diff={consensus.p95_diff:.1f}"
        return None


class PrimersField(FastaField):
    def __init__(self):
        super().__init__('primers', 'Detected primers')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.primers:
            return f"primers={','.join(consensus.primers)}"
        return None


class GroupField(FastaField):
    def __init__(self):
        super().__init__('group', 'Variant group number')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Extract from sample_name (e.g., "...-1" or "...-2.v1")
        match = re.search(r'-(\d+)(?:\.v\d+)?$', consensus.sample_name)
        if match:
            return f"group={match.group(1)}"
        return None


class VariantField(FastaField):
    def __init__(self):
        super().__init__('variant', 'Variant within group')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Extract from sample_name (e.g., "...-1.v1" -> "v1")
        match = re.search(r'\.(v\d+)$', consensus.sample_name)
        if match:
            return f"variant={match.group(1)}"
        return None  # Not a variant


# Field registry - field name is the key (codes = names)
FASTA_FIELDS = {
    'size': SizeField(),
    'ric': RicField(),
    'length': LengthField(),
    'rawric': RawRicField(),
    'snp': SnpField(),
    'p50diff': P50DiffField(),
    'p95diff': P95DiffField(),
    'primers': PrimersField(),
    'group': GroupField(),
    'variant': VariantField(),
}

# Preset definitions
FASTA_FIELD_PRESETS = {
    'default': ['size', 'ric', 'rawric', 'snp', 'primers'],
    'minimal': ['size', 'ric'],
    'qc': ['size', 'ric', 'length', 'p50diff', 'p95diff'],
    'full': ['size', 'ric', 'length', 'rawric', 'snp', 'p50diff', 'p95diff', 'primers'],
    'id-only': [],
}
```

### 5. Add CLI Argument

**Location:** `parse_arguments()` function

```python
parser.add_argument("--fasta-fields", type=str, default="default",
                    help="FASTA header fields to output. Can be: "
                         "(1) a preset name (default, minimal, qc, full, id-only), "
                         "(2) comma-separated field names (size, ric, length, rawric, "
                         "snp, p50diff, p95diff, primers, group, variant), or "
                         "(3) a combination of presets and fields (e.g., minimal,qc or "
                         "minimal,p50diff,p95diff). Duplicates are removed, order preserved "
                         "left to right. Default: default")
```

### 6. Parse Field Specification with Preset Composition

**New function:** `parse_fasta_fields(spec: str) -> List[FastaField]`

```python
def parse_fasta_fields(spec: str) -> List[FastaField]:
    """
    Parse --fasta-fields specification into list of field objects.
    Supports preset composition with union semantics.

    Args:
        spec: Comma-separated list of preset names and/or field names
              Examples:
                - "default" (single preset)
                - "minimal,qc" (preset union)
                - "size,ric,primers" (field list)
                - "minimal,p50diff,p95diff" (preset + fields)

    Returns:
        List of FastaField objects in specified order, duplicates removed

    Raises:
        ValueError: If spec contains unknown preset or field names
    """
    spec = spec.strip().lower()
    if not spec:
        raise ValueError("Empty field specification")

    # Parse comma-separated items (can be presets or field names)
    items = [item.strip() for item in spec.split(',')]

    # Expand presets and collect all field names, preserving order
    all_field_names = []
    seen = set()  # Track duplicates

    for item in items:
        # Check if it's a preset
        if item in FASTA_FIELD_PRESETS:
            # Expand preset
            for field_name in FASTA_FIELD_PRESETS[item]:
                if field_name not in seen:
                    all_field_names.append(field_name)
                    seen.add(field_name)
        elif item in FASTA_FIELDS:
            # It's a field name
            if item not in seen:
                all_field_names.append(item)
                seen.add(item)
        else:
            # Unknown item - provide helpful error
            available_fields = ', '.join(sorted(FASTA_FIELDS.keys()))
            available_presets = ', '.join(sorted(FASTA_FIELD_PRESETS.keys()))
            raise ValueError(
                f"Unknown preset or field name: '{item}'\n"
                f"Available presets: {available_presets}\n"
                f"Available fields: {available_fields}"
            )

    # Convert field names to field objects
    fields = [FASTA_FIELDS[name] for name in all_field_names]

    return fields
```

### 7. Create Header Formatting Function

**New function:** `format_fasta_header(consensus: ConsensusInfo, fields: List[FastaField]) -> str`

```python
def format_fasta_header(consensus: ConsensusInfo, fields: List[FastaField]) -> str:
    """
    Format FASTA header with specified fields.

    Args:
        consensus: Consensus information
        fields: List of fields to include (in order)

    Returns:
        Formatted header line (without leading '>')
    """
    parts = [consensus.sample_name]

    for field in fields:
        value = field.format_value(consensus)
        if value is not None:  # Skip fields that aren't applicable
            parts.append(value)

    return ' '.join(parts)
```

### 8. Update File Writing Functions

**Modify:** `write_specimen_data_files()` and `write_output_files()`

Replace current header construction (appears in multiple places):
```python
# OLD (example from write_specimen_data_files, line ~1176):
header_parts = [f"size={consensus.size}", f"ric={consensus.ric}"]
if consensus.merged_ric and len(consensus.merged_ric) > 0:
    ric_values = sorted(consensus.merged_ric, reverse=True)
    header_parts.append(f"merged_ric={'+'.join(str(r) for r in ric_values)}")
if consensus.snp_count is not None and consensus.snp_count > 0:
    header_parts.append(f"snp={consensus.snp_count}")
if consensus.primers:
    header_parts.append(f"primers={','.join(consensus.primers)}")
f.write(f">{consensus.sample_name} {' '.join(header_parts)}\n")
```

With:
```python
# NEW:
header = format_fasta_header(consensus, fasta_fields)
f.write(f">{header}\n")
```

**Note:** This replacement occurs in 4 places:
1. Individual FASTA files (final consensus) - line ~1185
2. .raw FASTA files - line ~1200
3. summary.fasta (final consensus) - line ~1313
4. summary.fasta (.raw files) - line ~1321

### 9. Thread Fields Specification Through Call Chain

Add `fasta_fields` parameter to:
- `write_specimen_data_files()` - receives from main()
- `write_output_files()` - receives from main()
- `main()` - parses at startup, passes to write functions

**In main():**
```python
# Parse field specification early
try:
    fasta_fields = parse_fasta_fields(args.fasta_fields)
except ValueError as e:
    logging.error(str(e))
    sys.exit(1)

# Later, in processing loop:
specimen_raw_consensuses = write_specimen_data_files(
    final_consensus, merge_traceability, naming_info,
    args.summary_dir, fastq_dir,
    fastq_lookup, original_consensus_lookup,
    fasta_fields  # NEW parameter
)

# And at end:
write_output_files(
    all_final_consensus, all_raw_consensuses,
    args.summary_dir, temp_log_file.name,
    fasta_fields  # NEW parameter
)
```

## New Fields Details

### `length` - Sequence Length

**Rationale:** Common metadata in FASTA files, useful for filtering/QC

**Implementation:**
```python
class LengthField(FastaField):
    def __init__(self):
        super().__init__('length', 'length', 'Sequence length', True)

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"length={len(consensus.sequence)}"
```

### `group` - Group Number

**Rationale:** Expose the variant group number explicitly (currently embedded in sample_name)

**Implementation:**
```python
class GroupField(FastaField):
    def __init__(self):
        super().__init__('group', 'group', 'Variant group number', True)

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Extract from sample_name (e.g., "...-1" or "...-2.v1")
        match = re.search(r'-(\d+)(?:\.v\d+)?$', consensus.sample_name)
        if match:
            return f"group={match.group(1)}"
        return None
```

### `variant` - Variant Identifier

**Rationale:** Distinguish between merged consensus and within-group variants

**Implementation:**
```python
class VariantField(FastaField):
    def __init__(self):
        super().__init__('variant', 'variant', 'Variant within group', False)

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Extract from sample_name (e.g., "...-1.v1" -> "v1")
        match = re.search(r'\.(v\d+)$', consensus.sample_name)
        if match:
            return f"variant={match.group(1)}"
        return None  # Not a variant, it's the main merged consensus
```

## Testing Strategy

### Unit Tests
1. **Field parsing:**
   - Valid presets resolve correctly
   - Custom field lists parse correctly
   - Invalid field codes raise clear errors
   - Field order is preserved

2. **Header formatting:**
   - Each field formats correctly
   - Conditional fields omitted when not applicable
   - Empty field list produces ID-only output

3. **Integration:**
   - All presets produce valid output
   - Custom combinations work correctly

### Functional Tests

**Test dataset:** Use existing test data (ont37/demux1020)

```bash
# Test all presets
for preset in default minimal full id-only blast-like; do
    speconsense-summarize --source clusters --summary-dir test_${preset} \
        --fasta-fields ${preset} --min-ric 3
    # Verify output format
    head test_${preset}/summary.fasta
done

# Test custom combinations
speconsense-summarize --fasta-fields size,ric,length --min-ric 3
speconsense-summarize --fasta-fields primers,snp --min-ric 3
```

### Validation Checks
1. Output files parseable by BioPython SeqIO
2. Header format matches specification
3. Backward compatibility: default preset matches current behavior
4. File size comparison: minimal < default < full

## Backward Compatibility

**Default behavior unchanged:**
- `--fasta-fields default` produces current output
- Omitting `--fasta-fields` uses "default" preset
- Existing scripts continue to work

**Migration path:**
- Users wanting old behavior: no changes needed
- Users wanting custom output: add `--fasta-fields <spec>`

## Future Extensibility

**Adding new fields is straightforward:**
1. Define new `FastaField` subclass
2. Add to `FASTA_FIELDS` registry
3. Update `full` preset if appropriate
4. Document in help text

**Example future fields:**
- `cluster_ids`: Original cluster IDs that were merged
- `identity`: Sequence identity metrics
- `coverage`: Read coverage statistics
- `taxonomy`: Taxonomic assignment (if integrated with classifier)

## Estimated Effort

| Task | Time |
|------|------|
| Update ConsensusInfo and parse_consensus_header for stability fields | 0.5 hours |
| Define field registry and classes | 1.5 hours |
| Implement parsing and validation | 1 hour |
| Update write functions | 1 hour |
| Implement new fields (length, group, variant) | 1 hour |
| Testing and validation | 1.5 hours |
| Documentation updates | 0.5 hours |
| **Total** | **~7 hours** |

## Alternative Designs Considered

### 1. BLAST-style numbered formats
```bash
--fasta-fields 0  # default
--fasta-fields 1  # minimal
--fasta-fields "6 size ric primers"  # custom
```
**Rejected:** Numbers are opaque; keyword presets are clearer

### 2. Separate flags per field
```bash
--include-size --include-ric --no-primers
```
**Rejected:** Too many flags; less flexible; harder to specify order

### 3. Template strings
```bash
--fasta-header-template "{id} size={size} ric={ric}"
```
**Rejected:** More complex to implement; harder to validate; overkill for this use case

## Future Extensions (Optional)

These features are not part of the initial implementation but could be added later if needed:

1. **Custom field separators:**
   ```bash
   --fasta-field-sep " "  # space (default)
   --fasta-field-sep "\t"  # tab-delimited
   ```

2. **Value-only mode:** For tools expecting positional values
   ```bash
   --fasta-fields-format "key=value"  # default: size=638
   --fasta-fields-format "value-only"  # just: 638
   ```

3. **Additional computed fields:**
   - `gc_content`: GC percentage of sequence
   - `complexity`: Sequence complexity metric
   - `identity`: Sequence identity to reference (if available)

## Output Format Changes

### Before (current):
```
>ONT01.01-A01-...-1 size=638 ric=638 merged_ric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC
```

### After (with default preset):
```
>ONT01.01-A01-...-1 size=638 ric=638 rawric=333+305 snp=1 primers=5'-ITS1F,3'-ITS4_RC
```

**Changes:**
- `merged_ric` → `rawric` (emphasizes these are .raw file RiC values)

### New Capabilities:

```bash
# Include stability metrics from speconsense
--fasta-fields qc
>ONT01.01-A01-...-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=0.0

# Combine presets
--fasta-fields minimal,qc
>ONT01.01-A01-...-1 size=638 ric=638 length=589 p50diff=0.0 p95diff=0.0

# Add specific fields to a preset
--fasta-fields default,group,variant
>ONT01.01-A01-...-1 size=638 ric=638 rawric=333+305 snp=1 primers=... group=1

# Just the essentials
--fasta-fields minimal
>ONT01.01-A01-...-1 size=638 ric=638
```

## Summary

This design provides:
- ✅ Flexible field selection (presets + custom lists)
- ✅ Clean, intuitive syntax (field codes = field names)
- ✅ Backward compatible (default preset unchanged)
- ✅ Stability fields available but optional
- ✅ Extensible for future fields (taxonomy, coverage, etc.)
- ✅ Similar to familiar tools (BLAST, vsearch)
- ✅ Clear error messages with helpful suggestions
- ✅ Moderate implementation effort (~7 hours including stability field extraction)

**Key Design Decisions:**
1. **Field codes = field names** - No translation needed, maximum clarity
2. **Concise, consistent naming:**
   - `rawric` (not `merged_ric`) - emphasizes .raw file origin
   - `p50diff`, `p95diff` (not `median_diff`, `p95_diff`) - consistent percentile notation
   - `ric` unchanged (established convention)
3. **Upstream consistency** - speconsense (core.py) writes `p50diff`/`p95diff` to match
4. **Preset composition** - Combine presets and fields: `minimal,qc` or `default,group,variant`
5. **Stability fields available** - Optional QC data, not in default preset
6. **Conditional output** - Omit null/empty fields for clean headers
7. **Manual-only fields** - `length`, `group`, `variant` available but not in any preset
