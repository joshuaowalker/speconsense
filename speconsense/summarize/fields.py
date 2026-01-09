"""FASTA field customization support for speconsense-summarize.

Allows users to control which metadata fields appear in FASTA headers.
"""

import re
from typing import List, Optional

from speconsense.types import ConsensusInfo


class FastaField:
    """Base class for FASTA header field definitions."""

    def __init__(self, name: str, description: str):
        self.name = name  # Field name (matches field code for clarity)
        self.description = description

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        """Format field value for this consensus. Returns None if not applicable."""
        raise NotImplementedError


class SizeField(FastaField):
    def __init__(self):
        super().__init__('size', 'Total reads across merged variants')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"size={consensus.size}"


class RicField(FastaField):
    def __init__(self):
        super().__init__('ric', 'Reads in consensus')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"ric={consensus.ric}"


class LengthField(FastaField):
    def __init__(self):
        super().__init__('length', 'Sequence length in bases')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        return f"length={len(consensus.sequence)}"


class RawRicField(FastaField):
    def __init__(self):
        super().__init__('rawric', 'RiC values of .raw source variants')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.raw_ric and len(consensus.raw_ric) > 0:
            ric_values = sorted(consensus.raw_ric, reverse=True)
            return f"rawric={'+'.join(str(r) for r in ric_values)}"
        return None


class RawLenField(FastaField):
    def __init__(self):
        super().__init__('rawlen', 'Lengths of merged source sequences')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.raw_len and len(consensus.raw_len) > 0:
            len_values = sorted(consensus.raw_len, reverse=True)
            return f"rawlen={'+'.join(str(l) for l in len_values)}"
        return None


class SnpField(FastaField):
    def __init__(self):
        super().__init__('snp', 'Number of IUPAC ambiguity positions from merging')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.snp_count is not None and consensus.snp_count > 0:
            return f"snp={consensus.snp_count}"
        return None


class AmbigField(FastaField):
    def __init__(self):
        super().__init__('ambig', 'Count of IUPAC ambiguity codes in consensus')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Count non-ACGT characters in the sequence
        ambig_count = sum(1 for c in consensus.sequence if c.upper() not in 'ACGT')
        if ambig_count > 0:
            return f"ambig={ambig_count}"
        return None


class PrimersField(FastaField):
    def __init__(self):
        super().__init__('primers', 'Detected primer names')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.primers:
            return f"primers={','.join(consensus.primers)}"
        return None


class RidField(FastaField):
    def __init__(self):
        super().__init__('rid', 'Mean read identity (percentage)')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.rid is not None:
            return f"rid={consensus.rid*100:.1f}"
        return None


class RidMinField(FastaField):
    def __init__(self):
        super().__init__('rid_min', 'Minimum read identity (percentage)')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.rid_min is not None:
            return f"rid_min={consensus.rid_min*100:.1f}"
        return None


class GroupField(FastaField):
    def __init__(self):
        super().__init__('group', 'Variant group number')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Extract from sample_name (e.g., "...-1.v1" or "...-2.v1.raw1")
        match = re.search(r'-(\d+)\.v\d+(?:\.raw\d+)?$', consensus.sample_name)
        if match:
            return f"group={match.group(1)}"
        return None


class VariantField(FastaField):
    def __init__(self):
        super().__init__('variant', 'Variant identifier within group')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        # Extract from sample_name (e.g., "...-1.v1" -> "v1" or "...-1.v1.raw1" -> "v1")
        match = re.search(r'\.(v\d+)(?:\.raw\d+)?$', consensus.sample_name)
        if match:
            return f"variant={match.group(1)}"
        return None


# Field registry - field name is the key (codes = names)
FASTA_FIELDS = {
    'size': SizeField(),
    'ric': RicField(),
    'length': LengthField(),
    'rawric': RawRicField(),
    'rawlen': RawLenField(),
    'snp': SnpField(),
    'ambig': AmbigField(),
    'rid': RidField(),
    'rid_min': RidMinField(),
    'primers': PrimersField(),
    'group': GroupField(),
    'variant': VariantField(),
}

# Preset definitions
FASTA_FIELD_PRESETS = {
    'default': ['size', 'ric', 'rawric', 'rawlen', 'snp', 'ambig', 'primers'],
    'minimal': ['size', 'ric'],
    'qc': ['size', 'ric', 'length', 'rid', 'ambig'],
    'full': ['size', 'ric', 'length', 'rawric', 'rawlen', 'snp', 'ambig', 'rid', 'primers'],
    'id-only': [],
}


def validate_field_registry():
    """Validate that all preset fields exist in registry."""
    for preset_name, field_names in FASTA_FIELD_PRESETS.items():
        for field_name in field_names:
            if field_name not in FASTA_FIELDS:
                raise ValueError(f"Preset '{preset_name}' references unknown field '{field_name}'")


# Validate at module load
validate_field_registry()


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
                - "minimal,rid" (preset + fields)

    Returns:
        List of FastaField objects in specified order, duplicates removed

    Raises:
        ValueError: If spec contains unknown preset or field names
    """
    spec = spec.strip().lower()
    if not spec:
        # Default to "default" preset if empty
        spec = "default"

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
                f"  Available presets: {available_presets}\n"
                f"  Available fields: {available_fields}"
            )

    # Convert field names to field objects
    fields = [FASTA_FIELDS[name] for name in all_field_names]

    return fields


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
