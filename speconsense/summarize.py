#!/usr/bin/env python3

import os
import re
import glob
import csv
import json
import shutil
import argparse
import itertools
import logging
import datetime
from typing import List, Dict, Tuple, Optional, NamedTuple
from collections import defaultdict

import edlib
from Bio import SeqIO
from adjusted_identity import score_alignment, AdjustmentParams
import tempfile
import subprocess
from tqdm import tqdm
import numpy as np
from io import StringIO

# Import homopolymer-aware alignment functions from core
from speconsense.core import extract_alignments_from_msa, ReadAlignment, analyze_positional_variation


# IUPAC nucleotide ambiguity codes mapping
IUPAC_CODES = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N',
}

# IUPAC equivalencies for edlib alignment
# This allows edlib to treat IUPAC ambiguity codes as matching their constituent bases
IUPAC_EQUIV = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"),
               ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
               ("W", "A"), ("W", "T"), ("M", "A"), ("M", "C"),
               ("S", "C"), ("S", "G"), ("K", "G"), ("K", "T"),
               ("B", "C"), ("B", "G"), ("B", "T"),
               ("D", "A"), ("D", "G"), ("D", "T"),
               ("H", "A"), ("H", "C"), ("H", "T"),
               ("V", "A"), ("V", "C"), ("V", "G"), ]

# Standard adjustment parameters for consistent sequence comparison
# Used by both substitution distance calculation and adjusted identity distance
STANDARD_ADJUSTMENT_PARAMS = AdjustmentParams(
    normalize_homopolymers=True,    # Enable homopolymer normalization
    handle_iupac_overlap=False,     # Disable IUPAC overlap - use standard IUPAC semantics (Y≠M)
    normalize_indels=False,         # Disable indel normalization
    end_skip_distance=0,            # No end trimming - sequences must match end-to-end
    max_repeat_motif_length=1       # Single-base repeats for homopolymer normalization
)

# Maximum number of variants to evaluate for MSA-based merging
# Beyond this limit, subset evaluation becomes computationally impractical
# (2^N subsets to evaluate - 2^10 = 1024 is manageable, 2^20 = 1M+ is not)
MAX_MSA_MERGE_VARIANTS = 10


class ConsensusInfo(NamedTuple):
    """Information about a consensus sequence from speconsense output."""
    sample_name: str
    cluster_id: str
    sequence: str
    ric: int
    size: int
    file_path: str
    snp_count: Optional[int] = None  # Number of SNPs from IUPAC consensus generation
    primers: Optional[List[str]] = None  # List of detected primer names
    raw_ric: Optional[List[int]] = None  # RiC values of .raw source variants
    rid: Optional[float] = None  # Mean read identity (internal consistency metric)
    rid_min: Optional[float] = None  # Minimum read identity (worst-case read)


class ClusterQualityData(NamedTuple):
    """Quality metrics for a cluster (no visualization matrix)."""
    consensus_seq: str
    position_error_rates: List[float]  # Per-position error rates (0-1) in consensus space
    position_error_counts: List[int]  # Per-position error counts in consensus space
    read_identities: List[float]  # Per-read identity scores (0-1)
    position_stats: Optional[List] = None  # Detailed PositionStats for debugging (optional)


# FASTA field customization support
# Allows users to control which metadata fields appear in FASTA headers

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


class SnpField(FastaField):
    def __init__(self):
        super().__init__('snp', 'Number of IUPAC ambiguity positions')

    def format_value(self, consensus: ConsensusInfo) -> Optional[str]:
        if consensus.snp_count is not None and consensus.snp_count > 0:
            return f"snp={consensus.snp_count}"
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
    'snp': SnpField(),
    'rid': RidField(),
    'rid_min': RidMinField(),
    'primers': PrimersField(),
    'group': GroupField(),
    'variant': VariantField(),
}

# Preset definitions
FASTA_FIELD_PRESETS = {
    'default': ['size', 'ric', 'rawric', 'snp', 'primers'],
    'minimal': ['size', 'ric'],
    'qc': ['size', 'ric', 'length', 'rid'],
    'full': ['size', 'ric', 'length', 'rawric', 'snp', 'rid', 'primers'],
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


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Speconsense output with advanced variant handling.")
    parser.add_argument("--min-ric", type=int, default=3,
                        help="Minimum Reads in Consensus (RiC) threshold (default: 3)")
    parser.add_argument("--source", type=str, default="clusters",
                        help="Source directory containing Speconsense output (default: clusters)")
    parser.add_argument("--summary-dir", type=str, default="__Summary__",
                        help="Output directory for summary files (default: __Summary__)")
    parser.add_argument("--fasta-fields", type=str, default="default",
                        help="FASTA header fields to output. Can be: "
                             "(1) a preset name (default, minimal, qc, full, id-only), "
                             "(2) comma-separated field names (size, ric, length, rawric, "
                             "snp, rid, rid_min, primers, group, variant), or "
                             "(3) a combination of presets and fields (e.g., minimal,qc or "
                             "minimal,rid). Duplicates removed, order preserved "
                             "left to right. Default: default")

    # Merge phase parameters
    parser.add_argument("--merge-snp", action="store_true", default=True,
                        help="Enable SNP-based merging (default: True)")
    parser.add_argument("--merge-indel-length", type=int, default=0,
                        help="Maximum length of individual indels allowed in merging (default: 0 = disabled)")
    parser.add_argument("--merge-position-count", type=int, default=2,
                        help="Maximum total SNP+indel positions allowed in merging (default: 2)")
    parser.add_argument("--merge-min-size-ratio", type=float, default=0.0,
                        help="Minimum size ratio (smaller/larger) for merging clusters (default: 0.0 = disabled)")
    parser.add_argument("--disable-homopolymer-equivalence", action="store_true",
                        help="Disable homopolymer equivalence in merging (treat AAA vs AAAA as different)")

    # Backward compatibility: support old --snp-merge-limit parameter
    parser.add_argument("--snp-merge-limit", type=int, dest="_snp_merge_limit_deprecated",
                        help=argparse.SUPPRESS)  # Hidden but functional

    # Group and selection phase parameters
    parser.add_argument("--group-identity", "--variant-group-identity",
                        dest="group_identity", type=float, default=0.9,
                        help="Identity threshold for variant grouping using HAC (default: 0.9)")
    parser.add_argument("--select-max-variants", "--max-variants",
                        dest="select_max_variants", type=int, default=-1,
                        help="Maximum total variants to output per group (default: -1 = no limit, 0 also means no limit)")
    parser.add_argument("--select-max-groups", "--max-groups",
                        dest="select_max_groups", type=int, default=-1,
                        help="Maximum number of groups to output per specimen (default: -1 = all groups)")
    parser.add_argument("--select-strategy", "--variant-selection",
                        dest="select_strategy", choices=["size", "diversity"], default="size",
                        help="Variant selection strategy: size or diversity (default: size)")

    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level")

    args = parser.parse_args()

    # Handle backward compatibility for deprecated parameters
    import sys
    if args._snp_merge_limit_deprecated is not None:
        if '--snp-merge-limit' in sys.argv:
            logging.warning("--snp-merge-limit is deprecated, use --merge-position-count instead")
        args.merge_position_count = args._snp_merge_limit_deprecated

    if '--variant-group-identity' in sys.argv:
        logging.warning("--variant-group-identity is deprecated, use --group-identity instead")

    if '--max-variants' in sys.argv:
        logging.warning("--max-variants is deprecated, use --select-max-variants instead")

    if '--max-groups' in sys.argv:
        logging.warning("--max-groups is deprecated, use --select-max-groups instead")

    if '--variant-selection' in sys.argv:
        logging.warning("--variant-selection is deprecated, use --select-strategy instead")

    return args


def setup_logging(log_level: str, log_file: str = None):
    """Setup logging configuration with optional file output."""
    # Clear any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    # Set up root logger
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, log_level))
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return log_file
    
    return None


def parse_consensus_header(header: str) -> Tuple[Optional[str], Optional[int], Optional[int],
                                                   Optional[List[str]], Optional[float], Optional[float],
                                                   Optional[bool], Optional[int]]:
    """
    Extract information from Speconsense consensus FASTA header.

    Parses read identity metrics.

    Returns:
        Tuple of (sample_name, ric, size, primers, rid, rid_min)
    """
    sample_match = re.match(r'>([^ ]+) (.+)', header)
    if not sample_match:
        return None, None, None, None, None, None

    sample_name = sample_match.group(1)
    info_string = sample_match.group(2)

    # Extract RiC value
    ric_match = re.search(r'ric=(\d+)', info_string)
    ric = int(ric_match.group(1)) if ric_match else 0

    # Extract size value
    size_match = re.search(r'size=(\d+)', info_string)
    size = int(size_match.group(1)) if size_match else 0

    # Extract primers value
    primers_match = re.search(r'primers=([^,\s]+(?:,[^,\s]+)*)', info_string)
    primers = primers_match.group(1).split(',') if primers_match else None

    # Extract read identity metrics (percentages in headers, convert to fractions)
    rid_match = re.search(r'rid=([\d.]+)', info_string)
    rid = float(rid_match.group(1)) / 100.0 if rid_match else None

    rid_min_match = re.search(r'rid_min=([\d.]+)', info_string)
    rid_min = float(rid_min_match.group(1)) / 100.0 if rid_min_match else None

    return sample_name, ric, size, primers, rid, rid_min


def load_consensus_sequences(source_folder: str, min_ric: int) -> List[ConsensusInfo]:
    """Load all consensus sequences from speconsense output files."""
    consensus_list = []

    # Find all consensus FASTA files matching the new naming pattern
    fasta_pattern = os.path.join(source_folder, "*-all.fasta")
    fasta_files = sorted(glob.glob(fasta_pattern))

    for fasta_file in fasta_files:
        logging.debug(f"Processing consensus file: {fasta_file}")

        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                sample_name, ric, size, primers, rid, rid_min = \
                    parse_consensus_header(f">{record.description}")

                if sample_name and ric >= min_ric:
                    # Extract cluster ID from sample name (e.g., "sample-c1" -> "c1")
                    cluster_match = re.search(r'-c(\d+)$', sample_name)
                    cluster_id = cluster_match.group(0) if cluster_match else sample_name

                    consensus_info = ConsensusInfo(
                        sample_name=sample_name,
                        cluster_id=cluster_id,
                        sequence=str(record.seq),
                        ric=ric,
                        size=size,
                        file_path=fasta_file,
                        snp_count=None,  # No SNP info from original speconsense output
                        primers=primers,
                        raw_ric=None,  # Not available in original speconsense output
                        rid=rid,  # Mean read identity if available
                        rid_min=rid_min,  # Minimum read identity if available
                    )
                    consensus_list.append(consensus_info)

    logging.info(f"Loaded {len(consensus_list)} consensus sequences from {len(fasta_files)} files")
    return consensus_list


def load_metadata_from_json(source_folder: str, sample_name: str) -> Optional[Dict]:
    """Load metadata JSON file for a consensus sequence.

    Args:
        source_folder: Source directory containing cluster_debug folder
        sample_name: Sample name (e.g., "sample-c1")

    Returns:
        Dictionary with metadata, or None if file not found or error
    """
    # Construct path to metadata file
    debug_dir = os.path.join(source_folder, "cluster_debug")
    metadata_file = os.path.join(debug_dir, f"{sample_name}-metadata.json")

    if not os.path.exists(metadata_file):
        logging.debug(f"Metadata file not found: {metadata_file}")
        return None

    try:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        return metadata
    except Exception as e:
        logging.warning(f"Failed to load metadata from {metadata_file}: {e}")
        return None


def identify_outliers(final_consensus: List, all_raw_consensuses: List, source_folder: str) -> Dict:
    """Identify sequences with low read identity using statistical outlier detection.

    Flags sequences with mean read identity (rid) below (mean - 2×std) for the dataset.
    This identifies the ~2.5% lowest values that may warrant review.

    Note: rid_min (minimum read identity) is not used because single outlier reads
    don't significantly impact consensus quality. Positional analysis better captures
    systematic issues like mixed clusters or variants.

    Args:
        final_consensus: List of final consensus sequences
        all_raw_consensuses: List of all raw consensus sequences (unused, kept for API compatibility)
        source_folder: Source directory (unused, kept for API compatibility)

    Returns:
        Dictionary with:
        {
            'statistical_outliers': List of (cons, rid),
            'no_issues': List of consensus sequences with good quality,
            'global_stats': {'mean_rid', 'std_rid', 'stat_threshold_rid'}
        }
    """
    # Calculate global statistics for all sequences with identity metrics
    all_rids = []

    for cons in final_consensus:
        if cons.rid is not None:
            all_rids.append(cons.rid)

    # Calculate mean and std for statistical outlier detection
    mean_rid = np.mean(all_rids) if all_rids else 1.0
    std_rid = np.std(all_rids) if len(all_rids) > 1 else 0.0

    # Threshold for statistical outliers (2 standard deviations below mean)
    stat_threshold_rid = mean_rid - 2 * std_rid

    # Categorize sequences
    statistical = []
    no_issues = []

    for cons in final_consensus:
        rid = cons.rid if cons.rid is not None else 1.0

        if rid < stat_threshold_rid:
            statistical.append((cons, rid))
        else:
            no_issues.append(cons)

    return {
        'statistical_outliers': statistical,
        'no_issues': no_issues,
        'global_stats': {
            'mean_rid': mean_rid,
            'std_rid': std_rid,
            'stat_threshold_rid': stat_threshold_rid
        }
    }


def analyze_positional_identity_outliers(
    consensus_info,
    source_folder: str,
    min_variant_frequency: float,
    min_variant_count: int
) -> Optional[Dict]:
    """Analyze positional error rates and identify high-error positions.

    Args:
        consensus_info: ConsensusInfo object for the sequence
        source_folder: Source directory containing cluster_debug folder
        min_variant_frequency: Global threshold for flagging positions (from metadata)
        min_variant_count: Minimum variant count for phasing (from metadata)

    Returns:
        Dictionary with positional analysis:
        {
            'num_outlier_positions': int,
            'mean_outlier_error_rate': float,  # Mean error rate across outlier positions only
            'total_nucleotide_errors': int,    # Sum of error counts at outlier positions
            'outlier_threshold': float,
            'outlier_positions': List of (position, error_rate, error_count) tuples
        }
        Returns None if MSA file not found or analysis fails

        Note: Error rates already exclude homopolymer length differences due to
        homopolymer normalization in analyze_positional_variation()
    """
    # Skip analysis for low-RiC sequences (insufficient data for meaningful statistics)
    # Need at least 2 * min_variant_count to confidently phase two variants
    min_ric_threshold = 2 * min_variant_count
    if consensus_info.ric < min_ric_threshold:
        logging.debug(f"Skipping positional analysis for {consensus_info.sample_name}: "
                     f"RiC {consensus_info.ric} < {min_ric_threshold}")
        return None

    # Construct path to MSA file
    debug_dir = os.path.join(source_folder, "cluster_debug")

    # Try to find the MSA file
    # MSA files use the original cluster naming (e.g., "specimen-c1")
    # not the summarized naming (e.g., "specimen-1.v1")
    msa_file = None

    # Extract specimen name and cluster ID
    # consensus_info.sample_name might be "specimen-1.v1" (summarized)
    # consensus_info.cluster_id should be "-c1" (original cluster)

    # Build the base name from specimen + cluster_id
    # If sample_name is "ONT01.23-...-1.v1" and cluster_id is "-c1"
    # we need to reconstruct "ONT01.23-...-c1"

    sample_name = consensus_info.sample_name
    cluster_id = consensus_info.cluster_id

    # Remove any HAC group/variant suffix from sample_name to get specimen base
    # Pattern: "-\d+\.v\d+" (e.g., "-1.v1")
    specimen_base = re.sub(r'-\d+\.v\d+$', '', sample_name)

    # Reconstruct original cluster name
    original_cluster_name = f"{specimen_base}{cluster_id}"

    # Look for the MSA file with correct extension
    msa_fasta = os.path.join(debug_dir, f"{original_cluster_name}-RiC{consensus_info.ric}-msa.fasta")
    if os.path.exists(msa_fasta):
        msa_file = msa_fasta

    if not msa_file:
        logging.debug(f"No MSA file found for {original_cluster_name}")
        return None

    # Analyze cluster quality using core.py's positional analysis
    quality_data = analyze_cluster_quality(msa_file, consensus_info.sequence)

    if not quality_data or not quality_data.position_error_rates:
        logging.debug(f"Failed to analyze cluster quality for {original_cluster_name}")
        return None

    position_error_rates = quality_data.position_error_rates
    position_error_counts = quality_data.position_error_counts
    position_stats = quality_data.position_stats

    # Use global min_variant_frequency as threshold
    # Positions above this could be undetected/unphased variants
    threshold = min_variant_frequency
    outlier_positions = [
        (i, rate, count)
        for i, (rate, count) in enumerate(zip(position_error_rates, position_error_counts))
        if rate > threshold
    ]

    # Build detailed outlier info including base composition
    outlier_details = []
    if position_stats:
        for i, rate, count in outlier_positions:
            if i < len(position_stats):
                ps = position_stats[i]
                outlier_details.append({
                    'consensus_position': ps.consensus_position,
                    'msa_position': ps.msa_position,
                    'error_rate': rate,
                    'error_count': count,
                    'coverage': ps.coverage,
                    'consensus_nucleotide': ps.consensus_nucleotide,
                    'base_composition': dict(ps.base_composition),
                    'homopolymer_composition': dict(ps.homopolymer_composition) if ps.homopolymer_composition else {},
                    'sub_count': ps.sub_count,
                    'ins_count': ps.ins_count,
                    'del_count': ps.del_count,
                })

    # Calculate statistics for outlier positions only
    if outlier_positions:
        mean_outlier_error = np.mean([rate for _, rate, _ in outlier_positions])
        total_nucleotide_errors = sum(count for _, _, count in outlier_positions)
    else:
        mean_outlier_error = 0.0
        total_nucleotide_errors = 0

    return {
        'num_outlier_positions': len(outlier_positions),
        'mean_outlier_error_rate': mean_outlier_error,
        'total_nucleotide_errors': total_nucleotide_errors,
        'outlier_threshold': threshold,
        'outlier_positions': outlier_positions,
        'outlier_details': outlier_details,
        'consensus_seq': quality_data.consensus_seq,
        'ric': consensus_info.ric,
    }


def run_spoa_msa(sequences: List[str]) -> List:
    """
    Run SPOA to create multiple sequence alignment.

    Args:
        sequences: List of DNA sequence strings

    Returns:
        List of SeqRecord objects with aligned sequences (including gaps)
    """
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

            # Run SPOA with alignment output (-r 2) and global alignment mode (-l 1)
            result = subprocess.run(
                ['spoa', temp_input.name, '-r', '2', '-l', '1'],
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


def identify_indel_events(aligned_seqs: List, alignment_length: int) -> List[Tuple[int, int]]:
    """
    Identify consecutive runs of indel columns (events).

    An indel event is a maximal consecutive run of columns containing gaps.
    Each event represents a single biological insertion or deletion.

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        alignment_length: Length of the alignment

    Returns:
        List of (start_col, end_col) tuples, where end_col is inclusive
    """
    events = []
    in_event = False
    start_col = None

    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        has_gap = '-' in column
        has_bases = any(c != '-' for c in column)

        # Indel column: mix of gaps and bases
        if has_gap and has_bases:
            if not in_event:
                # Start new event
                in_event = True
                start_col = col_idx
        else:
            # Not an indel column (either all gaps or all bases)
            if in_event:
                # End current event
                events.append((start_col, col_idx - 1))
                in_event = False

    # Handle event that extends to end of alignment
    if in_event:
        events.append((start_col, alignment_length - 1))

    return events


def is_homopolymer_event(aligned_seqs: List, start_col: int, end_col: int) -> bool:
    """
    Classify a complete indel event as homopolymer or structural.

    An event is homopolymer if:
    1. All bases in the event region (across all sequences, all columns) are identical
    2. At least one flanking solid column has all sequences showing the same base

    This matches adjusted-identity semantics where AAA ≈ AAAA.

    Examples:
        Homopolymer:  ATAAA--GC vs ATAAAAGC  (event has all A's, flanked by A)
        Structural:   ATAA-GC vs ATG-AGC     (event has A, flanked by A vs G)
        Structural:   ATC--GC vs ATCATGC     (event has A and T - not homopolymer)

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        start_col: First column of the indel event (inclusive)
        end_col: Last column of the indel event (inclusive)

    Returns:
        True if homopolymer event, False if structural
    """
    # Extract all bases from the event region (excluding gaps)
    bases_in_event = set()
    for col_idx in range(start_col, end_col + 1):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        bases_in_event.update(c for c in column if c != '-')

    # Must have exactly one base type across the entire event
    if len(bases_in_event) != 1:
        return False

    event_base = list(bases_in_event)[0]
    alignment_length = len(aligned_seqs[0].seq)

    # Check flanking columns for matching homopolymer context
    # A valid flanking column must:
    # 1. Not be an indel column (all sequences have bases, no gaps)
    # 2. All bases match the event base

    # Check left flank
    if start_col > 0:
        left_col = start_col - 1
        left_column = [str(seq.seq[left_col]) for seq in aligned_seqs]
        left_bases = set(c for c in left_column if c != '-')
        left_has_gap = '-' in left_column

        if not left_has_gap and left_bases == {event_base}:
            return True

    # Check right flank
    if end_col < alignment_length - 1:
        right_col = end_col + 1
        right_column = [str(seq.seq[right_col]) for seq in aligned_seqs]
        right_bases = set(c for c in right_column if c != '-')
        right_has_gap = '-' in right_column

        if not right_has_gap and right_bases == {event_base}:
            return True

    # No valid homopolymer flanking found
    return False


def analyze_msa_columns(aligned_seqs: List) -> dict:
    """
    Analyze aligned sequences to count SNPs and indels.

    Distinguishes between structural indels (real insertions/deletions) and
    homopolymer indels (length differences in homopolymer runs like AAA vs AAAA).

    Uses event-based classification: consecutive indel columns are grouped into
    events, and each complete event is classified as homopolymer or structural.

    Important: All gaps (including terminal gaps) count as variant positions
    since variants within a group share the same primers.

    Returns dict with:
        'snp_count': number of positions with >1 non-gap base
        'structural_indel_count': number of structural indel events
        'structural_indel_length': length of longest structural indel event
        'homopolymer_indel_count': number of homopolymer indel events
        'homopolymer_indel_length': length of longest homopolymer indel event
        'indel_count': total indel events (for backward compatibility)
        'max_indel_length': max indel event length (for backward compatibility)
    """
    alignment_length = len(aligned_seqs[0].seq)

    # Step 1: Count SNPs
    snp_count = 0
    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        unique_bases = set(c for c in column if c != '-')
        has_gap = '-' in column

        # SNP position: multiple different bases with NO gaps
        # Columns with gaps are indels, not SNPs
        if len(unique_bases) > 1 and not has_gap:
            snp_count += 1

    # Step 2: Identify indel events (consecutive runs of indel columns)
    indel_events = identify_indel_events(aligned_seqs, alignment_length)

    # Step 3: Classify each event as homopolymer or structural
    structural_events = []
    homopolymer_events = []

    for start_col, end_col in indel_events:
        if is_homopolymer_event(aligned_seqs, start_col, end_col):
            homopolymer_events.append((start_col, end_col))
        else:
            structural_events.append((start_col, end_col))

    # Step 4: Calculate statistics
    # Count is number of events (not columns)
    structural_indel_count = len(structural_events)
    homopolymer_indel_count = len(homopolymer_events)

    # Length is the size of the longest event
    structural_indel_length = max((end - start + 1 for start, end in structural_events), default=0)
    homopolymer_indel_length = max((end - start + 1 for start, end in homopolymer_events), default=0)

    # Backward compatibility: total events and max length
    total_indel_count = structural_indel_count + homopolymer_indel_count
    max_indel_length = max(structural_indel_length, homopolymer_indel_length)

    return {
        'snp_count': snp_count,
        'structural_indel_count': structural_indel_count,
        'structural_indel_length': structural_indel_length,
        'homopolymer_indel_count': homopolymer_indel_count,
        'homopolymer_indel_length': homopolymer_indel_length,
        'indel_count': total_indel_count,  # Backward compatibility
        'max_indel_length': max_indel_length  # Backward compatibility
    }


def generate_all_subsets_by_size(variants: List[ConsensusInfo]) -> List[Tuple[int, ...]]:
    """
    Generate all possible non-empty subsets of variant indices.
    Returns subsets in descending order by total cluster size.

    This exhaustive approach guarantees finding the globally optimal merge
    when the number of variants is small (≤ MAX_MSA_MERGE_VARIANTS).

    Args:
        variants: List of variants to generate subsets from

    Returns:
        List of tuples of indices, sorted by total size descending
    """
    n = len(variants)
    sizes = [v.size for v in variants]

    # Build list of (total_size, subset_indices) tuples
    candidates = []

    # Generate all non-empty subsets
    for r in range(n, 0, -1):  # From largest to smallest subset size
        for indices in itertools.combinations(range(n), r):
            total_size = sum(sizes[i] for i in indices)
            candidates.append((total_size, indices))

    # Sort by total size descending
    candidates.sort(reverse=True, key=lambda x: x[0])

    # Return just the subset indices
    return [subset for _, subset in candidates]


def is_compatible_subset(variant_stats: dict, args) -> bool:
    """
    Check if variant statistics are within merge limits.

    By default, homopolymer indels are ignored (treated as compatible) to match
    adjusted-identity homopolymer normalization semantics where AAA ≈ AAAA.
    Only structural indels count against the limits.

    When --disable-homopolymer-equivalence is set, homopolymer indels are treated
    the same as structural indels and count against merge limits.
    """

    # Check SNP limit
    if variant_stats['snp_count'] > 0 and not args.merge_snp:
        return False

    # Determine which indels to count based on homopolymer equivalence setting
    if args.disable_homopolymer_equivalence:
        # Count both structural and homopolymer indels
        indel_count = variant_stats['structural_indel_count'] + variant_stats['homopolymer_indel_count']
        indel_length = max(variant_stats['structural_indel_length'],
                          variant_stats['homopolymer_indel_length'])
    else:
        # Only count structural indels (homopolymer indels ignored)
        indel_count = variant_stats['structural_indel_count']
        indel_length = variant_stats['structural_indel_length']

    # Check indel limits
    if indel_count > 0:
        if args.merge_indel_length == 0:
            return False
        if indel_length > args.merge_indel_length:
            return False

    # Check total position count
    total_positions = variant_stats['snp_count'] + indel_count
    if total_positions > args.merge_position_count:
        return False

    return True


def create_consensus_from_msa(aligned_seqs: List, variants: List[ConsensusInfo]) -> ConsensusInfo:
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
    raw_ric_values = sorted([v.ric for v in variants], reverse=True) if len(variants) > 1 else None

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
        raw_ric=raw_ric_values,
        rid=largest_variant.rid,  # Preserve identity metrics from largest variant
        rid_min=largest_variant.rid_min,
    )


def merge_group_with_msa(variants: List[ConsensusInfo], args) -> Tuple[List[ConsensusInfo], Dict, int]:
    """
    Find largest mergeable subset of variants using MSA-based evaluation with exhaustive search.

    Algorithm:
    1. Process variants in batches of up to MAX_MSA_MERGE_VARIANTS
    2. For each batch, run SPOA MSA once
    3. Exhaustively evaluate ALL subsets by total size (descending)
    4. Merge the best compatible subset found
    5. Remove merged variants and repeat with remaining

    This approach guarantees optimal results when N ≤ MAX_MSA_MERGE_VARIANTS.
    For N > MAX, processes top MAX per round (potentially suboptimal globally).

    Args:
        variants: List of ConsensusInfo from HAC group
        args: Command-line arguments with merge parameters

    Returns:
        (merged_variants, merge_traceability, potentially_suboptimal) where:
        - merged_variants is list of merged ConsensusInfo objects
        - traceability maps merged names to original cluster names
        - potentially_suboptimal is 1 if group had >MAX variants, 0 otherwise
    """
    if len(variants) == 1:
        return variants, {}, 0

    # Track if this group is potentially suboptimal (too many variants for global optimum)
    potentially_suboptimal = 1 if len(variants) > MAX_MSA_MERGE_VARIANTS else 0

    # Sort variants by size (largest first)
    remaining_variants = sorted(variants, key=lambda v: v.size, reverse=True)
    merged_results = []
    all_traceability = {}

    while remaining_variants:
        # Take up to MAX_MSA_MERGE_VARIANTS candidates
        candidates = remaining_variants[:MAX_MSA_MERGE_VARIANTS]

        # Apply size ratio filter if enabled (relative to largest in batch)
        if args.merge_min_size_ratio > 0:
            largest_size = candidates[0].size
            filtered_candidates = [v for v in candidates
                                  if (v.size / largest_size) >= args.merge_min_size_ratio]
            if len(filtered_candidates) < len(candidates):
                filtered_count = len(candidates) - len(filtered_candidates)
                logging.debug(f"Filtered out {filtered_count} variants with size ratio < {args.merge_min_size_ratio} relative to largest (size={largest_size})")
                candidates = filtered_candidates

        # Single candidate - just pass through
        if len(candidates) == 1:
            merged_results.append(candidates[0])
            remaining_variants.remove(candidates[0])
            continue

        logging.debug(f"Evaluating {len(candidates)} variants for merging (exhaustive subset search)")

        # Run SPOA MSA on candidates
        sequences = [v.sequence for v in candidates]
        aligned_seqs = run_spoa_msa(sequences)

        logging.debug(f"Generated MSA with length {len(aligned_seqs[0].seq)}")

        # Generate ALL subsets sorted by total size (exhaustive search)
        all_subsets = generate_all_subsets_by_size(candidates)

        logging.debug(f"Evaluating {len(all_subsets)} candidate subsets")

        # Find first (largest) compatible subset
        merged_this_round = False
        for subset_indices in all_subsets:
            subset_variants = [candidates[i] for i in subset_indices]
            subset_aligned = [aligned_seqs[i] for i in subset_indices]

            # Analyze MSA for this subset
            variant_stats = analyze_msa_columns(subset_aligned)

            # Check compatibility against merge limits
            if is_compatible_subset(variant_stats, args):
                # Only log "mergeable subset" message for actual merges (>1 variant)
                if len(subset_indices) > 1:
                    # Build detailed variant description
                    parts = []
                    if variant_stats['snp_count'] > 0:
                        parts.append(f"{variant_stats['snp_count']} SNPs")
                    if variant_stats['structural_indel_count'] > 0:
                        parts.append(f"{variant_stats['structural_indel_count']} structural indels")
                    if variant_stats['homopolymer_indel_count'] > 0:
                        parts.append(f"{variant_stats['homopolymer_indel_count']} homopolymer indels")

                    variant_desc = ", ".join(parts) if parts else "identical sequences"
                    logging.info(f"Found mergeable subset of {len(subset_indices)} variants: {variant_desc}")

                # Create merged consensus
                merged_consensus = create_consensus_from_msa(
                    subset_aligned, subset_variants
                )

                # Track merge provenance
                traceability = {
                    merged_consensus.sample_name: [v.sample_name for v in subset_variants]
                }
                all_traceability.update(traceability)

                # Add merged consensus to results
                merged_results.append(merged_consensus)

                # Remove merged variants from remaining pool
                for v in subset_variants:
                    if v in remaining_variants:
                        remaining_variants.remove(v)

                merged_this_round = True
                break

        # If no merge found, keep largest variant as-is and continue
        if not merged_this_round:
            logging.debug(f"No compatible merge found for largest variant (size={candidates[0].size})")
            merged_results.append(candidates[0])
            remaining_variants.remove(candidates[0])

    return merged_results, all_traceability, potentially_suboptimal


def bases_match_with_iupac(base1: str, base2: str) -> bool:
    """
    Check if two bases match, considering IUPAC ambiguity codes.
    Two bases match if their IUPAC expansions have any nucleotides in common.
    """
    if base1 == base2:
        return True
    
    # Handle gap characters
    if base1 == '-' or base2 == '-':
        return base1 == base2
    
    # Expand IUPAC codes and check for overlap
    expansion1 = expand_iupac_code(base1)
    expansion2 = expand_iupac_code(base2)
    
    # Bases match if their expansions have any nucleotides in common
    return bool(expansion1.intersection(expansion2))


def expand_iupac_code(base: str) -> set:
    """
    Expand an IUPAC code to its constituent nucleotides.
    Returns a set of nucleotides that the code represents.
    """
    iupac_expansion = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'},
    }
    
    return iupac_expansion.get(base.upper(), {'N'})


def create_variant_summary(primary_seq: str, variant_seq: str) -> str:
    """
    Compare a variant sequence to the primary sequence and create a summary string
    describing the differences. Returns a summary like:
    "3 substitutions, 1 single-nt indel, 1 short (<= 3nt) indel, 2 long indels"
    """
    if not primary_seq or not variant_seq:
        return "sequences empty - cannot compare"
    
    if primary_seq == variant_seq:
        return "identical sequences"
    
    try:
        # Get alignment from edlib with IUPAC awareness
        result = edlib.align(primary_seq, variant_seq, task="path", additionalEqualities=IUPAC_EQUIV)
        if result["editDistance"] == -1:
            return "alignment failed"
        
        # Get nice alignment to examine differences
        alignment = edlib.getNiceAlignment(result, primary_seq, variant_seq)
        if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
            return f"alignment parsing failed - edit distance {result['editDistance']}"
        
        query_aligned = alignment['query_aligned']
        target_aligned = alignment['target_aligned']
        
        # Categorize differences
        substitutions = 0
        single_nt_indels = 0  # Single nucleotide indels
        short_indels = 0      # 2-3 nt indels  
        long_indels = 0       # 4+ nt indels
        
        i = 0
        while i < len(query_aligned):
            query_char = query_aligned[i]
            target_char = target_aligned[i]
            
            # Check if characters are different, considering IUPAC codes
            if not bases_match_with_iupac(query_char, target_char):
                if query_char == '-' or target_char == '-':
                    # This is an indel - determine its length
                    indel_length = 1
                    
                    # Count consecutive indels
                    j = i + 1
                    while j < len(query_aligned) and (query_aligned[j] == '-' or target_aligned[j] == '-'):
                        indel_length += 1
                        j += 1
                    
                    # Categorize by length
                    if indel_length == 1:
                        single_nt_indels += 1
                    elif indel_length <= 3:
                        short_indels += 1
                    else:
                        long_indels += 1
                    
                    # Skip the rest of this indel
                    i = j
                    continue
                else:
                    # This is a substitution
                    substitutions += 1
            
            i += 1
        
        # Build summary string
        parts = []
        if substitutions > 0:
            parts.append(f"{substitutions} substitution{'s' if substitutions != 1 else ''}")
        if single_nt_indels > 0:
            parts.append(f"{single_nt_indels} single-nt indel{'s' if single_nt_indels != 1 else ''}")
        if short_indels > 0:
            parts.append(f"{short_indels} short (<= 3nt) indel{'s' if short_indels != 1 else ''}")
        if long_indels > 0:
            parts.append(f"{long_indels} long indel{'s' if long_indels != 1 else ''}")
        
        if not parts:
            return "identical sequences (IUPAC-compatible)"
        
        return ", ".join(parts)
        
    except Exception as e:
        return f"comparison failed: {str(e)}"


def calculate_adjusted_identity_distance(seq1: str, seq2: str) -> float:
    """Calculate adjusted identity distance between two sequences."""
    if not seq1 or not seq2:
        return 1.0  # Maximum distance

    if seq1 == seq2:
        return 0.0

    # Get alignment from edlib with IUPAC awareness
    result = edlib.align(seq1, seq2, task="path", additionalEqualities=IUPAC_EQUIV)
    if result["editDistance"] == -1:
        return 1.0
        
    # Get nice alignment for adjusted identity scoring
    alignment = edlib.getNiceAlignment(result, seq1, seq2)
    if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
        return 1.0

    # Calculate adjusted identity using standard parameters
    score_result = score_alignment(
        alignment['query_aligned'],
        alignment['target_aligned'],
        adjustment_params=STANDARD_ADJUSTMENT_PARAMS
    )

    # Convert adjusted identity to distance
    return 1.0 - score_result.identity


# ============================================================================
# Alignment Matrix Visualization Functions
# ============================================================================

def analyze_cluster_quality(
    msa_file: str,
    consensus_seq: str,
    max_reads: Optional[int] = None
) -> Optional[ClusterQualityData]:
    """
    Analyze cluster quality using core.py's analyze_positional_variation().

    Uses the canonical positional analysis from core.py to ensure consistent
    treatment of homopolymer length differences across the pipeline.

    Args:
        msa_file: Path to MSA FASTA file
        consensus_seq: Ungapped consensus sequence
        max_reads: Maximum reads to include (for downsampling large clusters)

    Returns:
        ClusterQualityData with position error rates and read identities, or None if failed
    """
    if not os.path.exists(msa_file):
        logging.debug(f"MSA file not found: {msa_file}")
        return None

    # Load MSA file content
    try:
        with open(msa_file, 'r') as f:
            msa_string = f.read()
    except Exception as e:
        logging.debug(f"Failed to read MSA file {msa_file}: {e}")
        return None

    # Extract alignments from MSA using core.py function with homopolymer normalization
    # This returns ReadAlignment objects with score_aligned field
    alignments, msa_consensus, msa_to_consensus_pos = extract_alignments_from_msa(
        msa_string,
        enable_homopolymer_normalization=True
    )

    if not alignments:
        logging.debug(f"No alignments found in MSA: {msa_file}")
        return None

    # Verify consensus matches: the passed-in consensus_seq may be trimmed (shorter) with IUPAC codes
    # The MSA consensus is untrimmed (longer) without IUPAC codes
    # Use edlib in HW mode to check if trimmed consensus is contained within MSA consensus
    if msa_consensus and msa_consensus != consensus_seq:
        # Use edlib HW mode (semi-global) to find consensus_seq within msa_consensus
        # This handles primer trimming (length difference) and IUPAC codes (via equivalencies)
        result = edlib.align(consensus_seq, msa_consensus, mode="HW", task="distance",
                             additionalEqualities=IUPAC_EQUIV)
        edit_distance = result["editDistance"]
        if edit_distance > 0:  # Any edits indicate a real mismatch
            logging.warning(f"Consensus mismatch in MSA file: {msa_file}")
            logging.warning(f"  MSA length: {len(msa_consensus)}, consensus length: {len(consensus_seq)}, edit distance: {edit_distance}")

    # Use the passed-in consensus (with IUPAC codes) as authoritative for quality analysis
    # This reflects the actual output sequence
    consensus_length = len(consensus_seq)

    if consensus_length == 0:
        logging.debug(f"Empty consensus sequence: {msa_file}")
        return None

    # Get consensus aligned sequence by parsing MSA string
    msa_handle = StringIO(msa_string)
    records = list(SeqIO.parse(msa_handle, 'fasta'))
    consensus_aligned = None
    for record in records:
        if 'Consensus' in record.description or 'Consensus' in record.id:
            consensus_aligned = str(record.seq).upper()
            break

    if consensus_aligned is None:
        logging.debug(f"No consensus found in MSA: {msa_file}")
        return None

    # Downsample reads if needed
    if max_reads and len(alignments) > max_reads:
        # Sort by read identity (using normalized edit distance) and take worst reads first, then best
        # This gives us a representative sample showing the quality range
        read_identities_temp = []
        for alignment in alignments:
            # Use normalized edit distance for identity calculation
            identity = 1.0 - (alignment.normalized_edit_distance / consensus_length) if consensus_length > 0 else 0.0
            read_identities_temp.append((identity, alignment))

        # Sort by identity
        read_identities_temp.sort(key=lambda x: x[0])

        # Take worst half and best half
        n_worst = max_reads // 2
        n_best = max_reads - n_worst
        sampled = read_identities_temp[:n_worst] + read_identities_temp[-n_best:]

        alignments = [alignment for _, alignment in sampled]
        logging.debug(f"Downsampled {len(read_identities_temp)} reads to {len(alignments)} for analysis")

    # Use core.py's canonical positional analysis
    position_stats = analyze_positional_variation(alignments, consensus_aligned, msa_to_consensus_pos)

    # Extract position error rates and counts for consensus positions only (skip insertion columns)
    consensus_position_stats = [ps for ps in position_stats if ps.consensus_position is not None]
    # Sort by consensus position to ensure correct order
    consensus_position_stats.sort(key=lambda ps: ps.consensus_position)
    position_error_rates = [ps.error_rate for ps in consensus_position_stats]
    position_error_counts = [ps.error_count for ps in consensus_position_stats]

    # Calculate per-read identities from alignments
    read_identities = []
    for alignment in alignments:
        # Use normalized edit distance for identity calculation
        identity = 1.0 - (alignment.normalized_edit_distance / consensus_length) if consensus_length > 0 else 0.0
        read_identities.append(identity)

    return ClusterQualityData(
        consensus_seq=consensus_seq,
        position_error_rates=position_error_rates,
        position_error_counts=position_error_counts,
        read_identities=read_identities,
        position_stats=consensus_position_stats
    )


def perform_hac_clustering(consensus_list: List[ConsensusInfo], 
                          variant_group_identity: float) -> Dict[int, List[ConsensusInfo]]:
    """
    Perform Hierarchical Agglomerative Clustering using complete linkage.
    Separates specimens from variants based on identity threshold.
    Returns groups of consensus sequences.
    """
    if len(consensus_list) <= 1:
        return {0: consensus_list}

    logging.debug(f"Performing HAC clustering with {variant_group_identity} identity threshold")

    n = len(consensus_list)
    distance_threshold = 1.0 - variant_group_identity
    
    # Initialize each sequence as its own cluster
    clusters = [[i] for i in range(n)]
    
    # Build initial distance matrix between individual sequences
    seq_distances = {}
    for i, j in itertools.combinations(range(n), 2):
        dist = calculate_adjusted_identity_distance(
            consensus_list[i].sequence,
            consensus_list[j].sequence
        )
        seq_distances[(i, j)] = dist
        seq_distances[(j, i)] = dist
    
    def cluster_distance(cluster1: List[int], cluster2: List[int]) -> float:
        """Calculate complete linkage distance between two clusters."""
        max_dist = 0.0
        for i in cluster1:
            for j in cluster2:
                if i < j:
                    dist = seq_distances[(i, j)]
                elif i > j:
                    dist = seq_distances[(j, i)]
                else:
                    dist = 0.0
                max_dist = max(max_dist, dist)
        return max_dist
    
    # Perform hierarchical clustering
    while len(clusters) > 1:
        # Find closest pair of clusters
        min_distance = float('inf')
        merge_pair = None
        
        for i, j in itertools.combinations(range(len(clusters)), 2):
            dist = cluster_distance(clusters[i], clusters[j])
            if dist < min_distance:
                min_distance = dist
                merge_pair = (i, j)
        
        # If minimum distance exceeds threshold, stop clustering
        if min_distance >= distance_threshold:
            break
            
        # Merge the closest clusters
        i, j = merge_pair
        merged_cluster = clusters[i] + clusters[j]
        
        # Remove old clusters (in reverse order to maintain indices)
        new_clusters = []
        for idx, cluster in enumerate(clusters):
            if idx != i and idx != j:
                new_clusters.append(cluster)
        new_clusters.append(merged_cluster)
        clusters = new_clusters
        
        logging.debug(f"Merged clusters with distance {min_distance:.3f}, now have {len(clusters)} clusters")
    
    # Convert clusters to groups of ConsensusInfo
    groups = {}
    for group_id, cluster_indices in enumerate(clusters):
        group_members = [consensus_list[idx] for idx in cluster_indices]
        groups[group_id] = group_members

    logging.debug(f"HAC clustering created {len(groups)} groups")
    for group_id, group_members in groups.items():
        member_names = [m.sample_name for m in group_members]
        # Convert group_id to final naming (group 0 -> 1, group 1 -> 2, etc.)
        final_group_name = group_id + 1
        logging.debug(f"Group {final_group_name}: {member_names}")

    return groups


def select_variants(group: List[ConsensusInfo],
                   max_variants: int,
                   variant_selection: str,
                   group_number: int = None) -> List[ConsensusInfo]:
    """
    Select variants from a group based on the specified strategy.
    Always includes the largest variant first.
    max_variants of 0 or -1 means no limit (return all variants).

    Logs variant summaries for ALL variants in the group, including those
    that will be skipped in the final output.

    Args:
        group: List of ConsensusInfo to select from
        max_variants: Maximum total variants per group (0 or -1 for no limit)
        variant_selection: Selection strategy ("size" or "diversity")
        group_number: Group number for logging prefix (optional)
    """
    # Sort by size, largest first
    sorted_group = sorted(group, key=lambda x: x.size, reverse=True)
    
    if not sorted_group:
        return []
    
    # The primary variant is always the largest
    primary_variant = sorted_group[0]

    # Build prefix for logging
    prefix = f"Group {group_number}: " if group_number is not None else ""

    # Only log Primary when there are other variants to compare against
    if len(sorted_group) > 1:
        logging.info(f"{prefix}Primary: {primary_variant.sample_name} (size={primary_variant.size}, ric={primary_variant.ric})")

    # Handle no limit case (0 or -1 means unlimited)
    if max_variants <= 0:
        selected = sorted_group
    elif len(group) <= max_variants:
        selected = sorted_group
    else:
        # Always include the largest (main) variant
        selected = [primary_variant]
        candidates = sorted_group[1:]

        if variant_selection == "size":
            # Select by size (max_variants - 1 because we already have primary)
            selected.extend(candidates[:max_variants - 1])
        else:  # diversity
            # Select by diversity (maximum distance from already selected)
            while len(selected) < max_variants and candidates:
                best_candidate = None
                best_min_distance = -1
                
                for candidate in candidates:
                    # Calculate minimum distance to all selected variants
                    min_distance = min(
                        calculate_adjusted_identity_distance(candidate.sequence, sel.sequence)
                        for sel in selected
                    )
                    
                    if min_distance > best_min_distance:
                        best_min_distance = min_distance
                        best_candidate = candidate
                
                if best_candidate:
                    selected.append(best_candidate)
                    candidates.remove(best_candidate)
    
    # Now generate variant summaries, showing selected variants first in their final order
    # Then show skipped variants
    
    # Log selected variants first (excluding primary, which is already logged)
    selected_secondary = selected[1:]  # Exclude primary variant
    for i, variant in enumerate(selected_secondary, 1):
        variant_summary = create_variant_summary(primary_variant.sequence, variant.sequence)
        logging.info(f"{prefix}Variant {i}: (size={variant.size}, ric={variant.ric}) - {variant_summary}")

    # Log skipped variants
    selected_names = {variant.sample_name for variant in selected}
    skipped_variants = [v for v in sorted_group[1:] if v.sample_name not in selected_names]

    for i, variant in enumerate(skipped_variants):
        # Calculate what the variant number would have been in the original sorted order
        original_position = next(j for j, v in enumerate(sorted_group) if v.sample_name == variant.sample_name)
        variant_summary = create_variant_summary(primary_variant.sequence, variant.sequence)
        logging.info(f"{prefix}Variant {original_position}: (size={variant.size}, ric={variant.ric}) - {variant_summary} - skipping")
    
    return selected


def create_output_structure(groups: Dict[int, List[ConsensusInfo]], 
                           max_variants: int,
                           variant_selection: str,
                           summary_folder: str) -> Tuple[List[ConsensusInfo], Dict]:
    """
    Create the final output structure with proper naming.
    Returns final consensus list and naming information.
    """
    os.makedirs(summary_folder, exist_ok=True)
    os.makedirs(os.path.join(summary_folder, 'FASTQ Files'), exist_ok=True)
    os.makedirs(os.path.join(summary_folder, 'variants'), exist_ok=True)
    os.makedirs(os.path.join(summary_folder, 'variants', 'FASTQ Files'), exist_ok=True)

    final_consensus = []
    naming_info = {}
    
    # Sort groups by size of largest member (descending)
    sorted_groups = sorted(groups.items(), 
                          key=lambda x: max(m.size for m in x[1]), 
                          reverse=True)
    
    for group_idx, (_, group_members) in enumerate(sorted_groups, 1):
        # Select variants for this group
        selected_variants = select_variants(group_members, max_variants, variant_selection, group_number=group_idx)
        
        # Create naming for this group
        group_naming = []
        
        for variant_idx, variant in enumerate(selected_variants):
            # All variants get .v suffix (primary is .v1, additional are .v2, .v3, etc.)
            new_name = f"{variant.sample_name.split('-c')[0]}-{group_idx}.v{variant_idx + 1}"
            
            # Create new ConsensusInfo with updated name
            renamed_variant = ConsensusInfo(
                sample_name=new_name,
                cluster_id=variant.cluster_id,
                sequence=variant.sequence,
                ric=variant.ric,
                size=variant.size,
                file_path=variant.file_path,
                snp_count=variant.snp_count,  # Preserve SNP count from original
                primers=variant.primers,  # Preserve primers
                raw_ric=variant.raw_ric,  # Preserve raw_ric
                rid=variant.rid,  # Preserve identity metrics
                rid_min=variant.rid_min,
            )

            final_consensus.append(renamed_variant)
            group_naming.append((variant.sample_name, new_name))
        
        naming_info[group_idx] = group_naming
    
    return final_consensus, naming_info


def write_consensus_fastq(consensus: ConsensusInfo,
                         merge_traceability: Dict[str, List[str]],
                         naming_info: Dict,
                         fastq_dir: str,
                         fastq_lookup: Dict[str, List[str]],
                         original_consensus_lookup: Dict[str, ConsensusInfo]):
    """Write FASTQ file for a consensus by concatenating existing FASTQ files."""
    import shutil
    
    # Find the original cluster name(s) by looking through naming_info
    original_clusters = []
    for group_naming in naming_info.values():
        for original_name, final_name in group_naming:
            if final_name == consensus.sample_name:
                # This original cluster contributed to our final consensus
                if original_name in merge_traceability:
                    # This was a merged cluster, get all original contributors
                    original_clusters.extend(merge_traceability[original_name])
                else:
                    # This was not merged, just add it directly
                    original_clusters.append(original_name)
                break
    
    if not original_clusters:
        logging.warning(f"Could not find contributing clusters for {consensus.sample_name}")
        return
    
    # Find FASTQ files for these clusters using lookup table
    fastq_output_path = os.path.join(fastq_dir, f"{consensus.sample_name}-RiC{consensus.ric}.fastq")
    input_files = []
    
    for cluster_name in original_clusters:
        # Look for specimen name from cluster name (e.g., "sample-c1" -> "sample")
        if '-c' in cluster_name:
            specimen_name = cluster_name.rsplit('-c', 1)[0]
            debug_files = fastq_lookup.get(specimen_name, [])

            # Get the original RiC value for this cluster
            original_ric = original_consensus_lookup.get(cluster_name)
            if not original_ric:
                logging.warning(f"Could not find original consensus info for {cluster_name}")
                continue

            # Filter files that match this specific cluster with exact RiC value
            # Match the full pattern: {specimen}-c{cluster}-RiC{exact_ric}-{stage}.fastq
            # This prevents matching multiple RiC values for the same cluster
            cluster_ric_pattern = f"{cluster_name}-RiC{original_ric.ric}-"
            matching_files = [f for f in debug_files if cluster_ric_pattern in f]

            # Validate that matched files exist and log any issues
            for mf in matching_files:
                if not os.path.exists(mf):
                    logging.warning(f"Matched file does not exist: {mf}")
                elif os.path.getsize(mf) == 0:
                    logging.warning(f"Matched file is empty: {mf}")

            input_files.extend(matching_files)

    if not input_files:
        logging.warning(f"No FASTQ files found for {consensus.sample_name} from clusters: {original_clusters}")
        return
    
    # Concatenate files directly without parsing
    files_processed = 0
    try:
        with open(fastq_output_path, 'wb') as outf:
            for input_file in input_files:
                try:
                    file_size = os.path.getsize(input_file)
                    if file_size > 0:
                        with open(input_file, 'rb') as inf:
                            shutil.copyfileobj(inf, outf)
                        files_processed += 1
                    else:
                        logging.debug(f"Skipping empty file: {input_file}")
                except Exception as e:
                    logging.debug(f"Could not concatenate {input_file}: {e}")
        
        # Check if the output file has content
        output_size = os.path.getsize(fastq_output_path)
        if output_size > 0:
            # Count reads for logging by quickly counting lines and dividing by 4
            with open(fastq_output_path, 'r') as f:
                line_count = sum(1 for line in f)
            read_count = line_count // 4
            logging.debug(f"Concatenated {files_processed}/{len(input_files)} files ({output_size:,} bytes) with ~{read_count} reads to {fastq_output_path}")
        else:
            # Debug: check what files were supposed to be concatenated
            file_info = []
            for input_file in input_files:
                if os.path.exists(input_file):
                    size = os.path.getsize(input_file)
                    file_info.append(f"{os.path.basename(input_file)}:{size}B")
                else:
                    file_info.append(f"{os.path.basename(input_file)}:missing")
            
            logging.warning(f"No data written for {consensus.sample_name} - input files: {', '.join(file_info)}")
            # Remove empty output file
            try:
                os.unlink(fastq_output_path)
            except:
                pass
            
    except Exception as e:
        logging.error(f"Failed to write concatenated FASTQ file {fastq_output_path}: {e}")


def write_specimen_data_files(specimen_consensus: List[ConsensusInfo],
                               merge_traceability: Dict[str, List[str]],
                               naming_info: Dict,
                               summary_folder: str,
                               fastq_dir: str,
                               fastq_lookup: Dict[str, List[str]],
                               original_consensus_lookup: Dict[str, ConsensusInfo],
                               fasta_fields: List[FastaField]
                               ) -> List[Tuple[ConsensusInfo, str]]:
    """
    Write individual FASTA and FASTQ files for a single specimen.
    Does NOT write summary files (summary.fasta, summary.txt).

    Args:
        fasta_fields: List of FastaField objects defining header format

    Returns:
        List of (raw_consensus, original_cluster_name) tuples for later use in summary.fasta
    """
    # Generate .raw file consensuses for merged variants
    raw_file_consensuses = []
    for consensus in specimen_consensus:
        # Only create .raw files if this consensus was actually merged
        if consensus.raw_ric and len(consensus.raw_ric) > 1:
            # Find the original cluster name from naming_info
            original_cluster_name = None
            for group_naming in naming_info.values():
                for orig_name, final_name in group_naming:
                    if final_name == consensus.sample_name:
                        original_cluster_name = orig_name
                        break
                if original_cluster_name:
                    break

            # Get contributing clusters from merge_traceability
            if original_cluster_name and original_cluster_name in merge_traceability:
                contributing_clusters = merge_traceability[original_cluster_name]

                # Sort by size (descending) to match .raw1, .raw2 ordering
                contributing_infos = []
                for cluster_name in contributing_clusters:
                    if cluster_name in original_consensus_lookup:
                        contributing_infos.append(original_consensus_lookup[cluster_name])

                contributing_infos.sort(key=lambda x: x.size, reverse=True)

                # Create .raw file entries
                for raw_idx, raw_info in enumerate(contributing_infos, 1):
                    raw_name = f"{consensus.sample_name}.raw{raw_idx}"

                    # Create new ConsensusInfo with .raw name but original sequence/metadata
                    raw_consensus = ConsensusInfo(
                        sample_name=raw_name,
                        cluster_id=raw_info.cluster_id,
                        sequence=raw_info.sequence,
                        ric=raw_info.ric,
                        size=raw_info.size,
                        file_path=raw_info.file_path,
                        snp_count=None,  # Pre-merge, no SNPs from merging
                        primers=raw_info.primers,
                        raw_ric=None,  # Pre-merge, not merged
                        rid=raw_info.rid,  # Preserve read identity metrics
                        rid_min=raw_info.rid_min,
                    )
                    raw_file_consensuses.append((raw_consensus, raw_info.sample_name))

    # Write individual FASTA files with custom field formatting
    for consensus in specimen_consensus:
        output_file = os.path.join(summary_folder, f"{consensus.sample_name}-RiC{consensus.ric}.fasta")
        with open(output_file, 'w') as f:
            header = format_fasta_header(consensus, fasta_fields)
            f.write(f">{header}\n")
            f.write(f"{consensus.sequence}\n")

    # Write FASTQ files for each final consensus containing all contributing reads
    for consensus in specimen_consensus:
        write_consensus_fastq(consensus, merge_traceability, naming_info, fastq_dir, fastq_lookup, original_consensus_lookup)

    # Write .raw files (individual FASTA and FASTQ for pre-merge variants)
    for raw_consensus, original_cluster_name in raw_file_consensuses:
        # Write individual FASTA file with custom field formatting
        output_file = os.path.join(summary_folder, 'variants', f"{raw_consensus.sample_name}-RiC{raw_consensus.ric}.fasta")
        with open(output_file, 'w') as f:
            header = format_fasta_header(raw_consensus, fasta_fields)
            f.write(f">{header}\n")
            f.write(f"{raw_consensus.sequence}\n")

        # Write FASTQ file by finding the original cluster's FASTQ
        # Look for specimen name from original cluster name
        if '-c' in original_cluster_name:
            specimen_name = original_cluster_name.rsplit('-c', 1)[0]
            debug_files = fastq_lookup.get(specimen_name, []) if fastq_lookup else []

            # Filter files that match this specific cluster with exact RiC value
            # Use the raw_consensus.ric which came from the original cluster
            cluster_ric_pattern = f"{original_cluster_name}-RiC{raw_consensus.ric}-"
            matching_files = [f for f in debug_files if cluster_ric_pattern in f]

            if matching_files:
                fastq_output_path = os.path.join(summary_folder, 'variants', 'FASTQ Files', f"{raw_consensus.sample_name}-RiC{raw_consensus.ric}.fastq")
                try:
                    with open(fastq_output_path, 'wb') as outf:
                        for input_file in matching_files:
                            if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
                                with open(input_file, 'rb') as inf:
                                    shutil.copyfileobj(inf, outf)
                    logging.debug(f"Wrote .raw FASTQ: {os.path.basename(fastq_output_path)}")
                except Exception as e:
                    logging.debug(f"Could not write .raw FASTQ for {raw_consensus.sample_name}: {e}")

    return raw_file_consensuses


def build_fastq_lookup_table(source_dir: str = ".") -> Dict[str, List[str]]:
    """
    Build a lookup table mapping specimen base names to their cluster FASTQ files.
    This avoids repeated directory scanning during file copying.
    """
    import glob
    import re

    lookup = defaultdict(list)

    # Initialize variables before conditional block
    debug_files = []
    selected_stage = None

    # Scan cluster_debug directory once to build lookup table
    cluster_debug_path = os.path.join(source_dir, "cluster_debug")
    if os.path.exists(cluster_debug_path):
        # Define priority order for stage types (first match wins)
        # This prevents including multiple versions of the same cluster
        stage_priority = ['sampled', 'reads', 'untrimmed']

        # Try each stage type in priority order until we find files
        for stage in stage_priority:
            debug_files = glob.glob(os.path.join(cluster_debug_path, f"*-{stage}.fastq"))
            if debug_files:
                selected_stage = stage
                break

        # If no files found with known stage types, try generic pattern
        if not debug_files:
            debug_files = glob.glob(os.path.join(cluster_debug_path, "*.fastq"))
            selected_stage = "unknown"

        # Use regex to robustly parse the filename pattern
        # Pattern: {specimen}-c{cluster}-RiC{size}-{stage}.fastq
        # Where stage can be: sampled, reads, untrimmed, or other variants
        pattern = re.compile(r'^(.+)-c(\d+)-RiC(\d+)-([a-z]+)\.fastq$')

        for fastq_path in debug_files:
            filename = os.path.basename(fastq_path)
            match = pattern.match(filename)
            if match:
                specimen_name = match.group(1)  # Extract specimen name
                # cluster_num = match.group(2)  # Available if needed
                # ric_value = match.group(3)    # Available if needed
                # stage = match.group(4)        # Stage: sampled, reads, untrimmed, etc.
                lookup[specimen_name].append(fastq_path)
            else:
                logging.warning(f"Skipping file with unexpected name pattern: {filename}")

    if debug_files:
        logging.debug(f"Built FASTQ lookup table for {len(lookup)} specimens with {sum(len(files) for files in lookup.values())} {selected_stage} files")
    else:
        logging.debug("No FASTQ files found in cluster_debug directory")
    return dict(lookup)


def write_position_debug_file(
    sequences_with_pos_outliers: List[Tuple],
    summary_folder: str,
    threshold: float
):
    """Write detailed debug information about high-error positions.

    Creates a separate file with per-position base composition and error details
    to help validate positional phasing and quality analysis.

    Args:
        sequences_with_pos_outliers: List of (ConsensusInfo, result_dict) tuples
        summary_folder: Output directory for the debug file
        threshold: Error rate threshold used for flagging positions
    """
    debug_path = os.path.join(summary_folder, 'position_errors_debug.txt')

    with open(debug_path, 'w') as f:
        f.write("POSITION ERROR DETAILED DEBUG REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Threshold: {threshold:.1%} (positions with error rate above this are flagged)\n")
        f.write(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        if not sequences_with_pos_outliers:
            f.write("No sequences with high-error positions found.\n")
            return

        # Sort by total nucleotide errors descending
        sorted_seqs = sorted(
            sequences_with_pos_outliers,
            key=lambda x: x[1].get('total_nucleotide_errors', 0),
            reverse=True
        )

        for cons, result in sorted_seqs:
            # Handle merged sequences (component_name in result)
            if 'component_name' in result:
                display_name = f"{cons.sample_name} (component: {result['component_name']})"
                ric = result.get('component_ric', cons.ric)
            else:
                display_name = cons.sample_name
                ric = result.get('ric', cons.ric)

            f.write("=" * 80 + "\n")
            f.write(f"SEQUENCE: {display_name}\n")
            f.write(f"RiC: {ric}\n")
            f.write(f"High-error positions: {result['num_outlier_positions']}\n")
            f.write(f"Mean error rate at flagged positions: {result['mean_outlier_error_rate']:.1%}\n")
            f.write(f"Total nucleotide errors: {result['total_nucleotide_errors']}\n")
            f.write("-" * 80 + "\n\n")

            outlier_details = result.get('outlier_details', [])
            if not outlier_details:
                # Fall back to basic info if detailed stats not available
                for pos, rate, count in result.get('outlier_positions', []):
                    f.write(f"  Position {pos+1}: error_rate={rate:.1%}, error_count={count}\n")
                f.write("\n")
                continue

            for detail in outlier_details:
                cons_pos = detail['consensus_position']
                msa_pos = detail.get('msa_position')
                # Display as 1-indexed for user-friendliness
                cons_pos_display = cons_pos + 1 if cons_pos is not None else "?"
                msa_pos_display = msa_pos + 1 if msa_pos is not None else "?"

                f.write(f"Position {cons_pos_display} (MSA: {msa_pos_display}):\n")
                f.write(f"  Consensus base: {detail['consensus_nucleotide']}\n")
                f.write(f"  Coverage: {detail['coverage']}\n")
                f.write(f"  Error rate: {detail['error_rate']:.1%}\n")
                f.write(f"  Error count: {detail['error_count']}\n")
                f.write(f"  Substitutions: {detail['sub_count']}, Insertions: {detail['ins_count']}, Deletions: {detail['del_count']}\n")

                # Format base composition (raw counts from MSA)
                base_comp = detail['base_composition']
                hp_comp = detail.get('homopolymer_composition', {})

                if base_comp:
                    total = sum(base_comp.values())
                    comp_str = ", ".join(
                        f"{base}:{count}({count/total*100:.0f}%)"
                        for base, count in sorted(base_comp.items(), key=lambda x: -x[1])
                        if count > 0
                    )
                    f.write(f"  Raw base composition: {comp_str}\n")

                # Format homopolymer composition if present
                if hp_comp and any(v > 0 for v in hp_comp.values()):
                    hp_str = ", ".join(
                        f"{base}:{count}"
                        for base, count in sorted(hp_comp.items(), key=lambda x: -x[1])
                        if count > 0
                    )
                    f.write(f"  Homopolymer length variants: {hp_str}\n")

                    # Calculate and show effective composition (raw - HP adjustments)
                    # HP variants are normalized away in error calculation
                    if base_comp:
                        effective_comp = {}
                        for base in base_comp:
                            raw = base_comp.get(base, 0)
                            hp_adj = hp_comp.get(base, 0)
                            effective = raw - hp_adj
                            if effective > 0:
                                effective_comp[base] = effective

                        if effective_comp:
                            eff_total = sum(effective_comp.values())
                            eff_str = ", ".join(
                                f"{base}:{count}({count/eff_total*100:.0f}%)"
                                for base, count in sorted(effective_comp.items(), key=lambda x: -x[1])
                                if count > 0
                            )
                            f.write(f"  Effective composition (HP-normalized): {eff_str}\n")

                f.write("\n")

            # Show context: consensus sequence around flagged positions
            consensus_seq = result.get('consensus_seq', '')
            if consensus_seq and outlier_details:
                f.write("Consensus sequence context (flagged positions marked with *):\n")
                # Mark positions in the sequence
                marked_positions = set()
                for detail in outlier_details:
                    if detail['consensus_position'] is not None:
                        marked_positions.add(detail['consensus_position'])

                # Show sequence in chunks of 60 with position markers
                chunk_size = 60
                for chunk_start in range(0, len(consensus_seq), chunk_size):
                    chunk_end = min(chunk_start + chunk_size, len(consensus_seq))
                    chunk = consensus_seq[chunk_start:chunk_end]

                    # Position line
                    f.write(f"  {chunk_start+1:>5}  ")
                    f.write(chunk)
                    f.write(f"  {chunk_end}\n")

                    # Marker line
                    f.write("         ")
                    for i in range(chunk_start, chunk_end):
                        if i in marked_positions:
                            f.write("*")
                        else:
                            f.write(" ")
                    f.write("\n")

                f.write("\n")

    logging.info(f"Position error debug file written to: {debug_path}")


def write_quality_report(final_consensus: List[ConsensusInfo],
                        all_raw_consensuses: List[Tuple[ConsensusInfo, str]],
                        summary_folder: str,
                        source_folder: str):
    """
    Write quality report with rid-based dual outlier detection.

    Uses mean read identity (rid) for outlier detection. rid_min is not used
    because single outlier reads don't significantly impact consensus quality;
    positional analysis better captures systematic issues.

    Structure:
    1. Executive Summary - High-level overview with attention flags
    2. Read Identity Analysis - Dual outlier detection (clustering threshold + statistical)
    3. Positional Identity Analysis - Sequences with problematic positions
    4. Variant Detection - Clusters with detected variants
    5. Merged Sequence Analysis - Quality issues in merged components
    6. Interpretation Guide - Actionable guidance with neutral tone

    Args:
        final_consensus: List of final consensus sequences
        all_raw_consensuses: List of (raw_consensus, original_name) tuples
        summary_folder: Output directory for report
        source_folder: Source directory containing cluster_debug with MSA files
    """
    from datetime import datetime

    quality_report_path = os.path.join(summary_folder, 'quality_report.txt')

    # Build .raw lookup: map merged sequence names to their .raw components
    raw_lookup = {}
    for raw_cons, original_name in all_raw_consensuses:
        base_match = re.match(r'(.+?)\.raw\d+$', raw_cons.sample_name)
        if base_match:
            base_name = base_match.group(1)
            if base_name not in raw_lookup:
                raw_lookup[base_name] = []
            raw_lookup[base_name].append(raw_cons)

    # Identify outliers using dual detection
    outlier_results = identify_outliers(final_consensus, all_raw_consensuses, source_folder)

    # Load min_variant_frequency and min_variant_count from metadata
    # These parameters are used as thresholds for positional outlier detection
    min_variant_frequency = None
    min_variant_count = None

    # Try to load from first available consensus sequence metadata
    for cons in final_consensus:
        # Extract specimen base name (metadata files are named by specimen, not cluster)
        sample_name = cons.sample_name
        specimen_base = re.sub(r'-\d+\.v\d+$', '', sample_name)

        metadata = load_metadata_from_json(source_folder, specimen_base)
        if metadata and 'parameters' in metadata:
            params = metadata['parameters']
            min_variant_frequency = params.get('min_variant_frequency', 0.2)
            min_variant_count = params.get('min_variant_count', 5)
            break

    # Fallback to defaults if not found
    if min_variant_frequency is None:
        min_variant_frequency = 0.2
        logging.warning("Could not load min_variant_frequency from metadata, using default: 0.2")
    if min_variant_count is None:
        min_variant_count = 5
        logging.warning("Could not load min_variant_count from metadata, using default: 5")

    # Analyze positional identity for all sequences
    # For merged sequences, analyze their .raw components instead (which have MSA files)
    pos_identity_results = {}
    sequences_with_pos_outliers = []

    # Collect all sequences to analyze
    # Use dict to deduplicate by sample_name (ConsensusInfo is unhashable due to list fields)
    sequences_to_analyze = {}
    for cons in final_consensus:
        sequences_to_analyze[cons.sample_name] = cons

    # For merged sequences, analyze their .raw components (which have MSA files)
    # For unmerged sequences, analyze them directly
    logging.info("Analyzing positional identity for quality report...")
    for cons in tqdm(sequences_to_analyze.values(), desc="Analyzing positional identity", unit="seq"):
        is_merged = cons.snp_count is not None and cons.snp_count > 0

        if is_merged:
            # Analyze the raw components for this merged sequence
            raw_components = raw_lookup.get(cons.sample_name, [])
            worst_result = None
            worst_outliers = 0

            for raw_cons in raw_components:
                result = analyze_positional_identity_outliers(
                    raw_cons, source_folder, min_variant_frequency, min_variant_count
                )
                if result:
                    result['component_name'] = raw_cons.sample_name
                    result['component_ric'] = raw_cons.ric
                    if result['num_outlier_positions'] > worst_outliers:
                        worst_outliers = result['num_outlier_positions']
                        worst_result = result

            # Report the worst component for this merged sequence
            if worst_result and worst_result['num_outlier_positions'] > 0:
                sequences_with_pos_outliers.append((cons, worst_result))
        else:
            # Unmerged sequence - analyze directly
            result = analyze_positional_identity_outliers(
                cons, source_folder, min_variant_frequency, min_variant_count
            )
            if result and result['num_outlier_positions'] > 0:
                sequences_with_pos_outliers.append((cons, result))

    # Sort positional outliers by total nucleotide errors (desc)
    sequences_with_pos_outliers.sort(key=lambda x: x[1].get('total_nucleotide_errors', 0), reverse=True)

    # Write detailed position debug file
    if sequences_with_pos_outliers:
        write_position_debug_file(sequences_with_pos_outliers, summary_folder, min_variant_frequency)

    # Identify merged sequences with quality issues in components
    merged_with_issues = []
    threshold_rid = outlier_results['global_stats']['stat_threshold_rid']

    for cons in final_consensus:
        is_merged = cons.snp_count is not None and cons.snp_count > 0
        if is_merged:
            raw_components = raw_lookup.get(cons.sample_name, [])
            if not raw_components:
                continue

            # Collect all component info and calculate weighted average rid
            components_info = []
            worst_rid = 1.0
            total_ric = 0
            weighted_rid_sum = 0

            for raw in raw_components:
                rid = raw.rid if raw.rid is not None else 1.0
                ric = raw.ric if raw.ric else 0
                components_info.append((raw, rid, ric))
                if rid < worst_rid:
                    worst_rid = rid
                if ric > 0:
                    total_ric += ric
                    weighted_rid_sum += rid * ric

            # Calculate weighted average rid
            weighted_avg_rid = weighted_rid_sum / total_ric if total_ric > 0 else 1.0

            # Sort components by rid ascending
            components_info.sort(key=lambda x: x[1])

            # Flag if worst component has concerning quality
            if worst_rid < threshold_rid:
                merged_with_issues.append((cons, components_info, worst_rid, weighted_avg_rid))

    # Sort merged issues by worst component quality
    merged_with_issues.sort(key=lambda x: x[2])

    # Write the report
    with open(quality_report_path, 'w') as f:
        # ====================================================================
        # SECTION 1: HEADER
        # ====================================================================
        f.write("=" * 80 + "\n")
        f.write("QUALITY REPORT - speconsense-summarize\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Source: {source_folder}\n")
        f.write("=" * 80 + "\n\n")

        # ====================================================================
        # SECTION 2: EXECUTIVE SUMMARY
        # ====================================================================
        f.write("EXECUTIVE SUMMARY\n")
        f.write("-" * 80 + "\n\n")

        total_seqs = len(final_consensus)
        total_merged = sum(1 for cons in final_consensus
                          if cons.snp_count is not None and cons.snp_count > 0)

        f.write(f"Total sequences: {total_seqs}\n")
        f.write(f"Merged sequences: {total_merged} ({100*total_merged/total_seqs:.1f}%)\n\n")

        # Global statistics
        stats = outlier_results['global_stats']
        f.write(f"Global Read Identity Statistics:\n")
        f.write(f"  Mean rid: {stats['mean_rid']:.1%} ± {stats['std_rid']:.1%}\n\n")

        # Sequences requiring attention
        f.write("Sequences Requiring Attention:\n")

        n_stat = len(outlier_results['statistical_outliers'])
        n_pos = len(sequences_with_pos_outliers)
        n_merged = len(merged_with_issues)

        # Count unique sequences flagged (use sample_name for deduplication)
        flagged_names = set()
        for c, _ in outlier_results['statistical_outliers']:
            flagged_names.add(c.sample_name)
        for c, _ in sequences_with_pos_outliers:
            flagged_names.add(c.sample_name)
        for c, _, _, _ in merged_with_issues:
            flagged_names.add(c.sample_name)
        total_flagged = len(flagged_names)

        # Count combined read identity issues (statistical outliers + merged with component issues)
        n_rid_issues = n_stat + n_merged

        f.write(f"  Total flagged: {total_flagged} ({100*total_flagged/total_seqs:.1f}%)\n")
        f.write(f"    - Low read identity: {n_rid_issues} ({n_stat} sequences + {n_merged} merged)\n")
        f.write(f"    - High-error positions: {n_pos}\n\n")

        # ====================================================================
        # SECTION 3: READ IDENTITY ANALYSIS
        # ====================================================================
        if n_rid_issues > 0:
            f.write("=" * 80 + "\n")
            f.write("READ IDENTITY ANALYSIS\n")
            f.write("=" * 80 + "\n\n")

            f.write("Sequences with mean read identity (rid) below mean - 2×std.\n")
            f.write(f"Threshold: {stats['stat_threshold_rid']:.1%}\n\n")

            f.write(f"{'Sequence':<50} {'RiC':<6} {'rid':<8}\n")
            f.write("-" * 64 + "\n")

            # Build combined list: non-merged outliers + merged with issues
            # Each entry is: (display_rid_for_sorting, is_merged, data)
            # For non-merged: data = (cons, rid)
            # For merged: data = (cons, components_info, worst_rid, weighted_avg_rid)
            combined_entries = []

            # Add non-merged statistical outliers
            for cons, rid in outlier_results['statistical_outliers']:
                is_merged = cons.snp_count is not None and cons.snp_count > 0
                if not is_merged:
                    combined_entries.append((rid, False, (cons, rid)))

            # Add merged sequences with issues
            for entry in merged_with_issues:
                cons, components_info, worst_rid, weighted_avg_rid = entry
                combined_entries.append((worst_rid, True, entry))

            # Sort by rid ascending (using worst component rid for merged)
            combined_entries.sort(key=lambda x: x[0])

            # Display entries
            for _, is_merged, data in combined_entries:
                if is_merged:
                    cons, components_info, worst_rid, weighted_avg_rid = data
                    # Parent row with [merged] tag and weighted average rid
                    name = cons.sample_name
                    name_with_tag = f"{name} [merged]"
                    name_truncated = name_with_tag[:49] if len(name_with_tag) > 49 else name_with_tag
                    rid_str = f"{weighted_avg_rid:.1%}"
                    f.write(f"{name_truncated:<50} {cons.ric:<6} {rid_str:<8}\n")

                    # Component rows (indented)
                    for raw, comp_rid, comp_ric in components_info:
                        # Extract cluster number from component name (e.g., "sample-c1" -> "c1")
                        comp_name = raw.sample_name
                        comp_display = f"  └─ {comp_name}"
                        comp_truncated = comp_display[:49] if len(comp_display) > 49 else comp_display
                        comp_rid_str = f"{comp_rid:.1%}"
                        f.write(f"{comp_truncated:<50} {comp_ric:<6} {comp_rid_str:<8}\n")
                else:
                    cons, rid = data
                    name_truncated = cons.sample_name[:49] if len(cons.sample_name) > 49 else cons.sample_name
                    rid_str = f"{rid:.1%}"
                    f.write(f"{name_truncated:<50} {cons.ric:<6} {rid_str:<8}\n")
            f.write("\n")

        # ====================================================================
        # SECTION 4: POSITIONAL IDENTITY ANALYSIS
        # ====================================================================
        if n_pos > 0:
            f.write("=" * 80 + "\n")
            f.write("POSITIONAL IDENTITY ANALYSIS\n")
            f.write("=" * 80 + "\n\n")

            f.write("Sequences with high-error positions (error rate > threshold at specific positions):\n")
            f.write(f"Threshold: {min_variant_frequency:.1%} (--min-variant-frequency from metadata)\n")
            f.write(f"Min RiC: {2 * min_variant_count} (2 × --min-variant-count)\n")
            f.write("Positions above threshold may indicate undetected/unphased variants.\n")
            f.write("For merged sequences, shows worst component.\n\n")

            # Sort by total nucleotide errors (descending)
            sorted_pos_outliers = sorted(
                sequences_with_pos_outliers,
                key=lambda x: x[1].get('total_nucleotide_errors', 0),
                reverse=True
            )

            # Calculate display names and find max length for dynamic column width
            display_data = []
            for cons, result in sorted_pos_outliers:
                if 'component_name' in result:
                    component_suffix = result['component_name'].split('.')[-1] if '.' in result['component_name'] else ''
                    display_name = f"{cons.sample_name} ({component_suffix})"
                    ric_val = result.get('component_ric', cons.ric)
                else:
                    display_name = cons.sample_name
                    ric_val = cons.ric
                display_data.append((display_name, ric_val, result))

            # Calculate column width based on longest name (minimum 40, cap at 70)
            max_name_len = max(len(name) for name, _, _ in display_data) if display_data else 40
            name_col_width = min(max(max_name_len + 2, 40), 70)

            f.write(f"{'Sequence':<{name_col_width}} {'RiC':<6} {'#Pos':<6} {'MeanErr':<8} {'TotalErr':<10}\n")
            f.write("-" * (name_col_width + 32) + "\n")

            for display_name, ric_val, result in display_data:
                mean_err = result.get('mean_outlier_error_rate', 0.0)
                total_err = result.get('total_nucleotide_errors', 0)
                num_pos = result['num_outlier_positions']

                f.write(f"{display_name:<{name_col_width}} {ric_val:<6} {num_pos:<6} "
                       f"{mean_err:<8.1%} {total_err:<10}\n")

            f.write("\n")

        # ====================================================================
        # SECTION 5: INTERPRETATION GUIDE
        # ====================================================================
        f.write("=" * 80 + "\n")
        f.write("INTERPRETATION GUIDE\n")
        f.write("=" * 80 + "\n\n")

        f.write("Read Identity Analysis:\n")
        f.write("-" * 40 + "\n")
        f.write("  Threshold: mean - 2×std (statistical outliers)\n")
        f.write("  RiC: Read-in-Cluster count\n")
        f.write("  rid: Mean read identity to consensus\n")
        f.write("  [merged]: Weighted average rid; components shown below\n\n")

        f.write("Positional Identity Analysis:\n")
        f.write("-" * 40 + "\n")
        f.write("  Threshold: --min-variant-frequency from metadata\n")
        f.write("  Min RiC: 2 × --min-variant-count\n")
        f.write("  #Pos: Count of positions exceeding error threshold\n")
        f.write("  MeanErr: Average error rate at flagged positions\n")
        f.write("  TotalErr: Sum of errors at flagged positions\n\n")


    # Log summary
    logging.info(f"Quality report written to: {quality_report_path}")
    if n_rid_issues > 0:
        logging.info(f"  {n_rid_issues} sequence(s) flagged for read identity ({n_stat} direct + {n_merged} merged)")
    else:
        logging.info("  All sequences show good read identity")

    if n_pos > 0:
        logging.info(f"  {n_pos} sequence(s) with high-error positions")
    if n_merged > 0:
        logging.info(f"  {n_merged} merged sequence(s) with component quality issues")



def write_output_files(final_consensus: List[ConsensusInfo],
                      all_raw_consensuses: List[Tuple[ConsensusInfo, str]],
                      summary_folder: str,
                      temp_log_file: str,
                      fasta_fields: List[FastaField]):
    """
    Write summary files only. Individual data files already written per-specimen.

    Args:
        fasta_fields: List of FastaField objects defining header format

    Writes:
    - summary.fasta: Combined index of all sequences
    - summary.txt: Statistics and totals
    - summarize_log.txt: Copy of processing log
    """

    # Write combined summary.fasta with custom field formatting
    # Include only final consensus sequences (not .raw pre-merge variants)
    summary_fasta_path = os.path.join(summary_folder, 'summary.fasta')
    with open(summary_fasta_path, 'w') as f:
        # Write final consensus sequences
        for consensus in final_consensus:
            header = format_fasta_header(consensus, fasta_fields)
            f.write(f">{header}\n")
            f.write(f"{consensus.sequence}\n")
    
    # Write summary statistics
    summary_txt_path = os.path.join(summary_folder, 'summary.txt')
    with open(summary_txt_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['Filename', 'Length', 'Reads in Consensus', 'Multiple'])
        
        unique_samples = set()
        total_ric = 0
        specimen_counters = {}
        
        for consensus in final_consensus:
            base_name = consensus.sample_name.split('-')[0]
            
            # Initialize counter for new specimen
            if base_name not in specimen_counters:
                specimen_counters[base_name] = 1
            else:
                specimen_counters[base_name] += 1
            
            multiple_id = specimen_counters[base_name]
            writer.writerow([consensus.sample_name, len(consensus.sequence), consensus.ric, multiple_id])
            unique_samples.add(base_name)
            total_ric += consensus.ric
        
        writer.writerow([])
        writer.writerow(['Total Unique Samples', len(unique_samples)])
        writer.writerow(['Total Consensus Sequences', len(final_consensus)])
        writer.writerow(['Total Reads in Consensus Sequences', total_ric])
    
    # Copy log file to summary directory as summarize_log.txt
    if temp_log_file:
        import shutil
        summarize_log_path = os.path.join(summary_folder, 'summarize_log.txt')
        try:
            # Flush any remaining log entries before copying
            logging.getLogger().handlers[1].flush() if len(logging.getLogger().handlers) > 1 else None
            shutil.copy2(temp_log_file, summarize_log_path)
            logging.info(f"Created log file: {summarize_log_path}")
        except Exception as e:
            logging.warning(f"Could not copy log file: {e}")


def process_single_specimen(file_consensuses: List[ConsensusInfo],
                           args) -> Tuple[List[ConsensusInfo], Dict[str, List[str]], Dict, int]:
    """
    Process a single specimen file: HAC cluster, MSA-based merge per group, and select final variants.
    Returns final consensus list, merge traceability, naming info, and limited_count for this specimen.

    Architecture (Phase 3):
    1. HAC clustering to separate variant groups (primary vs contaminants)
    2. MSA-based merging within each group
    3. Select representative variants per group
    """
    if not file_consensuses:
        return [], {}, {}, 0

    file_name = os.path.basename(file_consensuses[0].file_path)
    logging.info(f"Processing specimen from file: {file_name}")

    # Phase 1: HAC clustering to separate variant groups (moved before merging!)
    variant_groups = perform_hac_clustering(file_consensuses, args.group_identity)

    # Filter to max groups if specified
    if args.select_max_groups > 0 and len(variant_groups) > args.select_max_groups:
        # Sort groups by size of largest member
        sorted_for_filtering = sorted(
            variant_groups.items(),
            key=lambda x: max(m.size for m in x[1]),
            reverse=True
        )
        # Keep only top N groups
        variant_groups = dict(sorted_for_filtering[:args.select_max_groups])
        logging.info(f"Filtered to top {args.select_max_groups} groups by size (from {len(sorted_for_filtering)} total groups)")

    # Phase 2: MSA-based merging within each group
    merged_groups = {}
    all_merge_traceability = {}
    total_limited_count = 0

    for group_id, group_members in variant_groups.items():
        merged, traceability, limited_count = merge_group_with_msa(group_members, args)
        merged_groups[group_id] = merged
        all_merge_traceability.update(traceability)
        total_limited_count += limited_count

    # Phase 3: Select representative variants for each group in this specimen
    final_consensus = []
    naming_info = {}

    # Sort variant groups by size of largest member (descending)
    sorted_groups = sorted(merged_groups.items(),
                          key=lambda x: max(m.size for m in x[1]),
                          reverse=True)

    for group_idx, (_, group_members) in enumerate(sorted_groups):
        final_group_name = group_idx + 1

        # Select variants for this group
        selected_variants = select_variants(group_members, args.select_max_variants, args.select_strategy, group_number=final_group_name)

        # Create naming for this group within this specimen
        group_naming = []

        for variant_idx, variant in enumerate(selected_variants):
            # All variants get .v suffix (primary is .v1, additional are .v2, .v3, etc.)
            new_name = f"{variant.sample_name.split('-c')[0]}-{group_idx + 1}.v{variant_idx + 1}"

            # Create new ConsensusInfo with updated name
            renamed_variant = ConsensusInfo(
                sample_name=new_name,
                cluster_id=variant.cluster_id,
                sequence=variant.sequence,
                ric=variant.ric,
                size=variant.size,
                file_path=variant.file_path,
                snp_count=variant.snp_count,  # Preserve SNP count from merging
                primers=variant.primers,  # Preserve primers
                raw_ric=variant.raw_ric,  # Preserve raw_ric
                rid=variant.rid,  # Preserve identity metrics
                rid_min=variant.rid_min,
            )

            final_consensus.append(renamed_variant)
            group_naming.append((variant.sample_name, new_name))

        naming_info[group_idx + 1] = group_naming

    logging.info(f"Processed {file_name}: {len(final_consensus)} final variants across {len(merged_groups)} groups")
    logging.info("")  # Empty line for readability between specimens

    return final_consensus, all_merge_traceability, naming_info, total_limited_count


def main():
    """Main function to process command line arguments and run the summarization."""
    args = parse_arguments()

    # Parse FASTA field specification early
    try:
        fasta_fields = parse_fasta_fields(args.fasta_fields)
    except ValueError as e:
        logging.error(f"Invalid --fasta-fields specification: {e}")
        import sys
        sys.exit(1)

    # Set up logging with temporary log file
    import tempfile
    temp_log_file = tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False)
    temp_log_file.close()

    setup_logging(args.log_level, temp_log_file.name)

    logging.info("Starting enhanced speconsense summarization")
    logging.info(f"FASTA fields: {args.fasta_fields}")
    logging.info("Processing each specimen file independently to organize variants within specimens")
    
    # Load all consensus sequences
    consensus_list = load_consensus_sequences(args.source, args.min_ric)
    if not consensus_list:
        logging.error("No consensus sequences found")
        return
    
    # Group consensus sequences by input file (one file per specimen)
    file_groups = defaultdict(list)
    for cons in consensus_list:
        file_groups[cons.file_path].append(cons)

    # Create output directories before processing
    os.makedirs(args.summary_dir, exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'FASTQ Files'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'variants'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'variants', 'FASTQ Files'), exist_ok=True)

    # Build lookup tables once before processing loop
    fastq_lookup = build_fastq_lookup_table(args.source)
    original_consensus_lookup = {cons.sample_name: cons for cons in consensus_list}

    # Process each specimen file independently
    all_final_consensus = []
    all_merge_traceability = {}
    all_naming_info = {}
    all_raw_consensuses = []  # Collect .raw files from all specimens
    total_limited_merges = 0

    sorted_file_paths = sorted(file_groups.keys())
    for file_path in tqdm(sorted_file_paths, desc="Processing specimens", unit="specimen"):
        file_consensuses = file_groups[file_path]

        # Process specimen
        final_consensus, merge_traceability, naming_info, limited_count = process_single_specimen(
            file_consensuses, args
        )

        # Write individual data files immediately
        specimen_raw_consensuses = write_specimen_data_files(
            final_consensus,
            merge_traceability,
            naming_info,
            args.summary_dir,
            os.path.join(args.summary_dir, 'FASTQ Files'),
            fastq_lookup,
            original_consensus_lookup,
            fasta_fields
        )

        # Accumulate results for summary files
        all_final_consensus.extend(final_consensus)
        all_merge_traceability.update(merge_traceability)
        all_raw_consensuses.extend(specimen_raw_consensuses)
        total_limited_merges += limited_count

        # Update naming info with unique keys per specimen
        file_name = os.path.basename(file_path)
        for group_id, group_naming in naming_info.items():
            unique_key = f"{file_name}_{group_id}"
            all_naming_info[unique_key] = group_naming

    # Write summary files at end (after all processing)
    write_output_files(
        all_final_consensus,
        all_raw_consensuses,
        args.summary_dir,
        temp_log_file.name,
        fasta_fields
    )

    # Write quality report
    write_quality_report(
        all_final_consensus,
        all_raw_consensuses,
        args.summary_dir,
        args.source
    )

    logging.info(f"Enhanced summarization completed successfully")
    logging.info(f"Final output: {len(all_final_consensus)} consensus sequences in {args.summary_dir}")

    # Report if any variant groups were potentially suboptimal due to size
    if total_limited_merges > 0:
        logging.info(f"Note: {total_limited_merges} variant group(s) had >{MAX_MSA_MERGE_VARIANTS} variants (results potentially suboptimal)")

    # Clean up temporary log file
    try:
        os.unlink(temp_log_file.name)
    except Exception as e:
        logging.debug(f"Could not clean up temporary log file: {e}")


if __name__ == "__main__":
    main()