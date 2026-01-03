"""
Summarize subpackage for speconsense.

Provides post-processing of speconsense output: HAC variant grouping,
MSA-based merging with IUPAC ambiguity codes, and variant selection.
"""

# CLI and entry point
from .cli import main, parse_arguments, setup_logging, process_single_specimen

# IUPAC utilities and distance functions
from .iupac import (
    IUPAC_EQUIV,
    STANDARD_ADJUSTMENT_PARAMS,
    bases_match_with_iupac,
    expand_iupac_code,
    merge_bases_to_iupac,
    calculate_adjusted_identity_distance,
    calculate_overlap_aware_distance,
    create_variant_summary,
    primers_are_same,
)

# FASTA field classes
from .fields import (
    FastaField,
    FASTA_FIELDS,
    FASTA_FIELD_PRESETS,
    validate_field_registry,
    parse_fasta_fields,
    format_fasta_header,
)

# MSA analysis and quality assessment
from .analysis import (
    ClusterQualityData,
    MAX_MSA_MERGE_VARIANTS,
    run_spoa_msa,
    identify_indel_events,
    is_homopolymer_event,
    analyze_msa_columns,
    analyze_msa_columns_overlap_aware,
    analyze_cluster_quality,
    identify_outliers,
    analyze_positional_identity_outliers,
)

# MSA-based variant merging
from .merging import (
    generate_all_subsets_by_size,
    is_compatible_subset,
    create_consensus_from_msa,
    create_overlap_consensus_from_msa,
    merge_group_with_msa,
)

# HAC clustering and variant selection
from .clustering import (
    perform_hac_clustering,
    select_variants,
)

# File I/O operations
from .io import (
    parse_consensus_header,
    load_consensus_sequences,
    load_metadata_from_json,
    build_fastq_lookup_table,
    create_output_structure,
    write_consensus_fastq,
    write_specimen_data_files,
    write_position_debug_file,
    write_output_files,
)

__all__ = [
    # CLI
    "main",
    "parse_arguments",
    "setup_logging",
    "process_single_specimen",
    # IUPAC
    "IUPAC_EQUIV",
    "STANDARD_ADJUSTMENT_PARAMS",
    "bases_match_with_iupac",
    "expand_iupac_code",
    "merge_bases_to_iupac",
    "calculate_adjusted_identity_distance",
    "calculate_overlap_aware_distance",
    "create_variant_summary",
    "primers_are_same",
    # Fields
    "FastaField",
    "FASTA_FIELDS",
    "FASTA_FIELD_PRESETS",
    "validate_field_registry",
    "parse_fasta_fields",
    "format_fasta_header",
    # Analysis
    "ClusterQualityData",
    "MAX_MSA_MERGE_VARIANTS",
    "run_spoa_msa",
    "identify_indel_events",
    "is_homopolymer_event",
    "analyze_msa_columns",
    "analyze_msa_columns_overlap_aware",
    "analyze_cluster_quality",
    "identify_outliers",
    "analyze_positional_identity_outliers",
    # Merging
    "generate_all_subsets_by_size",
    "is_compatible_subset",
    "create_consensus_from_msa",
    "create_overlap_consensus_from_msa",
    "merge_group_with_msa",
    # Clustering
    "perform_hac_clustering",
    "select_variants",
    # I/O
    "parse_consensus_header",
    "load_consensus_sequences",
    "load_metadata_from_json",
    "build_fastq_lookup_table",
    "create_output_structure",
    "write_consensus_fastq",
    "write_specimen_data_files",
    "write_position_debug_file",
    "write_output_files",
]
