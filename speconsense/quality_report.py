"""
Quality report generation for speconsense-summarize.

This module handles the generation of quality reports with multiple analysis sections:
- Executive Summary
- Read Identity Analysis
- Positional Identity Analysis
- Overlap Merge Analysis
- Interpretation Guide
"""

import logging
import os
import re
from datetime import datetime
from typing import Dict, List, Tuple, TextIO

from tqdm import tqdm

# Import shared types
from speconsense.types import ConsensusInfo, OverlapMergeInfo

# Import helper functions from summarize (safe because summarize uses deferred import for this module)
from speconsense.summarize import (
    identify_outliers,
    analyze_positional_identity_outliers,
    load_metadata_from_json,
    write_position_debug_file,
)


def write_header_section(f: TextIO, source_folder: str):
    """Write the report header with timestamp and source info."""
    f.write("=" * 80 + "\n")
    f.write("QUALITY REPORT - speconsense-summarize\n")
    f.write("=" * 80 + "\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Source: {source_folder}\n")
    f.write("=" * 80 + "\n\n")


def write_executive_summary_section(
    f: TextIO,
    total_seqs: int,
    total_merged: int,
    stats: Dict,
    n_stat: int,
    n_pos: int,
    n_merged: int,
    total_flagged: int
):
    """Write the executive summary section with high-level statistics."""
    f.write("EXECUTIVE SUMMARY\n")
    f.write("-" * 80 + "\n\n")

    f.write(f"Total sequences: {total_seqs}\n")
    f.write(f"Merged sequences: {total_merged} ({100*total_merged/total_seqs:.1f}%)\n\n")

    f.write("Global Read Identity Statistics:\n")
    f.write(f"  Mean rid: {stats['mean_rid']:.1%} ± {stats['std_rid']:.1%}\n\n")

    f.write("Sequences Requiring Attention:\n")
    n_rid_issues = n_stat + n_merged
    f.write(f"  Total flagged: {total_flagged} ({100*total_flagged/total_seqs:.1f}%)\n")
    f.write(f"    - Low read identity: {n_rid_issues} ({n_stat} sequences + {n_merged} merged)\n")
    f.write(f"    - High-error positions: {n_pos}\n\n")


def write_read_identity_section(
    f: TextIO,
    outlier_results: Dict,
    merged_with_issues: List[Tuple[ConsensusInfo, List, float, float]],
    stats: Dict
):
    """Write the read identity analysis section."""
    n_stat = len(outlier_results['statistical_outliers'])
    n_merged = len(merged_with_issues)
    n_rid_issues = n_stat + n_merged

    if n_rid_issues == 0:
        return

    f.write("=" * 80 + "\n")
    f.write("READ IDENTITY ANALYSIS\n")
    f.write("=" * 80 + "\n\n")

    f.write("Sequences with mean read identity (rid) below mean - 2×std.\n")
    f.write(f"Threshold: {stats['stat_threshold_rid']:.1%}\n\n")

    f.write(f"{'Sequence':<50} {'RiC':<6} {'rid':<8}\n")
    f.write("-" * 64 + "\n")

    # Build combined list for sorting
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

    # Sort by rid ascending
    combined_entries.sort(key=lambda x: x[0])

    # Display entries
    for _, is_merged, data in combined_entries:
        if is_merged:
            cons, components_info, worst_rid, weighted_avg_rid = data
            name = cons.sample_name
            name_with_tag = f"{name} [merged]"
            name_truncated = name_with_tag[:49] if len(name_with_tag) > 49 else name_with_tag
            rid_str = f"{weighted_avg_rid:.1%}"
            f.write(f"{name_truncated:<50} {cons.ric:<6} {rid_str:<8}\n")

            # Component rows (indented)
            for raw, comp_rid, comp_ric in components_info:
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


def write_positional_identity_section(
    f: TextIO,
    sequences_with_pos_outliers: List[Tuple[ConsensusInfo, Dict]],
    min_variant_frequency: float,
    min_variant_count: int
):
    """Write the positional identity analysis section."""
    if not sequences_with_pos_outliers:
        return

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
        display_data.append((display_name, ric_val, cons, result))

    # Calculate column width based on longest name (minimum 40, cap at 70)
    max_name_len = max(len(name) for name, _, _, _ in display_data) if display_data else 40
    name_col_width = min(max(max_name_len + 2, 40), 70)

    f.write(f"{'Sequence':<{name_col_width}} {'RiC':<6} {'Ambig':<6} {'#Pos':<6} {'MeanErr':<8} {'TotalErr':<10}\n")
    f.write("-" * (name_col_width + 38) + "\n")

    for display_name, ric_val, cons, result in display_data:
        mean_err = result.get('mean_outlier_error_rate', 0.0)
        total_err = result.get('total_nucleotide_errors', 0)
        num_pos = result['num_outlier_positions']
        # Count IUPAC ambiguity codes in the consensus sequence (non-ACGT characters)
        ambig_count = sum(1 for c in cons.sequence if c.upper() not in 'ACGT')

        f.write(f"{display_name:<{name_col_width}} {ric_val:<6} {ambig_count:<6} {num_pos:<6} "
               f"{mean_err:<8.1%} {total_err:<10}\n")

    f.write("\n")


def write_overlap_merge_section(
    f: TextIO,
    overlap_merges: List[OverlapMergeInfo],
    min_merge_overlap: int
):
    """Write the overlap merge analysis section."""
    # Only include merges that extended beyond full overlap
    true_overlap_merges = [m for m in overlap_merges if m.prefix_bp > 0 or m.suffix_bp > 0]

    if not true_overlap_merges:
        return

    f.write("=" * 80 + "\n")
    f.write("OVERLAP MERGE ANALYSIS\n")
    f.write("=" * 80 + "\n\n")

    # Group merges by specimen
    specimen_merges: Dict[str, List[OverlapMergeInfo]] = {}
    for merge_info in true_overlap_merges:
        if merge_info.specimen not in specimen_merges:
            specimen_merges[merge_info.specimen] = []
        specimen_merges[merge_info.specimen].append(merge_info)

    f.write(f"{len(specimen_merges)} specimen(s) had overlap merges:\n\n")

    # Sort specimens by name
    for specimen in sorted(specimen_merges.keys()):
        merges = specimen_merges[specimen]
        merge_count = len(merges)
        max_iteration = max(m.iteration for m in merges)

        if max_iteration > 1:
            f.write(f"{specimen} ({merge_count} merge(s), iterative):\n")
        else:
            f.write(f"{specimen} ({merge_count} merge(s)):\n")

        # Sort by iteration
        for merge_info in sorted(merges, key=lambda m: m.iteration):
            iter_prefix = f"  Round {merge_info.iteration}: " if max_iteration > 1 else "  "

            # Format input clusters
            input_parts = []
            for cluster, length, ric in zip(
                merge_info.input_clusters,
                merge_info.input_lengths,
                merge_info.input_rics
            ):
                cluster_id = cluster.rsplit('-', 1)[-1] if '-' in cluster else cluster
                input_parts.append(f"{cluster_id} ({length}bp, RiC={ric})")

            f.write(f"{iter_prefix}Merged: {' + '.join(input_parts)} -> {merge_info.output_length}bp\n")

            # Calculate overlap as percentage of shorter sequence
            shorter_len = min(merge_info.input_lengths)
            overlap_pct = (merge_info.overlap_bp / shorter_len * 100) if shorter_len > 0 else 0
            f.write(f"    Overlap: {merge_info.overlap_bp}bp ({overlap_pct:.0f}% of shorter sequence)\n")
            f.write(f"    Extensions: prefix={merge_info.prefix_bp}bp, suffix={merge_info.suffix_bp}bp\n")

        f.write("\n")

    # Edge case warnings
    warnings = []
    for merge_info in true_overlap_merges:
        # Warn if overlap is within 10% of threshold
        if merge_info.overlap_bp < min_merge_overlap * 1.1:
            shorter_len = min(merge_info.input_lengths)
            if merge_info.overlap_bp < shorter_len:
                warnings.append(
                    f"{merge_info.specimen}: Small overlap relative to threshold "
                    f"({merge_info.overlap_bp}bp, threshold={min_merge_overlap}bp)"
                )

        # Warn if large length ratio (>3:1)
        max_len = max(merge_info.input_lengths)
        min_len = min(merge_info.input_lengths)
        if max_len > min_len * 3:
            warnings.append(
                f"{merge_info.specimen}: Large length ratio "
                f"({max_len}bp / {min_len}bp = {max_len/min_len:.1f}x)"
            )

    if warnings:
        f.write("Attention:\n")
        for warning in warnings:
            f.write(f"  * {warning}\n")
        f.write("\n")


def write_interpretation_guide_section(f: TextIO):
    """Write the interpretation guide section."""
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
    f.write("  Ambig: Count of IUPAC ambiguity codes in consensus\n")
    f.write("  #Pos: Count of positions exceeding error threshold\n")
    f.write("  MeanErr: Average error rate at flagged positions\n")
    f.write("  TotalErr: Sum of errors at flagged positions\n\n")


def write_quality_report(
    final_consensus: List[ConsensusInfo],
    all_raw_consensuses: List[Tuple[ConsensusInfo, str]],
    summary_folder: str,
    source_folder: str,
    overlap_merges: List[OverlapMergeInfo] = None,
    min_merge_overlap: int = 200
):
    """
    Write quality report with rid-based dual outlier detection.

    Uses mean read identity (rid) for outlier detection. rid_min is not used
    because single outlier reads don't significantly impact consensus quality;
    positional analysis better captures systematic issues.

    Structure:
    1. Executive Summary - High-level overview with attention flags
    2. Read Identity Analysis - Dual outlier detection (clustering threshold + statistical)
    3. Positional Identity Analysis - Sequences with problematic positions
    4. Overlap Merge Analysis - Details of overlap merges (when applicable)
    5. Interpretation Guide - Actionable guidance with neutral tone

    Args:
        final_consensus: List of final consensus sequences
        all_raw_consensuses: List of (raw_consensus, original_name) tuples
        summary_folder: Output directory for report
        source_folder: Source directory containing cluster_debug with MSA files
        overlap_merges: List of OverlapMergeInfo objects describing overlap merges
        min_merge_overlap: Threshold used for overlap merging (for edge case warnings)
    """
    if overlap_merges is None:
        overlap_merges = []

    quality_report_path = os.path.join(summary_folder, 'quality_report.txt')

    # Build .raw lookup: map merged sequence names to their .raw components
    raw_lookup: Dict[str, List[ConsensusInfo]] = {}
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
    min_variant_frequency = None
    min_variant_count = None

    for cons in final_consensus:
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
    sequences_with_pos_outliers: List[Tuple[ConsensusInfo, Dict]] = []
    sequences_to_analyze = {cons.sample_name: cons for cons in final_consensus}

    logging.info("Analyzing positional identity for quality report...")
    for cons in tqdm(sequences_to_analyze.values(), desc="Analyzing positional identity", unit="seq"):
        is_merged = cons.snp_count is not None and cons.snp_count > 0

        if is_merged:
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

            if worst_result and worst_result['num_outlier_positions'] > 0:
                sequences_with_pos_outliers.append((cons, worst_result))
        else:
            result = analyze_positional_identity_outliers(
                cons, source_folder, min_variant_frequency, min_variant_count
            )
            if result and result['num_outlier_positions'] > 0:
                sequences_with_pos_outliers.append((cons, result))

    sequences_with_pos_outliers.sort(key=lambda x: x[1].get('total_nucleotide_errors', 0), reverse=True)

    # Write detailed position debug file
    if sequences_with_pos_outliers:
        write_position_debug_file(sequences_with_pos_outliers, summary_folder, min_variant_frequency)

    # Identify merged sequences with quality issues
    merged_with_issues: List[Tuple[ConsensusInfo, List, float, float]] = []
    threshold_rid = outlier_results['global_stats']['stat_threshold_rid']

    for cons in final_consensus:
        is_merged = cons.snp_count is not None and cons.snp_count > 0
        if is_merged:
            raw_components = raw_lookup.get(cons.sample_name, [])
            if not raw_components:
                continue

            components_info = []
            worst_rid = 1.0
            total_ric = 0
            weighted_rid_sum = 0.0

            for raw in raw_components:
                rid = raw.rid if raw.rid is not None else 1.0
                ric = raw.ric if raw.ric else 0
                components_info.append((raw, rid, ric))
                if rid < worst_rid:
                    worst_rid = rid
                if ric > 0:
                    total_ric += ric
                    weighted_rid_sum += rid * ric

            weighted_avg_rid = weighted_rid_sum / total_ric if total_ric > 0 else 1.0
            components_info.sort(key=lambda x: x[1])

            if worst_rid < threshold_rid:
                merged_with_issues.append((cons, components_info, worst_rid, weighted_avg_rid))

    merged_with_issues.sort(key=lambda x: x[2])

    # Calculate summary statistics
    total_seqs = len(final_consensus)
    total_merged = sum(1 for cons in final_consensus if cons.snp_count is not None and cons.snp_count > 0)
    stats = outlier_results['global_stats']
    n_stat = len(outlier_results['statistical_outliers'])
    n_pos = len(sequences_with_pos_outliers)
    n_merged = len(merged_with_issues)
    n_rid_issues = n_stat + n_merged

    # Count unique flagged sequences
    flagged_names = set()
    for c, _ in outlier_results['statistical_outliers']:
        flagged_names.add(c.sample_name)
    for c, _ in sequences_with_pos_outliers:
        flagged_names.add(c.sample_name)
    for c, _, _, _ in merged_with_issues:
        flagged_names.add(c.sample_name)
    total_flagged = len(flagged_names)

    # Write the report
    with open(quality_report_path, 'w') as f:
        write_header_section(f, source_folder)

        write_executive_summary_section(
            f, total_seqs, total_merged, stats,
            n_stat, n_pos, n_merged, total_flagged
        )

        write_read_identity_section(f, outlier_results, merged_with_issues, stats)

        write_positional_identity_section(
            f, sequences_with_pos_outliers,
            min_variant_frequency, min_variant_count
        )

        write_overlap_merge_section(f, overlap_merges, min_merge_overlap)

        write_interpretation_guide_section(f)

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
