"""File input/output operations for speconsense-summarize.

Provides functions for loading consensus sequences, writing output files,
and managing the output directory structure.
"""

import os
import re
import glob
import csv
import json
import shutil
import logging
import datetime
from typing import List, Dict, Tuple, Optional
from collections import defaultdict

from Bio import SeqIO

from speconsense.types import ConsensusInfo

from .fields import FastaField, format_fasta_header
from .clustering import select_variants


def parse_consensus_header(header: str) -> Tuple[Optional[str], Optional[int], Optional[int],
                                                   Optional[List[str]], Optional[float], Optional[float]]:
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


def load_consensus_sequences(
    source_folder: str,
    min_ric: int,
    min_len: int = 0,
    max_len: int = 0
) -> List[ConsensusInfo]:
    """Load all consensus sequences from speconsense output files.

    Args:
        source_folder: Directory containing speconsense output files
        min_ric: Minimum Reads in Consensus threshold
        min_len: Minimum sequence length (0 = disabled)
        max_len: Maximum sequence length (0 = disabled)

    Returns:
        List of ConsensusInfo objects passing all filters
    """
    consensus_list = []
    filtered_by_ric = 0
    filtered_by_len = 0

    # Find all consensus FASTA files matching the new naming pattern
    fasta_pattern = os.path.join(source_folder, "*-all.fasta")
    fasta_files = sorted(glob.glob(fasta_pattern))

    for fasta_file in fasta_files:
        logging.debug(f"Processing consensus file: {fasta_file}")

        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                sample_name, ric, size, primers, rid, rid_min = \
                    parse_consensus_header(f">{record.description}")

                if not sample_name:
                    continue

                # RiC filter
                if ric < min_ric:
                    filtered_by_ric += 1
                    continue

                # Length filters (applied before merging to avoid chimeric contamination)
                seq_len = len(record.seq)
                if min_len > 0 and seq_len < min_len:
                    logging.debug(f"Filtered {sample_name}: length {seq_len} < min_len {min_len}")
                    filtered_by_len += 1
                    continue
                if max_len > 0 and seq_len > max_len:
                    logging.debug(f"Filtered {sample_name}: length {seq_len} > max_len {max_len}")
                    filtered_by_len += 1
                    continue

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

    # Log loading summary
    filter_parts = [f"Loaded {len(consensus_list)} consensus sequences from {len(fasta_files)} files"]
    if filtered_by_ric > 0:
        filter_parts.append(f"filtered {filtered_by_ric} by RiC")
    if filtered_by_len > 0:
        filter_parts.append(f"filtered {filtered_by_len} by length")
    logging.info(", ".join(filter_parts))

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
            # Use rsplit to split on the LAST '-c' (specimen names may contain '-c')
            specimen_base = variant.sample_name.rsplit('-c', 1)[0]
            new_name = f"{specimen_base}-{group_idx}.v{variant_idx + 1}"

            # Use _replace to preserve all fields while updating sample_name
            renamed_variant = variant._replace(sample_name=new_name)

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
    # Track cluster metadata alongside files: [(cluster_name, ric, [files]), ...]
    fastq_output_path = os.path.join(fastq_dir, f"{consensus.sample_name}-RiC{consensus.ric}.fastq")
    cluster_files = []

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
            valid_files = []
            for mf in matching_files:
                if not os.path.exists(mf):
                    logging.warning(f"Matched file does not exist: {mf}")
                elif os.path.getsize(mf) == 0:
                    logging.warning(f"Matched file is empty: {mf}")
                else:
                    valid_files.append(mf)

            if valid_files:
                cluster_files.append((cluster_name, original_ric.ric, valid_files))

    if not cluster_files:
        logging.warning(f"No FASTQ files found for {consensus.sample_name} from clusters: {original_clusters}")
        return

    # Concatenate files with cluster boundary delimiters
    # Each cluster gets a synthetic FASTQ record as a delimiter before its reads
    files_processed = 0
    try:
        with open(fastq_output_path, 'w') as outf:
            for idx, (cluster_name, ric, files) in enumerate(cluster_files, 1):
                # Count reads in this cluster's files
                cluster_reads = 0
                for f in files:
                    with open(f, 'r') as rf:
                        cluster_reads += sum(1 for _ in rf) // 4

                # Write cluster boundary delimiter
                outf.write(f"@CLUSTER_BOUNDARY_{idx}:{cluster_name}:RiC={ric}:reads={cluster_reads}\n")
                outf.write("NNNNNNNNNN\n")
                outf.write("+\n")
                outf.write("!!!!!!!!!!\n")

                # Write cluster reads
                for input_file in files:
                    try:
                        with open(input_file, 'r') as inf:
                            shutil.copyfileobj(inf, outf)
                        files_processed += 1
                    except Exception as e:
                        logging.debug(f"Could not concatenate {input_file}: {e}")

        # Check if the output file has content
        output_size = os.path.getsize(fastq_output_path)
        total_files = sum(len(files) for _, _, files in cluster_files)
        if output_size > 0:
            # Count reads for logging by quickly counting lines and dividing by 4
            with open(fastq_output_path, 'r') as f:
                line_count = sum(1 for line in f)
            read_count = line_count // 4
            logging.debug(f"Concatenated {files_processed}/{total_files} files from {len(cluster_files)} clusters ({output_size:,} bytes) with ~{read_count} reads to {fastq_output_path}")
        else:
            # Debug: check what files were supposed to be concatenated
            file_info = []
            for _, _, files in cluster_files:
                for input_file in files:
                    size = os.path.getsize(input_file) if os.path.exists(input_file) else 0
                    file_info.append(f"{os.path.basename(input_file)}:{size}B")

            logging.warning(f"No data written for {consensus.sample_name} - input files: {', '.join(file_info)}")
            # Remove empty output file
            try:
                os.unlink(fastq_output_path)
            except OSError:
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
        summarize_log_path = os.path.join(summary_folder, 'summarize_log.txt')
        try:
            # Flush any remaining log entries before copying
            logging.getLogger().handlers[1].flush() if len(logging.getLogger().handlers) > 1 else None
            shutil.copy2(temp_log_file, summarize_log_path)
            logging.info(f"Created log file: {summarize_log_path}")
        except Exception as e:
            logging.warning(f"Could not copy log file: {e}")
