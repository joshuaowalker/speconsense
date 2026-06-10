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
from speconsense.significance import DEFAULT_MIN_CER_FACTOR
from speconsense.msa import DEFAULT_MAX_ERR_FACTOR

from .fields import FastaField, format_fasta_header


_CLUSTER_SUFFIX_RE = re.compile(
    r'-\d+(?:\.v\d+|-full)'
    r'(?:\.raw(?:\d+|\.\d+\.v\d+))?'
    r'(?:\.(?:ns|lq|filtered))?$'
)


def strip_cluster_suffix(sample_name: str) -> str:
    """Return ``sample_name`` with the core/summarize cluster suffix removed.

    Handles the end-to-end gid.vid naming conventions:
    - ``specimen-1.v2`` → ``specimen``
    - ``specimen-1.v2.raw3`` → ``specimen`` (legacy)
    - ``specimen-1.v2.raw.1.v3`` → ``specimen`` (traceability)
    - ``specimen-1.v2.ns`` → ``specimen``
    - ``specimen-1.v2.lq`` → ``specimen``
    - ``specimen-1.v2.filtered`` → ``specimen``
    - ``specimen-1-full`` → ``specimen``

    Returns the input unchanged when no suffix is present.
    """
    return _CLUSTER_SUFFIX_RE.sub('', sample_name)


def parse_consensus_header(header: str) -> Tuple[Optional[str], Optional[int], Optional[int],
                                                   Optional[List[str]], Optional[float], Optional[float],
                                                   Optional[List[int]], Optional[List[int]], Optional[int],
                                                   Optional[float], Optional[float],
                                                   Optional[int], Optional[int]]:
    """
    Extract information from Speconsense consensus FASTA header.

    Parses read identity metrics, merge metadata, quality metrics
    (``cer_factor`` CER decision metric and ``err_factor`` cluster-homogeneity
    metric), and the core-assigned identity ranks (``group_rank`` / ``variant_rank``
    emitted as ``gid=`` / ``vid=``). Full per-position CER detail (p*, K, context
    tags, q_ctx values) lives only in the metadata JSON.

    Returns:
        Tuple of (sample_name, ric, size, primers, rid, rid_min, raw_ric,
        raw_len, snp_count, cer_factor, err_factor, group_rank, variant_rank).
    """
    sample_match = re.match(r'>([^ ]+) (.+)', header)
    if not sample_match:
        return (None, None, None, None, None, None, None, None, None,
                None, None, None, None)

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

    # Extract merge metadata (plus-delimited lists)
    rawric_match = re.search(r'rawric=(\d+(?:\+\d+)*)', info_string)
    raw_ric = [int(x) for x in rawric_match.group(1).split('+')] if rawric_match else None

    rawlen_match = re.search(r'rawlen=(\d+(?:\+\d+)*)', info_string)
    raw_len = [int(x) for x in rawlen_match.group(1).split('+')] if rawlen_match else None

    snp_match = re.search(r'snp=(\d+)', info_string)
    snp_count = int(snp_match.group(1)) if snp_match else None

    # CER fields
    cer_factor_match = re.search(r'\bcer_factor=(inf|[\d.]+)', info_string)
    if cer_factor_match:
        token = cer_factor_match.group(1)
        cer_factor = float('inf') if token == 'inf' else float(token)
    else:
        cer_factor = None

    err_factor_match = re.search(r'\berr_factor=([\d.]+)', info_string)
    err_factor = float(err_factor_match.group(1)) if err_factor_match else None

    gid_match = re.search(r'\bgid=(\d+)', info_string)
    group_rank = int(gid_match.group(1)) if gid_match else None

    vid_match = re.search(r'\bvid=(\d+)', info_string)
    variant_rank = int(vid_match.group(1)) if vid_match else None

    return (sample_name, ric, size, primers, rid, rid_min,
            raw_ric, raw_len, snp_count, cer_factor, err_factor,
            group_rank, variant_rank)


def load_consensus_sequences(
    source_folder: str,
    min_ric: int,
    min_len: int = 0,
    max_len: int = 0,
    specimen_id: str = None,
    min_cer_factor: float = DEFAULT_MIN_CER_FACTOR,
    max_err_factor: float = DEFAULT_MAX_ERR_FACTOR,
) -> Tuple[List[ConsensusInfo], List[ConsensusInfo], List[ConsensusInfo]]:
    """Load consensus sequences from speconsense output files.

    Args:
        source_folder: Directory containing speconsense output files
        min_ric: Minimum Reads in Consensus threshold
        min_len: Minimum sequence length (0 = disabled)
        max_len: Maximum sequence length (0 = disabled)
        specimen_id: If set, load only {specimen_id}-all.fasta instead of all
        min_cer_factor: Minimum per-position CER factor for primary output.
            Variants with cer_factor below this threshold are returned
            separately as ns records. Variants with cer_factor=None (anchors,
            no valid pairwise comparison, or legacy pre-CER output) always
            pass. Set to 0 to disable CER filtering.
        max_err_factor: Maximum cluster err_factor (observed/q_ctx-expected
            disagreement ratio). Variants above this threshold are returned
            separately as lq (low-quality) records. Variants with
            err_factor=None (legacy output missing the field) always pass.
            Set to 0 to disable err_factor filtering. When both filters fire
            on the same record, lq takes precedence over ns.

    Returns:
        Tuple of (passing, ns, lq) lists of ConsensusInfo. All lists have
        already passed the RiC and length filters.
    """
    consensus_list: List[ConsensusInfo] = []
    ns_list: List[ConsensusInfo] = []
    lq_list: List[ConsensusInfo] = []
    filtered_by_ric = 0
    filtered_by_len = 0

    cer_filter_enabled = min_cer_factor > 0
    err_filter_enabled = max_err_factor > 0

    # Find consensus FASTA files — narrow to single specimen if requested
    if specimen_id:
        fasta_pattern = os.path.join(source_folder, f"{specimen_id}-all.fasta")
    else:
        fasta_pattern = os.path.join(source_folder, "*-all.fasta")
    fasta_files = sorted(glob.glob(fasta_pattern))

    for fasta_file in fasta_files:
        logging.debug(f"Processing consensus file: {fasta_file}")

        # Pre-load per-specimen metadata JSON once (best effort) so each
        # ConsensusInfo can carry the raw err_factor sums needed for the
        # run-scale q_ctx calibration check in quality_report.
        specimen_base = os.path.basename(fasta_file)
        if specimen_base.endswith('-all.fasta'):
            specimen_base = specimen_base[:-len('-all.fasta')]
        metadata = load_metadata_from_json(source_folder, specimen_base)
        ef_lookup: Dict[str, Tuple[Optional[float], Optional[float], Optional[int]]] = {}
        if metadata:
            for variant in metadata.get('variants', []) or []:
                cid = variant.get('cluster_id')
                if cid:
                    ef_lookup[cid] = (
                        variant.get('err_factor_obs_sum'),
                        variant.get('err_factor_exp_sum'),
                        variant.get('err_factor_cols'),
                    )

        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                (sample_name, ric, size, primers, rid, rid_min,
                 _, _, _, cer_factor, err_factor,
                 group_rank, variant_rank) = \
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

                # Extract cluster ID from sample name (e.g., "sample-1.v2" -> "1.v2").
                # Core emits the gid.vid designator directly; preserve it verbatim
                # so summarize can round-trip the identifier.
                cluster_match = re.search(r'-(\d+\.v\d+)$', sample_name)
                cluster_id = cluster_match.group(1) if cluster_match else sample_name

                ef_obs, ef_exp, ef_cols = ef_lookup.get(cluster_id, (None, None, None))

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
                    cer_factor=cer_factor,
                    err_factor=err_factor,
                    err_factor_obs_sum=ef_obs,
                    err_factor_exp_sum=ef_exp,
                    err_factor_cols=ef_cols,
                    group_rank=group_rank,
                    variant_rank=variant_rank,
                )

                # Routing priority: lq (low quality) takes precedence over ns
                # (CER not significant). A cluster whose reads are internally
                # incoherent isn't worth evaluating for peer-artifact status.
                if (err_filter_enabled and err_factor is not None
                        and err_factor > max_err_factor):
                    lq_list.append(consensus_info)
                elif (cer_filter_enabled and cer_factor is not None
                        and cer_factor < min_cer_factor):
                    ns_list.append(consensus_info)
                else:
                    consensus_list.append(consensus_info)

    # Log loading summary
    filter_parts = [f"Loaded {len(consensus_list)} consensus sequences from {len(fasta_files)} files"]
    if lq_list:
        filter_parts.append(f"routed {len(lq_list)} to lq (err_factor > {max_err_factor})")
    if ns_list:
        filter_parts.append(f"routed {len(ns_list)} to ns (CER factor < {min_cer_factor})")
    if filtered_by_ric > 0:
        filter_parts.append(f"filtered {filtered_by_ric} by RiC")
    if filtered_by_len > 0:
        filter_parts.append(f"filtered {filtered_by_len} by length")
    logging.info(", ".join(filter_parts))

    return consensus_list, ns_list, lq_list


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
        # Extract specimen name from cluster name. Core emits the gid.vid suffix
        # format (e.g., "sample-1.v2"), so trim on the last "-{digits}.v{digits}".
        specimen_match = re.match(r'^(.+)-\d+\.v\d+$', cluster_name)
        if specimen_match:
            specimen_name = specimen_match.group(1)
            debug_files = fastq_lookup.get(specimen_name, [])

            # Get the original RiC value for this cluster
            original_ric = original_consensus_lookup.get(cluster_name)
            if not original_ric:
                logging.warning(f"Could not find original consensus info for {cluster_name}")
                continue

            # Filter files that match this specific cluster with exact RiC value
            # Match: {specimen}-{gid}.v{vid}-RiC{exact_ric}-{stage}.fastq.
            # The exact RiC filter prevents matching multiple RiC values for the
            # same cluster id.
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
                               fasta_fields: List[FastaField],
                               full_reads_by_name: Optional[Dict[str, List]] = None,
                               ) -> List[Tuple[ConsensusInfo, str]]:
    """
    Write individual FASTA and FASTQ files for a single specimen.
    Does NOT write summary files (summary.fasta, summary.txt).

    Args:
        fasta_fields: List of FastaField objects defining header format
        full_reads_by_name: Optional map of ``-full`` sample_name -> list
            of sampled SeqRecord reads. When present, the matching
            ``-full`` records emit a FASTQ of these reads in place of the
            usual cluster_debug lookup (which has no matching files for
            the synthetic ``-full`` record).

    Returns:
        List of (raw_consensus, original_cluster_name) tuples for later use in summary.fasta
    """
    full_reads_by_name = full_reads_by_name or {}
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

                # Create .raw file entries. Name by original core gid.vid for
                # direct traceability to the pre-merge core output
                # (e.g. specimen-1.v1.raw.2.v3 traces back to core's 2.v3).
                for raw_info in contributing_infos:
                    raw_name = f"{consensus.sample_name}.raw.{raw_info.group_rank}.v{raw_info.variant_rank}"

                    # Build .raw ConsensusInfo from the pre-merge cluster's full
                    # field set. Mirrors the .ns/.lq pass-through treatment so
                    # every per-cluster metric core emitted (cer_factor,
                    # err_factor, gid/vid, err_factor raw sums for the quality
                    # report's calibration check) is preserved on the .raw
                    # output. Only the merge-only fields (snp_count, raw_ric,
                    # raw_len, merge_indel_count) are reset since the .raw
                    # represents a single pre-merge cluster.
                    raw_consensus = raw_info._replace(
                        sample_name=raw_name,
                        snp_count=None,
                        raw_ric=None,
                        raw_len=None,
                        merge_indel_count=None,
                    )
                    raw_file_consensuses.append((raw_consensus, raw_info.sample_name))

    # Write individual FASTA files with custom field formatting
    for consensus in specimen_consensus:
        output_file = os.path.join(summary_folder, f"{consensus.sample_name}-RiC{consensus.ric}.fasta")
        with open(output_file, 'w') as f:
            header = format_fasta_header(consensus, fasta_fields)
            f.write(f">{header}\n")
            f.write(f"{consensus.sequence}\n")

    # Write FASTQ files for each final consensus containing all contributing reads.
    # -full records are synthetic (no single source cluster_debug FASTQ); their
    # sampled reads come in via full_reads_by_name and we write them directly.
    for consensus in specimen_consensus:
        if consensus.sample_name in full_reads_by_name:
            os.makedirs(fastq_dir, exist_ok=True)
            fastq_output_path = os.path.join(
                fastq_dir,
                f"{consensus.sample_name}-RiC{consensus.ric}.fastq",
            )
            try:
                with open(fastq_output_path, 'w') as outf:
                    SeqIO.write(full_reads_by_name[consensus.sample_name], outf, "fastq")
                logging.debug(f"Wrote -full FASTQ: {os.path.basename(fastq_output_path)}")
            except Exception as e:
                logging.warning(f"Could not write -full FASTQ for {consensus.sample_name}: {e}")
            continue
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
        specimen_name = strip_cluster_suffix(original_cluster_name)
        if specimen_name != original_cluster_name:
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


def _write_filtered_variant_files(
    consensuses: List[ConsensusInfo],
    summary_folder: str,
    fastq_lookup: Dict[str, List[str]],
    fasta_fields: List[FastaField],
    marker: str,
) -> None:
    """Emit filtered variants to ``variants/`` using the raw-variant layout.

    Writes one FASTA per record to
    ``variants/{sample_name}.{marker}-RiC{ric}.fasta`` and a matching FASTQ
    when the source cluster's debug FASTQ is available.
    """
    for consensus in consensuses:
        tagged_name = f"{consensus.sample_name}.{marker}"
        tagged_info = consensus._replace(sample_name=tagged_name)

        output_file = os.path.join(
            summary_folder, 'variants', f"{tagged_name}-RiC{consensus.ric}.fasta"
        )
        with open(output_file, 'w') as f:
            header = format_fasta_header(tagged_info, fasta_fields)
            f.write(f">{header}\n")
            f.write(f"{consensus.sequence}\n")

        original_cluster_name = consensus.sample_name
        specimen_name = strip_cluster_suffix(original_cluster_name)
        if specimen_name == original_cluster_name:
            continue
        debug_files = fastq_lookup.get(specimen_name, []) if fastq_lookup else []
        cluster_ric_pattern = f"{original_cluster_name}-RiC{consensus.ric}-"
        matching_files = [f for f in debug_files if cluster_ric_pattern in f]
        if not matching_files:
            continue

        fastq_output_path = os.path.join(
            summary_folder, 'variants', 'FASTQ Files',
            f"{tagged_name}-RiC{consensus.ric}.fastq"
        )
        try:
            with open(fastq_output_path, 'wb') as outf:
                for input_file in matching_files:
                    if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
                        with open(input_file, 'rb') as inf:
                            shutil.copyfileobj(inf, outf)
            logging.debug(f"Wrote {marker} FASTQ: {os.path.basename(fastq_output_path)}")
        except Exception as e:
            logging.debug(f"Could not write {marker} FASTQ for {tagged_name}: {e}")


def write_ns_variant_files(ns_consensuses: List[ConsensusInfo],
                           summary_folder: str,
                           fastq_lookup: Dict[str, List[str]],
                           fasta_fields: List[FastaField]) -> None:
    """Emit CER-filtered (ns) variants to variants/ using the raw-variant layout."""
    _write_filtered_variant_files(
        ns_consensuses, summary_folder, fastq_lookup, fasta_fields, marker='ns'
    )


def write_lq_variant_files(lq_consensuses: List[ConsensusInfo],
                           summary_folder: str,
                           fastq_lookup: Dict[str, List[str]],
                           fasta_fields: List[FastaField]) -> None:
    """Emit err_factor-filtered (lq) variants to variants/ using the raw-variant layout."""
    _write_filtered_variant_files(
        lq_consensuses, summary_folder, fastq_lookup, fasta_fields, marker='lq'
    )


def write_filtered_variant_files(filtered_consensuses: List[ConsensusInfo],
                                  summary_folder: str,
                                  fastq_lookup: Dict[str, List[str]],
                                  fasta_fields: List[FastaField]) -> None:
    """Emit selection-filtered variants to variants/ using the raw-variant layout.

    Selection-filtered variants passed quality gates (cer_factor, err_factor)
    but were excluded by selection or pruning parameters (--select-max-variants,
    --select-min-size-ratio, --select-max-groups, --prune-group-frac/abs).
    This distinguishes "quality problem" (.ns/.lq) from "selection decision"
    (.filtered).
    """
    _write_filtered_variant_files(
        filtered_consensuses, summary_folder, fastq_lookup, fasta_fields,
        marker='filtered',
    )


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
        # Pattern: {specimen}-{gid}.v{vid}-RiC{size}-{stage}.fastq
        # Where stage can be: sampled, reads, untrimmed, or other variants
        pattern = re.compile(r'^(.+)-(\d+\.v\d+)-RiC(\d+)-([a-z]+)\.fastq$')

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


def load_existing_specimen_outputs(summary_dir: str) -> List[ConsensusInfo]:
    """Load per-specimen FASTA outputs for aggregate-only mode.

    Scans summary_dir for individual specimen FASTA files (excluding summary.fasta),
    reconstructing ConsensusInfo objects from their headers.

    Returns:
        List of ConsensusInfo objects from existing per-specimen files
    """
    consensus_list = []
    fasta_pattern = os.path.join(summary_dir, "*.fasta")
    fasta_files = sorted(glob.glob(fasta_pattern))

    for fasta_file in fasta_files:
        basename = os.path.basename(fasta_file)
        if basename == "summary.fasta":
            continue

        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                (sample_name, ric, size, primers, rid, rid_min,
                 raw_ric, raw_len, snp_count, cer_factor, err_factor,
                 group_rank, variant_rank) = \
                    parse_consensus_header(f">{record.description}")

                if not sample_name:
                    continue

                cluster_match = re.search(r'-(\d+\.v\d+)$', sample_name)
                cluster_id = cluster_match.group(1) if cluster_match else sample_name

                consensus_info = ConsensusInfo(
                    sample_name=sample_name,
                    cluster_id=cluster_id,
                    sequence=str(record.seq),
                    ric=ric,
                    size=size,
                    file_path=fasta_file,
                    snp_count=snp_count,
                    primers=primers,
                    raw_ric=raw_ric,
                    raw_len=raw_len,
                    rid=rid,
                    rid_min=rid_min,
                    cer_factor=cer_factor,
                    err_factor=err_factor,
                    group_rank=group_rank,
                    variant_rank=variant_rank,
                )
                consensus_list.append(consensus_info)

    logging.info(f"Loaded {len(consensus_list)} existing sequences from {len(fasta_files)} files in {summary_dir}")
    return consensus_list


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
