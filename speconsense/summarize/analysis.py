"""MSA analysis support for speconsense-summarize.

Provides functions for running SPOA, identifying indel events, and analyzing
MSA columns (homopolymer / structural classification, overlap-aware spans).
"""

import logging
import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from speconsense.msa import MSAResult, extract_alignments_from_msa


# Maximum number of variants to evaluate for MSA-based merging (legacy constant)
# Batch size is now dynamically computed based on --merge-effort and group size.
# This constant is kept for backward compatibility and as the default MAX_MERGE_BATCH.
MAX_MSA_MERGE_VARIANTS = 8

# Merge effort batch size limits
MIN_MERGE_BATCH = 4
MAX_MERGE_BATCH = 8


def compute_merge_batch_size(group_size: int, effort: int) -> int:
    """Compute batch size for a group based on effort level.

    Uses formula: B = E + 1 - log2(V), clamped to [MIN_MERGE_BATCH, MAX_MERGE_BATCH]
    This keeps expected evaluations near 2^E per group.

    Args:
        group_size: Number of variants in the HAC group
        effort: Merge effort level (6-14, default 10)

    Returns:
        Batch size between MIN_MERGE_BATCH and MAX_MERGE_BATCH
    """
    import math

    if group_size <= 1:
        return 1

    log_v = int(math.log2(group_size))
    batch = effort + 1 - log_v

    return max(MIN_MERGE_BATCH, min(MAX_MERGE_BATCH, batch))


def run_spoa_for_cluster_metrics(
    sequences: Dict[str, str],
    disable_homopolymer_equivalence: bool = False,
    min_hp_length: int = 6,
) -> Optional[MSAResult]:
    """Run SPOA over a cluster's read set and return the parsed MSAResult.

    Mirrors core's ``_run_spoa_for_cluster_worker`` invocation (linear gap
    scoring ``-m1 -n-1 -g-1 -e-1`` which reduces ambiguous columns ~4× vs
    SPOA defaults) so that summarize-side metrics recomputed from the same
    reads are comparable to core's. Used to refresh ``rid``, ``rid_min``,
    and ``err_factor`` on merged clusters whose contributors come from
    multiple pre-merge clusters.

    Args:
        sequences: Mapping of read_id -> raw read sequence (ungapped).
        disable_homopolymer_equivalence: Pass through to
            ``extract_alignments_from_msa``. When True the per-read
            normalized edit distance is just the raw edit distance.
        min_hp_length: Minimum consensus HP run length for HP normalization
            in ``extract_alignments_from_msa``. Matches the
            ``--hp-normalization-length`` semantics shared across core and
            summarize.

    Returns:
        ``MSAResult`` with the SPOA consensus, raw MSA string (suitable for
        ``compute_cluster_err_factor``), parsed read alignments, and the
        MSA-to-consensus position map. Returns ``None`` if SPOA fails or
        the consensus can't be extracted.
    """
    if not sequences:
        return None

    temp_input = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
            for read_id, seq in sequences.items():
                f.write(f">{read_id}\n{seq}\n")
            temp_input = f.name

        cmd = [
            "spoa", temp_input,
            "-r", "2", "-l", "1", "-m", "1", "-n", "-1", "-g", "-1", "-e", "-1",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        enable_normalization = not disable_homopolymer_equivalence
        alignments, consensus, msa_to_consensus_pos = extract_alignments_from_msa(
            result.stdout,
            enable_homopolymer_normalization=enable_normalization,
            min_hp_length=min_hp_length,
        )
        if not consensus:
            return None

        return MSAResult(
            consensus=consensus,
            msa_string=result.stdout,
            alignments=alignments,
            msa_to_consensus_pos=msa_to_consensus_pos,
        )
    except subprocess.CalledProcessError as e:
        logging.debug(f"SPOA failed during merge-metric recompute: rc={e.returncode}")
        return None
    except Exception as e:
        logging.debug(f"Merge-metric recompute failed: {e}")
        return None
    finally:
        if temp_input and os.path.exists(temp_input):
            os.unlink(temp_input)


def run_spoa_msa(sequences: List[str], alignment_mode: int = 1) -> List:
    """
    Run SPOA to create multiple sequence alignment.

    Args:
        sequences: List of DNA sequence strings
        alignment_mode: SPOA alignment mode:
            0 = local (Smith-Waterman) - best for overlap merging
            1 = global (Needleman-Wunsch) - default, for same-length sequences
            2 = semi-global - alternative for overlap merging

    Returns:
        List of SeqRecord objects with aligned sequences (including gaps)
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
        try:
            # Write sequences to temporary file
            records = [
                SeqRecord(Seq(seq), id=f"seq{i}", description="")
                for i, seq in enumerate(sequences)
            ]
            SeqIO.write(records, temp_input, "fasta")
            temp_input.flush()

            # Run SPOA with alignment output (-r 2) and specified alignment mode
            result = subprocess.run(
                ['spoa', temp_input.name, '-r', '2', '-l', str(alignment_mode)],
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


def _aligned_hp_run_length(seq_str: str, start_col: int, end_col: int, base: str) -> int:
    """Measure the full HP run of ``base`` in an aligned sequence that touches
    the column span [start_col, end_col] (inclusive).

    Walks left from ``start_col-1``, across the span, and right from
    ``end_col+1``. Gap columns are skipped; non-gap columns matching ``base``
    are counted; any other non-gap base ends the walk on that side.

    Mirrors the measurement used by adjusted-identity's ``_hp_run_length``
    so that MSA-side HP classification uses the same MIN-of-runs semantics
    as the distance-side.
    """
    base_upper = base.upper()
    n = len(seq_str)
    count = 0

    # Walk left from just before start_col
    pos = start_col - 1
    while pos >= 0:
        c = seq_str[pos]
        if c == '-':
            pos -= 1
            continue
        if c.upper() == base_upper:
            count += 1
            pos -= 1
        else:
            break

    # Count base occurrences inside the span (gaps skipped)
    for pos in range(start_col, end_col + 1):
        c = seq_str[pos]
        if c == '-':
            continue
        if c.upper() == base_upper:
            count += 1

    # Walk right from just after end_col
    pos = end_col + 1
    while pos < n:
        c = seq_str[pos]
        if c == '-':
            pos += 1
            continue
        if c.upper() == base_upper:
            count += 1
            pos += 1
        else:
            break

    return count


def is_homopolymer_event(aligned_seqs: List, start_col: int, end_col: int,
                         min_hp_length: int = 1) -> bool:
    """
    Classify a complete indel event as homopolymer or structural.

    An event is homopolymer if:
    1. All bases in the event region (across all sequences, all columns) are identical
    2. At least one flanking solid column has all sequences showing the same base
    3. The shortest aligned HP run (min across sequences) is >= ``min_hp_length``.

    Default ``min_hp_length=1`` preserves the original binary classification.
    Higher thresholds demote short-HP length diffs back to structural indels,
    matching adjusted-identity's ``hp_normalize_min_length`` semantics.

    Examples:
        Homopolymer:  ATAAA--GC vs ATAAAAGC  (event has all A's, flanked by A)
        Structural:   ATAA-GC vs ATG-AGC     (event has A, flanked by A vs G)
        Structural:   ATC--GC vs ATCATGC     (event has A and T - not homopolymer)

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        start_col: First column of the indel event (inclusive)
        end_col: Last column of the indel event (inclusive)
        min_hp_length: Minimum HP run length (per-sequence, MIN across sequences)
            required to keep the event classified as homopolymer. 1 disables
            the check.

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
    has_matching_flank = False

    # Check left flank
    if start_col > 0:
        left_col = start_col - 1
        left_column = [str(seq.seq[left_col]) for seq in aligned_seqs]
        left_bases = set(c for c in left_column if c != '-')
        left_has_gap = '-' in left_column

        if not left_has_gap and left_bases == {event_base}:
            has_matching_flank = True

    # Check right flank
    if not has_matching_flank and end_col < alignment_length - 1:
        right_col = end_col + 1
        right_column = [str(seq.seq[right_col]) for seq in aligned_seqs]
        right_bases = set(c for c in right_column if c != '-')
        right_has_gap = '-' in right_column

        if not right_has_gap and right_bases == {event_base}:
            has_matching_flank = True

    if not has_matching_flank:
        return False

    # Short-HP demotion: if the shortest per-sequence HP run is below the
    # threshold, the length diff carries short-HP signal — treat as structural.
    if min_hp_length > 1:
        run_lengths = [
            _aligned_hp_run_length(str(seq.seq), start_col, end_col, event_base)
            for seq in aligned_seqs
        ]
        if run_lengths and min(run_lengths) < min_hp_length:
            return False

    return True


def analyze_msa_columns(aligned_seqs: List, min_hp_length: int = 1) -> dict:
    """
    Analyze aligned sequences to count SNPs and indels.

    Distinguishes between structural indels (real insertions/deletions) and
    homopolymer indels (length differences in homopolymer runs like AAA vs AAAA).

    Uses event-based classification: consecutive indel columns are grouped into
    events, and each complete event is classified as homopolymer or structural.

    Important: All gaps (including terminal gaps) count as variant positions
    since variants within a group share the same primers.

    Args:
        aligned_seqs: SPOA MSA sequences
        min_hp_length: Minimum HP run length to classify an event as
            homopolymer. Events whose shortest per-sequence run is below this
            threshold are demoted to structural. Default 1 disables the check
            and preserves pre-threshold behavior.

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
        if is_homopolymer_event(aligned_seqs, start_col, end_col, min_hp_length=min_hp_length):
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


def analyze_msa_columns_overlap_aware(aligned_seqs: List, min_overlap_bp: int,
                                       original_lengths: List[int],
                                       min_hp_length: int = 1) -> dict:
    """
    Analyze MSA columns, distinguishing terminal gaps from structural indels.

    Terminal gaps (from length differences at sequence ends) are NOT counted
    as structural indels when sequences have sufficient overlap in their
    shared region. This enables merging sequences from primer pools with
    different endpoints.

    Args:
        aligned_seqs: List of aligned sequences from SPOA
        min_overlap_bp: Minimum overlap required (0 to disable overlap mode)
        original_lengths: Original ungapped sequence lengths

    Returns dict with:
        'snp_count': SNPs in overlap region
        'structural_indel_count': Structural indels in overlap region only
        'structural_indel_length': Length of longest structural indel
        'homopolymer_indel_count': Homopolymer indels (anywhere)
        'homopolymer_indel_length': Length of longest homopolymer indel
        'terminal_gap_columns': Number of terminal gap columns (not counted as structural)
        'overlap_bp': Size of overlap region in base pairs
        'prefix_bp': Extension before overlap region (for logging)
        'suffix_bp': Extension after overlap region (for logging)
        'content_regions': List of (start, end) tuples per sequence (for span logging)
        'indel_count': Total events (backward compatibility)
        'max_indel_length': Max event length (backward compatibility)
    """
    alignment_length = len(aligned_seqs[0].seq)

    # Step 1: Find content region for each sequence (first non-gap to last non-gap)
    content_regions = []  # List of (start, end) tuples
    for seq in aligned_seqs:
        seq_str = str(seq.seq)
        # Find first and last non-gap positions
        first_base = next((i for i, c in enumerate(seq_str) if c != '-'), 0)
        last_base = alignment_length - 1 - next(
            (i for i, c in enumerate(reversed(seq_str)) if c != '-'), 0
        )
        content_regions.append((first_base, last_base))

    # Step 2: Calculate overlap region (intersection of all content regions)
    overlap_start = max(start for start, _ in content_regions)
    overlap_end = min(end for _, end in content_regions)

    # Calculate union region (for prefix/suffix extension reporting)
    union_start = min(start for start, _ in content_regions)
    union_end = max(end for _, end in content_regions)
    prefix_bp = overlap_start - union_start
    suffix_bp = union_end - overlap_end

    # Calculate actual overlap in base pairs (count only columns where all have bases)
    overlap_bp = 0
    if overlap_end >= overlap_start:
        for col_idx in range(overlap_start, overlap_end + 1):
            column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
            if all(c != '-' for c in column):
                overlap_bp += 1

    # Determine effective threshold for containment cases
    shorter_len = min(original_lengths)
    effective_threshold = min(min_overlap_bp, shorter_len)

    # Step 3: Count SNPs only within overlap region
    snp_count = 0
    for col_idx in range(overlap_start, overlap_end + 1):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
        unique_bases = set(c for c in column if c != '-')
        has_gap = '-' in column

        # SNP position: multiple different bases with NO gaps
        if len(unique_bases) > 1 and not has_gap:
            snp_count += 1

    # Step 4: Identify indel events, but only count those within overlap region
    indel_events = identify_indel_events(aligned_seqs, alignment_length)

    # Step 5: Classify each event and determine if it's in overlap region
    structural_events = []
    homopolymer_events = []
    terminal_gap_columns = 0

    for start_col, end_col in indel_events:
        # Check if this event is entirely within the overlap region
        is_in_overlap = (start_col >= overlap_start and end_col <= overlap_end)

        # Check if this is a terminal gap event (at the boundary of a content region)
        is_terminal = False
        for seq_start, seq_end in content_regions:
            # Terminal if event is adjacent to or outside a sequence's content region
            if end_col < seq_start or start_col > seq_end:
                is_terminal = True
                break
            # Also terminal if event is at the very edge of content
            if start_col == seq_start or end_col == seq_end:
                # Check if the gaps in this event are from this sequence's terminal
                for col_idx in range(start_col, end_col + 1):
                    column = [str(seq.seq[col_idx]) for seq in aligned_seqs]
                    for i, (s, e) in enumerate(content_regions):
                        if col_idx < s or col_idx > e:
                            if column[i] == '-':
                                is_terminal = True
                                break
                    if is_terminal:
                        break

        if is_terminal and overlap_bp >= effective_threshold:
            # Terminal gap from length difference - don't count as structural
            terminal_gap_columns += (end_col - start_col + 1)
        elif is_homopolymer_event(aligned_seqs, start_col, end_col, min_hp_length=min_hp_length):
            homopolymer_events.append((start_col, end_col))
        else:
            # Only count as structural if within overlap region
            if is_in_overlap:
                structural_events.append((start_col, end_col))
            else:
                # Outside overlap - this is a terminal gap
                terminal_gap_columns += (end_col - start_col + 1)

    # Step 6: Calculate statistics
    structural_indel_count = len(structural_events)
    homopolymer_indel_count = len(homopolymer_events)

    structural_indel_length = max((end - start + 1 for start, end in structural_events), default=0)
    homopolymer_indel_length = max((end - start + 1 for start, end in homopolymer_events), default=0)

    # Backward compatibility
    total_indel_count = structural_indel_count + homopolymer_indel_count
    max_indel_length = max(structural_indel_length, homopolymer_indel_length)

    return {
        'snp_count': snp_count,
        'structural_indel_count': structural_indel_count,
        'structural_indel_length': structural_indel_length,
        'homopolymer_indel_count': homopolymer_indel_count,
        'homopolymer_indel_length': homopolymer_indel_length,
        'terminal_gap_columns': terminal_gap_columns,
        'overlap_bp': overlap_bp,
        'prefix_bp': prefix_bp,
        'suffix_bp': suffix_bp,
        'content_regions': content_regions,
        'indel_count': total_indel_count,  # Backward compatibility
        'max_indel_length': max_indel_length  # Backward compatibility
    }
