"""Approach 1: per-cluster HP-length extraction (mode-based ground truth).

For each qualifying specimen's primary-anchor cluster, parse the SPOA MSA,
identify HP runs in the consensus, and extract the called HP length from
each read at each run. Aggregation and biological-variant filtering happen
in ``analyze.py``.

The MSA-walking helpers (``identify_hp_runs``, ``parse_msa_lightweight``,
``build_position_mappings``, ``find_anchor_columns``,
``extract_hp_length_from_read``, ``compute_mode``) are vendored from
``~/mm/analysis/hp_error_rate/hp_common.py``. They diverge intentionally
from ``speconsense.context.build_hp_runs`` (production), which defaults to
``min_length=2`` and operates on aligned-with-gaps coordinates — neither
of which fits the read-level HP-length-extraction use case.
"""

from __future__ import annotations

import glob
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import SeqIO

from .selection import SpecimenRecord


# ---------------------------------------------------------------------------
# MSA / HP-run utilities (vendored from hp_error_rate/hp_common.py)
# ---------------------------------------------------------------------------


def identify_hp_runs(sequence: str) -> List[Dict]:
    """Maximal HP runs in a sequence, including length-1 runs.

    Returns list of dicts: ``{base, start, end (exclusive), length,
    left_flank, right_flank}``.
    """
    if not sequence:
        return []

    runs = []
    run_start = 0
    run_base = sequence[0]
    for i in range(1, len(sequence) + 1):
        if i == len(sequence) or sequence[i] != run_base:
            left_flank = sequence[max(0, run_start - 2):run_start]
            right_flank = sequence[i:i + 2]
            runs.append({
                'base': run_base,
                'start': run_start,
                'end': i,
                'length': i - run_start,
                'left_flank': left_flank,
                'right_flank': right_flank,
            })
            if i < len(sequence):
                run_start = i
                run_base = sequence[i]
    return runs


def parse_msa_lightweight(msa_path: str) -> Tuple[str, str, List[str]]:
    """Parse SPOA MSA FASTA file.

    Returns ``(consensus_aligned, consensus_ungapped, read_aligned_sequences)``.
    Raises ``ValueError`` if no consensus record is found.
    """
    consensus_aligned: Optional[str] = None
    read_sequences: List[str] = []
    for record in SeqIO.parse(msa_path, 'fasta'):
        seq = str(record.seq).upper()
        if record.id == 'Consensus' or 'Consensus' in record.description:
            consensus_aligned = seq
        else:
            read_sequences.append(seq)
    if consensus_aligned is None:
        raise ValueError(f"No consensus sequence found in {msa_path}")
    return consensus_aligned, consensus_aligned.replace('-', ''), read_sequences


def build_position_mappings(consensus_aligned: str) -> Tuple[Dict[int, Optional[int]], Dict[int, int]]:
    """MSA column <-> consensus position bidirectional maps."""
    msa_to_cons: Dict[int, Optional[int]] = {}
    cons_to_msa: Dict[int, int] = {}
    cons_pos = 0
    for msa_col, ch in enumerate(consensus_aligned):
        if ch == '-':
            msa_to_cons[msa_col] = None
        else:
            msa_to_cons[msa_col] = cons_pos
            cons_to_msa[cons_pos] = msa_col
            cons_pos += 1
    return msa_to_cons, cons_to_msa


def find_anchor_columns(hp_run: Dict, cons_to_msa: Dict[int, int], msa_length: int) -> Tuple[int, int]:
    """MSA columns for the non-HP bases flanking ``hp_run``.

    Returns ``(left_anchor, right_anchor)``. ``-1`` denotes "no left anchor"
    (run at sequence start); ``msa_length`` denotes "no right anchor" (run
    at end).
    """
    start = hp_run['start']
    end = hp_run['end']
    left_anchor = cons_to_msa[start - 1] if start > 0 else -1
    right_anchor = cons_to_msa[end] if end in cons_to_msa else msa_length
    return left_anchor, right_anchor


def extract_hp_length_from_read(
    read_aligned: str,
    hp_base: str,
    left_anchor_col: int,
    right_anchor_col: int,
) -> Tuple[int, str]:
    """Called HP length for one read at one HP run.

    Status values: ``clean`` (all bases in region are hp_base), ``complex``
    (HP base interspersed with others — likely substitution), ``full_deletion``
    (no bases at all in region), ``no_coverage`` (read doesn't span this
    region).
    """
    region_start = left_anchor_col + 1
    region_end = right_anchor_col
    if region_start >= region_end:
        return 0, 'no_coverage'

    bases_in_region: List[str] = []
    all_gaps = True
    for col in range(region_start, region_end):
        base = read_aligned[col]
        if base != '-':
            bases_in_region.append(base)
            all_gaps = False

    left_ok = left_anchor_col < 0 or read_aligned[left_anchor_col] != '-'
    right_ok = right_anchor_col >= len(read_aligned) or read_aligned[right_anchor_col] != '-'
    if not left_ok and not right_ok and all_gaps:
        return 0, 'no_coverage'

    if not bases_in_region:
        return 0, 'full_deletion'

    hp_count = sum(1 for b in bases_in_region if b == hp_base)
    non_hp_count = len(bases_in_region) - hp_count
    if non_hp_count == 0:
        return hp_count, 'clean'
    return hp_count, 'complex'


def compute_mode(values: List[int]) -> int:
    """Mode of an integer list; ties resolve to the smaller value."""
    if not values:
        return 0
    counts: Dict[int, int] = {}
    for v in values:
        counts[v] = counts.get(v, 0) + 1
    max_count = max(counts.values())
    return min(v for v, c in counts.items() if c == max_count)


# ---------------------------------------------------------------------------
# Per-specimen extraction driver
# ---------------------------------------------------------------------------


@dataclass
class HPObservation:
    """One HP run's observations across all reads in its cluster."""
    specimen: str
    consensus_pos: int
    base: str
    spoa_length: int
    mode_length: int
    called_lengths: np.ndarray  # int16 array of read-level called lengths
    n_reads: int
    n_complex: int
    n_excluded: int


@dataclass
class NonHPCounts:
    """Per-cluster non-HP error counts (used by both approach 1's calibration
    output and approach 2's pooled aggregation)."""
    matches: int = 0
    substitutions: int = 0
    deletions: int = 0
    insertions: int = 0
    total_positions: int = 0


def compute_non_hp_error_counts(
    read_sequences: List[str],
    consensus_aligned: str,
    hp_positions_set: set,
) -> NonHPCounts:
    """Per-position error counts at non-HP consensus columns.

    HP positions (columns inside a length>1 run in ungapped consensus) are
    skipped. Insertion columns (consensus='-') contribute insertion counts
    only. Returned counts pool across all read sequences.
    """
    stats = NonHPCounts()
    _msa_to_cons, _ = build_position_mappings(consensus_aligned)

    for col, cons_base in enumerate(consensus_aligned):
        if cons_base == '-':
            for read_seq in read_sequences:
                if col < len(read_seq) and read_seq[col] != '-':
                    stats.insertions += 1
            continue

        cons_pos = _msa_to_cons[col]
        if cons_pos in hp_positions_set:
            continue

        stats.total_positions += 1
        for read_seq in read_sequences:
            if col >= len(read_seq):
                continue
            read_base = read_seq[col]
            if read_base == '-':
                stats.deletions += 1
            elif read_base == cons_base:
                stats.matches += 1
            else:
                stats.substitutions += 1

    return stats


def hp_positions_of_runs(hp_runs: List[Dict]) -> set:
    """Set of ungapped-consensus positions inside any HP run of length > 1."""
    out = set()
    for run in hp_runs:
        if run['length'] > 1:
            for pos in range(run['start'], run['end']):
                out.add(pos)
    return out


def _find_primary_msa(specimen_name: str, debug_dir: str) -> Optional[str]:
    """Locate the primary-anchor MSA file for a specimen.

    Tries the current schema-2.0 naming (``-1.v1-RiC*-msa.fasta``) first,
    falls back to the legacy ``-c1-RiC*-msa.fasta``. When multiple matches
    exist (cluster ID reused across pipeline passes), picks the highest
    RiC.
    """
    base = os.path.join(debug_dir, specimen_name)
    matches = glob.glob(f"{base}-1.v1-RiC*-msa.fasta")
    if not matches:
        matches = glob.glob(f"{base}-c1-RiC*-msa.fasta")
    if not matches:
        return None
    if len(matches) > 1:
        matches.sort(
            key=lambda p: int(p.split('-RiC')[1].split('-')[0]),
            reverse=True,
        )
    return matches[0]


def process_specimen(
    specimen_name: str,
    debug_dir: str,
) -> Tuple[Optional[List[HPObservation]], Optional[NonHPCounts]]:
    """Extract HP observations + non-HP error counts for one specimen.

    Returns ``(None, None)`` when the primary anchor's MSA cannot be located
    or has no reads.
    """
    msa_path = _find_primary_msa(specimen_name, debug_dir)
    if msa_path is None:
        return None, None

    try:
        consensus_aligned, consensus_ungapped, read_sequences = parse_msa_lightweight(msa_path)
    except Exception as e:  # noqa: BLE001
        logging.warning(f"Failed to parse {os.path.basename(msa_path)}: {e}")
        return None, None
    if not read_sequences:
        return None, None

    msa_length = len(consensus_aligned)
    _, cons_to_msa = build_position_mappings(consensus_aligned)
    hp_runs = identify_hp_runs(consensus_ungapped)
    hp_pos_set = hp_positions_of_runs(hp_runs)

    observations: List[HPObservation] = []
    for run in hp_runs:
        left, right = find_anchor_columns(run, cons_to_msa, msa_length)
        called: List[int] = []
        n_complex = 0
        n_excluded = 0
        for read_seq in read_sequences:
            length, status = extract_hp_length_from_read(read_seq, run['base'], left, right)
            if status == 'clean':
                called.append(length)
            elif status == 'complex':
                n_complex += 1
                called.append(length)
            else:  # full_deletion / no_coverage
                n_excluded += 1
        if not called:
            continue
        observations.append(HPObservation(
            specimen=specimen_name,
            consensus_pos=run['start'],
            base=run['base'],
            spoa_length=run['length'],
            mode_length=compute_mode(called),
            called_lengths=np.array(called, dtype=np.int16),
            n_reads=len(called),
            n_complex=n_complex,
            n_excluded=n_excluded,
        ))

    non_hp = compute_non_hp_error_counts(read_sequences, consensus_aligned, hp_pos_set)
    return observations, non_hp


# ---------------------------------------------------------------------------
# Parallel driver
# ---------------------------------------------------------------------------


def _worker(args: Tuple[str, str]):
    name, debug_dir = args
    try:
        return name, process_specimen(name, debug_dir)
    except Exception as e:  # noqa: BLE001
        return name, (None, None, str(e))  # third element flags worker error


def extract_all(
    qualifying: List[SpecimenRecord],
    input_dir: str,
    threads: int = 1,
    progress: bool = True,
) -> Tuple[List[HPObservation], NonHPCounts, int]:
    """Run approach-1 extraction over all qualifying specimens.

    Returns ``(all_observations, calibration_non_hp_counts, n_skipped)``.
    The calibration non-HP totals (pooled over primary-anchor clusters
    only) are a useful sanity-check sidecar; the canonical non-HP rates
    used in the output model come from ``nonhp.py`` (approach 2, all
    clusters).

    Parallelizes via ``ProcessPoolExecutor`` when ``threads > 1``.
    """
    debug_dir = os.path.join(input_dir, "cluster_debug")
    all_obs: List[HPObservation] = []
    calibration = NonHPCounts()
    n_skipped = 0

    try:
        from tqdm import tqdm
    except ImportError:
        def tqdm(it, **_kw):
            return it

    work = [(r.name, debug_dir) for r in qualifying]
    iterator = tqdm(work, desc="Approach 1 (HP runs)", disable=not progress)

    if threads <= 1:
        for name, debug in iterator:
            obs, non_hp = process_specimen(name, debug)
            if obs is None:
                n_skipped += 1
                continue
            all_obs.extend(obs)
            calibration.matches += non_hp.matches
            calibration.substitutions += non_hp.substitutions
            calibration.deletions += non_hp.deletions
            calibration.insertions += non_hp.insertions
            calibration.total_positions += non_hp.total_positions
        return all_obs, calibration, n_skipped

    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(process_specimen, name, debug) for name, debug in work]
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc="Approach 1 (HP runs)", disable=not progress):
            obs, non_hp = fut.result()
            if obs is None:
                n_skipped += 1
                continue
            all_obs.extend(obs)
            calibration.matches += non_hp.matches
            calibration.substitutions += non_hp.substitutions
            calibration.deletions += non_hp.deletions
            calibration.insertions += non_hp.insertions
            calibration.total_positions += non_hp.total_positions
    return all_obs, calibration, n_skipped
