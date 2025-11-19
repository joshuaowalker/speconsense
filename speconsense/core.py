#!/usr/bin/env python3

from collections import defaultdict
import argparse
import json
import logging
import os
import random
import statistics
import subprocess
import sys
import tempfile
import math
import itertools
from typing import List, Set, Tuple, Optional, Dict, Any, NamedTuple
from datetime import datetime

import edlib
from adjusted_identity import score_alignment, AdjustmentParams
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
from tqdm import tqdm

try:
    from speconsense import __version__
except ImportError:
    # Fallback for when running as a script directly (e.g., in tests)
    __version__ = "dev"

# IUPAC nucleotide ambiguity codes mapping
# Maps sets of nucleotides to their corresponding IUPAC code
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


class ErrorPosition(NamedTuple):
    """An error at a specific position in the MSA."""
    msa_position: int  # 0-indexed position in MSA alignment
    error_type: str  # 'sub', 'ins', or 'del'


class ReadAlignment(NamedTuple):
    """Alignment result for a single read against consensus."""
    read_id: str
    aligned_sequence: str  # Gapped sequence from MSA
    read_length: int

    # Raw metrics (count all differences including homopolymer length)
    edit_distance: int
    num_insertions: int
    num_deletions: int
    num_substitutions: int
    error_positions: List[ErrorPosition]  # Detailed error information

    # Homopolymer-normalized metrics (exclude homopolymer extensions)
    normalized_edit_distance: int  # Edit distance excluding homopolymer length differences
    normalized_error_positions: List[ErrorPosition]  # Only non-homopolymer errors
    score_aligned: str  # Scoring string from adjusted-identity ('|'=match, '='=homopolymer, ' '=error)


class PositionStats(NamedTuple):
    """Statistics for a single position in the MSA."""
    msa_position: int  # Position in MSA (0-indexed)
    consensus_position: Optional[int]  # Position in consensus (None for insertion columns)
    coverage: int
    error_count: int
    error_rate: float
    sub_count: int
    ins_count: int
    del_count: int
    is_homopolymer: bool
    consensus_nucleotide: str  # Base in consensus at this MSA position (or '-' for insertion)
    base_composition: Dict[str, int]  # {A: 50, C: 3, G: 45, T: 2, '-': 0}


class MSAResult(NamedTuple):
    """Result from SPOA multiple sequence alignment.

    Attributes:
        consensus: Ungapped consensus sequence
        msa_string: Raw MSA in FASTA format (for file writing)
        alignments: Parsed read alignments with gapped sequences
        msa_to_consensus_pos: Mapping from MSA position to consensus position
    """
    consensus: str
    msa_string: str
    alignments: List[ReadAlignment]
    msa_to_consensus_pos: Dict[int, Optional[int]]


# ============================================================================
# MSA Analysis Functions (moved from analyze.py to break circular dependency)
# ============================================================================

def parse_score_aligned_for_errors(
    score_aligned: str,
    read_aligned: str,
    consensus_aligned: str
) -> List[ErrorPosition]:
    """
    Parse score_aligned string to extract non-homopolymer errors.

    The score_aligned string from adjusted-identity uses these codes:
    - '|' : Exact match (not an error)
    - '=' : Ambiguous match or homopolymer extension (not counted as error)
    - ' ' : Substitution or indel (IS an error)
    - '.' : End-trimmed position (not counted)

    Args:
        score_aligned: Scoring string from adjusted-identity
        read_aligned: Aligned read sequence with gaps
        consensus_aligned: Aligned consensus sequence with gaps

    Returns:
        List of ErrorPosition for positions marked as errors (excluding homopolymer extensions)
    """
    normalized_errors = []

    for msa_pos, (score_char, read_base, cons_base) in enumerate(
        zip(score_aligned, read_aligned, consensus_aligned)
    ):
        # Skip matches and homopolymer extensions
        if score_char in ('|', '=', '.'):
            continue

        # This is a real error (substitution or indel) - classify it
        if read_base == '-' and cons_base != '-':
            error_type = 'del'
        elif read_base != '-' and cons_base == '-':
            error_type = 'ins'
        elif read_base != cons_base:
            error_type = 'sub'
        else:
            # Both are gaps or identical - should not happen if score_char indicates error
            continue

        normalized_errors.append(ErrorPosition(msa_pos, error_type))

    return normalized_errors


def extract_alignments_from_msa(
    msa_string: str,
    enable_homopolymer_normalization: bool = True
) -> Tuple[List[ReadAlignment], str, Dict[int, Optional[int]]]:
    """
    Extract read alignments from an MSA string with optional homopolymer normalization.

    The MSA contains aligned sequences where the consensus has header containing "Consensus".
    This function compares each read to the consensus at each aligned position.

    Error classification (raw metrics):
    - Both '-': Not an error (read doesn't cover this position)
    - Read '-', consensus base: Deletion (missing base in read)
    - Read base, consensus '-': Insertion (extra base in read)
    - Different bases: Substitution
    - Same base: Match (not an error)

    When enable_homopolymer_normalization=True, also computes normalized metrics that
    exclude homopolymer length differences using adjusted-identity library.

    IMPORTANT: Errors are reported at MSA positions, not consensus positions.
    This avoids ambiguity when multiple insertion columns map to the same consensus position.

    Args:
        msa_string: MSA content in FASTA format
        enable_homopolymer_normalization: If True, compute homopolymer-normalized metrics

    Returns:
        Tuple of:
        - list of ReadAlignment objects (with both raw and normalized metrics)
        - consensus sequence without gaps
        - mapping from MSA position to consensus position (None for insertion columns)
    """
    from io import StringIO

    # Define adjustment parameters for homopolymer normalization
    # Only normalize homopolymers (single-base repeats), no other adjustments
    HOMOPOLYMER_ADJUSTMENT_PARAMS = AdjustmentParams(
        normalize_homopolymers=True,
        handle_iupac_overlap=False,
        normalize_indels=False,
        end_skip_distance=0,
        max_repeat_motif_length=1  # Single-base repeats only
    )

    # Parse MSA
    msa_handle = StringIO(msa_string)
    records = list(SeqIO.parse(msa_handle, 'fasta'))

    if not records:
        logging.warning("No sequences found in MSA string")
        return [], "", {}

    # Find consensus sequence
    consensus_record = None
    read_records = []

    for record in records:
        if 'Consensus' in record.description or 'Consensus' in record.id:
            consensus_record = record
        else:
            read_records.append(record)

    if consensus_record is None:
        logging.warning("No consensus sequence found in MSA string")
        return [], "", {}

    consensus_aligned = str(consensus_record.seq).upper()
    msa_length = len(consensus_aligned)

    # Build mapping from MSA position to consensus position (excluding gaps)
    # For insertion columns (consensus has '-'), maps to None
    msa_to_consensus_pos = {}
    consensus_pos = 0
    for msa_pos in range(msa_length):
        if consensus_aligned[msa_pos] != '-':
            msa_to_consensus_pos[msa_pos] = consensus_pos
            consensus_pos += 1
        else:
            # Insertion column - no consensus position
            msa_to_consensus_pos[msa_pos] = None

    # Get consensus without gaps for return value
    consensus_ungapped = consensus_aligned.replace('-', '')

    # Process each read
    alignments = []

    for read_record in read_records:
        read_aligned = str(read_record.seq).upper()

        if len(read_aligned) != msa_length:
            logging.warning(f"Read {read_record.id} length mismatch with MSA length")
            continue

        # Compare read to consensus at each position
        error_positions = []
        num_insertions = 0
        num_deletions = 0
        num_substitutions = 0

        for msa_pos in range(msa_length):
            read_base = read_aligned[msa_pos]
            cons_base = consensus_aligned[msa_pos]

            # Skip if both are gaps (read doesn't cover this position)
            if read_base == '-' and cons_base == '-':
                continue

            # Classify error type and record at MSA position
            if read_base == '-' and cons_base != '-':
                # Deletion (missing base in read)
                error_positions.append(ErrorPosition(msa_pos, 'del'))
                num_deletions += 1
            elif read_base != '-' and cons_base == '-':
                # Insertion (extra base in read)
                error_positions.append(ErrorPosition(msa_pos, 'ins'))
                num_insertions += 1
            elif read_base != cons_base:
                # Substitution (different bases)
                error_positions.append(ErrorPosition(msa_pos, 'sub'))
                num_substitutions += 1
            # else: match, no error

        # Calculate edit distance and read length
        edit_distance = num_insertions + num_deletions + num_substitutions
        read_length = len(read_aligned.replace('-', ''))  # Length without gaps

        # Compute homopolymer-normalized metrics if enabled
        if enable_homopolymer_normalization:
            try:
                # Use adjusted-identity to get homopolymer-normalized scoring
                result = score_alignment(
                    read_aligned,
                    consensus_aligned,
                    HOMOPOLYMER_ADJUSTMENT_PARAMS
                )

                # Parse score_aligned string to extract normalized errors
                normalized_error_positions = parse_score_aligned_for_errors(
                    result.score_aligned,
                    read_aligned,
                    consensus_aligned
                )

                normalized_edit_distance = result.mismatches
                score_aligned_str = result.score_aligned

            except Exception as e:
                # If normalization fails, fall back to raw metrics
                logging.warning(f"Homopolymer normalization failed for read {read_record.id}: {e}")
                normalized_edit_distance = edit_distance
                normalized_error_positions = error_positions
                score_aligned_str = ""
        else:
            # Homopolymer normalization disabled - use raw metrics
            normalized_edit_distance = edit_distance
            normalized_error_positions = error_positions
            score_aligned_str = ""

        # Create alignment object with both raw and normalized metrics
        alignment = ReadAlignment(
            read_id=read_record.id,
            aligned_sequence=read_aligned,  # Store gapped sequence
            read_length=read_length,
            # Raw metrics
            edit_distance=edit_distance,
            num_insertions=num_insertions,
            num_deletions=num_deletions,
            num_substitutions=num_substitutions,
            error_positions=error_positions,
            # Normalized metrics
            normalized_edit_distance=normalized_edit_distance,
            normalized_error_positions=normalized_error_positions,
            score_aligned=score_aligned_str
        )
        alignments.append(alignment)

    return alignments, consensus_ungapped, msa_to_consensus_pos


def detect_homopolymer_run(seq: str, position: int, min_length: int = 4) -> bool:
    """
    Check if a position is within a homopolymer run.

    Args:
        seq: DNA sequence
        position: Position to check (0-indexed)
        min_length: Minimum length of homopolymer to detect

    Returns:
        True if position is within a homopolymer run of min_length or longer
    """
    if position >= len(seq) or position < 0:
        return False

    base = seq[position]

    # Find start of run
    start = position
    while start > 0 and seq[start - 1] == base:
        start -= 1

    # Find end of run
    end = position
    while end < len(seq) - 1 and seq[end + 1] == base:
        end += 1

    run_length = end - start + 1
    return run_length >= min_length


def analyze_positional_variation(
    alignments: List[ReadAlignment],
    consensus_seq: str,
    consensus_aligned: str,
    msa_to_consensus_pos: Dict[int, Optional[int]],
    overall_error_rate: float
) -> List[PositionStats]:
    """
    Analyze error rates at each position in the MSA.

    Identifies positions with elevated error rates, which may indicate
    homopolymer regions, structural variants, or other error-prone
    sequence contexts.

    IMPORTANT: All analysis is performed in MSA space (not consensus space).
    This correctly handles insertion columns where multiple MSA positions
    don't correspond to any consensus position.

    Args:
        alignments: List of read alignments (with errors at MSA positions)
        consensus_seq: Consensus sequence (ungapped)
        consensus_aligned: Consensus sequence (gapped, from MSA)
        msa_to_consensus_pos: Mapping from MSA position to consensus position
        overall_error_rate: Overall error rate (currently unused, kept for compatibility)

    Returns:
        List of PositionStats for each MSA position
    """
    msa_length = len(consensus_aligned)

    # Build error frequency matrix in MSA space
    # For each MSA position: [sub_count, ins_count, del_count, total_coverage]
    error_matrix = np.zeros((msa_length, 4), dtype=int)

    # Build base composition matrix in MSA space
    base_composition_matrix = [
        {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
        for _ in range(msa_length)
    ]

    # Process alignments to count errors at MSA positions
    for read_idx, alignment in enumerate(alignments):
        # Count this read as coverage for all MSA positions
        # Note: alignments span the full MSA
        for msa_pos in range(msa_length):
            error_matrix[msa_pos, 3] += 1  # coverage

        # Add errors at specific MSA positions
        for error_pos in alignment.error_positions:
            msa_pos = error_pos.msa_position
            if 0 <= msa_pos < msa_length:
                if error_pos.error_type == 'sub':
                    error_matrix[msa_pos, 0] += 1
                elif error_pos.error_type == 'ins':
                    error_matrix[msa_pos, 1] += 1
                elif error_pos.error_type == 'del':
                    error_matrix[msa_pos, 2] += 1

        # Extract base composition from aligned sequence
        read_aligned = alignment.aligned_sequence
        if len(read_aligned) != msa_length:
            continue

        # Track what base each read has at each MSA position
        for msa_pos in range(msa_length):
            read_base = read_aligned[msa_pos]

            # Normalize base
            if read_base in ['A', 'C', 'G', 'T', '-']:
                base_composition_matrix[msa_pos][read_base] += 1
            else:
                # Treat N or other ambiguous as gap
                base_composition_matrix[msa_pos]['-'] += 1

    # Calculate statistics for each MSA position
    position_stats = []

    for msa_pos in range(msa_length):
        sub_count = error_matrix[msa_pos, 0]
        ins_count = error_matrix[msa_pos, 1]
        del_count = error_matrix[msa_pos, 2]
        coverage = error_matrix[msa_pos, 3]

        # Total error events
        error_count = sub_count + ins_count + del_count
        error_rate = error_count / coverage if coverage > 0 else 0.0

        # Get consensus position (None for insertion columns)
        cons_pos = msa_to_consensus_pos[msa_pos]

        # Get consensus nucleotide at this MSA position
        cons_nucleotide = consensus_aligned[msa_pos]

        # Sequence context - calculate based on consensus position if available
        if cons_pos is not None:
            is_homopolymer = detect_homopolymer_run(consensus_seq, cons_pos, min_length=4)
        else:
            # Insertion column - no consensus context
            is_homopolymer = False

        # Get base composition for this MSA position
        base_comp = base_composition_matrix[msa_pos].copy()

        position_stats.append(PositionStats(
            msa_position=msa_pos,
            consensus_position=cons_pos,
            coverage=coverage,
            error_count=error_count,
            error_rate=error_rate,
            sub_count=sub_count,
            ins_count=ins_count,
            del_count=del_count,
            is_homopolymer=is_homopolymer,
            consensus_nucleotide=cons_nucleotide,
            base_composition=base_comp
        ))

    return position_stats


def is_variant_position_with_composition(
    position_stats: PositionStats,
    min_variant_frequency: float = 0.20,
    min_variant_count: int = 5
) -> Tuple[bool, List[str], str]:
    """
    Identify variant positions using simple frequency and count thresholds.

    This function determines if a position shows systematic variation (true biological
    variant) rather than scattered sequencing errors.

    Criteria for variant detection:
    1. At least one alternative allele must have frequency ≥ min_variant_frequency
    2. That allele must have count ≥ min_variant_count

    Args:
        position_stats: Position statistics including base composition
        min_variant_frequency: Minimum alternative allele frequency (default: 0.20 for 20%)
        min_variant_count: Minimum alternative allele read count (default: 5 reads)

    Returns:
        Tuple of (is_variant, variant_bases, reason)
        - is_variant: True if this position requires cluster separation
        - variant_bases: List of alternative bases meeting criteria (e.g., ['G', 'T'])
        - reason: Explanation of decision for logging/debugging
    """
    n = position_stats.coverage
    base_composition = position_stats.base_composition

    # Check we have composition data
    if not base_composition or sum(base_composition.values()) == 0:
        return False, [], "No composition data available"

    sorted_bases = sorted(
        base_composition.items(),
        key=lambda x: x[1],
        reverse=True
    )

    if len(sorted_bases) < 2:
        return False, [], "No alternative alleles observed"

    # Check each alternative allele (skip consensus base at index 0)
    variant_bases = []
    variant_details = []

    for base, count in sorted_bases[1:]:
        freq = count / n if n > 0 else 0

        # Must meet both frequency and count thresholds
        if freq >= min_variant_frequency and count >= min_variant_count:
            variant_bases.append(base)
            variant_details.append(f"{base}:{count}/{n}({freq:.1%})")

    if variant_bases:
        return True, variant_bases, f"Variant alleles: {', '.join(variant_details)}"

    return False, [], "No variants detected"


class SpecimenClusterer:
    def __init__(self, min_identity: float = 0.9,
                 inflation: float = 4.0,
                 min_size: int = 5,
                 min_cluster_ratio: float = 0.2,
                 max_sample_size: int = 500,
                 presample_size: int = 1000,
                 k_nearest_neighbors: int = 20,
                 sample_name: str = "sample",
                 disable_homopolymer_equivalence: bool = False,
                 output_dir: str = "clusters",
                 outlier_identity_threshold: Optional[float] = None,
                 enable_secondpass_phasing: bool = True,
                 min_variant_frequency: float = 0.20,
                 min_variant_count: int = 5):
        self.min_identity = min_identity
        self.inflation = inflation
        self.min_size = min_size
        self.min_cluster_ratio = min_cluster_ratio
        self.max_sample_size = max_sample_size
        self.presample_size = presample_size
        self.k_nearest_neighbors = k_nearest_neighbors
        self.sample_name = sample_name
        self.disable_homopolymer_equivalence = disable_homopolymer_equivalence
        self.output_dir = output_dir

        # Auto-calculate outlier identity threshold if not provided
        # Logic: min_identity accounts for 2×error (read-to-read comparison)
        # outlier_identity_threshold accounts for 1×error (read-to-consensus)
        # Therefore: outlier_identity_threshold = (1 + min_identity) / 2
        if outlier_identity_threshold is None:
            self.outlier_identity_threshold = (1.0 + min_identity) / 2.0
        else:
            self.outlier_identity_threshold = outlier_identity_threshold

        self.enable_secondpass_phasing = enable_secondpass_phasing
        self.min_variant_frequency = min_variant_frequency
        self.min_variant_count = min_variant_count
        self.sequences = {}  # id -> sequence string
        self.records = {}  # id -> SeqRecord object
        self.id_map = {}  # short_id -> original_id
        self.rev_id_map = {}  # original_id -> short_id

        # Create output directory and debug subdirectory
        os.makedirs(self.output_dir, exist_ok=True)
        self.debug_dir = os.path.join(self.output_dir, "cluster_debug")
        os.makedirs(self.debug_dir, exist_ok=True)

        # Initialize attributes that may be set later
        self.input_file = None
        self.augment_input = None
        self.algorithm = None
        self.orient_mode = None
        self.primers_file = None

    def write_metadata(self) -> None:
        """Write run metadata to JSON file for use by post-processing tools."""
        metadata = {
            "version": __version__,
            "timestamp": datetime.now().isoformat(),
            "sample_name": self.sample_name,
            "parameters": {
                "algorithm": self.algorithm,
                "min_identity": self.min_identity,
                "inflation": self.inflation,
                "min_size": self.min_size,
                "min_cluster_ratio": self.min_cluster_ratio,
                "max_sample_size": self.max_sample_size,
                "presample_size": self.presample_size,
                "k_nearest_neighbors": self.k_nearest_neighbors,
                "disable_homopolymer_equivalence": self.disable_homopolymer_equivalence,
                "outlier_identity_threshold": self.outlier_identity_threshold,
                "enable_secondpass_phasing": self.enable_secondpass_phasing,
                "min_variant_frequency": self.min_variant_frequency,
                "min_variant_count": self.min_variant_count,
                "orient_mode": self.orient_mode,
            },
            "input_file": self.input_file,
            "augment_input": self.augment_input,
        }

        # Add primer information if loaded
        if hasattr(self, 'primers') and self.primers:
            metadata["primers_file"] = self.primers_file
            metadata["primers"] = {}

            # Store primer sequences (avoid duplicates from RC versions)
            seen_primers = set()
            for primer_name, primer_seq in self.primers:
                # Skip RC versions (they end with _RC)
                if not primer_name.endswith('_RC') and primer_name not in seen_primers:
                    metadata["primers"][primer_name] = primer_seq
                    seen_primers.add(primer_name)

        # Write metadata file
        metadata_file = os.path.join(self.debug_dir, f"{self.sample_name}-metadata.json")
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logging.info(f"Wrote run metadata to {metadata_file}")

    def write_phasing_stats(self, initial_clusters_count: int, subclusters_count: int,
                           merged_count: int, final_count: int) -> None:
        """Write phasing statistics to JSON file after clustering completes.

        Args:
            initial_clusters_count: Number of clusters from initial clustering
            subclusters_count: Number of sub-clusters after phasing
            merged_count: Number of clusters after merging
            final_count: Number of final clusters after filtering
        """
        phasing_stats = {
            "phasing_enabled": self.enable_secondpass_phasing,
            "initial_clusters": initial_clusters_count,
            "phased_subclusters": subclusters_count,
            "after_merging": merged_count,
            "after_filtering": final_count,
            "clusters_split": subclusters_count > initial_clusters_count,
            "clusters_merged": merged_count < subclusters_count,
            "net_change": final_count - initial_clusters_count
        }

        # Write phasing stats to separate JSON file
        stats_file = os.path.join(self.debug_dir, f"{self.sample_name}-phasing_stats.json")
        with open(stats_file, 'w') as f:
            json.dump(phasing_stats, f, indent=2)

        logging.info(f"Wrote phasing statistics to {stats_file}")

    def add_sequences(self, records: List[SeqIO.SeqRecord],
                      augment_records: Optional[List[SeqIO.SeqRecord]] = None) -> None:
        """Add sequences to be clustered, with optional presampling."""
        all_records = records.copy()  # Start with primary records

        # Track the source of each record for potential logging/debugging
        primary_count = len(records)
        augment_count = 0

        # Add augmented records if provided
        if augment_records:
            augment_count = len(augment_records)
            all_records.extend(augment_records)

        if self.presample_size and len(all_records) > self.presample_size:
            logging.info(f"Presampling {self.presample_size} sequences from {len(all_records)} total "
                         f"({primary_count} primary, {augment_count} augmented)")

            # First, sort primary sequences by quality and take as many as possible
            primary_sorted = sorted(
                records,
                key=lambda r: -statistics.mean(r.letter_annotations["phred_quality"])
            )

            # Determine how many primary sequences to include (all if possible)
            primary_to_include = min(len(primary_sorted), self.presample_size)
            presampled = primary_sorted[:primary_to_include]

            # If we still have room, add augmented sequences sorted by quality
            remaining_slots = self.presample_size - primary_to_include
            if remaining_slots > 0 and augment_records:
                augment_sorted = sorted(
                    augment_records,
                    key=lambda r: -statistics.mean(r.letter_annotations["phred_quality"])
                )
                presampled.extend(augment_sorted[:remaining_slots])

            logging.info(f"Presampled {len(presampled)} sequences "
                         f"({primary_to_include} primary, {len(presampled) - primary_to_include} augmented)")
            all_records = presampled

        # Add all selected records to internal storage
        for record in all_records:
            self.sequences[record.id] = str(record.seq)
            self.records[record.id] = record

    def write_mcl_input(self, output_file: str) -> None:
        """Write similarity matrix in MCL input format using k-nearest neighbors approach."""
        self._create_id_mapping()

        # Calculate all pairwise similarities
        similarities = {}
        n = len(self.sequences)
        total_comparisons = (n * (n - 1)) // 2

        with tqdm(total=total_comparisons, desc="Calculating pairwise sequence similarities") as pbar:
            for id1 in self.sequences:
                similarities[id1] = {}
                for id2 in self.sequences:
                    if id1 >= id2:  # Only calculate upper triangle
                        continue
                    sim = self.calculate_similarity(self.sequences[id1], self.sequences[id2])
                    similarities[id1][id2] = sim
                    similarities.setdefault(id2, {})[id1] = sim  # Mirror for easy lookup
                    pbar.update(1)

        # Use k-nearest neighbors for each sequence to create sparse graph
        k = min(self.k_nearest_neighbors, n - 1)  # Connect to at most k neighbors
        min_edges_per_node = 3  # Ensure at least some connections per node

        with open(output_file, 'w') as f:
            for id1 in self.sequences:
                short_id1 = self.rev_id_map[id1]

                # Sort neighbors by similarity
                neighbors = sorted(
                    [(id2, sim) for id2, sim in similarities[id1].items()],
                    key=lambda x: x[1], reverse=True
                )

                # Take top k neighbors with sufficient similarity
                top_neighbors = [
                    (id2, sim) for id2, sim in neighbors[:k]
                    if sim >= self.min_identity
                ]

                # Ensure at least min_edges connections if possible
                if len(top_neighbors) < min_edges_per_node and len(neighbors) >= min_edges_per_node:
                    additional_needed = min_edges_per_node - len(top_neighbors)
                    for id2, sim in neighbors[k:k + additional_needed]:
                        if sim < self.min_identity * 0.9:  # Allow slightly lower threshold
                            break
                        top_neighbors.append((id2, sim))

                # Write edges
                for id2, sim in top_neighbors:
                    short_id2 = self.rev_id_map[id2]
                    # Transform similarity to emphasize differences
                    transformed_sim = sim ** 2  # Square the similarity
                    f.write(f"{short_id1}\t{short_id2}\t{transformed_sim:.6f}\n")

    def run_mcl(self, input_file: str, output_file: str) -> None:
        """Run MCL clustering algorithm with optimized parameters."""
        cmd = [
            "mcl",
            input_file,
            "--abc",  # Input is in ABC format (node1 node2 weight)
            "-I", str(self.inflation),  # Inflation parameter
            "-scheme", "7",  # More advanced flow simulation
            "-pct", "50",  # Prune weakest 50% of connections during iterations
            "-o", output_file  # Output file
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            logging.debug(f"MCL stdout: {result.stdout}")
            logging.debug(f"MCL stderr: {result.stderr}")

        except subprocess.CalledProcessError as e:
            logging.error(f"MCL failed with return code {e.returncode}")
            logging.error(f"Command: {' '.join(cmd)}")
            logging.error(f"Stderr: {e.stderr}")
            raise

    def merge_similar_clusters(self, clusters: List[Dict]) -> List[Dict]:
        """
        Merge clusters whose consensus sequences are identical or homopolymer-equivalent.
        Preserves provenance metadata through the merging process.

        Note: Primer trimming is performed before comparison to ensure clusters that differ
        only in primer regions are properly merged. Trimmed consensuses are used only for
        comparison and are discarded after merging.

        Args:
            clusters: List of cluster dictionaries with 'read_ids' and provenance fields

        Returns:
            List of merged cluster dictionaries with combined provenance
        """
        if not clusters:
            return []

        # Sort clusters by size, largest first
        clusters = sorted(clusters, key=lambda c: len(c['read_ids']), reverse=True)

        # Generate a consensus sequence for each cluster
        logging.info("Generating consensus sequences for cluster merging...")
        consensuses = []
        cluster_to_consensus = {}  # Map from cluster index to its consensus

        with tqdm(total=len(clusters), desc="Generating cluster consensuses") as pbar:
            for i, cluster_dict in enumerate(clusters):
                cluster_reads = cluster_dict['read_ids']

                # Sample from larger clusters to speed up consensus generation
                if len(cluster_reads) > self.max_sample_size:
                    # Sample by quality
                    qualities = []
                    for seq_id in cluster_reads:
                        record = self.records[seq_id]
                        mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                        qualities.append((mean_quality, seq_id))

                    sampled_ids = [seq_id for _, seq_id in
                                   sorted(qualities, reverse=True)[:self.max_sample_size]]
                    sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sampled_ids}
                else:
                    sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in cluster_reads}

                result = self.run_spoa(sampled_seqs)

                # Skip empty clusters (can occur when all sequences are filtered out)
                if result is None:
                    logging.warning(f"Cluster {i} produced no consensus (empty cluster), skipping")
                    pbar.update(1)
                    continue

                consensus = result.consensus

                # Trim primers before comparison to merge clusters that differ only in primer regions
                # The trimmed consensus is used only for comparison and discarded after merging
                if hasattr(self, 'primers'):
                    consensus, _ = self.trim_primers(consensus)  # Discard found_primers

                consensuses.append(consensus)
                cluster_to_consensus[i] = consensus
                pbar.update(1)

        consensus_to_clusters = defaultdict(list)
        
        if self.disable_homopolymer_equivalence:
            # Only merge exactly identical sequences
            for i, consensus in enumerate(consensuses):
                consensus_to_clusters[consensus].append(i)
        else:
            # Group by homopolymer-equivalent sequences
            for i, consensus in enumerate(consensuses):
                # Find if this consensus is homopolymer-equivalent to any existing group
                found_group = False
                for existing_consensus in consensus_to_clusters.keys():
                    if self.are_homopolymer_equivalent(consensus, existing_consensus):
                        consensus_to_clusters[existing_consensus].append(i)
                        found_group = True
                        break
                
                if not found_group:
                    consensus_to_clusters[consensus].append(i)

        merged = []
        merged_indices = set()

        # Handle clusters with equivalent consensus sequences
        for equivalent_clusters in consensus_to_clusters.values():
            if len(equivalent_clusters) > 1:
                # Merge clusters with equivalent consensus
                merged_read_ids = set()
                merged_from_list = []

                for idx in equivalent_clusters:
                    merged_read_ids.update(clusters[idx]['read_ids'])
                    merged_indices.add(idx)

                    # Track what we're merging from
                    cluster_info = {
                        'initial_cluster_num': clusters[idx]['initial_cluster_num'],
                        'allele_combo': clusters[idx].get('allele_combo'),
                        'size': len(clusters[idx]['read_ids'])
                    }
                    merged_from_list.append(cluster_info)

                # Create merged cluster with provenance
                merged_cluster = {
                    'read_ids': merged_read_ids,
                    'initial_cluster_num': None,  # Multiple sources
                    'allele_combo': None,  # Multiple alleles merged
                    'variant_positions': None,
                    'num_variants': 0,
                    'merged_from': merged_from_list  # Track merge provenance
                }
                merged.append(merged_cluster)

        # Add remaining unmerged clusters
        for i, cluster_dict in enumerate(clusters):
            if i not in merged_indices:
                merged.append(cluster_dict)

        merge_type = "identical" if self.disable_homopolymer_equivalence else "homopolymer-equivalent"
        if len(merged) < len(clusters):
            logging.info(f"Merged {len(clusters) - len(merged)} clusters with {merge_type} consensus sequences")
        else:
            logging.info(f"No clusters were merged (no {merge_type} consensus sequences found)")

        return merged

    def _find_root(self, merged_to: List[int], i: int) -> int:
        """Find the root index of a merged cluster using path compression."""
        if merged_to[i] != i:
            merged_to[i] = self._find_root(merged_to, merged_to[i])
        return merged_to[i]

    def write_cluster_files(self, cluster_num: int, cluster: Set[str],
                            consensus: str, trimmed_consensus: Optional[str] = None,
                            found_primers: Optional[List[str]] = None,
                            rid: Optional[float] = None,
                            rid_min: Optional[float] = None,
                            variant_positions: Optional[List[Dict]] = None,
                            actual_size: Optional[int] = None,
                            consensus_fasta_handle = None,
                            sampled_ids: Optional[Set[str]] = None,
                            msa: Optional[str] = None) -> None:
        """Write cluster files: reads FASTQ, MSA, consensus FASTA, and variant debug files.

        Read identity metrics measure internal cluster consistency (not accuracy vs. ground truth):
        - rid: Mean read identity - measures average agreement between reads and consensus
        - rid_min: Minimum read identity - captures worst-case outlier reads

        High identity values indicate homogeneous clusters with consistent reads.
        Low values may indicate heterogeneity, outliers, or poor consensus (especially at low RiC).

        Variant positions indicate potential unphased biological variation that may
        require cluster subdivision in future processing.
        """
        cluster_size = len(cluster)
        ric_size = min(actual_size or cluster_size, self.max_sample_size)

        # Create info string with size first
        info_parts = [f"size={cluster_size}", f"ric={ric_size}"]

        # Add read identity metrics (as percentages for readability)
        if rid is not None:
            info_parts.append(f"rid={rid*100:.1f}")
        if rid_min is not None:
            info_parts.append(f"rid_min={rid_min*100:.1f}")

        # Add variant flags
        if variant_positions is not None:
            num_variants = len(variant_positions)
            if num_variants > 0:
                info_parts.append(f"has_variants=true")
                info_parts.append(f"num_variants={num_variants}")
            else:
                info_parts.append(f"has_variants=false")

        if found_primers:
            info_parts.append(f"primers={','.join(found_primers)}")
        info_str = " ".join(info_parts)

        # Write reads FASTQ to debug directory with new naming convention
        reads_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-reads.fastq")
        with open(reads_file, 'w') as f:
            for seq_id in cluster:
                SeqIO.write(self.records[seq_id], f, "fastq")

        # Write sampled reads FASTQ (only sequences used for consensus generation)
        if sampled_ids is not None:
            sampled_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-sampled.fastq")
            with open(sampled_file, 'w') as f:
                for seq_id in sampled_ids:
                    SeqIO.write(self.records[seq_id], f, "fastq")

        # Write MSA (multiple sequence alignment) to debug directory
        if msa is not None:
            msa_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-msa.fasta")
            with open(msa_file, 'w') as f:
                f.write(msa)

        # Write variant positions to debug directory
        if variant_positions and len(variant_positions) > 0:
            variant_file = os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-variants.txt")
            with open(variant_file, 'w') as f:
                f.write(f"Variant positions for {self.sample_name}-c{cluster_num}\n")
                f.write(f"Total variants: {len(variant_positions)}\n\n")

                for var in variant_positions:
                    # Format position info
                    if var['consensus_position'] is not None:
                        pos_str = f"MSA:{var['msa_position']}/Consensus:{var['consensus_position']}"
                    else:
                        pos_str = f"MSA:{var['msa_position']}/Consensus:insertion"

                    # Format base composition
                    sorted_bases = sorted(var['base_composition'].items(), key=lambda x: x[1], reverse=True)
                    comp_str = ', '.join([f"{base}:{count}" for base, count in sorted_bases[:4]])

                    # Write variant info
                    f.write(f"Position {pos_str}\n")
                    f.write(f"  Coverage: {var['coverage']}\n")
                    f.write(f"  Error rate: {var['error_rate']*100:.1f}%\n")
                    f.write(f"  Base composition: {comp_str}\n")
                    f.write(f"  Variant bases: {var['variant_bases']}\n")
                    f.write(f"  Reason: {var['reason']}\n\n")

        # Write untrimmed consensus to debug directory
        with open(os.path.join(self.debug_dir, f"{self.sample_name}-c{cluster_num}-RiC{ric_size}-untrimmed.fasta"),
                  'w') as f:
            f.write(f">{self.sample_name}-c{cluster_num} {info_str}\n")
            f.write(consensus + "\n")

        # Write consensus to main output file if handle is provided
        if consensus_fasta_handle:
            final_consensus = trimmed_consensus if trimmed_consensus else consensus
            consensus_fasta_handle.write(f">{self.sample_name}-c{cluster_num} {info_str}\n")
            consensus_fasta_handle.write(final_consensus + "\n")

    def run_mcl_clustering(self, temp_dir: str) -> List[Set[str]]:
        """Run MCL clustering algorithm and return the clusters.

        Args:
            temp_dir: Path to temporary directory for intermediate files

        Returns:
            List of clusters, where each cluster is a set of sequence IDs
        """
        mcl_input = os.path.join(temp_dir, "input.abc")
        mcl_output = os.path.join(temp_dir, "output.cls")

        self.write_mcl_input(mcl_input)

        logging.info(f"Running MCL algorithm with inflation {self.inflation}...")
        self.run_mcl(mcl_input, mcl_output)
        return self.parse_mcl_output(mcl_output)

    def run_greedy_clustering(self, temp_dir: str) -> List[Set[str]]:
        """Run greedy clustering algorithm and return the clusters.

        This algorithm iteratively finds the sequence with the most connections above
        the similarity threshold and forms a cluster around it.

        Args:
            temp_dir: Path to temporary directory for intermediate files

        Returns:
            List of clusters, where each cluster is a set of sequence IDs
        """
        logging.info("Running greedy clustering algorithm...")

        # Build similarity matrix if not already built
        if not hasattr(self, 'alignments'):
            self.alignments = defaultdict(dict)
            self.build_similarity_matrix()

        # Initial clustering
        clusters = []
        available_ids = set(self.sequences.keys())

        cluster_count = 0

        while available_ids:
            center, members = self.find_cluster_center(available_ids)
            available_ids -= members

            clusters.append(members)
            cluster_count += 1

        return clusters

    def build_similarity_matrix(self) -> None:
        """Calculate all pairwise similarities between sequences."""
        logging.info("Calculating pairwise sequence similarities...")

        seq_ids = list(self.sequences.keys())
        total = len(seq_ids) * (len(seq_ids) - 1) // 2

        with tqdm(total=total, desc="Building similarity matrix") as pbar:
            for i, id1 in enumerate(seq_ids):
                for id2 in seq_ids[i + 1:]:
                    sim = self.calculate_similarity(
                        self.sequences[id1],
                        self.sequences[id2]
                    )

                    if sim >= self.min_identity:
                        self.alignments[id1][id2] = sim
                        self.alignments[id2][id1] = sim

                    pbar.update(1)

    def find_cluster_center(self, available_ids: Set[str]) -> Tuple[str, Set[str]]:
        """
        Find the sequence with most similar sequences above threshold,
        and return its ID and the IDs of its cluster members.
        """
        best_center = None
        best_members = set()
        best_count = -1

        for seq_id in available_ids:
            # Get all sequences that align with current sequence
            members = {other_id for other_id in self.alignments.get(seq_id, {})
                       if other_id in available_ids}

            if len(members) > best_count:
                best_count = len(members)
                best_center = seq_id
                best_members = members

        if best_center is None:
            # No alignments found, create singleton cluster with first available sequence
            singleton_id = next(iter(available_ids))
            return singleton_id, {singleton_id}

        best_members.add(best_center)  # Include center in cluster
        return best_center, best_members

    def cluster(self, algorithm: str = "graph") -> None:
        """Perform complete clustering process with variant phasing and write output files.

        Pipeline: Splitting (clustering + outlier removal + phasing) → Merging → Filtering → Output

        Args:
            algorithm: Clustering algorithm to use ('graph' for MCL or 'greedy')
        """
        # Create temporary directory outside the clustering method to allow
        # flexibility in choosing different clustering algorithms
        with tempfile.TemporaryDirectory() as temp_dir:
            # ============================================================
            # PHASE 1: INITIAL CLUSTERING
            # ============================================================
            if algorithm == "graph":
                try:
                    initial_clusters = self.run_mcl_clustering(temp_dir)
                except (subprocess.SubprocessError, FileNotFoundError) as e:
                    logging.error(f"MCL clustering failed: {str(e)}")
                    logging.error("You may need to install MCL: https://micans.org/mcl/")
                    logging.error("Falling back to greedy clustering algorithm...")
                    initial_clusters = self.run_greedy_clustering(temp_dir)
            elif algorithm == "greedy":
                initial_clusters = self.run_greedy_clustering(temp_dir)
            else:
                raise ValueError(f"Unknown clustering algorithm: {algorithm}")

            # Sort initial clusters by size (largest first)
            initial_clusters.sort(key=lambda c: len(c), reverse=True)

            logging.info(f"Initial clustering produced {len(initial_clusters)} clusters")

            # ============================================================
            # PHASE 2: DETECTION + PHASING (Splitting Phase)
            # ============================================================
            # Process each initial cluster to detect variants and phase reads
            all_subclusters = []  # Will hold all phased sub-clusters with provenance

            logging.info("Processing clusters for variant detection and phasing...")

            for initial_idx, cluster in enumerate(initial_clusters, 1):
                cluster_size = len(cluster)

                # Sample sequences for consensus generation if needed
                if cluster_size > self.max_sample_size:
                    logging.debug(f"Initial cluster {initial_idx}: Sampling {self.max_sample_size} from {cluster_size} reads")
                    qualities = []
                    for seq_id in cluster:
                        record = self.records[seq_id]
                        mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                        qualities.append((mean_quality, seq_id))
                    sampled_ids = {seq_id for _, seq_id in
                                   sorted(qualities, reverse=True)[:self.max_sample_size]}
                else:
                    sampled_ids = cluster

                # Generate consensus and MSA
                sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sampled_ids}
                result = self.run_spoa(sampled_seqs)

                if result is None:
                    logging.warning(f"Initial cluster {initial_idx}: Failed to generate consensus, skipping")
                    continue

                consensus = result.consensus
                msa = result.msa_string
                alignments = result.alignments
                msa_to_consensus_pos = result.msa_to_consensus_pos

                # Calculate per-read identity from MSA
                rid, rid_min = self.calculate_read_identity(alignments, consensus)

                # Optional: Remove outlier reads and regenerate consensus
                if self.outlier_identity_threshold is not None:
                    keep_ids, outlier_ids = self.identify_outlier_reads(
                        alignments, consensus, sampled_ids, self.outlier_identity_threshold
                    )

                    if outlier_ids:
                        logging.info(f"Initial cluster {initial_idx}: Removing {len(outlier_ids)}/{len(sampled_ids)} outlier reads, "
                                   f"regenerating consensus")

                        # Update sampled_ids to exclude outliers
                        sampled_ids = keep_ids

                        # CRITICAL: Also remove outliers from the full cluster
                        # This ensures outliers don't reappear in phased haplotypes
                        cluster = cluster - outlier_ids

                        # Regenerate consensus with filtered reads
                        sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sampled_ids}
                        result = self.run_spoa(sampled_seqs)

                        # Recalculate identity metrics
                        if result is not None:
                            consensus = result.consensus
                            msa = result.msa_string
                            alignments = result.alignments
                            msa_to_consensus_pos = result.msa_to_consensus_pos
                            rid, rid_min = self.calculate_read_identity(alignments, consensus)

                # Detect variant positions (if second-pass phasing enabled)
                variant_positions = []
                if consensus and alignments and self.enable_secondpass_phasing:
                    variant_positions = self.detect_variant_positions(
                        alignments, consensus, msa_to_consensus_pos
                    )

                    if variant_positions:
                        logging.info(f"Initial cluster {initial_idx}: Detected {len(variant_positions)} variant positions")

                # Phase reads into haplotypes
                # Note: Uses the FULL cluster (not just sampled_ids) for phasing
                phased_haplotypes = self.phase_reads_by_variants(
                    msa, consensus, cluster, variant_positions, alignments
                )

                # Store each haplotype as a sub-cluster with provenance
                for haplotype_idx, (allele_combo, haplotype_reads) in enumerate(phased_haplotypes):
                    subcluster = {
                        'read_ids': haplotype_reads,
                        'initial_cluster_num': initial_idx,
                        'allele_combo': allele_combo,
                        'variant_positions': variant_positions if allele_combo else None,
                        'num_variants': len(variant_positions) if variant_positions else 0
                    }
                    all_subclusters.append(subcluster)

            logging.info(f"After phasing, created {len(all_subclusters)} sub-clusters from {len(initial_clusters)} initial clusters")

            # ============================================================
            # PHASE 3: MERGING
            # ============================================================
            # Merge sub-clusters with identical/homopolymer-equivalent consensus sequences
            logging.info("Merging homopolymer-equivalent sub-clusters...")
            merged_subclusters = self.merge_similar_clusters(all_subclusters)

            logging.info(f"After merging, have {len(merged_subclusters)} clusters")

            # ============================================================
            # PHASE 4: FILTERING
            # ============================================================
            # Filter by absolute size
            large_clusters = [c for c in merged_subclusters if len(c['read_ids']) >= self.min_size]

            if len(large_clusters) < len(merged_subclusters):
                filtered_count = len(merged_subclusters) - len(large_clusters)
                logging.info(f"Filtered {filtered_count} clusters below minimum size ({self.min_size})")

            # Filter by relative size ratio
            if large_clusters and self.min_cluster_ratio > 0:
                largest_size = max(len(c['read_ids']) for c in large_clusters)
                before_ratio_filter = len(large_clusters)
                large_clusters = [c for c in large_clusters
                                 if len(c['read_ids']) / largest_size >= self.min_cluster_ratio]

                if len(large_clusters) < before_ratio_filter:
                    filtered_count = before_ratio_filter - len(large_clusters)
                    logging.info(f"Filtered {filtered_count} clusters below minimum ratio ({self.min_cluster_ratio})")

            # Sort by size and renumber as c1, c2, c3...
            large_clusters.sort(key=lambda c: len(c['read_ids']), reverse=True)

            total_sequences = len(self.sequences)
            sequences_covered = sum(len(c['read_ids']) for c in large_clusters)

            logging.info(f"Final: {len(large_clusters)} clusters covering {sequences_covered} sequences "
                        f"({sequences_covered / total_sequences:.1%} of total)")

            # ============================================================
            # PHASE 5: OUTPUT
            # ============================================================
            # Generate final consensus and write output files
            consensus_output_file = os.path.join(self.output_dir, f"{self.sample_name}-all.fasta")
            with open(consensus_output_file, 'w') as consensus_fasta_handle:
                for final_idx, cluster_dict in enumerate(large_clusters, 1):
                    cluster = cluster_dict['read_ids']
                    actual_size = len(cluster)

                    # Sample sequences for final consensus generation if needed
                    if len(cluster) > self.max_sample_size:
                        logging.info(f"Cluster {final_idx}: Sampling {self.max_sample_size} from {len(cluster)} reads for final consensus")
                        qualities = []
                        for seq_id in cluster:
                            record = self.records[seq_id]
                            mean_quality = statistics.mean(record.letter_annotations["phred_quality"])
                            qualities.append((mean_quality, seq_id))
                        sampled_ids = {seq_id for _, seq_id in
                                       sorted(qualities, reverse=True)[:self.max_sample_size]}
                    else:
                        sampled_ids = cluster

                    # Generate final consensus and MSA
                    sampled_seqs = {seq_id: self.sequences[seq_id] for seq_id in sampled_ids}
                    result = self.run_spoa(sampled_seqs)

                    # Calculate final identity metrics
                    rid, rid_min = None, None
                    consensus = None
                    msa = None

                    if result is not None:
                        consensus = result.consensus
                        msa = result.msa_string
                        alignments = result.alignments
                        rid, rid_min = self.calculate_read_identity(alignments, consensus)

                    # Note: No outlier removal or variant detection in output phase
                    # (already done in detection phase)

                    if consensus:
                        # Perform primer trimming
                        trimmed_consensus = None
                        found_primers = None
                        if hasattr(self, 'primers'):
                            trimmed_consensus, found_primers = self.trim_primers(consensus)

                        # Write output files (no variant positions in output phase)
                        self.write_cluster_files(
                            cluster_num=final_idx,
                            cluster=cluster,
                            consensus=consensus,
                            trimmed_consensus=trimmed_consensus,
                            found_primers=found_primers,
                            rid=rid,
                            rid_min=rid_min,
                            variant_positions=None,  # Don't report variants in final output
                            actual_size=actual_size,
                            consensus_fasta_handle=consensus_fasta_handle,
                            sampled_ids=sampled_ids,
                            msa=msa
                        )

            # ============================================================
            # WRITE PHASING STATISTICS
            # ============================================================
            self.write_phasing_stats(
                initial_clusters_count=len(initial_clusters),
                subclusters_count=len(all_subclusters),
                merged_count=len(merged_subclusters),
                final_count=len(large_clusters)
            )

    def _create_id_mapping(self) -> None:
        """Create short numeric IDs for all sequences."""
        for i, seq_id in enumerate(self.sequences.keys()):
            short_id = str(i)
            self.id_map[short_id] = seq_id
            self.rev_id_map[seq_id] = short_id

    def calculate_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity using edlib alignment."""
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0

        max_dist = int((1 - self.min_identity) * max(len(seq1), len(seq2)))
        result = edlib.align(seq1, seq2, task="distance", k=max_dist)

        if result["editDistance"] == -1:
            return 0.0

        return 1.0 - (result["editDistance"] / max(len(seq1), len(seq2)))

    def run_spoa(self, read_sequences: Dict[str, str]) -> Optional[MSAResult]:
        """Run SPOA to generate consensus sequence and MSA.

        Args:
            read_sequences: Dictionary mapping read IDs to sequence strings

        Returns:
            MSAResult containing consensus sequence, raw MSA string, and parsed records,
            or None if SPOA fails or input is empty
        """
        if not read_sequences:
            return None

        try:
            # Create temporary input file with actual read IDs
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
                for read_id, seq in read_sequences.items():
                    f.write(f">{read_id}\n{seq}\n")
                temp_input = f.name

            # Construct SPOA command with parameters
            cmd = [
                "spoa",
                temp_input,
                "-r", "2",  # Result mode 2: MSA + consensus (consensus is last)
                "-l", "1",  # Global alignment mode
                "-m", "5",  # Match score
                "-n", "-4",  # Mismatch penalty
                "-g", "-8",  # Gap opening penalty
                "-e", "-6",  # Gap extension penalty
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # Clean up input file
            os.unlink(temp_input)

            # Parse SPOA output - capture both MSA and consensus
            msa_output = result.stdout

            # Parse MSA into ReadAlignment objects
            # Homopolymer normalization is enabled unless explicitly disabled
            enable_normalization = not self.disable_homopolymer_equivalence
            alignments, consensus, msa_to_consensus_pos = extract_alignments_from_msa(
                msa_output,
                enable_homopolymer_normalization=enable_normalization
            )

            if not consensus:
                raise RuntimeError("SPOA did not generate consensus sequence")

            return MSAResult(
                consensus=consensus,
                msa_string=msa_output,
                alignments=alignments,
                msa_to_consensus_pos=msa_to_consensus_pos
            )

        except subprocess.CalledProcessError as e:
            logging.error(f"SPOA failed with return code {e.returncode}")
            logging.error(f"Command: {' '.join(cmd)}")
            logging.error(f"Stderr: {e.stderr}")
            return None

        except RuntimeError:
            # Re-raise RuntimeError (e.g., consensus generation failure) as hard error
            raise

        except Exception as e:
            logging.error(f"Error running SPOA: {str(e)}")
            return None

    def identify_outlier_reads(self, alignments: List[ReadAlignment], consensus_seq: str,
                              sampled_ids: Set[str], threshold: float) -> Tuple[Set[str], Set[str]]:
        """Identify outlier reads below identity threshold using homopolymer-normalized metrics.

        Uses normalized_edit_distance which excludes homopolymer length differences when
        homopolymer normalization is enabled (default behavior).

        Args:
            alignments: List of ReadAlignment objects from MSA
            consensus_seq: Ungapped consensus sequence
            sampled_ids: Set of read IDs that were sampled
            threshold: Identity threshold (e.g., 0.95 for 95% identity)
                      Reads with identity < threshold are considered outliers

        Returns:
            Tuple of (keep_ids, outlier_ids) where both are sets of read IDs
        """
        if not alignments or not consensus_seq:
            return sampled_ids, set()

        try:
            keep_ids = set()
            outlier_ids = set()

            consensus_length = len(consensus_seq)
            if consensus_length == 0:
                return sampled_ids, set()

            for alignment in alignments:
                # Calculate read identity using normalized edit distance
                # (excludes homopolymer length differences when normalization is enabled)
                error_rate = alignment.normalized_edit_distance / consensus_length
                identity = 1.0 - error_rate

                if identity >= threshold:
                    keep_ids.add(alignment.read_id)
                else:
                    outlier_ids.add(alignment.read_id)

            return keep_ids, outlier_ids

        except Exception as e:
            logging.warning(f"Failed to identify outlier reads: {e}")
            return sampled_ids, set()

    def calculate_read_identity(self, alignments: List[ReadAlignment], consensus_seq: str) -> Tuple[Optional[float], Optional[float]]:
        """Calculate read identity metrics from MSA alignments using homopolymer-normalized metrics.

        Read identity measures how well individual reads agree with the consensus sequence.
        This is an internal consistency metric - it does NOT measure accuracy against ground truth.
        High values indicate homogeneous clustering and/or low read error rates.
        Low values may indicate heterogeneous clusters, outliers, or poor consensus quality (esp. at low RiC).

        Uses normalized_edit_distance which excludes homopolymer length differences when
        homopolymer normalization is enabled (default behavior).

        Args:
            alignments: List of ReadAlignment objects from MSA
            consensus_seq: Ungapped consensus sequence

        Returns:
            Tuple of (mean_rid, min_rid) where:
            - mean_rid: Mean read identity across all reads (0.0-1.0)
            - min_rid: Minimum read identity (worst-case read) (0.0-1.0)
            Returns (None, None) if no valid alignments or consensus
        """
        if not alignments or not consensus_seq:
            return None, None

        try:
            # Calculate read identity (1 - error_rate) for each read using normalized metrics
            consensus_length = len(consensus_seq)
            if consensus_length == 0:
                return None, None

            identities = []
            for alignment in alignments:
                # Use normalized edit distance (excludes homopolymer differences)
                error_rate = alignment.normalized_edit_distance / consensus_length
                identity = 1.0 - error_rate
                identities.append(identity)

            if not identities:
                return None, None

            # Calculate per-read statistics
            mean_rid = np.mean(identities)
            min_rid = np.min(identities)

            logging.debug(f"Calculated identity metrics: rid={mean_rid:.3f}, rid_min={min_rid:.3f}")

            return mean_rid, min_rid

        except Exception as e:
            logging.warning(f"Failed to calculate read identity from MSA: {e}")
            return None, None

    def detect_variant_positions(self, alignments: List[ReadAlignment], consensus_seq: str,
                                msa_to_consensus_pos: Dict[int, Optional[int]]) -> List[Dict]:
        """Detect variant positions in MSA for future phasing.

        Identifies positions with significant alternative alleles that may indicate
        unphased biological variation requiring cluster subdivision.

        Args:
            alignments: List of ReadAlignment objects from MSA
            consensus_seq: Ungapped consensus sequence
            msa_to_consensus_pos: Mapping from MSA position to consensus position

        Returns:
            List of variant position dictionaries with metadata, or empty list if none found
        """
        if not alignments or not consensus_seq:
            return []

        try:
            # Extract consensus_aligned from first alignment (they all have the same length)
            if not alignments:
                logging.debug(f"No alignments found in MSA for variant detection")
                return []

            # Reconstruct consensus_aligned from consensus_seq and msa_to_consensus_pos
            msa_length = len(alignments[0].aligned_sequence)
            consensus_aligned = []
            for msa_pos in range(msa_length):
                cons_pos = msa_to_consensus_pos.get(msa_pos)
                if cons_pos is not None:
                    consensus_aligned.append(consensus_seq[cons_pos])
                else:
                    consensus_aligned.append('-')
            consensus_aligned = ''.join(consensus_aligned)

            # Analyze positional variation (error_threshold parameter unused but kept for compatibility)
            position_stats = analyze_positional_variation(
                alignments,
                consensus_seq,
                consensus_aligned,
                msa_to_consensus_pos,
                overall_error_rate=0.0  # Unused, kept for compatibility
            )

            # Identify variant positions
            variant_positions = []
            for pos_stat in position_stats:
                is_variant, variant_bases, reason = is_variant_position_with_composition(
                    pos_stat,
                    min_variant_frequency=self.min_variant_frequency,
                    min_variant_count=self.min_variant_count
                )

                if is_variant:
                    variant_positions.append({
                        'msa_position': pos_stat.msa_position,
                        'consensus_position': pos_stat.consensus_position,
                        'coverage': pos_stat.coverage,
                        'variant_bases': variant_bases,
                        'base_composition': pos_stat.base_composition,
                        'error_rate': pos_stat.error_rate,
                        'reason': reason
                    })

            return variant_positions

        except Exception as e:
            logging.warning(f"Failed to detect variant positions: {e}")
            return []

    def phase_reads_by_variants(
        self,
        msa_string: str,
        consensus_seq: str,
        cluster_read_ids: Set[str],
        variant_positions: List[Dict],
        alignments: Optional[List[ReadAlignment]] = None
    ) -> List[Tuple[str, Set[str]]]:
        """Phase reads into haplotypes based on variant positions.

        Splits a cluster into sub-clusters where each sub-cluster represents reads
        sharing the same alleles at all variant positions. Every unique combination
        of alleles observed in the data becomes a separate haplotype.

        Args:
            msa_string: MSA in FASTA format from SPOA (used if alignments not provided)
            consensus_seq: Ungapped consensus sequence
            cluster_read_ids: Set of read IDs in this cluster (for validation)
            variant_positions: List of variant position dicts from detect_variant_positions()
            alignments: Optional pre-parsed alignments (avoids reparsing)

        Returns:
            List of (allele_combo_string, read_id_set) tuples
            e.g., [("C-T-A", {id1, id2}), ("T-C-A", {id3, id4})]

            Allele combination format: bases separated by hyphens, in MSA position order
            Gap positions are represented as "-" (the gap character itself)
        """
        if not variant_positions:
            # No variants to phase - return single group with all reads
            return [(None, cluster_read_ids)]

        try:
            # Use pre-parsed alignments if provided, otherwise parse msa_string
            if not alignments:
                if not msa_string:
                    return [(None, cluster_read_ids)]
                # Homopolymer normalization is enabled unless explicitly disabled
                enable_normalization = not self.disable_homopolymer_equivalence
                alignments, _, _ = extract_alignments_from_msa(
                    msa_string,
                    enable_homopolymer_normalization=enable_normalization
                )

            if not alignments:
                logging.warning("No alignments found in MSA for phasing")
                return [(None, cluster_read_ids)]

            # Build read_sequences dict from alignments (using gapped sequences)
            read_sequences = {
                alignment.read_id: alignment.aligned_sequence
                for alignment in alignments
            }

            # Extract MSA positions for variants (sorted by position)
            variant_msa_positions = sorted([v['msa_position'] for v in variant_positions])

            # For each read, extract alleles at variant positions
            # Note: read_id is now the actual read ID (not seq{i})
            read_to_alleles = {}
            for read_id, aligned_seq in read_sequences.items():
                # Extract allele at each variant position
                alleles = []
                for msa_pos in variant_msa_positions:
                    if msa_pos < len(aligned_seq):
                        allele = aligned_seq[msa_pos]
                        alleles.append(allele)
                    else:
                        # Read doesn't cover this position - treat as gap
                        alleles.append('-')

                # Create allele combination string
                allele_combo = '-'.join(alleles)
                read_to_alleles[read_id] = allele_combo

            # Group reads by allele combination
            combo_to_reads = defaultdict(set)
            for read_id, allele_combo in read_to_alleles.items():
                combo_to_reads[allele_combo].add(read_id)

            # Convert to list of tuples
            result = [(combo, reads) for combo, reads in combo_to_reads.items()]

            # Sort by size (largest first) for consistency
            result.sort(key=lambda x: len(x[1]), reverse=True)

            logging.debug(f"Phased {len(cluster_read_ids)} reads into {len(result)} haplotypes "
                         f"based on {len(variant_positions)} variant positions")

            return result

        except Exception as e:
            logging.warning(f"Failed to phase reads by variants: {e}")
            # On error, return all reads as single group
            return [(None, cluster_read_ids)]

    def load_primers(self, primer_file: str) -> None:
        """Load primers from FASTA file with position awareness."""
        # Store primers in separate lists by position
        self.forward_primers = []
        self.reverse_primers = []
        self.forward_primers_rc = []  # RC of forward primers
        self.reverse_primers_rc = []  # RC of reverse primers
        
        # For backward compatibility with trim_primers
        self.primers = []  # Will be populated with all primers for existing code
        
        try:
            primer_count = {'forward': 0, 'reverse': 0, 'unknown': 0}
            
            for record in SeqIO.parse(primer_file, "fasta"):
                sequence = str(record.seq)
                sequence_rc = str(reverse_complement(sequence))
                
                # Parse position from header
                if "position=forward" in record.description:
                    self.forward_primers.append((record.id, sequence))
                    self.forward_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    primer_count['forward'] += 1
                elif "position=reverse" in record.description:
                    self.reverse_primers.append((record.id, sequence))
                    self.reverse_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    primer_count['reverse'] += 1
                else:
                    # For primers without position info, add to both lists
                    logging.warning(f"Primer {record.id} has no position annotation, treating as bidirectional")
                    self.forward_primers.append((record.id, sequence))
                    self.forward_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    self.reverse_primers.append((record.id, sequence))
                    self.reverse_primers_rc.append((f"{record.id}_RC", sequence_rc))
                    primer_count['unknown'] += 1
                
                # Maintain backward compatibility
                self.primers.append((record.id, sequence))
                self.primers.append((f"{record.id}_RC", sequence_rc))

            total_primers = sum(primer_count.values())
            if total_primers == 0:
                logging.warning("No primers were loaded. Primer trimming will be disabled.")
            else:
                logging.info(f"Loaded {total_primers} primers: {primer_count['forward']} forward, "
                           f"{primer_count['reverse']} reverse, {primer_count['unknown']} unknown")
        except Exception as e:
            logging.error(f"Error loading primers: {str(e)}")
            raise

    def orient_sequences(self) -> set:
        """Normalize sequence orientations based on primer matches.
        
        Scoring system:
        - +1 point if a forward primer is found at the expected position
        - +1 point if a reverse primer is found at the expected position
        - Maximum score: 2 (both primers found)
        
        Decision logic:
        - If one orientation scores >0 and the other scores 0: use the non-zero orientation
        - If both score 0 or both score >0: keep original orientation (ambiguous/failed)
        """
        if not hasattr(self, 'forward_primers') or not hasattr(self, 'reverse_primers'):
            logging.warning("No positioned primers loaded, skipping orientation")
            return set()
            
        if len(self.forward_primers) == 0 and len(self.reverse_primers) == 0:
            logging.warning("No positioned primers available, skipping orientation")
            return set()
        
        logging.info("Starting sequence orientation based on primer positions...")
        
        oriented_count = 0
        already_correct = 0
        failed_count = 0
        failed_sequences = set()  # Track which sequences failed orientation
        
        # Process each sequence
        for seq_id in tqdm(self.sequences, desc="Orienting sequences"):
            sequence = self.sequences[seq_id]
            
            # Test both orientations (scores will be 0, 1, or 2)
            forward_score = self._score_orientation(sequence, "forward")
            reverse_score = self._score_orientation(sequence, "reverse")
            
            # Decision logic
            if forward_score > 0 and reverse_score == 0:
                # Clear forward orientation
                already_correct += 1
                logging.debug(f"Kept {seq_id} as-is: forward_score={forward_score}, reverse_score={reverse_score}")
            elif reverse_score > 0 and forward_score == 0:
                # Clear reverse orientation - needs to be flipped
                self.sequences[seq_id] = str(reverse_complement(sequence))
                
                # Also update the record if it exists
                if seq_id in self.records:
                    record = self.records[seq_id]
                    record.seq = reverse_complement(record.seq)
                    # Reverse quality scores too if they exist
                    if 'phred_quality' in record.letter_annotations:
                        record.letter_annotations['phred_quality'] = \
                            record.letter_annotations['phred_quality'][::-1]
                
                oriented_count += 1
                logging.debug(f"Reoriented {seq_id}: forward_score={forward_score}, reverse_score={reverse_score}")
            else:
                # Both zero (no primers) or both non-zero (ambiguous) - orientation failed
                failed_count += 1
                failed_sequences.add(seq_id)  # Track this sequence as failed
                if forward_score == 0 and reverse_score == 0:
                    logging.debug(f"No primer matches for {seq_id}")
                else:
                    logging.debug(f"Ambiguous orientation for {seq_id}: forward_score={forward_score}, reverse_score={reverse_score}")
        
        logging.info(f"Orientation complete: {already_correct} kept as-is, "
                    f"{oriented_count} reverse-complemented, {failed_count} orientation failed")
        
        # Return set of failed sequence IDs for potential filtering
        return failed_sequences
    
    def _score_orientation(self, sequence: str, orientation: str) -> int:
        """Score how well primers match in the given orientation.
        
        Simple binary scoring:
        - +1 if a forward primer is found at the expected position
        - +1 if a reverse primer is found at the expected position
        
        Args:
            sequence: The sequence to test
            orientation: Either "forward" or "reverse"
            
        Returns:
            Score from 0-2 (integer)
        """
        score = 0
        
        if orientation == "forward":
            # Forward orientation:
            # - Check for forward primers at 5' end (as-is)
            # - Check for RC of reverse primers at 3' end
            if self._has_primer_match(sequence, self.forward_primers, "start"):
                score += 1
            if self._has_primer_match(sequence, self.reverse_primers_rc, "end"):
                score += 1
        else:
            # Reverse orientation:
            # - Check for reverse primers at 5' end (as-is)
            # - Check for RC of forward primers at 3' end
            if self._has_primer_match(sequence, self.reverse_primers, "start"):
                score += 1
            if self._has_primer_match(sequence, self.forward_primers_rc, "end"):
                score += 1
        
        return score
    
    def _has_primer_match(self, sequence: str, primers: List[Tuple[str, str]], end: str) -> bool:
        """Check if any primer matches at the specified end of sequence.
        
        Args:
            sequence: The sequence to search in
            primers: List of (name, sequence) tuples to search for
            end: Either "start" or "end"
            
        Returns:
            True if any primer has a good match, False otherwise
        """
        if not primers or not sequence:
            return False
        
        # Determine search region
        max_primer_len = max(len(p[1]) for p in primers) if primers else 50
        if end == "start":
            search_region = sequence[:min(max_primer_len * 2, len(sequence))]
        else:
            search_region = sequence[-min(max_primer_len * 2, len(sequence)):]
        
        for primer_name, primer_seq in primers:
            # Allow up to 25% errors
            k = len(primer_seq) // 4
            
            # Use edlib to find best match
            result = edlib.align(primer_seq, search_region, task="distance", mode="HW", k=k)
            
            if result["editDistance"] != -1:
                # Consider it a match if identity is >= 75%
                identity = 1.0 - (result["editDistance"] / len(primer_seq))
                if identity >= 0.75:
                    logging.debug(f"Found {primer_name} at {end} with identity {identity:.2%} "
                                f"(edit_dist={result['editDistance']}, len={len(primer_seq)})")
                    return True
        
        return False
    
    def trim_primers(self, sequence: str) -> Tuple[str, List[str]]:
        """
        Trim primers from start and end of sequence.
        Returns (trimmed_sequence, [found_primers])
        """
        if not hasattr(self, 'primers') or not self.primers:
            return sequence, []

        found_primers = []
        trimmed_seq = sequence

        logging.debug(f"Starting primer trimming on sequence of length {len(sequence)}")

        # Look for primer at 5' end
        best_start_dist = float('inf')
        best_start_primer = None
        best_start_end = None

        for primer_name, primer_seq in self.primers:
            k = len(primer_seq) // 4  # Allow ~25% errors
            search_region = sequence[:len(primer_seq) * 2]

            result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)
            logging.debug(f"Testing {primer_name} at 5' end (k={k})")

            if result["editDistance"] != -1:
                dist = result["editDistance"]
                if dist < best_start_dist:
                    best_start_dist = dist
                    best_start_primer = primer_name
                    best_start_end = result["locations"][0][1] + 1

        if best_start_primer:
            logging.debug(f"Found {best_start_primer} at 5' end (k={best_start_dist})")
            found_primers.append(f"5'-{best_start_primer}")
            trimmed_seq = trimmed_seq[best_start_end:]

        # Look for primer at 3' end
        best_end_dist = float('inf')
        best_end_primer = None
        best_end_start = None

        for primer_name, primer_seq in self.primers:
            k = len(primer_seq) // 4  # Allow ~25% errors
            search_region = sequence[-len(primer_seq) * 2:]

            result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)

            if result["editDistance"] != -1:
                dist = result["editDistance"]
                if dist < best_end_dist:
                    best_end_dist = dist
                    best_end_primer = primer_name
                    base_pos = len(trimmed_seq) - len(search_region)
                    best_end_start = base_pos + result["locations"][0][0]

        if best_end_primer:
            logging.debug(f"Found {best_end_primer} at 3' end (k={best_end_dist})")
            found_primers.append(f"3'-{best_end_primer}")
            trimmed_seq = trimmed_seq[:best_end_start]

        return trimmed_seq, found_primers

    def calculate_consensus_distance(self, seq1: str, seq2: str, require_merge_compatible: bool = False) -> int:
        """Calculate distance between two consensus sequences using adjusted identity.
        
        Uses custom adjustment parameters that enable only homopolymer normalization:
        - Homopolymer differences (e.g., AAA vs AAAAA) are treated as identical
        - Regular substitutions count as mismatches
        - Non-homopolymer indels optionally prevent merging
        
        Args:
            seq1: First consensus sequence
            seq2: Second consensus sequence
            require_merge_compatible: If True, return -1 when sequences have variations
                                     that cannot be represented in IUPAC consensus (indels)
        
        Returns:
            Distance between sequences (substitutions only), or -1 if require_merge_compatible=True
            and sequences contain non-homopolymer indels
        """
        if not seq1 or not seq2:
            return max(len(seq1), len(seq2))

        # Get alignment from edlib (uses global NW alignment by default)
        result = edlib.align(seq1, seq2, task="path")
        if result["editDistance"] == -1:
            # Alignment failed, return maximum possible distance
            return max(len(seq1), len(seq2))
            
        # Get nice alignment for adjusted identity scoring
        alignment = edlib.getNiceAlignment(result, seq1, seq2)
        if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
            # Fall back to edit distance if alignment extraction fails
            return result["editDistance"]
            
        # Configure custom adjustment parameters for homopolymer normalization only
        custom_params = AdjustmentParams(
            normalize_homopolymers=True,    # Enable homopolymer normalization
            handle_iupac_overlap=False,     # Disable IUPAC overlap handling
            normalize_indels=False,         # Disable indel normalization
            end_skip_distance=0,            # Disable end trimming
            max_repeat_motif_length=2       # Keep default motif length
        )
        
        # Create custom scoring format to distinguish indels from substitutions
        from adjusted_identity import ScoringFormat
        custom_format = ScoringFormat(
            match='|',
            substitution='X',     # Distinct code for substitutions
            indel_start='I',      # Distinct code for indels
            indel_extension='-',
            homopolymer_extension='=',
            end_trimmed='.'
        )
        
        # Calculate adjusted identity with custom format
        score_result = score_alignment(
            alignment['query_aligned'], 
            alignment['target_aligned'],
            adjustment_params=custom_params,
            scoring_format=custom_format
        )
        
        # Check for merge compatibility if requested
        # Both non-homopolymer indels ('I') and terminal overhangs ('.') prevent merging
        if require_merge_compatible:
            if 'I' in score_result.score_aligned:
                logging.debug(f"Non-homopolymer indel detected, sequences not merge-compatible")
                return -1  # Signal that merging should not occur
            if '.' in score_result.score_aligned:
                logging.debug(f"Terminal overhang detected, sequences not merge-compatible")
                return -1  # Signal that merging should not occur
        
        # Count only substitutions (not homopolymer adjustments or indels)
        # Note: mismatches includes both substitutions and non-homopolymer indels
        # For accurate distance when indels are present, we use the mismatches count
        distance = score_result.mismatches
        
        # Log details about the variations found
        substitutions = score_result.score_aligned.count('X')
        indels = score_result.score_aligned.count('I')
        homopolymers = score_result.score_aligned.count('=')
        
        logging.debug(f"Consensus distance: {distance} total mismatches "
                     f"({substitutions} substitutions, {indels} indels, "
                     f"{homopolymers} homopolymer adjustments)")
        
        return distance

    def are_homopolymer_equivalent(self, seq1: str, seq2: str) -> bool:
        """Check if two sequences are equivalent when considering only homopolymer differences.

        Uses adjusted-identity scoring with global alignment. Terminal overhangs (marked as '.')
        and non-homopolymer indels (marked as 'I') prevent merging, ensuring truncated sequences
        don't merge with full-length sequences.
        """
        if not seq1 or not seq2:
            return seq1 == seq2

        # Use calculate_consensus_distance with merge compatibility check
        # Global alignment ensures terminal gaps are counted as indels
        # Returns: -1 (non-homopolymer indels), 0 (homopolymer-equivalent), >0 (substitutions)
        # Only distance == 0 means truly homopolymer-equivalent
        distance = self.calculate_consensus_distance(seq1, seq2, require_merge_compatible=True)
        return distance == 0

    def parse_mcl_output(self, mcl_output_file: str) -> List[Set[str]]:
        """Parse MCL output file into clusters of original sequence IDs."""
        clusters = []
        with open(mcl_output_file) as f:
            for line in f:
                # Each line is a tab-separated list of cluster members
                short_ids = line.strip().split('\t')
                # Map short IDs back to original sequence IDs
                cluster = {self.id_map[short_id] for short_id in short_ids}
                clusters.append(cluster)
        return clusters



def main():
    parser = argparse.ArgumentParser(
        description="MCL-based clustering of nanopore amplicon reads"
    )
    parser.add_argument("input_file", help="Input FASTQ file")
    parser.add_argument("--augment-input", help="Additional FASTQ/FASTA file with sequences recovered after primary demultiplexing (e.g., from specimine)")
    parser.add_argument("--algorithm", type=str, default="graph", choices=["graph", "greedy"],
                        help="Clustering algorithm to use (default: graph)")
    parser.add_argument("--min-identity", type=float, default=0.9,
                        help="Minimum sequence identity threshold for clustering (default: 0.9)")
    parser.add_argument("--inflation", type=float, default=4.0,
                        help="MCL inflation parameter (default: 4.0)")
    parser.add_argument("--min-size", type=int, default=5,
                        help="Minimum cluster size (default: 5, 0 to disable)")
    parser.add_argument("--min-cluster-ratio", type=float, default=0.2,
                        help="Minimum size ratio between a cluster and the largest cluster (default: 0.2, 0 to disable)")
    parser.add_argument("--max-sample-size", type=int, default=500,
                        help="Maximum cluster size for consensus (default: 500)")
    parser.add_argument("--outlier-identity-threshold", type=float, default=None,
                        help="Identity threshold for outlier removal (e.g., 0.95 for 95%% identity). "
                             "Reads with identity below this threshold will be removed and consensus regenerated. "
                             "Auto-calculated as (1 + min_identity) / 2 if not specified. "
                             "For min_identity=0.9, default is 0.95.")
    parser.add_argument("--disable-secondpass-phasing", action="store_true",
                        help="Disable second-pass variant phasing (enabled by default). "
                             "First-pass MCL clustering already handles most variants.")
    parser.add_argument("--min-variant-frequency", type=float, default=0.20,
                        help="Minimum alternative allele frequency to call variant (default: 0.20 for 20%%)")
    parser.add_argument("--min-variant-count", type=int, default=5,
                        help="Minimum alternative allele read count to call variant (default: 5)")
    parser.add_argument("--presample", type=int, default=1000,
                        help="Presample size for initial reads (default: 1000, 0 to disable)")
    parser.add_argument("--k-nearest-neighbors", type=int, default=5,
                        help="Number of nearest neighbors for graph construction (default: 5)")
    parser.add_argument("--primers", help="FASTA file containing primer sequences (default: looks for primers.fasta in input file directory)")
    parser.add_argument("-O", "--output-dir", default="clusters",
                        help="Output directory for all files (default: clusters)")
    parser.add_argument("--disable-homopolymer-equivalence", action="store_true",
                        help="Disable homopolymer equivalence in cluster merging (only merge identical sequences)")
    parser.add_argument("--orient-mode", choices=["skip", "keep-all", "filter-failed"], default="skip",
                        help="Sequence orientation mode: skip (default, no orientation), keep-all (orient but keep failed), or filter-failed (orient and remove failed)")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    parser.add_argument("--version", action="version",
                        version=f"Speconsense {__version__}",
                        help="Show program's version number and exit")

    args = parser.parse_args()

    # Setup standard logging
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format=log_format
    )

    sample = os.path.splitext(os.path.basename(args.input_file))[0]
    clusterer = SpecimenClusterer(
        min_identity=args.min_identity,
        inflation=args.inflation,
        min_size=args.min_size,
        min_cluster_ratio=args.min_cluster_ratio,
        max_sample_size=args.max_sample_size,
        presample_size=args.presample,
        k_nearest_neighbors=args.k_nearest_neighbors,
        sample_name=sample,
        disable_homopolymer_equivalence=args.disable_homopolymer_equivalence,
        output_dir=args.output_dir,
        outlier_identity_threshold=args.outlier_identity_threshold,
        enable_secondpass_phasing=not args.disable_secondpass_phasing,
        min_variant_frequency=args.min_variant_frequency,
        min_variant_count=args.min_variant_count
    )

    # Log configuration
    if args.outlier_identity_threshold is not None:
        logging.info(f"Outlier removal enabled: min_identity_threshold={args.outlier_identity_threshold*100:.1f}% (user-specified)")
    else:
        # Auto-calculated threshold
        auto_threshold = (1.0 + args.min_identity) / 2.0
        logging.info(f"Outlier removal enabled: min_identity_threshold={auto_threshold*100:.1f}% (auto-calculated from min_identity={args.min_identity*100:.1f}%)")

    if not args.disable_secondpass_phasing:
        logging.info(f"Second-pass variant phasing enabled: min_freq={args.min_variant_frequency:.0%}, "
                    f"min_count={args.min_variant_count}")

    # Set additional attributes for metadata
    clusterer.input_file = os.path.abspath(args.input_file)
    clusterer.augment_input = os.path.abspath(args.augment_input) if args.augment_input else None
    clusterer.algorithm = args.algorithm
    clusterer.orient_mode = args.orient_mode

    # Read primary sequences
    logging.info(f"Reading sequences from {args.input_file}")
    format = "fasta" if args.input_file.endswith(".fasta") else "fastq"
    records = list(SeqIO.parse(args.input_file, format))
    logging.info(f"Loaded {len(records)} primary sequences")

    # Load augmented sequences if specified
    augment_records = None
    if args.augment_input:
        # Check if augment input file exists
        if not os.path.exists(args.augment_input):
            logging.error(f"Augment input file not found: {args.augment_input}")
            sys.exit(1)
            
        logging.info(f"Reading augmented sequences from {args.augment_input}")
        
        # Auto-detect format like main input
        augment_format = "fasta" if args.augment_input.endswith(".fasta") else "fastq"
        
        try:
            augment_records = list(SeqIO.parse(args.augment_input, augment_format))
            logging.info(f"Loaded {len(augment_records)} augmented sequences")
            
            if len(augment_records) == 0:
                logging.warning(f"No sequences found in augment input file: {args.augment_input}")
            
            # Add dummy quality scores to FASTA sequences so they can be written as FASTQ later
            if augment_format == "fasta":
                for record in augment_records:
                    if not hasattr(record, 'letter_annotations') or 'phred_quality' not in record.letter_annotations:
                        # Add dummy quality scores (quality 30 = '?' in FASTQ)
                        record.letter_annotations = {'phred_quality': [30] * len(record.seq)}
                logging.debug(f"Added quality scores to {len(augment_records)} FASTA sequences for downstream compatibility")
                
        except Exception as e:
            logging.error(f"Failed to read augment input file '{args.augment_input}': {e}")
            sys.exit(1)

    # Add sequences to clusterer (both primary and augmented)
    clusterer.add_sequences(records, augment_records)

    if args.primers:
        clusterer.primers_file = os.path.abspath(args.primers)
        clusterer.load_primers(args.primers)
    else:
        # Look for primers.fasta in the same directory as the input file
        input_dir = os.path.dirname(os.path.abspath(args.input_file))
        auto_primer_path = os.path.join(input_dir, "primers.fasta")

        if os.path.exists(auto_primer_path):
            logging.info(f"Found primers.fasta in input directory: {auto_primer_path}")
            clusterer.primers_file = os.path.abspath(auto_primer_path)
            clusterer.load_primers(auto_primer_path)
        else:
            logging.warning("No primer file specified and primers.fasta not found in input directory. Primer trimming will be disabled.")
            clusterer.primers_file = None
    
    # Handle sequence orientation based on mode
    if args.orient_mode != "skip":
        if hasattr(clusterer, 'forward_primers') and hasattr(clusterer, 'reverse_primers'):
            failed_sequences = clusterer.orient_sequences()

            # Filter failed sequences if requested
            if args.orient_mode == "filter-failed" and failed_sequences:
                logging.info(f"Filtering out {len(failed_sequences)} sequences with failed orientation")

                # Remove failed sequences from clusterer
                for seq_id in failed_sequences:
                    del clusterer.sequences[seq_id]
                    del clusterer.records[seq_id]

                remaining = len(clusterer.sequences)
                logging.info(f"Continuing with {remaining} successfully oriented sequences")
        else:
            logging.warning(f"--orient-mode={args.orient_mode} specified but no primers with position information loaded")

    # Write metadata file for use by post-processing tools
    clusterer.write_metadata()

    clusterer.cluster(algorithm=args.algorithm)
    print()

if __name__ == "__main__":
    main()

