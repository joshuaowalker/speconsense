"""Module-level worker functions for ProcessPoolExecutor.

These must be at module level to be picklable for multiprocessing.
Includes standalone versions of functions used by workers and config classes.
"""

import logging
import os
import subprocess
import tempfile
from io import StringIO
from typing import Dict, List, Optional, Set, Tuple

import edlib
import numpy as np
from Bio import SeqIO

from speconsense.msa import (
    MSAResult,
    ReadAlignment,
    analyze_positional_variation,
    call_iupac_ambiguities,
    calculate_within_cluster_error,
    extract_alignments_from_msa,
    filter_qualifying_haplotypes,
    group_reads_by_single_position,
    group_reads_by_position_tuple,
    is_variant_position_with_composition,
)
from speconsense.significance import is_variant_significant


# Configuration classes for parallel processing

class ClusterProcessingConfig:
    """Configuration for parallel cluster processing.

    Passed to worker processes to avoid needing to pickle the entire SpecimenClusterer.
    """
    __slots__ = ['outlier_identity_threshold', 'enable_secondpass_phasing',
                 'disable_homopolymer_equivalence', 'min_variant_frequency', 'min_variant_count',
                 'total_specimen_reads', 'assumed_error_rate', 'significance_level',
                 'min_k_position_gap', 'k_correlation_threshold']

    def __init__(self, outlier_identity_threshold: Optional[float],
                 enable_secondpass_phasing: bool,
                 disable_homopolymer_equivalence: bool,
                 min_variant_frequency: float,
                 min_variant_count: int,
                 total_specimen_reads: int = 0,
                 assumed_error_rate: float = 0.015,
                 significance_level: float = 1e-5,
                 min_k_position_gap: int = 10,
                 k_correlation_threshold: float = 0.9):
        self.outlier_identity_threshold = outlier_identity_threshold
        self.enable_secondpass_phasing = enable_secondpass_phasing
        self.disable_homopolymer_equivalence = disable_homopolymer_equivalence
        self.min_variant_frequency = min_variant_frequency
        self.min_variant_count = min_variant_count
        self.total_specimen_reads = total_specimen_reads
        self.assumed_error_rate = assumed_error_rate
        self.significance_level = significance_level
        self.min_k_position_gap = min_k_position_gap
        self.k_correlation_threshold = k_correlation_threshold


class ConsensusGenerationConfig:
    """Configuration for parallel final consensus generation.

    Passed to worker processes to avoid needing to pickle the entire SpecimenClusterer.
    """
    __slots__ = ['max_sample_size', 'enable_iupac_calling', 'min_ambiguity_frequency',
                 'min_ambiguity_count', 'disable_homopolymer_equivalence', 'primers']

    def __init__(self, max_sample_size: int,
                 enable_iupac_calling: bool,
                 min_ambiguity_frequency: float,
                 min_ambiguity_count: int,
                 disable_homopolymer_equivalence: bool,
                 primers: Optional[List[Tuple[str, str]]] = None):
        self.max_sample_size = max_sample_size
        self.enable_iupac_calling = enable_iupac_calling
        self.min_ambiguity_frequency = min_ambiguity_frequency
        self.min_ambiguity_count = min_ambiguity_count
        self.disable_homopolymer_equivalence = disable_homopolymer_equivalence
        self.primers = primers


# Worker functions

def _run_spoa_worker(args: Tuple[int, Dict[str, str], bool]) -> Tuple[int, Optional[MSAResult]]:
    """Worker function for parallel SPOA execution.

    Must be at module level for ProcessPoolExecutor pickling.

    Args:
        args: Tuple of (cluster_idx, sampled_seqs, disable_homopolymer_equivalence)

    Returns:
        Tuple of (cluster_idx, MSAResult or None)
    """
    cluster_idx, sampled_seqs, disable_homopolymer_equivalence = args

    if not sampled_seqs:
        return cluster_idx, None

    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
            for read_id, seq in sorted(sampled_seqs.items()):
                f.write(f">{read_id}\n{seq}\n")
            temp_input = f.name

        cmd = [
            "spoa", temp_input,
            "-r", "2", "-l", "1", "-m", "5", "-n", "-4", "-g", "-8", "-e", "-6",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        os.unlink(temp_input)

        enable_normalization = not disable_homopolymer_equivalence
        alignments, consensus, msa_to_consensus_pos = extract_alignments_from_msa(
            result.stdout, enable_homopolymer_normalization=enable_normalization
        )

        if not consensus:
            return cluster_idx, None

        return cluster_idx, MSAResult(
            consensus=consensus,
            msa_string=result.stdout,
            alignments=alignments,
            msa_to_consensus_pos=msa_to_consensus_pos
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"SPOA worker failed for cluster {cluster_idx}: return code {e.returncode}")
        return cluster_idx, None
    except Exception as e:
        logging.error(f"SPOA worker failed for cluster {cluster_idx}: {e}")
        return cluster_idx, None


def _run_spoa_for_cluster_worker(sequences: Dict[str, str],
                                  disable_homopolymer_equivalence: bool) -> Optional[MSAResult]:
    """Run SPOA for a set of sequences. Used by cluster processing worker."""
    if not sequences:
        return None

    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
            for read_id, seq in sorted(sequences.items()):
                f.write(f">{read_id}\n{seq}\n")
            temp_input = f.name

        cmd = [
            "spoa", temp_input,
            "-r", "2", "-l", "1", "-m", "5", "-n", "-4", "-g", "-8", "-e", "-6",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        os.unlink(temp_input)

        enable_normalization = not disable_homopolymer_equivalence
        alignments, consensus, msa_to_consensus_pos = extract_alignments_from_msa(
            result.stdout, enable_homopolymer_normalization=enable_normalization
        )

        if not consensus:
            return None

        return MSAResult(
            consensus=consensus,
            msa_string=result.stdout,
            alignments=alignments,
            msa_to_consensus_pos=msa_to_consensus_pos
        )
    except Exception:
        return None


def _identify_outlier_reads_standalone(alignments: List[ReadAlignment], consensus_seq: str,
                                        sampled_ids: Set[str], threshold: float) -> Tuple[Set[str], Set[str]]:
    """Identify outlier reads below identity threshold. Standalone version for workers."""
    if not alignments or not consensus_seq:
        return sampled_ids, set()

    try:
        keep_ids = set()
        outlier_ids = set()
        consensus_length = len(consensus_seq)
        if consensus_length == 0:
            return sampled_ids, set()

        for alignment in alignments:
            error_rate = alignment.normalized_edit_distance / consensus_length
            identity = 1.0 - error_rate
            if identity >= threshold:
                keep_ids.add(alignment.read_id)
            else:
                outlier_ids.add(alignment.read_id)

        return keep_ids, outlier_ids
    except Exception:
        return sampled_ids, set()


def _calculate_read_identity_standalone(alignments: List[ReadAlignment],
                                         consensus_seq: str) -> Tuple[Optional[float], Optional[float]]:
    """Calculate read identity metrics. Standalone version for workers."""
    if not alignments or not consensus_seq:
        return None, None

    try:
        consensus_length = len(consensus_seq)
        if consensus_length == 0:
            return None, None

        identities = []
        for alignment in alignments:
            error_rate = alignment.normalized_edit_distance / consensus_length
            identity = 1.0 - error_rate
            identities.append(identity)

        if not identities:
            return None, None

        return np.mean(identities), np.min(identities)
    except Exception:
        return None, None


def _detect_variant_positions_standalone(alignments: List[ReadAlignment], consensus_seq: str,
                                          msa_to_consensus_pos: Dict[int, Optional[int]],
                                          min_variant_frequency: float,
                                          min_variant_count: int) -> List[Dict]:
    """Detect variant positions in MSA. Standalone version for workers."""
    if not alignments or not consensus_seq:
        return []

    try:
        msa_length = len(alignments[0].aligned_sequence)
        consensus_aligned = []
        for msa_pos in range(msa_length):
            cons_pos = msa_to_consensus_pos.get(msa_pos)
            if cons_pos is not None:
                consensus_aligned.append(consensus_seq[cons_pos])
            else:
                consensus_aligned.append('-')
        consensus_aligned = ''.join(consensus_aligned)

        position_stats = analyze_positional_variation(alignments, consensus_aligned, msa_to_consensus_pos)

        variant_positions = []
        for pos_stat in position_stats:
            is_variant, variant_bases, reason = is_variant_position_with_composition(
                pos_stat,
                min_variant_frequency=min_variant_frequency,
                min_variant_count=min_variant_count
            )
            if is_variant:
                variant_positions.append({
                    'msa_position': pos_stat.msa_position,
                    'consensus_position': pos_stat.consensus_position,
                    'coverage': pos_stat.coverage,
                    'variant_bases': variant_bases,
                    'base_composition': pos_stat.base_composition,
                    'homopolymer_composition': pos_stat.homopolymer_composition,
                    'error_rate': pos_stat.error_rate,
                    'reason': reason
                })

        return variant_positions
    except Exception:
        return []


def _find_correlated_positions(
    k1_failed: List[Tuple[int, Set[str]]],
    msa_to_consensus_pos: Dict[int, Optional[int]],
    read_to_position_alleles: Dict[str, Dict[int, str]],
    all_read_ids: Set[str],
    total_reads: int,
    consensus: str,
    config: 'ClusterProcessingConfig'
) -> Optional[Tuple[int, Dict, Dict, float]]:
    """Find correlated variant positions for K>1 CER evaluation.

    For each seed position that failed K=1 CER, find other failed positions
    whose minority allele correlates with the seed. If the joint K-position
    pattern passes CER, return the best candidate split.

    Returns (best_pos, qualifying, non_qualifying, error) or None.
    The best_pos is the seed position; qualifying uses allele tuples as keys.
    """
    from speconsense.msa import calculate_within_cluster_error

    best_combo_error = float('inf')
    best_result = None
    all_positions = set(pos for pos, _ in k1_failed)

    for seed_idx, (seed_pos, seed_minority_reads) in enumerate(k1_failed):
        seed_cons_pos = msa_to_consensus_pos.get(seed_pos)

        correlated_positions = [seed_pos]
        joint_minority_reads = seed_minority_reads.copy()

        for other_idx, (other_pos, other_minority_reads) in enumerate(k1_failed):
            if other_idx == seed_idx:
                continue

            # Proximity filter in consensus coordinates
            other_cons_pos = msa_to_consensus_pos.get(other_pos)
            if seed_cons_pos is not None and other_cons_pos is not None:
                if abs(seed_cons_pos - other_cons_pos) < config.min_k_position_gap:
                    continue

            # Check correlation: what fraction of joint minority reads
            # also have the minority allele at the other position?
            overlap = joint_minority_reads & other_minority_reads
            if len(joint_minority_reads) > 0 and len(overlap) >= len(joint_minority_reads) * config.k_correlation_threshold:
                correlated_positions.append(other_pos)
                joint_minority_reads = overlap

        K = len(correlated_positions)
        if K < 2:
            continue

        # Group reads by allele tuple at all correlated positions
        allele_groups = group_reads_by_position_tuple(
            read_to_position_alleles, correlated_positions, all_read_ids
        )
        qualifying, non_qualifying = filter_qualifying_haplotypes(
            allele_groups, total_reads, config.min_variant_count, config.min_variant_frequency
        )
        if len(qualifying) < 2:
            continue

        min_M = min(len(reads) for reads in qualifying.values())
        is_sig, p_star = is_variant_significant(
            M=min_M, N=config.total_specimen_reads,
            L=len(consensus),
            assumed_error_rate=config.assumed_error_rate,
            alpha=config.significance_level,
            K=K
        )
        if not is_sig:
            continue

        logging.debug(f"K>1 CER pass: K={K} positions {correlated_positions}, "
                      f"M={min_M}, p*={p_star:.4f}")

        error = calculate_within_cluster_error(
            qualifying, read_to_position_alleles, set(correlated_positions), all_positions
        )
        if error < best_combo_error:
            best_combo_error = error
            best_result = (seed_pos, qualifying, non_qualifying, error)

    return best_result


def _recursive_phase_cluster_standalone(
    read_ids: Set[str],
    read_sequences: Dict[str, str],
    path: List[str],
    depth: int,
    config: ClusterProcessingConfig
) -> Tuple[List[Tuple[List[str], str, Set[str]]], Set[str]]:
    """Recursively phase a cluster. Standalone version for workers."""
    total_reads = len(read_ids)

    # Base case: cluster too small
    if total_reads < config.min_variant_count * 2:
        leaf_seqs = {rid: read_sequences[rid] for rid in read_ids}
        result = _run_spoa_for_cluster_worker(leaf_seqs, config.disable_homopolymer_equivalence)
        consensus = result.consensus if result else ""
        return [(path, consensus, read_ids)], set()

    # Generate MSA
    cluster_seqs = {rid: read_sequences[rid] for rid in read_ids}
    result = _run_spoa_for_cluster_worker(cluster_seqs, config.disable_homopolymer_equivalence)

    if result is None:
        return [(path, "", read_ids)], set()

    consensus = result.consensus
    alignments = result.alignments
    msa_to_consensus_pos = result.msa_to_consensus_pos

    # Detect variants
    variant_positions = _detect_variant_positions_standalone(
        alignments, consensus, msa_to_consensus_pos,
        config.min_variant_frequency, config.min_variant_count
    )

    if not variant_positions:
        return [(path, consensus, read_ids)], set()

    # Parse MSA for consensus_aligned
    msa_handle = StringIO(result.msa_string)
    records = list(SeqIO.parse(msa_handle, 'fasta'))

    consensus_aligned = None
    for record in records:
        if 'Consensus' in record.description or 'Consensus' in record.id:
            consensus_aligned = str(record.seq).upper()
            break

    if not consensus_aligned:
        return [(path, consensus, read_ids)], set()

    # Build read to alignment mapping
    read_to_alignment = {a.read_id: a for a in alignments}

    # Extract alleles at variant positions
    variant_msa_positions = sorted([v['msa_position'] for v in variant_positions])
    read_to_position_alleles = {}

    for read_id in read_ids:
        alignment = read_to_alignment.get(read_id)
        if not alignment:
            continue

        aligned_seq = alignment.aligned_sequence
        score_aligned = alignment.score_aligned

        position_alleles = {}
        for msa_pos in variant_msa_positions:
            if msa_pos < len(aligned_seq):
                allele = aligned_seq[msa_pos]
                if score_aligned and msa_pos < len(score_aligned):
                    if score_aligned[msa_pos] == '=':
                        allele = consensus_aligned[msa_pos]
                position_alleles[msa_pos] = allele
            else:
                position_alleles[msa_pos] = '-'

        read_to_position_alleles[read_id] = position_alleles

    # Find best split — Phase 1: K=1 evaluation
    best_pos = None
    best_error = float('inf')
    best_qualifying = None
    best_non_qualifying = None
    all_positions = set(variant_msa_positions)
    all_read_ids = set(read_to_position_alleles.keys())
    cer_active = config.assumed_error_rate > 0 and config.total_specimen_reads > 0

    # Track positions that fail K=1 CER for potential K>1 fallback
    k1_failed = []  # list of (pos, minority_reads) for K>1 candidates

    for pos in variant_msa_positions:
        allele_groups = group_reads_by_single_position(
            read_to_position_alleles, pos, all_read_ids
        )
        qualifying, non_qualifying = filter_qualifying_haplotypes(
            allele_groups, total_reads, config.min_variant_count, config.min_variant_frequency
        )
        if len(qualifying) < 2:
            continue

        # CER gate: check statistical significance before expensive error calculation
        if cer_active:
            min_M = min(len(reads) for reads in qualifying.values())
            is_sig, p_star = is_variant_significant(
                M=min_M, N=config.total_specimen_reads,
                L=len(consensus),
                assumed_error_rate=config.assumed_error_rate,
                alpha=config.significance_level
            )
            if not is_sig:
                logging.debug(f"CER suppression at MSA pos {pos}: p*={p_star:.4f} < {config.assumed_error_rate}")
                # Save for K>1 fallback: record minority allele reads
                minority_allele = min(qualifying, key=lambda a: len(qualifying[a]))
                k1_failed.append((pos, qualifying[minority_allele].copy()))
                continue

        error = calculate_within_cluster_error(qualifying, read_to_position_alleles, {pos}, all_positions)
        if error < best_error:
            best_error = error
            best_pos = pos
            best_qualifying = qualifying
            best_non_qualifying = non_qualifying

    # Phase 2: K>1 correlation-based fallback (only if no K=1 candidates found)
    if best_pos is None and len(k1_failed) >= 2 and cer_active:
        best_k_combo = _find_correlated_positions(
            k1_failed, msa_to_consensus_pos, read_to_position_alleles,
            all_read_ids, total_reads, consensus, config
        )
        if best_k_combo is not None:
            best_pos, best_qualifying, best_non_qualifying, best_error = best_k_combo

    if best_pos is None:
        return [(path, consensus, read_ids)], set()

    # Collect deferred reads
    all_deferred = set()
    for allele, reads in best_non_qualifying.items():
        all_deferred.update(reads)

    # Recurse
    all_leaves = []
    for allele, sub_read_ids in sorted(best_qualifying.items()):
        # K>1 alleles are tuples (e.g. ('A', 'G', 'T')); convert to string for path
        allele_str = ''.join(allele) if isinstance(allele, tuple) else allele
        new_path = path + [allele_str]
        sub_leaves, sub_deferred = _recursive_phase_cluster_standalone(
            sub_read_ids, read_sequences, new_path, depth + 1, config
        )
        all_leaves.extend(sub_leaves)
        all_deferred.update(sub_deferred)

    return all_leaves, all_deferred


def _phase_reads_by_variants_standalone(
    cluster_read_ids: Set[str],
    sequences: Dict[str, str],
    variant_positions: List[Dict],
    config: ClusterProcessingConfig
) -> List[Tuple[str, Set[str]]]:
    """Phase reads into haplotypes. Standalone version for workers."""
    if not variant_positions:
        return [(None, cluster_read_ids)]

    try:
        read_sequences = {rid: sequences[rid] for rid in cluster_read_ids if rid in sequences}
        if not read_sequences:
            return [(None, cluster_read_ids)]

        logging.debug(f"Recursive phasing with MSA regeneration: {len(variant_positions)} initial variants, {len(read_sequences)} reads")

        leaves, deferred = _recursive_phase_cluster_standalone(
            set(read_sequences.keys()), read_sequences, [], 0, config
        )

        if len(leaves) <= 1 and not deferred:
            if leaves:
                path, consensus, reads = leaves[0]
                if not path:
                    return [(None, cluster_read_ids)]
            return [(None, cluster_read_ids)]

        if len(leaves) == 1:
            return [(None, cluster_read_ids)]

        logging.debug(f"Recursive phasing: {len(leaves)} leaf haplotypes, {len(deferred)} deferred reads")

        # Reassign deferred reads
        if deferred:
            leaf_reads_updated = {tuple(path): set(reads) for path, consensus, reads in leaves}
            leaf_consensuses = {tuple(path): consensus for path, consensus, reads in leaves}

            for read_id in deferred:
                if read_id not in read_sequences:
                    continue

                read_seq = read_sequences[read_id]
                min_distance = float('inf')
                nearest_path = None

                for path_tuple, consensus in leaf_consensuses.items():
                    if not consensus:
                        continue
                    result = edlib.align(read_seq, consensus)
                    distance = result['editDistance']
                    if distance < min_distance:
                        min_distance = distance
                        nearest_path = path_tuple

                if nearest_path is not None:
                    leaf_reads_updated[nearest_path].add(read_id)

            # Log and decide on phasing
            total_phased = sum(len(reads) for reads in leaf_reads_updated.values())
            if total_phased < 2:
                return [(None, cluster_read_ids)]

            logging.debug(f"Phasing decision: SPLITTING cluster into {len(leaf_reads_updated)} haplotypes")

            return [('-'.join(path), reads) for path, reads in sorted(leaf_reads_updated.items(), key=lambda x: -len(x[1]))]
        else:
            logging.debug(f"Phasing decision: SPLITTING cluster into {len(leaves)} haplotypes")
            return [('-'.join(path), reads) for path, consensus, reads in sorted(leaves, key=lambda x: -len(x[2]))]

    except Exception as e:
        logging.warning(f"Failed to phase reads: {e}")
        return [(None, cluster_read_ids)]


def _process_cluster_worker(args) -> Tuple[List[Dict], Set[str]]:
    """Worker function for parallel cluster processing.

    Must be at module level for ProcessPoolExecutor pickling.

    Args:
        args: Tuple of (initial_idx, cluster_ids, sequences, qualities, config)

    Returns:
        Tuple of (subclusters, discarded_read_ids)
    """
    initial_idx, cluster_ids, sequences, qualities, config = args

    subclusters = []
    discarded_ids = set()

    # Sort by quality
    sorted_ids = sorted(cluster_ids, key=lambda x: (-qualities.get(x, 0), x))

    # Generate consensus and MSA
    cluster_seqs = {seq_id: sequences[seq_id] for seq_id in sorted_ids}
    result = _run_spoa_for_cluster_worker(cluster_seqs, config.disable_homopolymer_equivalence)

    if result is None:
        logging.warning(f"Initial cluster {initial_idx}: Failed to generate consensus, skipping")
        # Track these reads as discarded since we couldn't generate consensus
        discarded_ids.update(cluster_ids)
        return subclusters, discarded_ids

    consensus = result.consensus
    msa = result.msa_string
    alignments = result.alignments
    msa_to_consensus_pos = result.msa_to_consensus_pos

    _calculate_read_identity_standalone(alignments, consensus)

    cluster = set(cluster_ids)

    # Outlier removal
    if config.outlier_identity_threshold is not None:
        keep_ids, outlier_ids = _identify_outlier_reads_standalone(
            alignments, consensus, cluster, config.outlier_identity_threshold
        )

        if outlier_ids:
            if len(cluster) == 2 and len(outlier_ids) == 1:
                logging.debug(f"Initial cluster {initial_idx}: Split 2-read cluster due to 1 outlier")
                for read_id in cluster:
                    subclusters.append({
                        'read_ids': {read_id},
                        'initial_cluster_num': initial_idx,
                        'allele_combo': 'single-read-split'
                    })
                return subclusters, discarded_ids

            logging.debug(f"Initial cluster {initial_idx}: Removing {len(outlier_ids)}/{len(cluster)} outlier reads, "
                         f"regenerating consensus")

            discarded_ids.update(outlier_ids)
            cluster = cluster - outlier_ids

            # Regenerate consensus
            sorted_ids_filtered = sorted(cluster, key=lambda x: (-qualities.get(x, 0), x))
            cluster_seqs = {seq_id: sequences[seq_id] for seq_id in sorted_ids_filtered}
            result = _run_spoa_for_cluster_worker(cluster_seqs, config.disable_homopolymer_equivalence)

            if result is not None:
                consensus = result.consensus
                msa = result.msa_string
                alignments = result.alignments
                msa_to_consensus_pos = result.msa_to_consensus_pos
                _calculate_read_identity_standalone(alignments, consensus)

    # Detect variants
    variant_positions = []
    if consensus and alignments and config.enable_secondpass_phasing:
        variant_positions = _detect_variant_positions_standalone(
            alignments, consensus, msa_to_consensus_pos,
            config.min_variant_frequency, config.min_variant_count
        )

        if variant_positions:
            logging.debug(f"Initial cluster {initial_idx}: Detected {len(variant_positions)} variant positions")

    # Phase reads
    phased_haplotypes = _phase_reads_by_variants_standalone(
        cluster, sequences, variant_positions, config
    )

    for haplotype_idx, (allele_combo, haplotype_reads) in enumerate(phased_haplotypes):
        subclusters.append({
            'read_ids': haplotype_reads,
            'initial_cluster_num': initial_idx,
            'allele_combo': allele_combo
        })

    return subclusters, discarded_ids


def _trim_primers_standalone(sequence: str, primers: Optional[List[Tuple[str, str]]]) -> Tuple[str, List[str]]:
    """Trim primers from start and end of sequence. Standalone version for workers."""
    if not primers:
        return sequence, []

    found_primers = []
    trimmed_seq = sequence

    # Look for primer at 5' end
    best_start_dist = float('inf')
    best_start_primer = None
    best_start_end = None

    for primer_name, primer_seq in primers:
        k = len(primer_seq) // 4  # Allow ~25% errors
        search_region = sequence[:len(primer_seq) * 2]

        result = edlib.align(primer_seq, search_region, task="path", mode="HW", k=k)

        if result["editDistance"] != -1:
            dist = result["editDistance"]
            if dist < best_start_dist:
                best_start_dist = dist
                best_start_primer = primer_name
                best_start_end = result["locations"][0][1] + 1

    if best_start_primer:
        found_primers.append(f"5'-{best_start_primer}")
        trimmed_seq = trimmed_seq[best_start_end:]

    # Look for primer at 3' end
    best_end_dist = float('inf')
    best_end_primer = None
    best_end_start = None

    for primer_name, primer_seq in primers:
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
        found_primers.append(f"3'-{best_end_primer}")
        trimmed_seq = trimmed_seq[:best_end_start]

    return trimmed_seq, found_primers


def _generate_cluster_consensus_worker(args) -> Dict:
    """Worker function for parallel final consensus generation.

    Must be at module level for ProcessPoolExecutor pickling.

    Args:
        args: Tuple of (final_idx, cluster_read_ids, sequences, qualities, config)
              - final_idx: 1-based cluster index
              - cluster_read_ids: Set of read IDs in this cluster
              - sequences: Dict mapping read_id -> sequence
              - qualities: Dict mapping read_id -> mean quality score
              - config: ConsensusGenerationConfig

    Returns:
        Dict with all computed results for this cluster
    """
    final_idx, cluster_read_ids, sequences, qualities, config = args

    cluster = cluster_read_ids
    actual_size = len(cluster)

    # Sort all cluster reads by quality for consistent output ordering
    sorted_cluster_ids = sorted(
        cluster,
        key=lambda x: (-qualities.get(x, 0), x)
    )

    # Sample sequences for final consensus generation if needed
    if len(cluster) > config.max_sample_size:
        sorted_sampled_ids = sorted_cluster_ids[:config.max_sample_size]
        sampled_ids = set(sorted_sampled_ids)
    else:
        sorted_sampled_ids = sorted_cluster_ids
        sampled_ids = cluster

    # Generate final consensus and MSA
    sampled_seqs = {seq_id: sequences[seq_id] for seq_id in sorted_sampled_ids}
    result = _run_spoa_for_cluster_worker(sampled_seqs, config.disable_homopolymer_equivalence)

    # Calculate final identity metrics
    rid, rid_min = None, None
    consensus = None
    msa = None
    iupac_count = 0
    trimmed_consensus = None
    found_primers = None

    if result is not None:
        consensus = result.consensus
        msa = result.msa_string
        alignments = result.alignments
        rid, rid_min = _calculate_read_identity_standalone(alignments, consensus)

        if consensus:
            # Apply IUPAC ambiguity calling for unphased variant positions
            if config.enable_iupac_calling:
                consensus, iupac_count, iupac_details = call_iupac_ambiguities(
                    consensus=consensus,
                    alignments=result.alignments,
                    msa_to_consensus_pos=result.msa_to_consensus_pos,
                    min_variant_frequency=config.min_ambiguity_frequency,
                    min_variant_count=config.min_ambiguity_count
                )

            # Perform primer trimming
            if config.primers:
                trimmed_consensus, found_primers = _trim_primers_standalone(consensus, config.primers)

    return {
        'final_idx': final_idx,
        'cluster': cluster,
        'actual_size': actual_size,
        'consensus': consensus,
        'trimmed_consensus': trimmed_consensus,
        'found_primers': found_primers,
        'rid': rid,
        'rid_min': rid_min,
        'msa': msa,
        'sampled_ids': sampled_ids,
        'sorted_cluster_ids': sorted_cluster_ids,
        'sorted_sampled_ids': sorted_sampled_ids,
        'iupac_count': iupac_count
    }
