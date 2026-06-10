"""Group-level full consensus builder for speconsense-summarize.

Implements ``--enable-full-consensus``: per identity group, draw a
size-weighted, quality-sorted, noise-scrubbed read sample across the
pre-merge core variants in the bucket and build a single dominant-path
consensus from it. The output is named ``-{gid}-full`` and emitted
alongside the per-variant FASTAs. The artifact is meant as a query
substrate for BLAST against legacy unphased references — it reproduces
what NGSpeciesID would have produced on the input FASTQ if the original
reads had been pre-filtered for noise.

The deprecated ``--enable-full-consensus`` (removed in commit 4dc3797)
unioned indels across pre-merge consensuses, which produced a sequence
no real haplotype contained and regressed legacy BLAST matching by ~58
percentage points at raw-BLAST ≥99%. This implementation replaces that
policy with majority-wins gaps and IUPAC ambiguity only at columns that
clear ``min_ambiguity_frequency`` (inherited from the core run).
"""

from __future__ import annotations

import logging
import os
import statistics
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from speconsense.types import ConsensusInfo
from speconsense.msa import extract_alignments_from_msa

from .analysis import run_spoa_for_cluster_metrics
from .merging import (
    _resolve_contributor_fastq,
    build_full_consensus_from_msa,
)


def _gate_variants(
    variants: List[ConsensusInfo],
    min_ambiguity_frequency: float,
) -> List[ConsensusInfo]:
    """Sort variants by size descending and keep those whose size is at
    least ``min_ambiguity_frequency`` of the running total."""
    sorted_variants = sorted(variants, key=lambda v: v.size, reverse=True)
    if not sorted_variants:
        return []
    gated = [sorted_variants[0]]
    running_total = sorted_variants[0].size
    for variant in sorted_variants[1:]:
        if variant.size >= min_ambiguity_frequency * running_total:
            gated.append(variant)
            running_total += variant.size
        else:
            break
    return gated


def _allocate_reads(sizes: List[int], budget: int) -> List[int]:
    """Size-weighted allocation summing to ``budget`` (or to ``sum(sizes)``
    if smaller). Fractional remainders go to the largest-size variants."""
    total = sum(sizes)
    if total <= 0 or budget <= 0:
        return [0] * len(sizes)
    if total <= budget:
        return list(sizes)

    raw = [budget * size / total for size in sizes]
    floored = [int(x) for x in raw]
    assigned = sum(floored)
    remainder = budget - assigned
    if remainder > 0:
        fractional = sorted(
            range(len(sizes)),
            key=lambda i: (-(raw[i] - floored[i]), -sizes[i]),
        )
        for i in fractional[:remainder]:
            floored[i] += 1
    return floored


def _load_quality_sorted_reads(
    variant: ConsensusInfo,
    fastq_lookup: Optional[Dict[str, List[str]]],
    limit: int,
) -> List[SeqRecord]:
    """Return the top-``limit`` reads from this variant's cluster_debug
    FASTQ, sorted by mean Phred descending. Empty list if the FASTQ is
    unavailable.
    """
    if limit <= 0:
        return []
    path = _resolve_contributor_fastq(variant, fastq_lookup)
    if not path or not os.path.exists(path):
        return []
    try:
        records = list(SeqIO.parse(path, "fastq"))
    except Exception as exc:
        logging.debug(
            f"build_group_full_consensus: failed to parse {path}: {exc}"
        )
        return []

    def mean_phred(record: SeqRecord) -> float:
        qualities = record.letter_annotations.get("phred_quality") or []
        if not qualities:
            return 0.0
        return statistics.fmean(qualities)

    records.sort(key=mean_phred, reverse=True)
    return records[:limit]


def _is_cross_primer(variants: List[ConsensusInfo]) -> bool:
    primer_sets = {frozenset(v.primers or []) for v in variants}
    return len(primer_sets) > 1


def build_group_full_consensus(
    final_gid: int,
    group_members: List[ConsensusInfo],
    pass_track_count: int,
    fastq_lookup: Optional[Dict[str, List[str]]],
    min_ambiguity_frequency: float,
    max_sample_size: int,
    specimen_base: str,
    primary_file_path: str,
    group_size_total: Optional[int] = None,
    global_size_total: Optional[int] = None,
    min_position_frequency: float = 0.5,
    min_position_count: int = 3,
) -> Optional[Tuple[ConsensusInfo, List[SeqRecord]]]:
    """Build the ``-{gid}-full`` consensus and return it with the sampled
    reads, or ``None`` when the guards fire (single pass-track variant,
    single gated variant, no reads available, SPOA failure).

    Args:
        final_gid: The (post-conflation) identity group id this consensus
            represents. Becomes the ``-{gid}`` portion of the output name.
        group_members: Pre-MSA-merge core variants in this bucket. May
            include records from absorbed cross-primer-conflated groups.
        pass_track_count: Number of variants this group emits to the
            pass track after summarize's selection. The builder bails
            when ``< 2`` (single-variant case adds no information).
        fastq_lookup: Specimen -> list of cluster_debug FASTQ paths, used
            by ``_resolve_contributor_fastq`` to locate per-variant reads.
        min_ambiguity_frequency: Per-specimen IUPAC frequency threshold
            from the core run's metadata. Reused as the variant-gate
            threshold so a single number governs both "is this signal
            strong enough to include?" and "is this minor base strong
            enough to call IUPAC?".
        max_sample_size: Total read budget from the core run's metadata.
            Default fallback handled by the caller.
        specimen_base: Specimen name with any cluster suffix stripped.
        primary_file_path: File path for the emitted ConsensusInfo
            (matches the per-specimen ``-all.fasta`` source for downstream
            grouping; not used for FASTQ resolution).
    """
    if pass_track_count < 2:
        return None

    gated = _gate_variants(group_members, min_ambiguity_frequency)
    if len(gated) < 2:
        return None

    sizes = [v.size for v in gated]
    allocations = _allocate_reads(sizes, max_sample_size)

    sampled: List[SeqRecord] = []
    for variant, alloc in zip(gated, allocations):
        sampled.extend(_load_quality_sorted_reads(variant, fastq_lookup, alloc))

    if len(sampled) < 2:
        logging.debug(
            f"build_group_full_consensus: insufficient sampled reads for gid={final_gid} "
            f"(got {len(sampled)}); skipping -full"
        )
        return None

    sequences = {f"r{idx}": str(record.seq) for idx, record in enumerate(sampled)}

    alignment_mode = 0 if _is_cross_primer(gated) else 1
    msa_result = run_spoa_for_cluster_metrics(
        sequences,
        alignment_mode=alignment_mode,
    )
    if msa_result is None:
        logging.warning(
            f"build_group_full_consensus: SPOA failed for gid={final_gid}; skipping -full"
        )
        return None

    aligned_reads, _consensus_str, _pos_map = extract_alignments_from_msa(
        msa_result.msa_string,
        enable_homopolymer_normalization=False,
    )

    aligned_records = [
        SeqRecord(Seq(alignment.aligned_sequence), id=alignment.read_id, description="")
        for alignment in aligned_reads
    ]

    if len(aligned_records) < 2:
        logging.warning(
            f"build_group_full_consensus: parsed {len(aligned_records)} alignment rows for gid={final_gid}; "
            "skipping -full"
        )
        return None

    consensus_seq, snp_count = build_full_consensus_from_msa(
        aligned_records,
        min_ambiguity_frequency,
        min_position_frequency=min_position_frequency,
        min_position_count=min_position_count,
    )
    if not consensus_seq:
        return None

    primer_set = sorted({p for v in gated for p in (v.primers or [])}) or None
    total_size = sum(v.size for v in gated)
    total_ric = len(sampled)

    full_record = ConsensusInfo(
        sample_name=f"{specimen_base}-{final_gid}-full",
        cluster_id=f"{final_gid}-full",
        sequence=consensus_seq,
        ric=total_ric,
        size=total_size,
        file_path=primary_file_path,
        snp_count=snp_count if snp_count > 0 else None,
        primers=primer_set,
        raw_ric=None,
        raw_len=None,
        rid=None,
        rid_min=None,
        merge_indel_count=None,
        cer_factor=None,
        err_factor=None,
        err_factor_obs_sum=None,
        err_factor_exp_sum=None,
        err_factor_cols=None,
        group_rank=final_gid,
        variant_rank=None,
        group_size_total=group_size_total,
        global_size_total=global_size_total,
    )
    return full_record, sampled
