"""MSA-based variant merging for speconsense-summarize.

Provides functions for finding and merging compatible variants within HAC groups
using exhaustive subset evaluation with SPOA multiple sequence alignment.
"""

import itertools
import logging
import os
from typing import List, Tuple, Dict, Optional
from collections import defaultdict

from Bio import SeqIO

from speconsense.types import ConsensusInfo, OverlapMergeInfo
from speconsense.msa import compute_cluster_err_factor
from speconsense.context import classify_pairwise_differences
from speconsense.qctx import get_qctx
from speconsense.significance import compute_cer_factor

from .iupac import merge_bases_to_iupac, primers_are_same
from .analysis import (
    run_spoa_msa,
    run_spoa_for_cluster_metrics,
    analyze_msa_columns,
    analyze_msa_columns_overlap_aware,
    MAX_MSA_MERGE_VARIANTS,  # Kept for backward compatibility
    compute_merge_batch_size,
)


def generate_all_subsets_by_size(variants: List[ConsensusInfo]) -> List[Tuple[int, ...]]:
    """
    Generate all possible non-empty subsets of variant indices.
    Returns subsets in descending order by total cluster size.

    This exhaustive approach guarantees finding the globally optimal merge
    when the number of variants is small (<= MAX_MSA_MERGE_VARIANTS).

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


def is_compatible_subset(variant_stats: dict, args, prior_positions: dict = None) -> bool:
    """
    Check if variant statistics are within merge limits.

    By default, homopolymer indels are ignored (treated as compatible) to match
    adjusted-identity homopolymer normalization semantics where AAA ~ AAAA.
    Only structural indels count against the limits.

    When --disable-homopolymer-equivalence is set, homopolymer indels are treated
    the same as structural indels and count against merge limits.

    Args:
        variant_stats: Statistics from MSA analysis (snp_count, indel counts, etc.)
        args: Command-line arguments with merge parameters
        prior_positions: Optional dict with cumulative counts from prior merge rounds
                        {'snp_count': N, 'indel_count': M} - these are added to
                        current stats when checking limits for iterative merging
    """
    if prior_positions is None:
        prior_positions = {'snp_count': 0, 'indel_count': 0}

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

    # Check total position count (including prior merge rounds)
    total_positions = (variant_stats['snp_count'] + prior_positions['snp_count'] +
                      indel_count + prior_positions['indel_count'])
    if total_positions > args.merge_position_count:
        return False

    return True


def _build_merged_consensus_info(
    consensus_seq: list, snp_count: int, variants: List[ConsensusInfo]
) -> ConsensusInfo:
    """Assemble a ConsensusInfo from column-voting results and source variants.

    Field handling summary:
        - ``sequence``: column-voted, passed in as ``consensus_seq``.
        - ``size``, ``ric``: summed across contributors.
        - ``raw_ric``, ``raw_len``: flattened merge history, sorted desc.
        - ``snp_count``: count of new IUPAC positions from this column vote;
          the caller in ``merge_group_with_msa`` overwrites with the cumulative
          total across iterative rounds. Cumulation can overcount when the
          same physical position becomes ambiguous in multiple rounds; the
          final consensus's ``ambig`` (count of non-ACGT characters) is the
          canonical IUPAC site count.
        - ``primers``: union across contributors, sorted. Aligns with the
          set-based comparator in ``primers_are_same``. Cross-primer merges
          surface as multi-primer lists in the output header.
        - ``rid``, ``rid_min``, ``err_factor``, ``cer_factor``: inherited from
          the largest contributor here; ``merge_group_with_msa`` recomputes
          ``rid``, ``rid_min``, and ``err_factor`` from the merged cluster's
          SPOA MSA over the union of contributor reads, and recomputes
          ``cer_factor`` for same-primer merges (sets it to None for
          cross-primer merges).
        - ``group_rank``, ``variant_rank``, ``cluster_id``, ``file_path``,
          ``sample_name``: inherited from the largest contributor; the naming
          pass in ``process_single_specimen`` may overwrite ``sample_name``.

    Args:
        consensus_seq: List of consensus characters from column voting.
        snp_count: Number of ambiguous (multi-base) positions in this vote.
        variants: Source ConsensusInfo objects that were merged.

    Returns:
        ConsensusInfo with merged metadata. The merge-time recompute of
        ``rid``/``rid_min``/``err_factor``/``cer_factor`` happens in the
        caller after this helper returns.
    """
    consensus_sequence = ''.join(consensus_seq)
    total_size = sum(v.size for v in variants)
    total_ric = sum(v.ric for v in variants)

    # Collect RiC values, preserving any prior merge history
    raw_ric_values = []
    for v in variants:
        if v.raw_ric:
            raw_ric_values.extend(v.raw_ric)
        else:
            raw_ric_values.append(v.ric)
    raw_ric_values = sorted(raw_ric_values, reverse=True) if len(variants) > 1 else None

    # Collect lengths, preserving any prior merge history
    raw_len_values = []
    for v in variants:
        if v.raw_len:
            raw_len_values.extend(v.raw_len)
        else:
            raw_len_values.append(len(v.sequence))
    raw_len_values = sorted(raw_len_values, reverse=True) if len(variants) > 1 else None

    # Union primer names across contributors, sorted. Matches the set-based
    # semantics of primers_are_same and gives cross-primer merges an honest
    # multi-primer header (e.g., 5'-ITS1F,5'-ITS3,3'-ITS4_RC for a 2-pair
    # primer-pool conflation). Single-pair merges produce the same primer
    # list they would have inherited from the largest contributor.
    primer_set = {p for v in variants for p in (v.primers or [])}
    merged_primers = sorted(primer_set) if primer_set else None

    largest_variant = max(variants, key=lambda v: v.size)

    return ConsensusInfo(
        sample_name=largest_variant.sample_name,
        cluster_id=largest_variant.cluster_id,
        sequence=consensus_sequence,
        ric=total_ric,
        size=total_size,
        file_path=largest_variant.file_path,
        snp_count=snp_count if snp_count > 0 else None,
        primers=merged_primers,
        raw_ric=raw_ric_values,
        raw_len=raw_len_values,
        rid=largest_variant.rid,
        rid_min=largest_variant.rid_min,
        cer_factor=largest_variant.cer_factor,
        err_factor=largest_variant.err_factor,
        err_factor_obs_sum=largest_variant.err_factor_obs_sum,
        err_factor_exp_sum=largest_variant.err_factor_exp_sum,
        err_factor_cols=largest_variant.err_factor_cols,
        group_rank=largest_variant.group_rank,
        variant_rank=largest_variant.variant_rank,
    )


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
                # Multiple bases - generate IUPAC code (expanding any existing IUPAC codes)
                represented_bases = set(base_votes.keys())
                iupac_code = merge_bases_to_iupac(represented_bases)
                consensus_seq.append(iupac_code)
                snp_count += 1
        # else: majority wants gap, omit position

    return _build_merged_consensus_info(consensus_seq, snp_count, variants)


def create_overlap_consensus_from_msa(aligned_seqs: List, variants: List[ConsensusInfo]) -> ConsensusInfo:
    """
    Generate consensus from MSA where sequences may have different lengths.

    For overlap merging (primer pools with different endpoints):
    - In overlap region: Use size-weighted majority voting
    - In non-overlap regions: Keep content from whichever sequence(s) have it

    This produces a consensus spanning the union of all input sequences.

    Args:
        aligned_seqs: MSA sequences with gaps as '-'
        variants: Original ConsensusInfo objects (for size weighting)

    Returns:
        ConsensusInfo with merged consensus sequence spanning full length
    """
    consensus_seq = []
    snp_count = 0
    alignment_length = len(aligned_seqs[0].seq)

    # Find content region for each sequence
    content_regions = []
    for seq in aligned_seqs:
        seq_str = str(seq.seq)
        first_base = next((i for i, c in enumerate(seq_str) if c != '-'), 0)
        last_base = alignment_length - 1 - next(
            (i for i, c in enumerate(reversed(seq_str)) if c != '-'), 0
        )
        content_regions.append((first_base, last_base))

    # Calculate overlap region
    overlap_start = max(start for start, _ in content_regions)
    overlap_end = min(end for _, end in content_regions)

    # Process each column
    for col_idx in range(alignment_length):
        column = [str(seq.seq[col_idx]) for seq in aligned_seqs]

        # Determine which sequences have content at this position
        seqs_with_content = []
        for i, (start, end) in enumerate(content_regions):
            if start <= col_idx <= end:
                seqs_with_content.append(i)

        if not seqs_with_content:
            # No sequence has content here (shouldn't happen in valid MSA)
            continue

        # Check if we're in the overlap region
        in_overlap = overlap_start <= col_idx <= overlap_end

        if in_overlap:
            # Overlap region: use size-weighted majority voting (like original)
            votes_with_size = [(column[i], variants[i].size) for i in seqs_with_content]

            votes = defaultdict(int)
            for base, size in votes_with_size:
                votes[base.upper()] += size

            gap_votes = votes.get('-', 0)
            base_votes = {b: v for b, v in votes.items() if b != '-'}
            total_base_votes = sum(base_votes.values())

            if total_base_votes > gap_votes:
                if len(base_votes) == 1:
                    consensus_seq.append(list(base_votes.keys())[0])
                else:
                    represented_bases = set(base_votes.keys())
                    iupac_code = merge_bases_to_iupac(represented_bases)
                    consensus_seq.append(iupac_code)
                    snp_count += 1
            # else: majority wants gap in overlap, omit position
        else:
            # Non-overlap region: keep content from available sequences
            # (don't let gap votes from sequences that don't extend here remove content)
            bases_only = [column[i] for i in seqs_with_content if column[i] != '-']

            if bases_only:
                # Weight by size for consistency
                votes = defaultdict(int)
                for i in seqs_with_content:
                    if column[i] != '-':
                        votes[column[i].upper()] += variants[i].size

                if len(votes) == 1:
                    consensus_seq.append(list(votes.keys())[0])
                else:
                    represented_bases = set(votes.keys())
                    iupac_code = merge_bases_to_iupac(represented_bases)
                    consensus_seq.append(iupac_code)
                    snp_count += 1

    return _build_merged_consensus_info(consensus_seq, snp_count, variants)


def merge_group_with_msa(variants: List[ConsensusInfo], args) -> Tuple[List[ConsensusInfo], Dict, int, List[OverlapMergeInfo]]:
    """
    Find largest mergeable subset of variants using MSA-based evaluation with exhaustive search.

    Algorithm:
    1. Process variants in batches of up to MAX_MSA_MERGE_VARIANTS
    2. For each batch, run SPOA MSA once
    3. Exhaustively evaluate ALL subsets by total size (descending)
    4. Merge the best compatible subset found
    5. Remove merged variants and repeat with remaining
    6. When overlap mode is enabled, iterate the entire process on merged results
       until no more merges happen (handles prefix+suffix+full scenarios)

    This approach guarantees optimal results when N <= MAX_MSA_MERGE_VARIANTS.
    For N > MAX, processes top MAX per round (potentially suboptimal globally).

    Iterative merging (overlap mode only):
    - After first pass, merged results are fed back for another round
    - Cumulative SNP/indel counts are tracked across rounds
    - Continues until no merges occur in a round

    Args:
        variants: List of ConsensusInfo from HAC group
        args: Command-line arguments with merge parameters

    Returns:
        (merged_variants, merge_traceability, potentially_suboptimal, overlap_merges) where:
        - merged_variants is list of merged ConsensusInfo objects
        - traceability maps merged names to original cluster names
        - potentially_suboptimal is 1 if group had >MAX variants, 0 otherwise
        - overlap_merges is list of OverlapMergeInfo for quality reporting
    """
    if len(variants) == 1:
        return variants, {}, 0, []

    # Compute batch size based on effort and group size
    effort = getattr(args, 'merge_effort_value', 10)  # Default to balanced
    batch_size = compute_merge_batch_size(len(variants), effort)

    # Track if this group is potentially suboptimal (too many variants for global optimum)
    potentially_suboptimal = 1 if len(variants) > batch_size else 0

    all_traceability = {}
    overlap_merges = []  # Track overlap merge events for quality reporting

    # For iterative merging in overlap mode, we may need multiple rounds
    current_variants = variants
    iteration = 0
    max_iterations = 10  # Safety limit to prevent infinite loops

    while iteration < max_iterations:
        iteration += 1

        # Sort variants by size (largest first)
        remaining_variants = sorted(current_variants, key=lambda v: v.size, reverse=True)
        merged_results = []
        merges_this_iteration = 0

        while remaining_variants:
            # Take up to batch_size candidates (dynamically computed based on effort and group size)
            candidates = remaining_variants[:batch_size]

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

            if iteration > 1:
                logging.debug(f"Iteration {iteration}: Evaluating {len(candidates)} variants "
                              f"(batch_size={batch_size}) for merging")
            else:
                logging.debug(f"Evaluating {len(candidates)} variants (batch_size={batch_size}, "
                              f"effort={effort}) for merging (exhaustive subset search)")

            # Determine if overlap mode should be used for this merge batch
            # Same primers -> use global mode (chimeras have same primers but different lengths)
            # Different primers -> use overlap mode (legitimate primer pool variation)
            all_same_primers = all(
                primers_are_same(candidates[0].primers, v.primers)
                for v in candidates[1:]
            ) if len(candidates) > 1 else True
            use_overlap_mode = args.min_merge_overlap > 0 and not all_same_primers

            if args.min_merge_overlap > 0 and all_same_primers and len(candidates) > 1:
                # Log when primer constraint prevents overlap merging
                primer_str = ','.join(candidates[0].primers) if candidates[0].primers else 'unknown'
                logging.debug(f"Same primers [{primer_str}] detected - using global alignment instead of overlap")

            # Run SPOA MSA on candidates
            # Use local alignment mode (0) for overlap merging to get clean terminal gaps
            # Use global alignment mode (1) for standard same-length merging
            sequences = [v.sequence for v in candidates]
            spoa_mode = 0 if use_overlap_mode else 1
            aligned_seqs = run_spoa_msa(sequences, alignment_mode=spoa_mode)

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
                hp_min = getattr(args, 'hp_normalization_length', 1)
                if use_overlap_mode:
                    # Use overlap-aware analysis for primer pool scenarios
                    original_lengths = [len(v.sequence) for v in subset_variants]
                    variant_stats = analyze_msa_columns_overlap_aware(
                        subset_aligned, args.min_merge_overlap, original_lengths,
                        min_hp_length=hp_min,
                    )

                    # Check overlap requirement
                    shorter_len = min(original_lengths)
                    effective_threshold = min(args.min_merge_overlap, shorter_len)
                    if variant_stats['overlap_bp'] < effective_threshold:
                        # Insufficient overlap - skip this subset
                        continue
                else:
                    # Use standard analysis
                    variant_stats = analyze_msa_columns(subset_aligned, min_hp_length=hp_min)

                # Calculate cumulative positions from input sequences (for iterative merging)
                # Each sequence may carry positions from prior merges
                prior_snps = sum(v.snp_count or 0 for v in subset_variants)
                prior_indels = sum(v.merge_indel_count or 0 for v in subset_variants)
                prior_positions = {'snp_count': prior_snps, 'indel_count': prior_indels}

                # Check compatibility against merge limits (including cumulative positions)
                if is_compatible_subset(variant_stats, args, prior_positions):
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
                        iter_prefix = f"Iteration {iteration}: " if iteration > 1 else ""
                        if use_overlap_mode:
                            # Include prefix/suffix extension info for overlap merges
                            prefix_bp = variant_stats.get('prefix_bp', 0)
                            suffix_bp = variant_stats.get('suffix_bp', 0)
                            logging.info(f"{iter_prefix}Found mergeable subset of {len(subset_indices)} variants "
                                       f"(overlap={variant_stats.get('overlap_bp', 'N/A')}bp, "
                                       f"prefix={prefix_bp}bp, suffix={suffix_bp}bp): {variant_desc}")

                            # DEBUG: Show span details for each sequence in the merge
                            content_regions = variant_stats.get('content_regions', [])
                            if content_regions:
                                spans = [f"seq{i+1}=({s},{e})" for i, (s, e) in enumerate(content_regions)]
                                logging.debug(f"Merge spans: {', '.join(spans)}")
                        else:
                            logging.info(f"{iter_prefix}Found mergeable subset of {len(subset_indices)} variants: {variant_desc}")

                        # Calculate total positions for cumulative tracking
                        # Total = prior positions from input sequences + new positions from this merge
                        if args.disable_homopolymer_equivalence:
                            this_merge_indels = variant_stats['structural_indel_count'] + variant_stats['homopolymer_indel_count']
                        else:
                            this_merge_indels = variant_stats['structural_indel_count']
                        total_snps = prior_snps + variant_stats['snp_count']
                        total_indels = prior_indels + this_merge_indels

                    # Create merged consensus
                    if len(subset_indices) == 1:
                        # Single variant - use directly, preserving raw_ric and other metadata
                        merged_consensus = subset_variants[0]
                    elif use_overlap_mode:
                        # Use overlap-aware consensus generation
                        merged_consensus = create_overlap_consensus_from_msa(
                            subset_aligned, subset_variants
                        )
                    else:
                        merged_consensus = create_consensus_from_msa(
                            subset_aligned, subset_variants
                        )

                    # Update merged consensus with cumulative position counts for iterative tracking
                    if len(subset_indices) > 1:
                        merged_consensus = merged_consensus._replace(
                            snp_count=total_snps if total_snps > 0 else None,
                            merge_indel_count=total_indels if total_indels > 0 else None
                        )

                    # Track merge provenance - expand any intermediate merges
                    # so we always trace back to the original cluster names
                    original_clusters = []
                    for v in subset_variants:
                        if v.sample_name in all_traceability:
                            # This variant was itself merged, expand to its originals
                            original_clusters.extend(all_traceability[v.sample_name])
                        else:
                            original_clusters.append(v.sample_name)
                    traceability = {
                        merged_consensus.sample_name: original_clusters
                    }
                    all_traceability.update(traceability)

                    # Track overlap merge for quality reporting
                    if use_overlap_mode and len(subset_indices) > 1:
                        # Extract specimen name by stripping the gid.vid suffix
                        from .io import strip_cluster_suffix
                        specimen = strip_cluster_suffix(merged_consensus.sample_name)
                        overlap_merges.append(OverlapMergeInfo(
                            specimen=specimen,
                            iteration=iteration,
                            input_clusters=[v.sample_name for v in subset_variants],
                            input_lengths=[len(v.sequence) for v in subset_variants],
                            input_rics=[v.ric for v in subset_variants],
                            overlap_bp=variant_stats.get('overlap_bp', 0),
                            prefix_bp=variant_stats.get('prefix_bp', 0),
                            suffix_bp=variant_stats.get('suffix_bp', 0),
                            output_length=len(merged_consensus.sequence)
                        ))

                    # Add merged consensus to results
                    merged_results.append(merged_consensus)

                    # Remove merged variants from remaining pool
                    for v in subset_variants:
                        if v in remaining_variants:
                            remaining_variants.remove(v)

                    merged_this_round = True
                    if len(subset_indices) > 1:
                        merges_this_iteration += 1
                    break

            # If no merge found, keep largest variant as-is and continue
            if not merged_this_round:
                logging.debug(f"No compatible merge found for largest variant (size={candidates[0].size})")
                merged_results.append(candidates[0])
                remaining_variants.remove(candidates[0])

        # Check if we should do another iteration (overlap mode only)
        if args.min_merge_overlap > 0 and merges_this_iteration > 0 and len(merged_results) > 1:
            # More merges might be possible with the new merged sequences
            # Cumulative positions are tracked per-sequence via snp_count and merge_indel_count
            logging.debug(f"Iteration {iteration} complete: {merges_this_iteration} merges, "
                         f"{len(merged_results)} variants remaining, trying another round")
            current_variants = merged_results
        else:
            # No more iterations needed
            if iteration > 1:
                logging.debug(f"Iterative merging complete after {iteration} iterations")
            break

    return merged_results, all_traceability, potentially_suboptimal, overlap_merges


# ---------------------------------------------------------------------------
# Merge-time metric recomputation
#
# When ``merge_group_with_msa`` collapses ≥2 contributors into a single
# record, ``_build_merged_consensus_info`` inherits ``rid``, ``rid_min``,
# ``err_factor`` from the largest contributor. Those values describe the
# largest contributor's pre-merge cluster, not the merged cluster. The
# helpers below re-derive them from the union of contributor reads.
#
# Cross-primer merges (different-primer contributors conflated by overlap)
# remain meaningful for ``rid``/``rid_min``/``err_factor`` because each
# metric is per-position-coverage-aware: reads with terminal gaps simply
# don't contribute to columns they don't cover. ``cer_factor`` is more
# delicate — the noise model wasn't designed for merged-locus candidates —
# so the cross-primer branch sets it to None (treated as "no valid peer
# comparison" by summarize's routing).
# ---------------------------------------------------------------------------


def _resolve_contributor_fastq(
    contributor: ConsensusInfo,
    fastq_lookup: Optional[Dict[str, List[str]]],
) -> Optional[str]:
    """Return the cluster_debug FASTQ path for ``contributor``.

    Matches via the ``{sample_name}-RiC{ric}-`` substring that core uses
    for cluster_debug filenames. Returns None when the lookup is missing
    or no file matches.
    """
    if not fastq_lookup or not contributor.sample_name:
        return None
    # Lazy import to avoid a circular dependency with io.
    from .io import strip_cluster_suffix
    specimen_name = strip_cluster_suffix(contributor.sample_name)
    if specimen_name == contributor.sample_name:
        return None
    files = fastq_lookup.get(specimen_name, [])
    if not files:
        return None
    needle = f"{contributor.sample_name}-RiC{contributor.ric}-"
    for path in files:
        if needle in os.path.basename(path):
            return path
    return None


def _load_contributor_reads(
    contributors: List[ConsensusInfo],
    fastq_lookup: Optional[Dict[str, List[str]]],
) -> Dict[str, str]:
    """Load the union of reads from each contributor's debug FASTQ.

    Read IDs are namespaced ``{cluster_id}:{read_id}`` so cross-cluster
    duplicates don't collide in the SPOA input dict.
    """
    sequences: Dict[str, str] = {}
    for contributor in contributors:
        fastq_path = _resolve_contributor_fastq(contributor, fastq_lookup)
        if not fastq_path or not os.path.exists(fastq_path):
            continue
        try:
            for record in SeqIO.parse(fastq_path, "fastq"):
                key = f"{contributor.cluster_id}:{record.id}"
                sequences[key] = str(record.seq)
        except Exception as e:
            logging.debug(
                f"Could not parse {fastq_path} during merge recompute: {e}"
            )
    return sequences


def _refresh_cluster_metrics(
    merged_record: ConsensusInfo,
    contributors: List[ConsensusInfo],
    fastq_lookup: Optional[Dict[str, List[str]]],
    qctx_table: Optional[Dict[str, float]],
    hp_normalization_length: int,
    disable_homopolymer_equivalence: bool,
) -> ConsensusInfo:
    """Recompute rid, rid_min, err_factor for a merged cluster.

    Runs SPOA over the union of contributor reads (loaded from each
    contributor's ``cluster_debug`` FASTQ) and replaces the inherited
    metrics with values derived from the merged read set:

    - ``rid``, ``rid_min`` come from per-read homopolymer-normalized
      identity against SPOA's re-derived consensus, the same calculation
      core uses for single-cluster rid (``extract_alignments_from_msa``).
    - ``err_factor`` and its raw sums come from
      ``msa.compute_cluster_err_factor`` against the same MSA.

    When the reads can't be loaded (missing debug FASTQs, lookup failure)
    or SPOA fails, the merged record is returned unchanged — the inherited
    values stay in place as a best-effort fallback.

    Interpretation notes:
    - The recomputed metrics describe the merged cluster's homogeneity
      against an MSA built from its union reads. SPOA's MSA-derived
      consensus may differ slightly from the shipped column-voted
      ``sequence``; the metrics describe cluster internal structure rather
      than agreement with the shipped consensus per se.
    - Cross-primer reads with terminal gaps contribute zero to columns
      they don't cover, so err_factor remains well-defined under primer
      conflation. rid uses HP-normalized edit distance over the read's
      MSA length and is correspondingly unaffected by uncovered terminal
      regions of the merged-locus consensus.
    """
    if len(contributors) <= 1:
        return merged_record

    sequences = _load_contributor_reads(contributors, fastq_lookup)
    if len(sequences) < 2:
        logging.debug(
            f"Skipping merge-metric recompute for {merged_record.sample_name}: "
            f"could not load union reads"
        )
        return merged_record

    msa_result = run_spoa_for_cluster_metrics(
        sequences,
        disable_homopolymer_equivalence=disable_homopolymer_equivalence,
        min_hp_length=hp_normalization_length,
    )
    if msa_result is None:
        logging.debug(
            f"Skipping merge-metric recompute for {merged_record.sample_name}: "
            f"SPOA failed on {len(sequences)} reads"
        )
        return merged_record

    # rid / rid_min from per-read homopolymer-normalized identity against
    # SPOA's MSA consensus.
    new_rid: Optional[float] = merged_record.rid
    new_rid_min: Optional[float] = merged_record.rid_min
    if msa_result.alignments and msa_result.consensus:
        consensus_len = len(msa_result.consensus)
        if consensus_len > 0:
            identities = [
                1.0 - (a.normalized_edit_distance / consensus_len)
                for a in msa_result.alignments
            ]
            if identities:
                new_rid = sum(identities) / len(identities)
                new_rid_min = min(identities)

    # err_factor against the same MSA, weighted by q_ctx where available.
    new_err = merged_record.err_factor
    new_obs = merged_record.err_factor_obs_sum
    new_exp = merged_record.err_factor_exp_sum
    new_cols = merged_record.err_factor_cols
    if qctx_table is not None:
        ef, obs, exp, cols = compute_cluster_err_factor(
            msa_result.msa_string, qctx_table
        )
        # Preserve inherited sums when err_factor can't be computed; replace
        # otherwise. cols/obs/exp are 0 in the failure case so guarding on
        # ef avoids silently zeroing the quality report's calibration.
        if ef is not None:
            new_err = ef
            new_obs = obs
            new_exp = exp
            new_cols = cols

    return merged_record._replace(
        rid=new_rid,
        rid_min=new_rid_min,
        err_factor=new_err,
        err_factor_obs_sum=new_obs,
        err_factor_exp_sum=new_exp,
        err_factor_cols=new_cols,
    )


def _count_independent_sites(seq: str) -> int:
    """Bonferroni site count for context-aware CER (cer_in_practice §2.3).

    Replicated from ``speconsense.core.clusterer._count_independent_sites``
    to avoid a cross-package private-helper import. Each maximal run of
    identical bases counts as one site regardless of run length.
    """
    if not seq:
        return 0
    n = 1
    for i in range(1, len(seq)):
        if seq[i] != seq[i - 1]:
            n += 1
    return n


def _refresh_cer_factor(
    merged_record: ConsensusInfo,
    contributors: List[ConsensusInfo],
    bucket_members: List[ConsensusInfo],
    qctx_table: Optional[Dict[str, float]],
    alpha: float,
    hp_min_length: int,
) -> ConsensusInfo:
    """Recompute cer_factor for a merged cluster against post-merge peers.

    Walks the surviving records in the merged record's bucket and, for each
    peer with strictly larger size, classifies the pairwise differences via
    ``classify_pairwise_differences``, looks up each variant event's q_ctx,
    and runs ``compute_cer_factor`` with N = bucket-total reads (post-merge),
    M = candidate size. The minimum factor across peers is the candidate's
    new cer_factor (same convention core uses in
    ``SpecimenClusterer._compute_cer_for_candidate``).

    Cross-primer merges (contributors had different primer sets) bypass this
    pipeline and get ``cer_factor=None``. The CER noise model is per-locus
    by construction and isn't well-defined when the merged candidate spans
    multiple primer-pool loci stitched together.

    Returns ``merged_record`` unchanged when the recompute cannot proceed
    (no qctx table, no larger peer, etc.) — caller is responsible for
    setting cer_factor=None upstream if a stale inherited value is
    unacceptable.

    Interpretation:
    - The recomputed factor compares the merged candidate's shipped
      consensus (column-vote of contributor consensuses) to each larger
      peer's shipped consensus. K and the per-position q_ctx values are
      re-derived from the merged consensus, so the factor reflects the
      merged object's variant signature against post-merge peers.
    - For same-primer merges within a core identity group, this is the
      semantically correct CER call: the noise model was originally
      designed for exactly this comparison, just with a re-derived
      candidate state.
    - For cross-primer merges (None), summarize's routing treats the value
      identically to other "no valid peer comparison" cases — record always
      passes the cer_factor filter.
    """
    contributor_primer_sets = {frozenset(c.primers or []) for c in contributors}
    is_cross_primer = len(contributor_primer_sets) > 1
    if is_cross_primer:
        return merged_record._replace(cer_factor=None)

    if qctx_table is None:
        return merged_record

    larger_peers = [
        p for p in bucket_members
        if p is not merged_record and p.size > merged_record.size
    ]
    if not larger_peers:
        # Anchor of its bucket (or only peer) — no comparison available.
        return merged_record._replace(cer_factor=None)

    group_N = sum(m.size for m in bucket_members)

    best_factor: Optional[float] = None
    for peer in larger_peers:
        tags = classify_pairwise_differences(
            merged_record.sequence, peer.sequence,
            hp_min_length=hp_min_length,
        )
        if not tags:
            continue
        qctx_values: List[float] = []
        for tag in tags:
            q = get_qctx(tag, table=qctx_table)
            if q is None:
                continue
            qctx_values.append(q)
        if not qctx_values:
            continue
        longer = (
            peer.sequence
            if len(peer.sequence) >= len(merged_record.sequence)
            else merged_record.sequence
        )
        n_sites = _count_independent_sites(longer)
        factor = compute_cer_factor(
            N=group_N,
            M=merged_record.size,
            n_sites=n_sites,
            q_ctx_per_position=qctx_values,
            alpha=alpha,
        )
        if factor is None:
            continue
        if best_factor is None or factor < best_factor:
            best_factor = factor

    return merged_record._replace(cer_factor=best_factor)
