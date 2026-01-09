"""MSA-based variant merging for speconsense-summarize.

Provides functions for finding and merging compatible variants within HAC groups
using exhaustive subset evaluation with SPOA multiple sequence alignment.
"""

import itertools
import logging
from typing import List, Tuple, Dict
from collections import defaultdict

from speconsense.types import ConsensusInfo, OverlapMergeInfo

from .iupac import merge_bases_to_iupac, primers_are_same
from .analysis import (
    run_spoa_msa,
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

    # Create merged ConsensusInfo
    consensus_sequence = ''.join(consensus_seq)
    total_size = sum(v.size for v in variants)
    total_ric = sum(v.ric for v in variants)

    # Collect RiC values, preserving any prior merge history
    raw_ric_values = []
    for v in variants:
        if v.raw_ric:
            raw_ric_values.extend(v.raw_ric)  # Flatten prior merge history
        else:
            raw_ric_values.append(v.ric)
    raw_ric_values = sorted(raw_ric_values, reverse=True) if len(variants) > 1 else None

    # Collect lengths, preserving any prior merge history
    raw_len_values = []
    for v in variants:
        if v.raw_len:
            raw_len_values.extend(v.raw_len)  # Flatten prior merge history
        else:
            raw_len_values.append(len(v.sequence))
    raw_len_values = sorted(raw_len_values, reverse=True) if len(variants) > 1 else None

    # Use name from largest variant
    largest_variant = max(variants, key=lambda v: v.size)

    return ConsensusInfo(
        sample_name=largest_variant.sample_name,
        cluster_id=largest_variant.cluster_id,
        sequence=consensus_sequence,
        ric=total_ric,
        size=total_size,
        file_path=largest_variant.file_path,
        snp_count=snp_count if snp_count > 0 else None,
        primers=largest_variant.primers,
        raw_ric=raw_ric_values,
        raw_len=raw_len_values,
        rid=largest_variant.rid,  # Preserve identity metrics from largest variant
        rid_min=largest_variant.rid_min,
    )


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

    # Create merged ConsensusInfo
    consensus_sequence = ''.join(consensus_seq)
    total_size = sum(v.size for v in variants)
    total_ric = sum(v.ric for v in variants)

    # Collect RiC values, preserving any prior merge history
    raw_ric_values = []
    for v in variants:
        if v.raw_ric:
            raw_ric_values.extend(v.raw_ric)  # Flatten prior merge history
        else:
            raw_ric_values.append(v.ric)
    raw_ric_values = sorted(raw_ric_values, reverse=True) if len(variants) > 1 else None

    # Collect lengths, preserving any prior merge history
    raw_len_values = []
    for v in variants:
        if v.raw_len:
            raw_len_values.extend(v.raw_len)  # Flatten prior merge history
        else:
            raw_len_values.append(len(v.sequence))
    raw_len_values = sorted(raw_len_values, reverse=True) if len(variants) > 1 else None

    # Use name from largest variant
    largest_variant = max(variants, key=lambda v: v.size)

    return ConsensusInfo(
        sample_name=largest_variant.sample_name,
        cluster_id=largest_variant.cluster_id,
        sequence=consensus_sequence,
        ric=total_ric,
        size=total_size,
        file_path=largest_variant.file_path,
        snp_count=snp_count if snp_count > 0 else None,
        primers=largest_variant.primers,
        raw_ric=raw_ric_values,
        raw_len=raw_len_values,
        rid=largest_variant.rid,
        rid_min=largest_variant.rid_min,
    )


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
                if use_overlap_mode:
                    # Use overlap-aware analysis for primer pool scenarios
                    original_lengths = [len(v.sequence) for v in subset_variants]
                    variant_stats = analyze_msa_columns_overlap_aware(
                        subset_aligned, args.min_merge_overlap, original_lengths
                    )

                    # Check overlap requirement
                    shorter_len = min(original_lengths)
                    effective_threshold = min(args.min_merge_overlap, shorter_len)
                    if variant_stats['overlap_bp'] < effective_threshold:
                        # Insufficient overlap - skip this subset
                        continue
                else:
                    # Use standard analysis
                    variant_stats = analyze_msa_columns(subset_aligned)

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
                        # Extract specimen name (remove cluster suffix like -c1)
                        specimen = merged_consensus.sample_name.rsplit('-c', 1)[0] if '-c' in merged_consensus.sample_name else merged_consensus.sample_name
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
