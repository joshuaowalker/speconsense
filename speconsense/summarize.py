#!/usr/bin/env python3

import os
import re
import glob
import csv
import shutil
import argparse
import itertools
import logging
import statistics
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional, NamedTuple
from collections import defaultdict

import edlib
import numpy as np
from Bio import SeqIO
from adjusted_identity import score_alignment, AdjustmentParams, ScoringFormat


class ConsensusInfo(NamedTuple):
    """Information about a consensus sequence from speconsense output."""
    sample_name: str
    cluster_id: str
    sequence: str
    length: int
    ric: int
    size: int
    file_path: str


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Speconsense output with advanced variant handling.")
    parser.add_argument("--min-ric", type=int, default=3,
                        help="Minimum Reads in Consensus (RiC) threshold (default: 3)")
    parser.add_argument("--source", type=str, default=".",
                        help="Source directory containing Speconsense output (default: current directory)")
    parser.add_argument("--summary-dir", type=str, default="__Summary__",
                        help="Output directory for summary files (default: __Summary__)")
    parser.add_argument("--variant-merge-threshold", type=int, default=2,
                        help="Maximum number of SNPs allowed in merged variants (default: 2)")
    parser.add_argument("--variant-group-identity", type=float, default=0.9,
                        help="Identity threshold for variant grouping using HAC (default: 0.9)")
    parser.add_argument("--max-variants", type=int, default=2,
                        help="Maximum number of additional variants to output per group (default: 2)")
    parser.add_argument("--variant-selection", choices=["size", "diversity"], default="size",
                        help="Variant selection strategy: size or diversity (default: size)")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level")
    return parser.parse_args()


def setup_logging(log_level: str):
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def parse_consensus_header(header: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """
    Extract information from Speconsense consensus FASTA header.
    
    Note: Stability metrics (median_diff, p95_diff) are intentionally ignored
    as they will be dropped from final output - they're debug info only.
    """
    sample_match = re.match(r'>([^;]+);(.+)', header)
    if not sample_match:
        return None, None, None

    sample_name = sample_match.group(1)
    info_string = sample_match.group(2)

    # Extract RiC value
    ric_match = re.search(r'ric=(\d+)', info_string)
    ric = int(ric_match.group(1)) if ric_match else 0

    # Extract size value
    size_match = re.search(r'size=(\d+)', info_string)
    size = int(size_match.group(1)) if size_match else 0

    # Note: median_diff and p95_diff are available in original files but
    # are intentionally not extracted here as they will be dropped from final output
    
    return sample_name, ric, size


def load_consensus_sequences(source_folder: str, min_ric: int) -> List[ConsensusInfo]:
    """Load all consensus sequences from speconsense output files."""
    consensus_list = []
    
    # Find all consensus FASTA files matching the new naming pattern
    fasta_pattern = os.path.join(source_folder, "*_consensus_sequences.fasta")
    
    for fasta_file in sorted(glob.glob(fasta_pattern)):
        logging.info(f"Processing consensus file: {fasta_file}")
        
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                sample_name, ric, size = parse_consensus_header(f">{record.description}")
                
                if sample_name and ric >= min_ric:
                    # Extract cluster ID from sample name (e.g., "sample-c1" -> "c1")
                    cluster_match = re.search(r'-c(\d+)$', sample_name)
                    cluster_id = cluster_match.group(0) if cluster_match else sample_name
                    
                    consensus_info = ConsensusInfo(
                        sample_name=sample_name,
                        cluster_id=cluster_id,
                        sequence=str(record.seq),
                        length=len(record.seq),
                        ric=ric,
                        size=size,
                        file_path=fasta_file
                    )
                    consensus_list.append(consensus_info)
    
    logging.info(f"Loaded {len(consensus_list)} consensus sequences")
    return consensus_list


def calculate_substitution_distance(seq1: str, seq2: str) -> int:
    """
    Calculate distance between two sequences counting only substitutions.
    Uses adjusted identity with custom parameters to handle homopolymers
    but count only substitutions for the distance.
    """
    if not seq1 or not seq2:
        return max(len(seq1), len(seq2))
        
    # Get alignment from edlib
    result = edlib.align(seq1, seq2, task="path")
    if result["editDistance"] == -1:
        return max(len(seq1), len(seq2))
        
    # Get nice alignment for adjusted identity scoring
    alignment = edlib.getNiceAlignment(result, seq1, seq2)
    if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
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
    
    # Count only substitutions (not homopolymer adjustments or indels)
    substitutions = score_result.score_aligned.count('X')
    
    logging.debug(f"Substitution distance: {substitutions} substitutions")
    return substitutions


def merge_variants_per_file(consensus_list: List[ConsensusInfo], 
                           variant_merge_threshold: int) -> Tuple[List[ConsensusInfo], Dict[str, List[str]]]:
    """
    Merge consensus sequences that differ by at most variant_merge_threshold substitutions.
    Only merges variants within the same input file to prevent cross-specimen contamination.
    Returns merged consensus list and traceability information.
    """
    if variant_merge_threshold <= 0:
        logging.info("Variant merging disabled")
        return consensus_list, {}
    
    logging.info(f"Merging variants with threshold of {variant_merge_threshold} substitutions")
    
    # Group by input file to ensure we only merge within specimens
    file_groups = defaultdict(list)
    for cons in consensus_list:
        file_groups[cons.file_path].append(cons)
    
    merged_consensus = []
    merge_traceability = {}
    
    for file_path, file_consensuses in file_groups.items():
        if len(file_consensuses) <= 1:
            # No variants to merge in this file
            merged_consensus.extend(file_consensuses)
            continue
        
        # Extract file base name for logging
        file_name = os.path.basename(file_path)
        logging.info(f"Processing {len(file_consensuses)} variants from file {file_name}")
        
        # Sort by cluster size (largest first) for merging priority
        file_consensuses.sort(key=lambda x: x.size, reverse=True)
        
        # Build distance matrix
        n = len(file_consensuses)
        distances = np.zeros((n, n))
        
        for i, j in itertools.combinations(range(n), 2):
            dist = calculate_substitution_distance(
                file_consensuses[i].sequence,
                file_consensuses[j].sequence
            )
            distances[i, j] = dist
            distances[j, i] = dist
        
        # Greedy merging starting with largest cluster
        merged_indices = set()
        merge_groups = []
        
        for i in range(n):
            if i in merged_indices:
                continue
                
            # Start a new merge group with the current consensus
            current_group = [i]
            merged_indices.add(i)
            
            # Find all compatible consensuses
            for j in range(i + 1, n):
                if j in merged_indices:
                    continue
                    
                # Check if j is compatible with all members of current group
                compatible = True
                for group_member in current_group:
                    if distances[group_member, j] > variant_merge_threshold:
                        compatible = False
                        break
                
                if compatible:
                    current_group.append(j)
                    merged_indices.add(j)
            
            merge_groups.append(current_group)
        
        # Create merged consensus sequences
        for group_idx, group in enumerate(merge_groups):
            if len(group) == 1:
                # No merging needed
                merged_consensus.append(file_consensuses[group[0]])
            else:
                # Merge group - use the largest (first) as representative
                representative = file_consensuses[group[0]]
                
                # Create traceability entry
                original_names = [file_consensuses[idx].sample_name for idx in group]
                merge_key = representative.sample_name
                merge_traceability[merge_key] = original_names
                
                # Create new merged consensus info
                total_size = sum(file_consensuses[idx].size for idx in group)
                total_ric = sum(file_consensuses[idx].ric for idx in group)
                
                merged_info = ConsensusInfo(
                    sample_name=representative.sample_name,
                    cluster_id=representative.cluster_id,
                    sequence=representative.sequence,
                    length=representative.length,
                    ric=total_ric,
                    size=total_size,
                    file_path=representative.file_path
                )
                
                merged_consensus.append(merged_info)
                
                logging.info(f"Merged {len(group)} variants from {file_name}: {original_names}")
    
    logging.info(f"Variant merging completed: {len(consensus_list)} -> {len(merged_consensus)} sequences")
    return merged_consensus, merge_traceability


def calculate_adjusted_identity_distance(seq1: str, seq2: str) -> float:
    """Calculate adjusted identity distance between two sequences."""
    if not seq1 or not seq2:
        return 1.0  # Maximum distance
    
    if seq1 == seq2:
        return 0.0
        
    # Get alignment from edlib
    result = edlib.align(seq1, seq2, task="path")
    if result["editDistance"] == -1:
        return 1.0
        
    # Get nice alignment for adjusted identity scoring
    alignment = edlib.getNiceAlignment(result, seq1, seq2)
    if not alignment or not alignment.get('query_aligned') or not alignment.get('target_aligned'):
        return 1.0
        
    # Use standard adjusted identity parameters
    score_result = score_alignment(
        alignment['query_aligned'], 
        alignment['target_aligned']
    )
    
    # Convert adjusted identity to distance
    return 1.0 - score_result.identity


def perform_hac_clustering(consensus_list: List[ConsensusInfo], 
                          variant_group_identity: float) -> Dict[int, List[ConsensusInfo]]:
    """
    Perform Hierarchical Agglomerative Clustering using complete linkage.
    Separates specimens from variants based on identity threshold.
    Returns groups of consensus sequences.
    """
    if len(consensus_list) <= 1:
        return {0: consensus_list}
    
    logging.info(f"Performing HAC clustering with {variant_group_identity} identity threshold")
    
    n = len(consensus_list)
    distance_threshold = 1.0 - variant_group_identity
    
    # Initialize each sequence as its own cluster
    clusters = [[i] for i in range(n)]
    
    # Build initial distance matrix between individual sequences
    seq_distances = {}
    for i, j in itertools.combinations(range(n), 2):
        dist = calculate_adjusted_identity_distance(
            consensus_list[i].sequence,
            consensus_list[j].sequence
        )
        seq_distances[(i, j)] = dist
        seq_distances[(j, i)] = dist
    
    def cluster_distance(cluster1: List[int], cluster2: List[int]) -> float:
        """Calculate complete linkage distance between two clusters."""
        max_dist = 0.0
        for i in cluster1:
            for j in cluster2:
                if i < j:
                    dist = seq_distances[(i, j)]
                elif i > j:
                    dist = seq_distances[(j, i)]
                else:
                    dist = 0.0
                max_dist = max(max_dist, dist)
        return max_dist
    
    # Perform hierarchical clustering
    while len(clusters) > 1:
        # Find closest pair of clusters
        min_distance = float('inf')
        merge_pair = None
        
        for i, j in itertools.combinations(range(len(clusters)), 2):
            dist = cluster_distance(clusters[i], clusters[j])
            if dist < min_distance:
                min_distance = dist
                merge_pair = (i, j)
        
        # If minimum distance exceeds threshold, stop clustering
        if min_distance >= distance_threshold:
            break
            
        # Merge the closest clusters
        i, j = merge_pair
        merged_cluster = clusters[i] + clusters[j]
        
        # Remove old clusters (in reverse order to maintain indices)
        new_clusters = []
        for idx, cluster in enumerate(clusters):
            if idx != i and idx != j:
                new_clusters.append(cluster)
        new_clusters.append(merged_cluster)
        clusters = new_clusters
        
        logging.debug(f"Merged clusters with distance {min_distance:.3f}, now have {len(clusters)} clusters")
    
    # Convert clusters to groups of ConsensusInfo
    groups = {}
    for group_id, cluster_indices in enumerate(clusters):
        group_members = [consensus_list[idx] for idx in cluster_indices]
        groups[group_id] = group_members
    
    logging.info(f"HAC clustering created {len(groups)} groups")
    for group_id, group_members in groups.items():
        member_names = [m.sample_name for m in group_members]
        logging.info(f"Group {group_id}: {member_names}")
    
    return groups


def select_variants(group: List[ConsensusInfo], 
                   max_variants: int,
                   variant_selection: str) -> List[ConsensusInfo]:
    """
    Select variants from a group based on the specified strategy.
    Always includes the largest variant first.
    """
    if len(group) <= max_variants + 1:  # +1 for main variant
        return sorted(group, key=lambda x: x.size, reverse=True)
    
    # Sort by size, largest first
    sorted_group = sorted(group, key=lambda x: x.size, reverse=True)
    
    # Always include the largest (main) variant
    selected = [sorted_group[0]]
    candidates = sorted_group[1:]
    
    if variant_selection == "size":
        # Select by size
        selected.extend(candidates[:max_variants])
    else:  # diversity
        # Select by diversity (maximum distance from already selected)
        while len(selected) < max_variants + 1 and candidates:
            best_candidate = None
            best_min_distance = -1
            
            for candidate in candidates:
                # Calculate minimum distance to all selected variants
                min_distance = min(
                    calculate_adjusted_identity_distance(candidate.sequence, sel.sequence)
                    for sel in selected
                )
                
                if min_distance > best_min_distance:
                    best_min_distance = min_distance
                    best_candidate = candidate
            
            if best_candidate:
                selected.append(best_candidate)
                candidates.remove(best_candidate)
    
    return selected


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
    os.makedirs(os.path.join(summary_folder, 'original_clusters'), exist_ok=True)
    
    final_consensus = []
    naming_info = {}
    
    # Sort groups by size of largest member (descending)
    sorted_groups = sorted(groups.items(), 
                          key=lambda x: max(m.size for m in x[1]), 
                          reverse=True)
    
    for group_idx, (_, group_members) in enumerate(sorted_groups, 1):
        # Select variants for this group
        selected_variants = select_variants(group_members, max_variants, variant_selection)
        
        # Create naming for this group
        group_naming = []
        
        for variant_idx, variant in enumerate(selected_variants):
            if variant_idx == 0:
                # Main variant gets simple numeric suffix
                new_name = f"{variant.sample_name.split('-c')[0]}-{group_idx}"
            else:
                # Additional variants get .v suffix
                new_name = f"{variant.sample_name.split('-c')[0]}-{group_idx}.v{variant_idx}"
            
            # Create new ConsensusInfo with updated name
            renamed_variant = ConsensusInfo(
                sample_name=new_name,
                cluster_id=variant.cluster_id,
                sequence=variant.sequence,
                length=variant.length,
                ric=variant.ric,
                size=variant.size,
                file_path=variant.file_path
            )
            
            final_consensus.append(renamed_variant)
            group_naming.append((variant.sample_name, new_name))
        
        naming_info[group_idx] = group_naming
    
    return final_consensus, naming_info


def write_output_files(final_consensus: List[ConsensusInfo],
                      merge_traceability: Dict[str, List[str]],
                      naming_info: Dict,
                      summary_folder: str):
    """Write all output files."""
    
    # Write individual FASTA files (without stability metrics)
    for consensus in final_consensus:
        output_file = os.path.join(summary_folder, f"{consensus.sample_name}.fasta")
        with open(output_file, 'w') as f:
            # Clean header without stability metrics (p95_diff, median_diff)
            f.write(f">{consensus.sample_name};length={consensus.length} ric={consensus.ric} size={consensus.size}\n")
            f.write(f"{consensus.sequence}\n")
    
    # Write combined summary.fasta (without stability metrics)
    summary_fasta_path = os.path.join(summary_folder, 'summary.fasta')
    with open(summary_fasta_path, 'w') as f:
        for consensus in final_consensus:
            # Clean header without stability metrics (p95_diff, median_diff)
            f.write(f">{consensus.sample_name};length={consensus.length} ric={consensus.ric} size={consensus.size}\n")
            f.write(f"{consensus.sequence}\n")
    
    # Write summary statistics
    summary_txt_path = os.path.join(summary_folder, 'summary.txt')
    with open(summary_txt_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['Filename', 'Length', 'Reads in Consensus', 'Multiple'])
        
        unique_samples = set()
        total_ric = 0
        
        for consensus in final_consensus:
            writer.writerow([consensus.sample_name, consensus.length, consensus.ric, 1])
            base_name = consensus.sample_name.split('-')[0]
            unique_samples.add(base_name)
            total_ric += consensus.ric
        
        writer.writerow([])
        writer.writerow(['Total Unique Samples', len(unique_samples)])
        writer.writerow(['Total Consensus Sequences', len(final_consensus)])
        writer.writerow(['Total Reads in Consensus Sequences', total_ric])
    
    # Write variants_merged.txt for traceability
    variants_merged_path = os.path.join(summary_folder, 'variants_merged.txt')
    with open(variants_merged_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['Final_Name', 'Original_Clusters', 'Merged_From'])
        
        for consensus in final_consensus:
            original_clusters = []
            merged_from = []
            
            # Check if this was renamed
            for group_id, renamings in naming_info.items():
                for original_name, final_name in renamings:
                    if final_name == consensus.sample_name:
                        original_clusters.append(original_name)
                        
                        # Check if this was also merged
                        if original_name in merge_traceability:
                            merged_from.extend(merge_traceability[original_name])
                        else:
                            merged_from.append(original_name)
                        break
            
            writer.writerow([
                consensus.sample_name,
                ';'.join(original_clusters) if original_clusters else consensus.sample_name,
                ';'.join(merged_from) if merged_from else consensus.sample_name
            ])


def process_single_specimen(file_consensuses: List[ConsensusInfo], 
                           args) -> Tuple[List[ConsensusInfo], Dict[str, List[str]], Dict]:
    """
    Process a single specimen file: merge variants, perform HAC clustering, and select final variants.
    Returns final consensus list, merge traceability, and naming info for this specimen.
    """
    if not file_consensuses:
        return [], {}, {}
    
    file_name = os.path.basename(file_consensuses[0].file_path)
    logging.info(f"Processing specimen from file: {file_name}")
    
    # Phase 1: Merge variants within this specimen file
    merged_consensus, merge_traceability = merge_variants_within_specimen(file_consensuses, args.variant_merge_threshold)
    
    if not merged_consensus:
        return [], merge_traceability, {}
    
    # Phase 2: HAC clustering to separate organisms within this specimen (primary vs contaminants)
    organism_groups = perform_hac_clustering(merged_consensus, args.variant_group_identity)
    
    # Phase 3: Select representative variants for each organism in this specimen
    final_consensus = []
    naming_info = {}
    
    # Sort organism groups by size of largest member (descending)
    sorted_groups = sorted(organism_groups.items(), 
                          key=lambda x: max(m.size for m in x[1]), 
                          reverse=True)
    
    for organism_idx, (_, group_members) in enumerate(sorted_groups):
        # Select variants for this organism
        selected_variants = select_variants(group_members, args.max_variants, args.variant_selection)
        
        # Create naming for this organism within this specimen
        group_naming = []
        
        for variant_idx, variant in enumerate(selected_variants):
            if variant_idx == 0:
                # Main variant gets local organism number (1, 2, etc. within this specimen)
                new_name = f"{variant.sample_name.split('-c')[0]}-{organism_idx + 1}"
            else:
                # Additional variants get .v suffix  
                new_name = f"{variant.sample_name.split('-c')[0]}-{organism_idx + 1}.v{variant_idx}"
            
            # Create new ConsensusInfo with updated name
            renamed_variant = ConsensusInfo(
                sample_name=new_name,
                cluster_id=variant.cluster_id,
                sequence=variant.sequence,
                length=variant.length,
                ric=variant.ric,
                size=variant.size,
                file_path=variant.file_path
            )
            
            final_consensus.append(renamed_variant)
            group_naming.append((variant.sample_name, new_name))
        
        naming_info[organism_idx + 1] = group_naming
    
    logging.info(f"Processed {file_name}: {len(organism_groups)} organisms -> {len(final_consensus)} final variants")
    
    return final_consensus, merge_traceability, naming_info


def merge_variants_within_specimen(file_consensuses: List[ConsensusInfo], 
                                 variant_merge_threshold: int) -> Tuple[List[ConsensusInfo], Dict[str, List[str]]]:
    """
    Merge consensus sequences within a single specimen file that differ by at most variant_merge_threshold substitutions.
    Returns merged consensus list and traceability information.
    """
    if variant_merge_threshold <= 0 or len(file_consensuses) <= 1:
        if variant_merge_threshold <= 0:
            logging.debug("Variant merging disabled for this specimen")
        return file_consensuses, {}
    
    file_name = os.path.basename(file_consensuses[0].file_path)
    logging.info(f"Merging {len(file_consensuses)} variants within {file_name} (threshold: {variant_merge_threshold} substitutions)")
    
    # Sort by cluster size (largest first) for merging priority
    file_consensuses.sort(key=lambda x: x.size, reverse=True)
    
    # Build distance matrix
    n = len(file_consensuses)
    distances = np.zeros((n, n))
    
    for i, j in itertools.combinations(range(n), 2):
        dist = calculate_substitution_distance(
            file_consensuses[i].sequence,
            file_consensuses[j].sequence
        )
        distances[i, j] = dist
        distances[j, i] = dist
    
    # Greedy merging starting with largest cluster
    merged_indices = set()
    merge_groups = []
    
    for i in range(n):
        if i in merged_indices:
            continue
            
        # Start a new merge group with the current consensus
        current_group = [i]
        merged_indices.add(i)
        
        # Find all compatible consensuses
        for j in range(i + 1, n):
            if j in merged_indices:
                continue
                
            # Check if j is compatible with all members of current group
            compatible = True
            for group_member in current_group:
                if distances[group_member, j] > variant_merge_threshold:
                    compatible = False
                    break
            
            if compatible:
                current_group.append(j)
                merged_indices.add(j)
        
        merge_groups.append(current_group)
    
    # Create merged consensus sequences
    merged_consensus = []
    merge_traceability = {}
    
    for group in merge_groups:
        if len(group) == 1:
            # No merging needed
            merged_consensus.append(file_consensuses[group[0]])
        else:
            # Merge group - use the largest (first) as representative
            representative = file_consensuses[group[0]]
            
            # Create traceability entry
            original_names = [file_consensuses[idx].sample_name for idx in group]
            merge_key = representative.sample_name
            merge_traceability[merge_key] = original_names
            
            # Create new merged consensus info
            total_size = sum(file_consensuses[idx].size for idx in group)
            total_ric = sum(file_consensuses[idx].ric for idx in group)
            
            merged_info = ConsensusInfo(
                sample_name=representative.sample_name,
                cluster_id=representative.cluster_id,
                sequence=representative.sequence,
                length=representative.length,
                ric=total_ric,
                size=total_size,
                file_path=representative.file_path
            )
            
            merged_consensus.append(merged_info)
            
            logging.info(f"Merged {len(group)} variants in {file_name}: {original_names}")
    
    return merged_consensus, merge_traceability


def main():
    """Main function to process command line arguments and run the summarization."""
    args = parse_arguments()
    setup_logging(args.log_level)
    
    logging.info("Starting enhanced speconsense summarization")
    logging.info("Note: Stability metrics (median_diff, p95_diff) will be dropped from final output")
    logging.info("Original stability metrics are preserved in cluster_debug/ files for debugging")
    logging.info("Processing each specimen file independently to organize variants within specimens")
    
    # Load all consensus sequences
    consensus_list = load_consensus_sequences(args.source, args.min_ric)
    if not consensus_list:
        logging.error("No consensus sequences found")
        return
    
    # Group consensus sequences by input file (one file per specimen)
    file_groups = defaultdict(list)
    for cons in consensus_list:
        file_groups[cons.file_path].append(cons)
    
    # Process each specimen file independently
    all_final_consensus = []
    all_merge_traceability = {}
    all_naming_info = {}
    
    for file_path in sorted(file_groups.keys()):
        file_consensuses = file_groups[file_path]
        
        final_consensus, merge_traceability, naming_info = process_single_specimen(
            file_consensuses, args
        )
        
        all_final_consensus.extend(final_consensus)
        all_merge_traceability.update(merge_traceability)
        
        # Update naming info with unique keys per specimen
        file_name = os.path.basename(file_path)
        for group_id, group_naming in naming_info.items():
            unique_key = f"{file_name}_{group_id}"
            all_naming_info[unique_key] = group_naming
    
    # Create output directories
    os.makedirs(args.summary_dir, exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'FASTQ Files'), exist_ok=True)
    os.makedirs(os.path.join(args.summary_dir, 'original_clusters'), exist_ok=True)
    
    # Write output files
    write_output_files(all_final_consensus, all_merge_traceability, all_naming_info, args.summary_dir)
    
    logging.info(f"Enhanced summarization completed successfully")
    logging.info(f"Final output: {len(all_final_consensus)} consensus sequences in {args.summary_dir}")
    
    # Copy original consensus files for traceability and debugging
    # These files contain the original stability metrics (median_diff, p95_diff)
    original_clusters_dir = os.path.join(args.summary_dir, 'original_clusters')
    copied_files = set()  # Track unique files to avoid duplicates
    
    for consensus in consensus_list:
        if os.path.exists(consensus.file_path) and consensus.file_path not in copied_files:
            shutil.copy(consensus.file_path, original_clusters_dir)
            copied_files.add(consensus.file_path)
            logging.debug(f"Copied original consensus file: {os.path.basename(consensus.file_path)}")
    
    logging.info(f"Copied {len(copied_files)} original consensus files with stability metrics to original_clusters/")


if __name__ == "__main__":
    main()