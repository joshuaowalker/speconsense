"""
Shared type definitions for speconsense.

This module contains data classes used across multiple modules,
extracted to avoid circular imports.
"""

from typing import List, Optional, NamedTuple


class ConsensusInfo(NamedTuple):
    """Information about a consensus sequence from speconsense output."""
    sample_name: str
    cluster_id: str
    sequence: str
    ric: int
    size: int
    file_path: str
    snp_count: Optional[int] = None  # Number of SNPs from IUPAC consensus generation
    primers: Optional[List[str]] = None  # List of detected primer names
    raw_ric: Optional[List[int]] = None  # RiC values of .raw source variants
    raw_len: Optional[List[int]] = None  # Lengths of merged source sequences
    rid: Optional[float] = None  # Mean read identity (internal consistency metric)
    rid_min: Optional[float] = None  # Minimum read identity (worst-case read)
    merge_indel_count: Optional[int] = None  # Number of indels consumed by merging (for cumulative tracking)


class OverlapMergeInfo(NamedTuple):
    """Information about a single overlap merge event for quality reporting."""
    specimen: str           # Specimen name
    iteration: int          # Merge iteration (1 = first pass, 2+ = iterative)
    input_clusters: List[str]   # Cluster IDs involved in merge
    input_lengths: List[int]    # Original sequence lengths
    input_rics: List[int]       # RiC values of input sequences
    overlap_bp: int             # Overlap region size in bp
    prefix_bp: int              # Extension before overlap
    suffix_bp: int              # Extension after overlap
    output_length: int          # Final merged sequence length
