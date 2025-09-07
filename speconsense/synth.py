#!/usr/bin/env python3
"""
Speconsense-synth: Synthetic read generator for testing consensus algorithms.

Generates simulated reads from input sequences with controlled error rates
for testing clustering and consensus generation behavior.
"""

import argparse
import sys
import random
import logging
from typing import List, Tuple, Dict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_ratios(ratio_str: str, num_sequences: int) -> List[float]:
    """Parse ratio string and validate against number of sequences.
    
    Args:
        ratio_str: Comma-separated ratios (e.g., "50,30,20")
        num_sequences: Number of input sequences
        
    Returns:
        List of normalized ratios summing to 1.0
    """
    if not ratio_str:
        # Equal distribution by default
        return [1.0 / num_sequences] * num_sequences
    
    ratios = [float(r.strip()) for r in ratio_str.split(',')]
    
    if len(ratios) != num_sequences:
        raise ValueError(f"Number of ratios ({len(ratios)}) must match number of sequences ({num_sequences})")
    
    if any(r < 0 for r in ratios):
        raise ValueError("Ratios must be non-negative")
    
    total = sum(ratios)
    if total == 0:
        raise ValueError("At least one ratio must be positive")
    
    # Normalize to sum to 1.0
    return [r / total for r in ratios]


def error_rate_to_phred(error_rate: float) -> int:
    """Convert error rate to Phred quality score.
    
    Args:
        error_rate: Probability of error (0-1)
        
    Returns:
        Phred quality score (0-40, capped)
    """
    if error_rate <= 0:
        return 40  # Cap at Q40
    if error_rate >= 1:
        return 0
    
    import math
    phred = -10 * math.log10(error_rate)
    return min(40, max(0, int(round(phred))))


def introduce_errors(sequence: str, error_rate: float, rng: random.Random) -> str:
    """Introduce errors into a sequence at specified rate.
    
    Each position has error_rate chance of mutation.
    Error types (insertion, deletion, substitution) are equally likely.
    
    Args:
        sequence: Original sequence
        error_rate: Probability of error at each position
        rng: Random number generator for reproducibility
        
    Returns:
        Mutated sequence
    """
    if error_rate <= 0:
        return sequence
    
    result = []
    bases = ['A', 'C', 'G', 'T']
    
    for base in sequence:
        if rng.random() < error_rate:
            # Error occurs - choose type
            error_type = rng.choice(['insertion', 'deletion', 'substitution'])
            
            if error_type == 'deletion':
                # Skip this base
                continue
            elif error_type == 'insertion':
                # Add current base plus a random insertion
                result.append(base)
                result.append(rng.choice(bases))
            else:  # substitution
                # Replace with different base
                alternatives = [b for b in bases if b != base.upper()]
                result.append(rng.choice(alternatives))
        else:
            result.append(base)
    
    return ''.join(result)


def normalize_sequence(seq_str: str, seq_id: str) -> str:
    """Normalize a sequence by removing whitespace and converting to uppercase.
    
    Args:
        seq_str: Input sequence string
        seq_id: Sequence ID for error messages
        
    Returns:
        Normalized sequence (uppercase, no whitespace)
        
    Raises:
        ValueError: If sequence contains non-ACGT bases
    """
    # Remove all whitespace and convert to uppercase
    normalized = ''.join(seq_str.split()).upper()
    
    # Check for non-ACGT bases
    valid_bases = set('ACGT')
    invalid_bases = set(normalized) - valid_bases
    
    if invalid_bases:
        logging.warning(f"Sequence '{seq_id}' contains non-ACGT bases: {sorted(invalid_bases)}")
        logging.warning(f"These bases will be treated as-is but may cause unexpected behavior")
    
    return normalized


def generate_reads(sequences: List[SeqRecord], 
                  num_reads: int,
                  error_rate: float,
                  ratios: List[float],
                  seed: int = None) -> List[SeqRecord]:
    """Generate synthetic reads from input sequences.
    
    Args:
        sequences: Input sequences to generate reads from
        num_reads: Total number of reads to generate
        error_rate: Per-base error probability
        ratios: Relative abundance of each sequence
        seed: Random seed for reproducibility
        
    Returns:
        List of synthetic reads as SeqRecord objects
    """
    rng = random.Random(seed)
    reads = []
    
    # Calculate reads per sequence based on ratios
    reads_per_seq = []
    cumulative = 0
    for i, ratio in enumerate(ratios):
        if i == len(ratios) - 1:
            # Last sequence gets remaining reads to ensure exact count
            reads_per_seq.append(num_reads - cumulative)
        else:
            count = int(round(num_reads * ratio))
            reads_per_seq.append(count)
            cumulative += count
    
    # Generate reads for each sequence
    quality = error_rate_to_phred(error_rate)
    
    for seq_idx, (seq_record, seq_reads) in enumerate(zip(sequences, reads_per_seq)):
        # Normalize sequence
        sequence_str = normalize_sequence(str(seq_record.seq), seq_record.id)
        
        for read_idx in range(seq_reads):
            # Introduce errors
            mutated_seq = introduce_errors(sequence_str, error_rate, rng)
            
            # Create read ID with provenance (unique ID + source)
            read_id = f"read_{len(reads) + 1}_from_{seq_record.id}"
            
            # Create SeqRecord with quality scores
            read_record = SeqRecord(
                Seq(mutated_seq),
                id=read_id,
                description=f"source={seq_record.id} error_rate={error_rate:.3f}",
                letter_annotations={'phred_quality': [quality] * len(mutated_seq)}
            )
            
            reads.append(read_record)
    
    # Shuffle reads to mix sequences
    rng.shuffle(reads)
    
    return reads


def main():
    """Main entry point for speconsense-synth."""
    parser = argparse.ArgumentParser(
        description='Generate synthetic reads from reference sequences with controlled error rates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate 1000 reads with 10% error rate
  speconsense-synth reference.fasta -n 1000 -e 0.1 -o synthetic_reads.fastq
  
  # Generate reads from multiple sequences with specific ratios
  speconsense-synth variants.fasta -n 5000 -e 0.15 --ratios 70,30 -o mixed_reads.fastq
  
  # Set seed for reproducible results
  speconsense-synth reference.fasta -n 1000 -e 0.1 --seed 42 -o synthetic_reads.fastq
        """
    )
    
    parser.add_argument('input', help='Input FASTA file with reference sequence(s)')
    parser.add_argument('-n', '--num-reads', type=int, default=1000,
                       help='Number of reads to generate (default: 1000)')
    parser.add_argument('-e', '--error-rate', type=float, default=0.1,
                       help='Per-base error rate (default: 0.1)')
    parser.add_argument('-o', '--output', default='synthetic_reads.fastq',
                       help='Output FASTQ file (default: synthetic_reads.fastq)')
    parser.add_argument('--ratios', type=str,
                       help='Comma-separated ratios for multiple sequences (e.g., 70,30)')
    parser.add_argument('--seed', type=int,
                       help='Random seed for reproducibility')
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level (default: INFO)')
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Validate error rate
    if not 0 <= args.error_rate <= 1:
        parser.error("Error rate must be between 0 and 1")
    
    # Load input sequences
    try:
        sequences = list(SeqIO.parse(args.input, 'fasta'))
        if not sequences:
            parser.error(f"No sequences found in {args.input}")
        logging.info(f"Loaded {len(sequences)} sequence(s) from {args.input}")
    except Exception as e:
        parser.error(f"Failed to read input file: {e}")
    
    # Parse ratios
    try:
        ratios = parse_ratios(args.ratios, len(sequences))
        if len(sequences) > 1:
            ratio_str = ', '.join(f"{s.id}:{r:.1%}" for s, r in zip(sequences, ratios))
            logging.info(f"Sequence ratios: {ratio_str}")
    except ValueError as e:
        parser.error(str(e))
    
    # Generate synthetic reads
    logging.info(f"Generating {args.num_reads} reads with {args.error_rate:.1%} error rate")
    if args.seed is not None:
        logging.info(f"Using random seed: {args.seed}")
    
    reads = generate_reads(
        sequences=sequences,
        num_reads=args.num_reads,
        error_rate=args.error_rate,
        ratios=ratios,
        seed=args.seed
    )
    
    # Write output
    try:
        with open(args.output, 'w') as f:
            SeqIO.write(reads, f, 'fastq')
        logging.info(f"Wrote {len(reads)} reads to {args.output}")
        
        # Report statistics
        total_bases = sum(len(r.seq) for r in reads)
        avg_length = total_bases / len(reads) if reads else 0
        logging.info(f"Total bases: {total_bases:,}, Average read length: {avg_length:.1f}")
        
    except Exception as e:
        logging.error(f"Failed to write output: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()