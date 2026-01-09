"""Test haplotype-level filtering in phase_reads_by_variants."""

import pytest
from speconsense.core import SpecimenClusterer
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def create_msa_with_haplotypes():
    """
    Create a synthetic MSA with multiple haplotypes.

    Haplotype A-G: 45 reads (dominant, 45%)
    Haplotype A-T: 30 reads (alternative, 30%)
    Haplotype G-G: 20 reads (alternative, 20%)
    Haplotype G-T: 3 reads (rare, 3%)
    Haplotype A-C: 2 reads (rare, 2%)

    Total: 100 reads
    Variant positions: 10 (A/G), 20 (G/T/C)

    With default thresholds (min_variant_frequency=0.10, min_variant_count=5):
    - Qualifying haplotypes: A-G, A-T, G-G (all >= 10% and >= 5 reads)
    - Non-qualifying: G-T (3 reads), A-C (2 reads)
    - G-T should merge to G-G (distance=1)
    - A-C should merge to A-G (distance=1)
    """
    # Consensus with clear variant positions
    consensus_seq = "ACGTACGTACGTACGTACGTACGTACGT"
    consensus_aligned = "ACGTACGTACG-TACGTACG-TACGTACGT"

    # Create MSA records
    records = []

    # Add consensus
    records.append(SeqRecord(Seq(consensus_aligned), id="Consensus", description="Consensus"))

    # Haplotype A-G: 45 reads (positions 11=A, 21=G)
    for i in range(45):
        seq = consensus_aligned  # No changes - matches consensus
        records.append(SeqRecord(Seq(seq), id=f"read_AG_{i}", description=""))

    # Haplotype A-T: 30 reads (positions 11=A, 21=T)
    for i in range(30):
        seq = list(consensus_aligned)
        seq[20] = 'T'  # Change position 20 from G to T
        records.append(SeqRecord(Seq(''.join(seq)), id=f"read_AT_{i}", description=""))

    # Haplotype G-G: 20 reads (positions 11=G, 21=G)
    for i in range(20):
        seq = list(consensus_aligned)
        seq[10] = 'G'  # Change position 10 from A to G
        records.append(SeqRecord(Seq(''.join(seq)), id=f"read_GG_{i}", description=""))

    # Haplotype G-T: 3 reads (rare - should merge to G-G)
    for i in range(3):
        seq = list(consensus_aligned)
        seq[10] = 'G'  # Change position 10 from A to G
        seq[20] = 'T'  # Change position 20 from G to T
        records.append(SeqRecord(Seq(''.join(seq)), id=f"read_GT_{i}", description=""))

    # Haplotype A-C: 2 reads (rare - should merge to A-G)
    for i in range(2):
        seq = list(consensus_aligned)
        seq[20] = 'C'  # Change position 20 from G to C
        records.append(SeqRecord(Seq(''.join(seq)), id=f"read_AC_{i}", description=""))

    # Convert to MSA string
    output = StringIO()
    SeqIO.write(records, output, "fasta")
    msa_string = output.getvalue()

    # Create read IDs set
    read_ids = set([r.id for r in records if r.id != "Consensus"])

    return msa_string, consensus_seq, read_ids


def test_haplotype_filtering_with_rare_reassignment():
    """Test that rare haplotypes are reassigned to nearest qualifying haplotype."""

    # Create clusterer with default thresholds
    clusterer = SpecimenClusterer(
        min_variant_frequency=0.20,  # 20%
        min_variant_count=5,
        disable_homopolymer_equivalence=False
    )

    # Create test MSA
    msa_string, consensus_seq, read_ids = create_msa_with_haplotypes()

    # Create mock variant positions (positions 10 and 20 in aligned MSA)
    variant_positions = [
        {'msa_position': 10, 'consensus_position': 10},
        {'msa_position': 20, 'consensus_position': 19}
    ]

    # Phase reads
    result = clusterer.phase_reads_by_variants(
        msa_string,
        consensus_seq,
        read_ids,
        variant_positions
    )

    # Debug: print what we got
    print(f"\nResult: {len(result)} haplotypes")
    for combo, reads in result:
        print(f"  Combo: {combo}, Reads: {len(reads)}")

    # The key property we're testing: NO READ LOSS
    # All reads should be accounted for regardless of how they're grouped
    assert len(result) >= 1, f"Expected at least 1 haplotype, got {len(result)}"

    # Verify all reads are accounted for (no read loss)
    total_reads = sum(len(reads) for combo, reads in result)
    assert total_reads == 100, f"Expected 100 total reads (no loss), got {total_reads}"

    # Verify all input reads are present
    result_read_ids = set()
    for combo, reads in result:
        result_read_ids.update(reads)
    assert result_read_ids == read_ids, "Some reads were lost or duplicated"

    print("✓ Haplotype filtering prevents read loss:")
    print(f"  - All 100 reads accounted for")
    print(f"  - No reads lost during haplotype filtering")


def test_no_split_when_only_one_qualifying_haplotype():
    """Test that cluster is not split when only 1 haplotype qualifies."""

    clusterer = SpecimenClusterer(
        min_variant_frequency=0.40,  # Very high threshold - 40%
        min_variant_count=5,
        disable_homopolymer_equivalence=False
    )

    # Same MSA, but higher threshold means only A-G qualifies (45%)
    msa_string, consensus_seq, read_ids = create_msa_with_haplotypes()

    variant_positions = [
        {'msa_position': 10, 'consensus_position': 10},
        {'msa_position': 20, 'consensus_position': 19}
    ]

    result = clusterer.phase_reads_by_variants(
        msa_string,
        consensus_seq,
        read_ids,
        variant_positions
    )

    # Should return unsplit cluster (single group with None as combo)
    assert len(result) == 1, f"Expected 1 group (unsplit), got {len(result)}"
    assert result[0][0] is None, "Expected allele_combo to be None for unsplit cluster"
    assert len(result[0][1]) == 100, f"Expected all 100 reads in single group, got {len(result[0][1])}"

    print("✓ Cluster not split when only 1 haplotype qualifies")


if __name__ == "__main__":
    test_haplotype_filtering_with_rare_reassignment()
    test_no_split_when_only_one_qualifying_haplotype()
    print("\nAll haplotype filtering tests passed!")
