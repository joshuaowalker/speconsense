"""Lightweight read records for the clustering pipeline.

BioPython ``SeqRecord`` objects cost ~45KB per ONT read once the per-base
quality int list is materialized — prohibitive for eDNA-scale inputs
(``--presample 0`` on a 1M-read FASTQ would need ~45GB). ``ReadRecord``
stores the same information as three strings (~1.5KB/read): the full FASTQ
title line, the sequence, and the ASCII-33 quality string.

Parity contracts with the SeqRecord path they replace:

- ``mean_quality`` equals ``statistics.mean(letter_annotations["phred_quality"])``
  exactly: sum(ord(c) - 33) over the ASCII string is the same integer sum,
  and a single exact-int division produces the identical float.
- ``write_fastq`` emits ``@{title}\\n{seq}\\n+\\n{qual}\\n`` — byte-identical
  to ``SeqIO.write(record, handle, "fastq")`` for records parsed from FASTQ,
  where ``record.description`` is the full title line and BioPython writes
  ``@{description}`` with a bare ``+`` separator and unwrapped lines.
- ``reverse_complemented`` mirrors orientation's in-place SeqRecord mutation:
  reverse-complemented sequence, reversed quality string.
"""

from typing import Iterator, List, NamedTuple

from Bio.Seq import reverse_complement


class ReadRecord(NamedTuple):
    """One sequencing read: id (first title token), full title, seq, ASCII qual."""
    id: str
    title: str
    seq: str
    qual: str

    def mean_quality(self) -> float:
        """Mean Phred quality; exactly equals the SeqRecord-based computation."""
        if not self.qual:
            return 0.0
        qb = self.qual.encode('ascii')
        return (sum(qb) - 33 * len(qb)) / len(qb)

    def reverse_complemented(self) -> "ReadRecord":
        return self._replace(
            seq=str(reverse_complement(self.seq)),
            qual=self.qual[::-1],
        )


def from_seqrecord(record) -> ReadRecord:
    """Convert a BioPython SeqRecord (FASTQ- or FASTA-parsed, or test-built)."""
    quals = record.letter_annotations.get("phred_quality") if record.letter_annotations else None
    qual = ''.join(chr(q + 33) for q in quals) if quals else ""
    title = record.description if record.description else record.id
    return ReadRecord(id=record.id, title=title, seq=str(record.seq), qual=qual)


def parse_fastq(path: str) -> Iterator[ReadRecord]:
    """Fast 4-line FASTQ reader (~10x SeqIO.parse, no per-base int lists)."""
    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            header = header.rstrip('\n')
            if not header.startswith('@'):
                raise ValueError(f"Malformed FASTQ header in {path}: {header[:50]!r}")
            seq = handle.readline().rstrip('\n')
            plus = handle.readline()
            if not plus.startswith('+'):
                raise ValueError(f"Malformed FASTQ separator in {path} after {header[:50]!r}")
            qual = handle.readline().rstrip('\n')
            title = header[1:]
            read_id = title.split(None, 1)[0] if title else title
            yield ReadRecord(id=read_id, title=title, seq=seq, qual=qual)


def parse_input(path: str, fmt: str) -> List[ReadRecord]:
    """Parse an input file into ReadRecords. FASTQ uses the fast reader;
    FASTA falls back to SeqIO (quality strings empty)."""
    if fmt == "fastq":
        return list(parse_fastq(path))
    from Bio import SeqIO
    return [from_seqrecord(r) for r in SeqIO.parse(path, fmt)]


def write_fastq(handle, record: ReadRecord) -> None:
    """Write one FASTQ record, byte-identical to BioPython's writer for
    FASTQ-parsed records."""
    handle.write(f"@{record.title}\n{record.seq}\n+\n{record.qual}\n")
