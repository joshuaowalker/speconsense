"""Per-position variant context classification for context-aware CER.

Given two consensus sequences (a candidate cluster and a validated reference
cluster from the same identity group), classify each pairwise difference into
a ContextTag. The tags drive q_ctx lookup so that the unified CER equation can
weight per-position error probabilities by the context of each differing
position.

The classifier produces one ContextTag per "variant event" — a substitution or
contiguous indel block — to match the K count produced by
count_variant_differences.

HP context is determined from the reference consensus, since the artifact
hypothesis under test is that the candidate cluster's reads are miscalled
copies of the reference's sequence.
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple

import edlib

from speconsense.distances import IUPAC_EQUIV, bases_match_with_iupac


class ContextClass(Enum):
    """Coarse class of a variant event.

    NON_HP_SUB: substitution at a non-HP position OR within an HP run.
        HP-interior and HP-interface substitutions are not given special
        treatment because empirical analysis shows their rate is
        indistinguishable from non-HP substitutions.

    NON_HP_INDEL: indel at a non-HP position. Rare in amplicon data; treated
        as a degenerate context with a placeholder rate.

    HP_LENGTH_CHANGE: indel within an HP run that changes the run's apparent
        length. Tagged with the HP base, reference run length, and the change
        direction (negative for deletions, positive for insertions).
    """
    NON_HP_SUB = "non-hp-sub"
    NON_HP_INDEL = "non-hp-indel"
    HP_LENGTH_CHANGE = "hp-length-change"


@dataclass(frozen=True)
class ContextTag:
    """Per-position context for a single variant event.

    For ``NON_HP_SUB`` and ``NON_HP_INDEL``, only ``cls`` is meaningful.

    For ``HP_LENGTH_CHANGE``:
        base: the HP base ('A'/'C'/'G'/'T').
        length: the HP run length in the reference (number of consensus bases).
        direction: signed integer; negative = deletion of |direction| bases,
            positive = insertion of |direction| bases.
    """
    cls: ContextClass
    base: Optional[str] = None
    length: Optional[int] = None
    direction: Optional[int] = None

    def to_string(self) -> str:
        """Compact string form used in JSON metadata CER reproduction records."""
        if self.cls == ContextClass.NON_HP_SUB:
            return "non-hp-sub"
        if self.cls == ContextClass.NON_HP_INDEL:
            return "non-hp-indel"
        if self.cls == ContextClass.HP_LENGTH_CHANGE:
            assert self.base is not None and self.length is not None and self.direction is not None
            dir_token = f"del{abs(self.direction)}" if self.direction < 0 else f"ins{abs(self.direction)}"
            return f"hp-l{self.length}-{self.base}-{dir_token}"
        return self.cls.value


# ---------------------------------------------------------------------------
# HP run detection (reference-aligned coordinates)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HPRun:
    """A homopolymer run in an aligned reference sequence.

    Attributes:
        base: HP base ('A','C','G','T').
        length: Number of consensus (non-gap) bases in the run.
        start_col: First column (in the aligned reference) belonging to the run,
            including any leading gap-extension.
        end_col: Last column (inclusive) belonging to the run, including any
            trailing gap-extension.
    """
    base: str
    length: int
    start_col: int
    end_col: int


def build_hp_runs(consensus_aligned: str, min_length: int = 2) -> List[Optional[HPRun]]:
    """Return per-column HP run info for an aligned consensus.

    For each column in the alignment, the returned list contains the HPRun
    that column belongs to (including internal gaps and adjacent gap
    extensions), or None if the column is not part of an HP run.

    Args:
        consensus_aligned: Aligned consensus sequence (may contain '-').
        min_length: Minimum number of consensus bases (non-gap) for a run to
            qualify. Default 2; runs of length 1 are not HP for our purposes.

    Returns:
        List of length len(consensus_aligned) of Optional[HPRun].
    """
    n = len(consensus_aligned)
    runs_by_col: List[Optional[HPRun]] = [None] * n

    non_gap = [(i, consensus_aligned[i]) for i in range(n) if consensus_aligned[i] != '-']
    if not non_gap:
        return runs_by_col

    def _emit_run(base: str, positions: List[int]) -> None:
        if len(positions) < min_length:
            return
        # Tag all MSA columns from first to last position (inclusive of internal gaps)
        run = HPRun(base=base, length=len(positions),
                    start_col=positions[0], end_col=positions[-1])
        for col in range(run.start_col, run.end_col + 1):
            runs_by_col[col] = run

    cur_base = non_gap[0][1]
    cur_positions = [non_gap[0][0]]
    for i in range(1, len(non_gap)):
        pos, base = non_gap[i]
        if base == cur_base:
            cur_positions.append(pos)
        else:
            _emit_run(cur_base, cur_positions)
            cur_base = base
            cur_positions = [pos]
    _emit_run(cur_base, cur_positions)

    # Extend HP run tagging into adjacent gap columns whose nearest non-gap
    # neighbor is the run's base. Iterate to stability so chains of gap
    # columns propagate.
    changed = True
    while changed:
        changed = False
        for col in range(n):
            if consensus_aligned[col] != '-' or runs_by_col[col] is not None:
                continue
            left_run = runs_by_col[col - 1] if col > 0 else None
            right_run = runs_by_col[col + 1] if col < n - 1 else None
            chosen: Optional[HPRun] = None
            if left_run and right_run:
                if left_run is right_run:
                    chosen = left_run
            elif left_run:
                chosen = left_run
            elif right_run:
                chosen = right_run
            if chosen is not None:
                runs_by_col[col] = chosen
                changed = True

    return runs_by_col


# ---------------------------------------------------------------------------
# Pairwise classification
# ---------------------------------------------------------------------------


def classify_pairwise_differences(
    candidate_consensus: str,
    reference_consensus: str,
    hp_min_length: int = 2,
) -> Optional[List[ContextTag]]:
    """Classify each variant event between candidate and reference.

    Aligns candidate to reference with edlib (IUPAC-aware), then walks the
    alignment producing one ContextTag per substitution or contiguous indel
    block. The list length equals K as counted by count_variant_differences,
    which is the number of variant events used in CER.

    Args:
        candidate_consensus: Ungapped candidate cluster consensus.
        reference_consensus: Ungapped validated cluster consensus (artifact
            hypothesis source).
        hp_min_length: Minimum HP run length (in reference, in consensus bases)
            to qualify as an HP context. Default 2.

    Returns:
        List of ContextTag, one per variant event, in alignment order.
        Empty list if sequences are identical.
        None if alignment fails or inputs are invalid.
    """
    if not candidate_consensus or not reference_consensus:
        return None

    if candidate_consensus == reference_consensus:
        return []

    try:
        result = edlib.align(
            candidate_consensus, reference_consensus,
            task="path", additionalEqualities=IUPAC_EQUIV,
        )
        if result["editDistance"] == -1:
            return None
        alignment = edlib.getNiceAlignment(result, candidate_consensus, reference_consensus)
        if not alignment:
            return None
    except Exception:
        return None

    candidate_aligned = alignment.get('query_aligned')
    reference_aligned = alignment.get('target_aligned')
    if not candidate_aligned or not reference_aligned:
        return None

    return _classify_from_alignment(candidate_aligned, reference_aligned, hp_min_length)


def _classify_from_alignment(
    candidate_aligned: str,
    reference_aligned: str,
    hp_min_length: int,
) -> List[ContextTag]:
    """Walk a pre-aligned candidate/reference pair and emit ContextTags."""
    assert len(candidate_aligned) == len(reference_aligned)
    runs = build_hp_runs(reference_aligned, min_length=hp_min_length)

    tags: List[ContextTag] = []
    n = len(candidate_aligned)
    i = 0
    while i < n:
        c = candidate_aligned[i]
        r = reference_aligned[i]

        if bases_match_with_iupac(c, r):
            i += 1
            continue

        if c == '-' or r == '-':
            indel_is_deletion_in_candidate = (c == '-')
            block_start = i
            block_len = 0
            while i < n:
                this_c = candidate_aligned[i]
                this_r = reference_aligned[i]
                if not (this_c == '-' or this_r == '-'):
                    break
                # Direction must remain the same (don't merge a candidate-gap run
                # with a reference-gap run as a single block).
                this_is_del = (this_c == '-')
                if this_is_del != indel_is_deletion_in_candidate:
                    break
                block_len += 1
                i += 1

            tags.append(_classify_indel_block(
                runs=runs,
                block_start=block_start,
                block_len=block_len,
                deletion_in_candidate=indel_is_deletion_in_candidate,
                candidate_aligned=candidate_aligned,
                reference_aligned=reference_aligned,
            ))
            continue

        tags.append(ContextTag(cls=ContextClass.NON_HP_SUB))
        i += 1

    return tags


def _classify_indel_block(
    runs: List[Optional[HPRun]],
    block_start: int,
    block_len: int,
    deletion_in_candidate: bool,
    candidate_aligned: str,
    reference_aligned: str,
) -> ContextTag:
    """Classify a single contiguous indel block into NON_HP_INDEL or HP_LENGTH_CHANGE.

    A block is HP_LENGTH_CHANGE if all of the following hold:
      * Every column in the block falls within the same HP run (in the reference).
      * Every gapped base in the block matches that run's base.

    Otherwise the block is NON_HP_INDEL.
    """
    # Identify the run associated with the first column of the block.
    first_run = runs[block_start]
    if first_run is None:
        return ContextTag(cls=ContextClass.NON_HP_INDEL)

    run_base = first_run.base
    is_hp_block = True
    for col in range(block_start, block_start + block_len):
        col_run = runs[col]
        if col_run is None or col_run is not first_run:
            is_hp_block = False
            break
        # Confirm the gapped base matches the run base.
        if deletion_in_candidate:
            base_at = reference_aligned[col]
        else:
            base_at = candidate_aligned[col]
        if base_at != '-' and base_at.upper() != run_base.upper():
            is_hp_block = False
            break

    if not is_hp_block:
        return ContextTag(cls=ContextClass.NON_HP_INDEL)

    direction = -block_len if deletion_in_candidate else block_len
    return ContextTag(
        cls=ContextClass.HP_LENGTH_CHANGE,
        base=run_base,
        length=first_run.length,
        direction=direction,
    )
