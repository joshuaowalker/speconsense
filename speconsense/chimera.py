"""Two-parent recombinant (PCR chimera) detection.

The CER framework tests one artifact hypothesis: "this candidate's reads are
miscalled copies of a single reference" (see context.py). A PCR chimera — a
recombinant of two real haplotypes from the same identity group — is the
inverse case: it differs from parent A at many positions (all inherited from
parent B) and vice versa. Real linked differences make the joint q* collapse,
so a chimera receives a *high* cer_factor and looks like the strongest
variant. err_factor cannot catch it either: a chimeric cluster's reads are
internally homogeneous. This module closes that hole with a direct
recombinant test.

The test is uchime-style (Edgar 2011), simplified because our position is
stronger than de novo chimera detection over raw amplicons: the candidate
parents are already known (the larger peers in the candidate's identity
group), abundance is trustworthy (cluster read counts, not dereplicated read
multiplicities), and the candidates are denoised consensuses, so a true
chimera matches its two-parent model near-perfectly (the same reasoning that
led UCHIME2/3 to a perfect-model criterion with a strong abundance prior).

For a candidate C and each pair of eligible parents (A, B):

1. SPOA-align [A, B, C] with the pipeline's linear gap params.
2. Collect *diagnostic columns*: positions where A and B are distinct plain
   ACGT bases. Gap and IUPAC columns are excluded, so homopolymer wobble and
   ambiguity codes don't pollute the signal.
3. Best single-parent model: mismatches of C to A or to B over those columns,
   whichever is smaller.
4. Best chimera model: one breakpoint splitting the columns into a prefix
   explained by one parent and a suffix explained by the other (both
   orientations scanned, O(k)).
5. Flag when the chimera model is near-perfect, beats the single-parent model
   decisively, and has real support on both sides of the breakpoint.

Thresholds were validated empirically on the ont98 dataset (2026-08-12 probe:
10,480 candidates, 171 flagged, 65 pass-track flags with >=10 diagnostic
sites explained perfectly; flagged pass-track records had median
cer_factor 15 — confirming that CER actively promotes chimeras).

Abundance eligibility follows uchime's abskew: parents must have more than
DEFAULT_ABSKEW times the candidate's read count (strictly larger at the
default of 1.0). A chimera forms partway through PCR and undergoes fewer
doublings than either parent, so this guards against flagging a large
haplotype as a recombinant of rarer "parents".
"""

import logging
import os
import subprocess
import tempfile
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

# Parents must have MORE THAN this multiple of the candidate's read count
# (strict inequality — uchime --abskew). uchime defaults to 2.0 (16 for
# denoised data), but the ont98 validation that calibrated the gates below
# used strictly-larger parents, and real minor haplotypes routinely sit
# within 2x of the chimeras they parent (e.g. a 9-read haplotype parenting
# a 5-read chimera). The near-perfect-partition gates carry the
# specificity; abundance is only a weak prior here because rDNA
# copy-number variation decouples read counts from template abundance.
DEFAULT_ABSKEW = 1.0

# A candidate is tested against pairs drawn from its MAX_PARENTS largest
# eligible peers. Chimera parents are almost always abundant, so this caps
# the O(peers^2) pair loop at negligible recall cost.
MAX_PARENTS = 6

# Flagging gates (see module docstring for empirical provenance):
MIN_DIAGNOSTIC_SITES = 6   # min columns where the parents disagree
MAX_CHIMERA_MISMATCHES = 1  # sites the chimera model may leave unexplained
MIN_IMPROVEMENT = 4         # sites chimera model must explain beyond single
MIN_SIDE_SUPPORT = 3        # sites agreeing with each parent on its own side

_ACGT = frozenset("ACGT")


class ChimeraHit(NamedTuple):
    """Result of a positive chimera test for one candidate.

    parent_a / parent_b are caller-scoped keys (core passes subcluster
    indices) identifying the pair whose recombination best explains the
    candidate. Orientation "A|B" means the candidate's prefix matches
    parent_a and its suffix matches parent_b; "B|A" is the reverse.
    """
    parent_a: int
    parent_b: int
    score: int              # single_mismatches - chimera_mismatches
    diagnostic_sites: int
    single_mismatches: int
    chimera_mismatches: int
    left_support: int
    right_support: int
    orientation: str


def diagnostic_columns(a_row: str, b_row: str, c_row: str) -> List[Tuple[str, str, str]]:
    """Columns of a 3-row MSA where the two parents are distinct plain bases.

    Returned in MSA order as (a_base, b_base, c_base) triples. Gap and
    IUPAC-containing columns are excluded so indel placement, homopolymer
    wobble, and ambiguity codes cannot masquerade as diagnostic signal.
    """
    return [
        (a, b, c)
        for a, b, c in zip(a_row, b_row, c_row)
        if a in _ACGT and b in _ACGT and a != b
    ]


def scan_breakpoint(diag: Sequence[Tuple[str, str, str]]) -> Tuple[int, int, int, int, str]:
    """Best single-parent and best one-breakpoint chimera model over diag.

    Returns (single_mismatches, chimera_mismatches, left_support,
    right_support, orientation). left/right_support count sites agreeing
    with the assigned parent on each side of the best breakpoint.
    """
    k = len(diag)
    pref_mm_a = [0]
    pref_mm_b = [0]
    for a, b, c in diag:
        pref_mm_a.append(pref_mm_a[-1] + (c != a))
        pref_mm_b.append(pref_mm_b[-1] + (c != b))
    tot_mm_a = pref_mm_a[-1]
    tot_mm_b = pref_mm_b[-1]
    m_single = min(tot_mm_a, tot_mm_b)

    best_m = m_single
    best = (0, 0, "-")
    for j in range(k + 1):
        m_ab = pref_mm_a[j] + (tot_mm_b - pref_mm_b[j])
        if m_ab < best_m:
            best_m = m_ab
            best = (j - pref_mm_a[j], (k - j) - (tot_mm_b - pref_mm_b[j]), "A|B")
        m_ba = pref_mm_b[j] + (tot_mm_a - pref_mm_a[j])
        if m_ba < best_m:
            best_m = m_ba
            best = (j - pref_mm_b[j], (k - j) - (tot_mm_a - pref_mm_a[j]), "B|A")
    return m_single, best_m, best[0], best[1], best[2]


def evaluate_triple(a_row: str, b_row: str, c_row: str) -> Optional[Tuple[int, int, int, int, int, int, str]]:
    """Apply the flagging gates to one aligned (parentA, parentB, candidate).

    Returns (score, k, m_single, m_chimera, left, right, orientation) when
    the candidate looks chimeric, else None.
    """
    diag = diagnostic_columns(a_row, b_row, c_row)
    k = len(diag)
    if k < MIN_DIAGNOSTIC_SITES:
        return None
    m_single, m_chimera, left, right, orientation = scan_breakpoint(diag)
    if (m_chimera <= MAX_CHIMERA_MISMATCHES
            and (m_single - m_chimera) >= MIN_IMPROVEMENT
            and left >= MIN_SIDE_SUPPORT
            and right >= MIN_SIDE_SUPPORT):
        return m_single - m_chimera, k, m_single, m_chimera, left, right, orientation
    return None


def _spoa_msa_rows(seqs: Sequence[str]) -> Optional[List[str]]:
    """Align seqs with the pipeline's SPOA invocation; return aligned rows.

    Input order anchors the SPOA graph and therefore the coordinate frame —
    callers must pass sequences in a deterministic order.
    """
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
            for i, s in enumerate(seqs):
                f.write(f">s{i}\n{s}\n")
            temp_input = f.name
        try:
            result = subprocess.run(
                ["spoa", temp_input,
                 "-r", "2", "-l", "1", "-m", "1", "-n", "-1", "-g", "-1", "-e", "-1"],
                capture_output=True, text=True, check=True)
        finally:
            os.unlink(temp_input)
    except (subprocess.CalledProcessError, OSError) as e:
        logging.debug(f"SPOA failed during chimera check: {e}")
        return None

    rows: List[str] = []
    cur: List[str] = []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if cur:
                rows.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        rows.append("".join(cur))
    # spoa emits the input rows in order (a consensus record, if requested,
    # would follow them; -r 2 emits MSA rows only for our invocation).
    if len(rows) < len(seqs):
        return None
    return rows[:len(seqs)]


def detect_chimera(candidate_seq: str,
                   parents: Sequence[Tuple[int, str]]) -> Optional[ChimeraHit]:
    """Test candidate_seq against every pair of the given parents.

    parents: (key, consensus) pairs, already abundance-filtered and ordered
    deterministically by the caller (size descending); at most MAX_PARENTS
    are considered. Returns the strongest ChimeraHit, or None. Ties keep the
    first (largest-parent) pair, so results are deterministic.
    """
    parents = list(parents)[:MAX_PARENTS]
    best: Optional[ChimeraHit] = None
    for i in range(len(parents)):
        for j in range(i + 1, len(parents)):
            key_a, seq_a = parents[i]
            key_b, seq_b = parents[j]
            rows = _spoa_msa_rows([seq_a, seq_b, candidate_seq])
            if rows is None:
                continue
            verdict = evaluate_triple(rows[0], rows[1], rows[2])
            if verdict is None:
                continue
            score, k, m_single, m_chimera, left, right, orientation = verdict
            if best is None or score > best.score:
                best = ChimeraHit(
                    parent_a=key_a, parent_b=key_b, score=score,
                    diagnostic_sites=k, single_mismatches=m_single,
                    chimera_mismatches=m_chimera, left_support=left,
                    right_support=right, orientation=orientation)
    return best
