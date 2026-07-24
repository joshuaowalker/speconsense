# Architecture: Variant Significance, err_factor, and IUPAC Handling

How speconsense decides whether a variant is real, and the sequence-comparison
machinery that decision rests on.

## Three-stage architecture

Phasing does **not** gate on significance. The pipeline instead:

1. **Phases indiscriminately** — split by haplotype wherever MSA position analysis supports it
2. **Groups by identity** — complete linkage on `--group-identity` (default `0.85`) via
   `scipy.cluster.hierarchy`. Every pair in a group must meet the threshold, which prevents the
   transitive chains that would otherwise collapse closely related variants in eDNA-style mixtures
3. **Annotates pairwise via CER** — after the fact, per candidate

The same identity groups gate read reassignment and discard recovery: reads can only move
within their identity group.

## CER (critical error rate)

Key pieces:

- `significance.compute_critical_error_rate(N, M, L, alpha, K)` — p* under a uniform model
  (q=p/3), with combinatorial Bonferroni correction for `K>1` multi-position variants
- `context.classify_variant_context()` — produces one `ContextTag` per variant event
  (substitution or contiguous indel block). **HP context comes from the reference consensus**,
  because the artifact hypothesis under test is that the candidate's reads are miscalled copies
  of the reference
- `qctx.get_qctx(tag, table)` — per-position error rate for a context; HP runs longer than the
  table's max route to blanket homopolymer normalization

The shipped pipeline reports `cer_factor` (per-position multiplicative inflation) directly from
the joint q* solver. The uniform-p* form is retained for the variant-significance paper's tables.

### No gating in core

Every non-anchor candidate is pairwise-compared against all larger peers in its identity group,
annotated with a `cer_factor` (worst-case across peers), and flows through to the FASTA output
regardless of value. The reference pool accumulates all processed clusters regardless of factor —
`min_factor` is inherently conservative for artifact-vs-artifact cases.

**Summarize** applies the user-visible routing via `--min-cer-factor` (default `1.0`; `0`
disables). Records with `cer_factor=None` — anchors, or a failed pairwise comparison — always
pass the quality gates.

The FASTA header carries only `cer_factor=`. Full per-position detail (p*, K, context tags,
q_ctx values, reference idx) lives in the metadata JSON via `_build_variant_record`.

## err_factor (cluster homogeneity)

Complementary to CER and **peer-independent**. For each non-gap consensus column, `err_factor` is
the fraction of reads disagreeing with the consensus divided by the q_ctx rate predicted for that
column's context (HP run length, or non-HP), aggregated cluster-wide.

- Values near 1.0 — reads consistent with basecaller noise
- Values ≫ 1.0 — internal heterogeneity beyond what sequencing noise would produce

Because it needs no peer, it distinguishes novel-but-real sequences (low `err_factor`) from noise
combinations (high). Computed in `msa.compute_cluster_err_factor(msa_string, qctx_table)` during
final consensus generation; emitted as `err_factor=` in the FASTA header and stored with raw
`obs_sum`/`exp_sum`/`cols` in the metadata JSON for reproducibility.

Summarize filters on it via `--max-err-factor` (default `1.5`; `0` disables). The `1.5` default is
safe because MAD outlier removal at final consensus (`speconsense.outliers.detect_rid_outliers`)
strips the single-read outliers that would otherwise inflate `err_factor` on real clusters.

## Identity rank (gid/vid) ordering

Assigned in `_assign_identity_ranks` (core, Phase 13):

- **Groups** are ranked by anchor (largest-member) size descending, so `gid=1` is always the
  largest group.
- **Variants within a group** are ranked by a quality-aware tier
  `(expected_to_pass_summarize, size)` descending — clusters expected to survive summarize's
  default cer/err filters get the low vids, then size descending within each tier.

This keeps primary-track vids contiguous in the common case (no gap where a low-`cer_factor` or
high-`err_factor` variant gets routed to `.ns`/`.lq`) while preserving traceability: vids are
stamped once here, and summarize never renumbers them.

The tier boundaries are the shared constants `DEFAULT_MIN_CER_FACTOR` (`significance.py`) and
`DEFAULT_MAX_ERR_FACTOR` (`msa.py`), which also back summarize's `--min-cer-factor` /
`--max-err-factor` defaults, so the two cannot drift apart.

Because the largest cluster can be demoted below `vid=1` when its `err_factor` is high,
`fit_error_model` re-derives the primary anchor by max read count rather than trusting the
`1.v1` label.

## IUPAC ambiguity handling

- `IUPAC_CODES` (`msa.py`) maps nucleotide sets to codes (R=A/G, Y=C/T, …)
- `IUPAC_EQUIV` (`distances.py`) lets edlib treat ambiguity codes as matching their constituent bases
- `STANDARD_ADJUSTMENT_PARAMS` (`distances.py`) defines consistent comparison parameters:
  - homopolymer normalization enabled (treats `AAA` = `AAAAA`)
  - IUPAC overlap disabled — standard IUPAC semantics, so Y≠M
  - no end trimming (`end_skip_distance=0`)
  - single-base repeats for homopolymer normalization
- Safeguards `MIN_COVERAGE_THRESHOLD=0.5` and `MAX_ADJUSTMENT_RATIO=1.5` fall back to raw edlib
  identity when terminal-gap exclusion would inflate adjusted identity on length-mismatched
  sequences.

## MCL graph construction (weakly-held design decision)

The K-NN similarity graph for MCL (`scalability/base.py`) uses an **asymmetric edge storage
pattern**: `similarities[id1]` only contains entries for neighbors `id2 > id1` (lexicographically).
This contradicts the MCL documentation's recommendation of symmetric graphs, but produces good
clustering results.

The pattern emerged from the original implementation and affects tie-breaking when multiple
neighbors have identical similarity scores. Switching to symmetric storage would change which
neighbors are selected during K-NN construction, and therefore change clustering results.
Validate against the existing test cases before changing it.
