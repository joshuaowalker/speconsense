# Architecture: Summarize Naming, Tracks, and Merge Semantics

How `speconsense-summarize` names records, routes them to output tracks, and recomputes
fields when records collapse together. See `docs/customizing-fasta-headers.md` for the
user-facing field reference and `docs/understanding-ric-and-merging.md` for the conceptual
treatment of RiC and merging.

## Output structure

**Core (`speconsense`)**
- `{sample_name}-all.fasta` — all consensus sequences for a specimen
- `cluster_debug/`
  - `{sample_name}-c{num}-RiC{size}-reads.fastq` — original reads
  - `{sample_name}-c{num}-RiC{size}-sampled.fastq` — sampled reads for consensus
  - `{sample_name}-c{num}-RiC{size}-untrimmed.fasta` — untrimmed consensus

**Summarize** — `__Summary__/`, with a dual namespace:
- original clustering: `-c1`, `-c2`, `-c3`
- summarization: `-1.v1`, `-1.v2`, `-2.v1` (group and variant)

Four tracks in `__Summary__/variants/`:
- pass — `{name}-RiC{ric}.fasta`
- `.ns` — `{name}.ns-RiC{ric}.fasta`
- `.lq` — `{name}.lq-RiC{ric}.fasta`
- `.filtered` — `{name}.filtered-RiC{ric}.fasta`

Pre-merge contributors: `{name}.raw.{gid}.v{vid}-RiC{ric}.fasta`, named by the contributor's
core `gid.vid` for direct traceability.

## Four-track variant routing

Every variant goes to exactly one track:

- **pass** — quality gates met and selected for review. Written to `-all.fasta` and `variants/`.
- **`.ns`** (not significant) — `cer_factor` below `--min-cer-factor` (default `1.0`). Likely a
  sequencing artifact relative to a larger peer.
- **`.lq`** (low quality) — `err_factor` above `--max-err-factor` (default `1.5`). Internal
  cluster heterogeneity exceeds basecaller noise expectations. **`.lq` wins** when both `.ns` and
  `.lq` thresholds fire.
- **`.filtered`** — passes quality gates but excluded by a selection or pruning decision:
  `--select-max-variants`, `--select-min-size-ratio`, `--select-max-groups`, or
  `--prune-group-ratio`/`--prune-group-count`.

The `.ns`/`.lq` versus `.filtered` split is deliberate: it distinguishes a *quality problem*
from a *selection decision*.

## Naming policy (`process_single_specimen`)

Summarize preserves core's `gid`/`vid` on every record that is **not moved between groups by
cross-primer conflation**.

- Variants dropped by within-group MSA merging, `--select-max-variants`, `--select-min-size-ratio`,
  `--min-cer-factor`, or `--max-err-factor` **leave their vids as gaps**. Gaps are expected output,
  not a bug.
- Variants moved into a survivor group by cross-primer conflation adopt the survivor's `gid` and
  get a freshly-minted `vid` strictly above the highest vid core ever wrote under that gid.
  Collision avoidance considers passed + `.ns` + `.lq` + `.filtered` records under that gid, plus
  any vids already minted in the same pass (this covers 3-or-more-group conflation).
- Vids under absorbed-group gids stay on disk under their original gid in the `.ns`/`.lq`/
  `.filtered` outputs, and do not block the survivor's namespace.

## Post-merge re-check and merge unwind

Two distinct mechanisms, both preserving per-cluster detail rather than emitting a bad synthetic:

- **Re-check (after Phase 4)** — merging recomputes `err_factor`/`cer_factor` on merged records.
  Any merged record whose recomputed metrics cross filter thresholds has its **merge cancelled**:
  the original pre-merge contributors are restored on the pass track, rather than routing the
  synthetic merged record to `.lq`/`.ns`.
- **Unwind (Phase 6/7)** — when a merged record is dropped by selection filtering, the merge is
  unwound: original pre-merge contributors are emitted individually on the `.filtered` track
  instead of the synthetic merged record, preserving per-cluster metrics for reviewer detail.

## Merge-time field handling

When ≥2 contributors collapse into one record (within-group MSA merge or cross-primer overlap
conflation), summarize aggregates per-field via `_build_merged_consensus_info` plus a follow-up
recompute pass in `process_single_specimen`:

- `size` / `ric` — summed
- `length` / `ambig` — re-derived from the column-voted consensus
- `rawric` / `rawlen` — carry flattened merge history
- `primers` — union of contributor primer names, sorted
- `rid` / `rid_min` / `err_factor` — re-derived from a SPOA MSA over the union of contributor
  reads (loaded from each contributor's `cluster_debug` FASTQ)
- `cer_factor` — recomputed against larger peers in the post-merge bucket for same-primer merges;
  set to `None` for cross-primer merges, because the CER noise model is per-locus
- `snp_count` — cumulative across iterative rounds, and can **over-count by 0–2** when the same
  physical position becomes ambiguous in multiple rounds. `ambig` is the canonical IUPAC-site count.

`.raw` pre-merge files are named `.raw.{gid}.v{vid}` using the contributor's original core
`gid.vid`, and carry every field from their source cluster (matching `.ns`/`.lq` pass-through),
resetting only the merge-history fields.

Recompute requires the contributors' debug FASTQs and the per-specimen metadata JSON's
`parameters.error_model`. Missing inputs fall back to the values inherited from the largest
contributor.

## Frequency fields

`group_frequency=` and `global_frequency=` — opt-in via `--fasta-fields full` or explicit listing.
Each is the variant's `size` as a percentage of a per-specimen denominator:

- `group_frequency` uses the **conflation-aware bucket total** — the sum of `size` over every
  record (passed + `.ns` + `.lq` + `.filtered`) in the same post-cross-primer-conflation bucket,
  so a moved variant is measured against the merged-group total
- `global_frequency` uses `total_input_reads` from the metadata JSON — the post-presample-cap
  count fed into clustering, **not** the literal cap

Denominators are computed in `process_single_specimen` after Phase 1b and stashed on each
`ConsensusInfo` (`group_size_total`, `global_size_total`). Both fields are suppressed when
denominators are unavailable (legacy inputs without `gid=`; specimens with missing metadata),
and both propagate through within-group MSA merges, the `-{gid}-full` builder, and the
`.raw` / `.ns` / `.lq` / `.filtered` writers.

## Locus labeling

`locus=` — opt-in via `--fasta-fields full`, `--fasta-fields locus`, or explicit listing.
Summarize classifies each consensus via `pyitsx.classify()` and stamps `locus=ITS`, `locus=ITS1`,
or `locus=ITS2`. The organism group comes from core's metadata (`parameters.pyitsx_organism`,
always recorded, default `F`).

Initialization is lazy: `_init_locus_labeler` checks `pyitsx.is_available()` and creates a
`ProfileDB` once per run. Silently omitted when pyitsx is unavailable or the organism's profiles
are not installed. **Independent of `--orient-mode`** — locus labeling works even with
`--orient-mode none`.

## Group full consensus (`--enable-full-consensus`)

Implemented in `speconsense/summarize/full_consensus.py`. For each identity group with ≥2 selected
variants on the pass track, summarize emits an additional `-{gid}-full` record.

The read pool is a size-weighted, top-mean-Phred sample of the **pre-merge core variants** in the
bucket (sourced from `variant_groups[final_gid]`, so within-group MSA merging and selection filters
do not constrain it). Contributors are gated by a running-total threshold using
`parameters.min_ambiguity_frequency` from core's metadata JSON (default 0.10): each contributor
must satisfy `size ≥ min_ambiguity_frequency × running_total`.

Sampled reads (budget = `parameters.max_sample_size`, default 100) go through SPOA with linear gap
scoring, using `-l 0` (local SW) when the gated bucket spans multiple primer sets and `-l 1`
(global NW) otherwise. The MSA is collapsed via `build_full_consensus_from_msa` with
one-vote-per-read and majority-wins gaps; columns where ≥2 bases each clear
`min_ambiguity_frequency` get IUPAC codes.

Output: a FASTA entry in `-all.fasta` and a sampled-reads FASTQ in `FASTQ Files/`. No `.raw`
lineage, since the record is synthetic. Suppressed when fewer than 2 pass-track variants exist for
the group, fewer than 2 contributors clear the gate, or the read pool cannot be loaded. `.ns` and
`.lq` records are excluded from the pool by construction (filtered upstream by
`load_consensus_sequences`).

**Intended use is BLAST query against legacy unphased ITS references with adjusted-identity
scoring.** `-full` is IUPAC-bearing and will silently degrade under raw BLAST.
