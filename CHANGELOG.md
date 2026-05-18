# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.1] - 2026-05-18

Correctness pass on the pipeline + significant summarize-side improvements. Headlines:

- **MAD outlier removal moved before CER validation.** Previously the shipped FASTA described a post-MAD cluster while its `cer_factor` annotation described a pre-MAD cluster. All annotations now describe the cluster actually shipped to disk.
- **Post-refinement merge (Phase 10).** Collapses clusters that converge into identical (or HP-equivalent) consensuses after read reassignment / discard recovery / second phasing, catching cases the earlier merges (operating on pre-MAD consensuses) miss.
- **Summarize naming preserves core's `gid`/`vid`** except for cross-primer-conflated records, which get freshly-minted vids above the highest vid core ever produced under the survivor's gid.
- **Summarize merge-time field handling recomputes per-cluster metrics** (`rid`, `rid_min`, `err_factor`, `cer_factor`) on the union of contributor reads so merged records' annotations describe the merged cluster, not the largest pre-merge contributor.
- **`--enable-full-consensus` reintroduced** with a principled construction policy (noise-scrubbed size-weighted sample, linear-gap SPOA, frequency-gated IUPAC, majority-wins gaps) — empirically beats SNP-merged v1 on multi-haplotype specimens against legacy MycoMap, with no losses exceeding 0.15 pp on the validation cohort. Re-enabled in the `compressed` profile.
- **CLI help is tiered.** Default `--help` shows only user-facing flags; pre-tuned algorithm internals, MAD knobs, and integral-phase disable flags move to `--help-advanced`.
- **Two new opt-in FASTA fields: `group_frequency=` and `global_frequency=`.** Per-variant percentages against the conflation-aware identity bucket and against the specimen's presampled input total, respectively.

### Added

- **Post-refinement merge (Phase 10)** (core) — A new phase between MAD-cleaned consensus generation (Phase 9) and CER validation (Phase 11). It collapses clusters whose post-MAD `trimmed_consensus` values are identical (or HP-equivalent under `--hp-normalization-length`) into a single survivor, then re-runs the Phase 9 worker on the union of read_ids so the survivor's MSA, consensus, rid metrics, and IUPAC count reflect the merged read set. Catches convergence introduced by Phases 6–9 — read reassignment (6), discard recovery (7), second phasing (8), and MAD itself (9) can leave two clusters with the same shipped consensus, which the earlier Phase 2 and Phase 4 merges (operating on pre-MAD SPOA consensuses) cannot detect. Honors `--disable-cluster-merging` exactly as Phases 2 and 4 do. Affected metadata: `cer_group_N` and `M` reflect the merged cluster's read counts where merges occurred; identity-group composition shrinks by exactly the merge count.
- **MAD tuning knobs exposed as CLI flags** (core) — Four new flags — `--mad-z-threshold`, `--mad-gap-factor`, `--mad-min-mad`, `--mad-min-drop-from-median` — surface the per-knob tunables of `speconsense.outliers.detect_rid_outliers` for calibration. The same keys (with dashes) are accepted in profiles.
- **`speconsense-summarize --enable-full-consensus` reintroduced** with a principled construction policy that supersedes the all-indel union removed in 0.8.0 (`4dc3797`). Per identity group with ≥2 selected variants on the pass track, summarize draws a size-weighted, top-mean-Phred read sample (budget = `parameters.max_sample_size` from the core metadata JSON, default 100) across the pre-merge core variants that clear a running-total gate (each contributor must satisfy `size ≥ min_ambiguity_frequency * running_total`, with `min_ambiguity_frequency` also sourced from the core run, default 0.10). The sampled reads go through SPOA with linear gap scoring (`-m1 -n-1 -g-1 -e-1`, matching core's per-cluster consensus); the resulting MSA is collapsed to a single consensus with one-vote-per-read, majority-wins gaps, and IUPAC ambiguity codes at columns where ≥2 bases each clear the `min_ambiguity_frequency` threshold. Cross-primer-conflated groups switch to SPOA local mode (`-l 0`) so reads with mid-stream-aligned coverage from alternative primer pairs don't get penalized for end gaps. The output is named `-{gid}-full`, written to the main `-all.fasta` alongside variants, with a companion sampled-reads FASTQ in `FASTQ Files/`. Suppressed when the group has only 1 selected variant on the pass track or fewer than 2 contributors clear the gate (degenerate single-source consensus). `.ns` and `.lq` records are not eligible to contribute. Re-enabled in the `compressed` profile. Intent (per the `~/mm/analysis/blast_query` POSITION paper, 2026-04-30): a query substrate for BLAST against legacy unphased ITS references — what NGSpeciesID would have produced on the input FASTQ if the original reads had been pre-filtered for noise. Note: `-full` is IUPAC-bearing by construction and should be queried with adjusted-identity scoring (e.g., MycoBLAST); under raw BLAST it carries the IUPAC penalty documented in the POSITION paper.
- **`--disable-full-consensus` CLI override** (summarize) — Symmetric with `--enable-full-consensus` so users running with the `compressed` profile (which now sets `enable-full-consensus: true`) can suppress the `-full` artifact on a per-invocation basis. Mirrors the `--disable-merging` / `--enable-merging` pattern.
- **`group_frequency=` and `global_frequency=` FASTA fields** (summarize) — Two new `--fasta-fields` choices. `group_frequency` is the variant's `size` as a percentage of its conflation-aware identity-bucket total (sum of `size` across passed + `.ns` + `.lq` records in the same post-conflation bucket — so cross-primer-conflated variants are measured "in the context of the merged groups"). `global_frequency` is the variant's `size` as a percentage of the specimen's `total_input_reads` (the post-presample-cap count fed into clustering, from the metadata JSON). Both fields are **suppressed by default** — they appear only in the `full` preset, or when listed explicitly (`--fasta-fields default,group_frequency`). When the denominator is unknown (legacy inputs without `gid=` headers; specimen with missing metadata) the field is silently omitted. Both propagate through within-group MSA merges, cross-primer conflation, the `.raw` / `.ns` / `.lq` writers, and the `-{gid}-full` group consensus.

### Changed

#### Pipeline reorder (MAD before CER)

- **MAD outlier removal now runs before CER validation.** Previously MAD lived inside the final-consensus worker (Phase 11), after CER (Phase 9) had already annotated each cluster against the pre-MAD M, N, and consensus. The shipped FASTA therefore described a post-MAD cluster while its `cer_factor` annotation described a pre-MAD cluster. The pipeline now runs MAD-cleaned consensus generation as Phase 9, with the post-refinement merge (Phase 10), CER (Phase 11), size filtering (Phase 12), and output emission (Phase 13) all consuming the post-MAD state. The `err_factor` calculation already used the post-MAD MSA and continues to do so.
- **`--min-size` and `--min-cluster-ratio` now apply to post-MAD read counts.** A cluster whose pre-MAD size exceeds `--min-size` but whose post-MAD size falls below will be dropped where previously it would have survived. The semantic shift is small in practice because MAD typically removes single-read outliers, but it is observable on edge cases.
- **Metadata is now internally consistent on M and N.** Previously `M` was post-MAD (read at metadata-emission time, after MAD had subtracted outliers in Phase 11) but `N = cer_group_N` was pre-MAD (stamped during Phase 9 CER before MAD ran). Both are now post-MAD, so the CER decision recorded in metadata references the same population that's reported in the variant record.
- **Pipeline phases renumbered 1–14** to reflect the new ordering (was 1–12 in 0.8.0).

#### Summarize naming

- **`speconsense-summarize` now preserves core's `gid`/`vid` on output.** Previously summarize renumbered every passing variant positionally (`-{group_idx+1}.v{variant_idx+1}`), which could shift names whenever filters dropped variants or within-group MSA merging consolidated them — and could collide with the unrenamed `.ns`/`.lq` records that already carried core's original names. The new policy: if a cluster is not moved between groups by cross-primer overlap conflation, its name round-trips core's `gid.vid` exactly. Variants dropped by within-group MSA merging, `--select-max-variants`, `--select-min-size-ratio`, `--min-cer-factor`, or `--max-err-factor` leave their vids as gaps in the survivor's namespace. When cross-primer overlap conflation moves a variant into a different (larger) group, it adopts the survivor group's `gid` and is assigned a freshly-minted `vid` strictly above the highest vid core ever produced under that gid (across passed + `.ns` + `.lq` records, plus any vids already allocated earlier in the same pass). The collision-avoidance set is computed per-survivor-gid and is gid-local — vids that core wrote under absorbed groups remain on disk under their original gid in the `.ns`/`.lq` outputs and do not block the survivor's namespace. Any record whose summarize-emitted name differs from its core-emitted name is necessarily a record that crossed a group boundary via cross-primer conflation.

#### Summarize merge-time field handling

- **`speconsense-summarize` now recomputes `rid`, `rid_min`, and `err_factor` on merged records.** Previously a merged record inherited these metrics from its largest pre-merge contributor, so the displayed values described that contributor's old cluster, not the merged cluster. Each ≥2-contributor merge now runs SPOA over the union of contributor reads (loaded from each contributor's `cluster_debug` FASTQ, scored with the same linear-gap parameters core uses) and re-derives the metrics from the resulting MSA. `rid` / `rid_min` use the same homopolymer-normalized per-read identity calculation as core; `err_factor` and its raw `obs_sum` / `exp_sum` / `cols` come from `msa.compute_cluster_err_factor`. The recomputed values describe the merged cluster's internal homogeneity against an MSA built from its union reads — for cross-primer merges, reads with terminal gaps contribute zero to columns they don't cover, so the metrics remain well-defined under primer-pool conflation. If the contributor FASTQs aren't available or SPOA fails, the inherited values are kept as a best-effort fallback (debug-logged).
- **`cer_factor` is recomputed for same-primer merges and set to `None` for cross-primer merges.** Same-primer merges (all contributors share a primer set) re-classify the merged candidate's pairwise differences against every larger peer in its post-merge bucket via `classify_pairwise_differences`, look up each variant event's `q_ctx`, and run `compute_cer_factor` with `N = bucket-total reads`, `M = candidate size`. The minimum factor across peers is the candidate's new `cer_factor` — the same convention core uses. Cross-primer merges (contributors had different primer sets) bypass the pipeline and get `cer_factor=None`: the CER noise model is per-locus by construction and isn't well-defined for merged-locus candidates spanning multiple primer-pool amplicons. `None` is treated by summarize's routing identically to other "no valid peer comparison" cases — record always passes the `--min-cer-factor` filter. Summarize reads `parameters.error_model` and `parameters.significance_level` from each specimen's metadata JSON to pick the same q_ctx table and alpha core used; if metadata is missing, recompute is skipped and the inherited value is kept.
- **Cross-primer merges now report a union of primer names in `primers=`.** Previously the merged record inherited `primers` from its largest contributor, which silently dropped primers contributed by smaller contributors (e.g., a cross-primer ITS + ITS2 merge would show only the ITS primers). The merged record's `primers` field is now `sorted({primer for contributor in contributors for primer in contributor.primers})`, matching the set semantics of `primers_are_same`. Each primer name keeps its `5'-` / `3'-` direction prefix, so a 2-pair cross-primer merge emits e.g. `primers=5'-ITS1F,5'-ITS3,3'-ITS4_RC`.
- **`.raw` pre-merge files now carry every field from the source cluster.** Previously the `.raw` `ConsensusInfo` construction in `write_specimen_data_files` only copied a subset of fields (`size`, `ric`, `primers`, `rid`, `rid_min`), silently dropping `cer_factor`, `err_factor`, the err_factor raw sums (used by the quality report's calibration check), `group_rank`, and `variant_rank`. The `.raw` construction now uses `ConsensusInfo._replace` to inherit every field from the source cluster and only resets the merge-only fields (`snp_count`, `raw_ric`, `raw_len`, `merge_indel_count`). Matches the pass-through treatment that `.ns` and `.lq` records already received.
- **Documented merge-time `snp` cumulation caveat.** `_build_merged_consensus_info`'s docstring now notes that `snp_count` is cumulative across iterative merge rounds and can over-count when the same physical position becomes ambiguous in multiple rounds. The final consensus's `ambig` (non-ACGT character count) is the canonical IUPAC site count. Practical over-count is 0–2 for typical merges; full bookkeeping would require tracking MSA coordinates across rounds and isn't justified for the corner case.

#### CLI help surface

- **`speconsense` and `speconsense-summarize` help is now tiered.** Default `--help` shows only the user-facing flags (Common + Tuning). Pre-tuned algorithm internals, MAD knobs, and integral-phase disable flags move into a new "Advanced (pre-tuned — rarely needed)" group hidden from `--help` and shown under `--help-advanced`. Flags themselves are unchanged — parsing, semantics, and profile-loading all work identically. The trailing line `For pre-tuned/internal flags, use --help-advanced.` points users to the expanded view. Implementation in `speconsense/_help.py` is a ~70-line `argparse.HelpFormatter` subclass + a `_HelpAdvancedAction` that toggles per-action help strings, applied via `install_advanced_help(parser)`. Specifically moved to Advanced in `speconsense`: `--inflation`, `--k-nearest-neighbors`, the four MAD knobs (`--mad-z-threshold`, `--mad-gap-factor`, `--mad-min-mad`, `--mad-min-drop-from-median`), and every `--disable-*` / `--enable-*` flag (position-phasing, read-reassignment, discard-recovery, second-phasing, noise-filter, mad-outlier-removal, ambiguity-calling, cluster-merging, homopolymer-equivalence). In `speconsense-summarize`: `--disable-merging` / `--enable-merging`, `--disable-homopolymer-equivalence` / `--enable-homopolymer-equivalence`. `--presample` stays in the default Performance group.
- **Three previously CLI-only disable flags became valid profile keys** to align with the other `disable-*` flags that were already profile-loadable: `disable-second-phasing`, `disable-noise-filter`, `disable-mad-outlier-removal`. No semantic change — they were always functional on the CLI; the profile validator now accepts them too.

#### Defaults

- **`--mad-z-threshold` default lowered from 3.0 → 1.5.** Empirically tuned on the ont98 ITS cohort (955 specimens, Dorado v5.0 SUP). The modified-Z rule operates on cluster *shape* (one read pulling away from a tight median) rather than absolute magnitude, and primary clusters tend to have the "one read sticks out" pattern that this rule was designed for while noise clusters spread heterogeneity across all reads. Lowering z makes the rule more selective for that pattern. Effects on the validator-facing routing: err_factor routing AUC 0.630 → 0.664 (+5.4% relative); F1 of the .lq/pass routing decision at the existing threshold 0.148 → 0.193 (+30% relative); pass-track records 3,609 → 3,662 (+1.5%); separation between noise and primary `err_factor_p90` +23% relative. Costs: −1 specimen in primary recall, +0.8pp pass-track noise share, per-record IUPAC ambiguity density essentially flat (+1.4%). The other three MAD knobs (`gap-factor`, `min-mad`, `min-drop-from-median`) stay at their previous defaults; the `min-drop-from-median = 0.02` safety floor remains load-bearing and prevents the lower z from over-flagging tight clusters.

### Removed

- **`create_output_structure` removed from the summarize public API.** It was a dead duplicate of the naming logic that lived in `process_single_specimen`; it had no internal or test callers. Removed from `speconsense.summarize` exports.

### Migration notes (0.8.0 → 0.8.1)

- **Re-running summarize against a 0.8.0 `clusters/` directory works** (no header schema change). However, the new merge-time recompute requires `cluster_debug/` FASTQs and the per-specimen `*_metadata.json`; if either is missing for a specimen, the merge-time recompute silently falls back to inherited values from the largest contributor (debug-logged). For full benefit, the underlying core run should also be on 0.8.1.
- **Pre-existing `__Summary__/` directories with renumbered variant names** will be reorganized on re-run: variant names now round-trip core's `gid`/`vid` exactly. Tools that hard-coded post-summarize names (e.g., `-1.v1`, `-1.v2`) against pre-0.8.1 outputs should be re-verified.
- **`-{gid}-full` is a new output filename** when `--enable-full-consensus` is on (or via the `compressed` profile). Tools that parse the variant FASTA must handle the bare `-{gid}` suffix in addition to `-{gid}.v{vid}`.

## [0.8.0] - 2026-05-01

This is a major release. Headline themes:

- **Statistical variant validation (CER framework).** Every non-anchor cluster is annotated with a `cer_factor` measuring how implausible it is as basecaller-noise replicates of a larger peer under a per-position context-aware error model. Filtering happens in `speconsense-summarize`.
- **Cluster homogeneity (`err_factor`).** A peer-independent metric for whether a cluster's reads are internally consistent with the assumed sequencing-noise profile. Routes outlier clusters to a separate `.lq` track.
- **Post-phasing refinement pipeline.** Three new phases (read reassignment, discard recovery, second phasing pass) move reads between clusters and re-admit previously-discarded reads when consensus concordance supports it.
- **Identity groups (`gid`/`vid`) emitted by core, honored by summarize.** Complete-linkage groups are computed once, in core, and travel through the FASTA header into summarize, which no longer re-clusters.
- **Context-aware error models.** Per-basecaller q_ctx tables (HP run length × event type) ship as YAML and are user-overridable.
- **Pipeline phases renumbered 1–12** (replaces the nested `4a/4b/4b2/4b3/4c` scheme).

### Added

#### CER framework (statistical variant significance)

- **Critical Error Rate (CER) variant validation** (core) — Every non-anchor cluster is annotated with a `cer_factor`: the minimum across larger peers of the per-position multiplicative inflation needed for that peer's observed minor-allele counts to be plausible under the basecaller's noise model. Anchors and clusters that fail to find a comparison peer carry `cer_factor=None` and always pass. Implementation lives in `speconsense.significance` (binomial survival with combinatorial Bonferroni correction over candidate sites; uniform error model with `q=p/3`) and `speconsense.context` (per-position context classification).
- **`--significance-level / --min-cer-factor` user controls** — Core: `--significance-level` (default `1e-5`) sets the alpha used in p* solving. Summarize: `--min-cer-factor` (default `1.0`, `0` to disable) routes records below the threshold to `__Summary__/variants/{name}.ns-RiC{ric}.fasta` (and matching FASTQ) instead of dropping them
- **Multi-position significance (K>1)** — When two clusters disagree at K>1 candidate sites, the binomial test uses joint counts and a combinatorial Bonferroni correction over `C(L, K)` site pairs/triples. The site count `L` is now the **HP-compressed length** (each HP run counts once), not the raw sequence length, since multi-base HP errors are dominated by single-event mechanisms
- **Context-aware error model (q_ctx)** — `speconsense/error_models/*.yaml` ships per-basecaller error rate tables keyed by `(non-hp-sub, non-hp-indel, hp-l1 ... hp-lN)`. Each variant event is classified by HP run length (from the **reference** consensus, since the artifact hypothesis is that the candidate's reads are miscalled copies of the reference) and looked up in the table. HP runs longer than the table's max route to blanket HP normalization rather than CER evaluation.
- **`--error-model` flag** (core) — Replaces the earlier `--qctx-profile`. Resolution order: filesystem path → `~/.config/speconsense/error_models/` → bundled. `--list-error-models` lists available models. Default `dorado-v5.0`. Bundled models: `dorado-v5.0`, `dorado-v3.5`. Dorado v5.0 q_ctx values were re-estimated from the ont98 dataset under the current pipeline
- **Per-cluster CER reproduction data in metadata JSON** — `clusters/{specimen}_metadata.json` now includes per-cluster `cer_details` (peer id, K, p*, joint q, per-position context tags and q_ctx values, reference idx) so CER decisions can be replayed offline
- **`cer_factor=` and `cer_details=` FASTA fields** — Replaces the old `cer=` / `cer.a=` fields. The header carries only the scalar factor; full per-position detail lives in metadata JSON

#### Cluster homogeneity (err_factor)

- **`err_factor` cluster homogeneity metric** (core) — For each non-gap consensus column, the fraction of reads disagreeing with the consensus divided by the q_ctx rate predicted for that column's context (HP run length or non-HP). Values near 1.0 mean reads consistent with basecaller noise; values ≫ 1.0 indicate internal heterogeneity beyond what sequencing noise produces. Computed in `msa.compute_cluster_err_factor`; emitted as `err_factor=` in the FASTA header and stored with raw `obs_sum`/`exp_sum`/`cols` in metadata JSON. Unlike `cer_factor`, `err_factor` is peer-independent — it distinguishes novel-but-real sequences (low) from noise combinations (high)
- **`--max-err-factor`** (summarize) — Default `1.5`; `0` disables. Records above threshold route to `__Summary__/variants/{name}.lq-RiC{ric}.fasta`. `.lq` takes precedence over `.ns` when both fire. The 1.5 default is safe because MAD outlier removal at final consensus removes single-read outliers that would otherwise inflate err_factor on real clusters

#### Post-phasing refinement (new pipeline phases)

- **Cross-cluster read reassignment within identity groups** (core, phase 6) — After variant phasing, reads can move between clusters in the same complete-linkage identity group when the new cluster's consensus better explains them. Concordance scoring uses position-aware variant agreement; an edit-distance guard prevents drift across distinct sequences. New flags `--disable-read-reassignment` / `--enable-read-reassignment` (default enabled)
- **Discard recovery** (core, phase 7) — Reads previously dropped (failed phasing, MAD outliers, hard-floor or noise-filter dropouts) are re-tested against surviving cluster consensuses; concordant reads re-enter their best-fitting cluster. New flags `--disable-discard-recovery` / `--enable-discard-recovery` (default enabled). Auto-skipped if `--disable-read-reassignment` is set, since recovery uses the same concordance machinery. Reassigned discards no longer get a `d-` prefix on read IDs
- **Second phasing pass** (core, phase 8) — Re-phases any clusters whose membership changed via reassignment/recovery. Gated by `--enable-position-phasing` AND `--enable-read-reassignment` — if reassignment is off, nothing changed since the first pass and there's no useful work to do
- **Noise filter** (core, phase 5) — `_filter_noisy_clusters` runs SPOA on each small cluster and disbands those with no-majority columns (positions where no base reaches 50%). Reads from disbanded clusters fall to the discard pool; some may return via phase 7
- **No-majority ambiguity calling** — When phasing produces a single qualifying haplotype, IUPAC codes are generated for every variant position. The lower count threshold (`--min-ambiguity-count`, default 3) captures variation that doesn't have enough support to form its own cluster

#### Identity groups (gid/vid) end-to-end

- **Complete-linkage identity grouping** (core) — Replaces the previous chained/single-linkage scheme. Every pair within a group must meet `--group-identity` (default `0.85`), preventing transitive collapse of close-but-distinct variants in eDNA-style mixtures. Implementation uses `scipy.cluster.hierarchy.fcluster`. Identity groups gate read reassignment, discard recovery, and CER validation
- **`gid=` and `vid=` FASTA fields** (core) — Each cluster receives an identity-group rank (`gid`) and within-group variant rank (`vid`). Output filenames switch from `-c{cluster_num}` to `-{gid}.v{vid}` end-to-end (core's `-all.fasta`, summarize's outputs, debug fastqs)
- **Summarize honors core grouping** — `speconsense-summarize` parses `gid`/`vid` from the FASTA header and skips its previous independent HAC re-clustering. Only **cross-primer overlap merging** still runs across core groups, for the primer-pool use case (ITS/ITS1/ITS2). Hard-fails on inputs lacking `gid=`/`vid=` (i.e., refuses to process pre-0.8 outputs without re-running core)
- **Cross-primer overlap merger** (summarize) — Configured by `--group-identity` (default `0.85`, matches core's anchor identity) and `--min-merge-overlap`. Merges different-primer core groups whose anchors overlap above the identity threshold

#### Other additions

- **MAD-based outlier removal in cluster consensus generation** (core) — Replaces the earlier `--outlier-identity` static threshold. Computes per-read identity to consensus, drops reads beyond `n × MAD` from the median (robust to skewed distributions). The MAD pass runs before CER validation, so CER, size filtering, err_factor, and final FASTA emission all observe the post-MAD cluster. Removes single-read outliers that would otherwise inflate `err_factor` and skew CER decisions
- **Per-specimen variant tree** (summarize) — `__Summary__/trees/{specimen}.txt` renders an ASCII hierarchy of every variant in a specimen (passed, `.ns`, `.lq` together), grouped by core identity group. Each non-anchor variant branches from the larger-size peer with the highest pairwise identity, with a one-line edit summary (substitutions, single-nt indels, short ≤3 nt indels, long indels). Filtered variants carry the on-disk status marker (`-1.v4.ns`, `-1.v8.lq`) so they're visually distinct from passed variants that may share the same vid after summarize's renaming
- **`--hp-normalization-length` parameter** (core + summarize) — Unifies the HP threshold across `are_homopolymer_equivalent`, `_classify_subcluster_groups`, `is_homopolymer_event`, distance calculations, and MSA-based variant merging. Default 6 matches the HP error rate paper's recommendation of CER-evaluating L≤5 HP variants. Set to 1 for legacy blanket-normalize-all behavior
- **`build_adjustment_params(hp_normalization_length=N)` factory** (`speconsense.distances`) — Returns an `AdjustmentParams` variant with the requested HP threshold, sharing all other STANDARD settings
- **`--disable-read-reassignment` / `--enable-read-reassignment`** (core) — Gates phase 6
- **`--disable-discard-recovery` / `--enable-discard-recovery`** (core) — Gates phase 7
- **HP error rate analysis paper** — Two new papers under `docs/`: `cer_in_practice/` (the CER framework as shipped) and `hp_error_rate/` (per-context error rate estimation)

### Changed

- **Identity grouping moved entirely to core** (architecture) — Previously summarize ran its own HAC on consensus sequences. Now core computes one set of complete-linkage identity groups; summarize honors them via FASTA-header `gid`/`vid`. Within-group MSA-based variant merging in summarize is unchanged
- **CER filtering moved core → summarize** (architecture) — Core annotates every cluster with `cer_factor` and writes everything to disk. The pass/ns decision is applied by summarize via `--min-cer-factor`. This lets users adjust the threshold without re-running core
- **Bonferroni site count uses HP-compressed length** (core) — Multi-position CER significance was previously divided by the raw sequence length. Now uses HP-compressed length (each run counts once) since multi-base HP errors are correlated, not independent. Slightly reduces conservatism for K>1 variants in HP-rich amplicons
- **`--hp-min-length` renamed to `--hp-normalization-length`** (core) — Same semantics, same default (6); previous name was opaque about what the threshold controls. Profile key `hp-min-length` likewise renamed
- **`--qctx-profile` renamed to `--error-model`** (core) — Same resolution order, same YAML format; previous name overloaded "profile" with the parameter-profile system. The `error_models/` package directory replaces the earlier name
- **Default `--max-err-factor=1.5`** (summarize) — Was unset in 0.7.x. The 1.5 default is safe in conjunction with MAD outlier removal at final consensus
- **Internal cluster-merging and identity-grouping HP threshold honored** (core) — `are_homopolymer_equivalent` and `_classify_subcluster_groups` were silently using a hardcoded 6 instead of `self.min_hp_length`; now respect the user-configured value
- **Per-specimen INFO log rationalized** (core) — Each pipeline stage now emits one uniform `[N clusters, R/T reads, P%]  {stage}: {outcome}` line, with a fixed running-state column on the left. Stage labels reworded for user-facing clarity (`Pre/Post-phasing merge` → `Cluster equivalence merge [(round 2)]`, `After phasing` → `Variant phasing`, `Second phasing pass` → `Variant phasing (round 2)`, `CER annotation` → `Significance testing`). Configuration and procedural messages (e.g. `Position-based variant phasing enabled`, `Running MCL algorithm`) demoted to DEBUG. No-op stages (no merges, no reads moved, no splits, no recovered reads) demoted to DEBUG. Per-cluster MAD outlier drops aggregated into a single `Outlier filter` summary line. The `Final` line is now emitted *after* outlier removal, so its totals reflect what was actually written to disk
- **Quality report rewritten** (summarize) — Replaced the descriptive multi-table dump with an inspection-oriented report. New sections: executive summary (with a count of specimens that had input reads but produced no clusters), passed variants worth inspecting (run-relative outliers on `err_factor` and ambiguity, mean ± 2σ), possible rescues (low-yield specimens whose filtered variants are close to the threshold), overlap merge details (only when present), run-wide parameter signals (suggesting tuning), pipeline activity. The CER and err_factor filter decisions now appear in the report; the "Read Identity" and "Positional Identity" outlier tables are dropped (subsumed by the err_factor outlier ranking)
- **Quality report yield computed from cluster size, not RiC** (summarize) — RiC reflects sampling cap, not biological yield; cluster size is the right denominator for per-specimen yield estimation
- **Variants qualified for phasing on full cluster OR quality-biased sample** (core) — Either a position passes the variant frequency threshold across the full cluster, or it passes within a quality-biased sample (top-quality reads). Catches signal that was being diluted out by low-quality reads under the old all-or-nothing test
- **Bundled profiles aligned with new pipeline** — `herbarium`: disable CER and err_factor filters in summarize (high recall; degraded DNA inflates err_factor and weak-signal variants warrant review). `largedata`: tighten core `group-identity` to 0.95 (matches summarize; smaller identity groups reduce work in phases 6/7/9 on bloated groups). `nostalgia`: disable read reassignment, discard recovery, CER and err_factor filters (older pipelines lacked all four). `strict`: migrate `min-cluster-ratio` (global) to `select-min-size-ratio` (per-group; preserves precision intent without suppressing secondary primer pairs in multi-target work)

### Fixed

- **CER solver numerical bounds for high-K variants** — The joint-q solver could overflow or underflow at very small p*. Tightened bisection bounds and clamped intermediate values
- **Specimens with all variants filtered now emit `.ns`/`.lq` files** (summarize) — The per-specimen processing loop previously iterated only over file paths in the passing list. Specimens whose every variant was routed to `.ns` or `.lq` were skipped entirely and their filter outputs were lost. The loop now iterates the union of file paths across passing/`.ns`/`.lq` lists
- **Merged variants now retain `gid`/`vid`** (summarize) — `merge_variants_with_msa` was constructing the merged `ConsensusInfo` without copying `group_rank` / `variant_rank` from the largest input variant, leaving them unset on every post-merge record. Now propagated

### Removed

- **`--enable-full-consensus` removed** (summarize) — The all-indel IUPAC consensus output (`.full`) was dominated on every metric tested in the blast_query analysis: no measurable benefit under adjusted identity, severe regression under raw BLAST (~58 pp at ≥99%) because no real haplotype contains all the unioned indels at once. The companion SNP-merging path (`--merge-indel-length`, `--merge-position-count`) is the recommended query output and is unchanged. The `compressed` profile no longer sets `enable-full-consensus: true`. User profiles still setting `enable-full-consensus:` will fail validation. Existing `.full` records in pre-existing `__Summary__/` directories are no longer recognized by `--aggregate-only`; regenerate by re-running summarize
- **`--enable-early-filter` / `--disable-early-filter` removed** (core) — The early filter applied `--min-size` and `--min-cluster-ratio` after pre-phasing merge to skip variant phasing on doomed-small clusters. In practice the flag was disabled by default, no bundled profile enabled it, and the unconditional hard-floor filter (drop clusters < 3 reads) covers the most useful case. Filtering by `--min-size` / `--min-cluster-ratio` continues at the post-phasing size-filtering phase. User profiles still setting `enable-early-filter:` will fail validation
- **`--outlier-identity` removed** (core) — Static identity threshold for outlier removal replaced by MAD-based detection at final consensus. No knob — runs unconditionally
- **HAC clustering code path in summarize removed** — Summarize no longer re-groups consensus sequences; identity groups come from core via `gid`/`vid`. Cross-primer overlap merger is the only remaining cross-group operation
- **`cer=` / `cer.a=` FASTA fields removed** — Replaced by `cer_factor=` (header) and `cer_details=` (metadata JSON)

### Pipeline numbering

- **Phases renumbered sequentially 1–12** (core) — The previous nested numbering (`4a`/`4b`/`4b2`/`4b3`/`4c`) was historical accretion that obscured the linear flow. Phases now run in plain sequence: (1) initial clustering → (2) pre-phasing merge → (3) variant phasing → (4) post-phasing merge → (5) noise filter → (6) read reassignment → (7) discard recovery → (8) second phasing pass → (9) CER validation → (10) size filtering → (11) output generation → (12) discarded-reads writeout. User-facing log labels were already name-based and are unchanged. CLI flags, profile keys, and metadata field names are unchanged

### Performance

- **Noise-filter SPOA parallelized** (core) — `_filter_noisy_clusters` now dispatches per-cluster SPOA via `ProcessPoolExecutor` when `--threads N>1` and more than 10 small clusters are queued, mirroring the existing parallelization pattern in `merge_similar_clusters`. Sequential fallback preserved for `--threads 1`. At 100K-read inputs this phase used to hang for tens of minutes on serial SPOA; with `--threads 8` it now completes in low single-digit minutes
- **Identity grouping sparsified for large cluster counts** (core) — `_form_identity_groups` (called by read reassignment, discard recovery, and CER validation) now uses `ScalablePairwiseOperation.compute_distance_matrix` when the cluster count exceeds 50 AND the scalability threshold (default 1001 reads). vsearch finds candidate pairs at the relaxed identity threshold (`group_identity × 0.9 = 0.765`); each candidate is rescored under the same HP-normalized metric as the dense path; missing pairs default to distance 1.0, which under complete linkage forces the corresponding clusters into separate groups. Reduces ~7.8M edlib calls per invocation (at 2792 clusters) to a few thousand candidate refinements. Falls back to the dense O(n²) loop when scalability is disabled, vsearch is missing, or the cluster count is small. `_get_scalable_operation` now accepts an optional `scoring_function` kwarg so call sites can override the default similarity metric
- **Discard screening sparsified for large discard pools** (core) — `_run_discard_reassignment`'s screening pass replaces the per-discard scan over all cluster consensuses with a batched vsearch query (top-K=10) refined by edlib on the candidate slice — same metric as the dense path. Activates when the discard pool exceeds the scalability threshold. At 100K reads with 15K discards this collapses ~42M edlib calls to a single vsearch index build plus ~150K candidate refinements. Falls back to the dense per-discard loop when scalability is disabled
- **CER validation within-group top-K for large identity groups** (core) — `_validate_identity_group` now compares each candidate against its top-K=10 most similar larger peers rather than all larger peers, when the group exceeds 50 clusters. The reported `cer_factor` is the minimum across compared pairs, dominated by the most-similar peer (smaller pairwise distance → smaller K → tighter Bonferroni → smaller factor); top-K preserves this minimum with very high probability. Min over a subset is always ≥ min over the full set, so the approximation can only over-report factors, never under-report — i.e., approximation errors result in slightly less filtering, never over-filtering. vsearch returns top-30 most similar peers per query (3× headroom for the case where many of a candidate's most-similar peers are smaller-than-self); the loop intersects with the strict-larger-peer reference pool and takes the first 10. Falls back to the dense O(K²) pairwise comparison when group size ≤ 50 or scalability is disabled

### Dependencies

- `adjusted-identity` requirement bumped to `>=0.2.7` for the `hp_normalize_min_length` `AdjustmentParams` field

### Migration notes (0.7.x → 0.8.0)

- **Re-run core before re-running summarize.** Summarize hard-fails on FASTA inputs lacking `gid=`/`vid=` headers. Pre-0.8 `clusters/` directories cannot be summarized in place — re-run `speconsense` to regenerate
- **Drop `enable-full-consensus`, `enable-early-filter`, `--outlier-identity` from user profiles.** All three keys/flags now fail validation
- **Rename `--qctx-profile` → `--error-model` and `--hp-min-length` → `--hp-normalization-length`** in any wrapper scripts. Profile keys likewise renamed (`qctx-profile` → `error-model`, `hp-min-length` → `hp-normalization-length`)
- **Filename schema changed from `-c{N}` to `-{gid}.v{vid}` everywhere.** Tools that parse speconsense filenames need updating. The original-namespace `-c{N}` debug filenames are no longer emitted by core
- **`cer=`/`cer.a=` FASTA fields are gone.** Tools reading FASTA headers should switch to `cer_factor=` (and `err_factor=` for the new homogeneity metric)

## [0.7.7] - 2026-03-17

### Added
- **`--specimen` flag** (summarize) — Process only a single specimen, outputting per-specimen files and a JSON summary to stdout. Designed for incremental invocation by specimux-suite orchestration layer
- **`--aggregate-only` flag** (summarize) — Regenerate aggregate `summary.fasta` from existing per-specimen outputs without reprocessing. Parses metadata (`raw_ric`, `raw_len`, `snp_count`) back from FASTA headers for full-fidelity reconstruction

### Fixed
- **Stale output file cleanup** — Single-specimen mode removes previous output files before reprocessing to prevent accumulation when variant counts change between runs
- **Temp log file leak** — Early-return paths (e.g., no sequences found) now clean up the temporary log file
- **Glob escaping in specimen cleanup** — Specimen IDs containing glob metacharacters no longer match unintended files

## [0.7.6] - 2026-02-08

### Changed
- **`.full` consensus suppressed for single-variant groups** — When `--enable-full-consensus` is active but a group has only one surviving post-merge variant, the `.full` output is now omitted (it would be redundant). `.full` is still generated when multiple post-merge variants survive in the group, even if `--select-max-variants` limits the number actually output

### Fixed
- **Truncated sequences no longer distort HAC grouping** — A short consensus (e.g., 220bp vs 660bp) could cause `calculate_adjusted_identity_distance` to return near-zero distance due to terminal gap exclusion scoring only the overlap region. This caused complete-linkage HAC to split full-length sequences (>99% identical) into separate groups. Added coverage and adjustment-ratio guards that fall back to raw edlib identity when adjusted identity is unreliable

## [0.7.5] - 2026-02-05

### Fixed
- **Full consensus now includes only surviving variants' components** — `.full` consensus previously applied `merge-min-size-ratio` independently against pre-merge sizes, which used a different (looser) reference than the post-merge `select-min-size-ratio` filter. This caused `.full` RiC to exceed the sum of variant RiCs when merging inflated the denominator. Now `.full` is built from exactly the pre-merge variants that contributed to post-merge variants surviving `select-min-size-ratio`

## [0.7.4] - 2026-02-05

### Added
- **`--select-min-size-ratio` parameter** (summarize) — Filters out post-merge variants whose size ratio to the largest variant in their group falls below the threshold (default: 0 = disabled, e.g. 0.2 for 20% cutoff). Applied after merging, before variant selection

### Changed
- **`compressed` bundled profile** — Now includes `select-min-size-ratio: 0.2` to match 20% calling threshold theme

## [0.7.3] - 2026-02-04

### Added
- **`--enable-full-consensus` flag** (summarize) — Generates a full IUPAC consensus per variant group from all pre-merge variants (gaps never win). Output named `{specimen}-{group}.full`
- **`compressed` bundled profile** — Aggressive variant merging with indel support, 20% frequency thresholds, and full consensus enabled
- **Complement flags for all boolean CLI options** — Allows CLI override of profile-set boolean values (e.g., `--enable-merging` overrides `--disable-merging` set by a profile)

### Changed
- Extracted shared metadata assembly into `_build_merged_consensus_info()` helper in `merging.py` (reduces duplication across three consensus-from-MSA functions)

## [0.7.0] - 2026-01-05

### Added
- **Scalability mode** - Automatic vsearch-based acceleration for large datasets (1000+ sequences), reducing O(n²) pairwise comparisons to ~O(n log n)
  - `--scale-threshold N` sets activation threshold (default: 1001, 0 to disable)
  - `--threads N` controls internal parallelism (speconsense default: 1, summarize default: auto)
  - Requires vsearch (`conda install bioconda::vsearch`), falls back to brute-force if missing
- **Profile system** - YAML-based parameter presets for reproducible workflows
  - Use with `-p herbarium`, `-p strict`, `-p nostalgia`, or `-p largedata`
  - User profiles in `~/.config/speconsense/profiles/`
  - Bundled profiles: `herbarium` (high-recall), `strict` (high-precision), `nostalgia` (simulate older pipelines), `largedata` (large input files)
- **Merge control options** (summarize):
  - `--disable-merging` skips MSA-based merge evaluation entirely (fastest when merging not needed)
  - `--merge-effort` controls merge thoroughness: `fast`, `balanced` (default), `thorough`, or numeric 6-14
- **Length filtering** (summarize) - `--min-len` and `--max-len` filter sequences before processing
- **Primer-constrained overlap merging** - Overlap-aware distance only applies when sequences have different primer pairs, preventing chimeras from incorrectly merging with shorter amplicons
- **Hybrid linkage for HAC clustering** - Uses single-linkage when overlap merging is enabled (for transitive overlap relationships), complete linkage otherwise
- **Performance options** - `--enable-early-filter` skips phasing on clusters that will be filtered; `--collect-discards` writes discarded reads for inspection
- **Cluster boundary delimiters** - Merged FASTQ files include synthetic boundary records between clusters for easier inspection
- **Profile parameters reference** - New `docs/profile-parameters.md` with complete parameter tables for all profiles

### Changed
- **Default `--min-cluster-ratio`**: 0.20 → 0.01 (keep small clusters for downstream curation)
- **Default `--min-variant-frequency`**: 0.20 → 0.10 (unified with ambiguity threshold)
- **Default `--max-sample-size`**: 500 → 100 (sufficient for consensus quality)
- **Default `--merge-min-size-ratio`** (summarize): 0.0 → 0.1 (prevent small clusters from adding IUPAC to large ones)
- **Reduced MAX_MSA_MERGE_VARIANTS** from 10 to 8 for better performance in large variant groups
- **Code architecture** - Refactored `core.py` into `core/` subpackage and `summarize.py` into `summarize/` subpackage
- **CLI organization** - Options grouped into logical categories in `--help` output

### Fixed
- **`--merge-snp` option** - Now correctly supports `--no-merge-snp` to disable SNP merging from command line
- **`--collect-discards` option** - Fixed incomplete discard tracking that missed some filtered reads
- **README documentation** - Fixed incorrect defaults and examples (variant-frequency, merge-min-size-ratio, missing tools/profiles)

## [0.6.6] - 2026-01-01

### Fixed
- **`raw_ric` metadata loss during iterative merging** - Fixed two bugs causing `rawric` field and variant files to be missing
  - Bug 1: When iteration 2+ found no new merges, single-variant pass-through incorrectly called `create_*_from_msa()` which reset `raw_ric` to None
  - Bug 2: `raw_ric` was not flattening prior merge history (unlike `raw_len`), causing iterative merges to show summed values instead of original cluster RiCs
  - Result: `rawric` field missing from FASTA headers, `.raw` variant files not created in `__Summary__/variants/`
  - Fix: Skip `create_*_from_msa()` for single-variant subsets; flatten `raw_ric` like `raw_len` to preserve original values

## [0.6.5] - 2025-12-27

### Fixed
- **Merge traceability loss during iterative merging** - Fixed bug where FASTQ output files contained fewer reads than expected
  - When iterative merging occurred in overlap mode, the merge traceability dictionary was being overwritten rather than expanded
  - This caused loss of original cluster names when a merged result got merged again
  - Result: incomplete FASTQ file concatenation (fewer reads than RiC value indicated)
  - Fix: traceability now expands intermediate merges to preserve all original cluster names through the merge chain

## [0.6.4] - 2025-12-20

### Fixed
- **IUPAC code expansion during variant merging** - Fixed bug where merging a base with an existing IUPAC code would incorrectly produce 'N'
  - Example: C + Y previously produced N, now correctly produces Y (since Y = C|T)
  - Added `merge_bases_to_iupac()` helper that expands existing IUPAC codes before combining
  - Affects all three merge code paths in summarize.py

## [0.6.3] - 2025-12-18

### Added
- **Separate IUPAC ambiguity calling thresholds** - Independent parameters for ambiguity calling vs variant phasing
  - New `--min-ambiguity-frequency` parameter (default: 0.10 = 10%, lower than phasing threshold)
  - New `--min-ambiguity-count` parameter (default: 3, lower than phasing threshold)
  - Allows conservative phasing while capturing more residual variation as IUPAC codes
  - Parameters recorded in run metadata JSON

### Fixed
- **Empty input file handling** - Fixed crash (ZeroDivisionError) when input file contains no sequences
  - Now exits gracefully with warning message and exit code 0
  - Added defensive check in size filtering to prevent division by zero

## [0.6.1] - 2025-12-18

### Added
- **Overlap merge feature for primer pools** - Merge sequences of different lengths when they share sufficient overlap
  - New `--min-merge-overlap` parameter (default: 200bp, 0 to disable)
  - Enables primer pool workflows where reads have overlapping but different coverage
  - Uses single-linkage HAC clustering in overlap mode for better grouping
  - Supports iterative merging for 3+ overlapping sequences
  - Containment case handled: allows merge when shorter sequence is fully contained
- **`rawlen` FASTA field** - Shows original sequence lengths before overlap merging
  - Format: `rawlen=630+361` (sorted descending by length)
  - Accumulates through iterative merges: `rawlen=630+361+269`
  - Added to default and full field presets
- **Quality report overlap merge section** - New analysis section in quality_report.txt
  - Groups merges by specimen with iteration details
  - Shows overlap bp, prefix/suffix extensions
  - Warns on edge cases (overlap near threshold, large length ratios)
- **Improved overlap merge logging** - Enhanced log messages with extension details
  - Format: `(overlap=360bp, prefix=100bp, suffix=169bp)`
  - DEBUG level shows merge span positions

### Fixed
- **Specimen name truncation** - Fixed parsing of specimen names containing '-c'
  - Example: `test-core-overlap-c1` now correctly extracts `test-core-overlap`
  - Changed `split('-c')` to `rsplit('-c', 1)` throughout

## [0.6.0] - 2025-12-05

### Changed
- **Code architecture refactoring** - Extracted MSA analysis functions into dedicated `msa.py` module
  - Moved `IUPAC_CODES`, `ErrorPosition`, `ReadAlignment`, `PositionStats`, `MSAResult` data structures
  - Moved `parse_score_aligned_for_errors`, `extract_alignments_from_msa`, `analyze_positional_variation` functions
  - Moved variant phasing functions: `is_variant_position_with_composition`, `call_iupac_ambiguities`, `calculate_within_cluster_error`, `group_reads_by_allele_combo`, `filter_qualifying_haplotypes`, `reassign_to_nearest_haplotype`, `select_correlated_positions`
  - Improves code organization and testability
- **Refactored cluster() method** - Split 318-line monolithic method into 6 focused phase methods
  - `_run_initial_clustering()` - Phase 1: MCL or greedy clustering
  - `_run_prephasing_merge()` - Phase 2: Merge clusters before variant detection
  - `_run_variant_phasing()` - Phase 3: Detect variants and phase reads
  - `_run_postphasing_merge()` - Phase 4: Merge subclusters after phasing
  - `_run_size_filtering()` - Phase 5: Apply size/ratio filters
  - `_write_cluster_outputs()` - Phase 6: Generate output files
  - Main `cluster()` method now orchestrates these phases cleanly
- **Narrowed exception handling** - Replaced broad `except Exception` with specific exception types
  - `run_spoa()`: Catches `FileNotFoundError`, `OSError`, `subprocess.CalledProcessError`, `RuntimeError`
  - `identify_outlier_reads()`: Catches `ZeroDivisionError`, `AttributeError`, `TypeError`
  - `calculate_read_identity()`: Catches `ZeroDivisionError`, `ValueError`, `TypeError`

### Added
- **Unit tests for variant phasing** - Added 14 new tests in `tests/test_variant_phasing.py`
  - Tests for `is_variant_position_with_composition()` edge cases
  - Tests for `calculate_within_cluster_error()` with various haplotype configurations
  - Tests for `select_correlated_positions()` exhaustive and beam search modes
  - Tests for `call_iupac_ambiguities()` basic functionality
  - Tests for `IUPAC_CODES` mapping correctness

## [0.5.0] - 2025-11-06

### Added
- **Homopolymer-aware variant merging** - MSA-based merging now distinguishes structural indels from homopolymer length differences
  - Analyzes SPOA alignment to classify each indel event (consecutive run) as structural or homopolymer
  - Homopolymer indels (e.g., AAA vs AAAA) are ignored when checking merge compatibility by default
  - Structural indels (true insertions/deletions) count against merge limits
  - Matches adjusted-identity semantics used throughout the pipeline
  - Enables more aggressive merging of biologically equivalent sequences while preserving structural variation
- **`--disable-homopolymer-equivalence` option** - Allows strict identity merging when homopolymer length variation should be preserved
  - Treats all indels (both homopolymer and structural) as blocking merges
  - Useful for applications where homopolymer length variation is biologically significant

### Fixed
- **Indel event counting** - Corrected `--merge-position-count` to count indel events instead of indel columns
  - Previously: A 3bp deletion counted as 3 positions, incorrectly limiting merging
  - Now: Each consecutive indel run counts as 1 event, matching documented behavior
  - Example: `--merge-position-count 3` now correctly allows up to 3 indel events (of any length ≤ `--merge-indel-length`)
  - Event-based approach groups consecutive indel columns, then classifies each complete event as homopolymer or structural
- **SNP detection** - Fixed to exclude alignment columns containing gaps
  - Previously: Columns with gaps and multiple bases (e.g., `['-', 'G', 'A']`) were incorrectly counted as both SNPs and indels
  - Now: SNPs are only columns with multiple different bases and NO gaps
  - Columns with gaps are exclusively classified as indels

### Changed
- **Enhanced merge logging** - Merge messages now distinguish between structural indels and homopolymer indels
  - Example: "Found mergeable subset of 2 variants: 2 SNPs, 1 homopolymer indels"
  - Provides clearer insight into why sequences are being merged
- **Improved edge case handling** - Homopolymer detection now requires strict flanking base agreement
  - Adjacent indel columns cannot validate each other as homopolymer context
  - All sequences must agree on the flanking base (not just any sequence)
  - Correctly distinguishes SNP+indel combinations from true homopolymer indels

### Documentation
- **Comprehensive README updates** - Added detailed documentation of homopolymer-aware merging
  - How homopolymer vs structural indel classification works
  - Default behavior and strict mode examples
  - Position counting examples with various parameter combinations
  - Integration with overall workflow
- **Updated CLAUDE.md** - Reflects current post-processing architecture with homopolymer-aware merging
- **Command-line reference** - Added `--disable-homopolymer-equivalence` to usage documentation

## [0.4.3] - 2025-11-05

### Changed
- **Reorganized .raw variant file locations** - Improved output directory structure for better organization
  - `.raw` FASTA files now in `variants/` subdirectory instead of main `__Summary__/` directory
  - `.raw` FASTQ files now in `variants/FASTQ Files/` instead of main `FASTQ Files/` directory
  - `summary.fasta` now contains only final consensus sequences (excludes .raw pre-merge variants)
  - Final consensus FASTA files remain in main `__Summary__/` directory
  - Final consensus FASTQ files remain in main `FASTQ Files/` directory
  - Cleaner separation between final outputs and traceability files

### Documentation
- **Updated README** - All directory structure examples and descriptions reflect new `variants/` organization
- **Updated output file descriptions** - Clear documentation of where each file type is located

## [0.4.2] - 2025-11-05

### Fixed
- **MSA alignment mode in speconsense-summarize** - Fixed SPOA to use global alignment mode (`-l 1`) instead of default local alignment
  - Corrects handling of terminal SNPs which were incorrectly treated as indels under local alignment
  - Restores consistency with previous edlib-based global alignment behavior
  - Terminal position differences (e.g., A vs G at position 1) now properly merge as IUPAC ambiguity codes (e.g., 'R')

## [0.4.1] - 2025-11-05

### Changed
- **Variant naming convention** - Primary variants now receive `.v1` suffix for consistency with additional variants
  - Previous: `sample-1`, `sample-1.v1`, `sample-1.v2`
  - New: `sample-1.v1`, `sample-1.v2`, `sample-1.v3`
  - Applies to all output files (FASTA, FASTQ) and raw files (`.v1.raw1`, `.v1.raw2`)
- **`--select-max-variants` parameter semantics** - Now includes primary variant in total count
  - Previous behavior: specified number of *additional* variants beyond primary
  - New behavior: specified *total* number of variants including primary
  - Both 0 and -1 are treated as unlimited (output all variants)
- **Regex patterns updated** - GroupField and VariantField extractors now correctly handle `.raw` suffix in new naming scheme

### Documentation
- **Updated README** - All examples now reflect new `.v1` suffix for primary variants
- **Updated CLAUDE.md** - Documented new summarization namespace format: `-1.v1`, `-1.v2`, `-2.v1`

## [0.4.0] - 2025-10-23

### Added
- **Quality assessment reporting** - Automatic generation of `quality_report.txt` for prioritized review of sequences with potential quality concerns
  - Identifies sequences with elevated variation (p50diff > 0 or p95diff > 0)
  - Highlights merged sequences with small components (RiC < 21)
  - Prioritizes issues from high to low severity for efficient triage
- **Customizable FASTA header fields** - New `--fasta-fields` option in speconsense-summarize for controlling metadata in output headers
  - Presets: `default`, `minimal`, `qc`, `full`, `id-only`
  - Custom field selection: specify individual fields (e.g., `size,ric,primers`)
  - Support for combining presets and fields
- **MSA-based variant merging** - Exhaustive subset evaluation approach for order-independent variant merging
  - Finds largest compatible subsets of variants for merging
  - Creates IUPAC consensus sequences with size-weighted majority voting
  - Ensures reproducible results regardless of input order
- **Indel support in variant merging** - New `--merge-indel-length` parameter to allow merging variants with short indels
  - Separate control for SNP count (`--merge-position-count`) and indel length limits
  - Both constraints must be satisfied for merging to proceed
- **.raw files for merged variants** - Pre-merge variant sequences preserved with original sequences and reads
  - Full traceability from final merged outputs to original clusters
  - `rawric` header field tracks contribution from each .raw component
- **Synthetic testing documentation** - Comprehensive guide for testing with speconsense-synth
  - Empirical findings on consensus quality, variant detection, and contamination scenarios
  - Quick reference for common testing patterns
  - Citations and limitations sections
- **Incremental output writing** - Performance optimization that writes output files as soon as they're ready
  - Reduces memory footprint for large datasets
  - Provides faster feedback during long processing runs

### Changed
- **Field name updates** - Standardized naming conventions for FASTA header fields
  - `merged_ric` → `rawric`
  - `median_diff` → `p50diff`
  - `p95_diff` → `p95diff`
- **Enhanced variant merging architecture** - HAC grouping now occurs before merging to prevent inappropriate merging of dissimilar sequences
  - Merging applied independently within each HAC group
  - Prevents contaminants from merging with primary targets
- **Improved logging** - Enhanced clarity in speconsense-summarize processing logs
  - Complete variant summaries including skipped variants
  - Detailed difference categorization (substitutions, single-nt indels, short indels, long indels)
  - Clear group context and selection rationale
- **Documentation reorganization** - Converted proposal documents to permanent user documentation
  - IUPAC phasing limitations and merging controls documented
  - Enhanced README with comprehensive guidance on all features

### Removed
- **Legacy code cleanup** - Removed deprecated pairwise merging approach
- **Simplified output structure** - Removed `raw_clusters/` directory from summary output (replaced by .raw files)
- **Parameter consolidation** - Removed `--output-raw-variants` parameter (raw files now generated automatically)

### Fixed
- **FASTQ consistency** - Fixed bugs in FASTQ file handling for merged variants
  - Proper aggregation of reads from all merged components
  - Consistent file naming and metadata

### Documentation
- **Updated README** - Comprehensive documentation of new features and parameter changes
- **Field name migration** - All examples and documentation updated to use new standardized field names
- **Enhanced user guides** - Detailed explanations of quality reporting, variant merging, and customization options

## [0.3.5] - 2025-10-19

### Fixed
- Fixed cluster file matching in `speconsense-summarize` to prevent false positives (e.g., cluster "c1" incorrectly matching "c10", "c11")
- Fixed specimen name parsing in FASTQ lookup table using robust regex pattern instead of fragile string splitting
- Added validation for matched FASTQ files to detect missing or empty files

### Changed
- Improved FASTQ file lookup with stage priority system (sampled → reads → untrimmed) to prevent duplicate processing
- Elevated logging level to WARNING for unexpected file conditions (empty files, unmatched patterns) for better visibility
- Enhanced file matching to use pattern `{cluster_name}-RiC` for more precise cluster identification

## [0.3.4] - 2025-10-19

### Fixed
- Fixed cluster size filtering to apply after merging identical/homopolymer-equivalent clusters, allowing small clusters with identical consensus sequences to merge before being evaluated against size thresholds
- Fixed duplicate output sequences by trimming primers before cluster merging comparison, ensuring clusters differing only in primer regions are properly merged
- Fixed sequence comparison to use standard IUPAC semantics with `end_skip_distance=0` instead of previous `end_skip_distance=20`, preventing incorrect distance calculations

### Added
- Added IUPAC-aware alignment using `additionalEqualities` in edlib for proper handling of ambiguity codes (Y, R, N, W, M, S, K, B, D, H, V)
- Added centralized `STANDARD_ADJUSTMENT_PARAMS` for consistent sequence comparison across all alignment functions

### Changed
- Standardized all alignment functions (`calculate_substitution_distance`, `calculate_adjusted_identity_distance`, `create_iupac_consensus`, `create_variant_summary`) to use IUPAC-aware alignment
- Updated variant summary message from "identical sequences (after alignment)" to "identical sequences (IUPAC-compatible)" for clarity

## [0.3.3] - 2025-09-07

### Added
- **Sequence orientation detection** - Added `--orient-mode` parameter for automatic detection and correction of sequence orientation based on primer positions
- **Three orientation modes** - `skip` (default, no orientation), `keep-all` (orient but keep failed), `filter-failed` (orient and remove failed)
- **Binary scoring system** - Simple and robust orientation detection using forward and reverse primer matches
- **Comprehensive primer handling** - Improved primer loading with position awareness (forward/reverse) and backward compatibility

### Changed
- **Enhanced primer storage** - Primers now stored separately by position for more accurate orientation detection
- **Quality score handling** - Properly reverses quality scores when sequences are reverse-complemented
- **Better logging** - Reports orientation statistics (kept as-is, reverse-complemented, failed)

### Testing
- Added comprehensive test suite for orientation feature
- Fixed existing tests to work with new default output directory structure

## [0.3.2] - 2025-08-29

### Added
- **Sampled FASTQ files** - Added `-sampled.fastq` files to cluster_debug containing only sequences used for consensus generation
- **Improved summarize process** - Summary FASTQ files now contain only consensus-contributing reads rather than entire clusters
- **Better debugging workflow** - Separate full cluster files for complete debugging and sampled files for meaningful summary output

## [0.3.1] - 2025-08-28

### Added
- **Output directory control** - Added `--output-dir` (`-O`) option to specify output directory (default: `clusters`)
- **Automatic primer detection** - Automatically searches for `primers.fasta` in input file directory when `--primers` not specified
- **Enhanced speconsense-summarize defaults** - `--source` now defaults to `clusters` to match speconsense output directory

### Fixed
- **Fixed Multiple ID counter** - Now correctly counts sequences within each specimen (resets per specimen)
- **Fixed FASTQ file path resolution** - Resolved issue where speconsense-summarize couldn't find FASTQ files when run from outside the clusters directory

### Changed
- **Enhanced output file naming** - Added `-RiC{number}` suffix to individual output files for compatibility with downstream tools
- **Improved header format** - Changed FASTA header delimiter from semicolon to space for consistency across all fields

### Documentation
- **Updated all examples** - File naming patterns, directory structures, and header formats now reflect current implementation
- **Enhanced help text** - Improved argument descriptions with default values and automatic behaviors

## [0.3.0] - 2025-08-28

### Performance
- **Dramatically improved file I/O performance** - Replaced BioPython FASTQ parsing with direct file concatenation, achieving orders of magnitude speedup
- **Eliminated directory scanning bottleneck** - Single lookup table build replaces hundreds of glob operations
- **Optimized raw file copying** - Pre-built file lookup system scales efficiently with large datasets

### Added
- **Complete variant selection framework** - Added size-based and diversity-based selection strategies with configurable limits
- **Hierarchical variant grouping** - Implemented HAC clustering to separate distinct biological sequences
- **SNP-based variant merging** - Greedy merging approach with IUPAC consensus generation
- **Comprehensive variant analysis** - Detailed logging with difference categorization (substitutions, indels by length)

### Documentation
- **Expanded clustering documentation** - Clear guidance on when to use greedy vs. graph clustering
- **Technical algorithm descriptions** - Proper terminology and algorithmic details
- **User decision framework** - Specific use cases and parameter recommendations
- **Complete workflow documentation** - Step-by-step processing pipeline with all options explained
- **Enhanced logging features** - IUPAC-aware comparisons and full traceability through processing steps
- **Parameter reference** - Comprehensive guide to all summarize options with examples
- **Example outputs** - Real log examples and file structure illustrations

### Changed
- **Removed legacy code** - Eliminated all non-optimized function versions
- **Simplified interfaces** - Clean function signatures and reduced conditional logic
- **Maintained backward compatibility** - All existing functionality preserved

## [0.2.1] - 2025-08-27

### Core Functionality
- Stable clustering and consensus generation with Markov Clustering and greedy algorithms
- Comprehensive variant merging with adjusted identity scoring
- SPOA-based consensus generation with stability assessment
- Complete output structure with traceability features
