# Architecture: Module Map and Pipeline

Reference for how a run is structured. The authoritative per-phase detail lives in
the `SpecimenClusterer.cluster()` docstring and the per-phase `_run_*` method docstrings
(`speconsense/core/clusterer.py`), and the `# Phase N:` comments in
`speconsense/summarize/cli.py`; this file is the map that tells you which phase to go read.

## Module map

**`speconsense/core/`** — clustering and consensus generation:
- `clusterer.py`: `SpecimenClusterer` — orchestrates the pipeline
- `workers.py`: parallel workers (SPOA, cluster processing, phasing, primer trimming,
  Phase 6 group reassignment)
- `records.py`: `ReadRecord` — lightweight read container (id/title/seq/qual strings)
  replacing BioPython SeqRecords in core; fast 4-line FASTQ parse/write, byte-identical
  output to SeqIO
- `cli.py`: argument parsing and input loading (`records.parse_input`)

**`speconsense/summarize/`** — post-processing:
- `cli.py`: entry point and `process_single_specimen`
- `iupac.py`: IUPAC-aware distances (re-exports shared helpers from `distances.py`)
- `fields.py`: FASTA header field classes and formatting
- `analysis.py`: MSA analysis, cluster quality, outlier detection
- `merging.py`: MSA-based variant merging with IUPAC consensus
- `clustering.py`: bucketing by core-assigned gid, cross-primer anchor merger, variant selection
- `full_consensus.py`: the `-{gid}-full` group consensus builder
- `io.py`: loading sequences, writing consensus FASTA/FASTQ/debug outputs

**Shared top-level modules** (used by both subpackages):
- `types.py`: `ConsensusInfo`, `OverlapMergeInfo` NamedTuples — exists to avoid circular imports
- `msa.py`: SPOA MSA analysis, homopolymer-normalized error detection, IUPAC generation,
  variant position detection/phasing support. Defines `IUPAC_CODES`, `DEFAULT_MAX_ERR_FACTOR`
- `distances.py`: IUPAC-aware edlib alignment, adjusted-identity distance, variant difference
  counting. Defines `IUPAC_EQUIV`, `STANDARD_ADJUSTMENT_PARAMS`
- `context.py`: per-position variant context classification (`ContextClass`, `ContextTag`)
  driving CER q_ctx lookup
- `qctx.py`: loads error models from `error_models/*.yaml`. `DEFAULT_MODEL_NAME = "dorado-v5.0"`,
  `MAX_HP_LENGTH = 5`. Resolution order: filesystem path → `~/.config/speconsense/error_models/` → bundled
- `significance.py`: critical error rate (p*) via binomial survival with Bonferroni correction.
  Defines `DEFAULT_MIN_CER_FACTOR`
- `outliers.py`: `detect_rid_outliers` — MAD outlier removal
- `quality_report.py`: multi-section quality report for summarize
- `cli.py`: top-level entry-point stub re-exporting `core.main`

**`speconsense/scalability/`** — optional acceleration for O(n²) pairwise work:
- `base.py`: `CandidateFinder` protocol, `ScalablePairwiseOperation`
- `vsearch.py`: vsearch-backed candidate finder; `config.py`: `ScalabilityConfig`
- Two activation regimes, both requiring vsearch on PATH and `--scale-threshold > 0`:
  **read-level** operations (initial-clustering K-NN for the MCL graph) gate on
  `len(seqs) >= scale_threshold` (default 1001, calibrated for read counts);
  **cluster-level** operations pass `force_scalable=True` and activate at `n > 50` items
  regardless of the read-count threshold — cluster-equivalence merging (HP-equivalence
  union-find, phases 2/4), identity grouping (`_form_identity_groups` sparse distance
  matrix — used by read reassignment, discard recovery, and CER validation), discard
  screening (top-K cluster matches per discard), CER validation within-group top-K (when
  a group exceeds 50). `--scale-threshold 0` disables all of it.
- Independent of vsearch, `--threads N>1` parallelizes SPOA via `ProcessPoolExecutor`:
  noise-filter consensus (Phase 5), validation-consensus batches (phases 6/7,
  `_generate_validation_consensuses`), per-group read reassignment (Phase 6,
  `_reassign_group_worker`), and the second phasing pass (Phase 8).

**`speconsense/profiles/`** — parameter presets. See `docs/profile-parameters.md` for the
user-facing reference. Implementation notes:
- Bundled YAML profiles: `compressed`, `herbarium`, `largedata`, `nostalgia`, `strict`, `example`.
  User profiles in `~/.config/speconsense/profiles/` take precedence.
- Override order: defaults → profile (explicit or auto-detected) → explicit CLI arguments
- `VALID_SPECONSENSE_KEYS` / `VALID_SUMMARIZE_KEYS` in `profiles/__init__.py` are the source of
  truth for acceptable keys; keys are strictly validated
- Reserved name `default` (`RESERVED_PROFILE_NAMES`) selects CLI defaults without loading a YAML
  file — use `-p default` on summarize to override auto-detection, or on core to record an
  explicit "no profile" in metadata
- **Auto-detection**: summarize scans metadata JSONs for the `profile` field written by core.
  If all specimens reporting a profile agree, it is auto-applied (with a logged parameter listing);
  mixed profiles warn. Specimens without metadata or without the field are ignored. Any `-p`
  (including `-p default`) skips auto-detection.

**`speconsense/error_models/`** — bundled per-basecaller error models (YAML), loadable by name
(`--error-model dorado-v5.0`), from `~/.config/speconsense/error_models/`, or by filesystem path.

**`speconsense/synth.py`** — the `speconsense-synth` CLI: synthetic read generator for testing
consensus algorithms. See `docs/synthetic-testing.md`.

**`speconsense/fit_error_model/`** — the `speconsense-fit-error-model` CLI. Offline q_ctx
re-estimation from a finished output tree;
writes a user model to `~/.config/speconsense/error_models/{name}.yaml`. Implements approach-1
(HP rates from primary-anchor MSAs, mode-as-ground-truth + Bonferroni outlier + bimodal filter)
and approach-2 (non-HP rates pooled across all-cluster MSAs) from HP paper §8 / CER paper §4.2.
`hp-l{N}` is the implied per-position rate `p` such that `(1-p)^N = frac_correct`
(HP paper §3.5) — **not** raw `1-frac_correct`.

## Core pipeline (`SpecimenClusterer.cluster()`)

Fourteen sequential phases. Each has a `_run_*` (or `_filter_*` / `_write_*`) method with a
docstring covering its contract; the `cluster()` docstring lists them all.

1. **Initial clustering** — MCL graph-based or greedy
2. **Pre-phasing merge** — combine HP-equivalent initial clusters
3. **Variant phasing** — split clusters by haplotype via MSA position analysis
4. **Post-phasing merge** — combine HP-equivalent subclusters
5. **Noise filter** — drop small clusters with no-majority columns
6. **Read reassignment** (optional) — concordance-based moves within identity groups
7. **Discard recovery** (optional, coupled to 6) — re-admit previously-dropped reads
8. **Second phasing pass** — re-phase clusters that gained reads via 6/7. Gated on *three*
   flags: any of `--disable-position-phasing`, `--disable-read-reassignment`, or
   `--disable-second-phasing` skips it. Phase 7 is coupled to Phase 6, so if 6 did not run,
   nothing has changed since Phase 3 and there is no work for 8 to do
   (see `clusterer.py` around the Phase 8 call site)
9. **Cluster consensus generation** — SPOA → MAD outlier removal → re-SPOA → IUPAC ambiguity
   calling → primer trimming; stamps the post-MAD MSA and consensus on each `cluster_dict`
10. **Post-refinement merge** — combine clusters whose post-MAD consensuses are identical or
    HP-equivalent; reruns the Phase 9 worker on each merge survivor
11. **CER validation** — annotate each non-anchor candidate with its `cer_factor`
12. **Size filtering** — drop clusters below `--min-size`
13. **Output emission** — compute `err_factor` on each stamped MSA, assign identity ranks
    (`_assign_identity_ranks`), write FASTA/FASTQ/MSA
14. **Discard reads written** (optional, `--collect-discards`)

Reads that fail clustering or are dropped by any filter accumulate in
`self.discarded_read_ids`; phase 7 attempts to recover the concordant ones.

### Orientation and primer trimming

These run outside the numbered phases — orientation during input processing, primer trimming
during final consensus. Two backends:

- **`primer`** — edlib-based primer matching at read ends
- **`pyitsx`** — HMM profile-based ITS strand detection via pyitsx (optional dependency, ITS loci
  only; `--pyitsx-organism` selects the kingdom, default Fungi)

Both discard failed reads. `--orient-mode=pyitsx+primer` runs pyitsx first, then falls back to
primer-based orientation for reads it cannot classify (`strand=None`); chimeric reads are always
discarded. This matters for mixed runs where some amplicons target non-ITS loci not covered by
the ITSx HMMs. Plain `--orient-mode=pyitsx` discards all pyitsx failures with no fallback.

## Summarize pipeline (`process_single_specimen`)

1. Load sequences with RiC filtering
2. Bucket by core-assigned `gid=` (identity group rank) — **no re-grouping**
3. Cross-primer overlap conflation (`--min-merge-overlap`) merges different-primer core groups
   whose members overlap well (primer-pool use case: ITS/ITS1/ITS2)
4. Homopolymer-aware MSA-based merging within each (possibly-conflated) group
5. Post-merge re-check — cancel merges whose recomputed `err_factor`/`cer_factor` cross filter
   thresholds; restore the original contributors on the pass track
6. Selection size ratio filtering (`--select-min-size-ratio`) — removes tiny post-merge variants
7. Variant selection within each group (size-based or diversity-based)
8. Output generation with full traceability

Identity grouping is performed **once, in core**, via complete linkage on `--group-identity`
(default 0.85). Summarize honors those groups verbatim (and hard-fails on inputs lacking
`gid=`/`vid=`); it only merges *across* core groups when cross-primer overlap passes threshold.

## Algorithm selection

**Graph-based MCL (default)** — better at detecting sequence variants within a specimen; more
clusters, higher granularity; requires MCL on PATH (falls back to greedy if missing).

**Greedy (`--algorithm greedy`)** — faster, fewer clusters, focuses on well-separated sequences;
good for distinguishing distinct targets from contaminants.

## External dependencies

- **SPOA** (required, must be on PATH) — consensus generation
- **MCL** (optional, recommended) — graph clustering
- **vsearch** (optional) — enables the `scalability` candidate-finder backend
- **edlib** — edit distances, used with `IUPAC_EQUIV` for ambiguity-aware alignment
- **adjusted-identity** (`>=0.2.4`, from GitHub) — IUPAC-aware sequence alignment
- **BioPython**, **NumPy**, **SciPy**, **tqdm**, **PyYAML** — see `pyproject.toml`

## Integration context

Designed to replace NGSpeciesID in the ONT fungal barcoding pipeline from protocols.io.
Processes demultiplexed FASTQ (typically downstream of specimux) and generates consensus
sequences suitable for taxonomic identification. See
`~/mm/code/specimux-suite/INTEGRATION.md` for the profile system and subprocess invocation
contracts with specimux-suite.
