# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Speconsense is a Python tool for high-quality clustering and consensus generation from Oxford Nanopore amplicon reads. It's designed as an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.

## Key Commands

### Running the tools
```bash
# Main clustering tool (basic usage)
speconsense input.fastq

# With common options
speconsense input.fastq --algorithm greedy --min-size 10 --primers primers.fasta

# Post-processing tool
speconsense-summarize --min-ric 5 --source /path/to/output --summary-dir MyResults

# Synthetic data generator for testing
speconsense-synth reference.fasta --num-reads 1000 --error-rate 0.05 --output synthetic.fastq
```

### Development setup
```bash
# Install in development mode
pip install -e .

# External dependencies (must be in PATH)
conda install bioconda::spoa  # Required
conda install bioconda::mcl   # Optional but recommended
```

### Testing
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_orientation.py

# Run tests with specific markers
pytest -m "not slow"
```

## Code Architecture

### Main Components

**speconsense/core/** - Core clustering and consensus generation subpackage:
- `clusterer.py`: `SpecimenClusterer` class — main orchestrator for the clustering pipeline
- `workers.py`: Worker functions for parallel processing (SPOA, cluster processing, phasing, primer trimming)
- `cli.py`: Command-line interface and argument parsing
- Two clustering algorithms: graph-based MCL (default) and greedy clustering
- Consensus generation via external SPOA tool
- Cluster merging based on consensus sequence similarity
- Stability assessment through subsampling

**speconsense/summarize/** - Post-processing utility subpackage:
- `cli.py`: Command-line interface and main entry point
- `iupac.py`: IUPAC-aware distance calculations (re-exports shared helpers from `distances.py`)
- `fields.py`: FASTA header field classes and formatting
- `analysis.py`: MSA analysis, cluster quality, outlier detection
- `merging.py`: MSA-based variant merging with IUPAC consensus
- `clustering.py`: Bucketing by core-assigned gid, cross-primer anchor merger, variant selection
- `io.py`: File I/O (loading sequences, writing consensus FASTA/FASTQ/debug outputs)

**Shared top-level modules** (used by both subpackages):
- `types.py`: `ConsensusInfo`, `OverlapMergeInfo` NamedTuples — avoids circular imports. Use `._replace()` to create modified copies.
- `msa.py`: SPOA MSA analysis, homopolymer-normalized error detection, IUPAC generation, variant position detection/phasing support. Defines `IUPAC_CODES`.
- `distances.py`: IUPAC-aware edlib alignment, adjusted-identity distance, variant difference counting. Defines `IUPAC_EQUIV` and `STANDARD_ADJUSTMENT_PARAMS`.
- `context.py`: Per-position variant context classification (`ContextClass`, `ContextTag`) driving CER q_ctx lookup. HP context comes from the reference consensus.
- `qctx.py`: Loads error models (context-specific error-rate tables) from `error_models/*.yaml`. `DEFAULT_MODEL_NAME = "dorado-v5.0"`. HP runs beyond `MAX_HP_LENGTH=5` route to blanket normalization. Resolution order: filesystem path → `~/.config/speconsense/error_models/` → bundled.
- `significance.py`: Critical error rate (p*) via binomial survival with Bonferroni correction, uniform error model (q=p/3). The shipped pipeline reports `cer_factor` (per-position multiplicative inflation) directly from the joint q* solver; the uniform-p* form is retained for the variant-significance paper's tables.
- `quality_report.py`: Multi-section quality report for `speconsense-summarize`.
- `cli.py`: Top-level entry-point stub that re-exports `core.main`.

**speconsense/scalability/** - Optional acceleration for O(n²) pairwise work:
- `base.py`: `CandidateFinder` protocol, `ScalablePairwiseOperation` (top-K neighbors for MCL, sparse distance matrix for identity grouping, equivalence groups for cluster merging)
- `vsearch.py`: vsearch-backed candidate finder
- `config.py`: `ScalabilityConfig`
- Activates when `len(seqs) >= scale_threshold` AND vsearch is on PATH. Active uses: initial-clustering K-NN (MCL graph), cluster-equivalence merging (HP-equivalence union-find), identity grouping (`_form_identity_groups` sparse distance matrix, used by read reassignment, discard recovery, and CER validation), discard screening (top-K cluster matches per discard), CER validation within-group top-K (per-group most-similar peers when group exceeds 50). Noise-filter SPOA is parallelized via `ProcessPoolExecutor` when `--threads N>1` (no vsearch dependency).

**speconsense/profiles/** - Profile system for parameter presets:
- YAML profiles (`compressed`, `herbarium`, `largedata`, `nostalgia`, `strict`, `example`) bundled in the package
- User profiles in `~/.config/speconsense/profiles/` take precedence over bundled
- Override order: defaults → profile → explicit CLI arguments
- Valid keys are strictly validated; profile keys use dashes (e.g., `disable-merging`), argparse attrs use underscores
- `VALID_SPECONSENSE_KEYS` / `VALID_SUMMARIZE_KEYS` in `profiles/__init__.py` are the source of truth for acceptable keys

**speconsense/error_models/** - Bundled per-basecaller error models (YAML), loadable by name (`--error-model dorado-v5.0`), from `~/.config/speconsense/error_models/`, or by filesystem path.

**speconsense/synth.py** - Synthetic read generator for testing consensus algorithms.

### Key Processing Pipeline

Read in `SpecimenClusterer.cluster()` (`speconsense/core/clusterer.py`); 14 sequential phases:

1. **Initial clustering** — MCL graph-based or greedy
2. **Pre-phasing merge** — combine HP-equivalent initial clusters
3. **Variant phasing** — split clusters by haplotype via MSA position analysis
4. **Post-phasing merge** — combine HP-equivalent subclusters
5. **Noise filter** — drop small clusters with no-majority columns
6. **Read reassignment** (optional) — concordance-based moves within identity groups
7. **Discard recovery** (optional, coupled to 6) — re-admit previously-dropped reads
8. **Second phasing pass** — re-phase any clusters that gained reads via 6/7
9. **Cluster consensus generation** — SPOA → MAD outlier removal → re-SPOA → IUPAC ambiguity calling → primer trimming; stamps post-MAD MSA and consensus on each cluster_dict
10. **Post-refinement merge** — combine clusters whose post-MAD consensuses are identical or HP-equivalent; reruns Phase 9 worker on each merge survivor
11. **CER validation** — annotate each non-anchor candidate with its `cer_factor`
12. **Size filtering** — drop clusters below `--min-size` and `--min-cluster-ratio`
13. **Output emission** — write FASTA/FASTQ/MSA; compute `err_factor` on stamped MSA
14. **Discard reads written** (optional, `--collect-discards`)

Orientation (when `--orient-mode` ≠ skip) and primer trimming run during input processing and final consensus respectively, outside the numbered phases. Reads that fail clustering or are dropped by any filter accumulate in `self.discarded_read_ids`; phase 7 attempts to recover concordant ones.

### Post-Processing Pipeline (speconsense-summarize)

1. Load sequences with RiC filtering
2. Bucket by core-assigned `gid=` (identity group rank) — no re-grouping
3. Cross-primer overlap conflation (`--min-merge-overlap`) merges different-primer core groups whose members overlap well (primer-pool use case: ITS/ITS1/ITS2)
4. Homopolymer-aware MSA-based merging within each (possibly-conflated) group
5. Selection size ratio filtering (`--select-min-size-ratio`) — removes tiny post-merge variants
6. Variant selection within each group (size-based or diversity-based)
7. Output generation with full traceability

Note: Identity grouping is performed once, in core, via complete linkage on `--group-identity` (default 0.85). Summarize honors those groups verbatim (hard-fails on inputs lacking `gid=`/`vid=`) and only merges *across* core groups when cross-primer overlap passes threshold.

**Naming policy** (`process_single_specimen`): summarize preserves core's `gid`/`vid` on every record that is **not moved between groups by cross-primer conflation**. Variants dropped by within-group MSA merging, `--select-max-variants`, `--select-min-size-ratio`, `--min-cer-factor`, or `--max-err-factor` leave their vids as gaps. Variants moved into a survivor group by cross-primer conflation adopt the survivor's `gid` and get a freshly-minted `vid` strictly above the highest vid core ever wrote under that gid — collision-avoidance considers passed + `.ns` + `.lq` records under that gid plus any vids already minted in the same pass (covers 3+ group conflation). Vids under absorbed-group gids stay on disk under their original gid in `.ns`/`.lq` outputs and do not block the survivor's namespace.

**Merge-time field handling**: when ≥2 contributors collapse into one record (within-group MSA merge or cross-primer overlap conflation), summarize aggregates per-field via `_build_merged_consensus_info` and a follow-up recompute pass in `process_single_specimen`. `size`/`ric` are summed; `length`/`ambig` are re-derived from the column-voted consensus; `rawric`/`rawlen` carry flattened merge history; `primers` is the union of contributor primer names sorted; `rid`/`rid_min`/`err_factor` are re-derived from a SPOA MSA over the union of contributor reads (loaded from each contributor's `cluster_debug` FASTQ); `cer_factor` is recomputed against larger peers in the post-merge bucket for same-primer merges and set to `None` for cross-primer merges (CER noise model is per-locus). `snp_count` is cumulative across iterative rounds and can over-count by 0–2 when the same physical position becomes ambiguous in multiple rounds — `ambig` is the canonical IUPAC-site count. `.raw` pre-merge files carry every field from their source cluster (matching `.ns`/`.lq` pass-through), only resetting the merge-history fields. Recompute requires the contributors' debug FASTQs and the per-specimen metadata JSON's `parameters.error_model`; missing inputs fall back to the inherited values from the largest contributor.

**Group full consensus** (`--enable-full-consensus`, `speconsense/summarize/full_consensus.py`): per identity group with ≥2 selected variants on the pass track, summarize emits an additional `-{gid}-full` record built from a size-weighted, top-mean-Phred sample of the **pre-merge core variants** in the bucket (sourced from `variant_groups[final_gid]`, so within-group MSA merging and selection filters do not constrain the read pool). Variants are gated by a running-total threshold using `parameters.min_ambiguity_frequency` from core's metadata JSON (default 0.10): each contributor must satisfy `size ≥ min_ambiguity_frequency × running_total`. The sampled reads (budget = `parameters.max_sample_size`, default 100) go through SPOA with linear gap scoring; `-l 0` (local SW) when the gated bucket spans multiple primer sets, else `-l 1` (global NW). The MSA is collapsed via `build_full_consensus_from_msa` with one-vote-per-read and majority-wins gaps; columns where ≥2 bases each clear `min_ambiguity_frequency` get IUPAC codes. Output: FASTA entry in `-all.fasta`, sampled-reads FASTQ in `FASTQ Files/`, no `.raw` lineage (synthetic record). Suppressed when `<2` pass-track variants exist for the group, `<2` contributors clear the gate, or the read pool can't be loaded. `.ns` and `.lq` records are excluded from the pool by construction (filtered upstream by `load_consensus_sequences`). Intended use is BLAST query against legacy unphased ITS references with adjusted-identity scoring — `-full` is IUPAC-bearing and will silently degrade under raw BLAST.

### External Dependencies

- **SPOA**: Required for consensus generation, must be in PATH. When running SPOA, the candidate sequence must be the first input.
- **MCL**: Optional but recommended for graph clustering (falls back to greedy if missing)
- **vsearch**: Optional; enables the `scalability` candidate-finder backend for large inputs
- **edlib**: Edit distance calculations; used with `IUPAC_EQUIV` for ambiguity-aware alignment
- **adjusted-identity**: IUPAC-aware sequence alignment (from GitHub, `>=0.2.4`)
- **BioPython**, **NumPy**, **SciPy**, **tqdm**, **PyYAML**: see `pyproject.toml`

### Output Structure

**Main consensus files** (from speconsense):
- `{sample_name}-all.fasta` - All consensus sequences for a specimen
- `cluster_debug/` directory:
  - `{sample_name}-c{num}-RiC{size}-reads.fastq` - Original reads
  - `{sample_name}-c{num}-RiC{size}-sampled.fastq` - Sampled reads for consensus
  - `{sample_name}-c{num}-RiC{size}-untrimmed.fasta` - Untrimmed consensus

**Summary files** (from speconsense-summarize):
- `__Summary__/` directory with processed outputs
- Dual namespace system:
  - Original clustering: `-c1`, `-c2`, `-c3` format
  - Summarization: `-1.v1`, `-1.v2`, `-2.v1` format for groups and variants

### IUPAC Ambiguity Code Handling

The codebase uses IUPAC nucleotide ambiguity codes throughout:
- `IUPAC_CODES` (in `msa.py`) maps nucleotide sets to codes (R=A/G, Y=C/T, etc.)
- `IUPAC_EQUIV` (in `distances.py`) enables edlib alignment to treat ambiguity codes as matching their constituent bases
- `STANDARD_ADJUSTMENT_PARAMS` (in `distances.py`) defines consistent sequence comparison parameters:
  - Homopolymer normalization enabled (treats "AAA" = "AAAAA")
  - IUPAC overlap disabled (uses standard IUPAC semantics: Y≠M)
  - No end trimming (`end_skip_distance=0`)
  - Single-base repeats for homopolymer normalization
- Safeguards `MIN_COVERAGE_THRESHOLD=0.5` and `MAX_ADJUSTMENT_RATIO=1.5` fall back to raw edlib identity when terminal-gap exclusion would inflate adjusted identity on length-mismatched sequences.

### Variant Significance and CER

Phasing uses a three-stage architecture: phase indiscriminately, group by identity, then annotate pairwise via CER. Identity grouping uses **complete linkage** (every pair in a group must meet `--group-identity`, default `0.85`) via `scipy.cluster.hierarchy`, preventing transitive chains that would otherwise collapse closely related variants in eDNA-style mixtures. The same groups gate read reassignment and discard recovery — reads can only move within their identity group. Key pieces:
- `significance.compute_critical_error_rate(N, M, L, alpha, K)` — p* under uniform model (q=p/3), with combinatorial Bonferroni for `K>1` multi-position variants.
- `context.classify_variant_context()` produces one `ContextTag` per variant event (substitution or contiguous indel block). HP context comes from the reference consensus — the artifact hypothesis under test is that the candidate's reads are miscalled copies of the reference.
- `qctx.get_qctx(tag, table)` returns a per-position error rate; HP runs longer than the table's max route to blanket homopolymer normalization.

**No gating in core.** Every non-anchor candidate is pairwise-compared against all larger peers in its identity group, annotated with a `cer_factor` (per-position multiplicative inflation; worst-case across peers), and flows through to the FASTA output. The reference pool accumulates all processed clusters regardless of factor — `min_factor` is inherently conservative for artifact-vs-artifact cases. **Summarize applies the user-visible pass/ns decision** via `--min-cer-factor` (default `1.0`; `0` disables). Records below the threshold are routed to `__Summary__/variants/{name}.ns-RiC{ric}.fasta` with matching FASTQ, mirroring the `.raw` layout. Records with `cer_factor=None` (anchors, failed pairwise comparison) always pass. The FASTA header carries only `cer_factor=`; full per-position detail (p*, K, context tags, q_ctx values, reference idx) lives in the metadata JSON via `_build_variant_record`.

### Cluster Homogeneity (err_factor)

Complementary to CER, `err_factor` is a cluster-wide observed-vs-q_ctx-expected disagreement ratio: for each non-gap consensus column, the fraction of reads disagreeing with the consensus divided by the q_ctx rate predicted for that column's context (HP run length or non-HP). Values near 1.0 indicate reads consistent with basecaller noise; values ≫ 1.0 indicate internal heterogeneity beyond what sequencing noise would produce. Unlike `cer_factor`, `err_factor` is peer-independent — it distinguishes novel-but-real sequences (low) from noise combinations (high).

Computed in `msa.compute_cluster_err_factor(msa_string, qctx_table)` during final consensus generation; emitted as `err_factor=` in the FASTA header and stored with raw `obs_sum`/`exp_sum`/`cols` in metadata JSON for reproducibility. **Summarize filters on it** via `--max-err-factor` (default `1.5`; `0` disables). Records above threshold route to `__Summary__/variants/{name}.lq-RiC{ric}.fasta`; `.lq` takes precedence over `.ns` when both fire. The 1.5 default is safe because MAD outlier removal at final consensus (see ``speconsense.outliers.detect_rid_outliers``) removes single-read outliers that would otherwise inflate err_factor on real clusters.

### Algorithm Selection

**Graph-based MCL (default)**:
- Better at detecting sequence variants within specimens
- More clusters, higher granularity
- Requires MCL tool in PATH (falls back to greedy if missing)

**Greedy clustering (`--algorithm greedy`)**:
- Faster, simpler output
- Fewer clusters, focuses on well-separated sequences
- Good for detecting distinct targets vs. contaminants

### MCL Graph Construction (Design Decision)

The K-NN similarity graph for MCL (in `scalability/base.py`) uses an **asymmetric edge storage pattern** where `similarities[id1]` only contains entries for neighbors `id2 > id1` (lexicographically). This is a weakly-held design decision that produces good clustering results despite the MCL documentation recommending symmetric graphs.

**Why asymmetric?** The pattern emerged from the original implementation and affects tie-breaking when multiple neighbors have identical similarity scores. Changing to symmetric storage would alter which neighbors are selected during K-NN construction, potentially changing clustering results. Validate against existing test cases before changing.

### Integration Context

Designed to replace NGSpeciesID in the ONT fungal barcoding pipeline from protocols.io. Processes demultiplexed FASTQ files and generates consensus sequences suitable for taxonomic identification. Typically used downstream of specimux for demultiplexing.

### Configuration

Parameters are controlled via CLI arguments, optionally pre-set via YAML profiles (`-p/--profile`, `--list-profiles`). Key parameters:
- Identity thresholds (`--min-identity`)
- Clustering algorithm choice (`--algorithm`)
- Sample size limits (`--max-sample-size`, `--presample`)
- Cluster size filtering (`--min-size`, `--min-cluster-ratio`)
- Primer handling (`--primers`, `--orient-mode`)
- Variant phasing (`--disable-position-phasing`, `--min-variant-frequency`, `--significance-level`)
- Post-phasing refinement: `--disable-read-reassignment` (concordance-based reassignment within identity groups), `--disable-discard-recovery` (re-admit dropped reads; auto-skipped if read reassignment is disabled). The second phasing pass is gated by `--disable-position-phasing` AND `--disable-read-reassignment` — disabling either skips it.
- Error model selection (`--error-model`, `--hp-normalization-length`)
- Summarize CER filter (`--min-cer-factor`, default `1.0`, `0` disables)
- Summarize err_factor filter (`--max-err-factor`, default `1.5`; `0` disables)
- Summarize HP threshold (`--hp-normalization-length`, default `6`; matches core; `1` restores legacy blanket-normalize-all behavior)

## Integration with specimux-suite

See `~/mm/code/specimux-suite/INTEGRATION.md` for integration contracts with specimux-suite, including:
- Profile system (`-p/--profile`, `--list-profiles`) — profiles in `~/.config/speconsense/profiles/` and bundled
- Subprocess invocation contract — how specimux-suite constructs speconsense commands
