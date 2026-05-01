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
- `significance.py`: Critical error rate (p*) via binomial survival with Bonferroni correction, uniform error model (q=p/3). CER only reported when `p* < SIGNAL_DESTRUCTION_THRESHOLD (0.75)`.
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

Read in `SpecimenClusterer.cluster()` (`speconsense/core/clusterer.py`); 12 sequential phases:

1. **Initial clustering** — MCL graph-based or greedy
2. **Pre-phasing merge** — combine HP-equivalent initial clusters
3. **Variant phasing** — split clusters by haplotype via MSA position analysis
4. **Post-phasing merge** — combine HP-equivalent subclusters
5. **Noise filter** — drop small clusters with no-majority columns
6. **Read reassignment** (optional) — concordance-based moves within identity groups
7. **Discard recovery** (optional, coupled to 6) — re-admit previously-dropped reads
8. **Second phasing pass** — re-phase any clusters that gained reads via 6/7
9. **CER validation** — annotate each non-anchor candidate with its `cer_factor`
10. **Size filtering** — drop clusters below `--min-size` and `--min-cluster-ratio`
11. **Output generation** — final consensus, MAD outlier removal, FASTA writing
12. **Discard reads written** (optional, `--collect-discards`)

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
- CER is reported only when `p* < 0.75` (signal destruction threshold: at 0.75 the implied reference population equals the variant count under the uniform model).

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
