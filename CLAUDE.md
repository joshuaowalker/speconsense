# CLAUDE.md

Guidance for Claude Code working in this repository.

## Project

Speconsense clusters Oxford Nanopore amplicon reads and generates high-quality consensus
sequences — an experimental alternative to NGSpeciesID in the fungal DNA barcoding pipeline.
Two entry points do the work: `speconsense` (clustering + consensus) and
`speconsense-summarize` (post-processing, quality routing, naming).

## Commands

```bash
pip install -e .                  # dev install
pytest                            # full suite (config in pytest.ini)
pytest tests/test_orientation.py  # single file

conda install bioconda::spoa      # required, must be on PATH
conda install bioconda::mcl       # optional, recommended
```

README.md covers installation, usage, options, and output layout for users. `--help` is
current for every flag; prefer it over any list written down here.

## Where things are

| Need | Read |
|---|---|
| Module map, the 14 core phases, summarize's 8, orientation backends | `docs/architecture/pipeline.md` |
| CER, `err_factor`, identity grouping, gid/vid ranks, IUPAC/distance machinery | `docs/architecture/variant-model.md` |
| Naming policy, four output tracks, merge-time field recompute, `-full` | `docs/architecture/summarize-outputs.md` |
| Profile parameters / FASTA header fields / RiC semantics (user-facing) | `docs/profile-parameters.md`, `docs/customizing-fasta-headers.md`, `docs/understanding-ric-and-merging.md` |

Phase-level detail lives in code: the `SpecimenClusterer.cluster()` docstring and the per-phase
`_run_*` method docstrings (`speconsense/core/clusterer.py`) and the `# Phase N:` comments in
`speconsense/summarize/cli.py`.
Those are authoritative — the docs above are the map, not the territory.

## Gotchas

Things that will bite you, that reading the nearby code will not reveal.

**The MCL K-NN graph is deliberately asymmetric.** `similarities[id1]` holds only neighbors
`id2 > id1` (`scalability/base.py`). This contradicts the MCL docs and looks like a bug; it is a
weakly-held decision that affects tie-breaking and changes clustering results if "fixed."
Validate against existing tests before touching it.

**`DEFAULT_MIN_CER_FACTOR` (`significance.py`) and `DEFAULT_MAX_ERR_FACTOR` (`msa.py`) are shared
on purpose.** They back both core's vid tiering and summarize's `--min-cer-factor` /
`--max-err-factor` defaults, so the two cannot drift. Don't hardcode `1.0` / `1.5` anywhere.

**Core stamps `gid`/`vid` once; summarize never renumbers.** Gaps in the vid sequence are expected
output — a variant routed to `.ns`/`.lq`/`.filtered` leaves its number behind. Summarize hard-fails
on inputs lacking `gid=`/`vid=`; identity grouping happens only in core.

**The largest cluster is not necessarily `1.v1`.** Variant ranking is
`(expected_to_pass_summarize, size)`, so a high-`err_factor` cluster gets demoted. This is why
`fit_error_model` re-derives the primary anchor by read count instead of trusting the label.

**`qctx.MAX_HP_LENGTH = 5` and `--hp-normalization-length` (default 6) are different knobs.** The
first is where the error-model table stops and blanket HP normalization takes over; the second is
the distance-comparison threshold shared by core and summarize.

**`snp_count` can over-count by 0–2** across iterative merge rounds when the same physical position
goes ambiguous more than once. `ambig` is the canonical IUPAC-site count.

**SPOA output depends on input order.** The first sequence anchors the graph and therefore the
alignment coordinate frame. `_run_spoa_for_cluster_worker` writes dict entries in insertion order,
so callers control the frame by construction order. Corollary: never build SPOA input by iterating
a set — hash-seed order leaks into the coordinate frame and makes results vary run-to-run (this
was a real bug, fixed in 0.8.6 by sorting on `(-mean_quality, read_id)`). Sort explicitly. All cluster-level SPOA calls use linear gap
scoring (`-m 1 -n -1 -g -1 -e -1`), which produces roughly 4x fewer ambiguities than SPOA's
defaults — keep them in sync across `core/workers.py` and `summarize/analysis.py`.

**`-full` records carry IUPAC codes.** They are meant for BLAST against legacy unphased ITS
references *with adjusted-identity scoring*, and degrade silently under raw BLAST.

**`ConsensusInfo` / `OverlapMergeInfo` are NamedTuples** (`types.py`, structured that way to avoid
circular imports). Use `._replace()` to derive modified copies.

**Profile keys use dashes, argparse attrs use underscores** (`disable-merging` →
`disable_merging`). `VALID_SPECONSENSE_KEYS` / `VALID_SUMMARIZE_KEYS` in `profiles/__init__.py`
are the source of truth, and validation is strict.

**In `summarize/fields.py`, the `.raw` portion of the variant regex must stay outside the capture
group**, or `VariantField` returns `v1.raw2` instead of `v1`.

## Conventions

Match the surrounding code — this codebase leans on substantial module and function docstrings
that explain *why*, especially in `clusterer.py`. Follow that where the reasoning is non-obvious.

Behavior changes that affect FASTA headers, naming, or output tracks need a test; `tests/` mixes
direct-import unit tests with subprocess-driven integration tests of the CLIs.

## Related

`~/mm/code/specimux-suite/INTEGRATION.md` — profile system and subprocess invocation contracts
with specimux-suite.
