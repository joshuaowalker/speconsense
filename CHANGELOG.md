# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
