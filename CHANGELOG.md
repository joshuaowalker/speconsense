# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
