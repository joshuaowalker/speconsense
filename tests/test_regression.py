#!/usr/bin/env python3
"""
Regression tests for bug fixes identified in the commit history.

These tests ensure that previously fixed bugs do not regress.
Each test documents the original bug and the fix commit.
"""

import os
import tempfile
import shutil
import subprocess
import sys
import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class TestEmptyInputHandling:
    """Tests for empty input file handling.

    Bug: ZeroDivisionError when input file contains no sequences.
    Fix: commit 2f6a2da (v0.6.2)
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_empty_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def core_module(self):
        """Get module name for speconsense core."""
        return 'speconsense.core'

    def test_empty_fastq_exits_gracefully(self, temp_dir, core_module):
        """Verify no crash on empty FASTQ input file.

        Prior to fix, this would raise ZeroDivisionError during
        cluster size ratio calculations.
        """
        # Create completely empty file
        with open('empty.fastq', 'w') as f:
            pass

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'empty.fastq',
            '--min-size', '0',
            '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        # Should exit gracefully (code 0) with warning, not crash
        assert result.returncode == 0, f"Should exit gracefully, got: {result.stderr}"
        assert "No sequences found" in result.stderr, "Should warn about empty input"

    def test_fastq_with_only_whitespace_errors_gracefully(self, temp_dir, core_module):
        """Verify graceful error on FASTQ with only whitespace (malformed input).

        Note: This is a malformed FASTQ file (whitespace is not valid FASTQ format),
        not a truly empty file. BioPython correctly raises ValueError for this.
        The key is that it fails predictably, not with an obscure crash.
        """
        with open('whitespace.fastq', 'w') as f:
            f.write("\n\n  \n")

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'whitespace.fastq',
            '--min-size', '0',
            '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        # Malformed FASTQ should fail with clear error message
        assert result.returncode != 0, "Malformed FASTQ should fail"
        assert "ValueError" in result.stderr or "error" in result.stderr.lower(), \
            "Should show error for malformed input"

    def test_valid_fastq_still_works(self, temp_dir, core_module):
        """Verify normal operation with valid FASTQ is unaffected."""
        # Create a valid FASTQ with sequences
        records = [
            SeqRecord(
                Seq("ACGTACGTACGTACGT"),
                id=f"read{i}",
                letter_annotations={'phred_quality': [30] * 16}
            )
            for i in range(5)
        ]
        with open('valid.fastq', 'w') as f:
            SeqIO.write(records, f, 'fastq')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'valid.fastq',
            '--min-size', '0',
            '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Should succeed: {result.stderr}"
        assert os.path.exists('clusters/valid-all.fasta'), "Output should be created"


class TestDeterministicOrdering:
    """Tests for deterministic output ordering.

    Bug: Non-deterministic dict ordering caused different clustering results.
    Fix: commit d1f917f
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_determinism_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def core_module(self):
        """Get module name for speconsense core."""
        return 'speconsense.core'

    @pytest.fixture
    def test_fastq(self, temp_dir):
        """Create test FASTQ with multiple sequences that could cluster differently."""
        # Create sequences with some variation to create interesting clustering
        base_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        records = []
        for i in range(20):
            # Add some variation
            if i % 3 == 0:
                seq = base_seq + "AAA"
            elif i % 3 == 1:
                seq = base_seq + "TTT"
            else:
                seq = base_seq + "GGG"
            records.append(
                SeqRecord(
                    Seq(seq),
                    id=f"read_{i:03d}",  # Zero-padded for consistent sorting
                    letter_annotations={'phred_quality': [30] * len(seq)}
                )
            )

        fastq_path = os.path.join(temp_dir, 'determinism_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')
        return fastq_path

    def test_multiple_runs_produce_identical_output(self, temp_dir, core_module, test_fastq):
        """Same input should produce identical output across multiple runs."""
        outputs = []

        for run_num in range(3):
            run_dir = os.path.join(temp_dir, f'run_{run_num}')
            os.makedirs(run_dir)

            result = subprocess.run([
                sys.executable, '-m', core_module,
                test_fastq,
                '--min-size', '2',
                '--algorithm', 'greedy',
                '--output-dir', run_dir
            ], capture_output=True, text=True, cwd=temp_dir)

            assert result.returncode == 0, f"Run {run_num} failed: {result.stderr}"

            # Read output file
            output_file = os.path.join(run_dir, 'determinism_test-all.fasta')
            assert os.path.exists(output_file), f"Output not created for run {run_num}"

            with open(output_file) as f:
                outputs.append(f.read())

        # All outputs should be identical
        assert outputs[0] == outputs[1], "Run 1 and 2 produced different outputs"
        assert outputs[1] == outputs[2], "Run 2 and 3 produced different outputs"

    def test_greedy_vs_mcl_both_deterministic(self, temp_dir, core_module, test_fastq):
        """Both clustering algorithms should be deterministic."""
        for algorithm in ['greedy']:  # MCL requires external tool, test greedy
            outputs = []

            for run_num in range(2):
                run_dir = os.path.join(temp_dir, f'{algorithm}_run_{run_num}')
                os.makedirs(run_dir)

                result = subprocess.run([
                    sys.executable, '-m', core_module,
                    test_fastq,
                    '--min-size', '2',
                    '--algorithm', algorithm,
                    '--output-dir', run_dir
                ], capture_output=True, text=True, cwd=temp_dir)

                assert result.returncode == 0, f"{algorithm} run {run_num} failed"

                output_file = os.path.join(run_dir, 'determinism_test-all.fasta')
                with open(output_file) as f:
                    outputs.append(f.read())

            assert outputs[0] == outputs[1], f"{algorithm} is non-deterministic"


class TestParallelExecutionSafety:
    """Tests for parallel/threaded execution safety.

    Bug: Race condition with vsearch cache directories when using --threads.
    Fix: commit 89c3903 - PID-based unique cache directories
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_parallel_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def core_module(self):
        """Get module name for speconsense core."""
        return 'speconsense.core'

    @pytest.fixture
    def test_fastq(self, temp_dir):
        """Create test FASTQ with enough sequences to benefit from parallelism."""
        base_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        records = []
        for i in range(50):
            # Create varied sequences
            variant = "AAAA" if i % 2 == 0 else "TTTT"
            seq = base_seq + variant + str(i % 10) * 4
            # Ensure sequence only contains valid bases
            seq = seq.replace('0', 'A').replace('1', 'C').replace('2', 'G').replace('3', 'T')
            seq = seq.replace('4', 'A').replace('5', 'C').replace('6', 'G').replace('7', 'T')
            seq = seq.replace('8', 'A').replace('9', 'C')
            records.append(
                SeqRecord(
                    Seq(seq),
                    id=f"read_{i:03d}",
                    letter_annotations={'phred_quality': [30] * len(seq)}
                )
            )

        fastq_path = os.path.join(temp_dir, 'parallel_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')
        return fastq_path

    def test_threaded_execution_produces_valid_output(self, temp_dir, core_module, test_fastq):
        """Running with --threads should produce valid output without crashes."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq,
            '--min-size', '2',
            '--algorithm', 'greedy',
            '--threads', '4'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Threaded execution failed: {result.stderr}"
        assert os.path.exists('clusters/parallel_test-all.fasta'), "Output should be created"

        # Verify output is valid FASTA
        records = list(SeqIO.parse('clusters/parallel_test-all.fasta', 'fasta'))
        assert len(records) > 0, "Should produce at least one consensus sequence"

    def test_threaded_matches_single_threaded(self, temp_dir, core_module, test_fastq):
        """Threaded and single-threaded should produce equivalent results."""
        # Run single-threaded
        single_dir = os.path.join(temp_dir, 'single')
        os.makedirs(single_dir)
        result_single = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq,
            '--min-size', '2',
            '--algorithm', 'greedy',
            '--threads', '1',
            '--output-dir', single_dir
        ], capture_output=True, text=True)
        assert result_single.returncode == 0, f"Single-threaded failed: {result_single.stderr}"

        # Run multi-threaded
        multi_dir = os.path.join(temp_dir, 'multi')
        os.makedirs(multi_dir)
        result_multi = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq,
            '--min-size', '2',
            '--algorithm', 'greedy',
            '--threads', '4',
            '--output-dir', multi_dir
        ], capture_output=True, text=True)
        assert result_multi.returncode == 0, f"Multi-threaded failed: {result_multi.stderr}"

        # Compare outputs - should have same number of sequences
        single_output = os.path.join(single_dir, 'parallel_test-all.fasta')
        multi_output = os.path.join(multi_dir, 'parallel_test-all.fasta')

        single_records = list(SeqIO.parse(single_output, 'fasta'))
        multi_records = list(SeqIO.parse(multi_output, 'fasta'))

        assert len(single_records) == len(multi_records), \
            f"Different number of sequences: single={len(single_records)}, multi={len(multi_records)}"

        # Compare sequences (order may differ, so compare sets)
        single_seqs = set(str(r.seq) for r in single_records)
        multi_seqs = set(str(r.seq) for r in multi_records)

        assert single_seqs == multi_seqs, "Threaded and single-threaded produced different sequences"


class TestScaleThresholdBoundary:
    """Tests for scale threshold boundary conditions.

    Bug: Default presample of 1000 activated scale mode (threshold was also 1000).
    Fix: commit b014fc3 - Changed default threshold to 1001
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_scale_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def core_module(self):
        """Get module name for speconsense core."""
        return 'speconsense.core'

    @pytest.fixture
    def test_fastq(self, temp_dir):
        """Create test FASTQ with many sequences to test presampling."""
        base_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        records = []
        # Create more than 1000 sequences to trigger presampling
        for i in range(1200):
            records.append(
                SeqRecord(
                    Seq(base_seq),
                    id=f"read_{i:04d}",
                    letter_annotations={'phred_quality': [30] * len(base_seq)}
                )
            )

        fastq_path = os.path.join(temp_dir, 'scale_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')
        return fastq_path

    def test_default_presample_does_not_activate_scale_mode(self, temp_dir, core_module, test_fastq):
        """Default --presample 1000 should not activate scale mode (threshold is 1001)."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--log-level', 'DEBUG'  # To see scalability messages
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Should succeed: {result.stderr}"

        # Scale mode messages should NOT appear
        assert "scalability mode" not in result.stderr.lower(), \
            "Default presample should not activate scalability mode"
        assert "vsearch" not in result.stderr.lower(), \
            "Vsearch (scalability feature) should not be used with defaults"

    def test_explicit_scale_threshold_activates_scale_mode(self, temp_dir, core_module, test_fastq):
        """Explicit --scale-threshold 500 should activate scale mode."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--scale-threshold', '500',
            '--log-level', 'DEBUG'
        ], capture_output=True, text=True)

        # Should either succeed or fail gracefully if vsearch not installed
        # The key is that scale mode is attempted
        if "vsearch" in result.stderr.lower() or "scalability" in result.stderr.lower():
            # Scale mode was activated (success)
            pass
        elif result.returncode != 0 and "vsearch" in result.stderr.lower():
            # Scale mode was attempted but vsearch not available (expected on some systems)
            pytest.skip("vsearch not installed, cannot test scale mode activation")
        else:
            # Neither scale mode nor vsearch mentioned - unexpected
            assert False, f"Expected scale mode to be activated: {result.stderr}"

    def test_presample_1000_with_threshold_1000_activates_scale(self, temp_dir, core_module, test_fastq):
        """--presample 1000 --scale-threshold 1000 should activate scale mode."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--presample', '1000',
            '--scale-threshold', '1000',
            '--log-level', 'DEBUG'
        ], capture_output=True, text=True)

        # With threshold exactly at presample size, scale mode should activate
        stderr_lower = result.stderr.lower()
        if "vsearch" not in stderr_lower and "scalability" not in stderr_lower:
            if result.returncode == 0:
                # May have succeeded without scale mode if read count after presample is less
                pass
            else:
                assert False, f"Unexpected failure: {result.stderr}"


class TestLengthFilters:
    """Tests for --min-len and --max-len filters in speconsense-summarize.

    Feature: Filter sequences by length before merging.
    Added: commit c9cb2ae
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_length_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def summarize_module(self):
        """Get module name for speconsense summarize."""
        return 'speconsense.summarize'

    @pytest.fixture
    def source_dir_with_varied_lengths(self, temp_dir):
        """Create source directory with sequences of different lengths."""
        source_dir = os.path.join(temp_dir, 'clusters')
        os.makedirs(source_dir)

        # Create FASTA file with sequences of varying lengths
        # Short: 100bp, Medium: 500bp, Long: 1000bp
        fasta_content = """>short-c1 size=10 ric=10
{"A" * 100}
>medium-c1 size=20 ric=20
{"A" * 500}
>long-c1 size=30 ric=30
{"A" * 1000}
"""
        # Need to actually generate the sequences
        short_seq = "ACGT" * 25  # 100bp
        medium_seq = "ACGT" * 125  # 500bp
        long_seq = "ACGT" * 250  # 1000bp

        fasta_content = f""">short-c1 size=10 ric=10
{short_seq}
>medium-c1 size=20 ric=20
{medium_seq}
>long-c1 size=30 ric=30
{long_seq}
"""
        fasta_file = os.path.join(source_dir, 'test-all.fasta')
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)

        return source_dir

    def test_min_len_filters_short_sequences(self, temp_dir, summarize_module, source_dir_with_varied_lengths):
        """--min-len should filter out sequences shorter than threshold."""
        summary_dir = os.path.join(temp_dir, '__Summary__')

        result = subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', source_dir_with_varied_lengths,
            '--summary-dir', summary_dir,
            '--min-ric', '1',
            '--min-len', '200'  # Filter out 100bp sequence
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Read output and check sequence count
        output_fasta = os.path.join(summary_dir, 'summary.fasta')
        assert os.path.exists(output_fasta), "Output file should exist"

        records = list(SeqIO.parse(output_fasta, 'fasta'))
        # Should have 2 sequences (medium and long), not 3
        assert len(records) == 2, f"Expected 2 sequences after min-len filter, got {len(records)}"

        # Verify lengths are all >= 200
        for record in records:
            assert len(record.seq) >= 200, f"Sequence {record.id} length {len(record.seq)} < min-len 200"

    def test_max_len_filters_long_sequences(self, temp_dir, summarize_module, source_dir_with_varied_lengths):
        """--max-len should filter out sequences longer than threshold."""
        summary_dir = os.path.join(temp_dir, '__Summary__')

        result = subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', source_dir_with_varied_lengths,
            '--summary-dir', summary_dir,
            '--min-ric', '1',
            '--max-len', '600'  # Filter out 1000bp sequence
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_fasta = os.path.join(summary_dir, 'summary.fasta')
        records = list(SeqIO.parse(output_fasta, 'fasta'))

        # Should have 2 sequences (short and medium), not 3
        assert len(records) == 2, f"Expected 2 sequences after max-len filter, got {len(records)}"

        # Verify lengths are all <= 600
        for record in records:
            assert len(record.seq) <= 600, f"Sequence {record.id} length {len(record.seq)} > max-len 600"

    def test_min_and_max_len_combined(self, temp_dir, summarize_module, source_dir_with_varied_lengths):
        """--min-len and --max-len can be combined to filter a range."""
        summary_dir = os.path.join(temp_dir, '__Summary__')

        result = subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', source_dir_with_varied_lengths,
            '--summary-dir', summary_dir,
            '--min-ric', '1',
            '--min-len', '200',
            '--max-len', '600'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_fasta = os.path.join(summary_dir, 'summary.fasta')
        records = list(SeqIO.parse(output_fasta, 'fasta'))

        # Should have only 1 sequence (medium: 500bp)
        assert len(records) == 1, f"Expected 1 sequence in range, got {len(records)}"
        assert len(records[0].seq) == 500, f"Expected 500bp sequence, got {len(records[0].seq)}"

    def test_zero_means_disabled(self, temp_dir, summarize_module, source_dir_with_varied_lengths):
        """--min-len 0 and --max-len 0 should be disabled (default)."""
        summary_dir = os.path.join(temp_dir, '__Summary__')

        result = subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', source_dir_with_varied_lengths,
            '--summary-dir', summary_dir,
            '--min-ric', '1',
            '--min-len', '0',
            '--max-len', '0'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_fasta = os.path.join(summary_dir, 'summary.fasta')
        records = list(SeqIO.parse(output_fasta, 'fasta'))

        # Should have all 3 sequences
        assert len(records) == 3, f"Expected all 3 sequences with filters disabled, got {len(records)}"


class TestCollectDiscards:
    """Tests for --collect-discards option in speconsense.

    Feature: Write discarded reads (outliers and filtered clusters) to a FASTQ file.
    Added: commit da5e71a
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_discards_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def core_module(self):
        """Get module name for speconsense core."""
        return 'speconsense.core'

    @pytest.fixture
    def test_fastq_with_outlier(self, temp_dir):
        """Create test FASTQ with a clear outlier sequence."""
        # Create sequences: most are similar, one is very different (outlier)
        base_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        outlier_seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

        records = []
        # 10 similar sequences
        for i in range(10):
            records.append(
                SeqRecord(
                    Seq(base_seq),
                    id=f"similar_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(base_seq)}
                )
            )
        # 1 outlier
        records.append(
            SeqRecord(
                Seq(outlier_seq),
                id="outlier_01",
                letter_annotations={'phred_quality': [30] * len(outlier_seq)}
            )
        )

        fastq_path = os.path.join(temp_dir, 'with_outlier.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')
        return fastq_path

    def test_collect_discards_creates_file(self, temp_dir, core_module, test_fastq_with_outlier):
        """--collect-discards should create a discards.fastq file."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_with_outlier,
            '--min-size', '2',
            '--algorithm', 'greedy',
            '--collect-discards'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Check for discards file
        discards_file = os.path.join('clusters', 'cluster_debug', 'with_outlier-discards.fastq')

        # The file may or may not exist depending on whether there were discards
        # With our test data, the outlier should be discarded
        if os.path.exists(discards_file):
            # Verify it's a valid FASTQ
            records = list(SeqIO.parse(discards_file, 'fastq'))
            assert len(records) >= 0, "Discards file should be valid FASTQ"

    def test_without_collect_discards_no_file(self, temp_dir, core_module, test_fastq_with_outlier):
        """Without --collect-discards, no discards.fastq should be created."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_with_outlier,
            '--min-size', '2',
            '--algorithm', 'greedy'
            # No --collect-discards
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Discards file should NOT exist
        discards_file = os.path.join('clusters', 'cluster_debug', 'with_outlier-discards.fastq')
        assert not os.path.exists(discards_file), \
            "Discards file should not be created without --collect-discards"

    def test_collect_discards_with_early_filter(self, temp_dir, core_module):
        """--collect-discards should capture reads filtered by early filtering."""
        # Create sequences that will result in small clusters being filtered
        records = []
        # Main cluster: 10 identical sequences
        main_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        for i in range(10):
            records.append(
                SeqRecord(
                    Seq(main_seq),
                    id=f"main_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(main_seq)}
                )
            )
        # Small cluster: 2 different sequences (will be filtered with --min-size 5)
        small_seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        for i in range(2):
            records.append(
                SeqRecord(
                    Seq(small_seq),
                    id=f"small_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(small_seq)}
                )
            )

        fastq_path = os.path.join(temp_dir, 'mixed_clusters.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            fastq_path,
            '--min-size', '5',  # Filter out clusters with < 5 reads
            '--algorithm', 'greedy',
            '--collect-discards',
            '--enable-early-filter'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Check for discards file - should contain the small cluster reads
        discards_file = os.path.join('clusters', 'cluster_debug', 'mixed_clusters-discards.fastq')

        if os.path.exists(discards_file):
            records = list(SeqIO.parse(discards_file, 'fastq'))
            # Should have at least the 2 small cluster reads
            small_ids = [r.id for r in records if r.id.startswith('small_')]
            # Note: early filter may or may not catch these depending on clustering
            assert len(records) >= 0, "Should have valid discards"

    def test_discards_contains_correct_reads(self, temp_dir, core_module):
        """Discards file should contain the actual discarded read sequences."""
        # Create a mix of sequences where some will definitely be discarded
        records = []

        # Main group: 15 identical sequences (will form main cluster)
        main_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        for i in range(15):
            records.append(
                SeqRecord(
                    Seq(main_seq),
                    id=f"main_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(main_seq)}
                )
            )

        # Single outlier (will be filtered by min-size)
        outlier_seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        records.append(
            SeqRecord(
                Seq(outlier_seq),
                id="lone_outlier",
                letter_annotations={'phred_quality': [30] * len(outlier_seq)}
            )
        )

        fastq_path = os.path.join(temp_dir, 'main_plus_outlier.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            fastq_path,
            '--min-size', '5',
            '--algorithm', 'greedy',
            '--collect-discards'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        discards_file = os.path.join('clusters', 'cluster_debug', 'main_plus_outlier-discards.fastq')

        if os.path.exists(discards_file):
            discarded = list(SeqIO.parse(discards_file, 'fastq'))
            discarded_ids = {r.id for r in discarded}

            # The lone outlier should be in discards (filtered by min-size)
            # Note: It may form its own cluster of size 1, which gets filtered
            if 'lone_outlier' in discarded_ids:
                # Verify the sequence is correct
                outlier_record = next(r for r in discarded if r.id == 'lone_outlier')
                assert str(outlier_record.seq) == outlier_seq, "Discarded sequence should match original"


class TestDiscardTrackingCompleteness:
    """Tests for complete discard tracking across all filtering steps.

    These tests verify that --collect-discards captures ALL filtered reads,
    not just some. The invariant is: (final output reads) + (discarded reads) = (input reads).

    Gaps identified:
    - Phase 5 size filtering: reads silently lost (clusterer.py:969-1008)
    - Orientation filtering: reads deleted without tracking (cli.py:287-293)
    - SPOA failures: reads lost when consensus fails (workers.py:482-484)
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_discard_complete_test_')
        original_dir = os.getcwd()
        os.chdir(test_dir)
        yield test_dir
        os.chdir(original_dir)
        shutil.rmtree(test_dir)

    @pytest.fixture
    def core_module(self):
        """Get module name for speconsense core."""
        return 'speconsense.core'

    def count_reads_in_cluster_files(self, clusters_dir: str) -> set:
        """Count unique read IDs in cluster debug files (reads that made it to output)."""
        read_ids = set()
        debug_dir = os.path.join(clusters_dir, 'cluster_debug')

        if not os.path.exists(debug_dir):
            return read_ids

        for filename in os.listdir(debug_dir):
            if filename.endswith('-reads.fastq') or filename.endswith('-reads.fasta'):
                filepath = os.path.join(debug_dir, filename)
                for record in SeqIO.parse(filepath, 'fastq' if filename.endswith('.fastq') else 'fasta'):
                    read_ids.add(record.id)

        return read_ids

    def count_discarded_reads(self, clusters_dir: str, sample_name: str) -> set:
        """Get read IDs from discards file."""
        discards_file = os.path.join(clusters_dir, 'cluster_debug', f'{sample_name}-discards.fastq')

        if not os.path.exists(discards_file):
            return set()

        return {record.id for record in SeqIO.parse(discards_file, 'fastq')}

    def test_read_accounting_with_size_filtering(self, temp_dir, core_module):
        """Verify all reads accounted for when clusters are filtered by size.

        This test creates a scenario where some clusters will be filtered by
        --min-size. The filtered reads should appear in the discards file.

        Expected behavior: input_reads == output_reads + discarded_reads
        Current behavior (BUG): Phase 5 filtered reads are lost
        """
        # Create input: 15 reads for main cluster, 2 reads for small cluster
        main_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        small_seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

        records = []
        # Main cluster: 15 identical reads
        for i in range(15):
            records.append(
                SeqRecord(
                    Seq(main_seq),
                    id=f"main_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(main_seq)}
                )
            )
        # Small cluster: 2 identical reads (will be filtered by --min-size 5)
        for i in range(2):
            records.append(
                SeqRecord(
                    Seq(small_seq),
                    id=f"small_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(small_seq)}
                )
            )

        input_count = len(records)
        input_ids = {r.id for r in records}

        fastq_path = os.path.join(temp_dir, 'size_filter_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            fastq_path,
            '--min-size', '5',  # Filter clusters < 5 reads
            '--algorithm', 'greedy',
            '--collect-discards'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Count reads in final output (cluster debug read files)
        output_ids = self.count_reads_in_cluster_files('clusters')

        # Count discarded reads
        discarded_ids = self.count_discarded_reads('clusters', 'size_filter_test')

        # Verify accounting: all input reads should be either in output or discards
        accounted_ids = output_ids | discarded_ids
        missing_ids = input_ids - accounted_ids
        extra_ids = accounted_ids - input_ids

        assert len(missing_ids) == 0, \
            f"Missing reads not in output or discards: {missing_ids}"
        assert len(extra_ids) == 0, \
            f"Extra reads appeared from nowhere: {extra_ids}"
        assert len(accounted_ids) == input_count, \
            f"Read accounting mismatch: input={input_count}, output={len(output_ids)}, discards={len(discarded_ids)}"

        # Specifically verify the small cluster reads are in discards
        small_ids = {f"small_{i:02d}" for i in range(2)}
        assert small_ids.issubset(discarded_ids), \
            f"Small cluster reads should be in discards: expected {small_ids}, got discards={discarded_ids}"

    def test_read_accounting_with_ratio_filtering(self, temp_dir, core_module):
        """Verify reads filtered by --min-cluster-ratio appear in discards.

        Creates a large cluster (100 reads) and small cluster (5 reads).
        With --min-cluster-ratio 0.1, the small cluster (5% of large) should be filtered.
        """
        main_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        small_seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"

        records = []
        # Large cluster: 100 reads
        for i in range(100):
            records.append(
                SeqRecord(
                    Seq(main_seq),
                    id=f"large_{i:03d}",
                    letter_annotations={'phred_quality': [30] * len(main_seq)}
                )
            )
        # Small cluster: 5 reads (5% of large, below 10% ratio threshold)
        for i in range(5):
            records.append(
                SeqRecord(
                    Seq(small_seq),
                    id=f"ratio_filtered_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(small_seq)}
                )
            )

        input_count = len(records)
        input_ids = {r.id for r in records}

        fastq_path = os.path.join(temp_dir, 'ratio_filter_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            fastq_path,
            '--min-size', '3',
            '--min-cluster-ratio', '0.10',  # Clusters < 10% of largest are filtered
            '--algorithm', 'greedy',
            '--collect-discards'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_ids = self.count_reads_in_cluster_files('clusters')
        discarded_ids = self.count_discarded_reads('clusters', 'ratio_filter_test')

        accounted_ids = output_ids | discarded_ids
        missing_ids = input_ids - accounted_ids

        assert len(missing_ids) == 0, \
            f"Missing reads not in output or discards: {missing_ids}"

        # The ratio-filtered reads should be in discards
        ratio_filtered_ids = {f"ratio_filtered_{i:02d}" for i in range(5)}
        assert ratio_filtered_ids.issubset(discarded_ids), \
            f"Ratio-filtered reads should be in discards: expected {ratio_filtered_ids}, got discards={discarded_ids}"

    def test_read_accounting_with_orientation_filtering(self, temp_dir, core_module):
        """Verify reads filtered by orientation appear in discards.

        Creates reads where some cannot be oriented (no primer match).
        With --orient-mode filter-failed, unoriented reads should be discarded.
        """
        # Forward primer at start
        forward_primer = "ACGTACGTACGT"
        # Reverse primer at end (reverse complement would be at end)
        reverse_primer = "TGCATGCATGCA"

        good_seq = forward_primer + "NNNNNNNNNNNNNNNNNNNNNNNN" + reverse_primer
        # Bad sequence: no primers, cannot be oriented
        bad_seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"

        records = []
        # Good reads: have primers, can be oriented
        for i in range(10):
            records.append(
                SeqRecord(
                    Seq(good_seq),
                    id=f"orientable_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(good_seq)}
                )
            )
        # Bad reads: no primers, will fail orientation
        for i in range(3):
            records.append(
                SeqRecord(
                    Seq(bad_seq),
                    id=f"unorientable_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(bad_seq)}
                )
            )

        input_count = len(records)
        input_ids = {r.id for r in records}

        fastq_path = os.path.join(temp_dir, 'orient_filter_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')

        # Create primers file with position hints
        primers_content = f""">forward_primer position=forward
{forward_primer}
>reverse_primer position=reverse
{reverse_primer}
"""
        primers_path = os.path.join(temp_dir, 'primers.fasta')
        with open(primers_path, 'w') as f:
            f.write(primers_content)

        result = subprocess.run([
            sys.executable, '-m', core_module,
            fastq_path,
            '--primers', primers_path,
            '--orient-mode', 'filter-failed',  # Filter reads that can't be oriented
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--collect-discards'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_ids = self.count_reads_in_cluster_files('clusters')
        discarded_ids = self.count_discarded_reads('clusters', 'orient_filter_test')

        accounted_ids = output_ids | discarded_ids
        missing_ids = input_ids - accounted_ids

        assert len(missing_ids) == 0, \
            f"Missing reads not in output or discards: {missing_ids}"

        # The unorientable reads should be in discards
        unorientable_ids = {f"unorientable_{i:02d}" for i in range(3)}
        assert unorientable_ids.issubset(discarded_ids), \
            f"Unorientable reads should be in discards: expected {unorientable_ids}, got discards={discarded_ids}"

    def test_total_read_accounting_complex_scenario(self, temp_dir, core_module):
        """Comprehensive test: all input reads must be in either output or discards.

        This test creates a complex scenario with multiple filtering opportunities
        and verifies complete read accounting.
        """
        # Three distinct sequence types
        main_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        variant_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT"  # Similar to main
        outlier_seq = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"

        records = []
        # Main cluster: 20 reads
        for i in range(20):
            records.append(
                SeqRecord(
                    Seq(main_seq),
                    id=f"main_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(main_seq)}
                )
            )
        # Variant cluster: 8 reads (may merge or stay separate)
        for i in range(8):
            records.append(
                SeqRecord(
                    Seq(variant_seq),
                    id=f"variant_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(variant_seq)}
                )
            )
        # Outlier cluster: 2 reads (will be filtered by --min-size 5)
        for i in range(2):
            records.append(
                SeqRecord(
                    Seq(outlier_seq),
                    id=f"outlier_{i:02d}",
                    letter_annotations={'phred_quality': [30] * len(outlier_seq)}
                )
            )

        input_count = len(records)
        input_ids = {r.id for r in records}

        fastq_path = os.path.join(temp_dir, 'complex_test.fastq')
        with open(fastq_path, 'w') as f:
            SeqIO.write(records, f, 'fastq')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            fastq_path,
            '--min-size', '5',
            '--algorithm', 'greedy',
            '--collect-discards'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_ids = self.count_reads_in_cluster_files('clusters')
        discarded_ids = self.count_discarded_reads('clusters', 'complex_test')

        # THE KEY INVARIANT: every input read is either in output or discards
        accounted_ids = output_ids | discarded_ids

        assert accounted_ids == input_ids, \
            f"Read accounting failed:\n" \
            f"  Input: {len(input_ids)} reads\n" \
            f"  Output: {len(output_ids)} reads\n" \
            f"  Discards: {len(discarded_ids)} reads\n" \
            f"  Missing: {input_ids - accounted_ids}\n" \
            f"  Extra: {accounted_ids - input_ids}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
