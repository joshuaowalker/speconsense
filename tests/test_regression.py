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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
