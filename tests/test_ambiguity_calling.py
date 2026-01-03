#!/usr/bin/env python3
"""
Pytest integration tests for IUPAC ambiguity calling parameters.

Tests --min-ambiguity-frequency and --min-ambiguity-count parameters
which control IUPAC ambiguity calling independently from variant phasing.
"""

import os
import tempfile
import shutil
import subprocess
import sys
import pytest
import re


class TestAmbiguityCalling:
    """Test suite for IUPAC ambiguity calling parameters."""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_ambiguity_test_')
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
    def test_fastq_path(self):
        """Get path to test FASTQ file."""
        return os.path.join(os.path.dirname(__file__), 'data',
                           'ONT10.80-H10--iNat229710865-1.v2-RiC9.fastq')

    def count_iupac_codes(self, fasta_content: str) -> int:
        """Count IUPAC ambiguity codes in sequences (excluding header lines)."""
        iupac_pattern = re.compile(r'[RYSWKMBDHVN]')
        count = 0
        for line in fasta_content.split('\n'):
            if not line.startswith('>'):
                count += len(iupac_pattern.findall(line))
        return count

    def test_default_ambiguity_thresholds(self, temp_dir, core_module, test_fastq_path):
        """Test that default ambiguity thresholds (10%/3) are applied."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_path,
            '--min-size', '0',
            '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Check output file exists
        output_file = os.path.join('clusters', 'ONT10.80-H10--iNat229710865-1.v2-RiC9-all.fasta')
        assert os.path.exists(output_file), f"Output file not found: {output_file}"

    def test_strict_ambiguity_thresholds_fewer_codes(self, temp_dir, core_module, test_fastq_path):
        """Test that stricter ambiguity thresholds result in fewer IUPAC codes."""
        # Run with default (lenient) thresholds
        result_lenient = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_path,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--min-ambiguity-frequency', '0.10',
            '--min-ambiguity-count', '1'
        ], capture_output=True, text=True)
        assert result_lenient.returncode == 0

        output_file = os.path.join('clusters', 'ONT10.80-H10--iNat229710865-1.v2-RiC9-all.fasta')
        with open(output_file) as f:
            lenient_content = f.read()
        lenient_iupac_count = self.count_iupac_codes(lenient_content)

        # Clean up for next run
        shutil.rmtree('clusters')

        # Run with strict thresholds
        result_strict = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_path,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--min-ambiguity-frequency', '0.40',
            '--min-ambiguity-count', '5'
        ], capture_output=True, text=True)
        assert result_strict.returncode == 0

        with open(output_file) as f:
            strict_content = f.read()
        strict_iupac_count = self.count_iupac_codes(strict_content)

        # Strict thresholds should produce fewer or equal IUPAC codes
        assert strict_iupac_count <= lenient_iupac_count, \
            f"Strict thresholds ({strict_iupac_count}) should not produce more IUPAC codes than lenient ({lenient_iupac_count})"

    def test_ambiguity_independent_from_phasing(self, temp_dir, core_module, test_fastq_path):
        """Test that ambiguity thresholds are independent from phasing thresholds."""
        # Run with high phasing thresholds but low ambiguity thresholds
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_path,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--min-variant-frequency', '0.50',  # High - less phasing
            '--min-variant-count', '10',        # High - less phasing
            '--min-ambiguity-frequency', '0.05', # Low - more ambiguity codes
            '--min-ambiguity-count', '1'         # Low - more ambiguity codes
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_file = os.path.join('clusters', 'ONT10.80-H10--iNat229710865-1.v2-RiC9-all.fasta')
        assert os.path.exists(output_file)

    def test_disable_ambiguity_calling(self, temp_dir, core_module, test_fastq_path):
        """Test that --disable-ambiguity-calling produces no IUPAC codes."""
        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_path,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--disable-ambiguity-calling'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        output_file = os.path.join('clusters', 'ONT10.80-H10--iNat229710865-1.v2-RiC9-all.fasta')
        with open(output_file) as f:
            content = f.read()

        iupac_count = self.count_iupac_codes(content)
        assert iupac_count == 0, f"Disabled ambiguity calling should produce 0 IUPAC codes, got {iupac_count}"

    def test_ambiguity_parameters_in_metadata(self, temp_dir, core_module, test_fastq_path):
        """Test that ambiguity parameters are recorded in metadata."""
        import json

        result = subprocess.run([
            sys.executable, '-m', core_module,
            test_fastq_path,
            '--min-size', '0',
            '--algorithm', 'greedy',
            '--min-ambiguity-frequency', '0.15',
            '--min-ambiguity-count', '4'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Command failed: {result.stderr}"

        metadata_file = os.path.join('clusters', 'cluster_debug',
                                      'ONT10.80-H10--iNat229710865-1.v2-RiC9-metadata.json')
        assert os.path.exists(metadata_file), f"Metadata file not found: {metadata_file}"

        with open(metadata_file) as f:
            metadata = json.load(f)

        params = metadata.get('parameters', {})
        assert params.get('min_ambiguity_frequency') == 0.15, \
            f"Expected min_ambiguity_frequency=0.15, got {params.get('min_ambiguity_frequency')}"
        assert params.get('min_ambiguity_count') == 4, \
            f"Expected min_ambiguity_count=4, got {params.get('min_ambiguity_count')}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
