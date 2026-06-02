#!/usr/bin/env python3
"""
Pytest integration tests for --augment-input functionality.
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


class TestAugmentInput:
    """Test suite for --augment-input functionality."""
    
    @pytest.fixture
    def test_data(self):
        """Create test FASTQ and FASTA files for testing."""
        # Create similar sequences that should cluster together
        seq1 = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
        seq2 = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'  # identical
        seq3 = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA'  # 1bp diff
        
        # Augmented sequence (identical to seq1/seq2)
        aug_seq = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
        
        # Different sequence for separate cluster
        different_seq = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        
        # Create main FASTQ file
        main_records = [
            SeqRecord(Seq(seq1), id='primary_1', letter_annotations={'phred_quality': [30]*len(seq1)}),
            SeqRecord(Seq(seq2), id='primary_2', letter_annotations={'phred_quality': [30]*len(seq2)}),
            SeqRecord(Seq(seq3), id='primary_3', letter_annotations={'phred_quality': [30]*len(seq3)}),
            SeqRecord(Seq(different_seq), id='primary_different', letter_annotations={'phred_quality': [30]*len(different_seq)})
        ]
        
        # Create augmented FASTA file (testing format auto-detection)
        augment_records = [
            SeqRecord(Seq(aug_seq), id='augmented_1', description='reference_sequence')
        ]
        
        return {
            'main_records': main_records,
            'augment_records': augment_records,
            'main_count': len(main_records),
            'augment_count': len(augment_records)
        }
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        test_dir = tempfile.mkdtemp(prefix='speconsense_test_')
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
    def summarize_module(self):
        """Get module name for speconsense summarize."""
        return 'speconsense.summarize'
    
    def create_test_files(self, test_data):
        """Create test files in current directory."""
        with open('test_main.fastq', 'w') as f:
            SeqIO.write(test_data['main_records'], f, 'fastq')
        
        with open('test_augment.fasta', 'w') as f:
            SeqIO.write(test_data['augment_records'], f, 'fasta')
    
    def test_error_handling_nonexistent_file(self, temp_dir, core_module, test_data):
        """Test error handling for nonexistent augment input file."""
        self.create_test_files(test_data)

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'nonexistent.fastq',
            '--min-size', '2', '--algorithm', 'greedy'
        ], capture_output=True, text=True)
        
        assert result.returncode == 1, "Should fail with exit code 1 for nonexistent file"
        assert "not found" in result.stderr, "Should show file not found error"
    
    def test_augmented_sequences_loaded(self, temp_dir, core_module, test_data):
        """Test that augmented sequences are properly loaded."""
        self.create_test_files(test_data)

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy',
            '--log-level', 'INFO'
        ], capture_output=True, text=True)
        
        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Loaded 1 augmented reads" in result.stderr, "Should load augmented sequences"
        assert "Loaded 4 primary reads" in result.stderr, "Should load primary sequences"
    
    def test_output_files_created(self, temp_dir, core_module, test_data):
        """Test that output files are created correctly."""
        self.create_test_files(test_data)

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        assert result.returncode == 0, "Speconsense should succeed"
        assert os.path.exists('clusters/test_main-all.fasta'), "Main output file should be created"
        assert os.path.exists('clusters/cluster_debug'), "Debug directory should be created"
    
    def test_augmented_sequence_in_cluster_output(self, temp_dir, core_module, test_data):
        """Test that augmented sequences appear in cluster debug files."""
        self.create_test_files(test_data)

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        assert result.returncode == 0, "Speconsense should succeed"

        # Check cluster debug files
        debug_files = [f for f in os.listdir('clusters/cluster_debug') 
                      if f.endswith('-reads.fastq') or f.endswith('-reads.fasta')]
        assert len(debug_files) > 0, "Should create cluster debug read files"
        
        # Check that augmented sequence is in cluster output
        with open(os.path.join('clusters/cluster_debug', debug_files[0]), 'r') as f:
            debug_content = f.read()
            assert 'augmented_1' in debug_content, "Augmented sequence should be in cluster output"
    
    def test_summarize_integration(self, temp_dir, core_module, summarize_module, test_data):
        """Test that summarize step handles augmented sequences correctly."""
        self.create_test_files(test_data)

        # Run speconsense
        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        assert result.returncode == 0, "Speconsense should succeed"

        # Run summarize
        result = subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', 'clusters', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Summarize should succeed: {result.stderr}"
        assert os.path.exists('__Summary__/FASTQ Files'), "Summary FASTQ Files directory should be created"
    
    def test_augmented_sequence_in_summary_output(self, temp_dir, core_module, summarize_module, test_data):
        """Test that augmented sequences appear in final summary output."""
        self.create_test_files(test_data)

        # Run speconsense
        subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        # Run summarize
        subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', 'clusters', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        # Check that augmented sequence is in final FASTQ output
        fastq_files = os.listdir('__Summary__/FASTQ Files')
        assert len(fastq_files) > 0, "Should create summary FASTQ files"
        
        with open(os.path.join('__Summary__/FASTQ Files', fastq_files[0]), 'r') as f:
            summary_content = f.read()
            assert 'augmented_1' in summary_content, "Augmented sequence should be in summary output"
    
    def test_final_consensus_counts(self, temp_dir, core_module, summarize_module, test_data):
        """Test that final consensus headers show correct sequence counts."""
        self.create_test_files(test_data)

        # Run speconsense
        subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        # Run summarize
        subprocess.run([
            sys.executable, '-m', summarize_module,
            '--source', 'clusters', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        # Check final FASTA header shows correct counts
        fasta_files = [f for f in os.listdir('__Summary__') 
                      if f.endswith('.fasta') and f != 'summary.fasta']
        assert len(fasta_files) > 0, "Should create final consensus FASTA files"
        
        with open(os.path.join('__Summary__', fasta_files[0]), 'r') as f:
            header = f.readline().strip()
            # Should show size=4 (3 similar primary + 1 augmented) or size=5 (all sequences)
            assert ('size=4' in header or 'size=5' in header), f"Unexpected sequence count in header: {header}"
    
    def test_multiple_augment_files(self, temp_dir, core_module, test_data):
        """Test that repeated --augment-input loads reads from every file."""
        self.create_test_files(test_data)

        # Second augment file with another read that clusters with the primaries
        aug_seq2 = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
        with open('test_augment2.fasta', 'w') as f:
            SeqIO.write([SeqRecord(Seq(aug_seq2), id='augmented_2', description='')], f, 'fasta')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq',
            '--augment-input', 'test_augment.fasta',
            '--augment-input', 'test_augment2.fasta',
            '--min-size', '2', '--algorithm', 'greedy',
            '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Loaded 2 augmented reads total from 2 files" in result.stderr

        # Both augmented reads should appear in cluster output
        debug_content = ""
        for fname in os.listdir('clusters/cluster_debug'):
            if fname.endswith('-reads.fastq'):
                with open(os.path.join('clusters/cluster_debug', fname)) as f:
                    debug_content += f.read()
        assert 'augmented_1' in debug_content, "First augment file's read should be clustered"
        assert 'augmented_2' in debug_content, "Second augment file's read should be clustered"

    def test_augmented_only_cluster_dropped(self, temp_dir, core_module, test_data):
        """Off-target augmented reads (no primary within --min-identity) never form a cluster.

        They are removed by the pre-clustering augment screen; the end result is the
        same as the Phase-12 augmented-only filter would produce (absent from output,
        present in discards).
        """
        self.create_test_files(test_data)

        # Augment file forming its own cluster (distinct from every primary)
        aug_only = 'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA'
        recs = [
            SeqRecord(Seq(aug_only), id='augonly_1', description=''),
            SeqRecord(Seq(aug_only), id='augonly_2', description=''),
            SeqRecord(Seq(aug_only), id='augonly_3', description=''),
        ]
        with open('test_augonly.fasta', 'w') as f:
            SeqIO.write(recs, f, 'fasta')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augonly.fasta',
            '--min-size', '2', '--algorithm', 'greedy',
            '--collect-discards', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Augment screen: kept 0/3" in result.stderr

        # The augmented-only consensus must not reach the final FASTA
        with open('clusters/test_main-all.fasta') as f:
            all_fasta = f.read()
        assert 'TGCATGCATGCA' not in all_fasta, "Augmented-only cluster should be absent from output"

        # Its reads should land in the discards file
        discards = os.path.join('clusters/cluster_debug', 'test_main-discards.fastq')
        assert os.path.exists(discards), "Discards file should be written with --collect-discards"
        with open(discards) as f:
            discard_content = f.read()
        assert 'augonly_1' in discard_content and 'augonly_2' in discard_content, \
            "Augmented-only reads should be discarded"

    def test_keep_augmented_only_clusters(self, temp_dir, core_module, test_data):
        """--keep-augmented-only-clusters retains clusters with no primary read."""
        self.create_test_files(test_data)

        aug_only = 'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA'
        recs = [
            SeqRecord(Seq(aug_only), id='augonly_1', description=''),
            SeqRecord(Seq(aug_only), id='augonly_2', description=''),
            SeqRecord(Seq(aug_only), id='augonly_3', description=''),
        ]
        with open('test_augonly.fasta', 'w') as f:
            SeqIO.write(recs, f, 'fasta')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augonly.fasta',
            '--min-size', '2', '--algorithm', 'greedy',
            '--keep-augmented-only-clusters', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Augmented-only filter" not in result.stderr

        with open('clusters/test_main-all.fasta') as f:
            all_fasta = f.read()
        assert 'TGCATGCATGCA' in all_fasta, "Augmented-only cluster should be retained with opt-out flag"

    def test_duplicate_id_collision_dedup(self, temp_dir, core_module, test_data):
        """A read ID present in both primary and augment inputs is deduped (primary wins)."""
        self.create_test_files(test_data)

        # Augment file reusing a primary read ID
        dup_seq = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
        with open('test_dup.fasta', 'w') as f:
            SeqIO.write([SeqRecord(Seq(dup_seq), id='primary_1', description='')], f, 'fasta')

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_dup.fasta',
            '--min-size', '2', '--algorithm', 'greedy',
            '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Skipped 1 read(s) with duplicate IDs" in result.stderr

    def test_augment_screen_keeps_matching(self, temp_dir, core_module, test_data):
        """Augmented reads within --min-identity of a primary survive the screen."""
        self.create_test_files(test_data)  # test_augment.fasta is identical to the primaries

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Augment screen: kept 1/1" in result.stderr, "Matching read should pass the screen"

    def test_augment_screen_skipped_when_keeping(self, temp_dir, core_module, test_data):
        """The screen is skipped under --keep-augmented-only-clusters."""
        self.create_test_files(test_data)

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '2', '--algorithm', 'greedy',
            '--keep-augmented-only-clusters', '--log-level', 'INFO'
        ], capture_output=True, text=True)

        assert result.returncode == 0, f"Speconsense should succeed: {result.stderr}"
        assert "Augment screen" not in result.stderr, "Screen must not run when keeping augmented-only clusters"

    def test_empty_augment_file_warning(self, temp_dir, core_module):
        """Test warning for empty augment input file."""
        # Create main file with one sequence, empty augment file
        with open('test_main.fastq', 'w') as f:
            f.write("@read1\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("IIIIIIIIIIII\n")
        with open('test_augment.fasta', 'w') as f:
            f.write("")

        result = subprocess.run([
            sys.executable, '-m', core_module,
            'test_main.fastq', '--augment-input', 'test_augment.fasta',
            '--min-size', '0', '--algorithm', 'greedy'
        ], capture_output=True, text=True)

        # Should warn about empty augment file
        assert "No reads found in augment input file" in result.stderr


if __name__ == "__main__":
    pytest.main([__file__, "-v"])