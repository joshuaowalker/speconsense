"""Tests for complement boolean CLI flags."""

import sys
from unittest.mock import patch

import pytest


class TestSummarizeComplementFlags:
    """Test complement flags in speconsense-summarize CLI."""

    def _parse(self, args_list):
        """Parse arguments for speconsense-summarize with mocked sys.argv."""
        with patch.object(sys, 'argv', ['speconsense-summarize'] + args_list):
            from speconsense.summarize.cli import parse_arguments
            return parse_arguments()

    def test_disable_merging_default(self):
        args = self._parse([])
        assert args.disable_merging is False

    def test_disable_merging_flag(self):
        args = self._parse(['--disable-merging'])
        assert args.disable_merging is True

    def test_enable_merging_complement(self):
        args = self._parse(['--enable-merging'])
        assert args.disable_merging is False

    def test_disable_homopolymer_equivalence_default(self):
        args = self._parse([])
        assert args.disable_homopolymer_equivalence is False

    def test_disable_homopolymer_equivalence_flag(self):
        args = self._parse(['--disable-homopolymer-equivalence'])
        assert args.disable_homopolymer_equivalence is True

    def test_enable_homopolymer_equivalence_complement(self):
        args = self._parse(['--enable-homopolymer-equivalence'])
        assert args.disable_homopolymer_equivalence is False

    def test_enable_full_consensus_default(self):
        args = self._parse([])
        assert args.enable_full_consensus is False

    def test_enable_full_consensus_flag(self):
        args = self._parse(['--enable-full-consensus'])
        assert args.enable_full_consensus is True

    def test_disable_full_consensus_complement(self):
        args = self._parse(['--disable-full-consensus'])
        assert args.enable_full_consensus is False

    def test_complement_overrides_profile_disable_merging(self):
        """--enable-merging should override profile that sets disable-merging: true."""
        args = self._parse(['-p', 'nostalgia', '--enable-merging'])
        # Regardless of what nostalgia sets, the explicit complement should win
        assert args.disable_merging is False

    def test_complement_overrides_profile_enable_full_consensus(self):
        """--disable-full-consensus should override profile that sets enable-full-consensus: true."""
        # Use a profile that sets enable-full-consensus: true
        # We test the mechanism: even if profile sets it, CLI complement wins
        args = self._parse(['-p', 'herbarium', '--disable-full-consensus'])
        assert args.enable_full_consensus is False


class TestCoreComplementFlags:
    """Test complement flags in speconsense core CLI.

    Tests argparse behavior by constructing the parser directly
    rather than going through main() which requires an input file.
    """

    @pytest.fixture(autouse=True)
    def _setup_parser(self):
        """Build the core parser for testing."""
        import argparse
        parser = argparse.ArgumentParser()

        # Replicate the boolean flags and their complements from core/cli.py
        parser.add_argument("--disable-position-phasing", action="store_true")
        parser.add_argument("--enable-position-phasing", action="store_false",
                            dest="disable_position_phasing")

        parser.add_argument("--disable-ambiguity-calling", action="store_true")
        parser.add_argument("--enable-ambiguity-calling", action="store_false",
                            dest="disable_ambiguity_calling")

        parser.add_argument("--disable-cluster-merging", action="store_true")
        parser.add_argument("--enable-cluster-merging", action="store_false",
                            dest="disable_cluster_merging")

        parser.add_argument("--disable-homopolymer-equivalence", action="store_true")
        parser.add_argument("--enable-homopolymer-equivalence", action="store_false",
                            dest="disable_homopolymer_equivalence")

        parser.add_argument("--enable-early-filter", action="store_true")
        parser.add_argument("--disable-early-filter", action="store_false",
                            dest="enable_early_filter")

        parser.add_argument("--collect-discards", action="store_true")
        parser.add_argument("--no-collect-discards", action="store_false",
                            dest="collect_discards")

        self.parser = parser

    def test_defaults_all_false(self):
        args = self.parser.parse_args([])
        assert args.disable_position_phasing is False
        assert args.disable_ambiguity_calling is False
        assert args.disable_cluster_merging is False
        assert args.disable_homopolymer_equivalence is False
        assert args.enable_early_filter is False
        assert args.collect_discards is False

    def test_primary_flags_set_true(self):
        args = self.parser.parse_args([
            '--disable-position-phasing',
            '--disable-ambiguity-calling',
            '--disable-cluster-merging',
            '--disable-homopolymer-equivalence',
            '--enable-early-filter',
            '--collect-discards',
        ])
        assert args.disable_position_phasing is True
        assert args.disable_ambiguity_calling is True
        assert args.disable_cluster_merging is True
        assert args.disable_homopolymer_equivalence is True
        assert args.enable_early_filter is True
        assert args.collect_discards is True

    def test_complement_flags_set_false(self):
        args = self.parser.parse_args([
            '--enable-position-phasing',
            '--enable-ambiguity-calling',
            '--enable-cluster-merging',
            '--enable-homopolymer-equivalence',
            '--disable-early-filter',
            '--no-collect-discards',
        ])
        assert args.disable_position_phasing is False
        assert args.disable_ambiguity_calling is False
        assert args.disable_cluster_merging is False
        assert args.disable_homopolymer_equivalence is False
        assert args.enable_early_filter is False
        assert args.collect_discards is False

    def test_complement_overrides_profile_default(self):
        """Complement flag should override a profile-set default (via set_defaults)."""
        # Simulate profile setting disable_position_phasing=True
        self.parser.set_defaults(disable_position_phasing=True)
        args = self.parser.parse_args(['--enable-position-phasing'])
        assert args.disable_position_phasing is False

    def test_complement_overrides_profile_enable_early_filter(self):
        """--disable-early-filter should override profile that sets enable_early_filter=True."""
        self.parser.set_defaults(enable_early_filter=True)
        args = self.parser.parse_args(['--disable-early-filter'])
        assert args.enable_early_filter is False

    def test_complement_overrides_profile_collect_discards(self):
        """--no-collect-discards should override profile that sets collect_discards=True."""
        self.parser.set_defaults(collect_discards=True)
        args = self.parser.parse_args(['--no-collect-discards'])
        assert args.collect_discards is False

    def test_last_flag_wins(self):
        """When both primary and complement are passed, last one wins."""
        args = self.parser.parse_args([
            '--disable-position-phasing',
            '--enable-position-phasing',
        ])
        assert args.disable_position_phasing is False

        args = self.parser.parse_args([
            '--enable-position-phasing',
            '--disable-position-phasing',
        ])
        assert args.disable_position_phasing is True
