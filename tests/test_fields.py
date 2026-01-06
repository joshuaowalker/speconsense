"""Tests for FASTA field formatting."""

import pytest

from speconsense.types import ConsensusInfo
from speconsense.summarize.fields import (
    # Field classes
    SizeField,
    RicField,
    LengthField,
    RawRicField,
    RawLenField,
    SnpField,
    AmbigField,
    PrimersField,
    RidField,
    RidMinField,
    GroupField,
    VariantField,
    # Registry and presets
    FASTA_FIELDS,
    FASTA_FIELD_PRESETS,
    # Functions
    parse_fasta_fields,
    format_fasta_header,
)


def make_consensus(
    sample_name="specimen-1.v1",
    sequence="ACGT",
    ric=10,
    size=15,
    raw_ric=None,
    raw_len=None,
    snp_count=None,
    primers=None,
    rid=None,
    rid_min=None,
):
    """Helper to create ConsensusInfo with sensible defaults."""
    return ConsensusInfo(
        sample_name=sample_name,
        cluster_id="c1",
        sequence=sequence,
        ric=ric,
        size=size,
        file_path="/test/path.fasta",
        snp_count=snp_count,
        primers=primers,
        raw_ric=raw_ric,
        raw_len=raw_len,
        rid=rid,
        rid_min=rid_min,
    )


class TestSizeField:
    """Tests for the size field."""

    def test_format_value(self):
        consensus = make_consensus(size=100)
        field = SizeField()
        assert field.format_value(consensus) == "size=100"

    def test_format_value_zero(self):
        consensus = make_consensus(size=0)
        field = SizeField()
        assert field.format_value(consensus) == "size=0"


class TestRicField:
    """Tests for the RiC field."""

    def test_format_value(self):
        consensus = make_consensus(ric=50)
        field = RicField()
        assert field.format_value(consensus) == "ric=50"

    def test_format_value_zero(self):
        consensus = make_consensus(ric=0)
        field = RicField()
        assert field.format_value(consensus) == "ric=0"


class TestLengthField:
    """Tests for the length field."""

    def test_format_value(self):
        consensus = make_consensus(sequence="ACGTACGTACGT")
        field = LengthField()
        assert field.format_value(consensus) == "length=12"

    def test_format_value_empty_sequence(self):
        consensus = make_consensus(sequence="")
        field = LengthField()
        assert field.format_value(consensus) == "length=0"


class TestRawRicField:
    """Tests for the rawric field - tracks RiC values of merged source variants."""

    def test_format_value_with_values(self):
        """rawric should format when raw_ric list is populated."""
        consensus = make_consensus(raw_ric=[50, 30, 10])
        field = RawRicField()
        result = field.format_value(consensus)
        # Values should be sorted descending
        assert result == "rawric=50+30+10"

    def test_format_value_unsorted_input(self):
        """rawric should sort values descending regardless of input order."""
        consensus = make_consensus(raw_ric=[10, 50, 30])
        field = RawRicField()
        result = field.format_value(consensus)
        assert result == "rawric=50+30+10"

    def test_format_value_single_value(self):
        """rawric with single value should still format."""
        consensus = make_consensus(raw_ric=[100])
        field = RawRicField()
        result = field.format_value(consensus)
        assert result == "rawric=100"

    def test_format_value_none(self):
        """rawric should return None when raw_ric is None."""
        consensus = make_consensus(raw_ric=None)
        field = RawRicField()
        assert field.format_value(consensus) is None

    def test_format_value_empty_list(self):
        """rawric should return None when raw_ric is empty list."""
        consensus = make_consensus(raw_ric=[])
        field = RawRicField()
        assert field.format_value(consensus) is None


class TestRawLenField:
    """Tests for the rawlen field - tracks lengths of merged source sequences."""

    def test_format_value_with_values(self):
        consensus = make_consensus(raw_len=[700, 650, 620])
        field = RawLenField()
        result = field.format_value(consensus)
        assert result == "rawlen=700+650+620"

    def test_format_value_unsorted_input(self):
        """rawlen should sort values descending."""
        consensus = make_consensus(raw_len=[620, 700, 650])
        field = RawLenField()
        result = field.format_value(consensus)
        assert result == "rawlen=700+650+620"

    def test_format_value_none(self):
        consensus = make_consensus(raw_len=None)
        field = RawLenField()
        assert field.format_value(consensus) is None

    def test_format_value_empty_list(self):
        consensus = make_consensus(raw_len=[])
        field = RawLenField()
        assert field.format_value(consensus) is None


class TestSnpField:
    """Tests for the SNP count field."""

    def test_format_value_with_count(self):
        consensus = make_consensus(snp_count=3)
        field = SnpField()
        assert field.format_value(consensus) == "snp=3"

    def test_format_value_zero(self):
        """Zero SNPs should return None (not displayed)."""
        consensus = make_consensus(snp_count=0)
        field = SnpField()
        assert field.format_value(consensus) is None

    def test_format_value_none(self):
        consensus = make_consensus(snp_count=None)
        field = SnpField()
        assert field.format_value(consensus) is None


class TestAmbigField:
    """Tests for the ambiguity code count field."""

    def test_format_value_with_ambiguity_codes(self):
        consensus = make_consensus(sequence="ACGTRYSWKM")  # 6 ambiguity codes
        field = AmbigField()
        assert field.format_value(consensus) == "ambig=6"

    def test_format_value_no_ambiguity(self):
        consensus = make_consensus(sequence="ACGTACGT")
        field = AmbigField()
        assert field.format_value(consensus) is None

    def test_format_value_lowercase(self):
        """Should count lowercase ambiguity codes too."""
        consensus = make_consensus(sequence="ACGTrysw")  # 4 ambiguity codes
        field = AmbigField()
        assert field.format_value(consensus) == "ambig=4"


class TestPrimersField:
    """Tests for the primers field."""

    def test_format_value_single_primer(self):
        consensus = make_consensus(primers=["ITS1"])
        field = PrimersField()
        assert field.format_value(consensus) == "primers=ITS1"

    def test_format_value_multiple_primers(self):
        consensus = make_consensus(primers=["ITS1", "ITS4"])
        field = PrimersField()
        assert field.format_value(consensus) == "primers=ITS1,ITS4"

    def test_format_value_none(self):
        consensus = make_consensus(primers=None)
        field = PrimersField()
        assert field.format_value(consensus) is None

    def test_format_value_empty_list(self):
        consensus = make_consensus(primers=[])
        field = PrimersField()
        assert field.format_value(consensus) is None


class TestRidField:
    """Tests for the mean read identity field."""

    def test_format_value(self):
        consensus = make_consensus(rid=0.95)
        field = RidField()
        assert field.format_value(consensus) == "rid=95.0"

    def test_format_value_decimal_precision(self):
        consensus = make_consensus(rid=0.9567)
        field = RidField()
        assert field.format_value(consensus) == "rid=95.7"

    def test_format_value_none(self):
        consensus = make_consensus(rid=None)
        field = RidField()
        assert field.format_value(consensus) is None


class TestRidMinField:
    """Tests for the minimum read identity field."""

    def test_format_value(self):
        consensus = make_consensus(rid_min=0.88)
        field = RidMinField()
        assert field.format_value(consensus) == "rid_min=88.0"

    def test_format_value_none(self):
        consensus = make_consensus(rid_min=None)
        field = RidMinField()
        assert field.format_value(consensus) is None


class TestGroupField:
    """Tests for the group number field."""

    def test_format_value_standard_name(self):
        consensus = make_consensus(sample_name="specimen-1.v1")
        field = GroupField()
        assert field.format_value(consensus) == "group=1"

    def test_format_value_multi_digit_group(self):
        consensus = make_consensus(sample_name="specimen-12.v3")
        field = GroupField()
        assert field.format_value(consensus) == "group=12"

    def test_format_value_raw_variant(self):
        """Should extract group from .raw suffixed names."""
        consensus = make_consensus(sample_name="specimen-2.v1.raw1")
        field = GroupField()
        assert field.format_value(consensus) == "group=2"

    def test_format_value_no_match(self):
        """Should return None for non-summarize naming."""
        consensus = make_consensus(sample_name="specimen-c1")
        field = GroupField()
        assert field.format_value(consensus) is None


class TestVariantField:
    """Tests for the variant identifier field."""

    def test_format_value_standard_name(self):
        consensus = make_consensus(sample_name="specimen-1.v1")
        field = VariantField()
        assert field.format_value(consensus) == "variant=v1"

    def test_format_value_higher_variant(self):
        consensus = make_consensus(sample_name="specimen-1.v5")
        field = VariantField()
        assert field.format_value(consensus) == "variant=v5"

    def test_format_value_raw_variant(self):
        """Should extract variant from .raw suffixed names."""
        consensus = make_consensus(sample_name="specimen-2.v1.raw1")
        field = VariantField()
        assert field.format_value(consensus) == "variant=v1"

    def test_format_value_no_match(self):
        consensus = make_consensus(sample_name="specimen-c1")
        field = VariantField()
        assert field.format_value(consensus) is None


class TestFieldRegistry:
    """Tests for the field registry."""

    def test_all_fields_registered(self):
        """All expected fields should be in the registry."""
        expected = {
            'size', 'ric', 'length', 'rawric', 'rawlen',
            'snp', 'ambig', 'rid', 'rid_min', 'primers',
            'group', 'variant'
        }
        assert set(FASTA_FIELDS.keys()) == expected

    def test_fields_have_correct_names(self):
        """Each field's name attribute should match its registry key."""
        for key, field in FASTA_FIELDS.items():
            assert field.name == key, f"Field '{key}' has name '{field.name}'"


class TestFieldPresets:
    """Tests for field presets."""

    def test_default_preset(self):
        assert 'size' in FASTA_FIELD_PRESETS['default']
        assert 'ric' in FASTA_FIELD_PRESETS['default']
        assert 'rawric' in FASTA_FIELD_PRESETS['default']

    def test_minimal_preset(self):
        assert FASTA_FIELD_PRESETS['minimal'] == ['size', 'ric']

    def test_id_only_preset(self):
        """id-only preset should be empty (just sample name)."""
        assert FASTA_FIELD_PRESETS['id-only'] == []

    def test_all_presets_reference_valid_fields(self):
        """All fields referenced in presets should exist in registry."""
        for preset_name, field_names in FASTA_FIELD_PRESETS.items():
            for field_name in field_names:
                assert field_name in FASTA_FIELDS, \
                    f"Preset '{preset_name}' references unknown field '{field_name}'"


class TestParseFastaFields:
    """Tests for parse_fasta_fields function."""

    def test_parse_preset(self):
        fields = parse_fasta_fields("minimal")
        assert len(fields) == 2
        assert fields[0].name == 'size'
        assert fields[1].name == 'ric'

    def test_parse_field_list(self):
        fields = parse_fasta_fields("size,length,primers")
        assert len(fields) == 3
        assert [f.name for f in fields] == ['size', 'length', 'primers']

    def test_parse_preset_union(self):
        """Multiple presets should union their fields."""
        fields = parse_fasta_fields("minimal,qc")
        names = [f.name for f in fields]
        # Should have minimal's fields plus qc's additions
        assert 'size' in names
        assert 'ric' in names
        assert 'length' in names
        assert 'rid' in names
        assert 'ambig' in names

    def test_parse_preset_plus_field(self):
        """Can add individual fields to a preset."""
        fields = parse_fasta_fields("minimal,primers")
        assert len(fields) == 3
        names = [f.name for f in fields]
        assert names == ['size', 'ric', 'primers']

    def test_parse_removes_duplicates(self):
        """Duplicate fields should only appear once."""
        fields = parse_fasta_fields("size,ric,size")
        assert len(fields) == 2

    def test_parse_empty_string_defaults(self):
        """Empty string should default to 'default' preset."""
        fields = parse_fasta_fields("")
        default_fields = parse_fasta_fields("default")
        assert [f.name for f in fields] == [f.name for f in default_fields]

    def test_parse_case_insensitive(self):
        fields = parse_fasta_fields("SIZE,RIC")
        assert len(fields) == 2
        assert fields[0].name == 'size'

    def test_parse_unknown_raises_error(self):
        with pytest.raises(ValueError) as exc_info:
            parse_fasta_fields("size,unknown_field")
        assert "unknown_field" in str(exc_info.value)
        assert "Available" in str(exc_info.value)

    def test_parse_id_only_returns_empty(self):
        fields = parse_fasta_fields("id-only")
        assert fields == []


class TestFormatFastaHeader:
    """Tests for format_fasta_header function."""

    def test_basic_format(self):
        consensus = make_consensus(sample_name="specimen-1.v1", size=100, ric=50)
        fields = parse_fasta_fields("size,ric")
        header = format_fasta_header(consensus, fields)
        assert header == "specimen-1.v1 size=100 ric=50"

    def test_skips_none_fields(self):
        """Fields that return None should not appear in header."""
        consensus = make_consensus(
            sample_name="specimen-1.v1",
            size=100,
            ric=50,
            raw_ric=None,  # Will return None
            snp_count=0,   # Will return None (zero = not displayed)
        )
        fields = parse_fasta_fields("size,ric,rawric,snp")
        header = format_fasta_header(consensus, fields)
        assert header == "specimen-1.v1 size=100 ric=50"
        assert "rawric" not in header
        assert "snp" not in header

    def test_includes_populated_optional_fields(self):
        """Fields with values should appear in header."""
        consensus = make_consensus(
            sample_name="specimen-1.v1",
            size=100,
            ric=50,
            raw_ric=[30, 20],
            snp_count=2,
        )
        fields = parse_fasta_fields("size,ric,rawric,snp")
        header = format_fasta_header(consensus, fields)
        assert header == "specimen-1.v1 size=100 ric=50 rawric=30+20 snp=2"

    def test_id_only(self):
        """Empty field list should just return sample name."""
        consensus = make_consensus(sample_name="specimen-1.v1")
        fields = parse_fasta_fields("id-only")
        header = format_fasta_header(consensus, fields)
        assert header == "specimen-1.v1"

    def test_full_preset(self):
        """Full preset should include all available fields."""
        consensus = make_consensus(
            sample_name="specimen-1.v1",
            sequence="ACGTRY",  # 2 ambig
            size=100,
            ric=50,
            raw_ric=[30, 20],
            raw_len=[600, 580],
            snp_count=1,
            primers=["ITS1", "ITS4"],
            rid=0.95,
        )
        fields = parse_fasta_fields("full")
        header = format_fasta_header(consensus, fields)

        assert "size=100" in header
        assert "ric=50" in header
        assert "length=6" in header
        assert "rawric=30+20" in header
        assert "rawlen=600+580" in header
        assert "snp=1" in header
        assert "ambig=2" in header
        assert "rid=95.0" in header
        assert "primers=ITS1,ITS4" in header


class TestRawRicIntegration:
    """Integration tests specifically for rawric field population.

    This addresses an issue where rawric was not populating correctly.
    """

    def test_rawric_from_single_merge(self):
        """rawric should capture RiC values when sequences are merged."""
        consensus = make_consensus(
            sample_name="specimen-1.v1",
            ric=80,  # Total after merge
            size=100,
            raw_ric=[50, 30],  # Two variants merged
        )
        field = RawRicField()
        result = field.format_value(consensus)
        assert result == "rawric=50+30"

    def test_rawric_in_header_after_merge(self):
        """rawric should appear in FASTA header for merged sequences."""
        consensus = make_consensus(
            sample_name="specimen-1.v1",
            ric=80,
            size=100,
            raw_ric=[50, 30],
        )
        fields = parse_fasta_fields("default")
        header = format_fasta_header(consensus, fields)
        assert "rawric=50+30" in header

    def test_rawric_not_in_header_for_unmerged(self):
        """rawric should not appear for sequences that weren't merged."""
        consensus = make_consensus(
            sample_name="specimen-1.v1",
            ric=50,
            size=60,
            raw_ric=None,
        )
        fields = parse_fasta_fields("default")
        header = format_fasta_header(consensus, fields)
        assert "rawric" not in header

    def test_rawric_many_values(self):
        """rawric should handle many merged variants."""
        consensus = make_consensus(
            raw_ric=[100, 80, 50, 30, 20, 10, 5],
        )
        field = RawRicField()
        result = field.format_value(consensus)
        assert result == "rawric=100+80+50+30+20+10+5"
