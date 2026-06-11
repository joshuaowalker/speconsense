"""Tests for the profile system."""

import argparse
import os
import tempfile
from pathlib import Path
from unittest.mock import patch
import pytest

from speconsense.profiles import (
    Profile,
    ProfileError,
    ProfileVersionError,
    ProfileValidationError,
    RESERVED_PROFILE_NAMES,
    check_version_compatible,
    list_profiles,
    list_bundled_profiles,
    apply_profile_to_args,
    ensure_user_profiles_dir,
    get_bundled_profile_path,
    PROFILES_DIR,
)


class TestVersionCompatibility:
    """Tests for version pattern matching."""

    def test_exact_match(self):
        assert check_version_compatible("0.7.0", "0.7.0") is True
        assert check_version_compatible("0.7.0", "0.7.1") is False

    def test_patch_wildcard(self):
        assert check_version_compatible("0.7.*", "0.7.0") is True
        assert check_version_compatible("0.7.*", "0.7.1") is True
        assert check_version_compatible("0.7.*", "0.7.99") is True
        assert check_version_compatible("0.7.*", "0.8.0") is False
        assert check_version_compatible("0.7.*", "1.7.0") is False

    def test_minor_wildcard(self):
        assert check_version_compatible("0.*", "0.7.0") is True
        assert check_version_compatible("0.*", "0.99.0") is True
        assert check_version_compatible("0.*", "1.0.0") is False

    def test_universal_wildcard(self):
        assert check_version_compatible("*", "0.7.0") is True
        assert check_version_compatible("*", "1.0.0") is True
        assert check_version_compatible("*", "99.99.99") is True

    def test_dev_version(self):
        assert check_version_compatible("*", "dev") is True
        assert check_version_compatible("0.*", "dev") is False


class TestBundledProfiles:
    """Tests for bundled profile access."""

    def test_list_bundled_profiles(self):
        profiles = list_bundled_profiles()
        assert "herbarium" in profiles
        assert "strict" in profiles
        assert "example" in profiles
        assert "nostalgia" in profiles

    def test_get_bundled_profile_path(self):
        path = get_bundled_profile_path("herbarium")
        assert path is not None
        assert path.exists()
        assert path.name == "herbarium.yaml"

    def test_get_bundled_profile_path_not_found(self):
        path = get_bundled_profile_path("nonexistent")
        assert path is None

    def test_load_bundled_profile(self):
        # Load without version check since test may run on different version
        profile = Profile.load("herbarium", check_version=False)
        assert profile.name == "herbarium"
        assert profile.description != ""
        assert "min-size" in profile.speconsense
        assert "min-ric" in profile.speconsense_summarize


class TestProfileLoading:
    """Tests for profile loading behavior."""

    def test_load_nonexistent_profile(self):
        with pytest.raises(ProfileError) as exc_info:
            Profile.load("nonexistent_profile_xyz", check_version=False)
        assert "not found" in str(exc_info.value)

    def test_load_with_version_check(self):
        # Bundled profiles' speconsense-version field tracks the current major
        # version line, so the version check should pass for the running build.
        profile = Profile.load("herbarium", check_version=True)
        assert profile.name == "herbarium"

    def test_user_profile_takes_precedence(self):
        """User profile should be loaded instead of bundled if both exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            # Create a user profile with same name as bundled
            user_profile = user_dir / "herbarium.yaml"
            user_profile.write_text("""
speconsense-version: "*"
description: "User override"
speconsense:
  min-size: 999
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                profile = Profile.load("herbarium", check_version=False)
                assert profile.description == "User override"
                assert profile.speconsense.get("min-size") == 999


class TestVersionValidation:
    """Tests for version compatibility errors."""

    def test_incompatible_version_raises_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            # Create profile requiring old version
            old_profile = user_dir / "old.yaml"
            old_profile.write_text("""
speconsense-version: "0.5.*"
description: "Old profile"
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                with patch("speconsense.profiles.__version__", "0.7.0"):
                    with pytest.raises(ProfileVersionError) as exc_info:
                        Profile.load("old", check_version=True)
                    assert "0.5.*" in str(exc_info.value)
                    assert "0.7.0" in str(exc_info.value)


class TestKeyValidation:
    """Tests for profile key validation."""

    def test_valid_keys_pass(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            valid_profile = user_dir / "valid.yaml"
            valid_profile.write_text("""
speconsense-version: "*"
description: "Valid profile"
speconsense:
  min-size: 10
speconsense-summarize:
  min-ric: 5
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                profile = Profile.load("valid", check_version=False)
                assert profile.speconsense["min-size"] == 10

    def test_unknown_speconsense_key_raises_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            bad_profile = user_dir / "bad.yaml"
            bad_profile.write_text("""
speconsense-version: "*"
description: "Bad profile"
speconsense:
  min-size: 10
  typo-option: 5
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                with pytest.raises(ProfileValidationError) as exc_info:
                    Profile.load("bad", check_version=False)
                assert "typo-option" in str(exc_info.value)

    def test_unknown_summarize_key_raises_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            bad_profile = user_dir / "bad.yaml"
            bad_profile.write_text("""
speconsense-version: "*"
description: "Bad profile"
speconsense-summarize:
  min-ric: 5
  invalid-key: 10
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                with pytest.raises(ProfileValidationError) as exc_info:
                    Profile.load("bad", check_version=False)
                assert "invalid-key" in str(exc_info.value)

    def test_hp_normalization_length_key_accepted_in_both_sections(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            (user_dir / "hp.yaml").write_text("""
speconsense-version: "*"
description: "HP threshold profile"
speconsense:
  hp-normalization-length: 6
speconsense-summarize:
  hp-normalization-length: 4
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                profile = Profile.load("hp", check_version=False)
                assert profile.speconsense["hp-normalization-length"] == 6
                assert profile.speconsense_summarize["hp-normalization-length"] == 4

    def test_legacy_hp_min_length_key_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            (user_dir / "legacy.yaml").write_text("""
speconsense-version: "*"
description: "Legacy HP key profile"
speconsense:
  hp-min-length: 6
""")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                with pytest.raises(ProfileValidationError) as exc_info:
                    Profile.load("legacy", check_version=False)
                assert "hp-min-length" in str(exc_info.value)


class TestDefaultProfile:
    """Tests for the reserved 'default' profile name."""

    def test_load_default_returns_empty_profile(self):
        profile = Profile.load("default")
        assert profile.name == "default"
        assert profile.speconsense == {}
        assert profile.speconsense_summarize == {}
        assert profile.version == "*"

    def test_load_default_does_not_require_yaml(self):
        """-p default works even without PyYAML installed."""
        with patch("speconsense.profiles.yaml", None):
            profile = Profile.load("default")
            assert profile.name == "default"
            assert profile.speconsense == {}

    def test_default_is_reserved(self):
        assert "default" in RESERVED_PROFILE_NAMES

    def test_default_yaml_file_is_shadowed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()
            (user_dir / "default.yaml").write_text(
                'speconsense-version: "*"\nspeconsense:\n  min-size: 99\n'
            )
            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                profile = Profile.load("default")
                assert profile.speconsense == {}


class TestProfileApplication:
    """Tests for applying profiles to argparse args."""

    def test_apply_profile_values(self):
        profile = Profile(
            name="test",
            version="*",
            description="Test",
            speconsense={"min-size": 10, "min-identity": 0.85},
            speconsense_summarize={},
        )

        args = argparse.Namespace(min_size=5, min_identity=0.9, other=True)
        explicit_args = set()  # No explicit args

        apply_profile_to_args(args, profile, "speconsense", explicit_args)

        assert args.min_size == 10
        assert args.min_identity == 0.85
        assert args.other is True  # Unchanged

    def test_explicit_args_override_profile(self):
        profile = Profile(
            name="test",
            version="*",
            description="Test",
            speconsense={"min-size": 10, "min-identity": 0.85},
            speconsense_summarize={},
        )

        args = argparse.Namespace(min_size=20, min_identity=0.9)
        explicit_args = {"min_size"}  # User explicitly set min_size

        apply_profile_to_args(args, profile, "speconsense", explicit_args)

        assert args.min_size == 20  # Not overridden
        assert args.min_identity == 0.85  # Overridden by profile

    def test_apply_summarize_profile(self):
        profile = Profile(
            name="test",
            version="*",
            description="Test",
            speconsense={},
            speconsense_summarize={"min-ric": 10, "select-max-variants": 3},
        )

        args = argparse.Namespace(min_ric=5, select_max_variants=-1)
        explicit_args = set()

        apply_profile_to_args(args, profile, "speconsense-summarize", explicit_args)

        assert args.min_ric == 10
        assert args.select_max_variants == 3


class TestListProfiles:
    """Tests for listing profiles."""

    def test_list_includes_bundled(self):
        profiles = list_profiles()
        assert "herbarium" in profiles
        assert "strict" in profiles
        assert "nostalgia" in profiles

    def test_list_includes_user_profiles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()

            (user_dir / "custom.yaml").write_text("speconsense-version: '*'\n")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                profiles = list_profiles()
                assert "custom" in profiles


class TestDetectProfileFromMetadata:
    """Tests for auto-detection of profile from core metadata."""

    def _write_metadata(self, base_dir, specimen, profile_value):
        """Write a minimal metadata JSON with the given profile value."""
        specimen_dir = os.path.join(base_dir, specimen, "cluster_debug")
        os.makedirs(specimen_dir, exist_ok=True)
        metadata = {"schema_version": "2.0", "profile": profile_value}
        with open(os.path.join(specimen_dir, f"{specimen}-metadata.json"), 'w') as f:
            import json
            json.dump(metadata, f)

    def test_unanimous_profile_detected(self):
        from speconsense.summarize.cli import _detect_profile_from_metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_metadata(tmpdir, "sp1", "compressed")
            self._write_metadata(tmpdir, "sp2", "compressed")
            assert _detect_profile_from_metadata(tmpdir) == "compressed"

    def test_mixed_profiles_returns_none(self, capsys):
        from speconsense.summarize.cli import _detect_profile_from_metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_metadata(tmpdir, "sp1", "compressed")
            self._write_metadata(tmpdir, "sp2", "strict")
            result = _detect_profile_from_metadata(tmpdir)
            assert result is None
            err = capsys.readouterr().err
            assert "different core profiles" in err

    def test_all_null_returns_none(self):
        from speconsense.summarize.cli import _detect_profile_from_metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_metadata(tmpdir, "sp1", None)
            self._write_metadata(tmpdir, "sp2", None)
            assert _detect_profile_from_metadata(tmpdir) is None

    def test_all_default_returns_none(self):
        from speconsense.summarize.cli import _detect_profile_from_metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_metadata(tmpdir, "sp1", "default")
            self._write_metadata(tmpdir, "sp2", "default")
            assert _detect_profile_from_metadata(tmpdir) is None

    def test_no_metadata_files_returns_none(self):
        from speconsense.summarize.cli import _detect_profile_from_metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            assert _detect_profile_from_metadata(tmpdir) is None

    def test_missing_profile_field_ignored(self):
        """Specimens without a profile field in metadata are ignored."""
        from speconsense.summarize.cli import _detect_profile_from_metadata
        import json
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_metadata(tmpdir, "sp1", "compressed")
            # Write a metadata JSON without a profile field
            sp2_dir = os.path.join(tmpdir, "sp2", "cluster_debug")
            os.makedirs(sp2_dir)
            with open(os.path.join(sp2_dir, "sp2-metadata.json"), 'w') as f:
                json.dump({"schema_version": "2.0"}, f)
            assert _detect_profile_from_metadata(tmpdir) == "compressed"

    def test_null_and_real_profile_detected(self):
        """Null profiles are ignored; unanimous real profiles are detected."""
        from speconsense.summarize.cli import _detect_profile_from_metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_metadata(tmpdir, "sp1", "compressed")
            self._write_metadata(tmpdir, "sp2", None)
            assert _detect_profile_from_metadata(tmpdir) == "compressed"


class TestUserDirectoryInit:
    """Tests for user profiles directory initialization."""

    def test_creates_directory_and_example(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "config" / "speconsense" / "profiles"

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                result = ensure_user_profiles_dir()

                assert result == user_dir
                assert user_dir.exists()
                assert (user_dir / "example.yaml").exists()

    def test_idempotent_when_example_exists(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()
            example = user_dir / "example.yaml"
            example.write_text("# Original content\n")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                ensure_user_profiles_dir()

                # Should not overwrite existing example
                assert example.read_text() == "# Original content\n"

    def test_skips_init_when_other_profiles_exist(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            user_dir = Path(tmpdir) / "profiles"
            user_dir.mkdir()
            (user_dir / "custom.yaml").write_text("speconsense-version: '*'\n")

            with patch("speconsense.profiles.PROFILES_DIR", user_dir):
                ensure_user_profiles_dir()

                # Should not create example when other profiles exist
                assert not (user_dir / "example.yaml").exists()
