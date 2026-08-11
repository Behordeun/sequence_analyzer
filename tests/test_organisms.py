"""Tests for organism-specific mode management."""

import json
from pathlib import Path

import pytest

from sequence_analyzer.core.organisms import (
    OrganismMode,
    get_organism_mode,
    list_organism_modes,
    load_organism_modes,
)


@pytest.fixture
def bundled_organisms_dir() -> Path:
    """Path to the bundled organism profiles directory."""
    return Path(__file__).resolve().parents[1] / "data" / "organisms"


@pytest.fixture
def tmp_organisms_dir(tmp_path: Path) -> Path:
    """Create a temporary organisms directory with minimal profiles for testing."""
    org_dir = tmp_path / "organisms"
    org_dir.mkdir()

    general = {
        "display_name": "Test General",
        "description": "Permissive test defaults.",
        "seq_type": "auto",
        "qc_ambiguity_threshold": 5.0,
        "qc_min_length": 10,
        "qc_max_gap_fraction": 50.0,
        "contamination_profile": "general",
        "motif_patterns": ["ATG"],
        "reference_accession": "",
        "reference_description": "",
    }
    (org_dir / "general.json").write_text(json.dumps(general), encoding="utf-8")

    strict = {
        "display_name": "Strict Organism",
        "description": "Strict thresholds for testing.",
        "seq_type": "DNA",
        "qc_ambiguity_threshold": 1.0,
        "qc_min_length": 200,
        "qc_max_gap_fraction": 10.0,
        "contamination_profile": "escherichia_coli",
        "motif_patterns": ["ATG", "TATAAT", "AGGAGG"],
        "reference_accession": "NC_000001",
        "reference_description": "Test reference chromosome",
    }
    (org_dir / "strict_org.json").write_text(json.dumps(strict), encoding="utf-8")

    return org_dir


@pytest.fixture
def loaded_tmp_modes(tmp_organisms_dir: Path) -> dict:
    """Pre-loaded modes from the temporary organisms directory."""
    return load_organism_modes(tmp_organisms_dir)


class TestLoadOrganismModes:
    def test_loads_bundled_profiles(self, bundled_organisms_dir: Path):
        modes = load_organism_modes(bundled_organisms_dir)
        assert "general" in modes
        assert "escherichia_coli" in modes
        assert "sars_cov_2" in modes
        assert "plasmodium_falciparum" in modes

    def test_loads_default_path_without_arguments(self):
        modes = load_organism_modes()
        assert "general" in modes
        assert "escherichia_coli" in modes
        assert len(modes) >= 7

    def test_loads_custom_directory(self, loaded_tmp_modes: dict):
        assert "general" in loaded_tmp_modes
        assert "strict_org" in loaded_tmp_modes
        assert len(loaded_tmp_modes) == 2

    def test_returns_mode_instances(self, loaded_tmp_modes: dict):
        for mode in loaded_tmp_modes.values():
            assert isinstance(mode, OrganismMode)

    def test_empty_directory_returns_empty_dict(self, tmp_path: Path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        modes = load_organism_modes(empty_dir)
        assert modes == {}

    def test_nonexistent_directory_returns_empty_dict(self, tmp_path: Path):
        fake_dir = tmp_path / "nonexistent"
        modes = load_organism_modes(fake_dir)
        assert modes == {}

    def test_skips_malformed_json(self, tmp_path: Path):
        org_dir = tmp_path / "organisms_bad"
        org_dir.mkdir()
        (org_dir / "good.json").write_text(
            json.dumps({"display_name": "Good", "seq_type": "DNA"}), encoding="utf-8"
        )
        (org_dir / "bad.json").write_text("not valid json {{{", encoding="utf-8")
        modes = load_organism_modes(org_dir)
        assert "good" in modes
        assert "bad" not in modes


class TestListOrganismModes:
    def test_returns_sorted_list(self, bundled_organisms_dir: Path):
        slugs = list_organism_modes(bundled_organisms_dir)
        assert slugs == sorted(slugs)
        assert "general" in slugs

    def test_custom_directory(self, tmp_organisms_dir: Path):
        slugs = list_organism_modes(tmp_organisms_dir)
        assert slugs == ["general", "strict_org"]


class TestGetOrganismMode:
    def test_returns_correct_mode(self, tmp_organisms_dir: Path):
        mode = get_organism_mode("strict_org", tmp_organisms_dir)
        assert mode.slug == "strict_org"
        assert mode.display_name == "Strict Organism"
        assert mode.seq_type == "DNA"

    def test_raises_on_unknown_slug(self, tmp_organisms_dir: Path):
        with pytest.raises(ValueError, match="Unknown organism mode"):
            get_organism_mode("nonexistent", tmp_organisms_dir)

    def test_qc_thresholds_loaded(self, tmp_organisms_dir: Path):
        mode = get_organism_mode("strict_org", tmp_organisms_dir)
        assert mode.qc_ambiguity_threshold == 1.0
        assert mode.qc_min_length == 200
        assert mode.qc_max_gap_fraction == 10.0

    def test_contamination_profile_loaded(self, tmp_organisms_dir: Path):
        mode = get_organism_mode("strict_org", tmp_organisms_dir)
        assert mode.contamination_profile == "escherichia_coli"

    def test_motif_patterns_loaded(self, tmp_organisms_dir: Path):
        mode = get_organism_mode("strict_org", tmp_organisms_dir)
        assert mode.motif_patterns == ["ATG", "TATAAT", "AGGAGG"]

    def test_reference_accession_loaded(self, tmp_organisms_dir: Path):
        mode = get_organism_mode("strict_org", tmp_organisms_dir)
        assert mode.reference_accession == "NC_000001"
        assert mode.reference_description == "Test reference chromosome"


class TestOrganismModeDataclass:
    def test_is_frozen(self, tmp_organisms_dir: Path):
        mode = get_organism_mode("general", tmp_organisms_dir)
        with pytest.raises(AttributeError):
            mode.slug = "modified"  # type: ignore[misc]

    def test_defaults_for_missing_fields(self, tmp_path: Path):
        org_dir = tmp_path / "minimal"
        org_dir.mkdir()
        # Only display_name, everything else should get defaults
        (org_dir / "minimal.json").write_text(
            json.dumps({"display_name": "Minimal"}), encoding="utf-8"
        )
        modes = load_organism_modes(org_dir)
        mode = modes["minimal"]
        assert mode.seq_type == "DNA"
        assert mode.qc_ambiguity_threshold == 5.0
        assert mode.qc_min_length == 10
        assert mode.qc_max_gap_fraction == 50.0
        assert mode.contamination_profile == "general"
        assert mode.motif_patterns == []
        assert mode.reference_accession == ""


class TestBundledProfiles:
    """Validate the structure of all bundled organism profiles."""

    def test_all_profiles_have_required_fields(self, bundled_organisms_dir: Path):
        modes = load_organism_modes(bundled_organisms_dir)
        assert len(modes) >= 7
        for slug, mode in modes.items():
            assert mode.display_name, f"{slug} missing display_name"
            assert mode.seq_type in ("DNA", "RNA", "auto"), f"{slug} bad seq_type"
            assert mode.qc_ambiguity_threshold > 0, f"{slug} bad ambiguity threshold"
            assert mode.qc_min_length >= 1, f"{slug} bad min_length"

    def test_ecoli_has_promoter_motifs(self, bundled_organisms_dir: Path):
        mode = get_organism_mode("escherichia_coli", bundled_organisms_dir)
        # -10 and -35 promoter elements
        assert "TATAAT" in mode.motif_patterns
        assert "TTGACA" in mode.motif_patterns

    def test_sars_cov_2_is_rna(self, bundled_organisms_dir: Path):
        mode = get_organism_mode("sars_cov_2", bundled_organisms_dir)
        assert mode.seq_type == "RNA"
        assert mode.reference_accession == "MN908947"

    def test_plasmodium_has_lenient_thresholds(self, bundled_organisms_dir: Path):
        mode = get_organism_mode("plasmodium_falciparum", bundled_organisms_dir)
        # AT-rich genome needs higher ambiguity tolerance
        assert mode.qc_ambiguity_threshold > 5.0
        assert mode.qc_min_length < 100
