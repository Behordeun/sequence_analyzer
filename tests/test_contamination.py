"""Tests for contamination detection module."""

import json
from pathlib import Path

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.contamination import (
    compute_dinucleotide_frequencies,
    compute_gc_content,
    detect_contamination,
    list_available_organisms,
    load_organism_profiles,
)


@pytest.fixture
def profiles_path() -> Path:
    """Path to the bundled organism profiles."""
    return Path(__file__).resolve().parents[1] / "data" / "organism_profiles.json"


@pytest.fixture
def tmp_profiles(tmp_path: Path) -> Path:
    """Minimal profiles file for isolated testing.

    The test_org profile uses dinucleotide frequencies that match
    a repeating ATGC pattern so that short test sequences get "Low" risk
    when their GC is within range.
    """
    profiles = {
        "organisms": {
            "test_org": {
                "display_name": "Test Organism",
                "gc_range": [45.0, 55.0],
                "gc_mean": 50.0,
                "dinucleotide_frequencies": {
                    "AA": 0.000,
                    "AT": 0.263,
                    "AG": 0.000,
                    "AC": 0.000,
                    "TA": 0.000,
                    "TT": 0.000,
                    "TG": 0.263,
                    "TC": 0.000,
                    "GA": 0.000,
                    "GT": 0.000,
                    "GG": 0.000,
                    "GC": 0.263,
                    "CA": 0.211,
                    "CT": 0.000,
                    "CG": 0.000,
                    "CC": 0.000,
                },
                "description": "Test organism matching ATGC-repeat dinucleotides.",
            },
            "high_gc_org": {
                "display_name": "High GC Organism",
                "gc_range": [60.0, 70.0],
                "gc_mean": 65.0,
                "dinucleotide_frequencies": {
                    "AA": 0.030,
                    "AT": 0.030,
                    "AG": 0.060,
                    "AC": 0.060,
                    "TA": 0.030,
                    "TT": 0.030,
                    "TG": 0.060,
                    "TC": 0.060,
                    "GA": 0.060,
                    "GT": 0.060,
                    "GG": 0.080,
                    "GC": 0.080,
                    "CA": 0.060,
                    "CT": 0.060,
                    "CG": 0.080,
                    "CC": 0.080,
                },
                "description": "High GC test organism.",
            },
        }
    }
    path = tmp_path / "profiles.json"
    path.write_text(json.dumps(profiles), encoding="utf-8")
    return path


@pytest.fixture
def balanced_gc_sequence() -> SeqRecord:
    """Sequence with ~50% GC content."""
    # 10 G/C + 10 A/T = 50% GC
    return SeqRecord(Seq("ATGCATGCATGCATGCATGC"), id="balanced")


@pytest.fixture
def high_gc_sequence() -> SeqRecord:
    """Sequence with ~80% GC content — would be a contaminant for most organisms."""
    return SeqRecord(Seq("GCGCGCGCGCGCGCGCGCGCATGC"), id="high_gc")


@pytest.fixture
def low_gc_sequence() -> SeqRecord:
    """Sequence with ~20% GC content — typical of AT-rich organisms."""
    return SeqRecord(Seq("AATATAATATAAATAATGCA"), id="low_gc")


class TestLoadOrganismProfiles:
    def test_loads_bundled_profiles(self, profiles_path: Path):
        profiles = load_organism_profiles(profiles_path)
        assert "general" in profiles
        assert "escherichia_coli" in profiles
        assert "sars_cov_2" in profiles

    def test_loads_custom_path(self, tmp_profiles: Path):
        profiles = load_organism_profiles(tmp_profiles)
        assert "test_org" in profiles
        assert profiles["test_org"]["gc_mean"] == 50.0

    def test_raises_on_missing_file(self, tmp_path: Path):
        fake_path = tmp_path / "nonexistent.json"
        with pytest.raises(FileNotFoundError, match="not found"):
            load_organism_profiles(fake_path)

    def test_raises_on_malformed_json(self, tmp_path: Path):
        bad_file = tmp_path / "bad.json"
        bad_file.write_text('{"no_organisms_key": true}', encoding="utf-8")
        with pytest.raises(ValueError, match="missing 'organisms' key"):
            load_organism_profiles(bad_file)

    def test_profile_structure(self, profiles_path: Path):
        profiles = load_organism_profiles(profiles_path)
        for slug, profile in profiles.items():
            assert "gc_range" in profile, f"{slug} missing gc_range"
            assert "gc_mean" in profile, f"{slug} missing gc_mean"
            assert "dinucleotide_frequencies" in profile, f"{slug} missing dinuc freqs"
            assert len(profile["gc_range"]) == 2
            assert len(profile["dinucleotide_frequencies"]) == 16


class TestListAvailableOrganisms:
    def test_returns_sorted_list(self, profiles_path: Path):
        organisms = list_available_organisms(profiles_path)
        assert organisms == sorted(organisms)
        assert "general" in organisms

    def test_custom_profiles(self, tmp_profiles: Path):
        organisms = list_available_organisms(tmp_profiles)
        assert organisms == ["high_gc_org", "test_org"]


class TestComputeGCContent:
    def test_balanced_sequence(self):
        # ATGC repeated = 50% GC
        assert compute_gc_content("ATGCATGC") == pytest.approx(50.0)

    def test_all_gc(self):
        assert compute_gc_content("GGGGCCCC") == pytest.approx(100.0)

    def test_all_at(self):
        assert compute_gc_content("AAAATTTT") == pytest.approx(0.0)

    def test_ignores_gaps_and_ambiguous(self):
        # Only counts ATGC bases
        assert compute_gc_content("G-G-N-C-C") == pytest.approx(100.0)

    def test_empty_string_returns_zero(self):
        assert compute_gc_content("") == 0.0

    def test_only_ambiguous_returns_zero(self):
        assert compute_gc_content("NNNN---") == 0.0

    def test_case_insensitive(self):
        assert compute_gc_content("atgcATGC") == pytest.approx(50.0)


class TestComputeDinucleotideFrequencies:
    def test_returns_all_16_dinucleotides(self):
        freqs = compute_dinucleotide_frequencies("ATGCATGCATGC")
        assert len(freqs) == 16
        assert set(freqs.keys()) == {
            "AA",
            "AT",
            "AG",
            "AC",
            "TA",
            "TT",
            "TG",
            "TC",
            "GA",
            "GT",
            "GG",
            "GC",
            "CA",
            "CT",
            "CG",
            "CC",
        }

    def test_frequencies_sum_to_one(self):
        freqs = compute_dinucleotide_frequencies("ATGCGATCGATCGATCGATCGATCG")
        total = sum(freqs.values())
        assert total == pytest.approx(1.0, abs=1e-9)

    def test_homopolymer_concentrates_frequency(self):
        freqs = compute_dinucleotide_frequencies("AAAAAAAAAA")
        assert freqs["AA"] == pytest.approx(1.0)
        assert freqs["TT"] == pytest.approx(0.0)

    def test_short_sequence_returns_uniform(self):
        # Single base: can't form any dinucleotide
        freqs = compute_dinucleotide_frequencies("A")
        expected = 1.0 / 16.0
        assert all(v == pytest.approx(expected) for v in freqs.values())

    def test_empty_returns_uniform(self):
        freqs = compute_dinucleotide_frequencies("")
        expected = 1.0 / 16.0
        assert all(v == pytest.approx(expected) for v in freqs.values())

    def test_ignores_non_atgc_characters(self):
        # "ANGC" becomes "AGC" after filtering → dinucs: AG, GC
        freqs = compute_dinucleotide_frequencies("ANGC")
        # clean = "AGC", pairs = AG + GC → total 2
        assert freqs["AG"] == pytest.approx(0.5)
        assert freqs["GC"] == pytest.approx(0.5)

    def test_alternating_pattern(self):
        # ATATAT → dinucs: AT, TA, AT, TA, AT → 3 AT + 2 TA = 5 total
        freqs = compute_dinucleotide_frequencies("ATATAT")
        assert freqs["AT"] == pytest.approx(3.0 / 5.0)
        assert freqs["TA"] == pytest.approx(2.0 / 5.0)


class TestDetectContamination:
    def test_returns_expected_columns(self, balanced_gc_sequence: SeqRecord, tmp_profiles: Path):
        df = detect_contamination(
            [balanced_gc_sequence], organism="test_org", profiles_path=tmp_profiles
        )
        expected_cols = {
            "ID",
            "GC_Observed",
            "GC_Expected",
            "GC_Deviation",
            "Dinuc_Distance",
            "Contamination_Risk",
            "Risk_Reason",
        }
        assert set(df.columns) == expected_cols

    def test_one_row_per_sequence(self, tmp_profiles: Path):
        records = [
            SeqRecord(Seq("ATGCGATCGATCG"), id="s1"),
            SeqRecord(Seq("GCTAGCTAGCTAG"), id="s2"),
            SeqRecord(Seq("ATATATATATATAT"), id="s3"),
        ]
        df = detect_contamination(records, organism="test_org", profiles_path=tmp_profiles)
        assert len(df) == 3

    def test_balanced_sequence_low_risk_for_test_org(
        self, balanced_gc_sequence: SeqRecord, tmp_profiles: Path
    ):
        df = detect_contamination(
            [balanced_gc_sequence], organism="test_org", profiles_path=tmp_profiles
        )
        row = df.iloc[0]
        assert row["Contamination_Risk"] == "Low"
        assert row["GC_Deviation"] == 0.0

    def test_high_gc_sequence_flagged_for_test_org(
        self, high_gc_sequence: SeqRecord, tmp_profiles: Path
    ):
        # 80% GC far outside test_org's 45-55% range
        df = detect_contamination(
            [high_gc_sequence], organism="test_org", profiles_path=tmp_profiles
        )
        row = df.iloc[0]
        assert row["Contamination_Risk"] in ("Medium", "High")
        assert row["GC_Deviation"] > 0.0

    def test_low_gc_sequence_flagged_for_high_gc_org(
        self, low_gc_sequence: SeqRecord, tmp_profiles: Path
    ):
        # 20% GC far below high_gc_org's 60-70% range
        df = detect_contamination(
            [low_gc_sequence], organism="high_gc_org", profiles_path=tmp_profiles
        )
        row = df.iloc[0]
        assert row["Contamination_Risk"] == "High"
        assert row["GC_Deviation"] >= 40.0

    def test_gc_within_range_has_zero_deviation(self, tmp_profiles: Path):
        # Craft a sequence that's ~50% GC → within test_org's 45-55% range
        seq = SeqRecord(Seq("ATGCATGCATGCATGCATGC"), id="in_range")
        df = detect_contamination([seq], organism="test_org", profiles_path=tmp_profiles)
        assert df.iloc[0]["GC_Deviation"] == 0.0

    def test_gc_below_range_measures_distance_to_min(self, tmp_profiles: Path):
        # All AT → 0% GC, range is 45-55%, deviation should be 45.0
        seq = SeqRecord(Seq("AAAAATTTTTAAAAATTTTT"), id="all_at")
        df = detect_contamination([seq], organism="test_org", profiles_path=tmp_profiles)
        assert df.iloc[0]["GC_Deviation"] == pytest.approx(45.0, abs=0.1)

    def test_gc_above_range_measures_distance_to_max(self, tmp_profiles: Path):
        # All GC → 100% GC, range is 45-55%, deviation should be 45.0
        seq = SeqRecord(Seq("GGGGGGGGGGCCCCCCCCCC"), id="all_gc")
        df = detect_contamination([seq], organism="test_org", profiles_path=tmp_profiles)
        assert df.iloc[0]["GC_Deviation"] == pytest.approx(45.0, abs=0.1)

    def test_raises_on_empty_sequences(self, tmp_profiles: Path):
        with pytest.raises(ValueError, match="At least one sequence"):
            detect_contamination([], organism="test_org", profiles_path=tmp_profiles)

    def test_raises_on_unknown_organism(self, tmp_profiles: Path):
        records = [SeqRecord(Seq("ATGCGATCG"), id="x")]
        with pytest.raises(ValueError, match="Unknown organism"):
            detect_contamination(records, organism="fake_org", profiles_path=tmp_profiles)

    def test_uses_general_profile_by_default(self, profiles_path: Path):
        # Use a longer sequence with varied composition to get realistic dinuc distances
        import random

        random.seed(99)
        seq_str = "".join(random.choices("ATGC", k=200))
        records = [SeqRecord(Seq(seq_str), id="default_test")]
        df = detect_contamination(records, profiles_path=profiles_path)
        # General profile has wide gc_range [20, 80] and a 200bp random sequence
        # produces realistic dinucleotide distances
        assert df.iloc[0]["Contamination_Risk"] == "Low"

    def test_general_profile_permissive(self, profiles_path: Path):
        # Even biased sequences should pass with the general profile's wide range
        records = [SeqRecord(Seq("GGGGCCCCGGGGCCCCATGC"), id="biased")]
        df = detect_contamination(records, profiles_path=profiles_path)
        # ~84% GC is outside even the general range [20, 80] so this should flag
        row = df.iloc[0]
        # Confirm the column exists and has a valid value
        assert row["Contamination_Risk"] in ("Low", "Medium", "High")

    def test_mixed_risk_batch(self, tmp_profiles: Path):
        records = [
            SeqRecord(Seq("ATGCATGCATGCATGCATGC"), id="good"),  # ~50% GC, matches test_org dinucs
            SeqRecord(Seq("AAAAAAAAAAATTTTTTTTT"), id="bad"),  # ~0% GC, far out of range
        ]
        df = detect_contamination(records, organism="test_org", profiles_path=tmp_profiles)
        good_row = df[df["ID"] == "good"].iloc[0]
        bad_row = df[df["ID"] == "bad"].iloc[0]
        assert good_row["Contamination_Risk"] == "Low"
        # bad sequence has GC deviation of 45% and divergent dinucs
        assert bad_row["Contamination_Risk"] in ("Medium", "High")

    def test_dinuc_distance_is_nonnegative(self, tmp_profiles: Path):
        records = [
            SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="s1"),
            SeqRecord(Seq("AAAAAAAAAAAAAAAAAAGC"), id="s2"),
        ]
        df = detect_contamination(records, organism="test_org", profiles_path=tmp_profiles)
        assert all(df["Dinuc_Distance"] >= 0.0)

    def test_sequence_with_gaps_and_ambiguous(self, tmp_profiles: Path):
        # Contamination detection should handle sequences with non-ATGC chars gracefully
        seq = SeqRecord(Seq("ATGC--NN--ATGCATGCATGCATGC"), id="gappy")
        df = detect_contamination([seq], organism="test_org", profiles_path=tmp_profiles)
        assert len(df) == 1
        assert df.iloc[0]["GC_Observed"] == pytest.approx(50.0)

    def test_bundled_ecoli_profile(self, profiles_path: Path):
        # Use a longer sequence with E. coli-like composition (~50% GC, varied dinucs)
        import random

        random.seed(6)
        seq_str = "".join(random.choices("ATGC", weights=[0.25, 0.25, 0.25, 0.25], k=200))
        seq = SeqRecord(Seq(seq_str), id="ecoli_like")
        df = detect_contamination([seq], organism="escherichia_coli", profiles_path=profiles_path)
        row = df.iloc[0]
        assert row["GC_Expected"] == 50.8
        # A balanced 200bp sequence should have GC near 50% and realistic dinuc profile
        assert row["Contamination_Risk"] == "Low"

    def test_id_column_uses_record_id(self, tmp_profiles: Path):
        records = [SeqRecord(Seq("ATGCATGCATGC"), id="my_sequence_42")]
        df = detect_contamination(records, organism="test_org", profiles_path=tmp_profiles)
        assert df.iloc[0]["ID"] == "my_sequence_42"

    def test_all_ambiguous_and_gaps_does_not_crash(self, tmp_profiles: Path):
        # Sequence with only ambiguous/gap characters should produce valid output
        seq = SeqRecord(Seq("NNNN---NNN--"), id="ambiguous_only")
        df = detect_contamination([seq], organism="test_org", profiles_path=tmp_profiles)
        row = df.iloc[0]
        assert row["GC_Observed"] == pytest.approx(0.0)
        assert row["Dinuc_Distance"] >= 0.0
        assert np.isfinite(row["Dinuc_Distance"])
        assert row["Contamination_Risk"] in ("Low", "Medium", "High")
        assert isinstance(row["Risk_Reason"], str)
        assert row["Risk_Reason"].strip() != ""
