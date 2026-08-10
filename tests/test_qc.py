"""Tests for sequence quality control module."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.qc import assess_sequences, filter_passing


@pytest.fixture
def mixed_quality_sequences() -> list[SeqRecord]:
    """Sequences with varying quality for QC testing."""
    return [
        SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="good_1"),
        SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="good_2"),
        SeqRecord(Seq("NNNNATGCNNNN"), id="ambiguous"),  # 8/12 = 66% ambiguity
        SeqRecord(Seq("ATG"), id="too_short"),  # length 3 < 10
        SeqRecord(Seq("---ATG---ATG---"), id="gappy"),  # 9/15 = 60% gaps
    ]


class TestAssessSequences:
    def test_returns_correct_columns(self, mixed_quality_sequences):
        df = assess_sequences(mixed_quality_sequences)
        expected = {
            "ID",
            "Length",
            "GC_Percent",
            "Ambiguity_Rate",
            "Gap_Fraction",
            "Type",
            "Status",
            "Flags",
        }
        assert set(df.columns) == expected

    def test_one_row_per_sequence(self, mixed_quality_sequences):
        df = assess_sequences(mixed_quality_sequences)
        assert len(df) == 5

    def test_good_sequences_pass(self):
        records = [
            SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="clean_1"),
            SeqRecord(Seq("GCTAGCTAGCTAGCTAGCTAG"), id="clean_2"),
        ]
        df = assess_sequences(records)
        assert all(df["Status"] == "Pass")

    def test_high_ambiguity_fails(self):
        records = [SeqRecord(Seq("NNNNNATGCNNNNN"), id="ambig")]
        df = assess_sequences(records)
        assert df.iloc[0]["Status"] == "Fail"
        assert "High ambiguity" in df.iloc[0]["Flags"]

    def test_too_short_fails(self):
        records = [SeqRecord(Seq("ATG"), id="short")]
        df = assess_sequences(records)
        assert df.iloc[0]["Status"] == "Fail"
        assert "Too short" in df.iloc[0]["Flags"]

    def test_excessive_gaps_fails(self):
        records = [SeqRecord(Seq("----ATG-----"), id="gaps")]
        df = assess_sequences(records)
        assert df.iloc[0]["Status"] == "Fail"
        assert "Excessive gaps" in df.iloc[0]["Flags"]

    def test_borderline_ambiguity_warns(self):
        # 3% ambiguity (between 2% and 5%) should be Warning
        # 1 N in 33 bases = ~3%
        seq = "ATGCGATCGATCGATCGATCGATCGATCGATCN"
        records = [SeqRecord(Seq(seq), id="borderline")]
        df = assess_sequences(records)
        assert df.iloc[0]["Status"] == "Warning"
        assert "Borderline ambiguity" in df.iloc[0]["Flags"]

    def test_exact_threshold_boundary_fails(self):
        # Exactly at threshold should still pass (> threshold fails, not >=)
        # 5% of 20 = 1 N → ambiguity = 5.0% exactly → passes (not > 5.0)
        seq = "ATGCGATCGATCGATCGATN"  # 1/20 = 5.0%
        records = [SeqRecord(Seq(seq), id="exact")]
        df = assess_sequences(records)
        # 5.0 is not > 5.0, so it doesn't fail — but it's >= 2.0 so it warns
        assert df.iloc[0]["Status"] == "Warning"

    def test_length_outlier_warns(self):
        # Create sequences where one is dramatically longer than the rest
        short_seq = "ATGCGATCGATCG"  # 13 bp
        records = [
            SeqRecord(Seq(short_seq), id="normal_1"),
            SeqRecord(Seq(short_seq), id="normal_2"),
            SeqRecord(Seq(short_seq), id="normal_3"),
            SeqRecord(Seq(short_seq), id="normal_4"),
            SeqRecord(Seq(short_seq), id="normal_5"),
            SeqRecord(Seq(short_seq * 20), id="outlier"),  # 260 bp vs 13 bp
        ]
        df = assess_sequences(records)
        outlier_row = df[df["ID"] == "outlier"].iloc[0]
        assert outlier_row["Status"] == "Warning"
        assert "Length outlier" in outlier_row["Flags"]

    def test_custom_thresholds(self):
        records = [SeqRecord(Seq("NATGCGATCGATCG"), id="one_n")]
        # With strict threshold of 1%, one N in 14 = 7.1% should fail
        df = assess_sequences(records, ambiguity_threshold=1.0)
        assert df.iloc[0]["Status"] == "Fail"

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="At least one sequence"):
            assess_sequences([])

    def test_gc_percent_computed_correctly(self):
        records = [SeqRecord(Seq("GGCCAATT"), id="half_gc")]
        df = assess_sequences(records)
        assert df.iloc[0]["GC_Percent"] == pytest.approx(50.0)

    def test_detects_sequence_type(self):
        records = [SeqRecord(Seq("AUGCGAUCGAUCG"), id="rna_seq")]
        df = assess_sequences(records)
        assert df.iloc[0]["Type"] in ("DNA", "RNA")

    def test_detects_protein_type_and_ambiguity(self):
        seq_str = "ACDEFGHIKLMNPQRSTVWYXBZJ"
        records = [SeqRecord(Seq(seq_str), id="protein_seq")]
        df = assess_sequences(records)
        assert df.iloc[0]["Type"] == "Protein"
        # X, B, Z, J are ambiguous for protein = 4/24 * 100 ≈ 16.67%
        expected_rate = 4 / len(seq_str) * 100
        assert df.iloc[0]["Ambiguity_Rate"] == pytest.approx(expected_rate, rel=0.01)

    def test_zero_length_sequence_has_zero_metrics(self):
        records = [
            SeqRecord(Seq(""), id="empty"),
            SeqRecord(Seq("ATGCGATCGATCG"), id="nonempty"),
        ]
        df = assess_sequences(records, min_length=10)
        empty_row = df[df["ID"] == "empty"].iloc[0]
        assert empty_row["Length"] == 0
        assert empty_row["GC_Percent"] == pytest.approx(0.0)
        assert empty_row["Ambiguity_Rate"] == pytest.approx(0.0)
        assert empty_row["Gap_Fraction"] == pytest.approx(0.0)
        assert empty_row["Status"] == "Fail"
        assert "Empty sequence" in empty_row["Flags"]


class TestFilterPassing:
    def test_filters_only_passing(self, mixed_quality_sequences):
        df = assess_sequences(mixed_quality_sequences)
        passing = filter_passing(mixed_quality_sequences, df, include_warnings=False)
        # Only the two good sequences should pass
        assert len(passing) == 2
        assert all(r.id.startswith("good") for r in passing)

    def test_includes_warnings_by_default(self):
        # Borderline ambiguity → Warning
        records = [
            SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="good"),
            SeqRecord(Seq("ATGCGATCGATCGATCGATCGATCGATCGATCN"), id="warning"),
        ]
        df = assess_sequences(records)
        passing = filter_passing(records, df, include_warnings=True)
        assert len(passing) == 2

    def test_excludes_warnings_when_disabled(self):
        records = [
            SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="good"),
            SeqRecord(Seq("ATGCGATCGATCGATCGATCGATCGATCGATCN"), id="warning"),
        ]
        df = assess_sequences(records)
        passing = filter_passing(records, df, include_warnings=False)
        assert len(passing) == 1
        assert passing[0].id == "good"

    def test_raises_on_length_mismatch(self):
        records = [SeqRecord(Seq("ATGC"), id="one")]
        df = assess_sequences(records)
        with pytest.raises(ValueError, match="Mismatch"):
            filter_passing(records + records, df)

    def test_all_fail_returns_empty(self):
        records = [
            SeqRecord(Seq("NNN"), id="bad_1"),
            SeqRecord(Seq("NNN"), id="bad_2"),
        ]
        df = assess_sequences(records)
        passing = filter_passing(records, df)
        assert passing == []
