"""Tests for sequence analysis functions."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.analysis import (
    analyze_sequences,
    compute_base_composition,
    compute_gc_skew,
)


class TestAnalyzeSequences:
    def test_returns_correct_columns_dna(self, sample_dna_records):
        df = analyze_sequences(sample_dna_records, is_rna=False)
        assert list(df.columns) == ["Sequence", "Length", "GC_Content", "A", "T", "G", "C"]

    def test_returns_correct_columns_rna(self, sample_rna_records):
        df = analyze_sequences(sample_rna_records, is_rna=True)
        assert list(df.columns) == ["Sequence", "Length", "GC_Content", "A", "U", "G", "C"]

    def test_one_row_per_sequence(self, sample_dna_records):
        df = analyze_sequences(sample_dna_records, is_rna=False)
        assert len(df) == 3

    def test_known_gc_content(self):
        # GGCC = 4 bases, all G or C, so GC = 100%
        records = [SeqRecord(Seq("GGCC"), id="full_gc")]
        df = analyze_sequences(records, is_rna=False)
        assert df.iloc[0]["GC_Content"] == pytest.approx(100.0)

    def test_zero_gc_content(self):
        # AAAA = no G or C
        records = [SeqRecord(Seq("AAAA"), id="no_gc")]
        df = analyze_sequences(records, is_rna=False)
        assert df.iloc[0]["GC_Content"] == pytest.approx(0.0)

    def test_nucleotide_counts_sum_to_length(self, sample_dna_records):
        df = analyze_sequences(sample_dna_records, is_rna=False)
        for _, row in df.iterrows():
            total = row["A"] + row["T"] + row["G"] + row["C"]
            assert total == row["Length"]

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="At least one sequence"):
            analyze_sequences([], is_rna=False)


class TestComputeGCSkew:
    def test_row_count_formula(self):
        # Sequence of length 10, window 3 → (10 - 3 + 1) = 8 rows
        records = [SeqRecord(Seq("ATGCGATCGA"), id="ten")]
        df = compute_gc_skew(records, window_size=3)
        assert len(df) == 8

    def test_window_larger_than_sequence_returns_empty(self):
        records = [SeqRecord(Seq("ATG"), id="short")]
        df = compute_gc_skew(records, window_size=10)
        assert len(df) == 0

    def test_window_equals_sequence_length(self):
        records = [SeqRecord(Seq("ATGC"), id="exact")]
        df = compute_gc_skew(records, window_size=4)
        assert len(df) == 1

    def test_multiple_sequences(self, sample_dna_records):
        df = compute_gc_skew(sample_dna_records, window_size=5)
        # Each record has length 13, so (13-5+1) = 9 rows per sequence
        assert len(df) == 9 * 3

    def test_known_skew_value(self):
        # "GGGG" → G=4, C=0, skew = (4-0)/(4+0) = 1.0
        records = [SeqRecord(Seq("GGGG"), id="all_g")]
        df = compute_gc_skew(records, window_size=4)
        assert df.iloc[0]["GC_Skew"] == pytest.approx(1.0)

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="At least one sequence"):
            compute_gc_skew([], window_size=10)

    def test_raises_on_invalid_window_size(self, sample_dna_records):
        with pytest.raises(ValueError, match="window_size must be >= 1"):
            compute_gc_skew(sample_dna_records, window_size=0)


class TestComputeBaseComposition:
    def test_returns_correct_columns(self, sample_dna_records):
        df = compute_base_composition(sample_dna_records)
        assert list(df.columns) == ["ID", "Length", "A%", "T%", "G%", "C%", "GC%"]

    def test_percentages_sum_correctly(self):
        # ATGC → A=25%, T=25%, G=25%, C=25%
        records = [SeqRecord(Seq("ATGC"), id="even")]
        df = compute_base_composition(records)
        row = df.iloc[0]
        assert row["A%"] == pytest.approx(25.0)
        assert row["T%"] == pytest.approx(25.0)
        assert row["G%"] == pytest.approx(25.0)
        assert row["C%"] == pytest.approx(25.0)
        assert row["GC%"] == pytest.approx(50.0)

    def test_gc_percent_matches_gc_bases(self):
        # GGCCAA → GC = 4/6 ≈ 66.67%
        records = [SeqRecord(Seq("GGCCAA"), id="gc_test")]
        df = compute_base_composition(records)
        assert df.iloc[0]["GC%"] == pytest.approx(4 / 6 * 100)

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="At least one sequence"):
            compute_base_composition([])

    def test_raises_on_zero_length_sequence(self):
        records = [SeqRecord(Seq(""), id="empty")]
        with pytest.raises(ValueError, match="zero length"):
            compute_base_composition(records)
