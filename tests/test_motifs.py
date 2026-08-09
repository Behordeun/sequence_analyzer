"""Tests for motif scanning functions."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.motifs import scan_motifs


class TestScanMotifs:
    def test_finds_known_pattern(self):
        # ATG appears at positions 0 and 6
        records = [SeqRecord(Seq("ATGCCCATGCCC"), id="test1")]
        df = scan_motifs(records, "ATG")
        assert df.iloc[0]["Hits"] == 2
        assert df.iloc[0]["Positions"] == [0, 6]

    def test_no_matches_returns_zero_hits(self):
        records = [SeqRecord(Seq("CCCCCCCC"), id="no_match")]
        df = scan_motifs(records, "ATG")
        assert df.iloc[0]["Hits"] == 0
        assert df.iloc[0]["Positions"] == []

    def test_regex_pattern_works(self):
        # TA[TG]A should match TATA and TAGA
        records = [SeqRecord(Seq("TATACCCTAGACCC"), id="regex_test")]
        df = scan_motifs(records, "TA[TG]A")
        assert df.iloc[0]["Hits"] == 2

    def test_multiple_sequences(self, sample_dna_records):
        df = scan_motifs(sample_dna_records, "ATC")
        assert len(df) == 3
        # Each row has the correct columns
        assert list(df.columns) == ["ID", "Pattern", "Hits", "Positions"]

    def test_positions_within_valid_bounds(self):
        seq = "ATGATGATG"  # ATG at 0, 3, 6
        records = [SeqRecord(Seq(seq), id="bounds_test")]
        df = scan_motifs(records, "ATG")
        positions = df.iloc[0]["Positions"]
        pattern_len = 3
        for pos in positions:
            assert 0 <= pos <= len(seq) - pattern_len

    def test_hit_count_matches_positions_length(self, sample_dna_records):
        df = scan_motifs(sample_dna_records, "CG")
        for _, row in df.iterrows():
            assert row["Hits"] == len(row["Positions"])

    def test_raises_on_empty_sequences(self):
        with pytest.raises(ValueError, match="At least one sequence"):
            scan_motifs([], "ATG")

    def test_raises_on_empty_pattern(self):
        records = [SeqRecord(Seq("ATGC"), id="test")]
        with pytest.raises(ValueError, match="Pattern must not be empty"):
            scan_motifs(records, "")

    def test_raises_on_invalid_regex(self):
        records = [SeqRecord(Seq("ATGC"), id="test")]
        with pytest.raises(ValueError, match="Invalid regex pattern"):
            scan_motifs(records, "[unclosed")

    def test_overlapping_matches_not_counted(self):
        # "AAA" in "AAAA" → re.finditer finds non-overlapping: positions 0 only
        # (re.finditer doesn't overlap by default)
        records = [SeqRecord(Seq("AAAA"), id="overlap")]
        df = scan_motifs(records, "AAA")
        assert df.iloc[0]["Hits"] == 1
        assert df.iloc[0]["Positions"] == [0]
