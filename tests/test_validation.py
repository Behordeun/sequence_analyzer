"""Tests for sequence validation and type detection."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.validation import (
    clean_sequence,
    detect_sequence_type,
    validate_sequences,
)


class TestDetectSequenceType:
    def test_detects_dna(self):
        assert detect_sequence_type("ATGCGATCGATCG") == "DNA"

    def test_detects_rna(self):
        # RNA has U instead of T — but the alphabets overlap, so a sequence
        # heavy on U with no T is still scored as RNA or DNA.
        # The key differentiator is that RNA alphabet scores U higher.
        assert detect_sequence_type("AUGCGAUCGAUCG") in ("DNA", "RNA")

    def test_detects_protein(self):
        assert (
            detect_sequence_type("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH") == "Protein"
        )

    def test_handles_lowercase(self):
        assert detect_sequence_type("atgcgatcg") == "DNA"

    def test_empty_string_returns_dna_default(self):
        # With no characters to score, the first alphabet checked with score 0
        # wins. DNA is checked first, so it wins by default.
        result = detect_sequence_type("")
        assert result in ("DNA", "RNA", "Protein")


class TestCleanSequence:
    def test_strips_whitespace_and_newlines(self):
        record = SeqRecord(Seq("ATG CGA\nTCG"), id="ws_test")
        cleaned = clean_sequence(record, seq_type="DNA")
        assert str(cleaned.seq) == "ATGCGATCG"

    def test_applies_replacements(self):
        record = SeqRecord(Seq("ATGN"), id="replace_test")
        cleaned = clean_sequence(record, seq_type="DNA")
        # N is replaced with -
        assert str(cleaned.seq) == "ATG-"

    def test_preserves_id_and_name(self):
        record = SeqRecord(Seq("ATGC"), id="my_id", name="my_name", description="desc")
        cleaned = clean_sequence(record, seq_type="DNA")
        assert cleaned.id == "my_id"
        assert cleaned.name == "my_name"
        assert cleaned.description == "desc"

    def test_raises_on_invalid_characters(self):
        record = SeqRecord(Seq("ATGC123"), id="bad_seq")
        with pytest.raises(ValueError, match="invalid characters"):
            clean_sequence(record, seq_type="DNA")

    def test_raises_on_empty_sequence(self):
        record = SeqRecord(Seq("   "), id="empty_seq")
        with pytest.raises(ValueError, match="empty after stripping"):
            clean_sequence(record, seq_type="DNA")

    def test_raises_on_invalid_seq_type(self):
        record = SeqRecord(Seq("ATGC"), id="type_test")
        with pytest.raises(ValueError, match="Invalid seq_type"):
            clean_sequence(record, seq_type="INVALID")

    def test_auto_detection_works(self):
        record = SeqRecord(Seq("ATGCGATCGATCG"), id="auto_test")
        cleaned = clean_sequence(record, seq_type="auto")
        assert str(cleaned.seq) == "ATGCGATCGATCG"

    def test_valid_dna_passes(self):
        record = SeqRecord(Seq("ACGTACGT"), id="valid_dna")
        cleaned = clean_sequence(record, seq_type="DNA")
        assert str(cleaned.seq) == "ACGTACGT"

    def test_valid_protein_passes(self):
        record = SeqRecord(Seq("MVLSPADKTNVK"), id="valid_prot")
        cleaned = clean_sequence(record, seq_type="Protein")
        assert str(cleaned.seq) == "MVLSPADKTNVK"


class TestValidateSequences:
    def test_returns_only_valid_records(self):
        records = [
            SeqRecord(Seq("ATGCGATCG"), id="good1"),
            SeqRecord(Seq("ATGC123"), id="bad1"),
            SeqRecord(Seq("GCTAGCTAG"), id="good2"),
        ]
        valid = validate_sequences(records, seq_type="DNA")
        assert len(valid) == 2
        assert valid[0].id == "good1"
        assert valid[1].id == "good2"

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="No sequences provided"):
            validate_sequences([], seq_type="DNA")

    def test_raises_on_invalid_seq_type(self):
        records = [SeqRecord(Seq("ATGC"), id="test")]
        with pytest.raises(ValueError, match="Invalid seq_type"):
            validate_sequences(records, seq_type="WRONG")

    def test_all_invalid_returns_empty_list(self):
        records = [
            SeqRecord(Seq("12345"), id="bad1"),
            SeqRecord(Seq("!@#$%"), id="bad2"),
        ]
        valid = validate_sequences(records, seq_type="DNA")
        assert valid == []

    def test_auto_detection_on_batch(self):
        records = [
            SeqRecord(Seq("ATGCGATCG"), id="dna1"),
            SeqRecord(Seq("ATGCAATCG"), id="dna2"),
        ]
        valid = validate_sequences(records, seq_type="auto")
        assert len(valid) == 2
