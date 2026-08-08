"""Tests for file parsing and format conversion."""

import pytest

from sequence_analyzer.io.parsers import convert_to_fasta, parse_sequence_file, parse_uploaded_file

SAMPLE_FASTA = """>seq1
ATGCGATCGATCG
>seq2
GCTAGCTAGCTAG
"""

SAMPLE_NEXUS = """#NEXUS
begin data;
dimensions ntax=2 nchar=13;
format datatype=dna missing=? gap=-;
matrix
seq1 ATGCGATCGATCG
seq2 GCTAGCTAGCTAG
;
end;
"""

SAMPLE_PHYLIP = """2 13
seq1      ATGCGATCGATCG
seq2      GCTAGCTAGCTAG
"""


class TestParseSequenceFile:
    def test_parses_fasta_content(self):
        records = parse_sequence_file(SAMPLE_FASTA, format_hint="fasta")
        assert len(records) == 2
        assert records[0].id == "seq1"
        assert str(records[0].seq) == "ATGCGATCGATCG"

    def test_auto_detects_fasta(self):
        records = parse_sequence_file(SAMPLE_FASTA, format_hint="auto")
        assert len(records) == 2

    def test_parses_nexus_content(self):
        records = parse_sequence_file(SAMPLE_NEXUS, format_hint="nexus")
        assert len(records) == 2

    def test_auto_detects_nexus(self):
        records = parse_sequence_file(SAMPLE_NEXUS, format_hint="auto")
        assert len(records) == 2

    def test_parses_phylip_content(self):
        records = parse_sequence_file(SAMPLE_PHYLIP, format_hint="phylip")
        assert len(records) == 2

    def test_auto_detects_phylip(self):
        records = parse_sequence_file(SAMPLE_PHYLIP, format_hint="auto")
        assert len(records) == 2

    def test_handles_bytes_input(self):
        records = parse_sequence_file(SAMPLE_FASTA.encode("utf-8"), format_hint="auto")
        assert len(records) == 2

    def test_empty_content_returns_empty_list(self):
        records = parse_sequence_file("", format_hint="auto")
        assert records == []

    def test_whitespace_only_returns_empty_list(self):
        records = parse_sequence_file("   \n\n  ", format_hint="auto")
        assert records == []

    def test_raises_on_unsupported_format(self):
        with pytest.raises(ValueError, match="Unsupported format"):
            parse_sequence_file(SAMPLE_FASTA, format_hint="genbank")

    def test_raises_on_bad_encoding(self):
        bad_bytes = b"\xff\xfe" + b"\x00" * 100
        with pytest.raises(ValueError, match="Cannot decode"):
            parse_sequence_file(bad_bytes, format_hint="auto")


class TestParseUploadedFile:
    def test_parses_fasta_bytes(self):
        records = parse_uploaded_file(SAMPLE_FASTA.encode("utf-8"))
        assert len(records) == 2
        assert records[0].id == "seq1"

    def test_empty_bytes_returns_empty(self):
        records = parse_uploaded_file(b"")
        assert records == []


class TestConvertToFasta:
    def test_fasta_input_roundtrips(self):
        result = convert_to_fasta(SAMPLE_FASTA)
        assert ">seq1" in result
        assert ">seq2" in result
        assert "ATGCGATCGATCG" in result

    def test_bytes_input_works(self):
        result = convert_to_fasta(SAMPLE_FASTA.encode("utf-8"))
        assert ">seq1" in result

    def test_nexus_converts_to_fasta(self):
        result = convert_to_fasta(SAMPLE_NEXUS)
        assert result.startswith(">")
