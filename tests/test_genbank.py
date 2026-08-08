"""Tests for GenBank fetching with mocked Entrez calls."""

from unittest.mock import MagicMock, patch

import pytest

from sequence_analyzer.core.genbank import fetch_sequence, fetch_sequences

MOCK_FASTA_RESPONSE = """>NM_001301717.1 Homo sapiens example
ATGCGATCGATCGATCGATCGATCGATCGATCG
"""


class TestFetchSequence:
    @patch("sequence_analyzer.core.genbank.Entrez.efetch")
    def test_successful_fetch_returns_seqrecord(self, mock_efetch):
        mock_handle = MagicMock()
        mock_handle.read.return_value = MOCK_FASTA_RESPONSE
        mock_handle.close = MagicMock()
        mock_efetch.return_value = mock_handle

        record = fetch_sequence("NM_001301717", email="test@example.com")
        assert record.id == "NM_001301717.1"
        assert "ATGCGATCGATCG" in str(record.seq)

    @patch("sequence_analyzer.core.genbank.Entrez.efetch")
    def test_network_failure_raises_connection_error(self, mock_efetch):
        mock_efetch.side_effect = Exception("Network unreachable")

        with pytest.raises(ConnectionError, match="Failed to reach NCBI"):
            fetch_sequence("NM_001301717", email="test@example.com")

    @patch("sequence_analyzer.core.genbank.Entrez.efetch")
    def test_not_found_raises_value_error(self, mock_efetch):
        mock_efetch.side_effect = Exception("HTTP Error 404: not found")

        with pytest.raises(ValueError, match="not found"):
            fetch_sequence("INVALID_ACC", email="test@example.com")

    def test_empty_accession_raises_value_error(self):
        with pytest.raises(ValueError, match="Accession must not be empty"):
            fetch_sequence("", email="test@example.com")

    def test_empty_email_raises_value_error(self):
        with pytest.raises(ValueError, match="Email is required"):
            fetch_sequence("NM_001301717", email="")


class TestFetchSequences:
    @patch("sequence_analyzer.core.genbank.fetch_sequence")
    def test_returns_successful_fetches(self, mock_fetch):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        mock_fetch.return_value = SeqRecord(Seq("ATGCGATCG"), id="NM_001")

        results = fetch_sequences(["NM_001", "NM_002"], email="test@example.com")
        assert len(results) == 2

    @patch("sequence_analyzer.core.genbank.fetch_sequence")
    def test_skips_failed_accessions(self, mock_fetch):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        def side_effect(acc, email):
            if acc == "BAD":
                raise ValueError("Not found")
            return SeqRecord(Seq("ATGC"), id=acc)

        mock_fetch.side_effect = side_effect

        results = fetch_sequences(["GOOD1", "BAD", "GOOD2"], email="test@example.com")
        assert len(results) == 2
        assert results[0].id == "GOOD1"
        assert results[1].id == "GOOD2"

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="At least one accession"):
            fetch_sequences([], email="test@example.com")

    def test_raises_on_empty_email(self):
        with pytest.raises(ValueError, match="Email is required"):
            fetch_sequences(["NM_001"], email="")
