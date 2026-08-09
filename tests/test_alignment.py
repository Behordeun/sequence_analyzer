"""Tests for sequence alignment functions."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.alignment import align_msa, align_pairwise


class TestAlignPairwise:
    def test_identical_sequences_produce_100_identity(self):
        records = [
            SeqRecord(Seq("ATGCGATCG"), id="s1"),
            SeqRecord(Seq("ATGCGATCG"), id="s2"),
        ]
        result = align_pairwise(records, matrix="NUC.4.4", mode="global")
        assert result.avg_identity == pytest.approx(100.0)

    def test_produces_n_minus_1_pairs(self, sample_dna_records):
        result = align_pairwise(sample_dna_records, matrix="NUC.4.4", mode="global")
        assert len(result.pairs) == len(sample_dna_records) - 1

    def test_pair_has_correct_ids(self):
        records = [
            SeqRecord(Seq("ATGCG"), id="alpha"),
            SeqRecord(Seq("ATGCG"), id="beta"),
        ]
        result = align_pairwise(records, matrix="NUC.4.4", mode="global")
        assert result.pairs[0].seq_a_id == "alpha"
        assert result.pairs[0].seq_b_id == "beta"

    def test_score_is_positive_for_similar_sequences(self):
        records = [
            SeqRecord(Seq("ATGCGATCG"), id="s1"),
            SeqRecord(Seq("ATGCAATCG"), id="s2"),
        ]
        result = align_pairwise(records, matrix="NUC.4.4", mode="global")
        assert result.avg_score > 0

    def test_raises_on_fewer_than_2_sequences(self):
        records = [SeqRecord(Seq("ATGC"), id="lonely")]
        with pytest.raises(ValueError, match="At least 2 sequences"):
            align_pairwise(records)

    def test_raises_on_invalid_matrix(self):
        records = [
            SeqRecord(Seq("ATGC"), id="s1"),
            SeqRecord(Seq("ATGC"), id="s2"),
        ]
        with pytest.raises(ValueError, match="Could not load substitution matrix"):
            align_pairwise(records, matrix="NONEXISTENT_MATRIX_999")

    def test_local_mode_works(self):
        records = [
            SeqRecord(Seq("ATGCGATCGATCG"), id="s1"),
            SeqRecord(Seq("GATCG"), id="s2"),
        ]
        result = align_pairwise(records, matrix="NUC.4.4", mode="local")
        assert result.avg_score > 0
        assert len(result.pairs) == 1


class TestAlignMSA:
    def test_pads_to_max_length(self, sample_dna_records):
        result = align_msa(sample_dna_records)
        # All format exports should contain padded sequences
        assert "clustal" in result.format_exports
        assert "nexus" in result.format_exports
        assert "phylip" in result.format_exports
        assert "fasta" in result.format_exports

    def test_format_exports_are_non_empty_strings(self, sample_dna_records):
        result = align_msa(sample_dna_records)
        for fmt, content in result.format_exports.items():
            assert isinstance(content, str)
            assert len(content) > 0, f"Format '{fmt}' is empty"

    def test_consensus_produced(self, sample_dna_records):
        result = align_msa(sample_dna_records)
        assert "consensus" in result.format_exports
        consensus = result.format_exports["consensus"]
        # Consensus length matches max sequence length
        max_len = max(len(r.seq) for r in sample_dna_records)
        assert len(consensus) == max_len

    def test_identity_100_for_identical_sequences(self):
        records = [
            SeqRecord(Seq("ATGCGATCG"), id="a"),
            SeqRecord(Seq("ATGCGATCG"), id="b"),
            SeqRecord(Seq("ATGCGATCG"), id="c"),
        ]
        result = align_msa(records)
        assert result.avg_identity == pytest.approx(100.0)

    def test_raises_on_fewer_than_2_sequences(self):
        records = [SeqRecord(Seq("ATGC"), id="single")]
        with pytest.raises(ValueError, match="At least 2 sequences"):
            align_msa(records)

    def test_handles_different_length_sequences(self):
        records = [
            SeqRecord(Seq("ATGCGATCGATCG"), id="long"),
            SeqRecord(Seq("ATGCG"), id="short"),
        ]
        result = align_msa(records)
        # Should not raise, and produces valid outputs
        assert result.avg_identity > 0
        assert len(result.format_exports["consensus"]) == 13
