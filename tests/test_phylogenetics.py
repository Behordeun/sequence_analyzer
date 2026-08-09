"""Tests for phylogenetic tree construction, bootstrapping, and entropy."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.phylogenetics import (
    bootstrap_tree,
    build_tree,
    compute_entropy,
)


@pytest.fixture
def divergent_sequences() -> list[SeqRecord]:
    """Four sequences with enough divergence to build meaningful trees."""
    return [
        SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="seq_a"),
        SeqRecord(Seq("ATGCAATCAATCAATCAATCA"), id="seq_b"),
        SeqRecord(Seq("GCTAGCTAGCTAGCTAGCTAG"), id="seq_c"),
        SeqRecord(Seq("GCTAGCAAGCAAGCAAGCAAG"), id="seq_d"),
    ]


class TestBuildTree:
    def test_tree_has_correct_terminal_count(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        terminals = list(result.tree.get_terminals())
        assert len(terminals) == 4

    def test_nj_produces_tree(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        assert result.tree is not None
        assert result.newick != ""

    def test_upgma_produces_tree(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="upgma")
        assert result.tree is not None
        assert result.newick != ""

    def test_newick_string_non_empty(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        assert len(result.newick) > 10

    def test_distance_matrix_returned(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        assert result.distance_matrix is not None

    def test_alignment_returned(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        assert result.alignment is not None

    def test_deduplicates_ids(self):
        records = [
            SeqRecord(Seq("ATGCGATCG"), id="dup"),
            SeqRecord(Seq("GCTAGCTAG"), id="dup"),
            SeqRecord(Seq("ATGCAATCG"), id="dup"),
        ]
        result = build_tree(records, method="upgma")
        terminal_names = [t.name for t in result.tree.get_terminals()]
        # All names should be unique
        assert len(terminal_names) == len(set(terminal_names))

    def test_raises_on_fewer_than_2_sequences(self):
        records = [SeqRecord(Seq("ATGC"), id="lonely")]
        with pytest.raises(ValueError, match="At least 2 sequences"):
            build_tree(records)

    def test_raises_on_invalid_method(self, divergent_sequences):
        with pytest.raises(ValueError, match="Method must be"):
            build_tree(divergent_sequences, method="invalid")


class TestBootstrapTree:
    def test_support_values_in_valid_range(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        support = bootstrap_tree(result.alignment, method="nj", replicates=20)
        for _tips, value in support.items():
            assert 0.0 <= value <= 100.0

    def test_returns_dict_of_frozensets(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        support = bootstrap_tree(result.alignment, method="nj", replicates=10)
        for key in support:
            assert isinstance(key, frozenset)
            assert all(isinstance(name, str) for name in key)

    def test_raises_on_low_replicates(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        with pytest.raises(ValueError, match="Replicates must be >= 10"):
            bootstrap_tree(result.alignment, replicates=5)

    def test_raises_on_empty_alignment(self):
        from Bio.Align import MultipleSeqAlignment

        empty_alignment = MultipleSeqAlignment([])
        with pytest.raises(ValueError, match="Alignment is empty"):
            bootstrap_tree(empty_alignment, method="nj", replicates=10)

    def test_raises_on_invalid_method(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        with pytest.raises(ValueError, match="Method must be"):
            bootstrap_tree(result.alignment, method="bad")


class TestComputeEntropy:
    def test_terminal_clades_have_zero_entropy(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        entropy_map = compute_entropy(result.tree)
        # At least some entries should be 0 (terminals)
        zero_entries = [v for v in entropy_map.values() if v == 0.0]
        assert len(zero_entries) >= 4  # 4 terminal clades

    def test_returns_dict_of_strings_to_floats(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="upgma")
        entropy_map = compute_entropy(result.tree)
        for key, value in entropy_map.items():
            assert isinstance(key, str)
            assert isinstance(value, float)

    def test_entropy_non_negative(self, divergent_sequences):
        result = build_tree(divergent_sequences, method="nj")
        entropy_map = compute_entropy(result.tree)
        for value in entropy_map.values():
            assert value >= 0.0
