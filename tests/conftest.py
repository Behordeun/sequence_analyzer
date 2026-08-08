"""Shared test fixtures for the sequence_analyzer test suite."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def sample_dna_records() -> list[SeqRecord]:
    """Three short DNA sequences for testing."""
    return [
        SeqRecord(Seq("ATGCGATCGATCG"), id="seq1"),
        SeqRecord(Seq("ATGCGATCAATCG"), id="seq2"),
        SeqRecord(Seq("GCTAGCTAGCTAG"), id="seq3"),
    ]


@pytest.fixture
def sample_rna_records() -> list[SeqRecord]:
    """Three short RNA sequences for testing."""
    return [
        SeqRecord(Seq("AUGCGAUCGAUCG"), id="rna1"),
        SeqRecord(Seq("AUGCGAUCAAUCG"), id="rna2"),
        SeqRecord(Seq("GCUAGCUAGCUAG"), id="rna3"),
    ]


@pytest.fixture
def single_dna_record() -> SeqRecord:
    """A single DNA sequence for testing."""
    return SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="single_dna")
