"""Tests for reference-based variant calling module."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.variants import (
    call_variants,
    compute_variant_density,
    summarize_variants,
)
from sequence_analyzer.models.sequences import VariantResult


@pytest.fixture
def reference() -> SeqRecord:
    """Reference sequence for variant calling tests."""
    return SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="reference")


@pytest.fixture
def identical_sample() -> SeqRecord:
    """Sample identical to the reference — produces zero variants."""
    return SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="identical")


@pytest.fixture
def snp_sample() -> SeqRecord:
    """Sample with a single SNP at position 5 (G->T)."""
    #                    01234567890
    return SeqRecord(Seq("ATGCTATCGATCGATCGATCG"), id="snp_sample")


@pytest.fixture
def multi_snp_sample() -> SeqRecord:
    """Sample with SNPs at positions 1, 10, 20 (A->C, A->T, C->A)."""
    #                         1         2
    #               1234567890123456789012
    return SeqRecord(Seq("CTGCGATCGTTCGATCGATAG"), id="multi_snp")


@pytest.fixture
def insertion_sample() -> SeqRecord:
    """Sample with an insertion (extra AAA between positions 5 and 6)."""
    return SeqRecord(Seq("ATGCGAAATCGATCGATCGATCG"), id="ins_sample")


@pytest.fixture
def deletion_sample() -> SeqRecord:
    """Sample with a deletion (missing bases 5-7: ATC removed)."""
    return SeqRecord(Seq("ATGCGGATCGATCGATCG"), id="del_sample")


class TestCallVariants:
    def test_identical_sequences_produce_no_variants(self, reference, identical_sample):
        result = call_variants(reference, [identical_sample])
        assert len(result.variants) == 0
        assert result.summary["total"] == 0
        assert result.reference_id == "reference"

    def test_single_snp_detected(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        snps = [v for v in result.variants if v.variant_type == "SNP"]
        assert len(snps) >= 1
        # Position 5: G->T
        snp = snps[0]
        assert snp.position == 5
        assert snp.ref_base == "G"
        assert snp.sample_base == "T"
        assert snp.sample_id == "snp_sample"

    def test_multiple_snps_detected(self, reference, multi_snp_sample):
        result = call_variants(reference, [multi_snp_sample])
        snps = [v for v in result.variants if v.variant_type == "SNP"]
        assert len(snps) >= 3
        assert result.summary["SNP"] >= 3

    def test_insertion_detected(self, reference, insertion_sample):
        result = call_variants(reference, [insertion_sample])
        insertions = [v for v in result.variants if v.variant_type == "insertion"]
        assert len(insertions) >= 1
        ins = insertions[0]
        assert ins.ref_base == "-"
        assert len(ins.sample_base) >= 1
        assert ins.sample_id == "ins_sample"

    def test_deletion_detected(self, reference, deletion_sample):
        result = call_variants(reference, [deletion_sample])
        deletions = [v for v in result.variants if v.variant_type == "deletion"]
        assert len(deletions) >= 1
        d = deletions[0]
        assert d.sample_base == "-"
        assert len(d.ref_base) >= 1
        assert d.sample_id == "del_sample"

    def test_multiple_samples_produces_variants_for_each(self, reference, snp_sample):
        sample2 = SeqRecord(Seq("ATGCGATCGATCGATCGATCC"), id="sample2")
        result = call_variants(reference, [snp_sample, sample2])
        sample_ids = {v.sample_id for v in result.variants}
        assert "snp_sample" in sample_ids
        assert "sample2" in sample_ids

    def test_returns_variant_result_type(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        assert isinstance(result, VariantResult)

    def test_summary_contains_all_types(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        assert "SNP" in result.summary
        assert "insertion" in result.summary
        assert "deletion" in result.summary
        assert "total" in result.summary

    def test_raises_on_empty_samples(self, reference):
        with pytest.raises(ValueError, match="At least one sample"):
            call_variants(reference, [])

    def test_raises_on_empty_reference(self):
        ref = SeqRecord(Seq(""), id="empty_ref")
        samples = [SeqRecord(Seq("ATGC"), id="s1")]
        with pytest.raises(ValueError, match="Reference sequence cannot be empty"):
            call_variants(ref, samples)

    def test_skips_empty_sample_without_crashing(self, reference):
        samples = [
            SeqRecord(Seq(""), id="empty"),
            SeqRecord(Seq("ATGCGATCGATCGATCGATCG"), id="good"),
        ]
        result = call_variants(reference, samples)
        # Should produce results only for the non-empty sample
        sample_ids = {v.sample_id for v in result.variants}
        assert "empty" not in sample_ids

    def test_variant_positions_are_positive(self, reference, multi_snp_sample):
        result = call_variants(reference, [multi_snp_sample])
        assert all(v.position >= 1 for v in result.variants)

    def test_variant_types_are_valid(self, reference, insertion_sample, deletion_sample):
        result = call_variants(reference, [insertion_sample, deletion_sample])
        valid_types = {"SNP", "insertion", "deletion"}
        assert all(v.variant_type in valid_types for v in result.variants)

    def test_reference_id_stored(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        assert result.reference_id == "reference"


class TestSummarizeVariants:
    def test_returns_dataframe_with_correct_columns(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        df = summarize_variants(result)
        expected_cols = {"Position", "Ref", "Alt", "Type", "Sample_ID"}
        assert set(df.columns) == expected_cols

    def test_empty_result_returns_empty_dataframe(self, reference, identical_sample):
        result = call_variants(reference, [identical_sample])
        df = summarize_variants(result)
        assert len(df) == 0
        assert set(df.columns) == {"Position", "Ref", "Alt", "Type", "Sample_ID"}

    def test_rows_sorted_by_position(self, reference, multi_snp_sample):
        result = call_variants(reference, [multi_snp_sample])
        df = summarize_variants(result)
        positions = df["Position"].tolist()
        assert positions == sorted(positions)

    def test_snp_row_values(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        df = summarize_variants(result)
        snp_rows = df[df["Type"] == "SNP"]
        assert len(snp_rows) >= 1
        row = snp_rows.iloc[0]
        assert row["Ref"] == "G"
        assert row["Alt"] == "T"
        assert row["Sample_ID"] == "snp_sample"

    def test_multiple_samples_in_table(self, reference, snp_sample):
        sample2 = SeqRecord(Seq("ATGCGATCGATCGATCGATCC"), id="s2")
        result = call_variants(reference, [snp_sample, sample2])
        df = summarize_variants(result)
        assert set(df["Sample_ID"].unique()).issuperset({"snp_sample", "s2"})


class TestComputeVariantDensity:
    def test_returns_correct_columns(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        df = compute_variant_density(result, reference_length=21, window_size=10)
        assert set(df.columns) == {"Window_Start", "Window_End", "Variant_Count"}

    def test_windows_cover_full_reference(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        df = compute_variant_density(result, reference_length=21, window_size=10)
        # Windows: 1-10, 11-20, 21-21
        assert df.iloc[0]["Window_Start"] == 1
        assert df.iloc[-1]["Window_End"] == 21

    def test_variant_counted_in_correct_window(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        df = compute_variant_density(result, reference_length=21, window_size=10)
        # SNP at position 5 should be in the first window (1-10)
        first_window = df.iloc[0]
        assert first_window["Variant_Count"] >= 1

    def test_no_variants_produces_zero_counts(self, reference, identical_sample):
        result = call_variants(reference, [identical_sample])
        df = compute_variant_density(result, reference_length=21, window_size=10)
        assert all(df["Variant_Count"] == 0)

    def test_raises_on_invalid_reference_length(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        with pytest.raises(ValueError, match="reference_length must be >= 1"):
            compute_variant_density(result, reference_length=0, window_size=10)

    def test_raises_on_invalid_window_size(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        with pytest.raises(ValueError, match="window_size must be >= 1"):
            compute_variant_density(result, reference_length=21, window_size=0)

    def test_window_size_larger_than_reference(self, reference, snp_sample):
        result = call_variants(reference, [snp_sample])
        df = compute_variant_density(result, reference_length=21, window_size=100)
        # Single window covering the entire reference
        assert len(df) == 1
        assert df.iloc[0]["Window_Start"] == 1
        assert df.iloc[0]["Window_End"] == 21
