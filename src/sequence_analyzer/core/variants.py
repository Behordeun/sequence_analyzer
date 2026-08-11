"""Reference-based variant calling from pairwise alignments.

Aligns each sample sequence against a reference and walks the alignment
to identify SNPs, insertions, and deletions. Produces a structured variant
table suitable for display or export.
"""

from __future__ import annotations

from collections import Counter
from enum import StrEnum

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.models.sequences import Variant, VariantResult


class VariantType(StrEnum):
    """Constrained variant type values."""

    SNP = "SNP"
    INSERTION = "insertion"
    DELETION = "deletion"


def _build_aligner(matrix: str = "NUC.4.4") -> PairwiseAligner:
    """Create a configured aligner instance for variant calling.

    Gap penalties are tuned so that single mismatches (score -4 in NUC.4.4)
    are always preferred over opening gap pairs. This ensures SNPs are
    reported as substitutions rather than deletion+insertion artifacts.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def _align_to_reference(
    reference: SeqRecord,
    sample: SeqRecord,
    aligner: PairwiseAligner,
) -> tuple[str, str]:
    """Align a single sample against the reference.

    Returns equal-length gapped strings built from the alignment coordinates.
    Biopython's Alignment.target/query attributes return ungapped strings in
    newer versions, so we reconstruct from coordinates for correctness.
    """
    best = aligner.align(reference.seq, sample.seq)[0]

    coords = best.coordinates
    ref_str = str(reference.seq)
    sample_str = str(sample.seq)

    aligned_ref_parts: list[str] = []
    aligned_sample_parts: list[str] = []

    for block_idx in range(len(coords[0]) - 1):
        ref_start, ref_end = coords[0][block_idx], coords[0][block_idx + 1]
        sample_start, sample_end = coords[1][block_idx], coords[1][block_idx + 1]

        ref_len = ref_end - ref_start
        sample_len = sample_end - sample_start

        if ref_len > 0 and sample_len > 0:
            aligned_ref_parts.append(ref_str[ref_start:ref_end])
            aligned_sample_parts.append(sample_str[sample_start:sample_end])
        elif ref_len == 0 and sample_len > 0:
            aligned_ref_parts.append("-" * sample_len)
            aligned_sample_parts.append(sample_str[sample_start:sample_end])
        elif ref_len > 0 and sample_len == 0:
            aligned_ref_parts.append(ref_str[ref_start:ref_end])
            aligned_sample_parts.append("-" * ref_len)

    return "".join(aligned_ref_parts), "".join(aligned_sample_parts)


def _handle_insertion(
    aligned_ref: str,
    aligned_sample: str,
    i: int,
    ref_pos: int,
    sample_id: str,
) -> tuple[Variant, int]:
    """Collect consecutive insertion bases and return the variant with updated index."""
    ins_bases: list[str] = []
    length = len(aligned_ref)
    while i < length and aligned_ref[i] == "-":
        ins_bases.append(aligned_sample[i])
        i += 1
    variant = Variant(
        position=max(ref_pos, 1),
        ref_base="-",
        sample_base="".join(ins_bases),
        variant_type=VariantType.INSERTION,
        sample_id=sample_id,
    )
    return variant, i


def _handle_deletion(
    aligned_ref: str,
    aligned_sample: str,
    i: int,
    ref_pos: int,
    sample_id: str,
) -> tuple[Variant, int, int]:
    """Collect consecutive deletion bases and return the variant with updated index and ref_pos."""
    del_bases: list[str] = []
    length = len(aligned_ref)
    start_pos = ref_pos + 1
    while i < length and aligned_sample[i] == "-" and aligned_ref[i] != "-":
        del_bases.append(aligned_ref[i])
        i += 1
        ref_pos += 1
    variant = Variant(
        position=start_pos,
        ref_base="".join(del_bases),
        sample_base="-",
        variant_type=VariantType.DELETION,
        sample_id=sample_id,
    )
    return variant, i, ref_pos


def _extract_variants(
    aligned_ref: str,
    aligned_sample: str,
    sample_id: str,
) -> list[Variant]:
    """Walk aligned strings to identify variants.

    Uses 1-based positioning relative to the reference. Groups consecutive
    gaps into single insertion/deletion events.
    """
    variants: list[Variant] = []
    ref_pos = 0
    i = 0
    length = len(aligned_ref)

    while i < length:
        ref_char = aligned_ref[i]
        sample_char = aligned_sample[i]

        if ref_char == "-":
            variant, i = _handle_insertion(aligned_ref, aligned_sample, i, ref_pos, sample_id)
            variants.append(variant)
        elif sample_char == "-":
            variant, i, ref_pos = _handle_deletion(
                aligned_ref, aligned_sample, i, ref_pos, sample_id
            )
            variants.append(variant)
        elif ref_char != sample_char:
            ref_pos += 1
            variants.append(
                Variant(
                    position=ref_pos,
                    ref_base=ref_char,
                    sample_base=sample_char,
                    variant_type=VariantType.SNP,
                    sample_id=sample_id,
                )
            )
            i += 1
        else:
            ref_pos += 1
            i += 1

    return variants


def _call_variants_for_sample(
    reference: SeqRecord,
    sample: SeqRecord,
    aligner: PairwiseAligner,
) -> list[Variant]:
    """Align and extract variants for a single sample."""
    if len(sample.seq) == 0:
        return []
    aligned_ref, aligned_sample = _align_to_reference(reference, sample, aligner)
    return _extract_variants(aligned_ref, aligned_sample, sample.id or "unnamed")


def call_variants(
    reference: SeqRecord,
    samples: list[SeqRecord],
    matrix: str = "NUC.4.4",
) -> VariantResult:
    """Identify variants in each sample relative to a reference sequence.

    Performs global pairwise alignment of each sample against the reference,
    then walks the alignment to detect SNPs, insertions, and deletions.

    Args:
        reference: The reference SeqRecord to compare against.
        samples: List of sample SeqRecords to call variants on.
        matrix: Substitution matrix for alignment (default NUC.4.4 for nucleotides).

    Returns:
        VariantResult containing all detected variants and a type summary.

    Raises:
        ValueError: If samples list is empty or reference sequence is empty.
    """
    if not samples:
        raise ValueError("At least one sample sequence is required for variant calling.")

    if len(reference.seq) == 0:
        raise ValueError("Reference sequence cannot be empty.")

    aligner = _build_aligner(matrix)
    all_variants: list[Variant] = []

    for sample in samples:
        all_variants.extend(_call_variants_for_sample(reference, sample, aligner))

    type_counts = Counter(v.variant_type for v in all_variants)
    summary = {
        "SNP": type_counts.get(VariantType.SNP, 0),
        "insertion": type_counts.get(VariantType.INSERTION, 0),
        "deletion": type_counts.get(VariantType.DELETION, 0),
        "total": len(all_variants),
    }

    return VariantResult(
        reference_id=reference.id or "unnamed",
        variants=all_variants,
        summary=summary,
    )


def summarize_variants(result: VariantResult) -> pd.DataFrame:
    """Convert a VariantResult into a positional variant table.

    Args:
        result: Output from call_variants().

    Returns:
        DataFrame with columns: Position, Ref, Alt, Type, Sample_ID.
        Empty DataFrame with correct columns if no variants found.
    """
    if not result.variants:
        return pd.DataFrame(columns=["Position", "Ref", "Alt", "Type", "Sample_ID"])

    rows = [
        {
            "Position": v.position,
            "Ref": v.ref_base,
            "Alt": v.sample_base,
            "Type": v.variant_type,
            "Sample_ID": v.sample_id,
        }
        for v in result.variants
    ]

    df = pd.DataFrame(rows)
    df = df.sort_values(["Position", "Sample_ID"]).reset_index(drop=True)
    return df


def compute_variant_density(
    result: VariantResult,
    reference_length: int,
    window_size: int = 50,
) -> pd.DataFrame:
    """Compute variant density across the reference in sliding windows.

    Useful for plotting variant hotspots along the reference sequence.

    Args:
        result: Output from call_variants().
        reference_length: Length of the reference sequence.
        window_size: Size of each window in bases.

    Returns:
        DataFrame with columns: Window_Start, Window_End, Variant_Count.

    Raises:
        ValueError: If reference_length < 1 or window_size < 1.
    """
    if reference_length < 1:
        raise ValueError("reference_length must be >= 1.")
    if window_size < 1:
        raise ValueError("window_size must be >= 1.")

    positions = [v.position for v in result.variants]
    rows: list[dict[str, int]] = []

    for start in range(1, reference_length + 1, window_size):
        end = min(start + window_size - 1, reference_length)
        count = sum(1 for p in positions if start <= p <= end)
        rows.append({"Window_Start": start, "Window_End": end, "Variant_Count": count})

    return pd.DataFrame(rows)
