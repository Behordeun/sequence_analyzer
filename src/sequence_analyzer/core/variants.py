"""Reference-based variant calling from pairwise alignments.

Aligns each sample sequence against a reference and walks the alignment
to identify SNPs, insertions, and deletions. Produces a structured variant
table suitable for display or export.
"""

from __future__ import annotations

from collections import Counter

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.models.sequences import Variant, VariantResult


def _align_to_reference(
    reference: SeqRecord,
    sample: SeqRecord,
    matrix: str = "NUC.4.4",
) -> tuple[str, str]:
    """Align a single sample against the reference using global pairwise alignment.

    Returns the aligned reference and sample strings (with gap characters)
    of equal length, suitable for position-by-position variant extraction.

    Gap penalties are set to strongly prefer mismatches over gap-open events,
    which ensures SNPs are reported as substitutions rather than
    deletion+insertion pairs.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(matrix)

    # Gap penalties tuned for variant calling: a single mismatch (score -4
    # in NUC.4.4) must always be preferred over opening a gap pair.
    # Open=-10 means the aligner only introduces gaps for genuine indels.
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(reference.seq, sample.seq)
    best = alignments[0]

    # Reconstruct gap-aligned strings from alignment coordinates.
    # coordinates is a 2xN array: row 0 = target (ref) positions,
    # row 1 = query (sample) positions. Gaps appear where one row
    # stays constant while the other advances.
    coords = best.coordinates
    ref_str = str(reference.seq)
    sample_str = str(sample.seq)

    aligned_ref_parts: list[str] = []
    aligned_sample_parts: list[str] = []

    for block_idx in range(len(coords[0]) - 1):
        ref_start, ref_end = coords[0][block_idx], coords[0][block_idx + 1]
        sample_start, sample_end = coords[1][block_idx], coords[1][block_idx + 1]

        ref_segment_len = ref_end - ref_start
        sample_segment_len = sample_end - sample_start

        if ref_segment_len > 0 and sample_segment_len > 0:
            # Aligned block: both advance
            aligned_ref_parts.append(ref_str[ref_start:ref_end])
            aligned_sample_parts.append(sample_str[sample_start:sample_end])
        elif ref_segment_len == 0 and sample_segment_len > 0:
            # Insertion in sample: reference has gaps
            aligned_ref_parts.append("-" * sample_segment_len)
            aligned_sample_parts.append(sample_str[sample_start:sample_end])
        elif ref_segment_len > 0 and sample_segment_len == 0:
            # Deletion in sample: sample has gaps
            aligned_ref_parts.append(ref_str[ref_start:ref_end])
            aligned_sample_parts.append("-" * ref_segment_len)

    return "".join(aligned_ref_parts), "".join(aligned_sample_parts)


def _extract_variants(
    aligned_ref: str,
    aligned_sample: str,
    sample_id: str,
) -> list[Variant]:
    """Walk aligned strings position-by-position to identify variants.

    Uses 1-based positioning relative to the reference (gaps in reference
    don't advance the reference position counter).

    Groups consecutive gaps into single insertion/deletion events rather
    than reporting each gap character independently.
    """
    variants: list[Variant] = []
    ref_pos = 0  # 1-based position in the ungapped reference

    i = 0
    length = len(aligned_ref)

    while i < length:
        ref_char = aligned_ref[i]
        sample_char = aligned_sample[i]

        if ref_char == "-":
            # Insertion in sample: reference has a gap
            # Collect consecutive insertion bases
            ins_bases: list[str] = []
            while i < length and aligned_ref[i] == "-":
                ins_bases.append(aligned_sample[i])
                i += 1
            # Position is the last reference base before the insertion
            variants.append(
                Variant(
                    position=max(ref_pos, 1),
                    ref_base="-",
                    sample_base="".join(ins_bases),
                    variant_type="insertion",
                    sample_id=sample_id,
                )
            )
        elif sample_char == "-":
            # Deletion in sample: sample has a gap
            # Collect consecutive deletion bases
            del_bases: list[str] = []
            while i < length and aligned_sample[i] == "-" and aligned_ref[i] != "-":
                del_bases.append(aligned_ref[i])
                i += 1
                ref_pos += 1
            # Position is where the deletion starts
            variants.append(
                Variant(
                    position=ref_pos - len(del_bases) + 1,
                    ref_base="".join(del_bases),
                    sample_base="-",
                    variant_type="deletion",
                    sample_id=sample_id,
                )
            )
        elif ref_char != sample_char:
            # SNP: both have a base but they differ
            ref_pos += 1
            variants.append(
                Variant(
                    position=ref_pos,
                    ref_base=ref_char,
                    sample_base=sample_char,
                    variant_type="SNP",
                    sample_id=sample_id,
                )
            )
            i += 1
        else:
            # Match: no variant
            ref_pos += 1
            i += 1

    return variants


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

    all_variants: list[Variant] = []

    for sample in samples:
        if len(sample.seq) == 0:
            continue

        aligned_ref, aligned_sample = _align_to_reference(reference, sample, matrix)
        sample_variants = _extract_variants(aligned_ref, aligned_sample, sample.id or "unnamed")
        all_variants.extend(sample_variants)

    type_counts = Counter(v.variant_type for v in all_variants)
    summary = {
        "SNP": type_counts.get("SNP", 0),
        "insertion": type_counts.get("insertion", 0),
        "deletion": type_counts.get("deletion", 0),
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
