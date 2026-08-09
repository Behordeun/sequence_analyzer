"""Sequence quality control: automated assessment and filtering.

Classifies sequences as Pass/Warning/Fail based on configurable thresholds
for ambiguity rate, length, gap fraction, and statistical outliers.
"""

import pandas as pd
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.validation import detect_sequence_type

# Characters treated as ambiguous per sequence type
_AMBIGUOUS_CHARS: dict[str, set[str]] = {
    "DNA": {"N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"},
    "RNA": {"N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"},
    "Protein": {"X", "B", "Z", "J"},
}


def assess_sequences(
    sequences: list[SeqRecord],
    ambiguity_threshold: float = 5.0,
    min_length: int = 10,
    max_gap_fraction: float = 50.0,
) -> pd.DataFrame:
    """Compute QC metrics and classify each sequence as Pass/Warning/Fail.

    Args:
        sequences: List of SeqRecords to assess.
        ambiguity_threshold: Max allowed ambiguity percentage before Fail (default 5%).
        min_length: Minimum sequence length before Fail (default 10).
        max_gap_fraction: Max allowed gap percentage before Fail (default 50%).

    Returns:
        DataFrame with columns: ID, Length, GC_Percent, Ambiguity_Rate,
        Gap_Fraction, Type, Status, Flags.

    Raises:
        ValueError: If sequences list is empty.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for QC assessment.")

    rows: list[dict[str, object]] = []

    for record in sequences:
        seq_str = str(record.seq).upper()
        length = len(seq_str)
        seq_type = detect_sequence_type(seq_str)

        # GC content
        g_count = seq_str.count("G")
        c_count = seq_str.count("C")
        gc_percent = (g_count + c_count) / length * 100 if length > 0 else 0.0

        # Ambiguity rate
        ambig_chars = _AMBIGUOUS_CHARS.get(seq_type, _AMBIGUOUS_CHARS["DNA"])
        ambig_count = sum(1 for ch in seq_str if ch in ambig_chars)
        ambiguity_rate = ambig_count / length * 100 if length > 0 else 0.0

        # Gap fraction
        gap_count = seq_str.count("-")
        gap_fraction = gap_count / length * 100 if length > 0 else 0.0

        rows.append(
            {
                "ID": record.id or "unnamed",
                "Length": length,
                "GC_Percent": round(gc_percent, 2),
                "Ambiguity_Rate": round(ambiguity_rate, 2),
                "Gap_Fraction": round(gap_fraction, 2),
                "Type": seq_type,
                "_raw_length": length,
                "_raw_gc": gc_percent,
            }
        )

    df = pd.DataFrame(rows)

    # Compute outlier boundaries for length and GC
    mean_length = df["_raw_length"].mean()
    std_length = df["_raw_length"].std(ddof=1) if len(df) > 1 else 0.0
    mean_gc = df["_raw_gc"].mean()
    std_gc = df["_raw_gc"].std(ddof=1) if len(df) > 1 else 0.0

    statuses: list[str] = []
    flags_list: list[str] = []

    for _, row in df.iterrows():
        flags: list[str] = []
        status = "Pass"

        # Fail conditions
        if row["Ambiguity_Rate"] > ambiguity_threshold:
            flags.append("High ambiguity")
            status = "Fail"

        if row["Length"] < min_length:
            flags.append("Too short")
            status = "Fail"

        if row["Gap_Fraction"] > max_gap_fraction:
            flags.append("Excessive gaps")
            status = "Fail"

        # Warning conditions (only if not already failed)
        if status != "Fail":
            if 2.0 <= row["Ambiguity_Rate"] <= ambiguity_threshold:
                flags.append("Borderline ambiguity")
                status = "Warning"

            if std_length > 0 and abs(row["_raw_length"] - mean_length) > 2 * std_length:
                flags.append("Length outlier")
                status = "Warning"

            if std_gc > 0 and abs(row["_raw_gc"] - mean_gc) > 2 * std_gc:
                flags.append("GC outlier")
                status = "Warning"

        statuses.append(status)
        flags_list.append("; ".join(flags) if flags else "")

    df["Status"] = statuses
    df["Flags"] = flags_list

    # Drop internal columns
    df = df.drop(columns=["_raw_length", "_raw_gc"])

    return df


def filter_passing(
    sequences: list[SeqRecord],
    qc_results: pd.DataFrame,
    include_warnings: bool = True,
) -> list[SeqRecord]:
    """Return sequences that passed QC.

    Args:
        sequences: Original list of SeqRecords (same order as qc_results).
        qc_results: DataFrame from assess_sequences().
        include_warnings: If True, include sequences with Warning status.

    Returns:
        Filtered list of SeqRecords.

    Raises:
        ValueError: If sequences and qc_results have different lengths.
    """
    if len(sequences) != len(qc_results):
        raise ValueError(f"Mismatch: {len(sequences)} sequences but {len(qc_results)} QC rows.")

    allowed_statuses = {"Pass"}
    if include_warnings:
        allowed_statuses.add("Warning")

    return [
        seq
        for seq, status in zip(sequences, qc_results["Status"], strict=True)
        if status in allowed_statuses
    ]
