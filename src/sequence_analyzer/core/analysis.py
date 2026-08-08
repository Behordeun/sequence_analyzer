"""Sequence analysis: nucleotide composition, GC content, GC skew, base composition."""

import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


def analyze_sequences(sequences: list[SeqRecord], is_rna: bool = False) -> pd.DataFrame:
    """Compute nucleotide composition, GC content, and length for each sequence.

    Args:
        sequences: List of validated SeqRecords.
        is_rna: If True, count U instead of T.

    Returns:
        DataFrame with columns: Sequence, Length, GC_Content, and individual
        nucleotide counts (A, T/U, G, C).

    Raises:
        ValueError: If sequences list is empty.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for analysis.")

    rows: list[dict[str, object]] = []
    for record in sequences:
        seq_str = str(record.seq).upper()
        gc_content = gc_fraction(seq_str) * 100

        second_base = "U" if is_rna else "T"
        nucleotide_counts = {
            "A": seq_str.count("A"),
            second_base: seq_str.count(second_base),
            "G": seq_str.count("G"),
            "C": seq_str.count("C"),
        }

        rows.append(
            {
                "Sequence": record.id,
                "Length": len(seq_str),
                "GC_Content": gc_content,
                **nucleotide_counts,
            }
        )

    return pd.DataFrame(rows)


def compute_gc_skew(sequences: list[SeqRecord], window_size: int = 100) -> pd.DataFrame:
    """Compute GC skew in sliding windows across each sequence.

    GC skew = (G - C) / (G + C) for each window. Useful for identifying
    replication origins in bacterial genomes.

    Args:
        sequences: List of validated SeqRecords.
        window_size: Number of bases per window.

    Returns:
        DataFrame with columns: ID, Position, GC_Skew.
        Returns empty DataFrame if no sequence is long enough for the window.

    Raises:
        ValueError: If sequences list is empty or window_size < 1.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for GC skew analysis.")
    if window_size < 1:
        raise ValueError(f"window_size must be >= 1, got {window_size}.")

    rows: list[dict[str, object]] = []
    for record in sequences:
        seq = str(record.seq).upper()
        seq_len = len(seq)

        for i in range(seq_len - window_size + 1):
            window = seq[i : i + window_size]
            g = window.count("G")
            c = window.count("C")
            skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
            rows.append({"ID": record.id, "Position": i, "GC_Skew": skew})

    return pd.DataFrame(rows)


def compute_base_composition(sequences: list[SeqRecord]) -> pd.DataFrame:
    """Compute percentage base composition for each sequence.

    Args:
        sequences: List of validated SeqRecords.

    Returns:
        DataFrame with columns: ID, Length, A%, T%, G%, C%, GC%.

    Raises:
        ValueError: If sequences list is empty or any sequence has zero length.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for base composition.")

    rows: list[dict[str, object]] = []
    for record in sequences:
        seq = record.seq
        length = len(seq)

        if length == 0:
            raise ValueError(f"Sequence '{record.id}' has zero length.")

        rows.append(
            {
                "ID": record.id,
                "Length": length,
                "A%": seq.count("A") / length * 100,
                "T%": seq.count("T") / length * 100,
                "G%": seq.count("G") / length * 100,
                "C%": seq.count("C") / length * 100,
                "GC%": (seq.count("G") + seq.count("C")) / length * 100,
            }
        )

    return pd.DataFrame(rows)
