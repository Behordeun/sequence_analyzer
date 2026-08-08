"""Motif scanning: find regex pattern occurrences across sequences."""

import re

import pandas as pd
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.models.sequences import MotifHit


def scan_motifs(sequences: list[SeqRecord], pattern: str) -> pd.DataFrame:
    """Scan each sequence for occurrences of a regex motif pattern.

    Args:
        sequences: List of validated SeqRecords.
        pattern: A regex pattern string (e.g., "ATG", "TA[TG]A").

    Returns:
        DataFrame with columns: ID, Pattern, Hits, Positions.
        Positions are 0-indexed start positions of each match.

    Raises:
        ValueError: If sequences is empty, pattern is empty, or pattern
            is an invalid regex.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for motif scanning.")
    if not pattern:
        raise ValueError("Pattern must not be empty.")

    try:
        compiled = re.compile(pattern)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern '{pattern}': {e}") from e

    hits: list[MotifHit] = []
    for record in sequences:
        seq_str = str(record.seq)
        matches = [m.start() for m in compiled.finditer(seq_str)]
        hits.append(
            MotifHit(
                sequence_id=record.id or '',
                pattern=pattern,
                hit_count=len(matches),
                positions=matches,
            )
        )

    rows = [
        {
            "ID": hit.sequence_id,
            "Pattern": hit.pattern,
            "Hits": hit.hit_count,
            "Positions": hit.positions,
        }
        for hit in hits
    ]
    return pd.DataFrame(rows)
