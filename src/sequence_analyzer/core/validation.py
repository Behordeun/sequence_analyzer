"""Sequence cleaning, validation, and type detection.

Accepts raw strings or SeqRecords. Returns cleaned SeqRecords or raises
ValueError when sequences contain invalid characters after cleaning.
"""

import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

ALPHABETS: dict[str, str] = {
    "DNA": "ACGTURYKMSWBDHVN-",
    "RNA": "ACGUURYKMSWBDHVN-",
    "Protein": "ACDEFGHIKLMNPQRSTVWYBXZJUO-",
}

REPLACEMENTS: dict[str, dict[str, str]] = {
    "DNA": {"N": "-", "X": "-", "J": "-"},
    "RNA": {"N": "-", "X": "-", "J": "-"},
    "Protein": {"X": "-", "J": "-", "O": "-"},
}

VALID_SEQ_TYPES = frozenset({"DNA", "RNA", "Protein", "auto"})


def detect_sequence_type(sequence: str) -> str:
    """Auto-detect whether a sequence is DNA, RNA, or Protein.

    Scores each alphabet by counting how many characters in the sequence
    belong to that alphabet. Returns the type with the highest match count.
    Returns "DNA" by default for empty sequences.
    """
    seq_upper = sequence.upper()
    best_type = "DNA"
    best_score = 0

    for seq_type, alphabet in ALPHABETS.items():
        score = sum(1 for char in seq_upper if char in alphabet)
        if score > best_score:
            best_score = score
            best_type = seq_type

    return best_type


def clean_sequence(record: SeqRecord, seq_type: str = "auto") -> SeqRecord:
    """Strip whitespace, apply character replacements, validate against alphabet.

    Args:
        record: A BioPython SeqRecord to clean.
        seq_type: One of "DNA", "RNA", "Protein", or "auto" for detection.

    Returns:
        A new SeqRecord with the cleaned sequence.

    Raises:
        ValueError: If seq_type is invalid, sequence is empty after cleaning,
            or invalid characters remain after replacement.
    """
    if seq_type not in VALID_SEQ_TYPES:
        raise ValueError(
            f"Invalid seq_type '{seq_type}'. Must be one of: {sorted(VALID_SEQ_TYPES)}"
        )

    raw = str(record.seq).upper().replace(" ", "").replace("\n", "").replace("\r", "")

    if not raw:
        raise ValueError(f"Sequence '{record.id}' is empty after stripping whitespace.")

    detected = detect_sequence_type(raw) if seq_type == "auto" else seq_type
    type_replacements = REPLACEMENTS.get(detected, {})
    cleaned = "".join(type_replacements.get(char, char) for char in raw)

    alphabet = ALPHABETS[detected]
    invalid_chars = sorted({char for char in cleaned if char not in alphabet})

    if invalid_chars:
        raise ValueError(
            f"Sequence '{record.id}' contains invalid characters for {detected}: {invalid_chars}"
        )

    new_record = SeqRecord(
        Seq(cleaned),
        id=record.id,
        name=record.name,
        description=record.description,
    )
    return new_record


def validate_sequences(records: list[SeqRecord], seq_type: str = "auto") -> list[SeqRecord]:
    """Clean and validate a batch of sequences.

    Returns only valid records. Logs warnings for sequences that fail validation
    rather than raising, so partial batches can still proceed.

    Args:
        records: List of SeqRecords to validate.
        seq_type: One of "DNA", "RNA", "Protein", or "auto".

    Returns:
        List of successfully validated SeqRecords.

    Raises:
        ValueError: If records is empty or seq_type is invalid.
    """
    if seq_type not in VALID_SEQ_TYPES:
        raise ValueError(
            f"Invalid seq_type '{seq_type}'. Must be one of: {sorted(VALID_SEQ_TYPES)}"
        )

    if not records:
        raise ValueError("No sequences provided for validation.")

    valid: list[SeqRecord] = []
    for record in records:
        try:
            cleaned = clean_sequence(record, seq_type)
            valid.append(cleaned)
        except ValueError as e:
            logger.warning("Skipping sequence: %s", e)

    return valid
