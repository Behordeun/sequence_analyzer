"""File parsing and format conversion for sequence files.

Handles FASTA, PHYLIP, NEXUS, and plain text formats. Provides utilities
for both file-path-based and in-memory (bytes) parsing.
"""

from io import StringIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _detect_format(content: str) -> str:
    """Guess file format from content header."""
    stripped = content.strip()

    if stripped.startswith(">"):
        return "fasta"

    if stripped.upper().startswith("#NEXUS"):
        return "nexus"

    if stripped.upper().startswith("CLUSTAL"):
        return "clustal"

    # PHYLIP starts with a line containing two integers (ntaxa nalign)
    first_line = stripped.split("\n")[0].strip()
    parts = first_line.split()
    if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
        return "phylip"

    # Default to FASTA as the most permissive
    return "fasta"


def parse_sequence_file(content: str | bytes, format_hint: str = "auto") -> list[SeqRecord]:
    """Parse sequence content from various formats.

    Supports FASTA, PHYLIP, NEXUS. Auto-detects format when format_hint is "auto".

    Args:
        content: String or bytes containing the sequence data.
        format_hint: One of "auto", "fasta", "phylip", "nexus".

    Returns:
        List of parsed SeqRecords. May be empty if no sequences found.

    Raises:
        ValueError: If format_hint is unsupported or content cannot be decoded.
    """
    valid_formats = {"auto", "fasta", "phylip", "nexus", "fasta-pearson", "clustal"}
    if format_hint not in valid_formats:
        raise ValueError(
            f"Unsupported format '{format_hint}'. Must be one of: {sorted(valid_formats)}"
        )

    if isinstance(content, bytes):
        try:
            content = content.decode("utf-8")
        except UnicodeDecodeError as e:
            raise ValueError(f"Cannot decode file content as UTF-8: {e}") from e

    if not content.strip():
        return []

    fmt = _detect_format(content) if format_hint == "auto" else format_hint

    try:
        records = list(SeqIO.parse(StringIO(content), fmt))
    except Exception as e:
        # Try fasta-pearson as fallback (handles some edge cases)
        try:
            records = list(SeqIO.parse(StringIO(content), "fasta-pearson"))
        except Exception:
            raise ValueError(f"Failed to parse content as '{fmt}': {e}") from e

    return records


def parse_uploaded_file(file_bytes: bytes) -> list[SeqRecord]:
    """Parse an uploaded file's raw bytes into SeqRecords.

    Handles UTF-8 decoding and auto-detects format.

    Args:
        file_bytes: Raw file content as bytes.

    Returns:
        List of parsed SeqRecords.

    Raises:
        ValueError: If content cannot be decoded or parsed.
    """
    return parse_sequence_file(file_bytes, format_hint="auto")


def convert_to_fasta(content: str | bytes) -> str:
    """Convert raw file content to FASTA-formatted string.

    Parses the input in whatever format it is, then re-serializes as FASTA.

    Args:
        content: Raw string or bytes of sequence data.

    Returns:
        FASTA-formatted string of all sequences.
    """
    if isinstance(content, bytes):
        try:
            content = content.decode("utf-8")
        except UnicodeDecodeError:
            content = content.decode("latin-1")

    records = parse_sequence_file(content, format_hint="auto")
    if not records:
        return content  # Return as-is if parsing found nothing

    output = StringIO()
    SeqIO.write(records, output, "fasta")
    return output.getvalue()
