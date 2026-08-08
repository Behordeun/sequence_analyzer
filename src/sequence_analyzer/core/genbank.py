"""GenBank sequence fetching via NCBI Entrez.

Wraps Biopython's Entrez API with proper error handling and configurable email.
"""

import logging

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def fetch_sequence(accession: str, email: str) -> SeqRecord:
    """Fetch a single sequence from NCBI GenBank by accession number.

    Args:
        accession: A valid NCBI accession (e.g., "NM_001301717").
        email: Email address for NCBI API compliance.

    Returns:
        A SeqRecord with the fetched sequence.

    Raises:
        ValueError: If accession is empty or the record is not found.
        ConnectionError: If NCBI is unreachable or returns an HTTP error.
    """
    if not accession or not accession.strip():
        raise ValueError("Accession must not be empty.")
    if not email or not email.strip():
        raise ValueError("Email is required for NCBI Entrez API usage.")

    accession = accession.strip()
    Entrez.email = email  # type: ignore[assignment]

    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        content = handle.read()
        handle.close()
    except Exception as e:
        error_msg = str(e).lower()
        if "404" in error_msg or "not found" in error_msg or "bad id" in error_msg:
            raise ValueError(f"Accession '{accession}' not found on NCBI.") from e
        raise ConnectionError(f"Failed to reach NCBI for accession '{accession}': {e}") from e

    if not content or content.startswith("Error") or "Nothing has been found" in content:
        raise ValueError(f"Accession '{accession}' not found on NCBI.")

    from io import StringIO

    record: SeqRecord = SeqIO.read(StringIO(content), "fasta")
    return record


def fetch_sequences(accessions: list[str], email: str) -> list[SeqRecord]:
    """Fetch multiple sequences from NCBI GenBank.

    Attempts to fetch each accession. Returns only successfully fetched records.
    Logs warnings for failed accessions rather than raising.

    Args:
        accessions: List of accession strings.
        email: Email address for NCBI API compliance.

    Returns:
        List of successfully fetched SeqRecords.

    Raises:
        ValueError: If accessions list is empty or email is missing.
    """
    if not accessions:
        raise ValueError("At least one accession is required.")
    if not email or not email.strip():
        raise ValueError("Email is required for NCBI Entrez API usage.")

    results: list[SeqRecord] = []
    for acc in accessions:
        try:
            record = fetch_sequence(acc, email)
            results.append(record)
            logger.info("Fetched %s successfully.", acc)
        except (ValueError, ConnectionError) as e:
            logger.warning("Failed to fetch %s: %s", acc, e)

    return results
