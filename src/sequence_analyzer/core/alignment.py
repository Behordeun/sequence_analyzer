"""Sequence alignment: pairwise and multiple sequence alignment (MSA).

Returns structured results. HTML rendering is the app layer's responsibility.
"""

from collections import Counter
from io import StringIO

from Bio.Align import MultipleSeqAlignment, PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.core.validation import detect_sequence_type
from sequence_analyzer.models.sequences import AlignmentPair, AlignmentResult


def _pad_sequences(sequences: list[SeqRecord], mol_type: str) -> list[SeqRecord]:
    """Pad sequences to equal length with gap characters."""
    max_len = max(len(r.seq) for r in sequences)
    return [
        SeqRecord(
            Seq(str(r.seq).ljust(max_len, "-")),
            id=r.id,
            name=r.id or "",
            description="",
            annotations={"molecule_type": mol_type},
        )
        for r in sequences
    ]


def _write_alignment_formats(alignment: MultipleSeqAlignment) -> dict[str, str]:
    """Serialize an alignment to multiple output formats."""
    from Bio import AlignIO

    format_exports: dict[str, str] = {}
    for fmt in ("clustal", "nexus", "phylip", "fasta"):
        buf = StringIO()
        AlignIO.write(alignment, buf, fmt)
        format_exports[fmt] = buf.getvalue()
    return format_exports


def align_pairwise(
    sequences: list[SeqRecord],
    matrix: str = "NUC.4.4",
    mode: str = "global",
) -> AlignmentResult:
    """Perform sequential pairwise alignment on a list of sequences.

    Aligns each consecutive pair (seq[0] vs seq[1], seq[1] vs seq[2], etc.)
    and returns an AlignmentResult with per-pair details plus averages.

    Args:
        sequences: At least 2 SeqRecords.
        matrix: Substitution matrix name (e.g., "NUC.4.4", "BLOSUM62").
        mode: Alignment mode — "global" or "local".

    Returns:
        AlignmentResult with pairs, avg_score, and avg_identity.

    Raises:
        ValueError: If fewer than 2 sequences, unsupported mode, or matrix cannot be loaded.
    """
    if len(sequences) < 2:
        raise ValueError("At least 2 sequences are required for pairwise alignment.")

    if mode not in {"global", "local"}:
        raise ValueError(f"Unsupported alignment mode: '{mode}'. Must be 'global' or 'local'.")

    aligner = PairwiseAligner()
    aligner.mode = mode

    try:
        aligner.substitution_matrix = substitution_matrices.load(matrix)
    except Exception as e:
        raise ValueError(f"Could not load substitution matrix '{matrix}': {e}") from e

    pairs: list[AlignmentPair] = []
    total_score = 0.0
    total_identity = 0.0

    for i in range(len(sequences) - 1):
        a, b = sequences[i], sequences[i + 1]
        aln = aligner.align(a.seq, b.seq)[0]
        aligned_a = str(aln.target)
        aligned_b = str(aln.query)

        length = len(aligned_a)
        identity = sum(x == y for x, y in zip(aligned_a, aligned_b, strict=False)) / length * 100

        pairs.append(
            AlignmentPair(
                seq_a_id=a.id or "",
                seq_b_id=b.id or "",
                aligned_a=aligned_a,
                aligned_b=aligned_b,
                score=float(aln.score),
                identity=identity,
            )
        )

        total_score += aln.score
        total_identity += identity

    count = len(pairs)
    return AlignmentResult(
        pairs=pairs,
        avg_score=total_score / count if count else 0.0,
        avg_identity=total_identity / count if count else 0.0,
    )


def align_msa(
    sequences: list[SeqRecord],
    seq_type: str = "auto",
) -> AlignmentResult:
    """Perform multiple sequence alignment by padding sequences to max length.

    Pads shorter sequences with gap characters to produce a simple MSA.
    Generates format exports (CLUSTAL, NEXUS, PHYLIP, FASTA) as strings.

    Args:
        sequences: At least 2 SeqRecords.
        seq_type: Sequence type for format annotation ("DNA", "RNA", "Protein", "auto").

    Returns:
        AlignmentResult with format_exports dict and identity percentage.

    Raises:
        ValueError: If fewer than 2 sequences provided.
    """
    if len(sequences) < 2:
        raise ValueError("At least 2 sequences are required for MSA.")

    mol_type = detect_sequence_type(str(sequences[0].seq)) if seq_type == "auto" else seq_type
    padded = _pad_sequences(sequences, mol_type)

    alignment = MultipleSeqAlignment(padded)
    alignment.annotations["molecule_type"] = mol_type

    format_exports = _write_alignment_formats(alignment)

    # Column-wise identity against first sequence as reference
    seq_matrix = [str(r.seq) for r in padded]
    match_count = 0
    total = 0

    for i in range(len(seq_matrix)):
        for j, base in enumerate(seq_matrix[i]):
            ref = seq_matrix[0][j]
            if base == ref:
                match_count += 1
            total += 1

    identity = (match_count / total * 100) if total else 0.0

    # Compute consensus (exclude gaps so actual bases win over gap characters)
    consensus = "".join(
        Counter(b for b in col if b != "-").most_common(1)[0][0]
        if any(b != "-" for b in col)
        else "-"
        for col in zip(*seq_matrix, strict=False)
    )
    format_exports["consensus"] = consensus

    return AlignmentResult(
        pairs=[],
        avg_score=0.0,
        avg_identity=identity,
        format_exports=format_exports,
    )
