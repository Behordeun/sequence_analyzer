"""Typed result containers for sequence analysis operations."""

from dataclasses import dataclass, field


@dataclass
class AlignmentPair:
    """A single pairwise alignment comparison."""

    seq_a_id: str
    seq_b_id: str
    aligned_a: str
    aligned_b: str
    score: float
    identity: float


@dataclass
class AlignmentResult:
    """Result container for alignment operations."""

    pairs: list[AlignmentPair] = field(default_factory=list)
    avg_score: float = 0.0
    avg_identity: float = 0.0
    format_exports: dict[str, str] = field(default_factory=dict)


@dataclass
class PhylogeneticResult:
    """Result container for phylogenetic analysis."""

    tree: object = None  # Bio.Phylo tree object
    distance_matrix: object = None  # DistanceMatrix object
    alignment: object = None  # MultipleSeqAlignment object
    newick: str = ""


@dataclass
class MotifHit:
    """Result for a single sequence's motif scan."""

    sequence_id: str
    pattern: str
    hit_count: int
    positions: list[int] = field(default_factory=list)
