"""Typed result containers for sequence analysis operations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.BaseTree import Tree
    from Bio.Phylo.TreeConstruction import DistanceMatrix


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

    tree: Tree | None = None
    distance_matrix: DistanceMatrix | None = None
    alignment: MultipleSeqAlignment | None = None
    newick: str = ""


@dataclass
class MotifHit:
    """Result for a single sequence's motif scan."""

    sequence_id: str
    pattern: str
    hit_count: int
    positions: list[int] = field(default_factory=list)


@dataclass
class Variant:
    """A single variant detected between a sample and a reference sequence."""

    position: int
    ref_base: str
    sample_base: str
    variant_type: str  # "SNP", "insertion", "deletion"
    sample_id: str


@dataclass
class VariantResult:
    """Result container for reference-based variant calling."""

    reference_id: str
    variants: list[Variant] = field(default_factory=list)
    summary: dict[str, int] = field(default_factory=dict)
