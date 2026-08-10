"""Phylogenetic tree construction, bootstrapping, and entropy analysis.

Handles distance matrix calculation, tree building via Neighbor Joining or UPGMA,
bootstrap support estimation, and Shannon entropy per clade.
"""

import math
import secrets
from collections import Counter
from io import StringIO

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.models.sequences import PhylogeneticResult

CladeKey = frozenset[str]


def _pad_and_deduplicate(sequences: list[SeqRecord]) -> MultipleSeqAlignment:
    """Pad sequences to equal length and deduplicate IDs."""
    max_len = max(len(s.seq) for s in sequences)

    seen: dict[str, int] = {}
    padded: list[SeqRecord] = []
    for s in sequences:
        sid = s.id or "unnamed"
        if sid in seen:
            seen[sid] += 1
            sid = f"{sid}_{seen[sid]}"
        else:
            seen[sid] = 1
        padded.append(SeqRecord(Seq(str(s.seq).ljust(max_len, "-")), id=sid))

    return MultipleSeqAlignment(padded)


def _build_tree_from_alignment(
    alignment: MultipleSeqAlignment,
    method_lower: str,
) -> Tree:
    """Build a tree from a MultipleSeqAlignment using the specified method."""
    dm = DistanceCalculator("identity").get_distance(alignment)
    constructor = DistanceTreeConstructor()
    return constructor.nj(dm) if method_lower == "nj" else constructor.upgma(dm)


def _resample_alignment(alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
    """Resample alignment columns with replacement for bootstrapping."""
    col_count = len(alignment[0])
    idx = [secrets.randbelow(col_count) for _ in range(col_count)]
    boot_records = [SeqRecord(Seq("".join(str(r.seq)[i] for i in idx)), id=r.id) for r in alignment]
    return MultipleSeqAlignment(boot_records)


def _collect_clade_tip_sets(tree: Tree) -> list[CladeKey]:
    """Collect frozensets of terminal names for each non-trivial clade."""
    clades: list[CladeKey] = []
    for clade in tree.find_clades():
        tips = frozenset(x.name for x in clade.get_terminals())
        if len(tips) > 1:
            clades.append(tips)
    return clades


def _shannon_entropy(counts: list[int]) -> float:
    """Compute Shannon entropy from a list of frequency counts."""
    total = sum(counts)
    if total == 0:
        return 0.0
    entropy = 0.0
    for c in counts:
        if c == 0:
            continue
        p = c / total
        entropy -= p * math.log2(p)
    return entropy


def build_tree(
    sequences: list[SeqRecord],
    method: str = "nj",
) -> PhylogeneticResult:
    """Construct a phylogenetic tree from sequences using NJ or UPGMA.

    Pads sequences to equal length, deduplicates IDs, computes a distance
    matrix via identity scoring, and builds the tree.

    Args:
        sequences: At least 2 SeqRecords.
        method: "nj" for Neighbor Joining, "upgma" for UPGMA.

    Returns:
        PhylogeneticResult with tree, distance_matrix, alignment, and newick string.

    Raises:
        ValueError: If fewer than 2 sequences or invalid method.
    """
    if len(sequences) < 2:
        raise ValueError("At least 2 sequences are required for tree construction.")

    method_lower = method.lower()
    if method_lower not in ("nj", "upgma"):
        raise ValueError(f"Method must be 'nj' or 'upgma', got '{method}'.")

    alignment = _pad_and_deduplicate(sequences)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm) if method_lower == "nj" else constructor.upgma(dm)

    newick_buf = StringIO()
    Phylo.write(tree, newick_buf, "newick")  # type: ignore[attr-defined]
    newick_str = newick_buf.getvalue().strip()

    return PhylogeneticResult(
        tree=tree,
        distance_matrix=dm,
        alignment=alignment,
        newick=newick_str,
    )


def bootstrap_tree(
    alignment: MultipleSeqAlignment,
    method: str = "nj",
    replicates: int = 50,
) -> dict[CladeKey, float]:
    """Run bootstrap analysis on an alignment.

    Resamples alignment columns with replacement, builds a tree from each
    replicate, and counts how often each clade appears. Returns support
    values as percentages.

    Args:
        alignment: A MultipleSeqAlignment (all sequences same length).
        method: "nj" or "upgma".
        replicates: Number of bootstrap replicates (minimum 10).

    Returns:
        Dict mapping frozenset of terminal names to support percentage (0-100).

    Raises:
        ValueError: If alignment is empty, replicates < 10, or method is invalid.
    """
    if not alignment:
        raise ValueError("Alignment is empty; cannot perform bootstrapping.")

    if replicates < 10:
        raise ValueError(f"Replicates must be >= 10, got {replicates}.")

    method_lower = method.lower()
    if method_lower not in ("nj", "upgma"):
        raise ValueError(f"Method must be 'nj' or 'upgma', got '{method}'.")

    support: Counter[CladeKey] = Counter()

    for _ in range(replicates):
        boot_aln = _resample_alignment(alignment)
        boot_tree = _build_tree_from_alignment(boot_aln, method_lower)

        for tips in _collect_clade_tip_sets(boot_tree):
            support[tips] += 1

    return {tips: (count / replicates * 100) for tips, count in support.items()}


def compute_entropy(tree: Tree) -> dict[str, float]:
    """Compute Shannon entropy for each clade in the tree.

    Terminal (leaf) clades have entropy 0. Internal clades have entropy based
    on the frequency distribution of their descendant terminal names.

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        Dict mapping clade name (or "unnamed_N") to entropy value.
    """
    entropy_map: dict[str, float] = {}

    for unnamed_counter, clade in enumerate(tree.find_clades()):
        if clade.is_terminal():
            name = clade.name or f"terminal_{unnamed_counter}"
            entropy_map[name] = 0.0
        else:
            tips = [x.name for x in clade.get_terminals()]
            counts = list(Counter(tips).values())
            entropy = _shannon_entropy(counts)
            name = clade.name or f"internal_{unnamed_counter}"
            entropy_map[name] = entropy

    return entropy_map
