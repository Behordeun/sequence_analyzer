"""Phylogenetic tree construction, bootstrapping, and entropy analysis.

Handles distance matrix calculation, tree building via Neighbor Joining or UPGMA,
bootstrap support estimation, and Shannon entropy per clade.
"""

import random
from collections import Counter
from io import StringIO

import numpy as np
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sequence_analyzer.models.sequences import PhylogeneticResult


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

    max_len = max(len(s.seq) for s in sequences)

    # Deduplicate IDs
    seen: dict[str, int] = {}
    padded: list[SeqRecord] = []
    for s in sequences:
        sid = s.id
        if sid in seen:
            seen[sid] += 1
            sid = f"{sid}_{seen[sid]}"
        else:
            seen[sid] = 1
        padded.append(SeqRecord(Seq(str(s.seq).ljust(max_len, "-")), id=sid))

    alignment = MultipleSeqAlignment(padded)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm) if method_lower == "nj" else constructor.upgma(dm)

    # Serialize to newick
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
) -> dict[frozenset[str], float]:
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
        ValueError: If replicates < 10 or method is invalid.
    """
    if replicates < 10:
        raise ValueError(f"Replicates must be >= 10, got {replicates}.")

    method_lower = method.lower()
    if method_lower not in ("nj", "upgma"):
        raise ValueError(f"Method must be 'nj' or 'upgma', got '{method}'.")

    aln_len = len(alignment[0])
    support: Counter[frozenset[str]] = Counter()

    for _ in range(replicates):
        # Resample columns with replacement
        idx = [random.randint(0, aln_len - 1) for _ in range(aln_len)]
        boot_records = [
            SeqRecord(Seq("".join(str(r.seq)[i] for i in idx)), id=r.id) for r in alignment
        ]
        boot_aln = MultipleSeqAlignment(boot_records)

        dm = DistanceCalculator("identity").get_distance(boot_aln)
        constructor = DistanceTreeConstructor()
        boot_tree = constructor.nj(dm) if method_lower == "nj" else constructor.upgma(dm)

        for clade in boot_tree.find_clades():
            tips = frozenset(x.name for x in clade.get_terminals())
            if len(tips) > 1:
                support[tips] += 1

    # Normalize to percentages
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
            freq = np.array(list(Counter(tips).values()), dtype=float)
            freq = freq / freq.sum()
            # Shannon entropy: -sum(p * log2(p))
            entropy = -float(np.sum(freq * np.log2(freq + 1e-12)))
            name = clade.name or f"internal_{unnamed_counter}"
            entropy_map[name] = entropy

    return entropy_map
