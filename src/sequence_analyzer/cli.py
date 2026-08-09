"""Command-line interface for the Sequence Analyzer toolkit.

Provides subcommands: analyze, align, tree, serve.
"""

import argparse
import sys
from pathlib import Path

from Bio.SeqRecord import SeqRecord


def _load_valid_sequences(
    input_path: Path, seq_type: str = "auto", min_count: int = 1
) -> list[SeqRecord]:
    """Read, parse, and validate sequences from a file.

    Shared across CLI subcommands to avoid duplicating the file→records→validation flow.
    """
    from sequence_analyzer.core.validation import validate_sequences
    from sequence_analyzer.io.parsers import parse_sequence_file

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    content = input_path.read_text(encoding="utf-8")
    records = parse_sequence_file(content, format_hint="auto")

    if not records:
        raise ValueError(f"No sequences found in {input_path}")

    valid = validate_sequences(records, seq_type=seq_type)

    if len(valid) < min_count:
        raise ValueError(
            f"Need at least {min_count} valid sequence(s), got {len(valid)}."
        )

    return valid


def main() -> None:
    """Top-level CLI entry point. Dispatches to subcommands."""
    parser = argparse.ArgumentParser(
        prog="sequence-analyzer",
        description="Bioinformatics toolkit for DNA/RNA sequence analysis, alignment, and phylogenetics.",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- analyze ---
    p_analyze = subparsers.add_parser("analyze", help="Run nucleotide composition analysis")
    p_analyze.add_argument("--input", "-i", required=True, help="Path to input FASTA file")
    p_analyze.add_argument("--output", "-o", help="Path to output CSV file (default: stdout)")
    p_analyze.add_argument(
        "--type",
        "-t",
        choices=["DNA", "RNA", "auto"],
        default="auto",
        help="Sequence type (default: auto)",
    )

    # --- align ---
    p_align = subparsers.add_parser("align", help="Align sequences")
    p_align.add_argument("--input", "-i", required=True, help="Path to input FASTA file")
    p_align.add_argument("--output", "-o", help="Path to output file (default: stdout)")
    p_align.add_argument(
        "--method",
        "-m",
        choices=["pairwise", "msa"],
        default="msa",
        help="Alignment method (default: msa)",
    )
    p_align.add_argument(
        "--format",
        "-f",
        choices=["fasta", "clustal", "nexus", "phylip"],
        default="fasta",
        help="Output format (default: fasta)",
    )

    # --- tree ---
    p_tree = subparsers.add_parser("tree", help="Build a phylogenetic tree")
    p_tree.add_argument("--input", "-i", required=True, help="Path to input FASTA file")
    p_tree.add_argument("--output", "-o", help="Path to output Newick file (default: stdout)")
    p_tree.add_argument(
        "--method",
        "-m",
        choices=["nj", "upgma"],
        default="nj",
        help="Tree construction method (default: nj)",
    )

    # --- serve ---
    subparsers.add_parser("serve", help="Launch the Streamlit web application")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    try:
        if args.command == "analyze":
            _cmd_analyze(args)
        elif args.command == "align":
            _cmd_align(args)
        elif args.command == "tree":
            _cmd_tree(args)
        elif args.command == "serve":
            _cmd_serve()
    except (ValueError, FileNotFoundError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _cmd_analyze(args: argparse.Namespace) -> None:
    """Run sequence analysis from the command line."""
    from sequence_analyzer.core.analysis import analyze_sequences
    from sequence_analyzer.core.validation import detect_sequence_type

    valid = _load_valid_sequences(Path(args.input), seq_type=args.type)

    # Derive is_rna from all validated sequences via majority vote
    if args.type == "auto":
        rna_count = sum(
            1 for r in valid if detect_sequence_type(str(r.seq)) == "RNA"
        )
        is_rna = rna_count > len(valid) // 2
    else:
        is_rna = args.type == "RNA"

    df = analyze_sequences(valid, is_rna=is_rna)

    output = df.to_csv(index=False)
    if args.output:
        Path(args.output).write_text(output, encoding="utf-8")
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(output)


def _cmd_align(args: argparse.Namespace) -> None:
    """Run sequence alignment from the command line."""
    from sequence_analyzer.core.alignment import align_msa, align_pairwise

    valid = _load_valid_sequences(Path(args.input), min_count=2)

    if args.method == "pairwise":
        result = align_pairwise(valid)
        lines = []
        for pair in result.pairs:
            lines.append(f">{pair.seq_a_id} vs {pair.seq_b_id}")
            lines.append(f"Score: {pair.score:.2f}, Identity: {pair.identity:.1f}%")
            lines.append(pair.aligned_a)
            lines.append(pair.aligned_b)
            lines.append("")
        output = "\n".join(lines)
    else:
        result = align_msa(valid)
        fmt = args.format
        output = result.format_exports.get(fmt, result.format_exports.get("fasta", ""))

    if args.output:
        Path(args.output).write_text(output, encoding="utf-8")
        print(f"Alignment written to {args.output}", file=sys.stderr)
    else:
        print(output)


def _cmd_tree(args: argparse.Namespace) -> None:
    """Build phylogenetic tree from the command line."""
    from sequence_analyzer.core.phylogenetics import build_tree

    valid = _load_valid_sequences(Path(args.input), min_count=2)
    result = build_tree(valid, method=args.method)

    if args.output:
        Path(args.output).write_text(result.newick, encoding="utf-8")
        print(f"Tree written to {args.output}", file=sys.stderr)
    else:
        print(result.newick)


def _cmd_serve() -> None:
    """Launch the Streamlit web UI."""
    from sequence_analyzer.app.main import run

    run()
