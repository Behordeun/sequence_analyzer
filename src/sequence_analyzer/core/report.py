"""Report generation: methods section, HTML report, and PDF export.

Produces publication-ready analysis reports from pipeline results.
"""

import datetime
import html as html_lib
from io import BytesIO

import pandas as pd

from sequence_analyzer import __version__


def generate_methods_section(
    seq_count: int,
    qc_pass_count: int,
    seq_type: str = "DNA",
    alignment_method: str = "msa",
    tree_method: str = "nj",
    bootstrap_replicates: int = 50,
    gc_window_size: int = 100,
) -> str:
    """Generate a methods paragraph in third-person past tense.

    Suitable for inclusion in a scientific publication's methods section.

    Args:
        seq_count: Total number of input sequences.
        qc_pass_count: Number of sequences passing QC.
        seq_type: Detected or specified sequence type.
        alignment_method: "msa" or "pairwise".
        tree_method: "nj" or "upgma".
        bootstrap_replicates: Number of bootstrap replicates performed.
        gc_window_size: Window size used for GC skew computation.

    Returns:
        A string containing the methods paragraph.
    """
    alignment_name = (
        "multiple sequence alignment (MSA) via gap-padding"
        if alignment_method == "msa"
        else "sequential pairwise alignment"
    )
    tree_name = "Neighbor-Joining (NJ)" if tree_method == "nj" else "UPGMA"

    return (
        f"A total of {seq_count} {seq_type} sequences were submitted for analysis. "
        f"Quality control filtering retained {qc_pass_count} sequences based on "
        f"ambiguity rate, minimum length, and gap fraction thresholds. "
        f"Nucleotide composition and GC content were computed for each sequence. "
        f"GC skew was calculated using a sliding window of {gc_window_size} bases. "
        f"Sequences were aligned using {alignment_name}. "
        f"A phylogenetic tree was constructed using the {tree_name} method "
        f"with {bootstrap_replicates} bootstrap replicates for branch support estimation. "
        f"All analyses were performed using Sequence Analyzer v{__version__}."
    )


def generate_html_report(
    title: str,
    methods_section: str,
    qc_table: pd.DataFrame | None = None,
    analysis_table: pd.DataFrame | None = None,
    figures_html: dict[str, str] | None = None,
    alignment_summary: str = "",
    tree_newick: str = "",
    notes: str = "",
) -> str:
    """Assemble a complete HTML report with methods, figures, and tables.

    Args:
        title: Report title (user-provided, will be HTML-escaped).
        methods_section: Generated methods paragraph (will be HTML-escaped).
        qc_table: QC assessment DataFrame (optional).
        analysis_table: Sequence analysis DataFrame (optional).
        figures_html: Dict of figure_name → Plotly HTML string (optional).
            TRUST ASSUMPTION: values must be from Plotly's to_html() output.
            Do not pass user-provided HTML without sanitization.
        alignment_summary: Text summary of alignment results (will be HTML-escaped).
        tree_newick: Newick format tree string (will be HTML-escaped).
        notes: User-provided notes (will be HTML-escaped).

    Returns:
        Complete HTML string suitable for saving as a standalone file.
    """
    date_str = datetime.date.today().isoformat()
    sections: list[str] = []

    # Header
    sections.append(f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{_escape_html(title)}</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 900px; margin: 0 auto; padding: 20px; color: #333; }}
h1 {{ color: #122C39; border-bottom: 2px solid #122C39; padding-bottom: 8px; }}
h2 {{ color: #2E7D32; margin-top: 30px; }}
table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; font-size: 14px; }}
th {{ background-color: #f4f4f9; }}
tr:nth-child(even) {{ background-color: #fafafa; }}
.metadata {{ color: #666; font-size: 12px; margin-top: 40px; border-top: 1px solid #ddd; padding-top: 10px; }}
pre {{ background: #f5f5f5; padding: 10px; border-radius: 4px; overflow-x: auto; }}
</style>
</head>
<body>
<h1>{_escape_html(title)}</h1>""")

    # Methods
    sections.append(f"<h2>Methods</h2>\n<p>{_escape_html(methods_section)}</p>")

    # QC Table
    if qc_table is not None and not qc_table.empty:
        sections.append("<h2>Quality Control</h2>")
        sections.append(qc_table.to_html(index=False, classes="qc-table"))

    # Analysis Table
    if analysis_table is not None and not analysis_table.empty:
        sections.append("<h2>Sequence Analysis</h2>")
        sections.append(analysis_table.to_html(index=False, classes="analysis-table"))

    # Figures — TRUST ASSUMPTION: figures_html values are Plotly-generated HTML
    # from fig.to_html(). Callers must not pass user-provided HTML here.
    # If reusing this function with untrusted input, sanitize figures_html
    # before calling (e.g., strip <script> tags or use a whitelist).
    if figures_html:
        for fig_name, fig_content in figures_html.items():
            sections.append(f"<h2>{_escape_html(fig_name)}</h2>")
            sections.append(f"<div class='figure'>{fig_content}</div>")

    # Alignment
    if alignment_summary:
        sections.append("<h2>Alignment</h2>")
        sections.append(f"<p>{_escape_html(alignment_summary)}</p>")

    # Tree
    if tree_newick:
        sections.append("<h2>Phylogenetic Tree</h2>")
        sections.append(f"<pre>{_escape_html(tree_newick)}</pre>")

    # Notes
    if notes:
        sections.append("<h2>Notes</h2>")
        sections.append(f"<p>{_escape_html(notes)}</p>")

    # Footer
    sections.extend([
        f'<div class="metadata">Generated by Sequence Analyzer v{__version__} on {date_str}</div>',
        "</body>\n</html>",
    ])

    return "\n".join(sections)


def generate_pdf_report(html_content: str) -> bytes:
    """Convert HTML report to PDF bytes using reportlab.

    Produces a simple PDF with the report content. For complex layouts,
    the HTML version is recommended.

    Args:
        html_content: The HTML report string (used for text extraction).

    Returns:
        PDF file content as bytes.
    """
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer

    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, topMargin=0.75 * inch)
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "ReportTitle",
        parent=styles["Heading1"],
        fontSize=18,
        spaceAfter=12,
    )
    heading_style = ParagraphStyle(
        "ReportHeading",
        parent=styles["Heading2"],
        fontSize=14,
        spaceAfter=8,
    )
    body_style = ParagraphStyle(
        "ReportBody",
        parent=styles["Normal"],
        fontSize=10,
        spaceAfter=6,
        leading=14,
    )

    story: list = []

    # Extract text content from HTML (simple approach for PDF)
    lines = _extract_text_from_html(html_content)

    for line in lines:
        if line.startswith("# "):
            story.append(Paragraph(line[2:], title_style))
        elif line.startswith("## "):
            story.extend([Spacer(1, 12), Paragraph(line[3:], heading_style)])
        elif line.strip():
            story.append(Paragraph(line, body_style))

    if not story:
        story.append(Paragraph("No content", body_style))

    doc.build(story)
    return buffer.getvalue()


def _escape_html(text: str) -> str:
    """Escape HTML special characters."""
    return (
        text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
    )


def _extract_text_from_html(html_content: str) -> list[str]:
    """Extract readable text lines from HTML for PDF generation."""
    import re

    text = html_content

    # Convert headings to markdown-style markers
    text = re.sub(r"<h1[^>]*>(.*?)</h1>", r"# \1", text)
    text = re.sub(r"<h2[^>]*>(.*?)</h2>", r"## \1", text)

    # Strip remaining tags
    text = re.sub(r"<[^>]+>", "", text)

    # Decode HTML entities using standard library
    text = html_lib.unescape(text)

    return [line.strip() for line in text.split("\n") if line.strip()]
