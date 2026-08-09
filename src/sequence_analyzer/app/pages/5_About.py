"""About page."""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles

apply_styles()
hide_streamlit_chrome()

st.title("About This Application")

st.markdown(
    """
## Why Was This Application Created?
The **Sequence Analyzer Application** was developed to provide a user-friendly, web-based tool for:
- **Biologists and researchers** needing an accessible sequence analysis platform.
- **Students and educators** who want to learn about sequence alignment and genetic analysis.
- **Data analysts in genomics** who need quick insights without complex installations.

## How Does It Work?
1. **Upload a sequence file** (FASTA, TXT, RTF, PHYLIP, or NEXUS) or enter a GenBank accession number.
2. **Select the type of analysis** (DNA, RNA, or Phylogenetic).
3. **Choose an alignment method** (Pairwise or Multiple Sequence Alignment).
4. **View results and download aligned sequences** for further study.

## What Makes It Unique?
- No installation required (web UI) + installable Python library + CLI
- Multiple sequence file format support
- Interactive visualizations with Plotly
- GenBank integration for direct sequence fetching
- Export in CSV, TXT, RTF, FASTA, CLUSTAL, NEXUS, PHYLIP formats
- GC Skew analysis and Motif scanning
- Phylogenetic tree construction with bootstrapping

## Technologies Used
- **Python** with **Biopython**, **pandas**, **scikit-bio**, **scikit-learn**
- **Streamlit** for the web interface
- **Plotly** for visualizations
- **uv** for package management

---

## Developers
- **Muhammad Abiodun SULAIMAN** -
  [GitHub](https://github.com/Behordeun) |
  [ORCID](https://orcid.org/0000-0001-9161-2608)
- **Bolaji Fatai OYEYEMI** -
  [GitHub](https://github.com/bollergene) |
  [ORCID](https://orcid.org/0000-0001-5564-6165)

Version: **0.2.0**
[GitHub Repository](https://github.com/Behordeun/sequence_analyzer)
"""
)

render_footer()
