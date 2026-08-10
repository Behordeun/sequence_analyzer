"""Homepage for the Sequence Analyzer Streamlit application."""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles

st.set_page_config(page_title="Sequence Analyzer App", page_icon="🔬", layout="wide")
apply_styles()
hide_streamlit_chrome()

st.title("🔬 Sequence Analyzer Application")

st.markdown(
    """
## Guided Pipeline

Start the **guided workflow** to walk through your analysis step by step:
**Ingest → QC → Analyze → Align → Tree → Report**

Your data carries forward through each stage, and you'll get a downloadable
report at the end.
""",
)

st.markdown(
    "\U0001f680 Use the **sidebar navigation** to access the Guided Pipeline "
    "or any individual tool."
)

st.markdown("---")

st.markdown(
    """
## Individual Tools

Or use the **navigation menu on the left** to access individual tools directly:
1. Sequence Alignment
2. DNA / RNA Analysis
3. Phylogenetic Trees

---

## What You Can Do
- **Align Sequences:** Upload one or multiple sequences and choose between Pairwise or
  Multiple Sequence Alignment.
- **Analyze DNA or RNA:** See nucleotide composition and GC content, with charts.
- **Retrieve Sequences from GenBank:** Enter accession numbers and fetch data directly.
- **Build Phylogenetic Trees:** Generate evolutionary trees using UPGMA or Neighbor-Joining
  and download them in Newick format.

---

## How It Works
1. **Upload** sequences in .fasta, .txt, or .rtf format, or input GenBank accession numbers.
2. **Choose** an analysis type (DNA or RNA).
3. **Select** alignment and tree-building methods.
4. **Visualize and Download** your results as CSV, RTF, TXT, or Newick files.
5. **Save** your session for future reference.
6. **Export** aligned sequences in FASTA, CLUSTAL, NEXUS, and PHYLIP formats.
""",
    unsafe_allow_html=True,
)

render_footer()
