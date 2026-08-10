"""Streamlit entry point for the Sequence Analyzer application.

Lives at the `app/` level so Streamlit's multipage discovery finds `pages/`
as a sibling directory and registers all sub-pages automatically.
"""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles

st.set_page_config(page_title="Sequence Analyzer App", page_icon="\U0001f52c", layout="wide")
apply_styles()
hide_streamlit_chrome()

st.title("\U0001f52c Sequence Analyzer Application")

st.markdown(
    """
## Guided Pipeline

Start the **guided workflow** to walk through your analysis step by step:
**Ingest \u2192 QC \u2192 Analyze \u2192 Align \u2192 Tree \u2192 Report**

Your data carries forward through each stage, and you'll get a downloadable
report at the end.
""",
)

st.markdown(
    "\U0001f680 Use the **sidebar navigation** to access individual analysis tools, "
    "or select a step below to start the guided pipeline."
)

col1, col2, col3 = st.columns(3)
with col1:
    st.page_link("pages/1_Sequence_Alignment.py", label="Sequence Alignment", icon="\U0001f9ec")
with col2:
    st.page_link("pages/2_DNA_Analysis.py", label="DNA Analysis", icon="\U0001f9ea")
with col3:
    st.page_link("pages/3_RNA_Analysis.py", label="RNA Analysis", icon="\U0001f52c")

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
