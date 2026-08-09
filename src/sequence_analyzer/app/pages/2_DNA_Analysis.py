"""DNA Sequence Analysis page."""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.pages.common_analysis import render_sequence_analysis_page
from sequence_analyzer.app.styles import apply_styles

st.set_page_config(layout="wide")
apply_styles()
hide_streamlit_chrome()

render_sequence_analysis_page(
    page_key="dna_analysis",
    title="🧬 DNA Sequence Analysis",
    seq_type="DNA",
    is_rna=False,
    upload_label="Upload DNA Sequence File",
    motif_placeholder="Enter DNA motif (e.g., ATG or TATA)",
    default_motif="ATG",
    download_prefix="dna_analysis",
)

render_footer()
