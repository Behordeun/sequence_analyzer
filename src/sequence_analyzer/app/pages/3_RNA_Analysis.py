"""RNA Sequence Analysis page."""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.pages._common_analysis import render_sequence_analysis_page
from sequence_analyzer.app.styles import apply_styles

st.set_page_config(layout="wide")
apply_styles()
hide_streamlit_chrome()

render_sequence_analysis_page(
    page_key="rna_analysis",
    title="🧬 RNA Sequence Analysis",
    seq_type="RNA",
    is_rna=True,
    upload_label="Upload RNA Sequence File",
    motif_placeholder="Enter RNA motif (e.g., AUG or UUU)",
    default_motif="AUG",
    download_prefix="rna_analysis",
)

render_footer()
