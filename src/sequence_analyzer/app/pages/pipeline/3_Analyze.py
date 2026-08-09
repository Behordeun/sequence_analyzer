"""Pipeline Stage 3: Analyze — nucleotide composition, GC skew, motifs."""

import plotly.express as px
import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.pages.pipeline.help_content import render_help_panel
from sequence_analyzer.app.pages.pipeline.navigation import (
    render_config_sidebar,
    render_progress_bar,
    render_stage_nav,
)
from sequence_analyzer.app.pages.pipeline.state import get_pipeline_state, mark_complete
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.analysis import (
    analyze_sequences,
    compute_base_composition,
    compute_gc_skew,
)
from sequence_analyzer.core.motifs import scan_motifs
from sequence_analyzer.core.validation import detect_sequence_type

GC_SKEW_LABEL = "GC Skew"

st.set_page_config(layout="wide", page_title="Pipeline: Analyze")
apply_styles()
hide_streamlit_chrome()
render_config_sidebar()

state = get_pipeline_state()
state.stage = "Analyze"

render_progress_bar()
st.title("📊 Stage 3: Sequence Analysis")


render_help_panel("Analyze")

sequences = state.valid_sequences or state.sequences
if not sequences:
    st.warning("No sequences available. Complete the Ingest or QC stage first.")
    render_stage_nav()
    render_footer()
    st.stop()

config = state.config

# Derive is_rna
if config.seq_type == "auto":
    rna_count = sum(1 for r in sequences if detect_sequence_type(str(r.seq)) == "RNA")
    is_rna = rna_count > len(sequences) // 2
else:
    is_rna = config.seq_type == "RNA"

# --- Run analysis ---
df = analyze_sequences(sequences, is_rna=is_rna)
state.analysis_results = df

st.subheader("Nucleotide Composition")
st.dataframe(df, use_container_width=True)

# --- Visualizations ---
col1, col2 = st.columns(2)

with col1:
    fig_gc = px.histogram(df, x="GC_Content", nbins=10, title="GC Content Distribution")
    st.plotly_chart(fig_gc, use_container_width=True)
    state.figures["GC Content Distribution"] = fig_gc.to_html(
        full_html=False, include_plotlyjs="cdn"
    )

with col2:
    fig_len = px.histogram(df, x="Length", nbins=10, title="Sequence Length Distribution")
    st.plotly_chart(fig_len, use_container_width=True)
    state.figures["Sequence Length Distribution"] = fig_len.to_html(
        full_html=False, include_plotlyjs="cdn"
    )

# --- GC Skew ---
st.subheader(GC_SKEW_LABEL)
df_skew = compute_gc_skew(sequences, window_size=config.gc_window_size)
if not df_skew.empty:
    fig_skew = px.line(df_skew, x="Position", y="GC_Skew", color="ID", title=GC_SKEW_LABEL)
    st.plotly_chart(fig_skew, use_container_width=True)
    state.figures[GC_SKEW_LABEL] = fig_skew.to_html(full_html=False, include_plotlyjs="cdn")
else:
    st.info("Sequences are shorter than the GC skew window size.")

# --- Motif scan ---
if config.motif_pattern:
    st.subheader("Motif Scan")
    try:
        df_motif = scan_motifs(sequences, config.motif_pattern)
        st.dataframe(df_motif, use_container_width=True)
    except ValueError as e:
        st.warning(f"Motif scan skipped: {e}")

# --- Base composition ---
st.subheader("Base Composition")
df_base = compute_base_composition(sequences)
st.dataframe(df_base, use_container_width=True)

mark_complete("Analyze")

st.markdown("---")
render_stage_nav()
render_footer()
