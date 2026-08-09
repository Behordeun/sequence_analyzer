"""Pipeline Stage 2: Quality Control — assess and filter sequences."""

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
from sequence_analyzer.core.qc import assess_sequences, filter_passing

st.set_page_config(layout="wide", page_title="Pipeline: QC")
apply_styles()
hide_streamlit_chrome()
render_config_sidebar()

state = get_pipeline_state()
state.stage = "QC"

render_progress_bar()
st.title("🔍 Stage 2: Quality Control")


render_help_panel("QC")

if not state.sequences:
    st.warning("No sequences loaded. Go back to the Ingest stage.")
    render_stage_nav()
    render_footer()
    st.stop()

# --- Run QC ---
config = state.config
qc_df = assess_sequences(
    state.sequences,
    ambiguity_threshold=config.qc_ambiguity_threshold,
    min_length=config.qc_min_length,
    max_gap_fraction=config.qc_max_gap_fraction,
)
state.qc_results = qc_df

# --- Summary ---
total = len(qc_df)
pass_count = (qc_df["Status"] == "Pass").sum()
warn_count = (qc_df["Status"] == "Warning").sum()
fail_count = (qc_df["Status"] == "Fail").sum()

col1, col2, col3, col4 = st.columns(4)
col1.metric("Total", total)
col2.metric("Pass", int(pass_count))
col3.metric("Warning", int(warn_count))
col4.metric("Fail", int(fail_count))

mean_length = qc_df["Length"].mean()
mean_gc = qc_df["GC_Percent"].mean()
st.markdown(f"**Mean Length:** {mean_length:.0f} bp | **Mean GC:** {mean_gc:.1f}%")

# --- Table ---
st.subheader("Per-Sequence Results")


# Color-code status
def _status_color(status: str) -> str:
    if status == "Pass":
        return "background-color: #c8e6c9"
    elif status == "Warning":
        return "background-color: #fff9c4"
    return "background-color: #ffcdd2"


styled_df = qc_df.style.applymap(_status_color, subset=["Status"])
st.dataframe(styled_df, use_container_width=True)

# --- Filter options ---
include_warnings = st.checkbox("Include sequences with warnings", value=True)

valid = filter_passing(state.sequences, qc_df, include_warnings=include_warnings)
state.valid_sequences = valid

st.markdown(f"**Sequences proceeding to analysis:** {len(valid)}")

if valid:
    mark_complete("QC")

# --- Download ---
st.download_button(
    "Download QC Report (CSV)",
    qc_df.to_csv(index=False),
    file_name="qc_report.csv",
    mime="text/csv",
)

st.markdown("---")
render_stage_nav()
render_footer()
