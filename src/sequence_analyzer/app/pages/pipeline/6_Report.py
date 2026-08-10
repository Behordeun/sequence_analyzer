"""Pipeline Stage 6: Report — generate and download analysis report."""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.pages.pipeline.help_content import render_help_panel
from sequence_analyzer.app.pages.pipeline.navigation import (
    render_config_sidebar,
    render_progress_bar,
)
from sequence_analyzer.app.pages.pipeline.state import get_pipeline_state, mark_complete
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.report import (
    generate_html_report,
    generate_methods_section,
    generate_pdf_report,
)
from sequence_analyzer.core.validation import detect_sequence_type


def _resolve_seq_type(config, state) -> str:
    """Derive the actual sequence type for the methods section."""
    if config.seq_type != "auto":
        return config.seq_type.upper()
    # Majority vote from validated sequences
    sequences = state.valid_sequences or state.sequences
    if not sequences:
        return "DNA"
    rna_count = sum(1 for r in sequences if detect_sequence_type(str(r.seq)) == "RNA")
    return "RNA" if rna_count > len(sequences) // 2 else "DNA"

st.set_page_config(layout="wide", page_title="Pipeline: Report")
apply_styles()
hide_streamlit_chrome()
render_config_sidebar()

state = get_pipeline_state()
state.stage = "Report"

render_progress_bar()
st.title("📄 Stage 6: Generate Report")


render_help_panel("Report")

# --- User inputs ---
state.report_title = st.text_input(
    "Report Title",
    value=state.report_title or "Sequence Analysis Report",
    key="report_title_input",
)
state.report_notes = st.text_area(
    "Notes (optional)",
    value=state.report_notes,
    height=100,
    key="report_notes_input",
)

# --- Determine what's available ---
config = state.config
has_qc = state.qc_results is not None
has_analysis = state.analysis_results is not None
has_alignment = state.alignment_result is not None
has_tree = state.tree_result is not None

completed = state.completed_stages
st.markdown("**Stages included in report:**")
for stage in completed:
    st.markdown(f"- {stage}")

if not completed:
    st.warning("Complete at least one stage before generating a report.")
    render_footer()
    st.stop()

# --- Generate ---
if st.button("Generate Report"):
    if not state.sequences:
        st.error("No sequences available. Complete the Ingest stage first.")
        st.stop()
    seq_count = len(state.sequences)
    qc_pass_count = len(state.valid_sequences) if state.valid_sequences else seq_count

    methods = generate_methods_section(
        seq_count=seq_count,
        qc_pass_count=qc_pass_count,
        seq_type=_resolve_seq_type(config, state),
        alignment_method=config.alignment_method,
        tree_method=config.tree_method,
        bootstrap_replicates=config.bootstrap_replicates,
        gc_window_size=config.gc_window_size,
    )

    alignment_summary = f"Method: {config.alignment_method}, Identity: {state.alignment_result.avg_identity:.1f}%" if has_alignment else ""
    tree_newick = state.tree_result.newick if has_tree else ""

    html = generate_html_report(
        title=state.report_title,
        methods_section=methods,
        qc_table=state.qc_results,
        analysis_table=state.analysis_results,
        figures_html=state.figures or None,
        alignment_summary=alignment_summary,
        tree_newick=tree_newick,
        notes=state.report_notes,
    )

    mark_complete("Report")
    st.success("Report generated!")

    # --- Preview ---
    with st.expander("Preview Report", expanded=False):
        import streamlit.components.v1 as components

        components.html(html, height=600, scrolling=True)

    # --- Downloads ---
    col1, col2 = st.columns(2)
    with col1:
        st.download_button(
            "Download HTML Report",
            html,
            file_name="analysis_report.html",
            mime="text/html",
        )
    with col2:
        with st.spinner("Generating PDF..."):
            pdf_bytes = generate_pdf_report(html)
        st.download_button(
            "Download PDF Report",
            pdf_bytes,
            file_name="analysis_report.pdf",
            mime="application/pdf",
        )

render_footer()
