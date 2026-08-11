"""Pipeline Stage 1: Ingest — upload files or fetch from GenBank."""

import os

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.pages.pipeline.help_content import render_help_panel
from sequence_analyzer.app.pages.pipeline.navigation import (
    render_config_sidebar,
    render_progress_bar,
    render_stage_nav,
)
from sequence_analyzer.app.pages.pipeline.state import (
    apply_organism_mode,
    get_pipeline_state,
    mark_complete,
)
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.genbank import fetch_sequence
from sequence_analyzer.core.organisms import load_organism_modes
from sequence_analyzer.io.parsers import parse_sequence_file

st.set_page_config(layout="wide", page_title="Pipeline: Ingest")
apply_styles()
hide_streamlit_chrome()
render_config_sidebar()

state = get_pipeline_state()
state.stage = "Ingest"

render_progress_bar()
st.title("📥 Stage 1: Ingest Sequences")


render_help_panel("Ingest")

# --- Organism Mode ---
st.subheader("Organism Mode")

modes = load_organism_modes()
mode_slugs = sorted(modes.keys())
mode_display = {slug: modes[slug].display_name for slug in mode_slugs}

current_mode = state.config.organism_mode
current_index = mode_slugs.index(current_mode) if current_mode in mode_slugs else 0

selected_mode = st.selectbox(
    "Select organism",
    options=mode_slugs,
    index=current_index,
    format_func=lambda slug: mode_display[slug],
    help="Pre-configures QC thresholds, contamination screening, and motif patterns for the selected organism.",
)

if selected_mode != state.config.organism_mode:
    apply_organism_mode(selected_mode)
    st.rerun()

# Show active mode details
active_mode = modes[state.config.organism_mode]
if active_mode.description:
    st.caption(active_mode.description)
if active_mode.reference_accession:
    st.caption(
        f"Suggested reference: {active_mode.reference_accession} ({active_mode.reference_description})"
    )

st.markdown("---")

# --- Input method ---
input_method = st.radio(
    "Input Method", ["Upload File(s)", "GenBank Accession Numbers"], horizontal=True
)

if input_method == "Upload File(s)":
    if files := st.file_uploader(
        "Upload sequence file(s)",
        type=["fasta", "txt", "rtf", "phy", "nex", "aln"],
        accept_multiple_files=True,
    ):
        all_records = []
        warnings = []
        for f in files:
            size_mb = f.size / (1024 * 1024)
            if size_mb > 50:
                st.error(f"{f.name} exceeds 50MB limit ({size_mb:.1f}MB). Skipped.")
                continue
            try:
                content = f.getvalue().decode("utf-8")
            except UnicodeDecodeError:
                st.error(f"{f.name} is not a valid UTF-8 text file. Skipped.")
                continue
            records = parse_sequence_file(content, format_hint="auto")
            if not records:
                warnings.append(f"No sequences found in {f.name}")
            all_records.extend(records)

        if all_records:
            state.sequences = all_records
            summary = f"{len(all_records)} sequences parsed from {len(files)} file(s)"
            state.ingest_summary = summary
            mark_complete("Ingest")
            st.success(summary)
        if warnings:
            for w in warnings:
                st.warning(w)

else:
    acc_input = st.text_area("Enter accession numbers (one per line)", height=120)
    genbank_email = st.text_input(
        "Email for NCBI API",
        value=os.environ.get("GENBANK_EMAIL", ""),
        help="Required by NCBI usage policy.",
    )
    if st.button("Fetch Sequences"):
        if not genbank_email.strip():
            st.error("Email is required for GenBank API access.")
        elif not acc_input.strip():
            st.error("Enter at least one accession number.")
        else:
            fetched = []
            for acc in acc_input.strip().splitlines():
                acc = acc.strip()
                if not acc:
                    continue
                try:
                    record = fetch_sequence(acc, email=genbank_email.strip())
                    fetched.append(record)
                except (ValueError, ConnectionError) as e:
                    st.warning(f"Could not fetch {acc}: {e}")

            if fetched:
                state.sequences = fetched
                summary = f"{len(fetched)} sequences fetched from GenBank"
                state.ingest_summary = summary
                mark_complete("Ingest")
                st.success(summary)

# --- Preview ---
if state.sequences:
    st.markdown(f"**Loaded:** {len(state.sequences)} sequences")
    if st.checkbox("Preview sequences", key="ingest_preview"):
        for rec in state.sequences[:20]:
            st.text(f">{rec.id} ({len(rec.seq)} bp)")

st.markdown("---")
render_stage_nav()
render_footer()
