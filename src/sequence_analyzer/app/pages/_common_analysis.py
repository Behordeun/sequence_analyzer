"""Shared logic for DNA/RNA sequence analysis pages."""

import plotly.express as px
import streamlit as st

from sequence_analyzer.core.analysis import analyze_sequences, compute_gc_skew
from sequence_analyzer.core.motifs import scan_motifs
from sequence_analyzer.core.validation import validate_sequences
from sequence_analyzer.io.parsers import parse_sequence_file


def render_sequence_analysis_page(
    page_key: str,
    title: str,
    seq_type: str,
    is_rna: bool,
    upload_label: str,
    motif_placeholder: str,
    default_motif: str,
    download_prefix: str,
) -> None:
    """Render a complete sequence analysis page with upload, analysis, and export."""
    seq_key = f"{page_key}_sequences"
    df_key = f"{page_key}_sequence_df"

    st.title(title)

    st.session_state.setdefault(seq_key, [])
    st.session_state.setdefault(df_key, None)

    uploaded_file = st.file_uploader(upload_label, type=["fasta", "txt", "rtf"])
    if uploaded_file:
        content = uploaded_file.getvalue().decode("utf-8")
        records = parse_sequence_file(content, format_hint="auto")
        valid = validate_sequences(records, seq_type=seq_type)
        st.session_state[seq_key] = valid
        st.success(f"{len(valid)} sequences loaded.")

    sequences = st.session_state[seq_key]

    if sequences and st.checkbox("Preview Sequences"):
        st.subheader("Uploaded Sequences")
        for record in sequences:
            st.write(f"{record.id} - {len(record.seq)} bp")
            st.code(str(record.seq), language="text")

    if sequences and st.button("Analyze"):
        st.session_state[df_key] = analyze_sequences(sequences, is_rna=is_rna)
        st.success(f"{seq_type} Analysis Completed!")

    df = st.session_state[df_key]
    if df is None:
        return

    if st.checkbox("Show Analysis Table"):
        st.dataframe(df)

    if st.checkbox("Show Visualizations"):
        st.subheader("GC Content Distribution")
        fig = px.histogram(df, x="GC_Content", nbins=10, title="GC Content Distribution")
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Sequence Length Distribution")
        fig = px.histogram(df, x="Length", nbins=10, title="Sequence Length Distribution")
        st.plotly_chart(fig, use_container_width=True)

    if st.checkbox("GC Skew Analysis"):
        window_size = st.slider("Select Window Size", 10, 200, 50, step=10)
        df_skew = compute_gc_skew(sequences, window_size=window_size)
        if not df_skew.empty:
            fig = px.line(
                df_skew,
                x="Position",
                y="GC_Skew",
                color="ID",
                title="GC Skew Across All Sequences",
            )
            st.plotly_chart(fig, use_container_width=True)

    if st.checkbox("Motif Scanner"):
        motif = st.text_input(motif_placeholder, value=default_motif)
        if motif.strip():
            df_motif = scan_motifs(sequences, motif.strip().upper())
            st.subheader("Motif Scan Results")
            st.dataframe(df_motif)

    file_format = st.selectbox("Select Download Format", [".csv", ".txt", ".rtf"])
    if file_format == ".csv":
        file_data = df.to_csv(index=False)
    elif file_format == ".txt":
        file_data = df.to_string(index=False)
    else:
        file_data = "{\\rtf1\\ansi\n" + df.to_string(index=False) + "\n}"

    st.download_button(
        f"Download ({file_format.upper()})",
        data=file_data,
        file_name=f"{download_prefix}{file_format}",
    )
