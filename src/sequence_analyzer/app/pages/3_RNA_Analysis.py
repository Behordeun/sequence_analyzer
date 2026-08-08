"""RNA Sequence Analysis page. Mirrors DNA page with is_rna=True."""

import plotly.express as px
import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.analysis import analyze_sequences, compute_gc_skew
from sequence_analyzer.core.motifs import scan_motifs
from sequence_analyzer.core.validation import validate_sequences
from sequence_analyzer.io.parsers import parse_sequence_file

st.set_page_config(layout="wide")
apply_styles()
hide_streamlit_chrome()

st.title("🧬 RNA Sequence Analysis")

if st.session_state.get("last_page") != "RNA_Analysis":
    st.session_state.clear()
    st.session_state["last_page"] = "RNA_Analysis"

st.session_state.setdefault("sequence_df", None)
st.session_state.setdefault("sequences", [])

uploaded_file = st.file_uploader("Upload RNA Sequence File", type=["fasta", "txt", "rtf"])
if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8")
    records = parse_sequence_file(content, format_hint="auto")
    valid = validate_sequences(records, seq_type="RNA")
    st.session_state["sequences"] = valid
    st.success(f"{len(valid)} sequences loaded.")

if st.session_state["sequences"] and st.checkbox("Preview Sequences"):
    st.subheader("Uploaded Sequences")
    for record in st.session_state["sequences"]:
        st.write(f"{record.id} - {len(record.seq)} bp")
        st.code(str(record.seq), language="text")

if st.session_state["sequences"] and st.button("Analyze"):
    st.session_state["sequence_df"] = analyze_sequences(st.session_state["sequences"], is_rna=True)
    st.success("RNA Analysis Completed!")

if st.session_state["sequence_df"] is not None:
    df = st.session_state["sequence_df"]

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
        df_skew = compute_gc_skew(st.session_state["sequences"], window_size=window_size)
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
        motif = st.text_input("Enter RNA motif (e.g., AUG or UUU)", value="AUG")
        if motif.strip():
            df_motif = scan_motifs(st.session_state["sequences"], motif.strip().upper())
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
        file_name=f"rna_analysis{file_format}",
    )

render_footer()
