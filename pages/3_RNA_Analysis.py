import json
import re
from io import BytesIO, StringIO

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq

from style_css import style
from utils import analyze_sequences, visualize_results

style()

st.markdown(
    """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
""",
    unsafe_allow_html=True,
)

st.title("🧬 RNA Sequence Analysis")

if (
    "last_page" not in st.session_state
    or st.session_state["last_page"] != "RNA_Analysis"
):
    st.session_state.clear()
    st.session_state["last_page"] = "RNA_Analysis"

if "sequence_df" not in st.session_state:
    st.session_state["sequence_df"] = None
if "sequences" not in st.session_state:
    st.session_state["sequences"] = []

VALID_RNA = re.compile(r"^[ACGUURYKMSWBDHVN-]+$")


def validate_sequence(record):
    seq = str(record.seq).upper().replace(" ", "").replace("\n", "")
    record.seq = Seq(seq)
    if not VALID_RNA.match(seq):
        st.warning(f"⚠️ Invalid characters found in {record.id}")
    return record


uploaded_file = st.file_uploader(
    "📁 Upload RNA Sequence File", type=["fasta", "txt", "rtf"]
)
if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8")
    records = list(SeqIO.parse(StringIO(content), "fasta-pearson"))
    st.session_state["sequences"] = [validate_sequence(r) for r in records]
    st.success(f"✅ {len(records)} sequences loaded.")

if st.session_state["sequences"] and st.checkbox("🔁 Show Reverse Complement"):
    for seq in st.session_state["sequences"]:
        seq.seq = seq.seq.reverse_complement()
    st.success("🔁 Reverse complements applied.")

if st.session_state["sequences"] and st.checkbox("👁️ Preview Sequences"):
    st.subheader("📄 Uploaded Sequences")
    for record in st.session_state["sequences"]:
        st.write(f"🔹 {record.id} - {len(record.seq)} bp")

if st.session_state["sequences"]:
    if st.button("🔬 Analyze"):
        st.session_state["sequence_df"] = analyze_sequences(
            st.session_state["sequences"], is_rna=True
        )
        st.success("✅ RNA Analysis Completed!")

if st.session_state["sequence_df"] is not None:
    df = st.session_state["sequence_df"]

    if st.checkbox("📊 Show Analysis Table"):
        st.dataframe(df)

    if st.checkbox("📈 Show Visualizations"):
        visualize_results(df)

    # Unified GC Skew Plot for all sequences
    if st.checkbox("🧪 GC Skew Analysis"):
        window_size = st.slider("Select Window Size", 10, 200, 50, step=10)
        rows = []

        for record in st.session_state["sequences"]:
            seq = str(record.seq).upper()
            for i in range(len(seq) - window_size + 1):
                win = seq[i : i + window_size]
                g, c = win.count("G"), win.count("C")
                skew = (g - c) / (g + c) if (g + c) else 0
                rows.append({"Sequence ID": record.id, "Position": i, "GC Skew": skew})

        skew_df = pd.DataFrame(rows)
        st.subheader("📉 GC Skew")
        fig = px.line(
            skew_df,
            x="Position",
            y="GC Skew",
            color="Sequence ID",
            labels={"GC Skew": "GC Skew Value"},
            title="GC Skew Across All Sequences",
        )
        st.plotly_chart(fig, use_container_width=True)

        # --- Multi-format download ---
        export_fmt = st.selectbox(
            "📥 Download GC Skew Plot As", ["PNG", "PDF", "SVG", "JPEG", "HTML"]
        )
        buffer = BytesIO()

        if export_fmt == "HTML":
            html = fig.to_html(full_html=False, include_plotlyjs="cdn")
            st.download_button(
                "⬇️ Download HTML",
                html,
                file_name="gc_skew.html",
                mime="text/html",
            )
        else:
            try:
                fig.write_image(buffer, format=export_fmt.lower())
                mime = (
                    "application/pdf"
                    if export_fmt == "PDF"
                    else f"image/{export_fmt.lower()}"
                )
                st.download_button(
                    f"⬇️ Download {export_fmt}",
                    buffer.getvalue(),
                    file_name=f"gc_skew.{export_fmt.lower()}",
                    mime=mime,
                )
            except Exception as e:
                st.error(f"Export failed: {e}")

    if st.checkbox("🧬 Motif Scanner"):
        motif = st.text_input("Enter RNA motif (e.g., AUG or UUU)", value="AUG")
        motif = motif.upper().strip()
        if motif:
            results = []
            for record in st.session_state["sequences"]:
                matches = [
                    i
                    for i in range(len(record.seq) - len(motif) + 1)
                    if record.seq[i : i + len(motif)] == motif
                ]
                results.append(
                    {
                        "Sequence": record.id,
                        "Matches": len(matches),
                        "Positions": matches,
                    }
                )
            st.subheader("🔎 Motif Scan Results")
            for r in results:
                st.write(
                    f"🧬 {r['Sequence']}: {r['Matches']} match(es) at positions {r['Positions']}"
                )

    file_format = st.selectbox("📂 Select Download Format", [".csv", ".txt", ".rtf"])
    file_name = f"rna_analysis{file_format}"
    if file_format == ".csv":
        file_data = df.to_csv(index=False)
    elif file_format == ".txt":
        file_data = df.to_string(index=False)
    elif file_format == ".rtf":
        file_data = "{\\rtf1\\ansi\n" + df.to_string(index=False) + "\n}"

    st.download_button(
        f"📥 Download ({file_format.upper()})", data=file_data, file_name=file_name
    )

    if st.checkbox("💾 Save Session"):
        session_name = st.text_input(
            "Enter session filename", value="rna_analysis_session"
        )
        session_obj = {
            "data": df.to_dict(),
            "sequences": [
                record.format("fasta") for record in st.session_state["sequences"]
            ],
        }
        st.download_button(
            "📥 Export Session (.json)",
            data=json.dumps(session_obj, indent=2),
            file_name=f"{session_name}.json",
        )

st.markdown("---")
st.markdown(
    """
<p style="text-align:center;font-size:14px">
    Developed by <a href="https://github.com/Behordeun">Behordeun</a> and 
    <a href="https://github.com/bollergene">Bollergene</a><br>
    📞 +2348108316393 | © Behordeun 2025
</p>
""",
    unsafe_allow_html=True,
)
