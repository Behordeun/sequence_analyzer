import datetime
import json
import re
from io import StringIO

import streamlit as st
from Bio import AlignIO, SeqIO, motifs
from Bio.Align import MultipleSeqAlignment, PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from style_css import style
from utils import convert_to_fasta, fetch_sequence

style()

st.markdown(
    """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    span.match { color: green; font-weight: bold; }
    span.mismatch { color: red; font-weight: bold; }
    </style>
""",
    unsafe_allow_html=True,
)

st.title("🔗 Sequence Alignment")

# Session state
if "sequences" not in st.session_state:
    st.session_state["sequences"] = []
if "aligned" not in st.session_state:
    st.session_state["aligned"] = ""
if "metadata" not in st.session_state:
    st.session_state["metadata"] = {}

# Input Controls
option = st.selectbox(
    "Choose Input Method",
    ["Upload Sequence File(s)", "Enter Assertion Numbers", "Reload Saved Session"],
)
sequence_type = st.radio("🔬 Select Sequence Type", ["DNA", "RNA", "Protein"])
alignment_method = st.radio(
    "🧬 Select Alignment Method", ["Pairwise Alignment", "Multiple Sequence Alignment"]
)

# Validators
VALID_DNA = re.compile(r"^[ACGTURYKMSWBDHVN-]+$")
VALID_PROTEIN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO-]+$")


def clean_and_validate(record, seq_type):
    seq = str(record.seq).upper().replace(" ", "").replace("\n", "")
    record.seq = Seq(seq)
    pattern = VALID_DNA if seq_type != "Protein" else VALID_PROTEIN
    if not pattern.match(seq):
        st.warning(f"⚠️ Sequence '{record.id}' contains invalid characters.")
    return record


# Upload sequences
if option == "Upload Sequence File(s)":
    uploaded_files = st.file_uploader(
        "Upload FASTA files", type=["fasta", "txt", "rtf"], accept_multiple_files=True
    )
    if uploaded_files:
        sequences = []
        for file in uploaded_files:
            content = convert_to_fasta(file)
            records = SeqIO.parse(StringIO(content), "fasta")
            for r in records:
                sequences.append(clean_and_validate(r, sequence_type))
        st.session_state["sequences"] = sequences
        st.success(f"✅ {len(sequences)} sequences uploaded.")

# GenBank input
elif option == "Enter Assertion Numbers":
    assertion_numbers = st.text_area(
        "Enter GenBank Accession Numbers (one per line)"
    ).splitlines()
    if assertion_numbers and st.button("🔍 Fetch Sequences"):
        records = []
        for acc in assertion_numbers:
            fetched = fetch_sequence(acc)
            try:
                rec = SeqIO.read(StringIO(fetched), "fasta")
                records.append(clean_and_validate(rec, sequence_type))
            except Exception:
                st.error(f"❌ Failed to parse {acc}")
        st.session_state["sequences"] = records
        st.success(f"✅ {len(records)} sequences retrieved.")

# Reload previous session
elif option == "Reload Saved Session":
    session_file = st.file_uploader(
        "📁 Upload Previous Session (.json, .nex, .phy)", type=["json", "nex", "phy"]
    )
    if session_file:
        name = session_file.name
        if name.endswith(".json"):
            session_data = json.load(session_file)
            if st.checkbox("👀 Preview Session Data"):
                st.code(session_data.get("result", "No content"), language="html")
            if st.button("🔄 Load Session"):
                st.session_state["sequences"] = list(
                    SeqIO.parse(StringIO("".join(session_data["sequences"])), "fasta")
                )
                st.session_state["aligned"] = session_data["result"]
                st.session_state["metadata"] = session_data.get("metadata", {})
                st.success("✅ Session Loaded.")
        else:
            try:
                fmt = "nexus" if name.endswith(".nex") else "phylip"
                alignment = AlignIO.read(session_file, fmt)
                st.session_state["sequences"] = [
                    clean_and_validate(r, sequence_type) for r in alignment
                ]
                st.success(f"✅ {fmt.upper()} alignment loaded.")
            except Exception:
                st.error("❌ Error loading alignment.")

# Sequence preview
if st.session_state["sequences"] and st.checkbox("👁️ Preview Sequences"):
    for seq in st.session_state["sequences"]:
        st.write(f"🔹 {seq.id} - {len(seq.seq)} bp")


# Highlight mismatches
def highlight_mismatches(seq1, seq2):
    html = []
    for a, b in zip(seq1, seq2):
        tag = "match" if a == b else "mismatch"
        html.append(f"<span class='{tag}'>{a}</span>")
    html.append("<br>")
    for a, b in zip(seq1, seq2):
        tag = "match" if a == b else "mismatch"
        html.append(f"<span class='{tag}'>{b}</span>")
    return "".join(html)


# Pairwise align
def align_pairwise(sequences, matrix, mode):
    aligner = PairwiseAligner()
    aligner.mode = mode
    try:
        aligner.substitution_matrix = substitution_matrices.load(matrix)
    except:
        st.warning("⚠️ Substitution matrix could not be loaded.")
    results = []
    for i in range(len(sequences) - 1):
        a = sequences[i]
        b = sequences[i + 1]
        alignment = aligner.align(a.seq, b.seq)[0]
        seqA = str(alignment.target)
        seqB = str(alignment.query)
        identity = sum(1 for x, y in zip(seqA, seqB) if x == y) / len(seqA) * 100
        score = alignment.score
        results.append(
            f"<b>{a.id} vs {b.id}</b><br>Score: {score:.2f}, Identity: {identity:.1f}%<br>"
        )
        results.append(highlight_mismatches(seqA, seqB))
        results.append("<br><br>")
    return "".join(results)


# MSA align
def align_msa(sequences):
    max_len = max(len(r.seq) for r in sequences)
    padded = [
        SeqRecord(Seq(str(r.seq).ljust(max_len, "-")), id=r.id) for r in sequences
    ]
    alignment = MultipleSeqAlignment(padded)
    motif = motifs.create([r.seq for r in alignment])
    consensus = motif.consensus
    clustal_io = StringIO()
    AlignIO.write(alignment, clustal_io, "clustal")
    st.session_state["clustal"] = clustal_io.getvalue()
    return (
        "<br>".join([f">{r.id}<br>{r.seq}" for r in padded])
        + f"<br>>Consensus<br>{consensus}"
    )


# Run alignment
if st.session_state["sequences"]:
    matrix_options = (
        ["BLOSUM62", "BLOSUM80"]
        if sequence_type == "Protein"
        else ["NUC.4.4", "IDENTITY"]
    )
    matrix_name = st.selectbox("Select Scoring Matrix", matrix_options)
    align_mode = st.selectbox(
        "🔧 Select Alignment Mode", ["global", "local", "semi-global"]
    )
    if st.button("🔄 Align Sequences"):
        if sequence_type == "DNA" and st.checkbox("🔁 Apply Reverse Complement"):
            for seq in st.session_state["sequences"]:
                seq.seq = seq.seq.reverse_complement()

        if alignment_method == "Pairwise Alignment":
            result_html = align_pairwise(
                st.session_state["sequences"], matrix_name, align_mode
            )
        else:
            result_html = align_msa(st.session_state["sequences"])

        st.session_state["aligned"] = result_html
        st.session_state["metadata"] = {
            "matrix": matrix_name,
            "method": alignment_method,
            "type": sequence_type,
            "mode": align_mode,
            "timestamp": str(datetime.datetime.now()),
        }
        st.success("✅ Alignment Complete!")

# Output tabs
if st.session_state["aligned"]:
    if st.checkbox("📄 Show Alignment Output"):
        st.tabs(["🧬 Alignment View"])

        #with tab1:
        with st.expander("🔬 Expand Alignment Output", expanded=True):
            st.markdown(st.session_state["aligned"], unsafe_allow_html=True)

        # with tab2:
        #     with st.expander("📖 View Raw HTML Code"):
        #         st.code(st.session_state["aligned"], language="html")

        if st.checkbox("📑 Show Metadata"):
            meta = st.session_state["metadata"]
            st.markdown(
                f"""
                <div style='background-color:#f0f0f0; padding:10px;'>
                <b>Type:</b> {meta.get("type")}<br>
                <b>Method:</b> {meta.get("method")}<br>
                <b>Matrix:</b> {meta.get("matrix")}<br>
                <b>Mode:</b> {meta.get("mode")}<br>
                <b>Time:</b> {meta.get("timestamp")}
                </div>
            """,
                unsafe_allow_html=True,
            )

        file_format = st.selectbox(
            "📥 Download Format", [".fasta", ".rtf", ".txt", ".clustal", ".nex", ".phy"]
        )
        data = st.session_state["aligned"].replace("<br>", "\n")
        if file_format == ".clustal":
            data = st.session_state.get("clustal", data)
        elif file_format == ".rtf":
            data = f"{{\\rtf1\\ansi\n{data}\n}}"
        st.download_button(
            "📥 Download Alignment", data, file_name=f"alignment{file_format}"
        )

        session_name = st.text_input("💾 Session Name", value="alignment_session")
        if st.button("💾 Save Alignment Session"):
            session = {
                "sequences": [r.format("fasta") for r in st.session_state["sequences"]],
                "result": st.session_state["aligned"],
                "metadata": st.session_state["metadata"],
            }
            json_data = json.dumps(session, indent=2)
            st.download_button(
                "📁 Download Session (.json)",
                data=json_data,
                file_name=f"{session_name}.json",
            )

# Footer
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
