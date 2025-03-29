import datetime
import json
import os
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

SESSION_PATH = "last_alignment_session.json"
if os.path.exists(SESSION_PATH):
    with open(SESSION_PATH) as f:
        try:
            session_data = json.load(f)
            st.session_state["sequences"] = list(
                SeqIO.parse(StringIO("".join(session_data["sequences"])), "fasta")
            )
            st.session_state["aligned"] = session_data["result"]
            st.session_state["metadata"] = session_data.get("metadata", {})
            st.success("✅ Auto-resumed last session")
        except Exception:
            pass

# User Input Options
option = st.selectbox(
    "Choose Input Method",
    ["Upload Sequence File(s)", "Enter Assertion Numbers", "Reload Saved Session"],
)
sequence_type = st.radio("🔬 Select Sequence Type", ["DNA", "RNA", "Protein"])
alignment_method = st.radio(
    "🧬 Select Alignment Method", ["Pairwise Alignment", "Multiple Sequence Alignment"]
)

# Input sequence validators
VALID_DNA = re.compile(r"^[ACGTURYKMSWBDHVN-]+$")
VALID_PROTEIN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO-]+$")


def clean_and_validate(record, seq_type):
    seq = str(record.seq).upper().replace("\n", "").replace(" ", "")
    record.seq = Seq(seq)
    pattern = VALID_DNA if seq_type != "Protein" else VALID_PROTEIN
    if not pattern.match(seq):
        st.warning(f"⚠️ Sequence '{record.id}' contains invalid characters: {seq}")
    return record


# Capture sequences into session
if "sequences" not in st.session_state:
    st.session_state["sequences"] = []

# Upload option
if option == "Upload Sequence File(s)":
    uploaded_files = st.file_uploader(
        "Upload Sequence Files (FASTA format)",
        type=["fasta", "txt", "rtf"],
        accept_multiple_files=True,
    )
    if uploaded_files:
        uploaded_sequences = []
        for uploaded_file in uploaded_files:
            fasta_content = convert_to_fasta(uploaded_file)
            fasta_records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
            for record in fasta_records:
                uploaded_sequences.append(clean_and_validate(record, sequence_type))
        st.session_state["sequences"] = uploaded_sequences
        st.success(f"{len(uploaded_sequences)} sequences uploaded!")

# GenBank option
elif option == "Enter Assertion Numbers":
    assertion_numbers = st.text_area(
        "Enter Assertion Numbers (one per line)"
    ).splitlines()
    if assertion_numbers and st.button("🔍 Fetch Sequences from GenBank"):
        genbank_sequences = [fetch_sequence(acc) for acc in assertion_numbers]
        fetched_sequences = [
            SeqIO.read(StringIO(seq), "fasta") for seq in genbank_sequences if seq
        ]
        for record in fetched_sequences:
            clean_and_validate(record, sequence_type)
        st.session_state["sequences"] = fetched_sequences
        st.success(f"{len(fetched_sequences)} sequences retrieved!")

# Reload saved .json or alignment
elif option == "Reload Saved Session":
    session_file = st.file_uploader(
        "📁 Upload Saved File (.json, .nex, .phy)", type=["json", "nex", "phy"]
    )
    if session_file:
        filename = session_file.name
        if filename.endswith(".json"):
            session_data = json.load(session_file)
            if st.checkbox("👀 Preview Session Before Loading"):
                st.code(session_data["result"], language="html")
            if st.button("🔁 Load Session"):
                st.session_state["sequences"] = list(
                    SeqIO.parse(StringIO("".join(session_data["sequences"])), "fasta")
                )
                st.session_state["aligned"] = session_data["result"]
                st.session_state["metadata"] = session_data.get("metadata", {})
                st.success("✅ Session reloaded successfully!")
        else:
            fmt = "nexus" if filename.endswith(".nex") else "phylip"
            try:
                alignment = AlignIO.read(session_file, fmt)
                records = [clean_and_validate(r, sequence_type) for r in alignment]
                st.session_state["sequences"] = records
                st.success(f"✅ {fmt.upper()} alignment loaded!")
                st.info("📛 Taxon Names:")
                for record in records:
                    st.write(f"🔹 {record.id} - {len(record.seq)} bp")
            except Exception as e:
                st.error(f"❌ Failed to load alignment: {e}")

# Show preview before alignment
if st.session_state["sequences"]:
    if st.checkbox("👁️ Preview Sequences"):
        st.subheader("📄 Identified Sequences:")
        for seq in st.session_state["sequences"]:
            st.write(f"✅ {seq.id} - {len(seq.seq)} bp")


# Alignment logic
def highlight_mismatches(seq1, seq2):
    return "<br>".join(
        [
            "".join(
                f"<span class='{'match' if a == b else 'mismatch'}'>{a}</span>"
                for a, b in zip(seq1, seq2)
            ),
            "".join(
                f"<span class='{'match' if a == b else 'mismatch'}'>{b}</span>"
                for a, b in zip(seq1, seq2)
            ),
        ]
    )


def align_pairwise(sequences, matrix):
    aligner = PairwiseAligner()
    aligner.mode = st.selectbox(
        "🔧 Select Alignment Mode", ["global", "local", "semi-global"]
    )
    try:
        aligner.substitution_matrix = substitution_matrices.load(matrix)
    except Exception:
        st.warning("⚠️ Could not load matrix.")
    results = []
    for i in range(len(sequences) - 1):
        alignment = aligner.align(sequences[i].seq, sequences[i + 1].seq)[0]
        score = alignment.score
        seqA, seqB = str(alignment.target), str(alignment.query)
        identity = sum(a == b for a, b in zip(seqA, seqB)) / len(seqA) * 100
        block = highlight_mismatches(seqA, seqB)
        results.append(
            f"<b>{sequences[i].id} vs {sequences[i+1].id}</b><br><b>Score:</b> {score:.2f}, <b>Identity:</b> {identity:.1f}%<br>{block}<br><br>"
        )
    return "<br>".join(results)


def align_msa(sequences):
    max_len = max(len(s.seq) for s in sequences)
    padded = [
        SeqRecord(Seq(str(s.seq).ljust(max_len, "-")), id=s.id) for s in sequences
    ]
    alignment = MultipleSeqAlignment(padded)
    motif = motifs.create([s.seq for s in alignment])
    consensus = motif.consensus
    clustal_io = StringIO()
    AlignIO.write(alignment, clustal_io, "clustal")
    st.session_state["clustal"] = clustal_io.getvalue()
    return (
        "<br>".join([f">{s.id}<br>{s.seq}" for s in padded])
        + f"<br>>Consensus<br>{consensus}"
    )


# Align Button
if st.session_state["sequences"]:
    if st.button("🔄 Align Sequence"):
        matrix_name = st.selectbox(
            "Choose Scoring Matrix",
            (
                ["BLOSUM62", "BLOSUM80"]
                if sequence_type == "Protein"
                else ["NUC.4.4", "IDENTITY"]
            ),
        )
        if sequence_type == "DNA":
            reverse_view = st.checkbox("Show reverse complement")
            if reverse_view:
                for seq in st.session_state["sequences"]:
                    seq.seq = seq.seq.reverse_complement()
                st.success("🔁 Reverse complement applied")

        if alignment_method == "Pairwise Alignment":
            aligned_html = align_pairwise(st.session_state["sequences"], matrix_name)
        else:
            aligned_html = align_msa(st.session_state["sequences"])

        st.session_state["aligned"] = aligned_html
        st.session_state["metadata"] = {
            "matrix": matrix_name,
            "method": alignment_method,
            "type": sequence_type,
            "timestamp": str(datetime.datetime.now()),
        }
        with open(SESSION_PATH, "w") as f:
            json.dump(
                {
                    "sequences": [
                        s.format("fasta") for s in st.session_state["sequences"]
                    ],
                    "result": aligned_html,
                    "metadata": st.session_state["metadata"],
                },
                f,
            )
        st.success("✅ Alignment complete!")

# Show results
if "aligned" in st.session_state:
    if st.checkbox("📄 Show Alignment Output"):
        if st.checkbox("📑 Show Alignment Metadata"):
            meta = st.session_state.get("metadata", {})
            st.markdown(
                f"""
                <div style='background-color:#f8f8f8; padding:10px;'>
                    <b>Method:</b> {meta.get("method", "N/A")}<br>
                    <b>Matrix:</b> {meta.get("matrix", "N/A")}<br>
                    <b>Type:</b> {meta.get("type", "N/A")}<br>
                    <b>Time:</b> {meta.get("timestamp", "N/A")}
                </div>
            """,
                unsafe_allow_html=True,
            )

        st.markdown(st.session_state["aligned"], unsafe_allow_html=True)

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
                "sequences": [s.format("fasta") for s in st.session_state["sequences"]],
                "result": st.session_state["aligned"],
                "metadata": st.session_state.get("metadata", {}),
            }
            json_str = json.dumps(session, indent=2)
            st.download_button(
                "📁 Download Session (.json)",
                json_str,
                file_name=f"{session_name}.json",
            )

# Footer
st.markdown("---")
st.markdown(
    """
<p style="text-align:center;font-size:14px">
Developed by <a href="https://github.com/Behordeun">Behordeun</a> and <a href="https://github.com/bollergene">Bollergene</a><br>
📞 +2348108316393 | © Behordeun 2025
</p>
""",
    unsafe_allow_html=True,
)
