import datetime
import json
import re
from collections import Counter
from io import StringIO

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

st.set_page_config(layout="wide")
st.markdown(
    """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    span.match { color: green; font-weight: bold; }
    span.mismatch { color: red; font-weight: bold; }
    span.motif { background-color: yellow; font-weight: bold; }
    </style>
""",
    unsafe_allow_html=True,
)

st.title("🔗 Sequence Alignment, GC Skew, and Motif Analysis")

st.session_state.setdefault("sequences", [])
st.session_state.setdefault("aligned", "")
st.session_state.setdefault("metadata", {})
st.session_state.setdefault("errors", [])

option = st.selectbox(
    "Choose Input Method",
    ["Upload Sequence File(s)", "Enter Assertion Numbers", "Reload Saved Session"],
)
seq_type_choice = st.radio("🔬 Sequence Type", ["Auto", "DNA", "RNA", "Protein"])
alignment_method = st.radio(
    "🧬 Alignment Method", ["Pairwise Alignment", "Multiple Sequence Alignment"]
)

ALPHABETS = {
    "DNA": "ACGTURYKMSWBDHVN-",
    "RNA": "ACGUURYKMSWBDHVN-",
    "Protein": "ACDEFGHIKLMNPQRSTVWYBXZJUO-",
}
REPLACEMENTS = {"N": "-", "X": "-", "J": "-", "U": "U", "O": "-"}


# --- Helper Functions ---
def auto_detect_seq_type(seq):
    seq = seq.upper()
    counts = {t: sum(1 for x in seq if x in ALPHABETS[t]) for t in ALPHABETS}
    return max(counts, key=counts.get)


def clean_and_validate(record):
    errors = []
    seq = str(record.seq).upper().replace(" ", "").replace("\n", "")
    detected = (
        auto_detect_seq_type(seq) if seq_type_choice == "Auto" else seq_type_choice
    )
    cleaned = "".join(REPLACEMENTS.get(n, n) for n in seq)
    if not all(n in ALPHABETS[detected] for n in cleaned):
        errors.append(f"❌ {record.id} skipped due to invalid characters")
        return None, errors, detected
    record.seq = Seq(cleaned)
    return record, errors, detected


def highlight_alignment(seq1, seq2, motif_regex):
    html1, html2 = [], []
    for a, b in zip(seq1, seq2):
        tag = "match" if a == b else "mismatch"
        html1.append(f"<span class='{tag}'>{a}</span>")
        html2.append(f"<span class='{tag}'>{b}</span>")
    return "".join(html1) + "<br>" + "".join(html2)


def motif_scan(sequences, pattern):
    results = []
    for rec in sequences:
        seq = str(rec.seq)
        matches = [m.start() + 1 for m in re.finditer(pattern, seq)]
        results.append(
            {
                "ID": rec.id,
                "Pattern": pattern,
                "Hits": len(matches),
                "Positions": matches,
            }
        )
    return pd.DataFrame(results)


def gc_skew_plot(sequences):
    data = []
    for rec in sequences:
        seq = str(rec.seq).upper()
        window = 100
        gc_skew = []
        for i in range(0, len(seq) - window + 1):
            win = seq[i : i + window]
            g = win.count("G")
            c = win.count("C")
            skew = (g - c) / (g + c) if (g + c) > 0 else 0
            gc_skew.append(skew)
        data.append((rec.id, gc_skew))
    fig = px.line(title="GC Skew", labels={"value": "GC Skew", "index": "Position"})
    for rec_id, skew in data:
        fig.add_scatter(x=list(range(len(skew))), y=skew, name=rec_id)
    return fig


def base_composition_summary(sequences):
    summary = []
    for rec in sequences:
        seq = str(rec.seq).upper()
        a, t, g, c = seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")
        total = len(seq.replace("-", ""))
        summary.append(
            {
                "ID": rec.id,
                "Length": total,
                "A%": round(a / total * 100, 2),
                "T%": round(t / total * 100, 2),
                "G%": round(g / total * 100, 2),
                "C%": round(c / total * 100, 2),
                "GC%": round((g + c) / total * 100, 2),
            }
        )
    return pd.DataFrame(summary)


def batch_export_metadata(sequences):
    data = [
        {"ID": rec.id, "Description": rec.description, "Length": len(rec.seq)}
        for rec in sequences
    ]
    return pd.DataFrame(data)


def align_pairwise(seqs, matrix, mode, motif_regex):
    aligner = PairwiseAligner()
    aligner.mode = mode
    try:
        aligner.substitution_matrix = substitution_matrices.load(matrix)
    except:
        st.warning("⚠️ Matrix not found")
    blocks = []
    for i in range(len(seqs) - 1):
        a, b = seqs[i], seqs[i + 1]
        aln = aligner.align(a.seq, b.seq)[0]
        sa, sb = str(aln.target), str(aln.query)
        identity = sum(x == y for x, y in zip(sa, sb)) / len(sa) * 100
        blocks.append(
            f"<b>{a.id} vs {b.id}</b><br>Score: {aln.score:.2f}, Identity: {identity:.1f}%<br>{highlight_alignment(sa, sb, motif_regex)}<br><br>"
        )
    return "".join(blocks)


def align_msa(seqs):
    max_len = max(len(r.seq) for r in seqs)
    padded = [SeqRecord(Seq(str(r.seq).ljust(max_len, "-")), id=r.id) for r in seqs]
    alignment = MultipleSeqAlignment(padded)
    consensus = "".join(
        "-" if all(b == "-" for b in col) else Counter(col).most_common(1)[0][0]
        for col in zip(*[str(rec.seq) for rec in padded])
    )
    clustal_io = StringIO()
    AlignIO.write(alignment, clustal_io, "clustal")
    st.session_state["clustal"] = clustal_io.getvalue()
    return (
        "<br>".join([f">{r.id}<br>{r.seq}" for r in padded])
        + f"<br>>Consensus<br>{consensus}"
    )


# --- File Upload Section ---
if option == "Upload Sequence File(s)":
    uploaded = st.file_uploader(
        "Upload FASTA", type=["fasta", "txt", "rtf"], accept_multiple_files=True
    )
    if uploaded:
        seqs, errors = [], []
        for f in uploaded:
            content = f.read().decode("utf-8")
            for rec in SeqIO.parse(StringIO(content), "fasta"):
                r, e, _ = clean_and_validate(rec)
                if r:
                    seqs.append(r)
                errors += e
        st.session_state.update(sequences=seqs, errors=errors)
        st.success(f"✅ {len(seqs)} valid sequences loaded")

# --- Alignment Section ---
if st.session_state["sequences"]:
    matrix = st.selectbox(
        "Scoring Matrix",
        (
            ["BLOSUM62", "BLOSUM80"]
            if seq_type_choice == "Protein"
            else ["NUC.4.4", "IDENTITY"]
        ),
    )
    mode = st.selectbox("Alignment Mode", ["global", "local", "semi-global"])
    motif_input = st.text_input("Motif Pattern (Regex)", value="ATG")
    motif_regex = re.compile(motif_input, re.IGNORECASE)

    if st.button("🔄 Align Sequences"):
        if seq_type_choice == "DNA" and st.checkbox("🔁 Reverse Complement"):
            for s in st.session_state["sequences"]:
                s.seq = s.seq.reverse_complement()

        aligned = (
            align_pairwise(st.session_state["sequences"], matrix, mode, motif_regex)
            if alignment_method == "Pairwise Alignment"
            else align_msa(st.session_state["sequences"])
        )

        st.session_state.update(
            aligned=aligned,
            metadata={
                "matrix": matrix,
                "method": alignment_method,
                "type": seq_type_choice,
                "mode": mode,
                "motif": motif_input,
                "timestamp": str(datetime.datetime.now()),
            },
        )
        st.success("✅ Alignment completed")

# --- Output Section ---
if st.session_state["aligned"]:
    st.subheader("📄 Aligned Sequences")
    st.markdown(st.session_state["aligned"], unsafe_allow_html=True)

    st.subheader("📊 GC Skew Visualization")
    st.plotly_chart(gc_skew_plot(st.session_state["sequences"]))

    st.subheader("🧬 Motif Scan Table")
    motif_df = motif_scan(st.session_state["sequences"], motif_input)
    st.dataframe(motif_df)
    st.download_button(
        "📥 Download Motif Table",
        motif_df.to_csv(index=False),
        file_name="motif_hits.csv",
    )

    st.subheader("🎯 Base Composition Summary")
    base_df = base_composition_summary(st.session_state["sequences"])
    st.dataframe(base_df)
    st.download_button(
        "📥 Download Composition",
        base_df.to_csv(index=False),
        file_name="base_composition.csv",
    )

    st.subheader("📂 Export Metadata")
    meta_df = batch_export_metadata(st.session_state["sequences"])
    st.dataframe(meta_df)
    st.download_button(
        "📤 Export Metadata (CSV)",
        meta_df.to_csv(index=False),
        file_name="sequence_metadata.csv",
    )

    fmt = st.selectbox(
        "\ud83d\udcc5 Download Alignment Format", [".fasta", ".txt", ".clustal", ".nex"]
    )
    raw = st.session_state["aligned"].replace("<br>", "\n")

    if fmt == ".clustal":
        out = st.session_state.get("clustal", raw)
    elif fmt == ".nex":
        # Convert aligned sequences into MultipleSeqAlignment object again
        max_len = max(len(r.seq) for r in st.session_state["sequences"])
        padded = [
            SeqRecord(Seq(str(r.seq).ljust(max_len, "-")), id=r.id)
            for r in st.session_state["sequences"]
        ]
        alignment = MultipleSeqAlignment(padded)
        nexus_io = StringIO()
        AlignIO.write(alignment, nexus_io, "nexus")
        out = nexus_io.getvalue()
    else:
        out = raw

    session_name = st.text_input("💾 Session Name", "alignment_session")
    if st.button("💾 Save Session"):
        json_obj = {
            "sequences": [s.format("fasta") for s in st.session_state["sequences"]],
            "result": st.session_state["aligned"],
            "metadata": st.session_state["metadata"],
        }
        st.download_button(
            "📁 Export Session (.json)",
            json.dumps(json_obj, indent=2),
            file_name=f"{session_name}.json",
        )

# --- Footer ---
st.markdown("---")
st.markdown(
    """
<p style="text-align:center;font-size:14px">
Developed by <a href="https://github.com/Behordeun">Behordeun</a> & <a href="https://github.com/bollergene">Bollergene</a><br>
📞 +2348108316393 | © Behordeun 2025
</p>
""",
    unsafe_allow_html=True,
)
