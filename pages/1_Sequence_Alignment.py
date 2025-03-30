import datetime
import re
from collections import Counter  # ✅ Add this line
from io import BytesIO, StringIO

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import convert_to_fasta, fetch_sequence

# --- Page Settings ---
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

st.title("🔗 Sequence Alignment, GC Skew, Motif Scan & Base Composition")

# --- Session State Initialization ---
for key in ["sequences", "aligned", "metadata"]:
    st.session_state.setdefault(key, [])
for key in ["show_alignment", "show_gc", "show_base", "show_motif", "show_chart"]:
    st.session_state.setdefault(key, False)

# --- Input Options ---
option = st.selectbox(
    "Choose Input Method", ["Upload Sequence File(s)", "Enter Accession Numbers"]
)
seq_type_choice = st.radio("🔬 Sequence Type", ["Auto", "DNA", "RNA", "Protein"])
alignment_method = st.radio(
    "🧬 Alignment Method", ["Pairwise Alignment", "Multiple Sequence Alignment"]
)

# --- Config ---
ALPHABETS = {
    "DNA": "ACGTURYKMSWBDHVN-",
    "RNA": "ACGUURYKMSWBDHVN-",
    "Protein": "ACDEFGHIKLMNPQRSTVWYBXZJUO-",
}
REPLACEMENTS = {"N": "-", "X": "-", "J": "-", "U": "U", "O": "-"}


# --- Helper Functions ---
def auto_detect_seq_type(seq):
    seq = seq.upper()
    return max(ALPHABETS, key=lambda t: sum(1 for x in seq if x in ALPHABETS[t]))


def clean_and_validate(record):
    seq = str(record.seq).upper().replace(" ", "").replace("\n", "")
    detected = (
        auto_detect_seq_type(seq) if seq_type_choice == "Auto" else seq_type_choice
    )
    cleaned = "".join(REPLACEMENTS.get(n, n) for n in seq)
    if not all(n in ALPHABETS[detected] for n in cleaned):
        st.warning(f"❌ {record.id} skipped due to invalid characters")
        return None
    record.seq = Seq(cleaned)
    return record


def highlight_alignment(seq1, seq2):
    html1 = "".join(
        f"<span class='{'match' if a == b else 'mismatch'}'>{a}</span>"
        for a, b in zip(seq1, seq2)
    )
    html2 = "".join(
        f"<span class='{'match' if a == b else 'mismatch'}'>{b}</span>"
        for a, b in zip(seq1, seq2)
    )
    return html1 + "<br>" + html2


def gc_skew_plot(sequences):
    rows = []
    window_size = 100
    for rec in sequences:
        seq = str(rec.seq).upper()
        for i in range(len(seq) - window_size + 1):
            win = seq[i : i + window_size]
            g, c = win.count("G"), win.count("C")
            skew = (g - c) / (g + c) if (g + c) else 0
            rows.append({"ID": rec.id, "Position": i, "GC Skew": skew})
    df = pd.DataFrame(rows)
    fig = px.line(df, x="Position", y="GC Skew", color="ID", title="GC Skew")
    return fig


def base_composition(sequences):
    return pd.DataFrame(
        [
            {
                "ID": s.id,
                "Length": len(s.seq),
                "A%": s.seq.count("A") / len(s.seq) * 100,
                "T%": s.seq.count("T") / len(s.seq) * 100,
                "G%": s.seq.count("G") / len(s.seq) * 100,
                "C%": s.seq.count("C") / len(s.seq) * 100,
                "GC%": (s.seq.count("G") + s.seq.count("C")) / len(s.seq) * 100,
            }
            for s in sequences
        ]
    )


def motif_scan_table(sequences, pattern):
    rows = []
    for rec in sequences:
        seq = str(rec.seq)
        matches = [m.start() + 1 for m in re.finditer(pattern, seq)]
        rows.append(
            {
                "ID": rec.id,
                "Pattern": pattern,
                "Hits": len(matches),
                "Positions": matches,
            }
        )
    return pd.DataFrame(rows)


def align_pairwise(sequences, matrix, mode):
    aligner = PairwiseAligner()
    aligner.mode = mode
    try:
        aligner.substitution_matrix = substitution_matrices.load(matrix)
    except:
        st.warning("⚠️ Matrix not loaded")
    html_blocks = []
    for i in range(len(sequences) - 1):
        a, b = sequences[i], sequences[i + 1]
        aln = aligner.align(a.seq, b.seq)[0]
        identity = (
            sum(x == y for x, y in zip(str(aln.target), str(aln.query)))
            / len(aln.target)
            * 100
        )
        html_blocks.append(
            f"<b>{a.id} vs {b.id}</b><br><b>Score:</b> {aln.score:.2f}, <b>Identity:</b> {identity:.1f}%<br>"
            + highlight_alignment(str(aln.target), str(aln.query))
            + "<br><br>"
        )
    return "".join(html_blocks)


def align_msa(sequences):
    max_len = max(len(r.seq) for r in sequences)
    mol_type = (
        auto_detect_seq_type(str(sequences[0].seq))
        if seq_type_choice == "Auto"
        else seq_type_choice.upper()
    )
    padded = [
        SeqRecord(
            Seq(str(r.seq).ljust(max_len, "-")),
            id=r.id,
            annotations={"molecule_type": mol_type},
        )
        for r in sequences
    ]
    alignment = MultipleSeqAlignment(padded)
    alignment.annotations["molecule_type"] = mol_type
    AlignIO.write(alignment, "sequences/clustal.aln", "clustal")
    AlignIO.write(alignment, "sequences/nexus.nex", "nexus")
    AlignIO.write(alignment, "sequences/phylip.phy", "phylip")
    st.session_state.update(
        {
            "clustal": open("sequences/clustal.aln").read(),
            "nexus": open("sequences/nexus.nex").read(),
            "phylip": open("sequences/phylip.phy").read(),
        }
    )
    seq_matrix = [str(r.seq) for r in padded]
    html_lines = []
    for i in range(len(seq_matrix)):
        line = f"<b>{padded[i].id}</b><br>"
        for j, base in enumerate(seq_matrix[i]):
            ref = seq_matrix[0][j]
            css = "match" if base == ref else "mismatch"
            line += f"<span class='{css}'>{base}</span>"
        html_lines.append(line + "<br>")
    consensus = "".join(
        Counter(col).most_common(1)[0][0] if any(b != "-" for b in col) else "-"
        for col in zip(*seq_matrix)
    )
    consensus_html = "<b>Consensus</b><br>" + "".join(
        (
            f"<span class='match'>{c}</span>"
            if c != "-"
            else "<span class='mismatch'>-</span>"
        )
        for c in consensus
    )
    return "<br>".join(html_lines + [consensus_html])


# --- File Input ---
if option == "Upload Sequence File(s)":
    uploaded = st.file_uploader(
        "Upload sequence file(s)",
        type=["fasta", "txt", "rtf", "phy", "nex"],
        accept_multiple_files=True,
    )
    if uploaded:
        valid = []
        for f in uploaded:
            content = convert_to_fasta(f)
            for rec in SeqIO.parse(StringIO(content), "fasta"):
                val = clean_and_validate(rec)
                if val:
                    valid.append(val)
        st.session_state["sequences"] = valid
        st.success(f"✅ {len(valid)} valid sequences loaded.")

elif option == "Enter Accession Numbers":
    acc_input = st.text_area("🔍 Enter Accession Numbers (one per line)", height=120)

    if st.button("📥 Retrieve Sequences"):
        fetched, raw = [], []
        for acc in acc_input.strip().splitlines():
            try:
                data = fetch_sequence(acc.strip())  # From your utils.py
                record = SeqIO.read(StringIO(data), "fasta")
                fetched.append(record)
                raw.append(data)
            except Exception as e:
                st.warning(f"⚠️ Could not fetch {acc}: {e}")
        st.session_state["sequences"] = fetched
        st.session_state["fetched_raw"] = raw
        st.success(f"✅ Retrieved {len(fetched)} sequence(s) from GenBank.")

    # After retrieval, allow export before alignment
    # Toggle view for retrieved sequences
    if st.session_state["sequences"]:
        st.session_state["show_fetched"] = st.checkbox(
            "👁️ Show Retrieved Sequences", value=True
        )

        if st.session_state["show_fetched"]:
            st.subheader("📄 Retrieved Sequences")
            for rec in st.session_state["sequences"]:
                st.text(f">{rec.id}\n{rec.seq}")

        # Export options
        export_fmt = st.selectbox(
            "📁 Download Fetched Sequences As", [".fasta", ".txt", ".json"]
        )
        if export_fmt == ".fasta":
            out = "\n".join(
                s.format("fasta").strip() for s in st.session_state["sequences"]
            )
        elif export_fmt == ".txt":
            out = "\n".join(f">{s.id}\n{s.seq}" for s in st.session_state["sequences"])
        elif export_fmt == ".json":
            out = {
                "fetched": [r.format("fasta") for r in st.session_state["sequences"]],
                "metadata": {
                    "source": "GenBank",
                    "count": len(st.session_state["sequences"]),
                },
            }
            import json

            out = json.dumps(out, indent=2)

        st.download_button(
            "⬇️ Download Retrieved Data", out, file_name=f"fetched_sequences{export_fmt}"
        )

# --- Alignment Options ---
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
    motif_pattern = st.text_input("🧬 Motif Pattern (Regex)", "ATG")

    # Toggles
    st.session_state["show_alignment"] = st.checkbox("📄 Show Alignment Result", True)
    st.session_state["show_gc"] = st.checkbox("📊 Show GC Skew")
    st.session_state["show_base"] = st.checkbox("🎯 Base Composition Summary")
    st.session_state["show_motif"] = st.checkbox("🧬 Motif Scan Table")

    if st.button("🔄 Align Sequences"):
        aligned = (
            align_pairwise(st.session_state["sequences"], matrix, mode)
            if alignment_method == "Pairwise Alignment"
            else align_msa(st.session_state["sequences"])
        )
        st.session_state["aligned"] = aligned
        st.session_state["metadata"] = {
            "matrix": matrix,
            "mode": mode,
            "type": seq_type_choice,
            "timestamp": str(datetime.datetime.now()),
        }
        st.success("✅ Alignment complete!")

# --- Output ---
if st.session_state["aligned"]:
    if st.session_state["show_alignment"]:
        st.subheader("📄 Alignment Output")
        st.markdown(st.session_state["aligned"], unsafe_allow_html=True)

        # Export styled HTML and RTF
        html_export = f"""
        <html><head>
        <style>
        span.match {{ color: green; font-weight: bold; }}
        span.mismatch {{ color: red; font-weight: bold; }}
        </style></head>
        <body><h3>Alignment</h3>{st.session_state['aligned']}</body></html>
        """
        st.download_button("⬇️ Export HTML", html_export, file_name="alignment.html")
        st.download_button(
            "⬇️ Export RTF",
            f"{{\\rtf1\\ansi\n{st.session_state['aligned'].replace('<br>', '\n')}\n}}",
            file_name="alignment.rtf",
        )

        fmt = st.selectbox("📁 Download Format", [".fasta", ".clustal", ".nex", ".phy"])
        file_map = {
            ".fasta": lambda: "\n".join(
                r.format("fasta") for r in st.session_state["sequences"]
            ),
            ".clustal": lambda: st.session_state.get("clustal", ""),
            ".nex": lambda: st.session_state.get("nexus", ""),
            ".phy": lambda: st.session_state.get("phylip", ""),
        }
        st.download_button(
            f"⬇️ Download {fmt}", file_map[fmt](), file_name=f"alignment{fmt}"
        )

    if st.session_state["show_gc"]:
        st.subheader("📊 GC Skew")
        fig_skew = gc_skew_plot(st.session_state["sequences"])
        st.plotly_chart(fig_skew)
        buffer = BytesIO()
        fig_skew.write_image(buffer, format="png")
        st.download_button(
            "⬇️ Download GC Skew Plot (PNG)", buffer.getvalue(), file_name="gc_skew.png"
        )

    if st.session_state["show_base"]:
        st.subheader("🎯 Base Composition")
        df_base = base_composition(st.session_state["sequences"])
        st.dataframe(df_base)
        st.download_button(
            "⬇️ Download Base Composition",
            df_base.to_csv(index=False),
            file_name="base_composition.csv",
        )

    if st.session_state["show_motif"]:
        st.subheader("🧬 Motif Scan Table")
        df_motif = motif_scan_table(st.session_state["sequences"], motif_pattern)
        st.dataframe(df_motif)
        st.download_button(
            "⬇️ Download Motif Table",
            df_motif.to_csv(index=False),
            file_name="motif_scan.csv",
        )

# --- Footer ---
st.markdown(
    """
---
<p style="text-align:center;font-size:14px">
Developed by <a href="https://github.com/Behordeun">Behordeun</a> & <a href="https://github.com/bollergene">Bollergene</a><br>
📞 +2348108316393 | © Behordeun 2025
</p>
""",
    unsafe_allow_html=True,
)
