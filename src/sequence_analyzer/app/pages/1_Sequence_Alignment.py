"""Sequence Alignment page: pairwise/MSA, GC skew, base composition, motif scanning."""

import datetime
import os
from io import StringIO

import plotly.express as px
import streamlit as st
from Bio import SeqIO

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.alignment import align_msa, align_pairwise
from sequence_analyzer.core.analysis import compute_base_composition, compute_gc_skew
from sequence_analyzer.core.genbank import fetch_sequence
from sequence_analyzer.core.motifs import scan_motifs
from sequence_analyzer.core.validation import clean_sequence
from sequence_analyzer.io.parsers import convert_to_fasta

st.set_page_config(layout="wide")
apply_styles()
hide_streamlit_chrome()

st.title("🔗 Sequence Alignment, GC Skew, Motif Scan & Base Composition")


# --- Helpers ---


def _get_seq_type(choice: str) -> str:
    """Normalize the UI sequence type choice to the internal format."""
    return "auto" if choice == "Auto" else choice.lower()


def _get_sequences() -> list:
    return st.session_state.get("sequences", [])


def _set_alignment_result(html: str, metadata: dict, format_exports: dict | None = None) -> None:
    st.session_state["aligned"] = html
    st.session_state["metadata"] = metadata
    if format_exports is not None:
        st.session_state["format_exports"] = format_exports
    else:
        st.session_state.pop("format_exports", None)


def _highlight_alignment(seq1: str, seq2: str) -> str:
    """Generate colored HTML showing matches/mismatches between two aligned sequences."""
    html1 = "".join(
        f"<span class='{'match' if a == b else 'mismatch'}'>{a}</span>"
        for a, b in zip(seq1, seq2, strict=False)
    )
    html2 = "".join(
        f"<span class='{'match' if a == b else 'mismatch'}'>{b}</span>"
        for a, b in zip(seq1, seq2, strict=False)
    )
    return html1 + "<br>" + html2


def _build_pairwise_html(result) -> str:
    html_blocks = []
    for pair in result.pairs:
        html_blocks.append(
            f"<b>{pair.seq_a_id} vs {pair.seq_b_id}</b><br>"
            f"<b>Score:</b> {pair.score:.2f}, "
            f"<b>Identity:</b> {pair.identity:.1f}%<br>"
            f"{_highlight_alignment(pair.aligned_a, pair.aligned_b)}<br><br>"
        )
    summary = (
        f"<b>Average Score:</b> {result.avg_score:.2f} | "
        f"<b>Average Identity:</b> {result.avg_identity:.1f}%<br><hr>"
    )
    return summary + "".join(html_blocks)


def _build_msa_html(result) -> str:
    return f"<b>MSA Identity:</b> {result.avg_identity:.2f}%<br><hr>"


def _run_alignment(seqs, seq_type_choice, matrix, mode, alignment_method):
    """Orchestrate alignment and return (html, metadata, format_exports or None)."""
    seq_type = _get_seq_type(seq_type_choice)
    metadata = {
        "matrix": matrix,
        "mode": mode,
        "type": seq_type_choice,
        "timestamp": str(datetime.datetime.now()),
    }

    if alignment_method == "Pairwise Alignment":
        result = align_pairwise(seqs, matrix=matrix, mode=mode)
        html = _build_pairwise_html(result)
        return html, metadata, None
    else:
        result = align_msa(seqs, seq_type=seq_type)
        html = _build_msa_html(result)
        return html, metadata, result.format_exports


def _load_from_upload(files, seq_type_choice):
    """Parse and validate sequences from uploaded files."""
    seq_type = _get_seq_type(seq_type_choice)
    valid = []
    for f in files:
        content = convert_to_fasta(f.read())
        for rec in SeqIO.parse(StringIO(content), "fasta"):
            try:
                valid.append(clean_sequence(rec, seq_type=seq_type))
            except ValueError as e:
                st.warning(f"Skipped {rec.id}: {e}")
    return valid


# --- Session state initialization ---

for key in ["sequences", "aligned", "metadata"]:
    st.session_state.setdefault(key, [])
for key in ["show_alignment", "show_gc", "show_base", "show_motif"]:
    st.session_state.setdefault(key, False)

# --- Input options ---

option = st.selectbox(
    "Choose Input Method",
    ["Upload Sequence File(s)", "Enter Accession Numbers"],
)
seq_type_choice = st.radio("Sequence Type", ["Auto", "DNA", "RNA", "Protein"])
alignment_method = st.radio(
    "Alignment Method", ["Pairwise Alignment", "Multiple Sequence Alignment"]
)

# --- File input ---

if option == "Upload Sequence File(s)":
    uploaded = st.file_uploader(
        "Upload sequence file(s)",
        type=["fasta", "txt", "rtf", "phy", "nex"],
        accept_multiple_files=True,
    )
    if uploaded:
        valid = _load_from_upload(uploaded, seq_type_choice)
        st.session_state["sequences"] = valid
        st.success(f"{len(valid)} valid sequences loaded.")

elif option == "Enter Accession Numbers":
    acc_input = st.text_area("Enter Accession Numbers (one per line)", height=120)
    genbank_email = st.text_input(
        "Email for NCBI API",
        value=os.environ.get("GENBANK_EMAIL", ""),
        help="Required by NCBI usage policy. Set GENBANK_EMAIL env var to avoid entering each time.",
    )
    if st.button("Retrieve Sequences"):
        if not genbank_email.strip():
            st.error("Email is required for GenBank API access.")
        else:
            fetched = []
            for acc in acc_input.strip().splitlines():
                try:
                    record = fetch_sequence(acc.strip(), email=genbank_email.strip())
                    fetched.append(record)
                except (ValueError, ConnectionError) as e:
                    st.warning(f"Could not fetch {acc}: {e}")
            st.session_state["sequences"] = fetched
            st.success(f"Retrieved {len(fetched)} sequence(s) from GenBank.")

# --- Alignment ---

sequences = _get_sequences()
if sequences:
    matrix = st.selectbox(
        "Scoring Matrix",
        ["BLOSUM62", "BLOSUM80"] if seq_type_choice == "Protein" else ["NUC.4.4", "IDENTITY"],
    )
    mode = st.selectbox("Alignment Mode", ["global", "local"])
    motif_pattern = st.text_input("Motif Pattern (Regex)", "ATG")

    st.session_state["show_alignment"] = st.checkbox("Show Alignment Result", True)
    st.session_state["show_gc"] = st.checkbox("Show GC Skew")
    st.session_state["show_base"] = st.checkbox("Base Composition Summary")
    st.session_state["show_motif"] = st.checkbox("Motif Scan Table")

    if st.button("Align Sequences"):
        html, metadata, format_exports = _run_alignment(
            sequences, seq_type_choice, matrix, mode, alignment_method
        )
        _set_alignment_result(html, metadata, format_exports)
        st.success("Alignment complete!")

# --- Output ---

if st.session_state.get("aligned"):
    if st.session_state["show_alignment"]:
        st.subheader("Alignment Output")
        st.markdown(st.session_state["aligned"], unsafe_allow_html=True)

        exports = st.session_state.get("format_exports", {})
        if exports:
            fmt = st.selectbox("Download Format", [".fasta", ".clustal", ".nex", ".phy"])
            fmt_map = {".fasta": "fasta", ".clustal": "clustal", ".nex": "nexus", ".phy": "phylip"}
            content = exports.get(fmt_map.get(fmt, "fasta"), "")
            if content:
                st.download_button(f"Download {fmt}", content, file_name=f"alignment{fmt}")
        else:
            st.info("Download formats are available for multiple sequence alignments.")

    if st.session_state["show_gc"]:
        st.subheader("GC Skew")
        df_skew = compute_gc_skew(sequences, window_size=100)
        if not df_skew.empty:
            fig = px.line(df_skew, x="Position", y="GC_Skew", color="ID", title="GC Skew")
            st.plotly_chart(fig, use_container_width=True)

    if st.session_state["show_base"]:
        st.subheader("Base Composition")
        df_base = compute_base_composition(sequences)
        st.dataframe(df_base)
        st.download_button(
            "Download Base Composition",
            df_base.to_csv(index=False),
            file_name="base_composition.csv",
        )

    if st.session_state["show_motif"]:
        st.subheader("Motif Scan Table")
        df_motif = scan_motifs(sequences, motif_pattern)
        st.dataframe(df_motif)
        st.download_button(
            "Download Motif Table",
            df_motif.to_csv(index=False),
            file_name="motif_scan.csv",
        )

render_footer()
