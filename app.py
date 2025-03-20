import tempfile
from io import StringIO

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqUtils import gc_fraction 

# Set Entrez email (Required by NCBI)
Entrez.email = "abiodun.msulaiman@gmail.com"


# Function to fetch sequence from GenBank
def fetch_sequence(accession):
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="fasta", retmode="text"
        )
        record = handle.read()
        handle.close()
        return record
    except Exception as e:
        return f"Error fetching {accession}: {str(e)}"


# Function to convert uploaded file to FASTA format
def convert_to_fasta(uploaded_file):
    try:
        content = uploaded_file.read().decode("utf-8")
        return content
    except Exception as e:
        return f"Error processing file: {str(e)}"


# Function to perform sequence alignment using Clustal Omega
def align_sequences(sequences):
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as fasta_file:
        fasta_file.write("\n".join(sequences).encode("utf-8"))
        fasta_file_path = fasta_file.name

    aligned_file_path = fasta_file_path.replace(".fasta", "_aligned.fasta")

    clustalomega_cline = ClustalOmegaCommandline(
        infile=fasta_file_path, outfile=aligned_file_path, verbose=True, auto=True
    )
    clustalomega_cline()

    return aligned_file_path


# Function to analyze sequences
def analyze_sequences(sequences, is_rna):
    seq_data = []
    for record in sequences:
        seq_str = str(record.seq)
        gc_content = gc_fraction(seq_str) * 100  # Convert fraction to percentage
        nucleotide_counts = {
            nuc: seq_str.count(nuc) for nuc in ("A", "U" if is_rna else "T", "G", "C")
        }

        seq_data.append(
            {
                "Sequence": record.id,
                "Length": len(seq_str),
                "GC_Content": gc_content,
                **nucleotide_counts,
            }
        )

    return pd.DataFrame(seq_data)


# Function to visualize analysis using Plotly
def visualize_results(df):
    st.subheader("GC Content Distribution")
    fig = px.histogram(
        df,
        x="GC_Content",
        nbins=10,
        title="GC Content Distribution",
        color_discrete_sequence=["blue"],
    )
    st.plotly_chart(fig)

    st.subheader("Sequence Length Distribution")
    fig = px.histogram(
        df,
        x="Length",
        nbins=10,
        title="Sequence Length Distribution",
        color_discrete_sequence=["green"],
    )
    st.plotly_chart(fig)

    st.subheader("GC Content vs Sequence Length")
    fig = px.scatter(
        df,
        x="Length",
        y="GC_Content",
        text="Sequence",
        title="GC Content vs Sequence Length",
        color_discrete_sequence=["red"],
    )
    fig.update_traces(textposition="top center")
    st.plotly_chart(fig)


# Streamlit UI
st.title("Open-Source Sequencing Application")
st.write(
    "Upload a sequence file or specify assertion numbers to fetch sequences from GenBank."
)

# Sidebar
st.sidebar.header("Options")
analysis_type = st.sidebar.radio(
    "Choose Analysis Type", ["DNA Analysis", "RNA Analysis"]
)
is_rna = True if analysis_type == "RNA Analysis" else False

# File Upload Section
uploaded_file = st.file_uploader(
    "Upload Sequence File (FASTA format)", type=["fasta", "txt", "rtf"]
)

# Assertion Number Input Section
assertion_numbers = st.text_area("Enter Assertion Numbers (one per line)").splitlines()

# Process File Upload
sequences = []
if uploaded_file:
    st.success("File uploaded successfully!")
    fasta_content = convert_to_fasta(uploaded_file)
    fasta_records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
    sequences.extend(fasta_records)

# Fetch Sequences from GenBank
if assertion_numbers:
    st.info("Fetching sequences from GenBank...")
    genbank_sequences = []
    for acc in assertion_numbers:
        seq_data = fetch_sequence(acc)
        if seq_data and "Error" not in seq_data:
            fasta_io = StringIO(seq_data)
            genbank_sequences.extend(list(SeqIO.parse(fasta_io, "fasta")))
        else:
            st.warning(f"Could not fetch: {acc}")

    sequences.extend(genbank_sequences)

# Proceed with Analysis if Sequences are Available
if sequences:
    st.success(f"{len(sequences)} sequences retrieved!")

    # Display sequence information
    for record in sequences:
        st.write(f"**{record.id}** - Length: {len(record.seq)} bp")

    # Perform Alignment
    st.subheader("Aligning Sequences...")
    fasta_sequences = [f">{rec.id}\n{rec.seq}" for rec in sequences]
    aligned_file_path = align_sequences(fasta_sequences)
    st.success("Sequence Alignment Completed!")

    # Load aligned sequences
    aligned_sequences = list(SeqIO.parse(aligned_file_path, "fasta"))

    # Sequence Analysis
    st.subheader("Sequence Analysis")
    sequence_df = analyze_sequences(aligned_sequences, is_rna)
    st.dataframe(sequence_df)

    # Save results
    sequence_df.to_csv("sequence_analysis_results.csv", index=False)
    st.download_button(
        "Download Analysis Results (CSV)",
        data=sequence_df.to_csv(index=False),
        file_name="sequence_analysis_results.csv",
    )

    # Visualize results
    visualize_results(sequence_df)

st.write("### Developed using Python, Streamlit & Plotly 🚀")
