import pandas as pd
import plotly.express as px
from Bio import Entrez
from Bio.SeqUtils import gc_fraction

# Set Entrez email (required for NCBI API)
Entrez.email = "abiodun.msulaiman@gmail.com"


# Function to fetch sequences from GenBank
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
    return uploaded_file.read().decode("utf-8")


# Function to analyze sequences
def analyze_sequences(sequences, is_rna):
    seq_data = []
    for record in sequences:
        seq_str = str(record.seq)
        gc_content = gc_fraction(seq_str) * 100
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


# Function to visualize analysis results using Plotly
def visualize_results(df):
    fig = px.histogram(
        df,
        x="GC_Content",
        nbins=10,
        title="GC Content Distribution",
        color_discrete_sequence=["blue"],
    )
    fig.show()

    fig = px.histogram(
        df,
        x="Length",
        nbins=10,
        title="Sequence Length Distribution",
        color_discrete_sequence=["green"],
    )
    fig.show()

    fig = px.scatter(
        df,
        x="Length",
        y="GC_Content",
        text="Sequence",
        title="GC Content vs Sequence Length",
        color_discrete_sequence=["red"],
    )
    fig.update_traces(textposition="top center")
    fig.show()
