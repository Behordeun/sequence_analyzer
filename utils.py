import subprocess
import tempfile

import pandas as pd
import plotly.express as px
from Bio import Entrez
from Bio.SeqUtils import gc_fraction

Entrez.email = "abc@example.com"


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


def convert_to_fasta(uploaded_file):
    return uploaded_file.read().decode("utf-8")


def align_sequences(sequences):
    # Save sequences to a temporary FASTA file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as fasta_file:
        fasta_file.write("\n".join(sequences).encode("utf-8"))
        fasta_file_path = fasta_file.name

    aligned_file_path = fasta_file_path.replace(".fasta", "_aligned.fasta")

    try:
        subprocess.run(
            ["clustalo", "-i", fasta_file_path, "-o", aligned_file_path, "--auto"],
            check=True,
        )
    except FileNotFoundError:
        return "Error: Clustal Omega (clustalo) not found. Please ensure it is installed and in your system PATH."

    return aligned_file_path


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
