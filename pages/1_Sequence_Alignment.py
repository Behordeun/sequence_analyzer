from io import StringIO

import streamlit as st
from Bio import SeqIO

from utils import align_sequences, convert_to_fasta, fetch_sequence

# Page Title
st.title("🔗 Sequence Alignment")


# User Selection: Upload or Fetch from GenBank
option = st.selectbox(
    "Choose Input Method", ["Upload Sequence File", "Enter Assertion Numbers"]
)


sequences = []


# Option 1: Upload Sequence File
if option == "Upload Sequence File":
    uploaded_file = st.file_uploader(
        "Upload Sequence File (FASTA format)", type=["fasta", "txt", "rtf"]
    )
    if uploaded_file:
        fasta_content = convert_to_fasta(uploaded_file)
        fasta_records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
        sequences.extend(fasta_records)
        st.success(f"{len(fasta_records)} sequences successfully uploaded!")


# Option 2: Enter Assertion Numbers
elif option == "Enter Assertion Numbers":
    assertion_numbers = st.text_area(
        "Enter Assertion Numbers (one per line)"
    ).splitlines()
    if assertion_numbers:
        st.info("Fetching sequences from GenBank...")
        genbank_sequences = [fetch_sequence(acc) for acc in assertion_numbers]
        sequences.extend(
            [SeqIO.read(StringIO(seq), "fasta") for seq in genbank_sequences if seq]
        )
        st.success(
            f"{len(genbank_sequences)} sequences successfully retrieved from GenBank!"
        )


# Perform Alignment if Sequences Are Available
if sequences:
    fasta_sequences = [f">{rec.id}\n{rec.seq}" for rec in sequences]
    aligned_file_path = align_sequences(fasta_sequences)
    st.success("✅ Sequence Alignment Completed!")

    # Let the user choose the file format for download
    file_format = st.selectbox("Select Download Format", [".fasta", ".rtf", ".txt"])

    # Determine file extension and format
    file_extension = file_format.lower()
    download_file_name = f"aligned_sequences{file_extension}"

    with open(aligned_file_path, "rb") as file:
        st.download_button(
            f"Download Alignment ({file_format.upper()})",
            data=file,
            file_name=download_file_name,
        )
