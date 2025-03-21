from io import StringIO

import streamlit as st
from Bio import SeqIO

from utils import analyze_sequences, visualize_results

st.title("🧬 DNA Sequence Analysis")

# File Upload Section
uploaded_file = st.file_uploader(
    "Upload DNA Sequence File (FASTA)", type=["fasta", "txt", "rtf"]
)

if uploaded_file:
    # Convert file to text mode using StringIO
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))

    # Read and analyze sequences
    sequences = list(SeqIO.parse(stringio, "fasta-pearson"))
    sequence_df = analyze_sequences(sequences, is_rna=False)

    # Display results
    st.dataframe(sequence_df)

    # User selects the download format
    file_format = st.selectbox("Select Download Format", [".csv", ".txt", ".rtf"])

    # Generate a dynamic file name
    download_file_name = f"dna_analysis{file_format}"

    # Convert DataFrame to selected format
    if file_format == ".csv":
        file_data = sequence_df.to_csv(index=False)
    elif file_format == ".txt":
        file_data = sequence_df.to_string(index=False)
    elif file_format == ".rtf":
        file_data = "{\\rtf1\\ansi\n" + sequence_df.to_string(index=False) + "\n}"

    # Provide download button
    st.download_button(
        f"Download DNA Analysis ({file_format.upper()})",
        data=file_data,
        file_name=download_file_name,
    )

    # Visualize results
    visualize_results(sequence_df)
