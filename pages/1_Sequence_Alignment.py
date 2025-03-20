import time
from io import StringIO

import requests
import streamlit as st
from Bio import SeqIO

from utils import convert_to_fasta, fetch_sequence

# Page Title
st.title("🔗 Sequence Alignment")

# User Selection: Upload or Fetch from GenBank
option = st.selectbox(
    "Choose Input Method", ["Upload Sequence File", "Enter Assertion Numbers"]
)

# User email input (EBI Clustal API requires email)
user_email = st.text_input("Enter your email (required for alignment)", "")

# Select sequence type (DNA or RNA)
sequence_type = st.radio("Select Sequence Type", ["DNA", "RNA"])
ebi_sequence_type = "dna" if sequence_type == "DNA" else "rna"

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


# Function to Perform MSA Using EBI Clustal Omega API
def align_sequences_ebi_clustal(sequences, user_email, seq_type):
    fasta_data = "\n".join([f">{seq.id}\n{seq.seq}" for seq in sequences])

    url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/"
    payload = {
        "email": user_email,
        "sequence": fasta_data,
        "stype": seq_type,  # Either "dna" or "rna"
        "format": "fasta",
    }

    response = requests.post(url, data=payload)

    if response.status_code == 200:
        job_id = response.text
        st.info("Waiting for alignment to complete...")

        # Check job status
        status_url = (
            f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
        )
        result_url = (
            f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/fasta"
        )

        while True:
            status = requests.get(status_url).text.strip()
            if status == "FINISHED":
                break
            time.sleep(5)  # Wait before checking again

        # Fetch aligned sequences
        aligned_sequences = requests.get(result_url).text
        return aligned_sequences

    else:
        return f"Error: {response.text}"


# Perform Alignment if Sequences Are Available
if sequences and user_email:
    aligned_sequences = align_sequences_ebi_clustal(
        sequences, user_email, ebi_sequence_type
    )
    st.success("✅ Sequence Alignment Completed!")

    # Let the user choose the file format for download
    file_format = st.selectbox("Select Download Format", [".fasta", ".rtf", ".txt"])

    # Determine file extension and format
    file_extension = file_format.lower()
    download_file_name = f"aligned_sequences{file_extension}"

    # Convert alignment to the selected format
    if file_format == ".fasta":
        file_data = aligned_sequences
    elif file_format == ".txt":
        file_data = aligned_sequences.replace(">", "\n>")  # Format for readability
    elif file_format == ".rtf":
        file_data = "{\\rtf1\\ansi\n" + aligned_sequences.replace(">", "\n>") + "\n}"

    # Provide download button
    st.download_button(
        f"Download Alignment ({file_format.upper()})",
        data=file_data,
        file_name=download_file_name,
    )
