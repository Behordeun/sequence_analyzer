from io import StringIO

import streamlit as st
from Bio import SeqIO, pairwise2
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import convert_to_fasta, fetch_sequence

# Page Title
st.title("🔗 Sequence Alignment")

# User Selection: Upload File(s) or Fetch from GenBank
option = st.selectbox(
    "Choose Input Method", ["Upload Sequence File(s)", "Enter Assertion Numbers"]
)

# Select sequence type (DNA or RNA)
sequence_type = st.radio("🔬 Select Sequence Type", ["DNA", "RNA"])

# Select alignment method (Pairwise or Multiple Sequence Alignment)
alignment_method = st.radio(
    "🧬 Select Alignment Method", ["Pairwise Alignment", "Multiple Sequence Alignment"]
)

# Initialize session state for sequences
if "sequences" not in st.session_state:
    st.session_state["sequences"] = []

# Option 1: Upload Multiple Sequence Files
if option == "Upload Sequence File(s)":
    uploaded_files = st.file_uploader(
        "Upload Sequence Files (FASTA format)",
        type=["fasta", "txt", "rtf"],
        accept_multiple_files=True,
    )

    if uploaded_files:
        total_sequences = 0
        uploaded_sequences = []
        for uploaded_file in uploaded_files:
            fasta_content = convert_to_fasta(uploaded_file)
            fasta_records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
            uploaded_sequences.extend(fasta_records)
            total_sequences += len(fasta_records)

        st.session_state["sequences"] = uploaded_sequences  # Store in session state
        st.success(f"{total_sequences} sequences successfully uploaded!")

        # Display detected sequences
        st.subheader("📄 Identified Sequences:")
        for seq in st.session_state["sequences"]:
            st.write(f"✅ {seq.id} - {len(seq.seq)} bp")

# Option 2: Enter Assertion Numbers
elif option == "Enter Assertion Numbers":
    assertion_numbers = st.text_area(
        "Enter Assertion Numbers (one per line)"
    ).splitlines()

    if assertion_numbers and st.button("🔍 Fetch Sequences from GenBank"):
        st.info("Fetching sequences from GenBank...")
        genbank_sequences = [fetch_sequence(acc) for acc in assertion_numbers]
        fetched_sequences = [
            SeqIO.read(StringIO(seq), "fasta") for seq in genbank_sequences if seq
        ]
        st.session_state["sequences"] = fetched_sequences  # Store in session state
        st.success(
            f"{len(fetched_sequences)} sequences successfully retrieved from GenBank!"
        )

        # Display detected sequences
        st.subheader("📄 Identified Sequences:")
        for seq in st.session_state["sequences"]:
            st.write(f"✅ {seq.id} - {len(seq.seq)} bp")


# Function for Pairwise Alignment using Biopython
def align_sequences_pairwise(sequences):
    if len(sequences) < 2:
        return "Error: Pairwise alignment requires at least two sequences."

    aligned_results = []
    for i in range(len(sequences) - 1):
        seq1 = sequences[i].seq
        seq2 = sequences[i + 1].seq
        alignments = pairwise2.align.globalxx(seq1, seq2)
        aligned_seq1, aligned_seq2, _, _, _ = alignments[0]
        aligned_results.append(f">{sequences[i].id}\n{aligned_seq1}")
        aligned_results.append(f">{sequences[i + 1].id}\n{aligned_seq2}")

    return "\n".join(aligned_results)


# Function for Multiple Sequence Alignment (MSA) using Biopython
def align_sequences_msa(sequences):
    if len(sequences) < 2:
        return "Error: MSA requires at least two sequences."

    # Find the longest sequence length
    max_length = max(len(seq.seq) for seq in sequences)

    # Pad shorter sequences with gaps (-)
    padded_sequences = [
        SeqRecord(Seq(str(seq.seq).ljust(max_length, "-")), id=seq.id, description="")
        for seq in sequences
    ]

    # Perform Multiple Sequence Alignment
    alignment = MultipleSeqAlignment(padded_sequences)

    # Generate consensus sequence
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()

    # Save the aligned sequences
    aligned_fasta = "\n".join([f">{seq.id}\n{seq.seq}" for seq in padded_sequences])
    aligned_fasta += f"\n>Consensus\n{consensus}"

    return aligned_fasta


# Align Sequence Button (Prevents automatic alignment)
if st.session_state["sequences"]:
    if st.button("🔄 Align Sequence"):
        if alignment_method == "Pairwise Alignment":
            st.session_state["aligned_sequences"] = align_sequences_pairwise(
                st.session_state["sequences"]
            )
        else:
            st.session_state["aligned_sequences"] = align_sequences_msa(
                st.session_state["sequences"]
            )

        st.success("✅ Sequence Alignment Completed!")

# Ensure the download button persists after alignment
if "aligned_sequences" in st.session_state:
    # Let the user choose the file format for download
    file_format = st.selectbox("📂 Select Download Format", [".fasta", ".rtf", ".txt"])

    # Generate a dynamic file name based on the number of sequences
    file_extension = file_format.lower()
    num_sequences = len(st.session_state["sequences"])
    download_file_name = f"aligned_sequences_{num_sequences}{file_extension}"

    # Convert alignment to the selected format
    if file_format == ".fasta":
        file_data = st.session_state["aligned_sequences"]
    elif file_format == ".txt":
        file_data = st.session_state["aligned_sequences"].replace(
            ">", "\n>"
        )  # Format for readability
    elif file_format == ".rtf":
        file_data = (
            "{\\rtf1\\ansi\n"
            + st.session_state["aligned_sequences"].replace(">", "\n>")
            + "\n}"
        )

    # Provide download button
    st.download_button(
        "📥 Download Alignment",
        data=file_data,
        file_name=download_file_name,
    )


# Footer
st.markdown("---")  # Divider
st.markdown(
    """
    **🔖 Open-Source Sequencing Application**  
    Developed with ❤️ using **Streamlit, Biopython, and Plotly by [Behordeun](https://github.com/Behordeun) and [Bollergene](https://github.com/bollergene)**  
    📅 Version: 1.0.0 | 🔗 [GitHub Repository](https://github.com/bioinformatics-project)
    
    """,
    unsafe_allow_html=True,
)
