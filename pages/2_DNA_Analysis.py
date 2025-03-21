from io import StringIO
import streamlit as st
from Bio import SeqIO

from style_css import style
from utils import analyze_sequences, visualize_results

# Apply custom styles
style()

# Hide Streamlit menu and footer
st.markdown(
    """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 DNA Sequence Analysis")

# **Reset session state when the user navigates to a different page**
if "last_page" not in st.session_state or st.session_state["last_page"] != "DNA_Analysis":
    st.session_state.clear()
    st.session_state["last_page"] = "DNA_Analysis"

# Initialize session state for sequence data and results
if "sequence_df" not in st.session_state:
    st.session_state["sequence_df"] = None

# File Upload Section
uploaded_file = st.file_uploader(
    "Upload DNA Sequence File (FASTA)", type=["fasta", "txt", "rtf"]
)

if uploaded_file:
    # Convert file to text mode using StringIO
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))

    # Read sequences (handling comments with 'fasta-pearson' format)
    sequences = list(SeqIO.parse(stringio, "fasta-pearson"))

    # Display detected sequences
    st.subheader("📄 Identified Sequences:")
    for seq in sequences:
        st.write(f"✅ {seq.id} - {len(seq.seq)} bp")

    # Analyze Button
    if st.button("🔬 Analyze"):
        st.session_state["sequence_df"] = analyze_sequences(sequences, is_rna=False)
        st.success("✅ Analysis Completed!")

# Display results if analysis has been performed
if st.session_state["sequence_df"] is not None:
    sequence_df = st.session_state["sequence_df"]

    # Show DataFrame
    st.dataframe(sequence_df)

    # User selects the download format
    file_format = st.selectbox("📂 Select Download Format", [".csv", ".txt", ".rtf"])

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
        f"📥 Download DNA Analysis ({file_format.upper()})",
        data=file_data,
        file_name=download_file_name,
    )

    # 📊 Visualize results
    visualize_results(sequence_df)

# Footer
st.markdown("---")  # Divider
st.markdown(
    """
    <p style="color: black; text-align: center; font-size: 15px;">
        Developed by 
        <a href="https://github.com/Behordeun" target="_blank" style="color: blue; text-decoration: none;">Behordeun</a> 
        and 
        <a href="https://github.com/bollergene" target="_blank" style="color: blue; text-decoration: none;">Bollergene</a>.
    </p>
""",
    unsafe_allow_html=True,
)

st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">
📞 +2348108316393
""",
    unsafe_allow_html=True,
)

st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">
Copyright | Behordeun 2025(c)
""",
    unsafe_allow_html=True,
)