import streamlit as st

st.set_page_config(page_title="Sequence Analyzer App", page_icon="🔬", layout="wide")

# Main Title
st.title("🔬 Sequence Analyzer Application")

st.markdown(
    """
## **About This Application**
The **Sequence Analyzer Application** is designed to assist researchers, students, and bioinformaticians in analyzing and visualizing DNA and RNA sequences.  

This tool provides a simple and accessible way to:
- **Align sequences** to identify similarities and variations.
- **Analyze DNA sequences** to determine GC content and nucleotide composition.
- **Analyze RNA sequences** with similar statistical insights.
- **Retrieve sequences from GenBank** using unique accession numbers.

## **How It Works**
1. **Upload a sequence file** (FASTA, TXT, or RTF) or enter GenBank accession numbers.
2. **Choose an analysis type** (DNA or RNA).
3. **Select a sequence alignment method** (Pairwise or Multiple Sequence Alignment).
4. **Download results** in various formats for further analysis or reporting.

This application is built using:
- **Streamlit** (for the interactive interface)
- **Biopython** (for sequence analysis and alignment)
- **Plotly** (for generating insightful visualizations)

---

## **🔍 Getting Started**
Use the **navigation menu on the left** to access different features.
""",
    unsafe_allow_html=True,
)

# Footer
st.markdown("---")  # Divider
st.markdown(
    """
**🔖 Sequence Analyzer Application**  
Developed with ❤️ by [Behordeun](https://github.com/Behordeun) and [Bollergene](https://github.com/bollergene).  

📅 Version: 1.0.0 | 🔗 [GitHub Repository](https://github.com/bioinformatics-project)
""",
    unsafe_allow_html=True,
)
