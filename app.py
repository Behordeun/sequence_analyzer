import streamlit as st

st.set_page_config(
    page_title="Open-Source Sequencing App", page_icon="🔬", layout="wide"
)

# Main Title
st.title("🔬 Open-Source Sequencing Application")
st.markdown(
    """
Welcome to the **Open-Source Sequencing Application**! This tool allows users to:
- **Align sequences** using Clustal Omega
- **Analyze DNA sequences** (GC content, nucleotide composition)
- **Analyze RNA sequences** (GC content, nucleotide composition)
- **Fetch sequences from GenBank** using accession numbers

---

### 🔍 **Get Started**:
Use the **navigation menu on the left** to explore different features.
"""
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
