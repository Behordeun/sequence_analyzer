import streamlit as st

from style_css import style

st.set_page_config(page_title="Sequence Analyzer App", page_icon="🔬", layout="wide")
style()

hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

# Main Title
st.title("🔬 Sequence Analyzer Application")

st.markdown(
    """
## **About This Application**
The **Sequence Analyzer Application** is built to help researchers, students, and bioinformaticians perform DNA, RNA, and Phylogenetic sequence analysis with ease, all from a web interface.

---

## ✅ **What You Can Do with This Tool**
- 🔡 **Align Sequences:** Upload one or multiple sequences and choose between **Pairwise or Multiple Sequence Alignment**.
- 🔬 **Analyze DNA or RNA:** See nucleotide composition and GC content, with charts.
- 🔗 **Retrieve Sequences from GenBank:** Just enter accession numbers and fetch data directly.
- 🌿 **Build Phylogenetic Trees:** Generate evolutionary trees using **UPGMA or Neighbor-Joining** and download them in Newick format.

---

## 🧪 **How It Works**
1. **Upload** sequences in **.fasta**, **.txt**, or **.rtf** format, or input GenBank accession numbers.
2. **Choose** an analysis type (DNA or RNA).
3. **Select** alignment and tree-building methods.
4. **Visualize and Download** your results as CSV, RTF, TXT, or Newick files.

---

## 🧭 **Getting Started**
Use the **navigation menu on the left** to:
1. Perform Sequence Alignment
2. Access DNA / RNA Analysis
3. Build Phylogenetic Trees
4. Learn more from the **"About"** section


""",
    unsafe_allow_html=True,
)

# Footer
st.markdown("---")
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
