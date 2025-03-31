import streamlit as st

from style_css import style

style()

hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.title("ℹ️ About This Application")

st.markdown(
    """
## **Why Was This Application Created?**
The **Sequence Analyzer Application** was developed to provide a user-friendly, web-based tool for:
- **Biologists and researchers** needing an accessible sequence analysis platform.
- **Students and educators** who want to learn about sequence alignment and genetic analysis.
- **Data analysts in genomics** who need quick insights without requiring complex installations.

## **How Does It Work?**
1. **Upload a sequence file** (FASTA, TXT, RTF, PHYLIP, or NEXUS) or enter a GenBank accession number.
2. **Select the type of analysis** (DNA, RNA, or Phylogenetic).
3. **Choose an alignment method** (Pairwise or Multiple Sequence Alignment).
4. **View results and download aligned sequences** for further study.

## **What Makes It Unique?**
- **No installation required** - Works in your browser.
- **Multiple sequence file support** - Handles FASTA, TXT, RTF, PHYLIP, and NEXUS formats.
- **Interactive visualizations** - Provides GC content, sequence statistics, and phylogenetic trees.
- **GenBank integration** - Fetches sequences directly from NCBI.
- **Export options** - Download aligned sequences and analysis results in various formats including CSV, TXT, RTF, FASTA, CLUSTAL, NEXUS, and PHYLIP.
- **Session management** - Save and export session data for future reference.

## **Technologies Used**
- **Python** - The backbone of the application.
- **Streamlit** - For an interactive, browser-based experience.
- **Biopython** - For sequence alignment, parsing, and analysis.
- **Plotly** - To generate detailed sequence visualizations.
- **Scikit-bio** - For advanced bioinformatics computations and analyses.

---

## **🚀 Future Improvements**
- **Support for protein sequences.**
- **Integration with multiple sequence alignment tools like MUSCLE or MAFFT.**
- **Cloud storage for saving alignment history.**
- **Real-time collaboration for bioinformatics teams.**

---

## **👨‍💻 Developers**
This application was developed by:
- **Muhammad Abiodun SULAIMAN: [ResearchGate](https://www.researchgate.net/profile/Muhammad-Sulaiman-19) | [Google Scholar](https://scholar.google.com/citations?user=0EqNhMQAAAAJ&hl=en) | [ORCID](https://orcid.org/0000-0001-9161-2608) | [LinkedIn](https://linkedin.com/in/muhammad_abiodun_sulaiman) | [GitHub](https://github.com/Behordeun) | [Email](mailto:abiodun.msulaiman@gmail.com)**
- **Bolaji Fatai OYEYEMI: [ResearchGate](https://www.researchgate.net/profile/Bolaji-Oyeyemi) | [Google Scholar](https://scholar.google.com/citations?user=D0LnYT0AAAAJ&hl=en) | [ORCID](https://orcid.org/0000-0001-5564-6165) | [LinkedIn](https://zoology.lifesciences.unilorin.edu.ng/staffmember/oyeyemi-bolaji-fatia/linkedin.com/in/bolaji-f-oyeyemi-phd-46a93363) | [GitHub](https://github.com/bollergene) | [Email](mailto:bolajioyeyemi@gmail.com)**  

📅 Version: **1.1.0**  
🔗 [GitHub Repository](https://github.com/Behordeun/sequence_analyzer)
"""
)

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
📞+2348108316393
""",
    unsafe_allow_html=True,
)

st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">
Copyright | Behordeun 2025(c)
""",
    unsafe_allow_html=True,
)
