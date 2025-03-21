import streamlit as st

st.title("ℹ️ About This Application")

st.markdown(
    """
## **Why Was This Application Created?**
The **Sequence Analyzer Application** was developed to provide a user-friendly, web-based tool for:
- **Biologists and researchers** needing an accessible sequence analysis platform.
- **Students and educators** who want to learn about sequence alignment and genetic analysis.
- **Data analysts in genomics** who need quick insights without requiring complex installations.

## **How Does It Work?**
1. **Upload a sequence file** (FASTA, TXT, or RTF) or enter a GenBank accession number.
2. **Select the type of analysis** (DNA or RNA).
3. **Choose an alignment method** (Pairwise or Multiple Sequence Alignment).
4. **View results and download aligned sequences** for further study.

## **What Makes It Unique?**
- **No installation required** – Works in your browser.
- **Multiple sequence file support** – Handles FASTA, TXT, and RTF formats.
- **Interactive visualizations** – Provides GC content and sequence statistics.
- **GenBank integration** – Fetches sequences directly from NCBI.
- **Export options** – Download aligned sequences in different formats.

## **Technologies Used**
- **Python** – The backbone of the application.
- **Streamlit** – For an interactive, browser-based experience.
- **Biopython** – For sequence alignment, parsing, and analysis.
- **Plotly** – To generate detailed sequence visualizations.

---

## **🚀 Future Improvements**
- **Support for protein sequences.**
- **Integration with multiple sequence alignment tools like MUSCLE or MAFFT.**
- **Cloud storage for saving alignment history.**
- **Real-time collaboration for bioinformatics teams.**

---

## **👨‍💻 Developers**
This application was developed by:
- **Muhammad Abiodun SULAIMAN: [ResearchGate](https://www.researchgate.net/profile/Muhammad-Sulaiman-19) | [Google Scholar](https://scholar.google.com/citations?user=0EqNhMQAAAAJ&hl=en) | [ORCID](https://orcid.org/0000-0001-9161-2608) | [LinkedIn](https://linkedin.com/in/muhammad-abiodun-sulaiman) | [GitHub](https://github.com/Behordeun) | [Email](mailto:abiodun.msulaiman@gmail.dom)**
- **Bolaji Fatai OYEYEMI: [ResearchGate](https://www.researchgate.net/profile/Bolaji-Oyeyemi) | [Google Scholar](https://scholar.google.com/citations?user=D0LnYT0AAAAJ&hl=en) | [ORCID](https://orcid.org/0000-0001-5564-6165) | [LinkedIn](https://zoology.lifesciences.unilorin.edu.ng/staffmember/oyeyemi-bolaji-fatia/linkedin.com/in/bolaji-f-oyeyemi-phd-46a93363) | [GitHub](https://github.com/bollergene) | [Email](mailto:bolajioyeyemi@gmail.com)**  

📅 Version: **1.0.0**  
🔗 [GitHub Repository](https://github.com/Behordeun/sequence_analyzer)
"""
)


# Footer
st.markdown("---")  # Divider
st.markdown(
    """
Developed by [Behordeun](https://github.com/Behordeun) and [Bollergene](https://github.com/bollergene).  

📅 Version: 1.0.0 | 🔗 [GitHub Repository](https://github.com/bioinformatics-project)

Copyright | Behordeun 2025(c)
""",
    unsafe_allow_html=True,
)
