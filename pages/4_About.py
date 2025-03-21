import streamlit as st

st.title("ℹ️ About This App")
st.markdown(
    """
This Open-Source Sequencing Application was developed using:
- **Python**
- **Streamlit**
- **Biopython**
- **Clustal Omega**
- **Plotly**

### 🛠 Features:
- **Sequence Alignment**
- **DNA Analysis**
- **RNA Analysis**
- **GenBank Sequence Retrieval**
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
