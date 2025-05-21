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


st.title("📞 Contact")

st.markdown(
    """

    ## **👨‍💻 Developers**
    This application was developed by: **[Muhammad Abiodun SULAIMAN](https://linkedin.com/in/muhammad_abiodun_sulaiman)** and **[Bolaji Fatai OYEYEMI](https://linkedin.com/in/bolaji-f-oyeyemi-phd-46a93363)**, who are passionate about bioinformatics and data analysis. They aim to provide a user-friendly platform for researchers and students to perform sequence analysis without the need for complex installations or configurations.
    ## **📧 Contact Information**
    
    - **Muhammad Abiodun SULAIMAN: [ResearchGate](https://www.researchgate.net/profile/Muhammad-Sulaiman-19) | [Google Scholar](https://scholar.google.com/citations?user=0EqNhMQAAAAJ&hl=en) | [ORCID](https://orcid.org/0000-0001-9161-2608) | [GitHub](https://github.com/Behordeun) | [Email](mailto:abiodun.msulaiman@gmail.com)**
    - **Bolaji Fatai OYEYEMI: [ResearchGate](https://www.researchgate.net/profile/Bolaji-Oyeyemi) | [Google Scholar](https://scholar.google.com/citations?user=D0LnYT0AAAAJ&hl=en) | [ORCID](https://orcid.org/0000-0001-5564-6165) | [GitHub](https://github.com/bollergene) | [Email](mailto:bolajioyeyemi@gmail.com)**  

    📅 Version: **1.1.2**  
    🔗 [GitHub Repository](https://github.com/Behordeun/sequence_analyzer)

"""
)


# Footer
st.markdown("---")  # Divider
st.markdown(
    """
    <p style="color: white; text-align: center; font-size: 15px;">
        Developed by 
        <a href="https://github.com/Behordeun" target="_blank" style="color: blue; text-decoration: none;">Behordeun</a> 
        and 
        <a href="https://github.com/bollergene" target="_blank" style="color: blue; text-decoration: none;">Bollergene</a>.
    </p>
""",
    unsafe_allow_html=True,
)

st.markdown(
    """<p style="color:white; text-align:center;font-size:15px;">
📞+2348108316393
""",
    unsafe_allow_html=True,
)

st.markdown(
    """<p style="color:whilte; text-align:center;font-size:15px;">
Copyright | Behordeun 2025(c)
""",
    unsafe_allow_html=True,
)
