"""Contact page."""

import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles

apply_styles()
hide_streamlit_chrome()

st.title("Contact")

st.markdown(
    """
## Developers
This application was developed by **Muhammad Abiodun SULAIMAN** and
**Bolaji Fatai OYEYEMI**, who are passionate about bioinformatics and
data analysis.

## Contact Information

- **Muhammad Abiodun SULAIMAN:**
  [GitHub](https://github.com/Behordeun) |
  [ORCID](https://orcid.org/0000-0001-9161-2608) |
  [Email](mailto:abiodun.msulaiman@gmail.com)
- **Bolaji Fatai OYEYEMI:**
  [GitHub](https://github.com/bollergene) |
  [ORCID](https://orcid.org/0000-0001-5564-6165) |
  [Email](mailto:bolajioyeyemi@gmail.com)

Version: **0.2.0**
[GitHub Repository](https://github.com/Behordeun/sequence_analyzer)
"""
)

render_footer()
