"""Shared UI components used across all pages."""

import streamlit as st


def hide_streamlit_chrome() -> None:
    """Hide the default Streamlit main menu and footer."""
    st.markdown(
        "<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>",
        unsafe_allow_html=True,
    )


def render_footer() -> None:
    """Render the shared application footer with developer credits."""
    st.markdown("---")
    st.markdown(
        """
        <p style="color: white; text-align: center; font-size: 15px;">
            Developed by
            <a href="https://github.com/Behordeun" target="_blank"
               style="color: blue; text-decoration: none;">Behordeun</a>
            and
            <a href="https://github.com/bollergene" target="_blank"
               style="color: blue; text-decoration: none;">Bollergene</a>.
        </p>
        <p style="color: white; text-align: center; font-size: 14px;">
            Copyright | Behordeun 2025(c)
        </p>
        """,
        unsafe_allow_html=True,
    )
