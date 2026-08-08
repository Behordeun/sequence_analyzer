"""Application CSS styles injected via st.markdown."""

import streamlit as st

CUSTOM_CSS = """
<style>
body {
    background-color: #f4f4f9;
    color: #333333;
    font-family: 'Arial', sans-serif;
}
.main-content {
    margin-top: 50px;
    padding: 20px;
    background-color: #ffffff;
    border-radius: 8px;
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
}
h1, h2, h3, h4, h5, h6 {
    color: #122C39;
}
.stButton>button {
    background-color: #122C39;
    color: #ffffff;
    border-radius: 5px;
    padding: 10px 20px;
    border: none;
    transition: background-color 0.3s ease;
}
.stButton>button:hover {
    background-color: #0f1e28;
}
.stDataFrame {
    border: 1px solid #ddd;
    border-radius: 5px;
    overflow: hidden;
}
footer {visibility: hidden;}
::-webkit-scrollbar {width: 8px;}
::-webkit-scrollbar-track {background: #f1f1f1;}
::-webkit-scrollbar-thumb {background: #888; border-radius: 10px;}
::-webkit-scrollbar-thumb:hover {background: #555;}
span.match { color: green; font-weight: bold; }
span.mismatch { color: red; font-weight: bold; }
span.motif { background-color: yellow; font-weight: bold; }
</style>
"""


def apply_styles() -> None:
    """Inject the application CSS styles."""
    st.markdown(CUSTOM_CSS, unsafe_allow_html=True)
