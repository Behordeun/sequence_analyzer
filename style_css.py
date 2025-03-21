import streamlit as st


def style():
    # Custom CSS for layout adjustments based on config.toml color scheme
    st.markdown(
        """
        <style>
        /* Apply theme colors */
        body {
            background-color: #FFFFFF; /* White background */
            color: black; /* Text color */
            font-family: 'sans-serif'; /* Font style */
        }

        .logo {
            position: absolute;
            top: 10px;
            left: 10px;
            width: 150px; /* Adjust the logo size */
        }

        .main-content {
            margin-top: 50px;
            padding: 20px;
            background-color: #122C39; /* Secondary background color */
            border-radius: 8px;
            box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.1);
        }

        .text-content {
            font-size: 18px;
            line-height: 1.8;
            color: black; /* Text color */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
