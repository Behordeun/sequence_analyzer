import streamlit as st


def style():
    # Custom CSS for layout adjustments based on config.toml color scheme
    st.markdown(
        """
        <style>
        /* Global styles */
        body {
            background-color: #f4f4f9; /* Light grey background for better readability */
            color: #333333; /* Dark grey text for contrast */
            font-family: 'Arial', sans-serif; /* Modern and clean font */
        }

        /* Logo styling */
        .logo {
            position: absolute;
            top: 10px;
            left: 10px;
            width: 150px; /* Adjust the logo size */
        }

        /* Main content styling */
        .main-content {
            margin-top: 50px;
            padding: 20px;
            background-color: #ffffff; /* White background for content */
            border-radius: 8px; /* Rounded corners for a softer look */
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1); /* Subtle shadow for depth */
        }

        /* Header styling */
        h1, h2, h3, h4, h5, h6 {
            color: #122C39; /* Consistent color for headers */
        }

        /* Button styling */
        .stButton>button {
            background-color: #122C39; /* Primary color for buttons */
            color: #ffffff; /* White text for contrast */
            border-radius: 5px; /* Rounded corners for buttons */
            padding: 10px 20px; /* Adequate padding for buttons */
            border: none; /* Remove default border */
            transition: background-color 0.3s ease; /* Smooth transition for hover effect */
        }

        .stButton>button:hover {
            background-color: #0f1e28; /* Darker shade on hover */
        }

        /* Table styling */
        .stDataFrame {
            border: 1px solid #ddd; /* Light border for tables */
            border-radius: 5px; /* Rounded corners for tables */
            overflow: hidden; /* Ensure content fits within the table */
        }

        /* Footer styling */
        footer {
            visibility: hidden; /* Hide default Streamlit footer */
        }

        /* Custom scrollbar */
        ::-webkit-scrollbar {
            width: 8px;
        }

        ::-webkit-scrollbar-track {
            background: #f1f1f1;
        }

        ::-webkit-scrollbar-thumb {
            background: #888;
            border-radius: 10px;
        }

        ::-webkit-scrollbar-thumb:hover {
            background: #555;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
