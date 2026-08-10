"""Streamlit application launcher. Called by the CLI `serve` subcommand."""

import sys
from pathlib import Path


def run() -> None:
    """Launch the Streamlit application."""
    from streamlit.web.cli import main_run

    app_path = str(Path(__file__).parent / "streamlit_app.py")
    original_argv = sys.argv.copy()
    sys.argv = ["streamlit", "run", app_path, "--server.headless=true"]
    try:
        main_run(app_path)
    finally:
        sys.argv = original_argv
