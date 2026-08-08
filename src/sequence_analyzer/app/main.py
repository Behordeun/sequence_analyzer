"""Streamlit application launcher. Called by the CLI `serve` subcommand."""

import sys
from pathlib import Path


def run() -> None:
    """Launch the Streamlit application."""
    from streamlit.web.cli import main_run

    app_path = str(Path(__file__).parent / "pages" / "homepage.py")
    sys.argv = ["streamlit", "run", app_path, "--server.headless=true"]
    main_run(app_path)
