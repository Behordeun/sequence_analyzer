"""Pipeline Stage 5: Tree — phylogenetic tree construction."""

import plotly.graph_objects as go
import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.pages.pipeline.help_content import render_help_panel
from sequence_analyzer.app.pages.pipeline.navigation import (
    render_config_sidebar,
    render_progress_bar,
    render_stage_nav,
)
from sequence_analyzer.app.pages.pipeline.state import get_pipeline_state, mark_complete
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.phylogenetics import bootstrap_tree, build_tree, compute_entropy

st.set_page_config(layout="wide", page_title="Pipeline: Tree")
apply_styles()
hide_streamlit_chrome()
render_config_sidebar()

state = get_pipeline_state()
state.stage = "Tree"

render_progress_bar()
st.title("🌿 Stage 5: Phylogenetic Tree")


render_help_panel("Tree")

sequences = state.valid_sequences or state.sequences
if len(sequences) < 2:
    st.warning("At least 2 sequences are required for tree construction.")
    render_stage_nav()
    render_footer()
    st.stop()

config = state.config

if st.button("Build Tree"):
    with st.spinner("Constructing tree..."):
        result = build_tree(sequences, method=config.tree_method)

        # Bootstrap if configured
        if config.bootstrap_replicates >= 10:
            support = bootstrap_tree(
                result.alignment,
                method=config.tree_method,
                replicates=config.bootstrap_replicates,
            )
            for clade in result.tree.find_clades():
                tips = frozenset(x.name for x in clade.get_terminals())
                clade.confidence = support.get(tips, 0.0)

        state.tree_result = result
        mark_complete("Tree")
        st.success("Tree built successfully!")

# --- Display ---
if state.tree_result:
    result = state.tree_result

    # Simple tree rendering
    def _layout_coords(tree):
        coords = {}
        y_counter = [0]

        def recurse(clade, x=0):
            if clade.is_terminal():
                coords[clade] = (x, y_counter[0])
                y_counter[0] += 1
            else:
                child_ys = []
                for child in clade.clades:
                    dx = child.branch_length or 0.1
                    recurse(child, x + dx)
                    child_ys.append(coords[child][1])
                coords[clade] = (x, sum(child_ys) / len(child_ys))

        recurse(tree.root)
        return coords

    coords = _layout_coords(result.tree)
    fig = go.Figure()

    for parent in coords:
        for child in parent.clades:
            if child in coords:
                x0, y0 = coords[parent]
                x1, y1 = coords[child]
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1, x1],
                        y=[y0, y0, y1],
                        mode="lines",
                        line={"color": "gray", "width": 1},
                        showlegend=False,
                    )
                )

    for clade, (x, y) in coords.items():
        if clade.is_terminal():
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[y],
                    mode="markers+text",
                    text=[clade.name or ""],
                    textposition="middle right",
                    marker={"size": 6, "color": "steelblue"},
                    showlegend=False,
                )
            )

    fig.update_layout(
        xaxis={"visible": False},
        yaxis={"visible": False},
        height=max(400, len(coords) * 20),
        margin={"l": 10, "r": 150, "t": 30, "b": 10},
    )
    st.plotly_chart(fig, use_container_width=True)

    # Entropy
    compute_entropy(result.tree)
    terminal_count = len(result.tree.get_terminals())
    st.markdown(f"**Terminals:** {terminal_count} | **Method:** {config.tree_method.upper()}")

    # Download
    st.download_button(
        "Download Tree (Newick)",
        result.newick,
        file_name="tree.nwk",
    )

st.markdown("---")
render_stage_nav()
render_footer()
