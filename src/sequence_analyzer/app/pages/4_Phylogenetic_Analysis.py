"""Phylogenetic Analysis page: tree construction, bootstrapping, entropy."""

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from sequence_analyzer.app.components import hide_streamlit_chrome, render_footer
from sequence_analyzer.app.styles import apply_styles
from sequence_analyzer.core.phylogenetics import bootstrap_tree, build_tree, compute_entropy
from sequence_analyzer.core.validation import validate_sequences
from sequence_analyzer.io.parsers import parse_sequence_file

st.set_page_config(layout="wide")
apply_styles()
hide_streamlit_chrome()

st.title("🌿 Phylogenetic Analysis")

for key in ["tree_result", "entropy_map"]:
    st.session_state.setdefault(key, None)

uploaded_file = st.file_uploader(
    "Upload Sequence Files", type=["fasta", "txt", "rtf", "phy", "nex"]
)

if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8")
    records = parse_sequence_file(content, format_hint="auto")
    valid = validate_sequences(records, seq_type="auto")
    st.success(f"{len(valid)} sequences uploaded")

    if st.checkbox("Preview Sequences"):
        for s in valid:
            st.write(f"{s.id} - {len(s.seq)} bp")
            st.code(str(s.seq), language="text")

    method = st.selectbox("Tree Method", ["nj", "upgma"])
    do_bootstrap = st.checkbox("Enable Bootstrapping")
    reps = st.slider("Bootstrap Replicates", 10, 100, 50) if do_bootstrap else 0

    if st.button("Build Tree"):
        result = build_tree(valid, method=method)
        st.session_state["tree_result"] = result

        if do_bootstrap:
            support = bootstrap_tree(result.alignment, method=method, replicates=reps)
            for clade in result.tree.find_clades():
                tips = frozenset(x.name for x in clade.get_terminals())
                clade.confidence = support.get(tips, 0.0)

        entropy_map = compute_entropy(result.tree)
        st.session_state["entropy_map"] = entropy_map
        st.success("Tree built successfully!")

if st.session_state["tree_result"]:
    result = st.session_state["tree_result"]

    # Entropy table
    if st.checkbox("Show Entropy Table"):
        entropy_map = st.session_state["entropy_map"] or {}
        df_entropy = pd.DataFrame(
            [{"Clade": k, "Entropy": round(v, 4)} for k, v in entropy_map.items()]
        )
        st.dataframe(df_entropy)

    # Simple tree visualization
    st.subheader("Phylogenetic Tree")

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
                # Horizontal line
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
        margin={"l": 10, "r": 10, "t": 30, "b": 10},
    )
    st.plotly_chart(fig, use_container_width=True)

    # Download
    export_fmt = st.selectbox("Export Format", ["Newick", "JSON"])
    if export_fmt == "Newick":
        st.download_button("Download Tree (Newick)", result.newick, file_name="tree.nwk")
    elif export_fmt == "JSON":
        import json

        data = {"newick": result.newick, "terminals": [t.name for t in result.tree.get_terminals()]}
        st.download_button("Download JSON", json.dumps(data, indent=2), file_name="tree.json")

render_footer()
