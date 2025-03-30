import base64
from io import StringIO

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import streamlit as st
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord

from utils import convert_to_fasta

# Apply custom styles
st.set_page_config(layout="wide")
st.markdown(
    """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🌿 Phylogenetic Analysis with Enhanced Visuals")

# Reset session on page change
if (
    "last_page" not in st.session_state
    or st.session_state["last_page"] != "Phylogenetic_Analysis"
):
    st.session_state.clear()
    st.session_state["last_page"] = "Phylogenetic_Analysis"

# Initialize alignment storage
st.session_state.setdefault("alignment", None)

# Upload file
uploaded_file = st.file_uploader(
    "Upload a DNA Sequence File (FASTA, TXT, RTF, PHYLIP, NEXUS)",
    type=["fasta", "txt", "rtf", "phy", "nex"],
)

# Select tree construction method
method = st.selectbox("🌳 Tree Construction Method", ["Neighbor Joining", "UPGMA"])
layout_style = st.selectbox(
    "🎨 Tree Layout Style", ["rectangular", "radial", "circular"]
)

# Optional toggles
show_branch_lengths = st.checkbox("📏 Show Branch Lengths", value=True)
show_support_values = st.checkbox("📊 Show Support Values", value=False)
show_distance_matrix = st.checkbox("📐 Show Distance Matrix Heatmap", value=True)

# Alignment process
if uploaded_file:
    content = convert_to_fasta(uploaded_file)
    sequences = list(SeqIO.parse(StringIO(content), "fasta"))

    st.subheader("📄 Identified Sequences:")
    for seq in sequences:
        st.write(f"✅ {seq.id} - {len(seq.seq)} bp")

    if st.button("🧬 Build Phylogenetic Tree"):
        max_len = max(len(s.seq) for s in sequences)
        padded = [Seq(str(s.seq).ljust(max_len, "-")) for s in sequences]
        aligned = MultipleSeqAlignment(
            [
                SeqRecord(seq, id=s.id, description="")
                for seq, s in zip(padded, sequences)
            ]
        )
        st.session_state["alignment"] = aligned
        st.success("✅ Alignment Completed!")

# Tree construction and visualization
if st.session_state["alignment"]:
    alignment = st.session_state["alignment"]

    st.subheader("📏 Distance Matrix")
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    df_matrix = pd.DataFrame(
        data=distance_matrix.matrix,
        index=distance_matrix.names,
        columns=distance_matrix.names
    )

    if show_distance_matrix:
        import plotly.express as px

        fig_heatmap = px.imshow(
            df_matrix, text_auto=".2f", color_continuous_scale="Viridis"
        )
        st.plotly_chart(fig_heatmap, use_container_width=True)

    constructor = DistanceTreeConstructor()
    tree = (
        constructor.nj(distance_matrix)
        if method == "Neighbor Joining"
        else constructor.upgma(distance_matrix)
    )

    st.subheader("🌳 Phylogenetic Tree (Interactive Plot)")

    def get_plotly_tree_data(tree):
        edges, labels, coords = [], {}, {}

        def traverse(clade, x=0, y=0, depth=0):
            label = str(clade.name) if clade.name else f"Node{depth}-{y}"
            labels[clade] = label
            coords[clade] = (x, y)
            if clade.clades:
                for i, child in enumerate(clade.clades):
                    new_x = x + (clade.branch_length or 0.1)
                    new_y = y + (i - len(clade.clades) / 2) * 0.5
                    edges.append((clade, child))
                    traverse(child, new_x, new_y, depth + 1)

        traverse(tree.root)
        return edges, labels, coords

    edges, labels, coords = get_plotly_tree_data(tree)
    fig = go.Figure()

    for parent, child in edges:
        x0, y0 = coords[parent]
        x1, y1 = coords[child]
        hover_text = ""
        if show_branch_lengths and child.branch_length is not None:
            hover_text += f"Length: {child.branch_length:.3f} "
        if show_support_values and hasattr(child, "confidence") and child.confidence:
            hover_text += f"Support: {child.confidence:.1f}"

        fig.add_trace(
            go.Scatter(
                x=[x0, x1],
                y=[-y0, -y1],
                mode="lines",
                line=dict(color="green", width=2),
                hoverinfo="text" if hover_text else "none",
                text=[hover_text],
            )
        )

    for clade, (x, y) in coords.items():
        fig.add_trace(
            go.Scatter(
                x=[x],
                y=[-y],
                mode="markers+text",
                text=[labels[clade]],
                textposition="middle right",
                marker=dict(color="blue", size=8),
                hoverinfo="text",
            )
        )

    fig.update_layout(
        showlegend=False,
        title="Phylogenetic Tree",
        xaxis=dict(showticklabels=False, zeroline=False),
        yaxis=dict(showticklabels=False, zeroline=False),
        height=600,
    )

    st.plotly_chart(fig)

    # 📥 Newick download
    tree_file = StringIO()
    Phylo.write(tree, tree_file, "newick")
    st.download_button(
        "📥 Download Tree (Newick)",
        tree_file.getvalue(),
        file_name="phylogenetic_tree.nwk",
    )

    # 📥 Image export
    formats = ["PNG", "JPEG", "SVG", "PDF"]
    selected_format = st.selectbox("📁 Select Tree Image Format", formats)

    if st.button("📤 Export Tree Image"):
        try:
            image_bytes = pio.to_image(
                fig, format=selected_format.lower(), width=1000, height=600
            )
            b64 = base64.b64encode(image_bytes).decode()
            mime_type = (
                "application/pdf"
                if selected_format == "PDF"
                else f"image/{selected_format.lower()}"
            )
            href = f'<a href="data:{mime_type};base64,{b64}" download="phylogenetic_tree.{selected_format.lower()}">📥 Download as {selected_format}</a>'
            st.markdown(href, unsafe_allow_html=True)
        except Exception as e:
            st.error(f"❌ Export failed: {e}")

# Footer
st.markdown("---")
st.markdown(
    """
    <p style="color: black; text-align: center; font-size: 15px;">
        Developed by 
        <a href="https://github.com/Behordeun" target="_blank" style="color: blue; text-decoration: none;">Behordeun</a> 
        and 
        <a href="https://github.com/bollergene" target="_blank" style="color: blue; text-decoration: none;">Bollergene</a>.
    </p>
""",
    unsafe_allow_html=True,
)
st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">
📞 +2348108316393
""",
    unsafe_allow_html=True,
)
st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">
Copyright | Behordeun 2025(c)
""",
    unsafe_allow_html=True,
)
