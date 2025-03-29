import base64
from io import StringIO

import plotly.graph_objects as go
import plotly.io as pio
import streamlit as st
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord

from style_css import style

# Apply custom styles
style()

st.markdown(
    """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
""",
    unsafe_allow_html=True,
)

st.title("🌿 Phylogenetic Analysis")

# Reset session on page load
if (
    "last_page" not in st.session_state
    or st.session_state["last_page"] != "Phylogenetic_Analysis"
):
    st.session_state.clear()
    st.session_state["last_page"] = "Phylogenetic_Analysis"

# Session state for key elements
if "alignment" not in st.session_state:
    st.session_state["alignment"] = None
if "tree" not in st.session_state:
    st.session_state["tree"] = None
if "distance_matrix" not in st.session_state:
    st.session_state["distance_matrix"] = None

# Upload file
uploaded_file = st.file_uploader(
    "Upload a DNA Sequence File (FASTA, TXT, RTF)", type=["fasta", "txt", "rtf"]
)

# Choose tree construction method and options
method = st.selectbox(
    "🌳 Select Tree Construction Method", ["Neighbor Joining", "UPGMA"]
)
show_branch_lengths = st.checkbox("📏 Show Branch Lengths", value=True)
show_support_values = st.checkbox("📊 Show Support Values", value=False)

# Run button
analyze_clicked = st.button("🔬 Analyze Phylogenetic Tree")

# Trigger analysis only when user clicks
if uploaded_file and analyze_clicked:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    sequences = list(SeqIO.parse(stringio, "fasta-pearson"))

    if sequences:
        st.subheader("📄 Identified Sequences:")
        for seq in sequences:
            st.write(f"✅ {seq.id} - {len(seq.seq)} bp")

        # Pad and align
        max_len = max(len(s.seq) for s in sequences)
        padded = [str(s.seq).ljust(max_len, "-") for s in sequences]
        aligned = MultipleSeqAlignment(
            [SeqRecord(seq, id=s.id, description="") for seq, s in zip(padded, sequences)]
        )

        st.session_state["alignment"] = aligned

        # Compute distance matrix and tree
        calculator = DistanceCalculator("identity")
        dist_matrix = calculator.get_distance(aligned)
        constructor = DistanceTreeConstructor()
        tree = (
            constructor.nj(dist_matrix)
            if method == "Neighbor Joining"
            else constructor.upgma(dist_matrix)
        )

        # Save to session
        st.session_state["distance_matrix"] = dist_matrix
        st.session_state["tree"] = tree
        st.success("✅ Phylogenetic Analysis Completed!")

# Show results
if st.session_state["alignment"] and st.session_state["tree"]:
    alignment = st.session_state["alignment"]
    tree = st.session_state["tree"]

    st.subheader("📏 Distance Matrix")
    st.dataframe(st.session_state["distance_matrix"])

    # Draw tree
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
                textposition="right center",
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

    # Download as Newick
    tree_file = StringIO()
    Phylo.write(tree, tree_file, "newick")
    st.download_button(
        "📥 Download Tree (Newick)",
        tree_file.getvalue(),
        file_name="phylogenetic_tree.nwk",
    )

    # Export as image
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
