import base64
from collections import Counter
from io import StringIO

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import streamlit as st
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord

from style_css import style

# 🎨 Apply styling
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

# Reset session
if (
    "last_page" not in st.session_state
    or st.session_state["last_page"] != "Phylogenetic_Analysis"
):
    st.session_state.clear()
    st.session_state["last_page"] = "Phylogenetic_Analysis"

if "alignment" not in st.session_state:
    st.session_state["alignment"] = None
if "tree" not in st.session_state:
    st.session_state["tree"] = None
if "distance_matrix" not in st.session_state:
    st.session_state["distance_matrix"] = None

uploaded_file = st.file_uploader(
    "Upload a DNA Sequence File (FASTA, TXT, RTF)", type=["fasta", "txt", "rtf"]
)

method = st.selectbox(
    "🌳 Select Tree Construction Method", ["Neighbor Joining", "UPGMA"]
)
show_branch_lengths = st.checkbox("📏 Show Branch Lengths", value=True)
show_support_values = st.checkbox("📊 Show Support Values", value=False)
analyze_clicked = st.button("🔬 Analyze Phylogenetic Tree")


# Convert DistanceMatrix to pandas DataFrame
def distance_matrix_to_dataframe(dist_matrix):
    names = dist_matrix.names
    matrix = dist_matrix.matrix

    full_matrix = []
    for i in range(len(names)):
        row = matrix[i] + [0.0] * (len(names) - len(matrix[i]))
        full_matrix.append(row)

    return pd.DataFrame(full_matrix, index=names, columns=names)


# Analysis
if uploaded_file and analyze_clicked:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    sequences = list(SeqIO.parse(stringio, "fasta-pearson"))

    if sequences:
        st.subheader("📄 Identified Sequences:")
        for seq in sequences:
            st.write(f"✅ {seq.id} - {len(seq.seq)} bp")

        # Ensure unique IDs
        name_counter = Counter()
        unique_ids = []
        for seq in sequences:
            base_id = seq.id
            name_counter[base_id] += 1
            new_id = (
                f"{base_id}_{name_counter[base_id]}"
                if name_counter[base_id] > 1
                else base_id
            )
            unique_ids.append(new_id)

        if any(c > 1 for c in name_counter.values()):
            st.warning("⚠️ Duplicate sequence IDs were found and automatically renamed.")

        max_len = max(len(s.seq) for s in sequences)
        padded = [str(s.seq).ljust(max_len, "-") for s in sequences]

        aligned = MultipleSeqAlignment(
            [
                SeqRecord(seq, id=new_id, description="")
                for seq, new_id in zip(padded, unique_ids)
            ]
        )

        st.session_state["alignment"] = aligned
        calculator = DistanceCalculator("identity")
        dist_matrix = calculator.get_distance(aligned)
        constructor = DistanceTreeConstructor()
        tree = (
            constructor.nj(dist_matrix)
            if method == "Neighbor Joining"
            else constructor.upgma(dist_matrix)
        )

        st.session_state["tree"] = tree
        st.session_state["distance_matrix"] = dist_matrix

        st.success("✅ Phylogenetic Analysis Completed!")

# Display output
if st.session_state["alignment"] and st.session_state["tree"]:
    tree = st.session_state["tree"]
    dist_matrix = st.session_state["distance_matrix"]
    df_matrix = distance_matrix_to_dataframe(dist_matrix)

    # Display matrix
    st.subheader("📏 Distance Matrix")
    st.dataframe(df_matrix)

    # Download as CSV
    csv_bytes = df_matrix.to_csv().encode("utf-8")
    st.download_button(
        "📥 Download Matrix as CSV", csv_bytes, file_name="distance_matrix.csv"
    )

    # 📊 Heatmap
    st.subheader("🌡️ Distance Matrix Heatmap")
    fig_heatmap = px.imshow(
        df_matrix,
        text_auto=".2f",
        color_continuous_scale="YlGnBu",
        labels=dict(color="Distance"),
        aspect="auto",
    )
    fig_heatmap.update_layout(
        xaxis_title="", yaxis_title="", xaxis=dict(tickangle=45), height=600
    )
    st.plotly_chart(fig_heatmap)

    # 🌳 Tree
    st.subheader("🌳 Phylogenetic Tree (Interactive)")

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
    fig_tree = go.Figure()

    for parent, child in edges:
        x0, y0 = coords[parent]
        x1, y1 = coords[child]
        hover_text = ""
        if show_branch_lengths and child.branch_length is not None:
            hover_text += f"Length: {child.branch_length:.3f} "
        if show_support_values and hasattr(child, "confidence") and child.confidence:
            hover_text += f"Support: {child.confidence:.1f}"

        fig_tree.add_trace(
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
        label = labels[clade]

        # Truncate label if it's too long
        visible_label = label if len(label) <= 20 else label[:17] + "..."
        
        fig_tree.add_trace(go.Scatter(
            x=[x],
            y=[-y],
            mode='markers+text',
            text=[visible_label],
            hovertext=[label],  # full label on hover
            #textangle=0,
            textfont=dict(size=10),
            textposition='middle right',
            marker=dict(color='blue', size=8),
            hoverinfo='text'
        ))

    fig_tree.update_layout(
        showlegend=False,
        title="Phylogenetic Tree",
        margin=dict(l=50, r=200, t=50, b=50),  # Extra space on right
        xaxis=dict(showticklabels=False),
        yaxis=dict(showticklabels=False),
        height=max(600, len(coords) * 30),  # dynamic height for more tips
    )

    st.plotly_chart(fig_tree)

    # Download tree in Newick
    tree_file = StringIO()
    Phylo.write(tree, tree_file, "newick")
    st.download_button(
        "📥 Download Tree (Newick)",
        tree_file.getvalue(),
        file_name="phylogenetic_tree.nwk",
    )

    # Export tree as image
    export_format = st.selectbox(
        "📁 Select Tree Image Format", ["PNG", "JPEG", "SVG", "PDF"]
    )
    if st.button("📤 Export Tree Image"):
        try:
            image_bytes = pio.to_image(
                fig_tree, format=export_format.lower(), width=1000, height=600
            )
            b64 = base64.b64encode(image_bytes).decode()
            mime_type = (
                "application/pdf"
                if export_format == "PDF"
                else f"image/{export_format.lower()}"
            )
            href = f'<a href="data:{mime_type};base64,{b64}" download="phylogenetic_tree.{export_format.lower()}">📥 Download as {export_format}</a>'
            st.markdown(href, unsafe_allow_html=True)
        except Exception as e:
            st.error(f"❌ Export failed: {e}")

# Footer
st.markdown("---")
st.markdown(
    """
    <p style="color: black; text-align: center; font-size: 15px;">
        Developed by 
        <a href="https://github.com/Behordeun" target="_blank" style="color: blue;">Behordeun</a> 
        and 
        <a href="https://github.com/bollergene" target="_blank" style="color: blue;">Bollergene</a>.
    </p>
""",
    unsafe_allow_html=True,
)
st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">📞 +2348108316393</p>""",
    unsafe_allow_html=True,
)
st.markdown(
    """<p style="color:black; text-align:center;font-size:15px;">Copyright | Behordeun 2025(c)</p>""",
    unsafe_allow_html=True,
)
