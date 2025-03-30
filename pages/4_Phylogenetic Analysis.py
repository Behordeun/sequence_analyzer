import base64
import json
from collections import Counter
from io import StringIO

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import skbio
import streamlit as st
from Bio import SeqIO
from skbio import DistanceMatrix
from skbio.sequence.distance import hamming
from skbio.tree import nj, upgma

st.set_page_config(layout="wide")
st.title("🌿 Bootstrap Phylogenetic Viewer")

for key in ["tree", "entropy", "layout", "records", "dm", "bootstraps"]:
    if key not in st.session_state:
        st.session_state[key] = None

with st.expander("📂 Upload Data", expanded=True):
    fasta_file = st.file_uploader("Upload Aligned Sequences (FASTA)", type="fasta")
    metadata_file = st.file_uploader(
        "📁 Optional Metadata File (CSV with ID, Group)", type="csv"
    )
    export_format = st.selectbox(
        "📄 Export Tree As", ["SVG", "PNG", "PDF", "JPEG", "HTML"]
    )
    save_combined = st.checkbox("📅 Save Complete Session")

# --- Config ---
distance_model = st.selectbox("📀 Distance Metric", ["identity", "hamming"])
tree_method = st.selectbox("🌲 Tree Method", ["Neighbor Joining", "UPGMA"])
tree_layout = st.selectbox("🎨 Layout", ["rectangular", "circular", "radial"])
color_by = st.selectbox("🧬 Color By", ["None", "Entropy", "Metadata"])
animate_tree = st.checkbox("🎥 Animate Transitions")
show_entropy_table = st.checkbox("🧪 Show Clade Entropy")
show_support_hist = st.checkbox("📈 Branch Length Histogram")
bootstrap = st.checkbox("🔁 Run Bootstraps")
bootstraps_n = st.slider("Bootstrap Replicates", 10, 100, 50) if bootstrap else 0
subtree_id = st.text_input("🌿 Extract Subtree (Node Name)", "")
run = st.button("🔬 Analyze")


# --- Helpers ---
def parse_alignment(fasta_bytes):
    records = list(SeqIO.parse(StringIO(fasta_bytes.decode("utf-8")), "fasta"))
    seqs = [skbio.DNA(str(r.seq), metadata={"id": r.id}) for r in records]
    return seqs, records


def compute_distance(seqs, model):
    base_len = len(seqs[0])
    filtered = [
        s for s in seqs if len(s) == base_len and set(str(s)).isdisjoint({"-", "N"})
    ]
    if len(filtered) < 2:
        st.error("❌ Not enough clean sequences.")
        st.stop()
    if len(filtered) != len(seqs):
        st.warning(f"⚠️ Skipped {len(seqs) - len(filtered)} invalid sequences.")

    if model == "identity":
        return DistanceMatrix.from_iterable(
            filtered,
            metric=lambda x, y: 1
            - sum(a == b for a, b in zip(str(x), str(y))) / len(x),
            key=lambda x: x.metadata["id"],
        )
    elif model == "hamming":
        return DistanceMatrix.from_iterable(
            filtered, metric=hamming, key=lambda x: x.metadata["id"]
        )


def build_tree(dm):
    return nj(dm) if tree_method == "Neighbor Joining" else upgma(dm)


def compute_entropy(tree):
    scores = []
    for node in tree.non_tips():
        tips = [t.name for t in node.tips()]
        freqs = pd.Series(tips).value_counts(normalize=True)
        entropy = -(freqs * freqs.apply(lambda p: np.log2(p))).sum()
        scores.append((node.name or "internal", entropy))
    return scores


def bootstrap_trees(seqs, reps):
    df = pd.DataFrame(
        [list(str(s)) for s in seqs], index=[s.metadata["id"] for s in seqs]
    )
    trees = []
    for _ in range(reps):
        resample = df.sample(n=df.shape[1], axis=1, replace=True)
        boot = [
            skbio.DNA("".join(str(b) if b else "-" for b in row), metadata={"id": i})
            for i, row in resample.iterrows()
        ]
        dm = compute_distance(boot, distance_model)
        trees.append(build_tree(dm))
    return trees


def annotate_support(main_tree, boots):
    support = Counter()
    for t in boots:
        for node in t.non_tips():
            key = frozenset(tip.name for tip in node.tips())
            support[key] += 1
    for node in main_tree.non_tips():
        key = frozenset(t.name for t in node.tips())
        node.name = f"{support[key]/len(boots)*100:.1f}%"
    return main_tree


def render_tree(tree, layout_style, color_by=None, metadata_map=None, entropy_map=None):
    coords, edges = {}, []

    def traverse(n, x=0, y=0):
        coords[n] = (x, y)
        for i, c in enumerate(n.children):
            nx, ny = (x + 1, y + (i - len(n.children) / 2))
            edges.append((n, c))
            traverse(c, nx, ny)

    traverse(tree.root)
    fig = go.Figure()
    for p, c in edges:
        x0, y0 = coords[p]
        x1, y1 = coords[c]
        fig.add_trace(
            go.Scatter(x=[x0, x1], y=[-y0, -y1], mode="lines", line=dict(color="gray"))
        )
    for n, (x, y) in coords.items():
        name = n.name or ""
        label = name
        color = "blue"
        if metadata_map and name in metadata_map:
            label += f" ({metadata_map[name]})"
            color = f"hsl({hash(metadata_map[name])%360},70%,60%)"
        elif entropy_map and name in entropy_map:
            color = f"rgba({int(entropy_map[name]*50)},0,255,0.8)"
        fig.add_trace(
            go.Scatter(
                x=[x],
                y=[-y],
                mode="markers+text",
                text=[label],
                textposition="middle right",
                marker=dict(size=8, color=color),
                hovertext=label,
                hoverinfo="text",
            )
        )
    fig.update_layout(
        title="Phylogenetic Tree",
        height=800,
        xaxis=dict(showticklabels=False),
        yaxis=dict(showticklabels=False),
        dragmode="pan",
    )
    return fig


# --- Main Logic ---
if run and fasta_file:
    seqs, records = parse_alignment(fasta_file.read())
    dm = compute_distance(seqs, distance_model)
    tree = build_tree(dm)
    if bootstrap:
        boots = bootstrap_trees(seqs, bootstraps_n)
        tree = annotate_support(tree, boots)
    st.session_state.tree = tree
    st.session_state.records = records
    st.session_state.dm = dm
    if color_by == "Entropy":
        entropy_list = compute_entropy(tree)
        st.session_state.entropy = dict(entropy_list)
    else:
        st.session_state.entropy = None
    st.success("✅ Tree Constructed!")

if st.session_state.tree:
    metadata_map = {}
    if metadata_file:
        df = pd.read_csv(metadata_file)
        if "ID" in df.columns and "Group" in df.columns:
            metadata_map = dict(zip(df["ID"], df["Group"]))

    fig = render_tree(
        st.session_state.tree,
        layout_style=tree_layout,
        color_by=color_by,
        metadata_map=metadata_map,
        entropy_map=st.session_state.entropy,
    )
    st.plotly_chart(fig, use_container_width=True)

    if show_entropy_table and st.session_state.entropy:
        df = pd.DataFrame(st.session_state.entropy.items(), columns=["Node", "Entropy"])
        st.dataframe(df.sort_values("Entropy", ascending=False))

    if show_support_hist:
        lengths = [
            c.length for c in st.session_state.tree.non_tips() if c.length is not None
        ]
        st.subheader("📈 Branch Length Distribution")
        st.bar_chart(pd.Series(lengths))

    if st.button("📥 Export Tree Image"):
        try:
            img = pio.to_image(fig, format=export_format.lower())
            b64 = base64.b64encode(img).decode()
            mime = (
                "application/pdf"
                if export_format == "PDF"
                else f"image/{export_format.lower()}"
            )
            st.markdown(
                f'<a href="data:{mime};base64,{b64}" download="tree.{export_format.lower()}">Download Tree</a>',
                unsafe_allow_html=True,
            )
        except Exception as e:
            st.error(f"Export failed: {e}")

    # Export tree as Nexus
    from io import StringIO

    from skbio.tree import TreeNode

    tree_str = st.session_state.tree.to_newick()
    tree_node = TreeNode.read(StringIO(tree_str))
    tree_nexus = StringIO()
    tree_node.write(tree_nexus, format="nexus")
    st.download_button(
        "📥 Download Tree (.nex)", tree_nexus.getvalue(), file_name="tree.nex"
    )

    if save_combined:
        session = {
            "tree": st.session_state.tree.to_newick(),
            "records": [r.format("fasta") for r in st.session_state.records],
            "matrix": pd.DataFrame(
                st.session_state.dm.data,
                index=st.session_state.dm.ids,
                columns=st.session_state.dm.ids,
            ).to_json(),
            "metadata": metadata_map,
            "layout": tree_layout,
        }
        st.download_button(
            "📁 Download Session (.json)",
            json.dumps(session, indent=2),
            file_name="session.json",
        )

    if subtree_id:
        try:
            subtree = st.session_state.tree.find(subtree_id)
            st.subheader(f"🌿 Subtree from {subtree_id}")
            sub_fig = render_tree(subtree, layout_style=tree_layout)
            st.plotly_chart(sub_fig, use_container_width=True)
            newick = subtree.to_newick()
            st.download_button(
                "📥 Download Subtree (Newick)",
                newick,
                file_name=f"subtree_{subtree_id}.nwk",
            )
        except Exception:
            st.error("❌ Could not extract specified subtree.")

st.markdown("---")
st.markdown(
    """
<p style="text-align:center;font-size:14px">
Developed by <a href="https://github.com/Behordeun">Behordeun</a> and <a href="https://github.com/bollergene">Bollergene</a><br>
📞 +2348108316393 | © Behordeun 2025
</p>
""",
    unsafe_allow_html=True,
)
