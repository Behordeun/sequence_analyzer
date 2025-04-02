import base64
import json
import random
from collections import Counter
from io import BytesIO, StringIO

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import streamlit as st
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

from utils import convert_to_fasta

st.set_page_config(layout="wide")
st.title("🌿 Phylogenetic Analysis")

# --- Initialize State ---
for key in ["alignment", "tree", "entropy_map", "annotations", "cluster_groups"]:
    st.session_state.setdefault(key, None)

# --- Uploads ---
uploaded_file = st.file_uploader(
    "📥 Upload Sequence Files", type=["fasta", "txt", "rtf", "phy", "nex"]
)

if uploaded_file:
    content = convert_to_fasta(uploaded_file)
    sequences = list(SeqIO.parse(StringIO(content), "fasta-blast"))
    st.success(f"✅ {len(sequences)} sequences uploaded")

    if st.checkbox("👁 Preview Sequences"):
        for s in sequences:
            st.write(f"🔹 {s.id} - {len(s.seq)} bp")
            st.code(str(s.seq), language="text")

global layout_style

# --- Tree Construction Section ---
if uploaded_file:
        if "tree_method" not in st.session_state:
            st.session_state.tree_method = "Neighbor Joining"

        if "layout_style" not in st.session_state:
            st.session_state.layout_style = "rectangular"

        st.session_state.tree_method = st.selectbox(
            "🌳 Tree Method", ["Neighbor Joining", "UPGMA"], key="select_tree_method"
        )

        st.session_state.layout_style = st.selectbox(
            "🎨 Layout", ["rectangular", "radial", "circular"], key="select_layout_style"
        )
        run_tree = st.button("🧬 Build Tree")

        method = st.session_state.tree_method
        layout_style = st.session_state.layout_style
        
        if run_tree:
            max_len = max(len(s.seq) for s in sequences)
            seen_ids = {}
            padded = []
            for s in sequences:
                sid = s.id
                if sid in seen_ids:
                    seen_ids[sid] += 1
                    sid = f"{sid}_{seen_ids[sid]}"
                else:
                    seen_ids[sid] = 1
                padded.append(SeqRecord(Seq(str(s.seq).ljust(max_len, "-")), id=sid))
    
            alignment = MultipleSeqAlignment(padded)
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor()
            tree = (
                constructor.nj(dm)
                if method == "Neighbor Joining"
                else constructor.upgma(dm)
            )
    
            st.session_state["alignment"] = alignment
            st.session_state["tree"] = tree


if st.session_state.tree:
    st.markdown("### ⚙️ Additional Options")

    # --- Tree Options ---
    #method = st.selectbox("🌳 Tree Method", ["Neighbor Joining", "UPGMA"])
    #layout_style = st.selectbox("🎨 Layout", ["rectangular", "radial", "circular"])
    bootstrap = st.checkbox("🔁 Enable Bootstrapping")
    reps = st.slider("Bootstrap Replicates", 10, 100, 50) if bootstrap else 0
    show_branch_lengths = st.checkbox("📏 Show Branch Lengths", True)
    show_labels = st.checkbox("🗣 Show Labels", True)
    show_support = st.checkbox("📊 Show Bootstrap Values", True)
    show_entropy = st.checkbox("🔥 Show Entropy", True)
    show_heatmap = st.checkbox("📐 Show Distance Matrix Heatmap")
    enable_annotation = st.checkbox("✏️ Enable Clade Annotation")
    enable_compare = st.checkbox("🧪 Enable Tree Comparison")
    enable_search = st.checkbox("🔍 Enable Clade Search")
    enable_metadata = st.checkbox("🧬 Enable Metadata")

    compare_tree = None
    if enable_compare:
        compare_tree = st.file_uploader(
            "📁 Upload Tree for Comparison (.nwk)", type="nwk"
        )

    # --- Metadata Handling (Optional Toggle) ---
    if enable_metadata:
        st.subheader("🧬 Metadata Generation & Export")

        grouping_method = st.selectbox(
            "Select Metadata Grouping Strategy",
            [
                "Group by Sequence Length",
                "Group by Taxonomy Prefix",
                "Cluster with K-Means",
            ], key = "grouping_method"
        )

        metadata = {}
        df_meta = None

        if st.session_state.get("alignment"):
            alignment = st.session_state.alignment
            records = alignment

            if grouping_method == "Group by Sequence Length":
                groups = ["Short" if len(r.seq) < 200 else "Long" for r in records]

            elif grouping_method == "Group by Taxonomy Prefix":
                groups = []
                for r in records:
                    # Try description first, fallback to ID
                    if r.description and r.description != r.id:
                        prefix = r.description.split()[0]
                    else:
                        prefix = r.id.split(".")[0].split("_")[0]
                    groups.append(prefix if prefix else "<unknown>")

            elif grouping_method == "Cluster with K-Means":
                try:
                    from sklearn.cluster import KMeans
                    from sklearn.feature_extraction.text import CountVectorizer

                    vectorizer = CountVectorizer(analyzer="char", ngram_range=(2, 2))
                    X = vectorizer.fit_transform([str(r.seq) for r in records])
                    kmeans = KMeans(n_clusters=2, random_state=42)
                    labels = kmeans.fit_predict(X)
                    groups = [f"Cluster {l+1}" for l in labels]
                except ImportError:
                    st.warning(
                        "⚠️ scikit-learn not installed. Cannot perform clustering."
                    )
                    groups = ["Unclustered"] * len(records)

            # Build DataFrame
            df_meta = pd.DataFrame({"ID": [r.id for r in records], "Group": groups})
            metadata = dict(zip(df_meta["ID"], df_meta["Group"]))
            st.session_state["metadata"] = metadata

            # Preview and download
            if st.checkbox("📋 Preview Generated Metadata"):
                st.dataframe(df_meta)

            st.download_button(
                "⬇️ Download Metadata CSV",
                df_meta.to_csv(index=False),
                file_name="generated_metadata.csv",
                mime="text/csv",
            )

    # --- Metadata Upload (Optional) ---
    metadata = {}
    df_meta = None

    enable_metadata_upload = st.checkbox("📎 Upload Metadata CSV Manually")
    if enable_metadata_upload:
        metadata_file = st.file_uploader("Upload Metadata CSV (ID, Group)", type="csv")
        if metadata_file:
            df_meta = pd.read_csv(metadata_file)
            if "ID" in df_meta.columns and "Group" in df_meta.columns:
                metadata = dict(zip(df_meta["ID"], df_meta["Group"]))
                st.session_state["metadata"] = metadata
                st.success("✅ Metadata file loaded")

                # Optional Preview
                if st.checkbox("📋 Preview Uploaded Metadata"):
                    st.dataframe(df_meta)

                st.download_button(
                    "⬇️ Download Uploaded Metadata",
                    df_meta.to_csv(index=False),
                    file_name="uploaded_metadata.csv",
                )
            else:
                st.error("❌ Metadata file must contain 'ID' and 'Group' columns.")

    # --- Bootstrap ---
    if st.session_state.tree and bootstrap:

        def bootstrap_trees(aligned, reps):
            trees = []
            for _ in range(reps):
                idx = [
                    random.randint(0, len(aligned[0]) - 1)
                    for _ in range(len(aligned[0]))
                ]
                boot = MultipleSeqAlignment(
                    [
                        SeqRecord(Seq("".join(str(r.seq)[i] for i in idx)), id=r.id)
                        for r in aligned
                    ]
                )
                dboot = DistanceCalculator("identity").get_distance(boot)
                t = (
                    DistanceTreeConstructor().nj(dboot)
                    if method == "Neighbor Joining"
                    else DistanceTreeConstructor().upgma(dboot)
                )
                trees.append(t)
            return trees

        boots = bootstrap_trees(st.session_state["alignment"], reps)
        support = Counter()
        for t in boots:
            for clade in t.find_clades():
                tips = frozenset([x.name for x in clade.get_terminals()])
                support[tips] += 1
        for clade in st.session_state.tree.find_clades():
            tips = frozenset([x.name for x in clade.get_terminals()])
            clade.confidence = support.get(tips, 0) / reps * 100

    # --- Entropy ---
    if show_entropy and st.session_state.tree:
        st.subheader("🔥 Entropy Analysis")
        entropy_map = {}
        entropy_rows = []

        for clade in st.session_state.tree.find_clades():
            if clade.is_terminal():
                entropy_map[clade] = 0
            else:
                tips = [x.name for x in clade.get_terminals()]
                freq = pd.Series(tips).value_counts(normalize=True)
                entropy = -(freq * np.log2(freq)).sum()
                entropy_map[clade] = entropy
            entropy_rows.append(
                {"Clade": clade.name or "-", "Entropy": round(entropy_map[clade], 4)}
            )

        st.session_state["entropy_map"] = entropy_map
        df_entropy = pd.DataFrame(entropy_rows)
        st.dataframe(df_entropy)

        st.download_button(
            "📥 Download Entropy Table (CSV)",
            data=df_entropy.to_csv(index=False),
            file_name="entropy_scores.csv",
        )

    # --- Distance Matrix ---
    if st.session_state["alignment"] and show_heatmap:
        st.subheader("📐 Distance Matrix Heatmap")
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(st.session_state["alignment"])
        df_dm = pd.DataFrame(dm.matrix, index=dm.names, columns=dm.names)
        st.dataframe(df_dm)
        st.download_button(
            "⬇️ Download Matrix", df_dm.to_csv(), file_name="distance_matrix.csv"
        )

        fig_heat = go.Figure(
            data=go.Heatmap(z=df_dm.values, x=df_dm.columns, y=df_dm.index)
        )
        fig_heat.update_layout(title="Pairwise Distance Heatmap")
        st.plotly_chart(fig_heat)

    # --- Clade Annotation ---
    if enable_annotation:
        st.subheader("✏️ Annotate Clades")
        clade_name = st.text_input("Enter Clade Name")
        clade_label = st.text_input("Label for Clade")

        if st.button("✅ Add Annotation"):
            if clade_name and clade_label:
                if (
                    "annotations" not in st.session_state
                    or st.session_state["annotations"] is None
                ):
                    st.session_state["annotations"] = {}
                st.session_state["annotations"][clade_name] = clade_label
                st.success(f"🔖 Tagged '{clade_name}' as '{clade_label}'")
            else:
                st.warning("⚠️ Please provide both clade name and label.")

    # --- Tree Rendering ---
    if st.session_state.tree:

        def layout_coords(tree, layout="rectangular"):
            coords = {}

            def recurse(clade, x=0, y=0):
                coords[clade] = (x, y)
                for i, child in enumerate(clade.clades):
                    dx = child.branch_length or 0.1
                    dy = (i - len(clade.clades) / 2) * 2
                    nx = x + dx if layout == "rectangular" else x + dx * np.cos(dy)
                    ny = y + dy if layout == "rectangular" else y + dx * np.sin(dy)
                    recurse(child, nx, ny)

            recurse(tree.root)
            return coords

        layout_style = st.session_state.layout_style
        
        coords = layout_coords(st.session_state.tree, layout_style)
        fig = go.Figure()

        for p in coords:
            for c in p.clades:
                if c in coords:
                    x0, y0 = coords[p]
                    x1, y1 = coords[c]
                    fig.add_trace(
                        go.Scatter(
                            x=[x0, x1],
                            y=[y0, y1],
                            mode="lines",
                            line=dict(color="gray"),
                        )
                    )

        for clade, (x, y) in coords.items():
            annotations = st.session_state.get("annotations") or {}
            label = annotations.get(clade.name, clade.name or "")
            tooltip = []

            # ✅ Add sequence name
            name = clade.name or "Unnamed"
            tooltip.append(f"🧬 Sequence: {name}")

            # ✅ Add taxonomy group from metadata if available
            if clade.name in metadata:
                tooltip.append(f"🧪 Group: {metadata[clade.name]}")

            # ✅ Add branch length
            if show_branch_lengths and clade.branch_length:
                tooltip.append(f"📏 Branch: {clade.branch_length:.2f}")

            # ✅ Add bootstrap support if available
            if show_support and getattr(clade, "confidence", None):
                tooltip.append(f"📊 Support: {clade.confidence:.1f}%")

            # ✅ Add entropy if enabled
            if show_entropy and st.session_state.entropy_map:
                entropy = st.session_state.entropy_map.get(clade, 0)
                tooltip.append(f"🔥 Entropy: {entropy:.3f}")
            else:
                entropy = 0

            # Color encoding
            color = "blue"
            if show_entropy:
                color = f"rgba({int(entropy*50)},0,255,0.9)"
            if clade.name in metadata:
                color = f"hsl({hash(metadata[clade.name]) % 360}, 70%, 60%)"

            # ✅ Click-to-focus via custom hovertemplate and selection highlighting
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[y],
                    mode="markers+text" if show_labels else "markers",
                    text=[label],
                    marker=dict(size=8, color=color),
                    hovertext="<br>".join(tooltip),
                    hoverinfo="text",
                    textposition="top center",
                    customdata=[clade.name],  # Used for future interactive selections
                )
            )

        fig.update_layout(
            title="🧬 Phylogenetic Tree",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            height=700,
        )
        st.plotly_chart(fig, use_container_width=True)

        # --- Optional Search ---
        if enable_search:
            search = st.text_input("🔎 Enter Clade Name to Search")
            if search:
                found = next(
                    (c for c in coords if c.name and search.lower() in c.name.lower()),
                    None,
                )
                if found:
                    st.success(f"✅ Found clade: {found.name}")
                else:
                    st.warning("❌ Clade not found")

        # --- Download Options ---
        export_fmt = st.selectbox(
            "📁 Export Format", ["Newick", "PNG", "SVG", "PDF", "JSON", "HTML"], key = "export format"
        )
        if export_fmt == "Newick":
            buf = StringIO()
            Phylo.write(st.session_state.tree, buf, "newick")
            st.download_button(
                "⬇️ Download Tree (Newick)", buf.getvalue(), file_name="phylogenetic_tree.nwk"
            )
        elif export_fmt == "JSON":
            data = {
                "nodes": [{"name": c.name} for c in coords],
                "annotations": st.session_state.annotations,
            }
            st.download_button(
                "⬇️ Download JSON", json.dumps(data, indent=2), file_name="phylogenetic_tree.json"
            )
        elif export_fmt == "HTML":
            html = pio.to_html(fig)
            st.download_button("⬇️ Download HTML", html, file_name="phylogenetic_tree.html")
        else:
            img = pio.to_image(fig, format=export_fmt.lower())
            b64 = base64.b64encode(img).decode()
            mime = (
                "application/pdf"
                if export_fmt == "PDF"
                else f"image/{export_fmt.lower()}"
            )
            st.markdown(
                f'<a href="data:{mime};base64,{b64}" download="phylogenetic_tree.{export_fmt.lower()}">📥 Download Tree</a>',
                unsafe_allow_html=True,
            )

        # --- Tree Comparison Rendering ---
        if compare_tree:
            st.subheader("🧪 Tree Comparison")

            # Load comparison tree
            compare_str = compare_tree.read().decode("utf-8")
            try:
                from Bio import Phylo
                from io import StringIO

                comp_tree = Phylo.read(StringIO(compare_str), "newick")
                st.success("✅ Comparison tree loaded")

                def render_tree(tree_obj, title):
                    coords = {}
                    def recurse(clade, x=0, y=0):
                        coords[clade] = (x, y)
                        for i, child in enumerate(clade.clades):
                            dx = child.branch_length or 0.1
                            dy = (i - len(clade.clades) / 2) * 2
                            recurse(child, x + dx, y + dy)

                    recurse(tree_obj.root)
                    fig = go.Figure()
                    for p in coords:
                        for c in p.clades:
                            if c in coords:
                                x0, y0 = coords[p]
                                x1, y1 = coords[c]
                                fig.add_trace(go.Scatter(
                                    x=[x0, x1], y=[y0, y1], mode="lines",
                                    line=dict(color="lightgray")))

                    for clade, (x, y) in coords.items():
                        name = clade.name or ""
                        fig.add_trace(go.Scatter(
                            x=[x], y=[y], mode="markers+text",
                            text=[name], marker=dict(size=6, color="green"),
                            textposition="top center"))

                    fig.update_layout(
                        title=title, height=600,
                        xaxis=dict(visible=False),
                        yaxis=dict(visible=False)
                    )
                    return fig

                col1, col2 = st.columns(2)
                with col1:
                    st.plotly_chart(render_tree(st.session_state.tree, "🔬 Your Tree"))
                with col2:
                    st.plotly_chart(render_tree(comp_tree, "🧪 Comparison Tree"))

            except Exception as e:
                st.error(f"❌ Error parsing comparison tree: {e}")

        # --- PDF Report ---
        if st.button("📄 Generate PDF Report"):
            pdf = BytesIO()
            c = canvas.Canvas(pdf, pagesize=letter)
            c.drawString(72, 750, "🧬 Phylogenetic Analysis Report")
            c.drawString(72, 730, f"Tree Method: {method}")
            c.drawString(72, 710, f"Layout: {layout_style}")
            c.drawString(72, 690, f"Bootstrap: {reps if bootstrap else 'None'}")
            c.save()
            st.download_button(
                "📄 Download Report", pdf.getvalue(), file_name="phylo_report.pdf"
            )

    # --- Footer ---
st.markdown("---")
st.markdown(
    """
<p style="text-align:center;font-size:14px">
Developed by <a href="https://github.com/Behordeun">Behordeun</a> & 
<a href="https://github.com/bollergene">Bollergene</a><br>
📞 +2348108316393 | © Behordeun 2025
</p>
""",
    unsafe_allow_html=True,
)
