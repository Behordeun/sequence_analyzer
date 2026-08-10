"""Contextual help content for each pipeline stage.

Each stage has an explanation of what it does and why it matters.
Key metrics have tooltip-style definitions. QC flags have explanations
of what triggered them and what the user should consider.
"""

import streamlit as st

STAGE_HELP: dict[str, dict[str, str]] = {
    "Ingest": {
        "title": "About Sequence Ingestion",
        "content": (
            "This stage loads your raw sequence data into the pipeline. "
            "You can upload files in common bioinformatics formats (FASTA, CLUSTAL, "
            "NEXUS, PHYLIP) or fetch sequences directly from NCBI GenBank using "
            "accession numbers.\n\n"
            "**Supported formats:**\n"
            "- **FASTA** — the most common; starts with `>`\n"
            "- **CLUSTAL** — alignment output from ClustalW/Omega\n"
            "- **NEXUS** — used by PAUP, MrBayes, and other phylogenetics tools\n"
            "- **PHYLIP** — compact format with a header line of `ntaxa nalign`\n\n"
            "**Tips:**\n"
            "- You can upload multiple files at once; they'll be concatenated\n"
            "- For GenBank, use accession numbers like `NM_001301717` or `MH156967`\n"
            "- Files over 50MB are rejected to keep the app responsive"
        ),
    },
    "QC": {
        "title": "About Quality Control",
        "content": (
            "Quality control checks each sequence against configurable thresholds "
            "before analysis begins. Catching problems here prevents misleading "
            "results downstream.\n\n"
            "**What gets checked:**\n"
            "- **Ambiguity rate** — percentage of ambiguous bases (N, R, Y, etc.). "
            "High ambiguity suggests sequencing errors or contamination.\n"
            "- **Length** — very short sequences may not carry enough information "
            "for meaningful alignment or phylogenetics.\n"
            "- **Gap fraction** — sequences dominated by gaps are likely artifacts "
            "of alignment or assembly problems.\n"
            "- **Outliers** — sequences whose length or GC content is >2 standard "
            "deviations from the batch mean may come from a different organism.\n\n"
            "**Statuses:**\n"
            "- **Pass** — meets all thresholds\n"
            "- **Warning** — borderline; review manually\n"
            "- **Fail** — exceeds a threshold; excluded by default"
        ),
    },
    "Analyze": {
        "title": "About Sequence Analysis",
        "content": (
            "This stage computes per-sequence statistics that help characterize "
            "your dataset before alignment.\n\n"
            "**Metrics computed:**\n"
            "- **GC Content** — the fraction of guanine + cytosine bases. "
            "Varies by organism (30-70% for most bacteria, ~40% for mammals). "
            "Unusual GC can indicate contamination.\n"
            "- **GC Skew** — (G-C)/(G+C) in sliding windows. Helps identify "
            "replication origins in bacterial genomes where the leading/lagging "
            "strands have different base compositions.\n"
            "- **Nucleotide composition** — raw counts of A, T/U, G, C per sequence.\n"
            "- **Motif scan** — finds occurrences of a regex pattern (e.g., start "
            "codons ATG, TATA boxes, restriction sites).\n\n"
            "**When is this useful?**\n"
            "- Screening for contamination (unexpected GC)\n"
            "- Identifying functional elements (motif scan)\n"
            "- Characterizing novel sequences before BLAST"
        ),
    },
    "Align": {
        "title": "About Sequence Alignment",
        "content": (
            "Alignment arranges sequences to identify regions of similarity that "
            "may reflect functional, structural, or evolutionary relationships.\n\n"
            "**Methods available:**\n"
            "- **MSA (Multiple Sequence Alignment)** — aligns all sequences "
            "simultaneously by padding shorter ones with gaps. Best for "
            "homologous sequences of similar length.\n"
            "- **Pairwise** — aligns consecutive pairs sequentially. Shows "
            "per-pair scores and identity percentages.\n\n"
            "**Identity percentage** tells you what fraction of aligned positions "
            "match between sequences. Higher identity suggests closer evolutionary "
            "relationship or stronger conservation.\n\n"
            "**Output formats:**\n"
            "- FASTA, CLUSTAL, NEXUS, PHYLIP — compatible with downstream tools "
            "like RAxML, MrBayes, MEGA, and IQ-TREE."
        ),
    },
    "Tree": {
        "title": "About Phylogenetic Trees",
        "content": (
            "A phylogenetic tree shows the inferred evolutionary relationships "
            "among your sequences based on their similarities and differences.\n\n"
            "**Methods:**\n"
            "- **Neighbor-Joining (NJ)** — fast, distance-based method. Good for "
            "exploratory analysis and large datasets.\n"
            "- **UPGMA** — assumes a molecular clock (equal rates of evolution). "
            "Use only when this assumption is reasonable.\n\n"
            "**Bootstrap support** tests how robust each branch is by resampling "
            "alignment columns and rebuilding the tree many times. Values >70% "
            "are generally considered well-supported.\n\n"
            "**Reading the tree:**\n"
            "- Branch length represents evolutionary distance (substitutions per site)\n"
            "- Closely grouped sequences share a recent common ancestor\n"
            "- The tree is unrooted unless UPGMA is used"
        ),
    },
    "Report": {
        "title": "About Report Generation",
        "content": (
            "This stage assembles all your results into a single downloadable "
            "document.\n\n"
            "**What's included:**\n"
            "- Auto-generated methods section (past tense, publication-ready)\n"
            "- QC summary table\n"
            "- Analysis metrics and charts\n"
            "- Alignment summary\n"
            "- Phylogenetic tree in Newick format\n"
            "- Any custom notes you add\n\n"
            "**Formats:**\n"
            "- **HTML** — best for viewing in a browser, includes interactive charts\n"
            "- **PDF** — best for printing or attaching to lab notebooks\n\n"
            "The methods section uses the exact parameters from your analysis, "
            "so you can paste it directly into a manuscript's Methods section."
        ),
    },
}

METRIC_TOOLTIPS: dict[str, str] = {
    "GC Content": (
        "Percentage of guanine + cytosine bases in the sequence. "
        "Ranges from 0-100%. Most organisms fall between 30-70%."
    ),
    "Ambiguity Rate": (
        "Percentage of ambiguous/degenerate bases (N, R, Y, K, M, S, W, etc.). "
        "High values suggest low-quality sequencing or unresolved bases."
    ),
    "Gap Fraction": (
        "Percentage of gap characters (-) in the sequence. "
        "High gap fractions in raw sequences indicate assembly issues."
    ),
    "Identity": (
        "Percentage of positions that match between aligned sequences. "
        "100% means identical; values >70% suggest homology."
    ),
    "Bootstrap Support": (
        "Percentage of bootstrap replicates that recovered a given branch. "
        "Values >70% indicate strong support; <50% is weak."
    ),
    "GC Skew": (
        "Calculated as (G-C)/(G+C) in a sliding window. "
        "Sign changes often mark replication origins in circular genomes."
    ),
}

QC_FLAG_EXPLANATIONS: dict[str, str] = {
    "High ambiguity": (
        "This sequence contains more ambiguous bases than the threshold allows. "
        "This typically means the sequencing quality was low or the base caller "
        "could not confidently assign nucleotides at many positions."
    ),
    "Too short": (
        "This sequence is shorter than the minimum length threshold. "
        "Very short sequences may not carry enough information for reliable "
        "alignment or phylogenetic inference."
    ),
    "Excessive gaps": (
        "More than half of this sequence consists of gap characters. "
        "This usually indicates an assembly artifact or a sequence that was "
        "already aligned and padded before upload."
    ),
    "Borderline ambiguity": (
        "The ambiguity rate is above 2% but below the fail threshold. "
        "The sequence may be usable but review it for potential quality issues."
    ),
    "Length outlier": (
        "This sequence's length is more than 2 standard deviations from the "
        "batch mean. It may be from a different gene or organism, or it may be "
        "a partial/truncated sequence."
    ),
    "GC outlier": (
        "This sequence's GC content is more than 2 standard deviations from "
        "the batch mean. It could indicate contamination from a different "
        "organism or a highly unusual genomic region."
    ),
}


def render_help_panel(stage: str) -> None:
    """Render a collapsible help panel for the given pipeline stage."""
    help_data = STAGE_HELP.get(stage)
    if not help_data:
        return

    # Check if user has dismissed help
    dismiss_key = f"help_dismissed_{stage}"
    if st.session_state.get(dismiss_key, False):
        if st.button("Show help", key=f"show_help_{stage}"):
            st.session_state[dismiss_key] = False
            st.rerun()
        return

    with st.expander(f"ℹ️ {help_data['title']}", expanded=False):
        st.markdown(help_data["content"])
        if st.button("Don't show again for this stage", key=f"dismiss_help_{stage}"):
            st.session_state[dismiss_key] = True
            st.rerun()
