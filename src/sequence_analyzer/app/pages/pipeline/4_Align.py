"""Pipeline Stage 4: Align — pairwise or multiple sequence alignment."""

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
from sequence_analyzer.core.alignment import align_msa, align_pairwise
from sequence_analyzer.core.variants import call_variants, summarize_variants

st.set_page_config(layout="wide", page_title="Pipeline: Align")
apply_styles()
hide_streamlit_chrome()
render_config_sidebar()

state = get_pipeline_state()
state.stage = "Align"

render_progress_bar()
st.title("🔗 Stage 4: Sequence Alignment")


render_help_panel("Align")

sequences = state.valid_sequences or state.sequences
if len(sequences) < 2:
    st.warning("At least 2 sequences are required for alignment. Go back and load more.")
    render_stage_nav()
    render_footer()
    st.stop()

config = state.config

if st.button("Run Alignment"):
    with st.spinner("Aligning sequences..."):
        if config.alignment_method == "pairwise":
            result = align_pairwise(sequences)
        else:
            result = align_msa(sequences, seq_type=config.seq_type)

        state.alignment_result = result
        mark_complete("Align")
        st.success("Alignment complete!")

# --- Display results ---
if state.alignment_result:
    result = state.alignment_result

    st.metric("Average Identity", f"{result.avg_identity:.1f}%")

    if result.pairs:
        st.subheader("Pairwise Results")
        for pair in result.pairs:
            with st.expander(f"{pair.seq_a_id} vs {pair.seq_b_id} — {pair.identity:.1f}%"):
                st.code(f"{pair.aligned_a}\n{pair.aligned_b}", language="text")
    elif result.format_exports:
        st.subheader("MSA Output")
        fmt = st.selectbox(
            "View format",
            ["clustal", "fasta", "phylip", "nexus"],
            key="align_view_fmt",
        )
        content = result.format_exports.get(fmt, "")
        if content:
            st.code(content, language="text")

        # Download
        dl_fmt = st.selectbox(
            "Download format",
            ["fasta", "clustal", "nexus", "phylip"],
            key="align_dl_fmt",
        )
        dl_content = result.format_exports.get(dl_fmt, "")
        if dl_content:
            st.download_button(
                f"Download ({dl_fmt})",
                dl_content,
                file_name=f"alignment.{dl_fmt}",
            )

st.markdown("---")

# --- Variant Calling ---
st.subheader("Variant Calling")

if len(sequences) >= 2:
    ref_options = {r.id: r for r in sequences}
    selected_ref_id = st.selectbox(
        "Reference sequence",
        options=list(ref_options.keys()),
        index=0,
        help="Select which sequence to use as the reference. All others become samples.",
    )

    if st.button("Call Variants"):
        reference = ref_options[selected_ref_id]
        samples = [s for s in sequences if s.id != selected_ref_id]

        if not samples:
            st.warning("Need at least one sample sequence besides the reference.")
        else:
            with st.spinner("Calling variants..."):
                try:
                    variant_result = call_variants(reference, samples)
                    state.variant_result = variant_result
                    st.success(
                        f"Found {variant_result.summary['total']} variant(s) "
                        f"across {len(samples)} sample(s)."
                    )
                except Exception as e:
                    st.error(f"Variant calling failed: {e}")

    # Display variant results
    if state.variant_result:
        vr = state.variant_result

        # Summary metrics
        v_col1, v_col2, v_col3, v_col4 = st.columns(4)
        v_col1.metric("Total", vr.summary.get("total", 0))
        v_col2.metric("SNPs", vr.summary.get("SNP", 0))
        v_col3.metric("Insertions", vr.summary.get("insertion", 0))
        v_col4.metric("Deletions", vr.summary.get("deletion", 0))

        # Variant table
        variant_df = summarize_variants(vr)
        if not variant_df.empty:
            st.dataframe(variant_df, use_container_width=True)

            # Download
            st.download_button(
                "Download Variant Table (CSV)",
                variant_df.to_csv(index=False),
                file_name="variants.csv",
                mime="text/csv",
            )
        else:
            st.info("No variants detected — samples are identical to the reference.")

st.markdown("---")
render_stage_nav()
render_footer()
