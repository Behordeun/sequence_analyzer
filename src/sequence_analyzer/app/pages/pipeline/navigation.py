"""Pipeline navigation: progress bar and stage navigation controls.

Renders a horizontal progress indicator and Next/Back buttons.
Handles stage transitions and disabling of incomplete stages.
"""

from __future__ import annotations

import streamlit as st

from sequence_analyzer.app.pages.pipeline.state import (
    STAGES,
    get_pipeline_state,
    is_stage_complete,
    reset_config,
    update_stage,
)


def render_progress_bar() -> None:
    """Render a horizontal progress indicator showing pipeline stages."""
    state = get_pipeline_state()
    current_idx = STAGES.index(state.stage) if state.stage in STAGES else 0

    cols = st.columns(len(STAGES))
    for i, (col, stage_name) in enumerate(zip(cols, STAGES, strict=True)):
        if stage_name in state.completed_stages:
            icon = "✅"
        elif i == current_idx:
            icon = "🔵"
        else:
            icon = "⚪"

        with col:
            st.markdown(
                f"<div style='text-align:center; font-size:12px;'>{icon}<br>{stage_name}</div>",
                unsafe_allow_html=True,
            )


def render_stage_nav() -> str | None:
    """Render Next/Back buttons and return the target stage if navigating.

    Returns:
        The stage name to navigate to, or None if no navigation occurred.
    """
    state = get_pipeline_state()
    current_idx = STAGES.index(state.stage) if state.stage in STAGES else 0

    col_back, _, col_next = st.columns([1, 4, 1])

    target: str | None = None

    with col_back:
        if current_idx > 0 and st.button("← Back"):
            target = STAGES[current_idx - 1]

    with col_next:
        if current_idx < len(STAGES) - 1:
            next_stage = STAGES[current_idx + 1]
            # Next is enabled only if current stage is complete
            can_proceed = is_stage_complete(state.stage)

            if st.button("Next →", disabled=not can_proceed):
                target = next_stage

    if target:
        update_stage(target)

    return target


def render_config_sidebar() -> None:
    """Render the pipeline configuration sidebar."""
    state = get_pipeline_state()
    config = state.config

    with st.sidebar:
        st.markdown("### Pipeline Config")

        config.seq_type = st.selectbox(
            "Sequence Type",
            ["auto", "DNA", "RNA", "Protein"],
            index=["auto", "DNA", "RNA", "Protein"].index(config.seq_type),
            key="pipeline_seq_type",
        )

        config.alignment_method = st.selectbox(
            "Alignment Method",
            ["msa", "pairwise"],
            index=["msa", "pairwise"].index(config.alignment_method),
            key="pipeline_align_method",
        )

        config.tree_method = st.selectbox(
            "Tree Method",
            ["nj", "upgma"],
            index=["nj", "upgma"].index(config.tree_method),
            key="pipeline_tree_method",
        )

        config.gc_window_size = st.slider(
            "GC Skew Window",
            min_value=10,
            max_value=500,
            value=config.gc_window_size,
            step=10,
            key="pipeline_gc_window",
        )

        config.motif_pattern = st.text_input(
            "Motif Pattern",
            value=config.motif_pattern,
            key="pipeline_motif",
        )

        config.bootstrap_replicates = st.slider(
            "Bootstrap Replicates",
            min_value=10,
            max_value=200,
            value=config.bootstrap_replicates,
            step=10,
            key="pipeline_bootstrap",
        )

        st.markdown("---")
        st.markdown("**QC Thresholds**")

        config.qc_ambiguity_threshold = st.number_input(
            "Max Ambiguity %",
            min_value=0.1,
            max_value=50.0,
            value=config.qc_ambiguity_threshold,
            step=0.5,
            key="pipeline_qc_ambig",
        )

        config.qc_min_length = st.number_input(
            "Min Length (bp)",
            min_value=1,
            max_value=1000,
            value=config.qc_min_length,
            step=1,
            key="pipeline_qc_min_len",
        )

        config.qc_max_gap_fraction = st.number_input(
            "Max Gap %",
            min_value=1.0,
            max_value=100.0,
            value=config.qc_max_gap_fraction,
            step=1.0,
            key="pipeline_qc_gaps",
        )

        st.markdown("---")
        if st.button("Reset to Defaults", key="pipeline_reset"):
            reset_config()
            st.rerun()
