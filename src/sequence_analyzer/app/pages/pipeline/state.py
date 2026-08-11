"""Pipeline state management: dataclasses and session state helpers.

All pipeline stages read from and write to a single PipelineState object
stored in st.session_state["pipeline"]. This keeps data flowing forward
without re-uploads or re-configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

    from sequence_analyzer.models.sequences import (
        AlignmentResult,
        PhylogeneticResult,
        VariantResult,
    )

STAGES = ["Ingest", "QC", "Analyze", "Align", "Tree", "Report"]


@dataclass
class PipelineConfig:
    """User-configurable parameters applied across all pipeline stages."""

    organism_mode: str = "general"
    seq_type: str = "auto"
    alignment_method: str = "msa"
    tree_method: str = "nj"
    gc_window_size: int = 100
    motif_pattern: str = "ATG"
    bootstrap_replicates: int = 50
    qc_ambiguity_threshold: float = 5.0
    qc_min_length: int = 10
    qc_max_gap_fraction: float = 50.0
    contamination_profile: str = "general"


@dataclass
class PipelineState:
    """Central state container for the guided pipeline.

    Each stage writes its output here so downstream stages can consume it.
    """

    stage: str = "Ingest"
    completed_stages: list[str] = field(default_factory=list)
    config: PipelineConfig = field(default_factory=PipelineConfig)

    # Ingest output
    sequences: list[SeqRecord] = field(default_factory=list)
    ingest_summary: str = ""

    # QC output
    qc_results: pd.DataFrame | None = None
    valid_sequences: list[SeqRecord] = field(default_factory=list)

    # Analysis output
    analysis_results: pd.DataFrame | None = None

    # Alignment output
    alignment_result: AlignmentResult | None = None

    # Variant calling output
    variant_result: VariantResult | None = None

    # Tree output
    tree_result: PhylogeneticResult | None = None

    # Report config
    report_title: str = ""
    report_notes: str = ""

    # Figures stored as Plotly HTML strings keyed by name
    figures: dict[str, str] = field(default_factory=dict)


def get_pipeline_state() -> PipelineState:
    """Retrieve or initialize the pipeline state from Streamlit session."""
    import streamlit as st

    if "pipeline" not in st.session_state:
        st.session_state["pipeline"] = PipelineState()
    return st.session_state["pipeline"]


def update_stage(stage: str) -> None:
    """Set the current pipeline stage."""
    state = get_pipeline_state()
    state.stage = stage


def mark_complete(stage: str) -> None:
    """Mark a stage as completed so downstream stages are unlocked."""
    state = get_pipeline_state()
    if stage not in state.completed_stages:
        state.completed_stages.append(stage)


def is_stage_complete(stage: str) -> bool:
    """Check if a stage has been completed."""
    state = get_pipeline_state()
    return stage in state.completed_stages


def reset_pipeline() -> None:
    """Reset the entire pipeline state to initial values."""
    import streamlit as st

    st.session_state["pipeline"] = PipelineState()


def reset_config() -> None:
    """Reset just the configuration to default values."""
    state = get_pipeline_state()
    state.config = PipelineConfig()


def apply_organism_mode(slug: str) -> None:
    """Apply an organism mode's settings to the current pipeline config.

    Loads the organism profile and overwrites the relevant PipelineConfig
    fields. Does not touch alignment or tree settings.

    Args:
        slug: Organism mode identifier (e.g., 'escherichia_coli').
    """
    from sequence_analyzer.core.organisms import get_organism_mode

    mode = get_organism_mode(slug)
    state = get_pipeline_state()
    config = state.config

    config.organism_mode = slug
    config.seq_type = mode.seq_type
    config.qc_ambiguity_threshold = mode.qc_ambiguity_threshold
    config.qc_min_length = mode.qc_min_length
    config.qc_max_gap_fraction = mode.qc_max_gap_fraction
    config.contamination_profile = mode.contamination_profile

    if mode.motif_patterns:
        config.motif_pattern = mode.motif_patterns[0]
