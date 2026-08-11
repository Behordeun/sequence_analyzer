"""Organism-specific mode management.

Loads pre-configured profiles that set QC thresholds, contamination
detection parameters, motif patterns, and sequence type defaults for
common research organisms. Eliminates repetitive manual configuration.
"""

from __future__ import annotations

import functools
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_ORGANISMS_DIR = Path(__file__).resolve().parents[3] / "data" / "organisms"


@dataclass(frozen=True)
class OrganismMode:
    """Pre-configured analysis parameters for a specific organism.

    Each field maps directly to a PipelineConfig setting. Profiles override
    values they specify; unspecified fields use the dataclass defaults.
    """

    slug: str
    display_name: str
    description: str = ""

    # Sequence type
    seq_type: str = "DNA"

    # QC thresholds
    qc_ambiguity_threshold: float = 5.0
    qc_min_length: int = 10
    qc_max_gap_fraction: float = 50.0

    # Contamination detection profile key (matches organism_profiles.json)
    contamination_profile: str = "general"

    # Motif patterns relevant to this organism
    motif_patterns: list[str] = field(default_factory=list)

    # Suggested reference accession (for variant calling guidance)
    reference_accession: str = ""
    reference_description: str = ""


def _parse_organism_mode(slug: str, data: dict[str, Any]) -> OrganismMode:
    """Build an OrganismMode from a parsed JSON dict."""
    return OrganismMode(
        slug=slug,
        display_name=data.get("display_name", slug),
        description=data.get("description", ""),
        seq_type=data.get("seq_type", "DNA"),
        qc_ambiguity_threshold=data.get("qc_ambiguity_threshold", 5.0),
        qc_min_length=data.get("qc_min_length", 10),
        qc_max_gap_fraction=data.get("qc_max_gap_fraction", 50.0),
        contamination_profile=data.get("contamination_profile", "general"),
        motif_patterns=data.get("motif_patterns", []),
        reference_accession=data.get("reference_accession", ""),
        reference_description=data.get("reference_description", ""),
    )


@functools.lru_cache(maxsize=4)
def _load_modes_cached(resolved_dir: str) -> dict[str, OrganismMode]:
    """Cache-backed loader keyed by resolved path string."""
    directory = Path(resolved_dir)
    modes: dict[str, OrganismMode] = {}

    if not directory.exists():
        return modes

    for path in sorted(directory.glob("*.json")):
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
            slug = path.stem
            modes[slug] = _parse_organism_mode(slug, data)
        except (json.JSONDecodeError, KeyError) as exc:
            logger.warning("Skipping malformed organism profile %s: %s", path.name, exc)

    return modes


def load_organism_modes(organisms_dir: Path | None = None) -> dict[str, OrganismMode]:
    """Load all available organism modes.

    Args:
        organisms_dir: Override directory for testing. Defaults to data/organisms/.

    Returns:
        Dictionary mapping organism slug to OrganismMode.
    """
    target = (organisms_dir or _ORGANISMS_DIR).resolve()
    return _load_modes_cached(str(target))


def list_organism_modes(organisms_dir: Path | None = None) -> list[str]:
    """Return sorted list of available organism mode slugs.

    Args:
        organisms_dir: Override directory for testing.

    Returns:
        Sorted list of mode identifiers.
    """
    modes = load_organism_modes(organisms_dir)
    return sorted(modes.keys())


def get_organism_mode(slug: str, organisms_dir: Path | None = None) -> OrganismMode:
    """Load a specific organism mode by slug.

    Args:
        slug: Organism identifier (e.g., 'sars_cov_2', 'escherichia_coli').
        organisms_dir: Override directory for testing.

    Returns:
        The matching OrganismMode.

    Raises:
        ValueError: If the slug is not found.
    """
    modes = load_organism_modes(organisms_dir)
    if slug not in modes:
        available = sorted(modes.keys())
        raise ValueError(f"Unknown organism mode '{slug}'. Available: {available}")
    return modes[slug]
