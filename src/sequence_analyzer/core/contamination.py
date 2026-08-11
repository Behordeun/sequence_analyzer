"""Contamination detection via GC content and dinucleotide frequency profiling.

Compares each sequence's composition against expected organism profiles to flag
sequences that likely came from a different source. Works as a post-QC step:
sequences that pass basic QC (length, ambiguity, gaps) get checked for
cross-species contamination before downstream analysis.
"""

from __future__ import annotations

import functools
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# Path traversal: contamination.py lives at src/sequence_analyzer/core/
# parents[0]=core, [1]=sequence_analyzer, [2]=src, [3]=project_root
_PROFILES_PATH = Path(__file__).resolve().parents[3] / "data" / "organism_profiles.json"

_DINUCLEOTIDES = [
    "AA",
    "AT",
    "AG",
    "AC",
    "TA",
    "TT",
    "TG",
    "TC",
    "GA",
    "GT",
    "GG",
    "GC",
    "CA",
    "CT",
    "CG",
    "CC",
]


@dataclass(frozen=True)
class RiskThresholds:
    """Configurable thresholds for contamination risk classification.

    Attributes:
        gc_norm_medium: GC deviation (normalized to range width) above which
            the sequence gets a medium GC score. Below this is low.
        gc_norm_high: GC deviation (normalized) above which the sequence
            gets a high GC score.
        dinuc_low_max: Dinucleotide distance below this is considered normal.
            Calibrated against random 200bp fragments from matching genomes.
        dinuc_medium_max: Dinucleotide distance below this is slightly atypical.
            Above this strongly indicates a different source organism.
    """

    gc_norm_medium: float = 0.5
    gc_norm_high: float = 0.5
    dinuc_low_max: float = 0.12
    dinuc_medium_max: float = 0.18


DEFAULT_THRESHOLDS = RiskThresholds()


@functools.lru_cache(maxsize=4)
def _load_profiles_cached(path_str: str) -> dict[str, Any]:
    """Cache-backed loader keyed by path string (hashable for lru_cache)."""
    target = Path(path_str)

    if not target.exists():
        raise FileNotFoundError(f"Organism profiles not found at {target}")

    try:
        with open(target, encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Organism profiles file is not valid JSON: {exc}") from exc

    if "organisms" not in data:
        raise ValueError("Organism profiles file missing 'organisms' key.")

    return data["organisms"]


def load_organism_profiles(path: Path | None = None) -> dict[str, Any]:
    """Load organism profiles from disk (cached after first read).

    Args:
        path: Override path for testing. Defaults to bundled organism_profiles.json.

    Returns:
        Dictionary keyed by organism slug with profile data.

    Raises:
        FileNotFoundError: If the profiles file does not exist.
        ValueError: If the file is malformed or missing the 'organisms' key.
    """
    target = path or _PROFILES_PATH
    return _load_profiles_cached(str(target))


def list_available_organisms(path: Path | None = None) -> list[str]:
    """Return sorted list of available organism profile slugs.

    Args:
        path: Override path for testing.

    Returns:
        Sorted list of organism identifiers (e.g., ['escherichia_coli', 'general', ...]).
    """
    profiles = load_organism_profiles(path)
    return sorted(profiles.keys())


def compute_dinucleotide_frequencies(sequence: str) -> dict[str, float]:
    """Compute observed dinucleotide frequencies for a sequence.

    Counts overlapping dinucleotides and normalizes to relative frequencies
    that sum to 1.0.

    Args:
        sequence: Uppercase nucleotide string (gaps and ambiguous chars excluded).

    Returns:
        Dict mapping each of the 16 canonical dinucleotides to its frequency.
        Returns uniform distribution if sequence is too short (<2 bases).
    """
    clean = "".join(ch for ch in sequence.upper() if ch in "ATGC")

    if len(clean) < 2:
        return {di: 1.0 / 16.0 for di in _DINUCLEOTIDES}

    counts = {di: 0 for di in _DINUCLEOTIDES}
    for i in range(len(clean) - 1):
        pair = clean[i : i + 2]
        if pair in counts:
            counts[pair] += 1

    total = sum(counts.values())
    if total == 0:
        return {di: 1.0 / 16.0 for di in _DINUCLEOTIDES}

    return {di: count / total for di, count in counts.items()}


def compute_gc_content(sequence: str) -> float:
    """Compute GC percentage from a sequence string.

    Only counts unambiguous ATGC bases in the denominator.

    Args:
        sequence: Uppercase nucleotide string.

    Returns:
        GC percentage (0-100). Returns 0.0 for sequences with no ATGC bases.
    """
    clean = "".join(ch for ch in sequence.upper() if ch in "ATGC")
    if not clean:
        return 0.0
    gc = clean.count("G") + clean.count("C")
    return (gc / len(clean)) * 100.0


def _euclidean_distance(observed: dict[str, float], expected: dict[str, float]) -> float:
    """Euclidean distance between two dinucleotide frequency vectors."""
    obs_vec = np.array([observed.get(di, 0.0) for di in _DINUCLEOTIDES])
    exp_vec = np.array([expected.get(di, 0.0) for di in _DINUCLEOTIDES])
    return float(np.sqrt(np.sum((obs_vec - exp_vec) ** 2)))


def _classify_risk(
    gc_deviation: float,
    dinuc_distance: float,
    gc_range_width: float,
    thresholds: RiskThresholds = DEFAULT_THRESHOLDS,
) -> tuple[str, str]:
    """Classify contamination risk based on GC deviation and dinucleotide distance.

    Uses a two-axis scoring system:
    - GC deviation as a fraction of the expected range width
    - Dinucleotide distance against empirically-calibrated thresholds

    Args:
        gc_deviation: How far the observed GC is from the nearest edge of the expected range.
            Zero means within range.
        dinuc_distance: Euclidean distance between observed and expected dinucleotide vectors.
        gc_range_width: Width of the expected GC range (max - min).
        thresholds: Configurable threshold values for risk classification.

    Returns:
        Tuple of (risk_level, reasoning) where risk_level is one of:
        "Low", "Medium", "High".
    """
    reasons: list[str] = []

    gc_score = 0
    dinuc_score = 0

    # Guard against zero-width ranges (would make normalization meaningless)
    effective_width = max(gc_range_width, 1.0)
    gc_norm = gc_deviation / effective_width

    if gc_deviation == 0.0:
        gc_score = 0
    elif gc_norm < thresholds.gc_norm_medium:
        gc_score = 1
        reasons.append(f"GC {gc_deviation:.1f}% outside expected range")
    else:
        gc_score = 2
        reasons.append(f"GC {gc_deviation:.1f}% far outside expected range")

    if dinuc_distance < thresholds.dinuc_low_max:
        dinuc_score = 0
    elif dinuc_distance < thresholds.dinuc_medium_max:
        dinuc_score = 1
        reasons.append(f"Dinucleotide profile slightly atypical (dist={dinuc_distance:.3f})")
    else:
        dinuc_score = 2
        reasons.append(f"Dinucleotide profile divergent (dist={dinuc_distance:.3f})")

    combined = gc_score + dinuc_score

    if combined == 0:
        return "Low", "Composition matches expected profile"
    elif combined <= 2:
        return "Medium", "; ".join(reasons)
    else:
        return "High", "; ".join(reasons)


def detect_contamination(
    sequences: list[SeqRecord],
    organism: str = "general",
    profiles_path: Path | None = None,
    thresholds: RiskThresholds | None = None,
) -> pd.DataFrame:
    """Screen sequences for potential contamination based on composition profiling.

    Compares each sequence's GC content and dinucleotide frequencies against
    the expected profile for the declared organism. Sequences that deviate
    significantly are flagged as potential contaminants.

    Args:
        sequences: List of SeqRecords to screen.
        organism: Organism profile key (e.g., 'escherichia_coli', 'sars_cov_2').
            Defaults to 'general' which uses permissive thresholds.
        profiles_path: Override path to profiles JSON for testing.
        thresholds: Optional RiskThresholds to override default classification bands.

    Returns:
        DataFrame with columns:
        - ID: Sequence identifier
        - GC_Observed: Measured GC percentage
        - GC_Expected: Mean expected GC for the organism
        - GC_Deviation: Distance from nearest edge of expected range (0 if within range)
        - Dinuc_Distance: Euclidean distance from expected dinucleotide profile
        - Contamination_Risk: "Low", "Medium", or "High"
        - Risk_Reason: Explanation of the risk classification

    Raises:
        ValueError: If sequences list is empty or organism profile not found.
    """
    if not sequences:
        raise ValueError("At least one sequence is required for contamination detection.")

    risk_thresholds = thresholds or DEFAULT_THRESHOLDS
    profiles = load_organism_profiles(profiles_path)

    if organism not in profiles:
        available = sorted(profiles.keys())
        raise ValueError(f"Unknown organism '{organism}'. Available: {available}")

    profile = profiles[organism]
    gc_min, gc_max = profile["gc_range"]
    gc_mean = profile["gc_mean"]
    expected_dinuc = profile["dinucleotide_frequencies"]
    gc_range_width = gc_max - gc_min

    rows: list[dict[str, object]] = []

    for record in sequences:
        seq_str = str(record.seq).upper()

        gc_observed = compute_gc_content(seq_str)
        observed_dinuc = compute_dinucleotide_frequencies(seq_str)
        dinuc_distance = _euclidean_distance(observed_dinuc, expected_dinuc)

        # GC deviation: 0 if within range, distance to nearest boundary otherwise
        if gc_observed < gc_min:
            gc_deviation = gc_min - gc_observed
        elif gc_observed > gc_max:
            gc_deviation = gc_observed - gc_max
        else:
            gc_deviation = 0.0

        risk_level, reason = _classify_risk(
            gc_deviation, dinuc_distance, gc_range_width, risk_thresholds
        )

        rows.append(
            {
                "ID": record.id or "unnamed",
                "GC_Observed": round(gc_observed, 2),
                "GC_Expected": gc_mean,
                "GC_Deviation": round(gc_deviation, 2),
                "Dinuc_Distance": round(dinuc_distance, 4),
                "Contamination_Risk": risk_level,
                "Risk_Reason": reason,
            }
        )

    return pd.DataFrame(rows)
