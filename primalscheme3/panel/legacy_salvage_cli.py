"""Native CLI parsing and capability validation for legacy dimer salvage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .legacy_salvage import LegacySalvageOptions


LEGACY_SALVAGE_OPTION_NAMES = frozenset(
    {
        "legacy_salvage",
        "legacy_salvage_thresholds",
        "legacy_salvage_floor",
        "legacy_salvage_max_edges_per_pool",
        "legacy_salvage_max_incident_species_per_pool",
        "legacy_salvage_min_reference_gain",
        "legacy_salvage_max_candidate_evaluations",
    }
)


def resolve_legacy_salvage_options(
    *,
    selection_algorithm: str,
    mode: Any,
    mapping: Any,
    region_bedfile: Path | None,
    input_bedfile: Path | None,
    config_input_bedfile: Path | None,
    strict_cutoff: float,
    explicit: set[str],
    legacy_salvage: str | None,
    legacy_salvage_thresholds: list[float] | None,
    legacy_salvage_floor: float | None,
    legacy_salvage_max_edges_per_pool: int | None,
    legacy_salvage_max_incident_species_per_pool: int | None,
    legacy_salvage_min_reference_gain: int | None,
    legacy_salvage_max_candidate_evaluations: int | None,
) -> LegacySalvageOptions | None:
    """Resolve the separate legacy salvage namespace, preserving old defaults."""
    values = {
        "legacy_salvage_thresholds": legacy_salvage_thresholds,
        "legacy_salvage_floor": legacy_salvage_floor,
        "legacy_salvage_max_edges_per_pool": legacy_salvage_max_edges_per_pool,
        "legacy_salvage_max_incident_species_per_pool": legacy_salvage_max_incident_species_per_pool,
        "legacy_salvage_min_reference_gain": legacy_salvage_min_reference_gain,
        "legacy_salvage_max_candidate_evaluations": legacy_salvage_max_candidate_evaluations,
    }
    requested_controls = sorted(name for name in LEGACY_SALVAGE_OPTION_NAMES if name in explicit and name != "legacy_salvage")
    requested = legacy_salvage is not None or bool(requested_controls)
    if not requested:
        return None
    if legacy_salvage != "bounded":
        if legacy_salvage not in (None, "off"):
            raise ValueError("--legacy-salvage must be off or bounded")
        if requested_controls:
            raise ValueError("legacy salvage controls require --legacy-salvage bounded")
        return None
    if selection_algorithm != "legacy":
        raise ValueError("--legacy-salvage bounded requires --selection-algorithm legacy")
    if str(getattr(mode, "value", mode)) != "equal":
        raise ValueError("legacy salvage requires --mode equal")
    if str(getattr(mapping, "value", mapping)) != "first":
        raise ValueError("legacy salvage requires --mapping first")
    if region_bedfile is not None:
        raise ValueError("legacy salvage does not support region inputs")
    if input_bedfile is not None or config_input_bedfile is not None:
        raise ValueError("legacy salvage does not support imported primer pairs")

    kwargs: dict[str, Any] = {"mode": "bounded"}
    if legacy_salvage_thresholds is not None:
        kwargs["thresholds"] = tuple(legacy_salvage_thresholds)
    if legacy_salvage_floor is not None:
        kwargs["floor"] = legacy_salvage_floor
    if legacy_salvage_max_edges_per_pool is not None:
        kwargs["max_edges_per_pool"] = legacy_salvage_max_edges_per_pool
    if legacy_salvage_max_incident_species_per_pool is not None:
        kwargs["max_incident_species_per_pool"] = legacy_salvage_max_incident_species_per_pool
    if legacy_salvage_min_reference_gain is not None:
        kwargs["min_reference_gain"] = legacy_salvage_min_reference_gain
    if legacy_salvage_max_candidate_evaluations is not None:
        kwargs["max_candidate_evaluations"] = legacy_salvage_max_candidate_evaluations
    options = LegacySalvageOptions(**kwargs)
    options.validate_for_strict_cutoff(float(strict_cutoff))
    return options
