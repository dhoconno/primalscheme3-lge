"""CLI validation for the opt-in legacy gap-completion workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any


GAP_EXPANSION_OPTION_NAMES = {
    "gap_expansion",
    "gap_expansion_max_anchors_per_msa",
    "gap_expansion_max_pairs_per_msa",
}


def resolve_gap_expansion_options(
    *,
    mode: str | None,
    parent: Path | None,
    max_anchors_per_msa: int | None,
    max_pairs_per_msa: int | None,
):
    """Resolve opt-in gap expansion without changing ordinary gap completion."""
    advanced = max_anchors_per_msa is not None or max_pairs_per_msa is not None
    normalized = None if mode is None else str(getattr(mode, "value", mode))
    if normalized not in (None, "off", "bounded"):
        raise ValueError("--gap-expansion must be off or bounded")
    if normalized in (None, "off") and advanced:
        raise ValueError(
            "gap expansion advanced bounds require --gap-expansion bounded"
        )
    if normalized != "bounded":
        return None
    if parent is None:
        raise ValueError("--gap-expansion bounded requires --gap-completion-parent")
    try:
        from primalscheme3.panel.gap_expansion import GapExpansionOptions
    except ImportError as error:  # pragma: no cover - core feature must be installed
        raise ValueError("gap expansion is unavailable in this PrimalScheme build") from error
    return GapExpansionOptions(
        mode="bounded",
        max_anchors_per_msa=(2000 if max_anchors_per_msa is None else max_anchors_per_msa),
        max_pairs_per_msa=(1000 if max_pairs_per_msa is None else max_pairs_per_msa),
    )


def resolve_gap_completion_parent(
    *,
    parent: Path | None,
    selection_algorithm: str,
    mode: Any,
    mapping: Any,
    terminal_gap_policy: Any,
    circular: bool,
    region_bedfile: Path | None,
    input_bedfile: Path | None,
    legacy_salvage: str | None,
) -> Path | None:
    """Validate and resolve the parent panel only when explicitly requested."""
    if parent is None:
        return None
    if selection_algorithm != "legacy":
        raise ValueError("--gap-completion-parent requires --selection-algorithm legacy")
    if str(getattr(mode, "value", mode)) != "equal":
        raise ValueError("gap completion requires --mode equal")
    if str(getattr(mapping, "value", mapping)) != "first":
        raise ValueError("gap completion requires --mapping first")
    if str(getattr(terminal_gap_policy, "value", terminal_gap_policy)) != "legacy":
        raise ValueError("gap completion requires --terminal-gap-policy legacy")
    if circular:
        raise ValueError("gap completion does not support circular references")
    if region_bedfile is not None:
        raise ValueError("gap completion does not support region inputs")
    if input_bedfile is not None:
        raise ValueError("gap completion does not support imported primer pairs")
    if legacy_salvage not in (None, "off"):
        raise ValueError("gap completion cannot be combined with --legacy-salvage bounded")
    resolved = parent.expanduser().resolve()
    if not resolved.is_dir():
        raise ValueError(f"gap completion parent must be an existing directory: {parent}")
    return resolved
