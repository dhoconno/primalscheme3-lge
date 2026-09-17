"""CLI validation for the opt-in legacy gap-completion workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any


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
