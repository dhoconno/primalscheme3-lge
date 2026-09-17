from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from primalscheme3.core.config import Config, MappingType, TerminalGapPolicy
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate


MSA = Path("tests/core/test_mismatch.fasta").resolve()


def _config() -> Config:
    config = Config(
        selection_algorithm="legacy",
        mapping=MappingType.FIRST,
        terminal_gap_policy=TerminalGapPolicy.LEGACY,
        n_pools=3,
    )
    config.use_matchdb = False
    return config


def _hashes(path: Path) -> dict[str, str]:
    return {
        name: hashlib.sha256((path / name).read_bytes()).hexdigest()
        for name in ("primer.bed", "amplicon.bed", "primertrim.amplicon.bed", "reference.fasta", "config.json")
    }


def test_gap_completion_workflow_preserves_parent_and_records_followup(tmp_path):
    parent = tmp_path / "parent"
    panelcreate(
        msa=[MSA], output_dir=parent, config=_config(), pm=ProgressManager(),
        force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
    )
    before = _hashes(parent)
    output = tmp_path / "followup"
    panelcreate(
        msa=[MSA], output_dir=output, config=_config(), pm=ProgressManager(),
        force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
        gap_completion_parent=parent,
    )
    assert _hashes(parent) == before
    report = json.loads((output / "gap-completion.json").read_text())
    coverage = json.loads((output / "gap-completion-coverage.json").read_text())
    assert report["followup"]["poolCount"] == 3
    assert report["parent"]["verification"]["referenceAndPrimerBed"] is True
    assert "0" in coverage["primary"] and "0" in coverage["combined"]
    assert (output / "panel-provenance.json").is_file()


def test_gap_completion_rejects_corrupt_parent_with_failure_receipt(tmp_path):
    parent = tmp_path / "parent"
    panelcreate(
        msa=[MSA], output_dir=parent, config=_config(), pm=ProgressManager(),
        force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
    )
    (parent / "reference.fasta").write_text((parent / "reference.fasta").read_text() + "A\n")
    output = tmp_path / "failed"
    with pytest.raises(ValueError, match="parent reference"):
        panelcreate(
            msa=[MSA], output_dir=output, config=_config(), pm=ProgressManager(),
            force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
            gap_completion_parent=parent,
        )
    receipt = output / "panel-provenance.json"
    assert receipt.is_file()
    assert json.loads(receipt.read_text())["status"] == "failure"
