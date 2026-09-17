from __future__ import annotations

import hashlib
import json
from pathlib import Path

from primalscheme3.core.config import Config, MappingType, TerminalGapPolicy
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.gap_expansion import GapExpansionOptions
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate


MSA = Path("tests/core/test_mismatch.fasta").resolve()
SCIENTIFIC = ("primer.bed", "amplicon.bed", "primertrim.amplicon.bed", "reference.fasta")


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
        for name in (*SCIENTIFIC, "config.json")
    }


def test_bounded_gap_expansion_workflow_preserves_parent_and_records_options(tmp_path):
    parent = tmp_path / "parent"
    panelcreate(
        msa=[MSA], output_dir=parent, config=_config(), pm=ProgressManager(),
        force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
    )
    before = _hashes(parent)
    output = tmp_path / "expanded"
    panelcreate(
        msa=[MSA], output_dir=output, config=_config(), pm=ProgressManager(),
        force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
        gap_completion_parent=parent,
        gap_expansion_options=GapExpansionOptions(
            mode="bounded", max_anchors_per_msa=100, max_pairs_per_msa=50
        ),
    )
    assert _hashes(parent) == before
    report = json.loads((output / "gap-completion.json").read_text())
    provenance = json.loads((output / "panel-provenance.json").read_text())
    assert report["followup"]["poolCount"] == 3
    assert report["parent"]["unchangedAfterRun"] is True
    assert report["gapExpansion"]["options"] == {
        "mode": "bounded", "maxAnchorsPerMsa": 100, "maxPairsPerMsa": 50
    }
    assert report["gapExpansion"]["perMsa"]
    assert provenance["resolvedOptions"]["gap_expansion_mode"] == "bounded"
    assert provenance["resolvedOptions"]["gap_expansion_max_anchors_per_msa"] == 100
    assert provenance["resolvedOptions"]["gap_expansion_max_pairs_per_msa"] == 50


def test_gap_expansion_off_matches_ordinary_gap_completion_scientific_outputs(tmp_path):
    parent = tmp_path / "parent"
    panelcreate(
        msa=[MSA], output_dir=parent, config=_config(), pm=ProgressManager(),
        force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
    )
    ordinary = tmp_path / "ordinary"
    off = tmp_path / "off"
    for output, options in ((ordinary, None), (off, GapExpansionOptions(mode="off"))):
        panelcreate(
            msa=[MSA], output_dir=output, config=_config(), pm=ProgressManager(),
            force=False, mode=PanelRunModes.EQUAL, offline_plots=False,
            gap_completion_parent=parent, gap_expansion_options=options,
        )
    assert _hashes(ordinary) == _hashes(off)
    report = json.loads((off / "gap-completion.json").read_text())
    assert report["gapExpansion"]["options"]["mode"] == "off"
