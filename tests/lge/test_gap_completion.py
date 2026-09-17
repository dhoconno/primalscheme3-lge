from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import dnaio
import pytest
from primalschemers import FKmer, RKmer

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config, MappingType, TerminalGapPolicy
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.gap_completion import (
    FollowupSelection,
    _parent_reference,
    gap_candidate_id,
    select_followup_candidates,
    trimmed_coverage,
)
from primalscheme3.panel.gap_expansion import (
    GapExpansionOptions,
    _anchor_tasks,
    _spatial_order,
)
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate


def pair(start: int, end: int, forward: str = "A" * 20, reverse: str = "C" * 20):
    return PrimerPair(
        FKmer([forward.encode()], start + len(forward)),
        RKmer([reverse.encode()], end - len(reverse)),
        0,
    )


def target(length: int = 200):
    return SimpleNamespace(
        msa_index=0,
        _chrom_name="target",
        _mapping_array=list(range(length)),
    )


def test_trimmed_coverage_uses_ungapped_half_open_union():
    msa = target(120)
    incumbent = pair(0, 80)
    overlapping = pair(40, 120)
    summary = trimmed_coverage({0: msa}, [incumbent, overlapping])
    assert summary["0"]["trimmedIntervals"] == [[20, 100]]
    assert summary["0"]["trimmedBases"] == 80
    assert summary["0"]["referenceLength"] == 120


def test_selector_accepts_one_base_gain_against_parent_union():
    msa = target(100)
    parent = pair(0, 60)
    candidate = pair(40, 81)
    result = select_followup_candidates(
        {0: msa},
        [parent],
        [candidate],
        pool_count=2,
        admission=lambda _candidate, _pool: None,
    )
    assert isinstance(result, FollowupSelection)
    assert result.accepted[0].candidate_id == gap_candidate_id(candidate, msa)
    assert result.accepted[0].trimmed_gain == 1
    assert result.followup_coverage["0"]["trimmedBases"] == 1
    assert result.combined_coverage["0"]["trimmedBases"] == 21
    assert result.candidate_union_coverage["0"]["trimmedBases"] == 21


def test_selector_uses_deterministic_candidate_and_pool_ties():
    msa = target(160)
    first = pair(0, 80, "A" * 20, "C" * 20)
    second = pair(80, 160, "G" * 20, "T" * 20)
    calls = []

    def admission(candidate, pool):
        calls.append((gap_candidate_id(candidate, msa), pool))
        return None

    result = select_followup_candidates(
        {0: msa}, [], [second, first], pool_count=2, admission=admission
    )
    assert [item.candidate_id for item in result.accepted] == sorted(
        item.candidate_id for item in result.accepted
    )
    assert {item.pool for item in result.accepted} == {0}
    assert all(pool in {0, 1} for _, pool in calls)


def test_selector_records_invalid_geometry_without_sorting_error():
    msa = target(100)
    invalid = pair(90, 130)
    valid = pair(0, 80)
    result = select_followup_candidates(
        {0: msa}, [], [invalid, valid], pool_count=1, admission=lambda *_: None
    )
    assert result.statuses[gap_candidate_id(invalid, msa)]["reason"] == "geometry"


def test_parent_reference_rejects_duplicate_target_identity(tmp_path):
    reference = tmp_path / "reference.fasta"
    with dnaio.FastaWriter(reference) as writer:
        writer.write(dnaio.SequenceRecord(name="target", sequence="AAAA"))
        writer.write(dnaio.SequenceRecord(name="target", sequence="CCCC"))
    with pytest.raises(ValueError, match="duplicate target"):
        _parent_reference(tmp_path)


def test_bounded_gap_expansion_records_anchor_and_pair_work(tmp_path):
    msa = Path("tests/core/test_mismatch.fasta").resolve()
    config = Config(
        selection_algorithm="legacy",
        mapping=MappingType.FIRST,
        terminal_gap_policy=TerminalGapPolicy.LEGACY,
        n_pools=2,
    )
    config.use_matchdb = False
    parent = tmp_path / "parent"
    panelcreate(
        msa=[msa], output_dir=parent, config=config, pm=ProgressManager(),
        mode=PanelRunModes.EQUAL, offline_plots=False,
    )
    output = tmp_path / "expanded"
    panelcreate(
        msa=[msa], output_dir=output, config=config, pm=ProgressManager(),
        mode=PanelRunModes.EQUAL, offline_plots=False,
        gap_completion_parent=parent,
        gap_expansion_options=GapExpansionOptions(
            mode="bounded", max_anchors_per_msa=20, max_pairs_per_msa=10
        ),
    )
    report = json.loads((output / "gap-completion.json").read_text())
    assert report["gapExpansion"]["options"]["mode"] == "bounded"
    assert report["gapExpansion"]["options"]["maxAnchorsPerMsa"] == 20


def test_gap_anchor_cap_interleaves_windows_directions_and_spatial_positions():
    msa = target(200)
    config = SimpleNamespace(amplicon_size_max=20)
    tasks = _anchor_tasks(msa, [(0, 10, 20), (0, 120, 130)], config)
    assert {(gap, direction) for gap, direction, _anchor in tasks[:4]} == {
        (0, "f"), (0, "r"), (1, "f"), (1, "r")
    }
    assert _spatial_order(list(range(0, 41, 10)))[:2] == [20, 0]
