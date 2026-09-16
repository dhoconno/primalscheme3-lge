from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from primalschemers import FKmer, RKmer, do_pool_interact

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.multiplex import PrimerPairCheck
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.legacy_salvage import (
    LegacySalvageOptions,
    legacy_candidate_id,
    native_dimer_score,
    run_legacy_salvage,
    validate_legacy_salvage,
)
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate


def pair(start: int, end: int, forward: str, reverse: str) -> PrimerPair:
    return PrimerPair(
        FKmer([forward.encode()], start + len(forward)),
        RKmer([reverse.encode()], end - len(reverse)),
        0,
    )


class FakePanel:
    def __init__(self, candidates, strict=(), strict_by_pool=None):
        self.n_pools = 2
        self._pools = [list(strict), []] if strict_by_pool is None else [list(items) for items in strict_by_pool]
        for pool, items in enumerate(self._pools):
            for item in items:
                item.pool = pool
        self._msa_dict = {}
        self._last_pp_added = list(strict)
        self._candidates = candidates

    def all_primerpairs(self):
        return [item for pool in self._pools for item in pool]

    def check_primerpair_can_be_added(self, primerpair, pool, otherseqs_bytes=None, *, dimer_score=None):
        if any(
            min(current.fprimer.starts()) < max(primerpair.rprimer.ends())
            and min(primerpair.fprimer.starts()) < max(current.rprimer.ends())
            and current.msa_index == primerpair.msa_index
            for current in self._pools[pool]
        ):
            return PrimerPairCheck.OVERLAP
        cutoff = -26 if dimer_score is None else dimer_score
        seqs = otherseqs_bytes or [seq for item in self._pools[pool] for seq in item.all_seq_bytes()]
        if do_pool_interact(primerpair.all_seq_bytes(), seqs, cutoff):
            return PrimerPairCheck.INTERACTING
        return PrimerPairCheck.OK

    def _add_primerpair(self, primerpair, pool, msa_index):
        primerpair.pool = pool
        self._pools[pool].append(primerpair)
        self._last_pp_added.append(primerpair)


def msa(candidates, *, length=120, index=0):
    return SimpleNamespace(
        msa_index=index,
        _chrom_name=f"target-{index}",
        _mapping_array=list(range(length)),
        array=[["A"] * length],
        primerpairs=list(candidates),
    )


def test_options_validate_threshold_floor_and_exposure_defaults():
    options = LegacySalvageOptions(mode="bounded")
    assert options.thresholds == (-28.0, -30.0, -32.0)
    assert options.floor == -32.0
    assert options.max_edges_per_pool == 8
    assert options.max_incident_species_per_pool == 4
    assert options.min_reference_gain == 1

    with pytest.raises(ValueError):
        LegacySalvageOptions(mode="bounded", thresholds=(-28, -33))
    with pytest.raises(ValueError):
        LegacySalvageOptions(mode="bounded", floor=-31, thresholds=(-32,))
    with pytest.raises(ValueError):
        LegacySalvageOptions(mode="bounded", min_reference_gain=0)
    assert LegacySalvageOptions(mode="off", thresholds=(-1,)).mode == "off"


def test_native_score_matches_boolean_kernel_boundary():
    score = native_dimer_score("CGATTCAAATGACGGCAGCA", "AGAACCGAGTGCTGACGTAA")
    assert score.score <= -29
    assert do_pool_interact(
        [score.sequences[0].encode()], [score.sequences[1].encode()], score.score
    )
    assert not do_pool_interact(
        [score.sequences[0].encode()], [score.sequences[1].encode()], score.score - 0.01
    )


def test_salvage_rescues_positive_trimmed_gap_at_next_cutoff():
    incumbent = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    incumbent_other_pool = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    rescued = pair(100, 200, "AGAACCGAGTGCTGACGTAA", "T" * 20)
    panel = FakePanel(
        [rescued],
        strict=(incumbent, incumbent_other_pool),
        strict_by_pool=((incumbent,), (incumbent_other_pool,)),
    )
    target = msa([rescued], length=220)

    result = run_legacy_salvage(
        panel,
        {0: target},
        LegacySalvageOptions(
            mode="bounded", thresholds=(-28, -30), max_edges_per_pool=8, max_incident_species_per_pool=4
        ),
    )

    assert result.stages[0].accepted == ()
    assert result.stages[1].accepted == (legacy_candidate_id(rescued, target),)
    assert panel._pools[0][-1] is rescued
    assert result.stages[0].coverage["0"]["trimmedBases"] == 60


def test_exposure_edges_are_deduplicated_and_cumulative_across_passes():
    incumbent = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    incumbent_other_pool = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    first = pair(100, 210, "AGAACCGAGTGCTGACGTAA", "T" * 20)
    second = pair(220, 300, "TATATCGCAATGCTGCGTGG", "C" * 20)
    panel = FakePanel(
        [first, second],
        strict=(incumbent, incumbent_other_pool),
        strict_by_pool=((incumbent,), (incumbent_other_pool,)),
    )
    panel.n_pools = 1
    panel._pools = [panel._pools[0]]
    target = msa([first, second], length=320)

    result = run_legacy_salvage(
        panel,
        {0: target},
        LegacySalvageOptions(
            mode="bounded",
            thresholds=(-28, -30),
            max_edges_per_pool=1,
            max_incident_species_per_pool=2,
        ),
    )

    assert result.stages[0].accepted == ()
    assert result.stages[1].accepted == (legacy_candidate_id(first, target),)
    rejection = result.stages[1].rejections[legacy_candidate_id(second, target)]
    assert rejection["reason"] == "pool-rejected"
    assert set(rejection["pools"].values()) == {"edge-budget"}


def test_disabled_mode_does_not_mutate_panel():
    candidate = pair(0, 40, "A" * 20, "C" * 20)
    panel = FakePanel([candidate])
    target = msa([candidate])
    result = run_legacy_salvage(panel, {0: target}, LegacySalvageOptions())
    assert result.stages == ()
    assert panel.all_primerpairs() == []


def test_fresh_validation_replays_strict_assignment_and_salvage():
    incumbent = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    incumbent_other_pool = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    rescued = pair(100, 200, "AGAACCGAGTGCTGACGTAA", "T" * 20)
    panel = FakePanel(
        [rescued],
        strict=(incumbent, incumbent_other_pool),
        strict_by_pool=((incumbent,), (incumbent_other_pool,)),
    )
    target = msa([rescued], length=220)
    run = run_legacy_salvage(
        panel,
        {0: target},
        LegacySalvageOptions(mode="bounded", thresholds=(-28, -30)),
    )
    validation = validate_legacy_salvage(panel, {0: target}, run)
    assert validation["valid"]
    assert validation["strictPoolAssignments"] == {
        legacy_candidate_id(incumbent, target): [0, 1],
    }

    panel._pools[0].remove(incumbent)
    panel._pools[1].append(incumbent)
    assert not validate_legacy_salvage(panel, {0: target}, run)["valid"]


def test_salvage_thresholds_use_configured_strict_cutoff():
    options = LegacySalvageOptions(mode="bounded", thresholds=(-30, -32), floor=-32)
    options.validate_for_strict_cutoff(-28)
    with pytest.raises(ValueError):
        options.validate_for_strict_cutoff(-31)


def test_salvage_prefers_pool_without_relaxed_exposure():
    incumbent = pair(0, 100, "CGATTCAAATGACGGCAGCA", "C" * 20)
    strict_pool_one = pair(300, 400, "G" * 20, "G" * 20)
    rescued = pair(100, 200, "AGAACCGAGTGCTGACGTAA", "T" * 20)
    panel = FakePanel(
        [rescued],
        strict_by_pool=((incumbent,), (strict_pool_one,)),
    )
    target = msa([rescued], length=420)
    result = run_legacy_salvage(
        panel,
        {0: target},
        LegacySalvageOptions(mode="bounded", thresholds=(-30,)),
    )
    assert result.stages[0].accepted == (legacy_candidate_id(rescued, target),)
    assert panel._pools[1][-1] is rescued
    evaluation = result.stages[0].evaluations[0]
    assert evaluation.cumulative_edges == 0


def test_bounded_panelcreate_writes_fresh_validation_and_raw_input_provenance(tmp_path):
    source = Path("tests/core/test_mismatch.fasta").absolute()
    output = tmp_path / "panel"
    config = Config()
    config.mapping = MappingType.FIRST
    config.use_matchdb = False
    panelcreate(
        [source],
        output,
        config,
        ProgressManager(),
        force=True,
        mode=PanelRunModes.EQUAL,
        offline_plots=False,
        legacy_salvage_options=LegacySalvageOptions(
            mode="bounded", thresholds=(-28,), max_candidate_evaluations=20
        ),
        executed_argv=["panel-create", "--legacy-salvage", "bounded"],
    )
    provenance = json.loads((output / "panel-provenance.json").read_text())
    assert provenance["status"] == "success"
    assert provenance["inputs"][0]["storedPath"] == "work/0000-test_mismatch.fasta"
    assert provenance["inputs"][0]["sourceAtStartMatchesStored"]
    salvage = json.loads((output / "legacy-salvage.json").read_text())
    assert salvage["schemaVersion"].endswith("legacy-dimer-salvage/v1")
    assert json.loads((output / "legacy-salvage-validation.json").read_text())["valid"]
    assert (output / "strict-primer.bed").is_file()
