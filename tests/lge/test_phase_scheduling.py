"""Cooperative phase reservations serve construction and actual exchanges."""

from dataclasses import replace

import pytest
from test_allele_search import instance, run

from primalscheme3.panel import allele_search as allele
from primalscheme3.panel.coverage_search import _Search


class Clock:
    def __init__(self):
        self.value = 0.0

    def __call__(self):
        return self.value

    def advance(self, amount=1):
        self.value += amount


def test_reserved_seed_and_repair_preparation_cutoffs_still_reach_exchange(monkeypatch):
    clock = Clock()
    original = _Search.fill

    def slow_seed(search, state, start, phase):
        if phase.startswith("seed:"):
            while True:
                clock.advance()
                search.tick()
        return original(search, state, start, phase)

    def slow_prepare(self, search, *args):
        while True:
            clock.advance()
            search.tick()

    def slow_cleanup(search, phase):
        while True:
            clock.advance()
            search.tick()

    monkeypatch.setattr(_Search, "fill", slow_seed)
    monkeypatch.setattr(allele._Proposals, "neighborhoods", slow_prepare)
    monkeypatch.setattr(_Search, "cleanup", slow_cleanup)
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11), "C": (13, 23)})
    edges = [(configs["A"].id, configs[x].id) for x in ("B", "C")]
    result = run(
        cat,
        configs,
        profile,
        edges=edges,
        clock=clock,
        options=allele.AlleleSearchOptions(
            starts=1, repair_rounds=1, time_limit=100, phase_scheduling="reserved"
        ),
    )
    assert result.validation["valid"]
    assert result.coverage.mean_coverage == 20 / 24
    assert result.metadata["work"]["repair_trials"] > 0
    progress = {p["phase"]: p for p in result.metadata["phase_progress"]}
    for name in (
        "seed:full",
        "seed:normal",
        "repair:0:0/preparation",
        "repair:0:0/cleanup",
    ):
        assert progress[name]["outcome"] == "phase-time-limit"
    assert progress["construction:0"]["work_delta"]["candidate_attempts"] > 0
    assert progress["repair:0:0/exchange"]["work_delta"]["repair_trials"] > 0
    assert not result.metadata["fixed_work_completed"]
    assert clock.value <= 100


def test_nonbinding_reserved_matches_serial_and_repeats_fixed_work():
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11), "C": (13, 23)})
    options = allele.AlleleSearchOptions(
        starts=2, repair_rounds=1, time_limit=100, phase_scheduling="serial"
    )
    serial = run(cat, configs, profile, clock=Clock(), options=options)
    reserved = [
        run(
            cat,
            configs,
            profile,
            clock=Clock(),
            options=replace(options, phase_scheduling="reserved"),
        )
        for _ in range(2)
    ]
    for result in reserved:
        assert result.assignments == serial.assignments
        assert result.objective == serial.objective
        assert result.metadata["work"] == serial.metadata["work"]
        assert result.metadata["fixed_work_completed"]
        assert (
            result.metadata["scheduling_policy"]["id"]
            == "initial-phase-reservations/v1"
        )
    assert reserved[0].ledger == reserved[1].ledger


def test_global_overshoot_stops_later_phases_and_preserves_incumbent(monkeypatch):
    from primalscheme3.panel.coverage_types import AlleleAssignment

    clock = Clock()

    def oversized(search, state, start, phase):
        clock.advance(150)
        search.tick()

    monkeypatch.setattr(_Search, "fill", oversized)
    cat, configs, profile = instance({"A": (3, 21)})
    baseline = (AlleleAssignment(configs["A"].id, 0, "strict"),)
    result = run(
        cat,
        configs,
        profile,
        baseline=baseline,
        clock=clock,
        options=allele.AlleleSearchOptions(
            starts=1, time_limit=100, phase_scheduling="reserved"
        ),
    )
    assert result.assignments == baseline
    assert result.metadata["stop_reason"] == "time-limit"
    assert result.metadata["active_phase_at_stop"] == "seed:full"
    assert result.metadata["deadline_overshoot_seconds"] == 50
    assert "construction:0" in result.metadata["phases_not_entered"]
    assert not result.metadata["fixed_work_completed"]


@pytest.mark.parametrize("policy", ["serial", "reserved"])
def test_work_cap_is_not_family_exhaustion_but_can_complete_fixed_recipe(policy):
    cat, configs, profile = instance({str(i): (i + 1, i + 4) for i in range(12)})
    result = run(
        cat,
        configs,
        profile,
        clock=Clock(),
        options=allele.AlleleSearchOptions(
            starts=1,
            repair_rounds=0,
            phase_scheduling=policy,
            work_limits=allele.AlleleWorkLimits(
                families_per_refresh=1, construction_candidate_attempts=1
            ),
        ),
    )
    assert result.metadata["fixed_work_completed"]
    assert result.metadata["exhausted_seed_modes"] == []
    assert all(p["outcome"] == "work-cap" for p in result.metadata["phase_progress"])
    assert not result.metadata["family_streams_exhausted"]["expanded"]


def test_native_option_roundtrip_and_capability_policy():
    from primalscheme3.core.config import Config
    from primalscheme3.panel.allele_options import AlleleOptions
    from primalscheme3.panel.coverage_provenance import capabilities_document

    config = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        phase_scheduling="reserved",
    )
    options = AlleleOptions.from_config(config)
    assert options.search_options().phase_scheduling == "reserved"
    assert Config(**config.to_dict()).to_dict() == config.to_dict()
    capability = capabilities_document()["alleleCoverage"]["phaseScheduling"]
    assert capability["default"] == "serial"
    assert capability["policies"]["reserved"]["repair_weights"] == {
        "preparation": 0.2,
        "cleanup": 0.2,
        "exchange": 0.6,
    }


@pytest.mark.parametrize("policy", ["serial", "reserved"])
def test_cancellation_stops_without_spending_reserved_work(policy):
    from primalscheme3.panel.coverage_types import AlleleAssignment

    cat, configs, profile = instance({"A": (3, 21)})
    baseline = (AlleleAssignment(configs["A"].id, 0, "strict"),)
    result = run(
        cat,
        configs,
        profile,
        baseline=baseline,
        clock=Clock(),
        cancelled=lambda: True,
        options=allele.AlleleSearchOptions(phase_scheduling=policy),
    )
    assert result.assignments == baseline
    assert result.metadata["stop_reason"] == "cancelled"
    assert result.metadata["active_phase_at_stop"] == "seed:full"
    assert result.metadata["work"]["candidate_attempts"] == 0
    assert not result.metadata["fixed_work_completed"]


def test_reserved_full_cloud_retains_ablation_contract():
    from test_allele_search import (
        test_full_cloud_ablation_never_selects_normal_subset_or_proposes_variants,
    )

    # Exercise the original adversarial full-cloud fixture with only the policy changed.
    original = allele.AlleleSearchOptions
    from unittest.mock import patch

    def reserved_options(**kwargs):
        return original(**kwargs, phase_scheduling="reserved")

    with patch.object(allele, "AlleleSearchOptions", reserved_options):
        test_full_cloud_ablation_never_selects_normal_subset_or_proposes_variants()


def test_reserved_salvage_inherits_policy_with_fresh_tier_budget():
    from test_allele_validation import dna, fixture, profile

    from primalscheme3.panel.coverage_salvage import SalvageOptions, run_salvage

    cat, configs, ledger = fixture(seq=dna(seed=91))
    strict = allele.search_allele_assignments(
        cat,
        profile(),
        allele.AlleleSearchOptions(
            starts=1, repair_rounds=1, time_limit=10, phase_scheduling="reserved"
        ),
        ledger,
    )
    result = run_salvage(
        strict,
        profile(),
        SalvageOptions(
            mode="bounded",
            thresholds=(-1000,),
            max_edges_per_pool=8,
            max_oligos_per_pool=4,
            time_limit=8,
        ),
    )
    assert result.tiers
    tier = result.tiers[0].result
    assert tier.metadata["options"]["phase_scheduling"] == "reserved"
    assert tier.metadata["options"]["time_limit"] == 8
    assert tier.metadata["scheduling_policy"] == strict.metadata["scheduling_policy"]
    assert tier.metadata["phase_progress"][0]["local_budget_seconds"] <= 0.8
    assert tier.validation["valid"]


def test_harness_rejects_local_interruption_even_with_global_completed():
    from test_cached_panel_benchmark import api

    metadata = {
        "stop_reason": "completed",
        "completed_starts": 1,
        "completed_repair_rounds": 2,
        "effective_seed_modes": ["full", "normal"],
        "completed_seed_modes": ["full", "normal"],
        "options": {"starts": 1, "repair_rounds": 2},
    }
    assert api().completed_fixed_work(metadata)
    assert not api().completed_fixed_work(metadata | {"fixed_work_completed": False})
    assert not api().completed_fixed_work(
        metadata | {"fixed_work_completed": True, "stop_reason": "cancelled"}
    )


def test_phase_lifecycle_is_queryable_without_synthetic_primer_entities():
    cat, configs, profile = instance({"A": (3, 21)})
    result = run(
        cat,
        configs,
        profile,
        clock=Clock(),
        options=allele.AlleleSearchOptions(
            starts=1, repair_rounds=1, phase_scheduling="reserved"
        ),
    )
    events = [
        e
        for e in result.history.events
        if e.kind in ("phase-started", "phase-finished")
    ]
    assert len(events) == 2 * len(result.metadata["phase_progress"])
    assert all(not e.entity_ids and e.stage_id == "strict" for e in events)
    assert {e.changes["phase"] for e in events} == {
        p["phase"] for p in result.metadata["phase_progress"]
    }


@pytest.mark.parametrize("value", ["unknown", True, None])
def test_invalid_policy_is_rejected(value):
    with pytest.raises(ValueError, match="phase_scheduling"):
        allele.AlleleSearchOptions(phase_scheduling=value)


def test_cli_exposes_policy_and_rejects_it_for_legacy():
    from typer.testing import CliRunner

    from primalscheme3.cli import app
    from primalscheme3.core.config import Config

    help_result = CliRunner().invoke(app, ["panel-create", "--help"])
    assert help_result.exit_code == 0
    assert "--phase-scheduling" in help_result.stdout
    with pytest.raises(ValueError):
        Config(phase_scheduling="reserved")
