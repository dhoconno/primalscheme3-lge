"""Bounded gain memoization preserves arithmetic and scientific search decisions."""

from dataclasses import replace
from types import SimpleNamespace

import pytest
from test_allele_search import GraphOracle, instance, run
from test_phase_scheduling import Clock

from primalscheme3.panel import allele_search as allele
from primalscheme3.panel.coverage_history import CoverageHistory
from primalscheme3.panel.coverage_types import ConfigurationLedger


def state_fixture():
    cat, configs, profile = instance(
        {"A": (3, 21), "B": (1, 11), "C": (13, 23)}, pools=2
    )
    views = {c.id: allele._ConfigurationView(c, cat) for c in configs.values()}
    pools = SimpleNamespace(
        diagnostics=lambda ids: {"violating_edge_count": 0, "incident_species_count": 0}
    )
    state = allele._AlleleState(
        cat,
        SimpleNamespace(candidate_by_id=views),
        profile,
        allele.AlleleSearchOptions(starts=1),
        pools,
    )
    return state, cat, configs, views


def test_lru_hits_zero_values_eviction_and_per_state_ownership(monkeypatch):
    monkeypatch.setattr(allele, "_GAIN_CACHE_LIMIT", 2)
    state, _, configs, views = state_fixture()
    a, b, c = [views[configs[k].id] for k in ("A", "B", "C")]
    calls = []
    original = state._fraction

    def fraction(*args):
        calls.append(args)
        return original(*args)

    monkeypatch.setattr(state, "_fraction", fraction)
    expected = {x.id: state.gain(x) for x in (a, b)}
    calls.clear()
    assert state.gain(a) == expected[a.id] and not calls
    state.gain(c)
    assert list(state._gain_cache) == [a.id, c.id]
    assert state.gain(b) == expected[b.id] and calls
    assert len(state._gain_cache) == 2
    state.add(a, 0)
    assert not state._gain_cache
    assert state.gain(a) == 0
    calls.clear()
    assert state.gain(a) == 0 and not calls
    other, *_ = state_fixture()
    assert not other._gain_cache


def test_every_successful_mutation_clears_and_failed_mutations_do_not():
    state, _, configs, views = state_fixture()
    a, b = [views[configs[k].id] for k in ("A", "B")]
    state.gain(a)
    state.add(a, 0)
    assert not state._gain_cache
    state.gain(b)
    retained = dict(state._gain_cache)
    with pytest.raises(ValueError, match="duplicate"):
        state.add(a, 1)
    assert dict(state._gain_cache) == retained
    with pytest.raises(KeyError):
        state.drop(b)
    assert dict(state._gain_cache) == retained
    state.add(b, 0)  # overlapping contributions still invalidate
    assert not state._gain_cache
    state.gain(a)
    state.drop(a)
    assert not state._gain_cache
    state.gain(a)
    state.add(a, 1)  # pool reassignment via drop/add
    assert not state._gain_cache
    assert state.assignments[a.id] == 1


def test_new_registration_keeps_existing_values_and_rejects_collision():
    state, cat, configs, views = state_fixture()
    history = CoverageHistory(run_id="gain-registration")
    ledger = ConfigurationLedger(cat.semantic_digest, (configs["A"],))
    proposals = allele._Proposals(
        cat,
        ledger,
        GraphOracle(configs.values()),
        state.profile,
        allele.StagePolicy(),
        state.options,
        history,
    )
    state.view = proposals.view
    a = proposals.view.candidate_by_id[configs["A"].id]
    old = state.gain(a)
    proposals._register(configs["B"])
    assert state._gain_cache[a.id] == old
    b = proposals.view.candidate_by_id[configs["B"].id]
    assert state.gain(b) == state._uncached_gain(b)
    bad = replace(configs["A"], full_interval=(0, 24))
    assert bad.id == configs["A"].id
    with pytest.raises(ValueError, match="collision"):
        proposals._register(bad)
    assert state.gain(a) == old


@pytest.mark.parametrize("scheduling", ["serial", "reserved"])
@pytest.mark.parametrize("variant_selection", ["subsets", "full-cloud"])
def test_fixed_work_history_objective_and_repair_parity(
    monkeypatch, scheduling, variant_selection
):
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11), "C": (13, 23)})
    edges = [(configs["A"].id, configs[k].id) for k in ("B", "C")]
    options = allele.AlleleSearchOptions(
        starts=2,
        repair_rounds=1,
        seed_modes=(),
        phase_scheduling=scheduling,
        variant_selection=variant_selection,
        time_limit=1000,
    )
    results = []
    for capacity in (0, 8192):
        monkeypatch.setattr(allele, "_GAIN_CACHE_LIMIT", capacity)
        results.append(
            run(
                cat,
                configs,
                profile,
                edges=edges,
                options=options,
                clock=Clock(),
                history=CoverageHistory(run_id="gain-parity"),
            )
        )
    plain, cached = results
    assert plain.assignments == cached.assignments
    assert plain.objective == cached.objective
    assert plain.coverage == cached.coverage
    assert plain.metadata == cached.metadata
    assert plain.ledger == cached.ledger
    for stream in ("evidence", "assessments", "events", "snapshots"):
        assert getattr(plain.history, stream) == getattr(cached.history, stream)
    assert cached.metadata["work"]["repair_trials"] > 0
    assert cached.coverage.mean_coverage == 20 / 24


def test_partial_distinct_class_weighting_multiple_targets_and_goal_transition():
    from allele_fixtures import target

    from primalscheme3.panel.allele_coverage import canonical_observations
    from primalscheme3.panel.coverage_types import VariantCatalog

    targets = (
        target("A" * 12, "A" * 12, "C" * 12, id="long"),
        target("A" * 8, "A" * 6 + "NN", id="short"),
        target("N" * 8, id="unassessable"),
    )
    observations = tuple(a for t in targets for a in canonical_observations(t))
    cat = VariantCatalog(targets, observations, (), ())
    long_a = next(
        a for a in observations if a.target_id == "long" and a.cells[0] == "A"
    )
    assert long_a.multiplicity == 2

    def candidate(cid, tid, positions):
        return SimpleNamespace(
            id=cid,
            target_id=tid,
            contribution=positions,
            forward_oligos=("A",),
            reverse_oligos=("T",),
            full_interval=(0, 8),
        )

    a = candidate("partial", "long", {long_a.id: frozenset(range(6))})
    b = candidate(
        "short",
        "short",
        {x.id: frozenset(range(4)) for x in observations if x.target_id == "short"},
    )
    c = candidate(
        "complete",
        "long",
        {
            x.id: frozenset(x.observed_positions)
            for x in observations
            if x.target_id == "long"
        },
    )
    unknown = candidate("unknown", "unassessable", {})
    view = SimpleNamespace(candidate_by_id={x.id: x for x in (a, b, c, unknown)})
    pools = SimpleNamespace(
        diagnostics=lambda ids: {"violating_edge_count": 0, "incident_species_count": 0}
    )
    state = allele._AlleleState(
        cat,
        view,
        allele.AlleleConstraintProfile(1, 12, 2),
        allele.AlleleSearchOptions(coverage_target=0.5),
        pools,
    )
    assert len(state.by_target) == 2  # unassessable target excluded from denominator
    assert state.gain(a) == pytest.approx(
        0.3125
    )  # equal classes, not duplicate row weight
    assert state.gain(b) == pytest.approx((7 / 12 + 0.5) / 2)
    assert state.gain(unknown) == 0
    for item in (a, b, c):
        assert state.gain(item) == state._uncached_gain(item)
        state.add(item, 0)
        assert not state._gain_cache
        for other in (a, b, c, unknown):
            assert state.gain(other) == state._uncached_gain(other)
    assert state.fraction("long") == 1  # crosses soft goal without gain rounding
    state.drop(c)
    assert state.fraction("long") == 0.25
    assert not state._gain_cache
    state.gain(a)
    state.add(unknown, 1)  # empty contribution still clears
    assert not state._gain_cache


def test_constructor_assignments_start_empty_cache():
    from primalscheme3.panel.coverage_search import Assignment

    state, cat, configs, views = state_fixture()
    a = views[configs["A"].id]
    initialized = allele._AlleleState(
        cat,
        state.view,
        state.profile,
        state.options,
        state.pool_evaluations,
        (Assignment(a.id, 0),),
    )
    assert initialized.assignments == {a.id: 0}
    assert not initialized._gain_cache
    assert initialized.gain(a) == 0


@pytest.mark.parametrize("cancelled", [False, True])
def test_timeout_and_cancellation_parity(monkeypatch, cancelled):
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11)})

    class TickClock:
        calls = 0

        def __call__(self):
            self.calls += 1
            return self.calls * 0.01

    results = []
    counts = []
    for capacity in (0, 8192):
        monkeypatch.setattr(allele, "_GAIN_CACHE_LIMIT", capacity)
        clock = TickClock()
        results.append(
            run(
                cat,
                configs,
                profile,
                clock=clock,
                cancelled=lambda: cancelled,
                options=allele.AlleleSearchOptions(time_limit=0.1),
                history=CoverageHistory(run_id="deadline-parity"),
            )
        )
        counts.append(clock.calls)
    assert counts[0] == counts[1]
    assert results[0].assignments == results[1].assignments
    assert results[0].metadata == results[1].metadata
    assert results[0].history.events == results[1].history.events
    assert results[0].metadata["stop_reason"] == (
        "cancelled" if cancelled else "time-limit"
    )


def test_generated_partial_subset_history_parity(monkeypatch):
    from allele_fixtures import target

    from primalscheme3.panel.allele_coverage import canonical_observations
    from primalscheme3.panel.coverage_types import (
        CandidateFamily,
        OligoSite,
        VariantCatalog,
    )
    from primalscheme3.panel.coverage_variants import make_configuration

    t = target("AACCCCTTAA", "AACCCCGGAA")
    f = OligoSite(t.id, "AA", "+", 2, (0, 2), accepting_profile_ids=("normal",))
    rs = tuple(
        OligoSite(t.id, s, "-", 6, (6, 8), accepting_profile_ids=("normal",))
        for s in ("AA", "CC")
    )
    family = CandidateFamily(t.id, (2, 6), (f.id,), tuple(r.id for r in rs))
    cat = VariantCatalog((t,), canonical_observations(t), (f,) + rs, (family,))
    full = make_configuration(cat, family.id, (f.id,), tuple(r.id for r in rs))
    results = []
    for capacity in (0, 8192):
        monkeypatch.setattr(allele, "_GAIN_CACHE_LIMIT", capacity)
        results.append(
            run(
                cat,
                {"full": full},
                allele.AlleleConstraintProfile(1, 10, 1, mismatch_kmersize=1),
                oracle=GraphOracle((full,), forbidden_sites=(rs[1].id,)),
                clock=Clock(),
                history=CoverageHistory(run_id="subset-parity"),
            )
        )
    a, b = results
    assert a.assignments == b.assignments
    assert a.objective == b.objective
    assert a.coverage == b.coverage
    assert a.metadata == b.metadata
    assert a.ledger == b.ledger
    for stream in ("evidence", "assessments", "events", "snapshots"):
        assert getattr(a.history, stream) == getattr(b.history, stream)
    assert b.coverage.mean_coverage == 0.2
    selected = b.ledger.configuration_by_id[b.assignments[0].configuration_id]
    assert selected.reverse_site_ids == (rs[0].id,)
