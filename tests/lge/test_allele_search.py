import pytest
from allele_fixtures import target
from allele_exhaustive_oracle import exhaustive
from primalscheme3.panel.coverage_types import (
    OligoSite,
    CandidateFamily,
    VariantCatalog,
    ConfigurationLedger,
    AlleleAssignment,
)
from primalscheme3.panel.allele_coverage import canonical_observations
from primalscheme3.panel.coverage_variants import make_configuration
from primalscheme3.panel.coverage_history import CoverageHistory
from primalscheme3.panel.allele_validation import AlleleConstraintProfile


def api():
    from primalscheme3.panel import allele_search

    return allele_search


def instance(intervals, rows=("A" * 24,), pools=1):
    t = target(*rows)
    sites = {}
    families = {}
    labels = {}
    for label, (a, b) in intervals.items():
        f = OligoSite(t.id, "A", "+", a, (a - 1, a), accepting_profile_ids=("normal",))
        r = OligoSite(t.id, "T", "-", b, (b, b + 1), accepting_profile_ids=("normal",))
        fam = CandidateFamily(t.id, (a, b), (f.id,), (r.id,))
        sites.update({f.id: f, r.id: r})
        families[fam.id] = fam
        labels[label] = (fam, f, r)
    cat = VariantCatalog(
        (t,), canonical_observations(t), tuple(sites.values()), tuple(families.values())
    )
    configurations = {
        label: make_configuration(cat, fam.id, (f.id,), (r.id,))
        for label, (fam, f, r) in labels.items()
    }
    return (
        cat,
        configurations,
        AlleleConstraintProfile(1, len(rows[0]), pools, mismatch_kmersize=1),
    )


class GraphOracle:
    def __init__(
        self, configurations, edges=(), invalid=(), forbidden_sites=(), pool_cap=None
    ):
        self.configurations = {c.id: c for c in configurations}
        self.edges = {frozenset(e) for e in edges}
        self.invalid = set(invalid)
        self.forbidden = set(forbidden_sites)
        self.pool_cap = pool_cap
        self.pool_calls = []

    def register(self, *configs):
        self.configurations.update({c.id: c for c in configs})

    def candidate_valid(self, cid):
        c = self.configurations[cid]
        return cid not in self.invalid and not self.forbidden.intersection(
            c.forward_site_ids + c.reverse_site_ids
        )

    def conflict(self, a, b):
        return frozenset((a, b)) in self.edges

    def pool_diagnostics(self, ids):
        ids = tuple(sorted(ids))
        self.pool_calls.append(ids)
        return {
            "reasons": (
                ["aggregate"]
                if self.pool_cap is not None and len(ids) > self.pool_cap
                else []
            ),
            "violating_edge_count": 0,
            "incident_species_count": 0,
        }

    def pool_valid(self, ids):
        return not self.pool_diagnostics(ids)["reasons"]


def run(
    cat,
    configs,
    profile,
    *,
    edges=(),
    options=None,
    baseline=(),
    history=None,
    clock=None,
    cancelled=None,
    oracle=None,
):
    a = api()
    ledger = ConfigurationLedger(cat.semantic_digest, tuple(configs.values()))
    return a.search_allele_configurations(
        cat,
        profile,
        options or a.AlleleSearchOptions(starts=1),
        ledger,
        oracle or GraphOracle(configs.values(), edges),
        baseline=baseline,
        history=history,
        **({} if clock is None else {"clock": clock}),
        cancelled=cancelled,
    )


def test_partial_support_and_unobserved_class_match_literal_oracle():
    cat, configs, profile = instance(
        {"part": (2, 18)}, rows=("A" * 24, "C" * 24, "N" * 24)
    )
    result = run(cat, configs, profile)
    expected, _ = exhaustive(cat, {c.id: c for c in configs.values()}, profile)
    assert result.objective == expected
    assert (
        result.validation["valid"] and result.validation["scope"] == "abstract-oracle"
    )
    assert sorted(c.covered_count for c in result.coverage.classes) == [0, 0, 16]
    assert result.coverage.mean_coverage == 1 / 3


def test_one_start_exchange_drops_blocker_and_adds_two():
    a = api()
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11), "C": (13, 23)})
    edges = [(configs["A"].id, configs[x].id) for x in ("B", "C")]
    greedy = run(
        cat,
        configs,
        profile,
        edges=edges,
        options=a.AlleleSearchOptions(starts=1, repair_rounds=0, seed_modes=()),
    )
    repaired = run(
        cat,
        configs,
        profile,
        edges=edges,
        options=a.AlleleSearchOptions(starts=1, repair_rounds=1, seed_modes=()),
    )
    expected, _ = exhaustive(
        cat, {c.id: c for c in configs.values()}, profile, edges=edges
    )
    assert greedy.coverage.mean_coverage == 18 / 24
    assert repaired.objective == expected
    assert repaired.coverage.mean_coverage == 20 / 24


def test_pool_move_and_aggregate_guard_match_tiny_oracle():
    cat, configs, profile = instance({"A": (1, 11), "B": (13, 23)}, pools=2)
    edge = (configs["A"].id, configs["B"].id)
    oracle = GraphOracle(configs.values(), (edge,), pool_cap=1)
    result = run(cat, configs, profile, edges=(edge,), oracle=oracle)
    expected, _ = exhaustive(
        cat, {c.id: c for c in configs.values()}, profile, edges=(edge,)
    )
    assert result.objective == expected
    assert len({x.pool for x in result.assignments}) == 2
    assert all(len(x) <= 2 for x in oracle.pool_calls)


def test_variant_deletion_rescues_partial_coverage():
    t = target("AACCCCTTAA", "AACCCCGGAA")
    f = OligoSite(t.id, "AA", "+", 2, (0, 2), accepting_profile_ids=("normal",))
    rs = tuple(
        OligoSite(t.id, s, "-", 6, (6, 8), accepting_profile_ids=("normal",))
        for s in ("AA", "CC")
    )
    fam = CandidateFamily(t.id, (2, 6), (f.id,), tuple(r.id for r in rs))
    cat = VariantCatalog((t,), canonical_observations(t), (f,) + rs, (fam,))
    full = make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs))
    profile = AlleleConstraintProfile(1, 10, 1, mismatch_kmersize=1)
    oracle = GraphOracle((full,), forbidden_sites=(rs[1].id,))
    result = run(
        cat,
        {"full": full},
        profile,
        oracle=oracle,
        history=CoverageHistory(run_id="rescue"),
    )
    assert result.coverage.mean_coverage == 0.2
    assert len(result.assignments) == 1
    selected = result.ledger.configuration_by_id[result.assignments[0].configuration_id]
    assert selected.reverse_site_ids == (rs[0].id,)
    assert result.metadata["proposal_work"]["expanded_states"] > 0


def test_seed_retention_cancellation_and_fixed_work():
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11), "C": (13, 23)})
    baseline = (
        AlleleAssignment(configs["B"].id, 0, "strict"),
        AlleleAssignment(configs["C"].id, 0, "strict"),
    )
    result = run(cat, configs, profile, baseline=baseline, cancelled=lambda: True)
    assert result.validation["valid"] and result.coverage.mean_coverage == 20 / 24
    assert result.metadata["stop_reason"] == "cancelled"
    assert result.metadata["proposal_work"]["unsearched_families"] == len(cat.families)
    one = run(cat, configs, profile)
    two = run(cat, configs, profile)
    assert one.assignments == two.assignments and one.objective == two.objective
    assert one.metadata["work"] == two.metadata["work"]


def test_invalid_work_controls_reject_before_search():
    a = api()
    for value in (True, 1.5, float("nan"), float("inf"), 0, -1):
        with pytest.raises(ValueError):
            a.AlleleWorkLimits(frontier_candidates=value)
    with pytest.raises(ValueError):
        a.AlleleSearchOptions(exchange_width=3)


def test_timeout_retains_seed_and_reports_unsearched_families():
    a = api()
    cat, configs, profile = instance({"A": (3, 21), "B": (1, 11), "C": (13, 23)})
    baseline = (AlleleAssignment(configs["A"].id, 0, "strict"),)

    class Clock:
        value = 0

        def __call__(self):
            self.value += 1
            return self.value

    result = run(
        cat,
        configs,
        profile,
        baseline=baseline,
        clock=Clock(),
        options=a.AlleleSearchOptions(time_limit=0.1, starts=1),
    )
    assert result.metadata["stop_reason"] == "time-limit"
    assert result.coverage.mean_coverage == 18 / 24
    assert result.metadata["proposal_work"]["unsearched_families"] == 3
    assert result.validation["valid"]


def test_reference_count_removal_preserves_overlapping_contributors():
    a = api()
    cat, configs, profile = instance({"A": (1, 15), "B": (8, 23)})
    history = CoverageHistory(run_id="counts")
    oracle = GraphOracle(configs.values())
    view = a._CatalogView(cat)
    view.candidate_by_id.update(
        {c.id: a._ConfigurationView(c, cat) for c in configs.values()}
    )
    pools = a._PoolEvaluations(oracle, cat, profile, a.StagePolicy(), history)
    state = a._AlleleState(cat, view, profile, a.AlleleSearchOptions(), pools)
    state.add(view.candidate_by_id[configs["A"].id], 0)
    state.add(view.candidate_by_id[configs["B"].id], 0)
    assert state.fraction(cat.targets[0].id) == 22 / 24
    state.drop(view.candidate_by_id[configs["A"].id])
    assert state.fraction(cat.targets[0].id) == 15 / 24
    assert state.oligos[0] == {"A", "T"}


def test_aggregate_failure_is_contextual_and_second_pool_remains_usable():
    cat, configs, profile = instance({"A": (1, 11), "B": (13, 23)}, pools=2)
    h = CoverageHistory(run_id="pools")
    result = run(
        cat,
        configs,
        profile,
        oracle=GraphOracle(configs.values(), pool_cap=1),
        history=h,
    )
    assert len(result.assignments) == 2
    failed = [x for x in h.assessments if x.outcome == "fail"]
    passed = [x for x in h.assessments if x.outcome == "pass"]
    assert failed and passed
    assert any(
        set(f.entity_ids).intersection(p.entity_ids) and f.pool != p.pool
        for f in failed
        for p in passed
    )
    assert all(
        result.metadata["candidate_status"][x.configuration_id]
        for x in result.assignments
    )
    assert len({x.context_digest for x in h.assessments}) > 1


def test_bounded_witness_deletion_precedes_family_loss():
    a = api()
    t = target("AACCCCTTAA", "AACCCCGGAA")
    f = OligoSite(t.id, "AA", "+", 2, (0, 2), accepting_profile_ids=("normal",))
    rs = tuple(
        OligoSite(t.id, s, "-", 6, (6, 8), accepting_profile_ids=("normal",))
        for s in ("AA", "CC")
    )
    fam = CandidateFamily(t.id, (2, 6), (f.id,), tuple(r.id for r in rs))
    cat = VariantCatalog((t,), canonical_observations(t), (f,) + rs, (fam,))
    full = make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs))

    class WitnessOracle(GraphOracle):
        def candidate_diagnostics(self, cid):
            if self.candidate_valid(cid):
                return {"reasons": [], "dimer_edges": []}
            return {
                "reasons": ["dimer-threshold"],
                "dimer_edges": [
                    {"id": "right-blocker", "site_ids": (rs[1].id,), "score": -30}
                ],
            }

    oracle = WitnessOracle((full,), forbidden_sites=(rs[1].id,))
    h = CoverageHistory(run_id="bounded-witness")
    result = run(
        cat,
        {"full": full},
        AlleleConstraintProfile(1, 10, 1, mismatch_kmersize=1),
        oracle=oracle,
        history=h,
        options=a.AlleleSearchOptions(
            starts=1, repair_rounds=1, subset_expansion_limit=1
        ),
    )
    assert result.coverage.mean_coverage == 0.2
    assert any("right-blocker" in e.changes.get("witness_ids", ()) for e in h.events)


def test_production_entry_freshly_validates_and_completes_history():
    from test_allele_validation import fixture, profile, dna

    a = api()
    cat, configs, ledger = fixture(seq=dna(seed=0))
    result = a.search_allele_assignments(
        cat, profile(), a.AlleleSearchOptions(starts=1, repair_rounds=0), ledger
    )
    assert result.validation["valid"], result.validation["violations"]
    assert result.assignments
    assert result.coverage.to_dict() == result.validation["coverage"]
    assert (
        result.history.last_complete_stage.id == result.metadata["history_snapshot_id"]
    )
    assert (
        result.validation["validation_scope"]
        == "catalog-authoritative-rows-and-in-memory-selection"
    )


def test_no_assessable_targets_skips_search():
    t = target("NNNN")
    cat = VariantCatalog((t,), canonical_observations(t), (), ())
    result = run(cat, {}, AlleleConstraintProfile(1, 4, 1, mismatch_kmersize=1))
    assert result.coverage.status == "no-assessable-targets"
    assert result.metadata["stop_reason"] == "no-assessable-targets"
    assert result.metadata["completed_starts"] == 0


def test_strict_hot_path_uses_incremental_guard_without_full_edge_materialization():
    class CompactOracle(GraphOracle):
        def pool_valid(self, ids):
            return True

        def pool_diagnostics(self, ids):
            raise AssertionError("strict valid hot path enumerated all pool edges")

    cat, configs, profile = instance({"A": (1, 11), "B": (13, 23)})
    result = run(cat, configs, profile, oracle=CompactOracle(configs.values()))
    assert result.coverage.mean_coverage == 20 / 24


def test_full_cloud_ablation_never_selects_normal_subset_or_proposes_variants():
    a = api()
    t = target("AACCCCTTAA", "AACCCCGGAA")
    f = OligoSite(t.id, "AA", "+", 2, (0, 2), accepting_profile_ids=("normal",))
    r = OligoSite(t.id, "AA", "-", 6, (6, 8), accepting_profile_ids=("normal",))
    bad = OligoSite(t.id, "CC", "-", 6, (6, 8), accepting_profile_ids=("high-gc",))
    fam = CandidateFamily(t.id, (2, 6), (f.id,), (r.id, bad.id))
    cat = VariantCatalog((t,), canonical_observations(t), (f, r, bad), (fam,))
    full = make_configuration(cat, fam.id, (f.id,), (r.id, bad.id))
    sub = make_configuration(cat, fam.id, (f.id,), (r.id,))
    profile = AlleleConstraintProfile(1, 10, 1, mismatch_kmersize=1)
    h = CoverageHistory(run_id="full-only")
    result = run(
        cat,
        {"full": full, "sub": sub},
        profile,
        oracle=GraphOracle((full, sub), forbidden_sites=(bad.id,)),
        history=h,
        baseline=(AlleleAssignment(sub.id, 0, "strict"),),
        options=a.AlleleSearchOptions(starts=1, variant_selection="full-cloud"),
    )
    assert not result.assignments
    assert result.metadata["effective_seed_modes"] == ["full"]
    assert not any(e.kind == "configuration-proposed" for e in h.events)
    assert result.metadata["proposal_work"].get("expanded_states", 0) == 0
    assert not result.metadata["baseline"]["supplied_valid"]


def test_parent_neighborhood_calls_obey_limit_with_many_same_family_configs():
    from test_coverage_variants import fixture

    a = api()
    cat, fam, f, rs = fixture()
    configs = tuple(make_configuration(cat, fam.id, (f.id,), (r.id,)) for r in rs)
    profile = AlleleConstraintProfile(1, 10, 1, mismatch_kmersize=1)
    oracle = GraphOracle(configs)
    options = a.AlleleSearchOptions(
        work_limits=a.AlleleWorkLimits(repair_neighborhoods_per_round=2)
    )
    h = CoverageHistory(run_id="neighborhood-bound")
    proposals = a._Proposals(
        cat,
        ConfigurationLedger(cat.semantic_digest, configs),
        oracle,
        profile,
        a.StagePolicy(),
        options,
        h,
    )
    pools = a._PoolEvaluations(oracle, cat, profile, a.StagePolicy(), h)
    state = a._AlleleState(cat, proposals.view, profile, options, pools)
    for c in configs:
        state.add(proposals.view.candidate_by_id[c.id], 0)
    search = a._Search(proposals.view, profile, options, oracle, lambda: 0)
    proposals.neighborhoods(search, state, 0, 0)
    assert proposals.work["neighborhood_calls"] <= 2


def test_rejected_candidate_retains_performed_measurements_and_unperformed_checks():
    from test_allele_validation import fixture, profile

    a = api()
    cat, configs, ledger = fixture()
    result = a.search_allele_assignments(
        cat,
        profile(),
        a.AlleleSearchOptions(starts=1, repair_rounds=0, seed_modes=()),
        ledger,
    )
    assert not result.assignments
    assert (
        result.history.last_complete_stage.dispositions[configs[0].id]
        == "rejected-stage-policy"
    )
    assessments = result.history.query(entity_id=configs[0].id)["assessments"]
    assert any(x.outcome == "fail" and "dimer" in x.check_name for x in assessments)
    assert any(x.evidence_ids for x in assessments)
    measured = result.history.query(entity_id=configs[0].id)["evidence"]
    assert any(e.values.get("dimer_edges") for e in measured)
