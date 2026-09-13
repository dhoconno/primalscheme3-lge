"""Independent interval/graph objectives plus real synthetic-science integration."""

import random
from dataclasses import replace
from itertools import product

import pytest

from primalscheme3.core.config import Config
from primalscheme3.panel.coverage_search import (
    SearchOptions,
    optimize_catalog,
    search_assignments,
)
from primalscheme3.panel.coverage_types import Assignment, Candidate, Catalog, Target
from primalscheme3.panel.coverage_validation import ConstraintProfile


def graph_catalog(intervals, lengths=None):
    lengths = lengths or {"t": 300}
    targets = tuple(
        Target(t, i, 0, (), (), "", length, (), ())
        for i, (t, length) in enumerate(lengths.items())
    )
    candidates = tuple(
        Candidate(
            name,
            target,
            (a, b),
            (a + 1, b - 1),
            (name + "F",),
            (name + "R",),
            (),
            (),
            (),
        )
        for name, target, a, b in intervals
    )
    return Catalog(targets, candidates, "abstract-graph", "{}", (), ())


def limits(pools=1, cap=None, per_target=None):
    return ConstraintProfile.from_config(
        Config(
            amplicon_size=150,
            amplicon_size_min=40,
            amplicon_size_max=300,
            n_pools=pools,
            mismatch_product_size=2000,
        ),
        max_amplicons=cap,
        max_amplicons_msa=per_target,
    )


class GraphOracle:
    """Explicit graph instance, independent of the production science oracle."""

    def __init__(self, cat, edges=(), invalid=()):
        self.cat = cat
        self.edges = {frozenset(e) for e in edges}
        self.invalid = set(invalid)
        self.counters = {"candidate_calls": 0, "pair_calls": 0}

    def candidate_valid(self, ident):
        self.counters["candidate_calls"] += 1
        return ident not in self.invalid

    def conflict(self, a, b):
        self.counters["pair_calls"] += 1
        return frozenset((a, b)) in self.edges


def run(
    cat, edges=(), *, options=None, profile=None, baseline=(), clock=None, invalid=()
):
    kwargs = {} if clock is None else {"clock": clock}
    return search_assignments(
        cat,
        profile or limits(),
        options or SearchOptions(starts=1),
        GraphOracle(cat, edges, invalid),
        baseline,
        **kwargs,
    )


def covered(cat, assignments, metric="full-span"):
    result = {t.id: set() for t in cat.targets}
    for a in assignments:
        c = cat.candidate_by_id[a.candidate_id]
        lo, hi = c.full_interval if metric == "full-span" else c.interior_interval
        result[c.target_id].update(range(lo, hi))
    return {t: len(values) for t, values in result.items()}


def independent_objective(cat, assignments, profile, metric="full-span", goal=0.9):
    # Direct discrete set union, no production interval/objective helpers.
    cov = covered(cat, assignments, metric)
    fractions = [
        cov[t.id] / t.reference_length if t.reference_length else 0
        for t in sorted(cat.targets, key=lambda t: t.id)
    ]
    deficits = [max(0, goal - f) / goal if goal else 0 for f in fractions]
    pools = [set() for _ in range(profile.n_pools)]
    distance = 0
    for a in assignments:
        c = cat.candidate_by_id[a.candidate_id]
        pools[a.pool].update(c.forward_oligos + c.reverse_oligos)
        distance += abs(c.full_interval[1] - c.full_interval[0] - profile.amplicon_size)
    burden = sum(len(pool) for pool in pools)
    imbalance = max(map(len, pools), default=0) - min(map(len, pools), default=0)
    return (
        max(deficits, default=0),
        sum(deficits),
        -sum(fractions) / len(fractions) if fractions else 0,
        burden,
        len(assignments),
        imbalance,
        distance,
    )


def test_one_start_actual_remove_refill_escapes_literal_260_trap():
    cat = graph_catalog([("A", "t", 20, 280), ("B", "t", 0, 150), ("C", "t", 150, 300)])
    edges = [("A", "B"), ("A", "C")]
    greedy = run(cat, edges, options=SearchOptions(starts=1, repair_rounds=0))
    repaired = run(cat, edges, options=SearchOptions(starts=1, repair_rounds=1))
    assert covered(cat, greedy.assignments) == {"t": 260}
    assert covered(cat, repaired.assignments) == {"t": 300}
    assert repaired.metadata["repairs_accepted"] >= 1
    assert repaired.metadata["completed_starts"] == 1


@pytest.mark.parametrize(
    "pools,cap,per_target,metric",
    [(1, None, None, "full-span"), (2, 3, 2, "full-span"), (3, 2, 1, "primer-trimmed")],
)
def test_tiny_exhaustive_objective_uses_independent_assignment_enumeration(
    pools, cap, per_target, metric
):
    cat = graph_catalog(
        [("A", "x", 0, 6), ("B", "x", 4, 10), ("C", "y", 0, 7), ("D", "y", 5, 10)],
        {"x": 10, "y": 10},
    )
    edges = {frozenset(pair) for pair in [("A", "B"), ("A", "C"), ("B", "D")]}
    p = limits(pools, cap, per_target)
    feasible = []
    for labels in product(range(-1, pools), repeat=4):
        assignments = tuple(
            Assignment(c.id, pool)
            for c, pool in zip(cat.candidates, labels, strict=False)
            if pool >= 0
        )
        if cap is not None and len(assignments) > cap:
            continue
        if per_target is not None and any(
            sum(
                cat.candidate_by_id[a.candidate_id].target_id == t.id
                for a in assignments
            )
            > per_target
            for t in cat.targets
        ):
            continue
        if any(
            a.pool == b.pool and frozenset((a.candidate_id, b.candidate_id)) in edges
            for i, a in enumerate(assignments)
            for b in assignments[i + 1 :]
        ):
            continue
        feasible.append(independent_objective(cat, assignments, p, metric))
    result = run(
        cat,
        edges,
        options=SearchOptions(metric=metric, starts=4, repair_rounds=2),
        profile=p,
    )
    assert independent_objective(cat, result.assignments, p, metric) == min(feasible)
    assert tuple(result.metadata["final_objective"]) == min(feasible)


def test_three_pool_repool_preserves_displaced_coverage():
    cat = graph_catalog(
        [
            ("A", "a", 0, 100),
            ("B", "b", 0, 100),
            ("C", "c", 0, 100),
            ("D", "d", 0, 100),
            ("E", "e", 0, 100),
        ],
        {t: 100 for t in "abcde"},
    )
    # D cannot enter any baseline pool. Move A to pool1 (B compatible),
    # placing D in pool0; moving C beside A/B then permits E in pool2.
    edges = [("D", x) for x in "ABC"] + [("E", x) for x in "ABCD"]
    base = tuple(Assignment(x, p) for x, p in zip("ABC", range(3), strict=False))
    result = run(
        cat,
        edges,
        profile=limits(3),
        baseline=base,
        options=SearchOptions(starts=1, repair_rounds=1),
    )
    assert sum(covered(cat, result.assignments).values()) == 500
    assert len({a.candidate_id for a in result.assignments}) == 5
    assert len({a.pool for a in result.assignments if a.candidate_id in "ABC"}) == 1


def test_seed_and_permutations_preserve_assignments_and_work_counts():
    cat = graph_catalog(
        [(str(i), "t", i * 10, i * 10 + 50) for i in range(10)], {"t": 140}
    )
    edges = [
        (str(i), str(j)) for i in range(10) for j in range(i + 1, 10) if abs(i - j) < 4
    ]
    options = SearchOptions(seed=902, starts=4, repair_rounds=2)
    first = run(cat, edges, options=options, profile=limits(3))
    shuffled = list(cat.candidates)
    random.Random(44).shuffle(shuffled)
    again = run(
        replace(cat, candidates=tuple(shuffled), targets=tuple(reversed(cat.targets))),
        reversed(edges),
        options=options,
        profile=limits(3),
    )
    assert first.assignments == again.assignments
    for key in (
        "final_objective",
        "objective_history",
        "completed_starts",
        "work",
        "repairs_accepted",
        "stop_reason",
    ):
        assert first.metadata[key] == again.metadata[key]


@pytest.mark.parametrize(
    "cap,per_target,expected", [(1, None, 1), (None, 1, 2), (0, None, 0), (None, 0, 0)]
)
def test_caps_remain_enforced_during_repairs(cap, per_target, expected):
    cat = graph_catalog(
        [
            ("A", "x", 0, 50),
            ("B", "x", 50, 100),
            ("C", "y", 0, 50),
            ("D", "y", 50, 100),
        ],
        {"x": 100, "y": 100},
    )
    result = run(cat, profile=limits(3, cap, per_target))
    assert len(result.assignments) == expected
    assert len({a.candidate_id for a in result.assignments}) == expected


def test_invalid_baseline_and_unreachable_or_empty_targets():
    cat = graph_catalog([("A", "t", 0, 100)])
    result = run(cat, invalid=("A",), baseline=(Assignment("A", 0),))
    assert result.assignments == ()
    assert result.metadata["baseline"]["supplied_valid"] is False
    assert result.metadata["baseline"]["violations"]
    assert result.metadata["per_target"]["t"]["covered_bases"] == 0
    empty = run(graph_catalog([], {"empty": 0}))
    assert empty.assignments == ()
    assert empty.metadata["final_objective"][:3] == [1.0, 1.0, 0.0]


def test_sparse_lazy_checks_and_full_catalogue_advancing_frontier():
    cat = graph_catalog([(f"{i:04}", "t", 0, 500) for i in range(500)], {"t": 500})
    oracle = GraphOracle(cat, invalid={c.id for c in cat.candidates[:-1]})
    result = search_assignments(
        cat, limits(), SearchOptions(starts=1, repair_rounds=0), oracle
    )
    assert covered(cat, result.assignments) == {"t": 500}
    assert oracle.counters["pair_calls"] < 500
    assert result.metadata["catalogue_candidates"] == 500


def test_timeout_returns_feasible_incumbent_after_completed_work_unit():
    cat = graph_catalog(
        [("A", "t", 0, 100), ("B", "t", 100, 200), ("C", "t", 200, 300)]
    )
    ticks = iter([0] * 12 + [1000] * 10000)
    result = run(
        cat, options=SearchOptions(starts=4, time_limit=1), clock=lambda: next(ticks)
    )
    assert result.metadata["stop_reason"] == "time-limit"
    assert len({a.candidate_id for a in result.assignments}) == len(result.assignments)
    assert result.metadata["validation"]["valid"]
    assert result.metadata["validation"]["scope"] == "abstract-oracle"


def test_options_reject_nonfinite_or_invalid_budgets():
    for kwargs in (
        {"time_limit": 0},
        {"time_limit": float("nan")},
        {"coverage_target": 1.1},
        {"coverage_target": float("nan")},
        {"starts": 0},
        {"repair_rounds": -1},
        {"metric": "unknown"},
    ):
        with pytest.raises(ValueError):
            SearchOptions(**kwargs)


def test_repair_removes_now_redundant_selected_candidate_without_replacement():
    cat = graph_catalog([("A", "t", 20, 280), ("B", "t", 0, 150), ("C", "t", 150, 300)])
    greedy = run(
        cat, profile=limits(3), options=SearchOptions(starts=1, repair_rounds=0)
    )
    repaired = run(
        cat, profile=limits(3), options=SearchOptions(starts=1, repair_rounds=1)
    )
    assert covered(cat, greedy.assignments) == {"t": 300}
    assert len(greedy.assignments) == 3
    assert covered(cat, repaired.assignments) == {"t": 300}
    assert {a.candidate_id for a in repaired.assignments} == {"B", "C"}


def real_trap():
    # Seed-selected random synthetic DNA, fixed independently of search. All
    # six 23-mers satisfy the unmodified production chemistry and specificity.
    from primalscheme3.core.seq_functions import reverse_complement

    rng = random.Random(1127)
    seq = "".join(rng.choice("ACGT") for _ in range(300))
    t = Target(
        "synthetic",
        0,
        0,
        ("row",),
        (tuple(seq),),
        seq,
        300,
        tuple(range(300)),
        tuple(range(301)),
    )
    cs = tuple(
        Candidate(
            ident,
            t.id,
            (a, b),
            (a + 23, b - 23),
            (seq[a : a + 23],),
            (reverse_complement(seq[b - 23 : b]),),
            (),
            (),
            (),
        )
        for ident, a, b in [("A", 20, 280), ("B", 0, 150), ("C", 150, 300)]
    )
    cat = Catalog((t,), cs, "synthetic-trap", "{}", (), ())
    p = replace(
        limits(),
        amplicon_size_min=150,
        amplicon_size_max=280,
        mismatch_kmersize=20,
        mismatch_fuzzy=False,
    )
    return cat, p


def test_production_oracle_and_fresh_validation_repair_real_synthetic_260_to_300():
    cat, p = real_trap()
    greedy = optimize_catalog(cat, p, SearchOptions(starts=1, repair_rounds=0))
    result = optimize_catalog(cat, p, SearchOptions(starts=1, repair_rounds=1))
    assert covered(cat, greedy.assignments) == {"synthetic": 260}
    assert covered(cat, result.assignments) == {"synthetic": 300}
    assert result.metadata["validation"]["valid"]
    assert (
        result.metadata["validation"]["schemaVersion"]
        == "primalscheme3.panel-validation/v1"
    )
    assert (
        result.metadata["validation"]["per_target"]["synthetic"]["full"][
            "covered_bases"
        ]
        == 300
    )
    assert result.metadata["oracle_counters"]["thermo_calls"] > 0
    assert result.metadata["validation"]["counters"]["thermo_calls"] > 0
    assert result.metadata["validation"]["allowed_secondary_products"]


def test_production_invalid_legacy_baseline_is_rejected_with_scientific_reasons():
    cat, p = real_trap()
    result = optimize_catalog(
        cat,
        p,
        SearchOptions(starts=1, repair_rounds=1),
        (Assignment("A", 0), Assignment("B", 0)),
    )
    assert result.metadata["baseline"]["supplied_valid"] is False
    assert result.metadata["baseline"]["violations"]
    assert result.metadata["validation"]["valid"]
    assert covered(cat, result.assignments) == {"synthetic": 300}


def test_production_timeout_preserves_independently_valid_fixed_baseline():
    cat, p = real_trap()
    ticks = iter([0] * 6 + [1000] * 10000)
    result = optimize_catalog(
        cat,
        p,
        SearchOptions(time_limit=1),
        (Assignment("B", 0),),
        clock=lambda: next(ticks),
    )
    assert result.assignments == (Assignment("B", 0),)
    assert result.metadata["stop_reason"] == "time-limit"
    assert result.metadata["validation"]["valid"]
    assert result.metadata["validation"]["support_diagnostics"]["B"]["joint_rows"] == [
        "row"
    ]
    assert set(result.metadata["phase_seconds"]) == {
        "baseline_validation",
        "search",
        "final_validation",
    }


def test_fresh_validator_rejects_poisoned_search_and_retains_valid_baseline(
    monkeypatch,
):
    from primalscheme3.panel.coverage_validation import CompatibilityOracle

    class CorruptPairCache(CompatibilityOracle):
        def conflict(self, first, second):
            return False

    monkeypatch.setattr(
        "primalscheme3.panel.coverage_search.CompatibilityOracle", CorruptPairCache
    )
    cat, p = real_trap()
    result = optimize_catalog(
        cat, p, SearchOptions(starts=1, repair_rounds=0), (Assignment("B", 0),)
    )
    assert result.metadata["rejected_final_validation"]["valid"] is False
    assert result.metadata["validation"]["valid"]
    assert result.assignments == (Assignment("B", 0),)
    assert result.metadata["per_target"]["synthetic"]["covered_bases"] == 150


def test_repair_budget_advances_after_improvement_to_other_deficient_targets():
    names = [f"t{i:02}" for i in range(30)]
    cat = graph_catalog(
        [
            (t + c, t, a, b)
            for t in names
            for c, a, b in [("A", 20, 280), ("B", 0, 150), ("C", 150, 300)]
        ],
        {t: 300 for t in names},
    )
    edges = [(t + "A", t + c) for t in names for c in "BC"]
    result = run(cat, edges, options=SearchOptions(starts=1, repair_rounds=1))
    assert sum(covered(cat, result.assignments).values()) == 9000


def test_two_candidate_blocking_neighborhood_replaces_both_selected_blockers():
    cat = graph_catalog(
        [
            ("A", "t", 0, 180),
            ("B", "t", 120, 300),
            ("C", "t", 0, 150),
            ("D", "t", 150, 300),
        ]
    )
    edges = [(old, new) for old in "AB" for new in "CD"]
    greedy = run(cat, edges, options=SearchOptions(starts=1, repair_rounds=0))
    repaired = run(cat, edges, options=SearchOptions(starts=1, repair_rounds=1))
    assert {a.candidate_id for a in greedy.assignments} == {"A", "B"}
    assert {a.candidate_id for a in repaired.assignments} == {"C", "D"}
    assert covered(cat, repaired.assignments) == {"t": 300}
    assert repaired.metadata["final_objective"][-1] == 0


def test_baseline_pool_label_permutation_is_semantically_canonical():
    cat = graph_catalog([("A", "x", 0, 100), ("B", "y", 0, 100)], {"x": 100, "y": 100})
    options = SearchOptions(starts=1, repair_rounds=0)
    first = run(
        cat,
        [("A", "B")],
        profile=limits(3),
        baseline=(Assignment("A", 2), Assignment("B", 0)),
        options=options,
    )
    second = run(
        replace(
            cat,
            targets=tuple(reversed(cat.targets)),
            candidates=tuple(reversed(cat.candidates)),
        ),
        [("B", "A")],
        profile=limits(3),
        baseline=(Assignment("B", 2), Assignment("A", 1)),
        options=options,
    )
    assert (
        first.assignments
        == second.assignments
        == (Assignment("A", 0), Assignment("B", 1))
    )
