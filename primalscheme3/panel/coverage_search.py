"""Bounded deterministic selection over a complete immutable catalogue.

This is a heuristic, not an optimality certificate. Candidate queues are never
pruned, but each construction/repair has explicit finite work limits. A wall
cutoff is checked between oracle predicates and search moves; a single chemistry
or specificity predicate cannot be interrupted. Discovery, baseline validation
and fresh final scientific validation are outside the search time budget.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from collections.abc import Callable
from dataclasses import asdict, dataclass
from itertools import combinations, product
from time import monotonic
from typing import Protocol

from .coverage_types import Assignment, Catalog
from .coverage_validation import (
    CompatibilityOracle,
    ConstraintProfile,
    validate_assignments,
)


@dataclass(frozen=True)
class SearchOptions:
    metric: str = "full-span"
    coverage_target: float = 0.9
    seed: int = 0
    starts: int = 4
    repair_rounds: int = 2
    time_limit: float = 120

    def __post_init__(self):
        if self.metric not in {"full-span", "primer-trimmed"}:
            raise ValueError("unknown coverage metric")
        if (
            not math.isfinite(self.coverage_target)
            or not 0 <= self.coverage_target <= 1
        ):
            raise ValueError("coverage_target must be finite and between zero and one")
        if not math.isfinite(self.time_limit) or self.time_limit <= 0:
            raise ValueError("time_limit must be positive and finite")
        for name, value, minimum in (
            ("starts", self.starts, 1),
            ("repair_rounds", self.repair_rounds, 0),
        ):
            if type(value) is not int or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if type(self.seed) is not int:
            raise ValueError("seed must be an integer")


@dataclass(frozen=True)
class SearchResult:
    assignments: tuple[Assignment, ...]
    metadata: dict


class SelectionOracle(Protocol):
    """Fixed unary/pair predicates for the ordinary combinatorial search core."""

    def candidate_valid(self, candidate_id: str) -> bool: ...
    def conflict(self, first: str, second: str) -> bool: ...


LIMITS = {
    "frontier_candidates": 64,
    "construction_candidate_attempts": 2048,
    "repair_candidate_probes_per_round": 128,
    "repair_neighborhoods_per_round": 256,
    "repair_trials_per_round": 256,
    "repair_removed_candidates": 2,
    "pool_lookahead_candidates": 4,
    "cleanup_moves_per_round": 64,
}
OBJECTIVE_TERMS = [
    "worst_normalized_shortfall",
    "summed_normalized_shortfall",
    "negative_mean_unique_coverage_fraction",
    "sum_per_pool_unique_oligos",
    "amplicon_count",
    "max_minus_min_pool_unique_oligos",
    "summed_absolute_requested_length_deviation",
]


def _merge(intervals):
    merged = []
    for lo, hi in sorted(intervals):
        if hi <= lo:
            continue
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    return tuple(merged)


def _length(intervals):
    return sum(hi - lo for lo, hi in intervals)


def _canonical(assignments):
    pools = {}
    for assignment in assignments:
        pools.setdefault(assignment.pool, []).append(assignment.candidate_id)
    groups = sorted(tuple(sorted(ids)) for ids in pools.values())
    return tuple(
        sorted(
            (
                Assignment(ident, pool)
                for pool, ids in enumerate(groups)
                for ident in ids
            ),
            key=lambda a: a.candidate_id,
        )
    )


def _signature(assignments):
    return tuple((a.candidate_id, a.pool) for a in _canonical(assignments))


class _State:
    """Coverage and unique oligos are rebuilt on removals, incremented on adds."""

    def __init__(self, catalog, profile, options, assignments=()):
        self.catalog, self.profile, self.options = catalog, profile, options
        self.assignments = {}
        self.pools = [set() for _ in range(profile.n_pools)]
        self.oligos = [set() for _ in range(profile.n_pools)]
        self.counts = Counter()
        self.coverage = {t.id: () for t in catalog.targets}
        self.distance = 0
        for assignment in assignments:
            self.add(catalog.candidate_by_id[assignment.candidate_id], assignment.pool)

    def interval(self, candidate):
        return (
            candidate.full_interval
            if self.options.metric == "full-span"
            else candidate.interior_interval
        )

    def gain(self, candidate):
        intervals = self.coverage[candidate.target_id]
        return _length(_merge((*intervals, self.interval(candidate)))) - _length(
            intervals
        )

    def add(self, candidate, pool):
        self.assignments[candidate.id] = pool
        self.pools[pool].add(candidate.id)
        self.oligos[pool].update(candidate.forward_oligos + candidate.reverse_oligos)
        self.counts[candidate.target_id] += 1
        self.coverage[candidate.target_id] = _merge(
            (*self.coverage[candidate.target_id], self.interval(candidate))
        )
        self.distance += abs(
            candidate.full_interval[1]
            - candidate.full_interval[0]
            - self.profile.amplicon_size
        )

    def vector(self):
        fractions = [
            self.fraction(t.id)
            for t in sorted(self.catalog.targets, key=lambda t: t.id)
        ]
        deficits = [
            (
                max(0.0, self.options.coverage_target - f)
                / self.options.coverage_target
                if self.options.coverage_target
                else 0.0
            )
            for f in fractions
        ]
        burdens = [len(pool) for pool in self.oligos]
        return (
            max(deficits, default=0.0),
            sum(deficits),
            -sum(fractions) / len(fractions) if fractions else 0.0,
            sum(burdens),
            len(self.assignments),
            max(burdens, default=0) - min(burdens, default=0),
            self.distance,
        )

    def fraction(self, target_id):
        length = self.catalog.target_by_id[target_id].reference_length
        return _length(self.coverage[target_id]) / length if length else 0.0

    def items(self):
        return tuple(
            Assignment(ident, pool) for ident, pool in sorted(self.assignments.items())
        )

    def key(self):
        return (*self.vector(), _signature(self.items()))

    def per_target(self):
        return {
            t.id: {
                "reference_length": t.reference_length,
                "selected_count": self.counts[t.id],
                "intervals": [list(interval) for interval in self.coverage[t.id]],
                "covered_bases": _length(self.coverage[t.id]),
                "coverage_fraction": self.fraction(t.id),
                "normalized_shortfall": (
                    max(0.0, self.options.coverage_target - self.fraction(t.id))
                    / self.options.coverage_target
                    if self.options.coverage_target
                    else 0.0
                ),
            }
            for t in sorted(self.catalog.targets, key=lambda t: t.id)
        }


class _TimeLimit(Exception):
    pass


class _PhaseLimit(Exception):
    pass


class _Cancelled(Exception):
    pass


class _Search:
    def __init__(
        self,
        catalog,
        profile,
        options,
        oracle,
        clock,
        *,
        state_factory=None,
        candidate_source=None,
        candidate_filter=None,
        aggregate_guard=None,
        record_hook=None,
        neighborhood_hook=None,
        cancelled=None,
        limits=None,
    ):
        self.catalog, self.profile, self.options = catalog, profile, options
        self.oracle, self.clock = oracle, clock
        self.state_factory = state_factory
        self.candidate_source = candidate_source
        self.candidate_filter = candidate_filter
        self.aggregate_guard = aggregate_guard
        self.record_hook = record_hook
        self.neighborhood_hook = neighborhood_hook
        self.cancelled = cancelled
        self.limits = dict(LIMITS) | (limits or {})
        self.ids = tuple(sorted(catalog.candidate_by_id))
        self.valid_cache = {}
        self.conflict_cache = {}
        self.work = Counter(
            {
                k: 0
                for k in (
                    "candidate_evaluations",
                    "pair_evaluations",
                    "candidate_attempts",
                    "frontier_evaluations",
                    "constructions",
                    "construction_limit_hits",
                    "repair_candidate_probes",
                    "repair_neighborhoods",
                    "repair_trials",
                    "repair_limit_hits",
                    "cleanup_moves",
                )
            }
        )
        self.deadline = math.inf
        self.phase_deadline = math.inf
        self.best = ()
        self.best_key = self.state().key()
        self.history = [{"phase": "empty", "objective": list(self.best_key[:-1])}]
        self.repairs_accepted = 0

    def state(self, assignments=()):
        if self.state_factory is not None:
            return self.state_factory(assignments)
        return _State(self.catalog, self.profile, self.options, assignments)

    def refresh(self, state, phase):
        if self.candidate_source is None:
            return False
        self.tick()
        changed = self.candidate_source(self, state, phase)
        self.ids = tuple(sorted(self.catalog.candidate_by_id))
        return changed

    def tick(self):
        if self.deadline < math.inf and self.cancelled is not None and self.cancelled():
            raise _Cancelled
        now = self.clock()
        if now >= self.deadline:
            raise _TimeLimit
        if now >= self.phase_deadline:
            raise _PhaseLimit

    def valid(self, ident):
        if ident not in self.valid_cache:
            self.tick()
            self.work["candidate_evaluations"] += 1
            self.valid_cache[ident] = bool(self.oracle.candidate_valid(ident))
        return self.valid_cache[ident]

    def conflict(self, first, second):
        key = tuple(sorted((first, second)))
        if key not in self.conflict_cache:
            self.tick()
            self.work["pair_evaluations"] += 1
            self.conflict_cache[key] = bool(self.oracle.conflict(*key))
        return self.conflict_cache[key]

    def record(self, state, phase):
        key = state.key()
        if key < self.best_key:
            previous = self.best
            self.best = _canonical(state.items())
            self.best_key = key
            self.history.append({"phase": phase, "objective": list(key[:-1])})
            if self.record_hook is not None:
                self.record_hook(previous, self.best, phase)
            return True
        return False

    def capped(self, state, candidate):
        p = self.profile
        return (
            p.max_amplicons is not None and len(state.assignments) >= p.max_amplicons
        ) or (
            p.max_amplicons_msa is not None
            and state.counts[candidate.target_id] >= p.max_amplicons_msa
        )

    def fits(self, state, candidate, pool):
        return (
            candidate.id not in state.assignments
            and not self.capped(state, candidate)
            and self.valid(candidate.id)
            and not any(
                self.conflict(candidate.id, other)
                for other in sorted(state.pools[pool])
            )
            and (
                self.aggregate_guard is None
                or self.aggregate_guard(state, candidate, pool)
            )
        )

    def check(self, assignments):
        state = self.state()
        violations = []
        for assignment in assignments:
            ident, pool = assignment.candidate_id, assignment.pool
            if ident not in self.catalog.candidate_by_id:
                violations.append(
                    {"candidate_id": ident, "reason": "unknown_candidate"}
                )
            elif type(pool) is not int or not 0 <= pool < self.profile.n_pools:
                violations.append({"candidate_id": ident, "reason": "invalid_pool"})
            elif ident in state.assignments:
                violations.append(
                    {"candidate_id": ident, "reason": "duplicate_candidate"}
                )
            else:
                candidate = self.catalog.candidate_by_id[ident]
                if self.capped(state, candidate):
                    violations.append({"candidate_id": ident, "reason": "count_cap"})
                if not self.valid(ident):
                    violations.append({"candidate_id": ident, "reason": "intrinsic"})
                for other in sorted(state.pools[pool]):
                    if self.conflict(ident, other):
                        violations.append(
                            {
                                "candidate_id": ident,
                                "other_candidate_id": other,
                                "pool": pool,
                                "reason": "pair_conflict",
                            }
                        )
                if self.aggregate_guard is not None and not self.aggregate_guard(
                    state, candidate, pool
                ):
                    violations.append(
                        {
                            "candidate_id": ident,
                            "pool": pool,
                            "reason": "aggregate-pool",
                        }
                    )
                state.add(candidate, pool)
        return {
            "scope": "abstract-oracle",
            "valid": not violations,
            "violations": violations,
        }

    def queues(self, state, start, attempted=()):
        queues = {t.id: [] for t in self.catalog.targets}
        rng = random.Random(self.options.seed + start)
        priority = {ident: rng.random() for ident in self.ids}
        for ident in self.ids:
            if (
                ident not in state.assignments
                and ident not in attempted
                and self.valid_cache.get(ident) is not False
                and (self.candidate_filter is None or self.candidate_filter(ident))
            ):
                queues[self.catalog.candidate_by_id[ident].target_id].append(ident)
        for values in queues.values():
            if start < 2:
                values.sort(
                    key=lambda ident: (
                        (
                            -state.queue_weight(self.catalog.candidate_by_id[ident])
                            if hasattr(state, "queue_weight")
                            else -_length(
                                (state.interval(self.catalog.candidate_by_id[ident]),)
                            )
                        ),
                        ident,
                    )
                )
            else:
                values.sort(key=lambda ident: (priority[ident], ident))
        return queues, priority

    def fill(self, state, start, phase):
        self.tick()
        self.work["constructions"] += 1
        # Built-in construction states only add assignments, so a consumed
        # zero-gain, capped, or incompatible candidate cannot become feasible
        # until a later fill starts from a different state, including repair.
        attempted = set()
        self.refresh(state, phase)
        queues, priority = self.queues(state, start, attempted)
        positions = {target: 0 for target in queues}
        for _ in range(self.limits["construction_candidate_attempts"]):
            self.tick()
            active = [
                target
                for target, values in queues.items()
                if positions[target] < len(values)
            ]
            if not active:
                if self.refresh(state, phase):
                    queues, priority = self.queues(state, start, attempted)
                    positions = {target: 0 for target in queues}
                    continue
                return "exhausted"
            active.sort(
                key=lambda target: (
                    state.fraction(target),
                    len(queues[target]) - positions[target],
                    target,
                )
            )
            # One head per target first, then deeper items: scarce loci do not
            # lose their frontier to a large catalogue from a single target.
            frontier = []
            depth = 0
            while len(frontier) < self.limits["frontier_candidates"]:
                before = len(frontier)
                for target in active:
                    index = positions[target] + depth
                    if index < len(queues[target]):
                        frontier.append(queues[target][index])
                        if len(frontier) == self.limits["frontier_candidates"]:
                            break
                if len(frontier) == before:
                    break
                depth += 1

            def rank(ident, queues=queues, positions=positions, priority=priority):
                candidate = self.catalog.candidate_by_id[ident]
                gain = state.gain(candidate)
                fraction = state.fraction(candidate.target_id)
                deficit = max(0.0, self.options.coverage_target - fraction)
                scarcity = (
                    len(queues[candidate.target_id]) - positions[candidate.target_id]
                )
                if start == 0:
                    return (-gain, -deficit, scarcity, ident)
                return (
                    -deficit,
                    scarcity,
                    -gain,
                    priority[ident] if start > 1 else 0.0,
                    ident,
                )

            self.work["frontier_evaluations"] += len(frontier)
            ident = min(frontier, key=rank)
            candidate = self.catalog.candidate_by_id[ident]
            values = queues[candidate.target_id]
            # Move consumed frontier entries to the prefix without deleting
            # or truncating the rest of the complete target queue.
            index = values.index(ident, positions[candidate.target_id])
            values[index], values[positions[candidate.target_id]] = (
                values[positions[candidate.target_id]],
                values[index],
            )
            positions[candidate.target_id] += 1
            self.work["candidate_attempts"] += 1
            attempted.add(ident)
            if (
                state.gain(candidate) <= 0
                or self.capped(state, candidate)
                or not self.valid(ident)
            ):
                continue
            feasible = [
                pool
                for pool in range(self.profile.n_pools)
                if self.fits(state, candidate, pool)
            ]
            if not feasible:
                continue
            # Preserve compatibility for a bounded frontier lookahead without
            # materializing a dense candidate conflict graph.
            lookahead = [
                other
                for other in frontier
                if other != ident and other not in state.assignments
            ][: self.limits["pool_lookahead_candidates"]]
            oligos = set(candidate.forward_oligos + candidate.reverse_oligos)

            def pool_rank(pool, ident=ident, lookahead=lookahead, oligos=oligos):
                lost = sum(
                    not any(
                        self.conflict(other, selected)
                        for selected in sorted(state.pools[pool])
                    )
                    and self.conflict(other, ident)
                    for other in lookahead
                )
                return (
                    len(oligos - state.oligos[pool]),
                    lost,
                    len(state.oligos[pool] | oligos),
                    len(state.pools[pool]),
                    pool,
                )

            pool = min(feasible, key=pool_rank)
            state.add(candidate, pool)
            self.record(state, phase)
        self.work["construction_limit_hits"] += 1

        return "work-cap"

    def cleanup(self, phase):
        moves = 0
        for ident in sorted(a.candidate_id for a in self.best):
            for destination in [-1, *range(self.profile.n_pools)]:
                self.tick()
                if moves >= self.limits["cleanup_moves_per_round"]:
                    return "work-cap"
                current = self.state(self.best)
                if ident not in current.assignments:
                    break
                if destination == current.assignments[ident]:
                    continue
                moves += 1
                self.work["cleanup_moves"] += 1
                trial = self.state(
                    a for a in current.items() if a.candidate_id != ident
                )
                if destination >= 0:
                    candidate = self.catalog.candidate_by_id[ident]
                    if not self.fits(trial, candidate, destination):
                        continue
                    trial.add(candidate, destination)
                if self.record(trial, phase):
                    self.repairs_accepted += 1

        return "exhausted" if moves else "no-eligible-work"

    def prepare_repair(self, start, round_index):
        if self.neighborhood_hook is not None:
            self.tick()
            try:
                return (
                    self.neighborhood_hook(
                        self, self.state(self.best), start, round_index
                    )
                    or "completed"
                )
            finally:
                self.ids = tuple(sorted(self.catalog.candidate_by_id))
        return "completed"

    def repair(self, start, round_index):
        self.prepare_repair(start, round_index)
        self.cleanup(f"cleanup:{start}:{round_index}")
        return self.exchange(start, round_index)

    def exchange(self, start, round_index):
        initial = self.state(self.best)
        # Consider deficits first, but also fully covered targets: burden,
        # count and pool balance remain objective terms after reaching goal.
        candidates = sorted(
            (
                ident
                for ident in self.ids
                if ident not in initial.assignments
                and self.valid_cache.get(ident) is not False
                and (self.candidate_filter is None or self.candidate_filter(ident))
            ),
            key=lambda ident: (
                initial.fraction(self.catalog.candidate_by_id[ident].target_id),
                -initial.gain(self.catalog.candidate_by_id[ident]),
                ident,
            ),
        )
        probes = neighborhoods = trials = 0
        for ident in candidates:
            self.tick()
            if any(a.candidate_id == ident for a in self.best):
                continue
            if probes >= self.limits["repair_candidate_probes_per_round"]:
                self.work["repair_limit_hits"] += 1
                return "work-cap"
            probes += 1
            self.work["repair_candidate_probes"] += 1
            if not self.valid(ident):
                continue
            candidate = self.catalog.candidate_by_id[ident]
            improved = False
            for pool in range(self.profile.n_pools):
                current = self.state(self.best)
                if ident in current.assignments:
                    break
                mandatory = tuple(
                    other
                    for other in sorted(current.pools[pool])
                    if self.conflict(ident, other)
                )
                if len(mandatory) > self.limits["repair_removed_candidates"]:
                    continue
                optional = sorted(set(current.assignments) - set(mandatory))
                for size in range(
                    self.limits["repair_removed_candidates"] - len(mandatory) + 1
                ):
                    for extra in combinations(optional, size):
                        self.tick()
                        if (
                            neighborhoods
                            >= self.limits["repair_neighborhoods_per_round"]
                            or trials >= self.limits["repair_trials_per_round"]
                        ):
                            self.work["repair_limit_hits"] += 1
                            return "work-cap"
                        neighborhoods += 1
                        self.work["repair_neighborhoods"] += 1
                        removed = tuple(sorted((*mandatory, *extra)))
                        retained = tuple(
                            a for a in current.items() if a.candidate_id not in removed
                        )
                        base = self.state(retained)
                        if not self.fits(base, candidate, pool):
                            continue
                        # Reinsert displaced candidates in another pool or
                        # drop them; two-candidate neighborhoods permit a
                        # second relocation blocker to move as well.
                        destinations = [
                            p for p in range(self.profile.n_pools) if p != pool
                        ] + [-1]
                        for locations in product(destinations, repeat=len(removed)):
                            self.tick()
                            if trials >= self.limits["repair_trials_per_round"]:
                                self.work["repair_limit_hits"] += 1
                                return "work-cap"
                            trials += 1
                            self.work["repair_trials"] += 1
                            trial = self.state(retained)
                            trial.add(candidate, pool)
                            valid = True
                            for other, destination in zip(
                                removed, locations, strict=False
                            ):
                                if destination < 0:
                                    continue
                                displaced = self.catalog.candidate_by_id[other]
                                if not self.fits(trial, displaced, destination):
                                    valid = False
                                    break
                                trial.add(displaced, destination)
                            if not valid:
                                continue
                            before = self.best_key
                            phase = f"repair:{start}:{round_index}"
                            self.record(trial, phase)
                            self.fill(trial, start, phase)
                            if self.best_key < before:
                                self.repairs_accepted += 1
                                improved = True
                                break
                        if improved:
                            break
                    if improved:
                        break
                if improved:
                    break

        return "exhausted" if probes else "no-eligible-work"


def search_assignments(
    catalog: Catalog,
    profile: ConstraintProfile,
    options: SearchOptions,
    oracle: SelectionOracle,
    baseline=(),
    *,
    clock: Callable[[], float] = monotonic,
) -> SearchResult:
    """Solve a fixed candidate/constraint instance with an abstract oracle.

    This general combinatorial core checks its final assignments against the
    supplied predicates; it does not assert scientific panel validity. Native
    workflows must use :func:`optimize_catalog`, whose oracle and fresh final
    validator cannot be replaced through this API.
    """
    search = _Search(catalog, profile, options, oracle, clock)
    baseline = tuple(baseline)
    baseline_begin = clock()
    report = search.check(baseline)
    if report["valid"]:
        search.record(search.state(baseline), "supplied-baseline")
    baseline_vector = list(search.best_key[:-1])
    baseline_seconds = max(0.0, clock() - baseline_begin)
    begin = clock()
    search.deadline = begin + options.time_limit
    completed_starts = completed_repairs = 0
    stop_reason = "completed"
    baseline_construction = None
    try:
        for start in range(options.starts):
            state = search.state()
            search.fill(state, start, f"construction:{start}")
            completed_starts += 1
            if start == 0:
                baseline_construction = list(state.vector())
            for round_index in range(options.repair_rounds):
                search.repair(start, round_index)
                completed_repairs += 1
    except _TimeLimit:
        stop_reason = "time-limit"
    elapsed = max(0.0, clock() - begin)
    search.deadline = math.inf
    final_begin = clock()
    validation = search.check(search.best)
    final_seconds = max(0.0, clock() - final_begin)
    state = search.state(search.best)
    metadata = {
        "schemaVersion": "primalscheme3.panel-optimizer/v1",
        "algorithm": "bounded-multistart-coverage/v1",
        "objective_version": "lexicographic-coverage/v1",
        "objective_terms": OBJECTIVE_TERMS.copy(),
        "objective_direction": "minimize-lexicographically",
        "tie_break": "candidate-id order with canonical symmetric pool labels",
        "options": asdict(options),
        "catalogue_candidates": len(catalog.candidates),
        "catalogue_targets": len(catalog.targets),
        "catalog_semantic_digest": catalog.semantic_digest,
        "search_limits": dict(LIMITS),
        "completed_starts": completed_starts,
        "completed_repair_rounds": completed_repairs,
        "repairs_accepted": search.repairs_accepted,
        "work": dict(search.work),
        "stop_reason": stop_reason,
        "objective_history": search.history,
        "baseline": {
            "supplied_count": len(baseline),
            "supplied_valid": report["valid"],
            "violations": report["violations"],
            "retained_objective": baseline_vector,
            "construction": "deterministic-gain-first-under-current-rules",
            "construction_objective": baseline_construction,
        },
        "baseline_objective": (
            baseline_construction
            if baseline_construction is not None
            else baseline_vector
        ),
        "final_objective": list(state.vector()),
        "per_target": state.per_target(),
        "validation": validation,
        "phase_seconds": {
            "baseline_validation": baseline_seconds,
            "search": elapsed,
            "final_validation": final_seconds,
        },
        "limitations": [
            "Heuristic bounded search; no global optimality guarantee.",
            "Complete catalogue retained; frontier, construction and repair limits can leave candidates or moves unexamined.",
            "Fixed-input/seed/work determinism excludes runs interrupted by wall time; one oracle predicate is indivisible.",
            "Discovery and baseline/final validation are outside the search time limit.",
            "Baseline construction uses current constraints and is not a rerun of historical unseeded legacy selection.",
            "Abstract oracle checks do not establish panel-v1 scientific validity.",
        ],
    }
    return SearchResult(search.best, metadata)


def optimize_catalog(
    catalog: Catalog,
    profile: ConstraintProfile,
    options: SearchOptions,
    baseline=(),
    *,
    clock: Callable[[], float] = monotonic,
) -> SearchResult:
    """Run production compatibility and independently validate the final panel."""
    baseline = tuple(baseline)
    begin = clock()
    baseline_report = validate_assignments(
        catalog, baseline, profile, options.metric, options.coverage_target
    )
    baseline_seconds = max(0.0, clock() - begin)
    valid_baseline = baseline if baseline_report["valid"] else ()
    oracle = CompatibilityOracle(catalog, profile)
    result = search_assignments(
        catalog, profile, options, oracle, valid_baseline, clock=clock
    )
    begin = clock()
    validation = validate_assignments(
        catalog, result.assignments, profile, options.metric, options.coverage_target
    )
    metadata = result.metadata
    assignments = result.assignments
    if not validation["valid"]:
        # A fresh-validator disagreement can never publish a failing panel.
        # Retain the independently valid fixed baseline, or the feasible empty
        # incumbent, and preserve the rejected proposal's diagnostics.
        metadata["rejected_final_validation"] = validation
        assignments = _canonical(valid_baseline)
        validation = validate_assignments(
            catalog, assignments, profile, options.metric, options.coverage_target
        )
        if not validation["valid"]:
            raise ValueError("independent validation rejected the fallback baseline")
        state = _State(catalog, profile, options, assignments)
        metadata["final_objective"] = list(state.vector())
        metadata["per_target"] = state.per_target()
        metadata["validation_fallback"] = True
    metadata["baseline"]["supplied_count"] = len(baseline)
    metadata["baseline"]["supplied_valid"] = baseline_report["valid"]
    metadata["baseline"]["violations"] = baseline_report["violations"]
    metadata["validation"] = validation
    metadata["profile"] = profile.to_dict()
    metadata["oracle_counters"] = dict(oracle.counters)
    metadata["phase_seconds"]["baseline_validation"] += baseline_seconds
    metadata["phase_seconds"]["final_validation"] += max(0.0, clock() - begin)
    metadata["limitations"] = [
        item
        for item in metadata["limitations"]
        if not item.startswith("Abstract oracle")
    ]
    return SearchResult(assignments, metadata)
