"""Allele-aware adapters for the shared bounded construction/repair engine.

The abstract entry point tests combinatorial behavior. Scientific callers use
search_allele_assignments, which owns real kernels and fresh final validation.
"""

from __future__ import annotations

import json
import math
from collections import Counter, OrderedDict, defaultdict, deque
from dataclasses import asdict, dataclass, field, replace
from functools import cached_property, partial
from importlib.metadata import version
from statistics import mean
from time import monotonic

from .allele_coverage import allele_summary, configuration_coverage, coverage_utility
from .allele_validation import (
    AlleleCompatibilityOracle,
    AlleleConstraintProfile,
    StagePolicy,
    validate_allele_assignments,
)
from .coverage_history import CoverageHistory
from .coverage_search import _Cancelled, _canonical, _Search, _signature, _TimeLimit
from .coverage_types import (
    AlleleAssignment,
    Assignment,
    ConfigurationLedger,
    CoverageSummary,
    VariantCatalog,
    semantic_id,
)
from .coverage_variants import (
    ConfigurationIneligible,
    ProposalContext,
    ProposalWitness,
    SubsetLimits,
    make_configuration,
    propose_configurations,
)
from .phase_scheduling import PhaseScheduler, Reservation, scheduling_policy


@dataclass(frozen=True)
class AlleleWorkLimits:
    frontier_candidates: int = 64
    construction_candidate_attempts: int = 2048
    repair_candidate_probes_per_round: int = 128
    repair_neighborhoods_per_round: int = 256
    repair_trials_per_round: int = 256
    pool_lookahead_candidates: int = 4
    cleanup_moves_per_round: int = 64
    families_per_refresh: int = 16

    def __post_init__(self):
        for name, value in asdict(self).items():
            if type(value) is not int or value < 1:
                raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class AlleleSearchOptions:
    coverage_target: float = 0.95
    seed: int = 0
    starts: int = 4
    repair_rounds: int = 2
    time_limit: float = 120
    subset_beam_width: int = 16
    subset_expansion_limit: int = 256
    exchange_width: int = 2
    requested_amplicon_size: int | None = None
    work_limits: AlleleWorkLimits = field(default_factory=AlleleWorkLimits)
    seed_modes: tuple[str, ...] = ("full", "normal")
    phase_scheduling: str = "serial"
    variant_selection: str = "subsets"

    def __post_init__(self):
        if (
            isinstance(self.coverage_target, bool)
            or not math.isfinite(self.coverage_target)
            or not 0 < self.coverage_target <= 1
        ):
            raise ValueError("coverage_target must be finite and in (0, 1]")
        if (
            isinstance(self.time_limit, bool)
            or not math.isfinite(self.time_limit)
            or self.time_limit <= 0
        ):
            raise ValueError("time_limit must be positive and finite")
        for name, minimum in (
            ("starts", 1),
            ("repair_rounds", 0),
            ("subset_beam_width", 1),
            ("subset_expansion_limit", 1),
            ("exchange_width", 0),
        ):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if self.exchange_width > 2:
            raise ValueError("exchange_width must be <= 2")
        if type(self.seed) is not int:
            raise ValueError("seed must be an integer")
        if self.requested_amplicon_size is not None and (
            type(self.requested_amplicon_size) is not int
            or self.requested_amplicon_size <= 0
        ):
            raise ValueError("requested_amplicon_size must be a positive integer")
        if not isinstance(self.work_limits, AlleleWorkLimits):
            raise ValueError("work_limits must be AlleleWorkLimits")
        if len(set(self.seed_modes)) != len(self.seed_modes) or any(
            mode not in ("full", "normal") for mode in self.seed_modes
        ):
            raise ValueError("seed_modes may contain full and normal once each")
        if self.phase_scheduling not in ("serial", "reserved"):
            raise ValueError("phase_scheduling must be serial or reserved")
        if self.variant_selection not in ("subsets", "full-cloud"):
            raise ValueError("variant_selection must be subsets or full-cloud")
        object.__setattr__(self, "seed_modes", tuple(self.seed_modes))


@dataclass(frozen=True)
class AlleleSearchResult:
    catalog: VariantCatalog
    ledger: ConfigurationLedger
    assignments: tuple[AlleleAssignment, ...]
    coverage: CoverageSummary
    validation: dict
    objective: tuple
    metadata: dict
    history: CoverageHistory = field(repr=False, compare=False)


@dataclass(frozen=True)
class _ConfigurationView:
    configuration: object
    catalog: VariantCatalog

    @property
    def id(self):
        return self.configuration.id

    @property
    def target_id(self):
        return self.configuration.target_id

    @property
    def full_interval(self):
        return self.configuration.full_interval

    @cached_property
    def forward_oligos(self):
        return tuple(
            self.catalog.site_by_id[s].sequence
            for s in self.configuration.forward_site_ids
        )

    @cached_property
    def reverse_oligos(self):
        return tuple(
            self.catalog.site_by_id[s].sequence
            for s in self.configuration.reverse_site_ids
        )

    @cached_property
    def is_full_family(self):
        family = self.catalog.family_by_id[self.configuration.family_id]

        def eligible(ids):
            return tuple(
                sorted(
                    sid
                    for sid in ids
                    if self.catalog.site_by_id[sid].accepting_profile_ids
                    and self.catalog.site_by_id[sid].reference_footprint is not None
                    and self.catalog.site_by_id[sid].mapping_failure is None
                )
            )

        return self.configuration.forward_site_ids == eligible(
            family.forward_site_ids
        ) and self.configuration.reverse_site_ids == eligible(family.reverse_site_ids)

    @cached_property
    def contribution(self):
        interiors = configuration_coverage(self.catalog, self.configuration)
        return {
            aid: positions
            & frozenset(self.catalog.observation_by_id[aid].observed_positions)
            for aid, positions in interiors.items()
        }


class _CatalogView:
    def __init__(self, catalog):
        self.targets = catalog.targets
        self.target_by_id = catalog.target_by_id
        self.candidate_by_id = {}


class _AssessmentOracle:
    """Persist performed checks while leaving numerical/policy caches with the oracle."""

    def __init__(self, oracle, catalog, profile, policy, history):
        self.oracle, self.catalog, self.profile, self.policy, self.history = (
            oracle,
            catalog,
            profile,
            policy,
            history,
        )
        self.kernels = {
            name: version(name) for name in ("primalschemers", "primer3-py")
        }

    def __getattr__(self, name):
        return getattr(self.oracle, name)

    def _record(self, ids, scope, result):
        dependency = {
            "catalog": self.catalog.semantic_digest,
            "profile": self.profile.cache_key,
            "policy": self.policy.cache_key,
            "configuration_ids": ids,
            "kernels": self.kernels,
        }
        evidence = self.history.record_evidence(
            entity_ids=ids,
            measurement=scope,
            dependency_key=dependency,
            values=result,
            status="evaluated",
        )
        checks = result.get("checks") or {
            scope: "fail" if result["reasons"] else "pass"
        }
        assessments = []
        for check, outcome in checks.items():
            assessment = self.history.assess(
                stage_id=self.policy.stage_id,
                entity_ids=ids,
                pool=None,
                context_digest=semantic_id(scope, dependency),
                profile_id="allele-panel-v1",
                kernel_versions=self.kernels,
                thresholds=self.profile.to_dict() | asdict(self.policy),
                check_name=check,
                outcome=outcome,
                reason=(
                    "upstream-check-failed"
                    if outcome == "not-evaluated"
                    else ",".join(result["reasons"]) or "within-policy"
                ),
                evidence_ids=(evidence.id,),
            )
            assessments.append(assessment.id)
        parents = tuple(
            sorted(
                {
                    event.id
                    for cid in ids
                    if (event := self.history.latest_event(cid)) is not None
                }
            )
        )
        self.history.emit(
            stage_id=self.policy.stage_id,
            kind="assessed",
            entity_ids=ids,
            parent_event_ids=parents,
            assessment_ids=tuple(assessments),
            changes={
                "scope": scope,
                "reasons": result["reasons"],
                "context_digest": semantic_id(scope, dependency),
            },
        )

    def candidate_valid(self, cid):
        if hasattr(self.oracle, "candidate_diagnostics"):
            result = self.oracle.candidate_diagnostics(cid)
            valid = not result["reasons"]
        else:
            valid = self.oracle.candidate_valid(cid)
            result = {"reasons": [] if valid else ["abstract-candidate-rejected"]}
        self._record((cid,), "candidate-compatibility", result)
        return valid

    def conflict(self, first, second):
        if hasattr(self.oracle, "pair_diagnostics"):
            result = self.oracle.pair_diagnostics(first, second)
            conflict = bool(result["reasons"])
        else:
            conflict = self.oracle.conflict(first, second)
            result = {"reasons": ["abstract-pair-conflict"] if conflict else []}
        self._record(tuple(sorted((first, second))), "pair-compatibility", result)
        return conflict


class _PoolEvaluations:
    """Full membership is the verdict key; raw kernel caching belongs to oracle."""

    def __init__(self, oracle, catalog, profile, policy, history):
        self.oracle, self.catalog, self.profile, self.policy, self.history = (
            oracle,
            catalog,
            profile,
            policy,
            history,
        )
        self.cache, self.assessed = OrderedDict(), set()

    def diagnostics(self, ids):
        key = tuple(sorted(ids))
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        valid = self.oracle.pool_valid(key)
        if valid and self.policy.active_cutoff == self.policy.strict_cutoff:
            result = {
                "reasons": [],
                "violating_edge_count": 0,
                "incident_species_count": 0,
                "diagnostic_scope": "incremental-strict-pool-guard",
            }
        elif valid and hasattr(self.oracle, "pool_exposure_counts"):
            edges, incident = self.oracle.pool_exposure_counts(key)
            result = {
                "reasons": [],
                "violating_edge_count": edges,
                "incident_species_count": incident,
                "diagnostic_scope": "incremental-pool-exposure",
            }
        else:
            result = self.oracle.pool_diagnostics(key)
        self.cache[key] = result
        if len(self.cache) > 2048:
            self.cache.popitem(last=False)
        return result

    def allows(self, state, candidate, pool):
        ids = tuple(sorted(state.pools[pool] | {candidate.id}))
        result = self.diagnostics(ids)
        context = semantic_id("pool-context", {"pool": pool, "configuration_ids": ids})
        if context not in self.assessed:
            self.assessed.add(context)
            evidence = self.history.record_evidence(
                entity_ids=ids,
                measurement="pool-evaluation",
                dependency_key={
                    "catalog": self.catalog.semantic_digest,
                    "profile": self.profile.cache_key,
                    "policy": self.policy.cache_key,
                    "configurations": ids,
                },
                values=result,
                status="evaluated",
            )
            self.history.assess(
                stage_id=self.policy.stage_id,
                entity_ids=ids,
                pool=pool,
                context_digest=context,
                profile_id="allele-panel-v1",
                kernel_versions={},
                thresholds=asdict(self.policy),
                check_name="aggregate-pool-exposure",
                outcome="fail" if result["reasons"] else "pass",
                reason=",".join(result["reasons"]) or "within-policy",
                evidence_ids=(evidence.id,),
            )
        return not result["reasons"]


class _AlleleState:
    """Exact observed-base contribution counts; removals preserve other contributors."""

    def __init__(self, catalog, view, profile, options, pools, assignments=()):
        self.catalog, self.view, self.profile, self.options, self.pool_evaluations = (
            catalog,
            view,
            profile,
            options,
            pools,
        )
        self.assignments = {}
        self.pools = [set() for _ in range(profile.n_pools)]
        self.oligos = [set() for _ in range(profile.n_pools)]
        self.species_counts = [Counter() for _ in range(profile.n_pools)]
        self.counts = Counter()
        self.coverage = {a.id: Counter() for a in catalog.observations}
        self.by_target = defaultdict(list)
        for a in catalog.observations:
            if a.observed_positions:
                self.by_target[a.target_id].append(a)
        self.distance = 0
        requested = options.requested_amplicon_size
        if requested is None:
            requested = json.loads(catalog.resolved_config_json).get("amplicon_size")
        self.requested_size = (
            requested
            if requested is not None
            else (profile.amplicon_size_min + profile.amplicon_size_max) / 2
        )
        if (
            isinstance(self.requested_size, bool)
            or not isinstance(self.requested_size, (int, float))
            or not math.isfinite(self.requested_size)
            or self.requested_size <= 0
        ):
            raise ValueError(
                "resolved requested amplicon size must be positive and finite"
            )
        for assignment in assignments:
            self.add(view.candidate_by_id[assignment.candidate_id], assignment.pool)

    def _fraction(self, target_id, extra=None):
        alleles = self.by_target.get(target_id, ())
        if not alleles:
            return None
        return mean(
            len(
                set(self.coverage[a.id])
                | (extra.get(a.id, frozenset()) if extra else set())
            )
            / len(a.observed_positions)
            for a in alleles
        )

    def fraction(self, target_id):
        value = self._fraction(target_id)
        return value if value is not None else 0.0

    def gain(self, candidate):
        before = self._fraction(candidate.target_id)
        if before is None:
            return 0.0
        after = self._fraction(candidate.target_id, candidate.contribution)
        return (
            coverage_utility((after,), goal=self.options.coverage_target)
            - coverage_utility((before,), goal=self.options.coverage_target)
        ) / len(self.by_target)

    def queue_weight(self, candidate):
        return self.gain(candidate)

    def add(self, candidate, pool):
        if candidate.id in self.assignments:
            raise ValueError("duplicate configuration in search state")
        self.assignments[candidate.id] = pool
        self.pools[pool].add(candidate.id)
        self.counts[candidate.target_id] += 1
        for aid, positions in candidate.contribution.items():
            self.coverage[aid].update(positions)
        self.species_counts[pool].update(
            set(candidate.forward_oligos + candidate.reverse_oligos)
        )
        self.oligos[pool] = set(self.species_counts[pool])
        self.distance += abs(
            candidate.full_interval[1]
            - candidate.full_interval[0]
            - self.requested_size
        )

    def drop(self, candidate):
        pool = self.assignments.pop(candidate.id)
        self.pools[pool].remove(candidate.id)
        self.counts[candidate.target_id] -= 1
        for aid, positions in candidate.contribution.items():
            self.coverage[aid].subtract(positions)
            self.coverage[aid] += Counter()
        self.species_counts[pool].subtract(
            set(candidate.forward_oligos + candidate.reverse_oligos)
        )
        self.species_counts[pool] += Counter()
        self.oligos[pool] = set(self.species_counts[pool])
        self.distance -= abs(
            candidate.full_interval[1]
            - candidate.full_interval[0]
            - self.requested_size
        )

    def vector(self):
        fractions = tuple(self.fraction(tid) for tid in sorted(self.by_target))
        burdens = [len(x) for x in self.oligos]
        exposure = [self.pool_evaluations.diagnostics(ids) for ids in self.pools if ids]
        values = (
            -(coverage_utility(fractions, goal=self.options.coverage_target) or 0),
            -mean(fractions) if fractions else 0,
            -min(fractions) if fractions else 0,
            sum(x["violating_edge_count"] for x in exposure),
            sum(x["incident_species_count"] for x in exposure),
            sum(burdens),
            len(self.assignments),
            max(burdens, default=0) - min(burdens, default=0),
            self.distance,
        )
        return tuple(round(value, 12) for value in values)

    def items(self):
        return tuple(
            Assignment(cid, pool) for cid, pool in sorted(self.assignments.items())
        )

    def key(self):
        return (*self.vector(), _signature(self.items()))


class _Proposals:
    def __init__(self, catalog, ledger, oracle, profile, policy, options, history):
        (
            self.catalog,
            self.oracle,
            self.profile,
            self.policy,
            self.options,
            self.history,
        ) = (catalog, oracle, profile, policy, options, history)
        self.records = dict(ledger.configuration_by_id)
        self.by_family = defaultdict(dict)
        self.view = _CatalogView(catalog)
        self.seed_ids = {"full": set(), "normal": set()}
        self.mode = "expanded"
        self.cursors = Counter()
        self.searched_families = set()
        self.work = Counter()
        self.omitted = set()
        self.pending_streams = 0
        groups = defaultdict(deque)
        for family in catalog.families:
            groups[family.target_id].append(family.id)
        ordered = []
        while groups:
            for tid in sorted(tuple(groups)):
                ordered.append(groups[tid].popleft())
                if not groups[tid]:
                    del groups[tid]
        self.families = tuple(ordered)
        for config in ledger.configurations:
            self._register(config)

    def _register(self, config):
        if config.id in self.view.candidate_by_id:
            if self.records[config.id] != config:
                raise ValueError("configuration collision")
            return
        self.records[config.id] = config
        self.by_family[config.family_id][config.id] = config
        self.view.candidate_by_id[config.id] = _ConfigurationView(config, self.catalog)
        self.oracle.register(config)

    def _seed(self, family_id, mode):
        family = self.catalog.family_by_id[family_id]

        def eligible(ids):
            return tuple(
                sid
                for sid in ids
                if self.catalog.site_by_id[sid].accepting_profile_ids
                and (
                    mode != "normal"
                    or "normal" in self.catalog.site_by_id[sid].accepting_profile_ids
                )
                and self.catalog.site_by_id[sid].reference_footprint is not None
                and self.catalog.site_by_id[sid].mapping_failure is None
            )

        forward = eligible(family.forward_site_ids)
        reverse = eligible(family.reverse_site_ids)
        try:
            config = make_configuration(
                self.catalog,
                family_id,
                forward,
                reverse,
                min_size=self.profile.amplicon_size_min,
                max_size=self.profile.amplicon_size_max,
            )
        except ConfigurationIneligible as exc:
            self.work["ineligible_seeds"] += 1
            context = {
                "catalog_digest": self.catalog.semantic_digest,
                "family_id": family_id,
                "seed_mode": mode,
                "forward_site_ids": forward,
                "reverse_site_ids": reverse,
                "profile": self.profile.to_dict(),
                "stage_policy": asdict(self.policy),
            }
            attempt = semantic_id("seed-attempt/v1", context)
            entities = (attempt, family_id, *forward, *reverse)
            evidence = self.history.record_evidence(
                entity_ids=entities,
                measurement="seed-construction",
                dependency_key=context,
                values=context | {"reason": exc.reason},
                status="evaluated",
            )
            assessment = self.history.assess(
                stage_id=self.policy.stage_id,
                entity_ids=entities,
                pool=None,
                context_digest=attempt,
                profile_id="allele-panel-v1",
                kernel_versions={},
                thresholds=self.profile.to_dict() | asdict(self.policy),
                check_name="seed-construction",
                outcome="fail",
                reason=exc.reason,
                evidence_ids=(evidence.id,),
            )
            self.history.emit(
                stage_id=self.policy.stage_id,
                kind="seed-rejected",
                entity_ids=entities,
                assessment_ids=(assessment.id,),
                changes=context | {"attempt_id": attempt, "reason": exc.reason},
            )
            return
        self._register(config)
        self.seed_ids[mode].add(config.id)
        self.history.emit(
            stage_id=self.policy.stage_id,
            kind="seed-configuration",
            entity_ids=(config.id,),
            changes={"seed_mode": mode, "family_id": family_id},
        )

    def _expand(
        self, family_id, *, parents=(), witnesses=(), pool=None, context_digest=""
    ):
        local = ConfigurationLedger(
            self.catalog.semantic_digest, tuple(self.by_family[family_id].values())
        )
        batch = propose_configurations(
            self.catalog,
            local,
            family_id,
            parent_ids=parents,
            witnesses=witnesses,
            context=ProposalContext(
                pool,
                context_digest,
                self.profile.amplicon_size_min,
                self.profile.amplicon_size_max,
            ),
            stage=self.policy.stage_id,
            limits=SubsetLimits(
                self.options.subset_beam_width, self.options.subset_expansion_limit
            ),
            history=self.history,
        )
        for cid in batch.new_ids:
            self._register(batch.ledger.configuration_by_id[cid])
        self.searched_families.add(family_id)
        self.work["neighborhood_calls"] += 1
        self.work["expanded_states"] += batch.expansion_count
        self.work["materialized_configurations"] += len(batch.new_ids)
        self.work["truncated_neighborhoods"] += int(batch.truncated)
        self.omitted.update(batch.omitted_materialized_ids)
        self.pending_streams += batch.unexpanded_frontier_count

    def refresh(self, search, state, phase):
        begin = self.cursors[self.mode]
        end = min(
            len(self.families), begin + self.options.work_limits.families_per_refresh
        )
        for index in range(begin, end):
            search.tick()
            family = self.families[index]
            if self.options.variant_selection == "full-cloud":
                self._seed(family, "full")
                self.searched_families.add(family)
            elif self.mode == "expanded":
                self._expand(family)
            else:
                self._seed(family, self.mode)
            self.cursors[self.mode] = index + 1
        return end > begin

    def allowed(self, cid):
        if self.options.variant_selection == "full-cloud":
            return self.view.candidate_by_id[cid].is_full_family
        return self.mode == "expanded" or cid in self.seed_ids[self.mode]

    def neighborhoods(self, search, state, start, round_index):
        """Prioritize retained-family variant edits, then numerical conflict witnesses."""
        if self.options.variant_selection == "full-cloud":
            return "no-eligible-work"
        proposed = set()
        context = semantic_id(
            "incumbent", [(a.candidate_id, a.pool) for a in state.items()]
        )
        limit = self.options.work_limits.repair_neighborhoods_per_round
        parent_families = defaultdict(list)
        for cid in sorted(state.assignments):
            parent_families[self.records[cid].family_id].append(cid)
        for family_id, parents in sorted(parent_families.items()):
            if len(proposed) >= limit:
                break
            search.tick()
            pool_ids = {state.assignments[cid] for cid in parents}
            self._expand(
                family_id,
                parents=tuple(parents),
                pool=next(iter(pool_ids)) if len(pool_ids) == 1 else None,
                context_digest=context,
            )
            proposed.add(family_id)
        probes = 0
        for cid in tuple(sorted(self.records)):
            if (
                probes >= self.options.work_limits.repair_candidate_probes_per_round
                or len(proposed) >= limit
            ):
                break
            if cid in state.assignments:
                continue
            search.tick()
            probes += 1
            self.work["contextual_probes"] += 1
            config = self.records[cid]
            if config.family_id in proposed:
                continue
            witnesses = []
            diagnostics = []
            if hasattr(self.oracle, "candidate_diagnostics") and not search.valid(cid):
                diagnostics.append(self.oracle.candidate_diagnostics(cid))
            if hasattr(self.oracle, "pair_diagnostics"):
                for other in sorted(state.assignments):
                    if (
                        probes
                        >= self.options.work_limits.repair_candidate_probes_per_round
                    ):
                        break
                    search.tick()
                    probes += 1
                    self.work["contextual_probes"] += 1
                    if search.conflict(cid, other):
                        diagnostics.append(self.oracle.pair_diagnostics(cid, other))
            for report in diagnostics:
                for key in ("dimer_edges", "rejected_products", "uncertainty_blocks"):
                    for witness in report.get(key, ()):
                        if key == "dimer_edges":
                            reasons = report.get("reasons", ())
                            score = witness["score"]
                            active = (
                                "dimer-threshold" in reasons
                                and self.policy.rejects(score)
                            )
                            exposure = (
                                "salvage-exposure-budget" in reasons
                                and score <= self.policy.strict_cutoff
                            )
                            if not (active or exposure):
                                continue
                        ids = tuple(witness.get("site_ids", ()))
                        if ids:
                            witnesses.append(
                                ProposalWitness(
                                    witness.get("id", semantic_id("witness", witness)),
                                    ids,
                                )
                            )
            if witnesses:
                self._expand(
                    config.family_id,
                    parents=(cid,),
                    witnesses=tuple(witnesses),
                    context_digest=context,
                )
                proposed.add(config.family_id)
        if (
            probes >= self.options.work_limits.repair_candidate_probes_per_round
            or len(proposed) >= limit
        ):
            return "work-cap"
        return "exhausted" if probes or proposed else "no-eligible-work"

    def ledger(self):
        return ConfigurationLedger(
            self.catalog.semantic_digest, tuple(self.records.values())
        )

    def metadata(self):
        return dict(self.work) | {
            "catalog_families": len(self.families),
            "searched_families": len(self.searched_families),
            "unsearched_families": len(self.families) - len(self.searched_families),
            "seed_families_examined": {
                mode: self.cursors[mode] for mode in ("full", "normal")
            },
            "omitted_materialized_ids": tuple(sorted(self.omitted)),
            "pending_neighborhood_streams_observed": self.pending_streams,
            "frontier_count_kind": "sum-of-pending-neighborhood-streams-at-call-return-not-exact-unseen-subsets",
        }


def _decisions(history, policy):
    def record(previous, current, phase):
        before = {a.candidate_id: a.pool for a in previous}
        after = {a.candidate_id: a.pool for a in current}
        context = semantic_id("incumbent", tuple(sorted(after.items())))
        for cid in sorted(set(before) | set(after)):
            if before.get(cid) == after.get(cid):
                continue
            prior = history.latest_event(cid)
            history.emit(
                stage_id=policy.stage_id,
                kind=(
                    "moved"
                    if cid in before and cid in after
                    else ("selected" if cid in after else "replaced")
                ),
                entity_ids=(cid,),
                parent_event_ids=(prior.id,) if prior else (),
                changes={
                    "phase": phase,
                    "old_pool": before.get(cid),
                    "new_pool": after.get(cid),
                    "incumbent_digest": context,
                },
            )

    return record


def _finish_history(result, *, complete):
    history = result.history
    selected = {a.configuration_id for a in result.assignments}
    stage = result.metadata["stage_policy"]["stage_id"]
    previous = history.last_complete_stage
    # A single checkpoint pass is intentional; never performed per family.
    touched = set()
    seed_context = set()
    failed_attempts = set()
    for event in history.iter_events(start=previous.events_count if previous else 0):
        touched.update(event.entity_ids)
        if event.kind == "seed-rejected":
            seed_context.update(event.entity_ids)
            failed_attempts.add(event.changes["attempt_id"])
    dispositions = {}
    validity = result.metadata["candidate_status"]
    for cid in touched:
        dispositions[cid] = (
            ("selected-strict" if stage == "strict" else "selected-salvage")
            if cid in selected
            else (
                "feasible-but-not-selected"
                if validity.get(cid) is True
                else (
                    "rejected-stage-policy"
                    if validity.get(cid) is False
                    else "search-limit-not-explored"
                )
            )
        )
    for entity in seed_context - set(validity) - selected:
        dispositions[entity] = "evaluated-seed-context"
    for attempt in failed_attempts:
        dispositions[attempt] = "rejected-seed-construction"
    snapshot = history.complete_stage(
        stage_id=stage,
        dispositions=dispositions,
        catalog_digest=result.catalog.semantic_digest,
        ledger_digest=result.ledger.semantic_digest,
        prior_snapshot_id=previous.id if previous else None,
        completeness="complete" if complete else "incomplete",
    )
    result.metadata["history_snapshot_id"] = snapshot.id


def search_allele_configurations(
    catalog: VariantCatalog,
    profile: AlleleConstraintProfile,
    options: AlleleSearchOptions,
    ledger: ConfigurationLedger,
    oracle,
    baseline=(),
    *,
    policy=StagePolicy(),
    history=None,
    clock=monotonic,
    cancelled=None,
    _complete_history=True,
) -> AlleleSearchResult:
    """Abstract combinatorial entry. A valid result here is not a scientific certificate."""
    options = options or AlleleSearchOptions()
    if ledger.catalog_digest != catalog.semantic_digest:
        raise ValueError("ledger/catalog mismatch")
    history = (
        history if history is not None else CoverageHistory(run_id="allele-search")
    )
    oracle = _AssessmentOracle(oracle, catalog, profile, policy, history)
    proposals = _Proposals(catalog, ledger, oracle, profile, policy, options, history)
    pools = _PoolEvaluations(oracle, catalog, profile, policy, history)
    factory = lambda assignments=(): _AlleleState(
        catalog, proposals.view, profile, options, pools, assignments
    )
    search = _Search(
        proposals.view,
        profile,
        options,
        oracle,
        clock,
        state_factory=factory,
        candidate_source=proposals.refresh,
        candidate_filter=proposals.allowed,
        aggregate_guard=pools.allows,
        record_hook=_decisions(history, policy),
        neighborhood_hook=proposals.neighborhoods,
        cancelled=cancelled,
        limits=asdict(options.work_limits)
        | {"repair_removed_candidates": options.exchange_width},
    )
    baseline = tuple(baseline)
    # Repeating the same exact configuration across pools cannot improve coverage
    # and weakly worsens exposure/burden; retain one explicitly recorded occurrence.
    normalized = {}
    for assignment in sorted(baseline, key=lambda a: (a.configuration_id, a.pool)):
        normalized.setdefault(assignment.configuration_id, assignment)
    duplicates = len(baseline) - len(normalized)
    internal = tuple(
        Assignment(a.configuration_id, a.pool) for a in normalized.values()
    )
    baseline_report = search.check(internal)
    if any(a.stage_id != policy.stage_id for a in baseline):
        baseline_report["valid"] = False
        baseline_report["violations"].append({"reason": "stage-mismatch"})
    if options.variant_selection == "full-cloud" and any(
        a.candidate_id in proposals.view.candidate_by_id
        and not proposals.view.candidate_by_id[a.candidate_id].is_full_family
        for a in internal
    ):
        baseline_report["valid"] = False
        baseline_report["violations"].append(
            {"reason": "baseline-incompatible-with-full-cloud-policy"}
        )
    if baseline_report["valid"]:
        search.record(search.state(internal), "supplied-baseline")
    retained_baseline = search.best
    empty_summary = allele_summary(catalog, (), ledger, goal=options.coverage_target)
    begin = clock()
    search.deadline = begin + options.time_limit
    stop = "completed"
    starts = repairs = 0
    seeds_completed = []
    seed_modes = tuple(
        mode
        for mode in options.seed_modes
        if options.variant_selection != "full-cloud" or mode == "full"
    )
    scheduler = PhaseScheduler(search, proposals, history, policy.stage_id)
    reserved = options.phase_scheduling == "reserved"
    initial_weight = (
        (0.2 if seed_modes else 0) + 0.4 + (0.4 if options.repair_rounds else 0)
    )
    initial = Reservation(search.deadline, initial_weight) if reserved else None
    planned = ["seed:" + mode for mode in seed_modes]
    for start in range(options.starts):
        planned.append(f"construction:{start}")
        for round_index in range(options.repair_rounds):
            planned.extend(
                f"repair:{start}:{round_index}/{part}"
                for part in ("preparation", "cleanup", "exchange")
            )
    exhausted_seeds = []
    try:
        if empty_summary.status == "no-assessable-targets":
            stop = "no-assessable-targets"
        else:
            for mode in seed_modes:
                proposals.mode = mode
                state = search.state()
                outcome = scheduler.run(
                    "seed:" + mode,
                    partial(search.fill, state, 0, "seed:" + mode),
                    initial,
                    0.2 / len(seed_modes),
                )
                if outcome != "phase-time-limit":
                    seeds_completed.append(mode)
                if outcome == "exhausted":
                    exhausted_seeds.append(mode)
            proposals.mode = "expanded"
            for start in range(options.starts):
                state = search.state()
                outcome = scheduler.run(
                    f"construction:{start}",
                    partial(search.fill, state, start, f"construction:{start}"),
                    initial if start == 0 else None,
                    0.4,
                )
                if outcome != "phase-time-limit":
                    starts += 1
                for round_index in range(options.repair_rounds):
                    # Reserve the first cycle only. Subsequent starts retain the
                    # same deterministic serial order and share the remainder.
                    if initial is not None and start == 0:
                        repair_end = initial.cutoff(
                            clock(), 0.4 / options.repair_rounds
                        )
                        subphase = Reservation(repair_end, 1.0)
                    else:
                        subphase = None
                    outcomes = []
                    for part, weight, action in (
                        (
                            "preparation",
                            0.2,
                            partial(search.prepare_repair, start, round_index),
                        ),
                        (
                            "cleanup",
                            0.2,
                            partial(search.cleanup, f"cleanup:{start}:{round_index}"),
                        ),
                        ("exchange", 0.6, partial(search.exchange, start, round_index)),
                    ):
                        outcomes.append(
                            scheduler.run(
                                f"repair:{start}:{round_index}/{part}",
                                action,
                                subphase,
                                weight,
                            )
                        )
                    if "phase-time-limit" not in outcomes:
                        repairs += 1
    except _TimeLimit:
        stop = "time-limit"
    except _Cancelled:
        stop = "cancelled"
    elapsed = max(0, clock() - begin)
    search.deadline = math.inf
    search.phase_deadline = math.inf
    validation = search.check(search.best)
    if not validation["valid"]:
        search.best = retained_baseline
        validation = search.check(search.best)
        if not validation["valid"]:
            raise ValueError("abstract validation rejected baseline")
    final_ledger = proposals.ledger()
    assignments = tuple(
        AlleleAssignment(a.candidate_id, a.pool, policy.stage_id) for a in search.best
    )
    summary = allele_summary(
        catalog, assignments, final_ledger, goal=options.coverage_target
    )
    metadata = {
        "schemaVersion": "primalscheme3.panel-optimizer/v2",
        "algorithm": "bounded-allele-coverage/v1",
        "objective_version": "soft-observed-allele-coverage/v1",
        "objective_direction": "minimize-lexicographically",
        "objective_terms": [
            "negative-mean-utility",
            "negative-mean-coverage",
            "negative-worst-coverage",
            "violating-edges",
            "incident-species",
            "sequence-pool-instances",
            "assignments",
            "pool-burden-range",
            "size-deviation",
            "canonical-assignments",
        ],
        "options": asdict(options),
        "stage_policy": asdict(policy),
        "work": dict(search.work),
        "proposal_work": proposals.metadata(),
        "completed_starts": starts,
        "completed_repair_rounds": repairs,
        "completed_seed_modes": seeds_completed,
        "exhausted_seed_modes": exhausted_seeds,
        "scheduling_policy": scheduling_policy(options.phase_scheduling),
        "phase_progress": scheduler.progress,
        "phases_not_entered": [
            name
            for name in planned
            if name not in {p["phase"] for p in scheduler.progress}
        ],
        "active_phase_at_stop": scheduler.active_at_stop,
        "deadline_overshoot_seconds": max(0.0, elapsed - options.time_limit),
        "fixed_work_completed": stop == "completed"
        and all(p["outcome"] != "phase-time-limit" for p in scheduler.progress),
        "family_streams_exhausted": {
            mode: cursor >= len(proposals.families)
            for mode, cursor in proposals.cursors.items()
        },
        "effective_seed_modes": list(seed_modes),
        "work_truncated": bool(
            search.work["construction_limit_hits"]
            or search.work["repair_limit_hits"]
            or proposals.work["truncated_neighborhoods"]
            or len(proposals.searched_families) < len(proposals.families)
        ),
        "repairs_accepted": search.repairs_accepted,
        "stop_reason": stop,
        "search_seconds": elapsed,
        "baseline": {
            "supplied_count": len(baseline),
            "supplied_valid": baseline_report["valid"],
            "violations": baseline_report["violations"],
            "redundant_configuration_placements_removed": duplicates,
        },
        "candidate_status": dict(search.valid_cache),
        "objective_history": search.history,
        "limitations": [
            "Bounded heuristic; no global optimality guarantee.",
            "Unsearched families and unexpanded subset streams may contain better configurations.",
            "Fixed-work reproducibility excludes cancellation and wall-time cutoffs.",
            "One proposal batch or oracle predicate is indivisible; wall-time can overrun by that work.",
            "Abstract oracle checks do not establish scientific panel validity.",
        ],
    }
    result = AlleleSearchResult(
        catalog,
        final_ledger,
        assignments,
        summary,
        validation,
        search.state(search.best).key(),
        metadata,
        history,
    )
    if _complete_history:
        _finish_history(result, complete=False)
    return result


def search_allele_assignments(
    catalog: VariantCatalog,
    profile: AlleleConstraintProfile,
    options: AlleleSearchOptions | None = None,
    ledger: ConfigurationLedger | None = None,
    baseline=(),
    *,
    policy=StagePolicy(),
    history=None,
    clock=monotonic,
    cancelled=None,
    score_cache=None,
) -> AlleleSearchResult:
    """Production search with independently reconstructed baseline and final validation."""
    options = options or AlleleSearchOptions()
    ledger = (
        ledger
        if ledger is not None
        else ConfigurationLedger(catalog.semantic_digest, ())
    )
    baseline = tuple(baseline)
    baseline_report = validate_allele_assignments(
        catalog, baseline, ledger, profile, policy=policy, goal=options.coverage_target
    )
    valid_baseline = baseline if baseline_report["valid"] else ()
    oracle = AlleleCompatibilityOracle(
        catalog, ledger, profile, policy, score_cache=score_cache
    )
    result = search_allele_configurations(
        catalog,
        profile,
        options,
        ledger,
        oracle,
        valid_baseline,
        policy=policy,
        history=history,
        clock=clock,
        cancelled=cancelled,
        _complete_history=False,
    )
    validation = validate_allele_assignments(
        catalog,
        result.assignments,
        result.ledger,
        profile,
        policy=policy,
        goal=options.coverage_target,
        expected_summary=result.coverage,
        history=result.history,
    )
    if not validation["valid"]:
        result.metadata["rejected_final_validation"] = validation
        eligible = (
            valid_baseline if result.metadata["baseline"]["supplied_valid"] else ()
        )
        unique = {}
        for assignment in sorted(eligible, key=lambda a: (a.configuration_id, a.pool)):
            unique.setdefault(assignment.configuration_id, assignment)
        fallback = tuple(
            AlleleAssignment(a.candidate_id, a.pool, policy.stage_id)
            for a in _canonical(
                Assignment(a.configuration_id, a.pool) for a in unique.values()
            )
        )
        summary = allele_summary(
            catalog, fallback, result.ledger, goal=options.coverage_target
        )
        validation = validate_allele_assignments(
            catalog,
            fallback,
            result.ledger,
            profile,
            policy=policy,
            goal=options.coverage_target,
            expected_summary=summary,
            history=result.history,
        )
        if not validation["valid"]:
            raise ValueError("independent validation rejected fallback baseline")
        _decisions(result.history, policy)(
            tuple(Assignment(a.configuration_id, a.pool) for a in result.assignments),
            tuple(Assignment(a.configuration_id, a.pool) for a in fallback),
            "independent-validation-fallback",
        )
        view = _CatalogView(catalog)
        view.candidate_by_id.update(
            {c.id: _ConfigurationView(c, catalog) for c in result.ledger.configurations}
        )
        pools = _PoolEvaluations(oracle, catalog, profile, policy, result.history)
        objective = _AlleleState(
            catalog,
            view,
            profile,
            options,
            pools,
            tuple(Assignment(a.configuration_id, a.pool) for a in fallback),
        ).key()
        result = replace(
            result, assignments=fallback, coverage=summary, objective=objective
        )
        result.metadata["validation_fallback"] = True
    result = replace(result, validation=validation)
    result.metadata["baseline"]["scientific_validation"] = baseline_report
    result.metadata["baseline"]["supplied_count"] = len(baseline)
    result.metadata["limitations"] = [
        x for x in result.metadata["limitations"] if not x.startswith("Abstract oracle")
    ]
    _finish_history(result, complete=True)
    return result
