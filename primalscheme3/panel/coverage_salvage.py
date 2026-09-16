"""Optional research salvage tiers with a retained independently valid strict result.

Exposure caps are computational search limits, not calibrated assay safety limits.
All tier interactions use the same chemistry, support and specificity contracts.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace
from time import monotonic

from .allele_search import AlleleSearchOptions, search_allele_assignments
from .allele_validation import StagePolicy, validate_allele_assignments


@dataclass(frozen=True)
class SalvageOptions:
    mode: str = "off"
    thresholds: tuple[float, ...] = (-28.0, -30.0, -32.0)
    max_stages: int = 3
    max_edges_per_pool: int = 8
    max_oligos_per_pool: int = 4
    time_limit: float = 60.0

    def __post_init__(self):
        if self.mode not in ("off", "bounded"):
            raise ValueError("salvage mode must be off or bounded")
        for name, minimum in (
            ("max_stages", 1),
            ("max_edges_per_pool", 0),
            ("max_oligos_per_pool", 0),
        ):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(name + " must be an integer >= " + str(minimum))
        if self.max_stages > 3:
            raise ValueError("at most three salvage stages supported")
        if (
            isinstance(self.time_limit, bool)
            or not math.isfinite(self.time_limit)
            or self.time_limit <= 0
        ):
            raise ValueError("salvage time_limit must be positive and finite")
        previous = -26.0
        for cutoff in self.thresholds:
            if (
                isinstance(cutoff, bool)
                or not math.isfinite(cutoff)
                or cutoff >= previous
            ):
                raise ValueError(
                    "thresholds must be finite and strictly decreasing from -26"
                )
            previous = cutoff
        if not self.thresholds or len(self.thresholds) > self.max_stages:
            raise ValueError("threshold ladder must contain one to max_stages entries")
        object.__setattr__(self, "thresholds", tuple(float(x) for x in self.thresholds))


@dataclass(frozen=True)
class SalvageTier:
    stage_id: str
    policy: StagePolicy
    status: str
    result: object | None
    comparison: dict
    error: str | None = None
    elapsed_seconds: float = 0.0


@dataclass(frozen=True)
class SalvageRun:
    strict: object
    tiers: tuple[SalvageTier, ...]
    options: SalvageOptions
    stop_reason: str = "completed"


def compare_to_strict(strict, tier):
    def sites(result):
        return {
            (sid, a.pool)
            for a in result.assignments
            for sid in (
                result.ledger.configuration_by_id[a.configuration_id].forward_site_ids
                + result.ledger.configuration_by_id[a.configuration_id].reverse_site_ids
            )
        }

    old, new = sites(strict), sites(tier)

    def delta(a, b):
        return None if a is None or b is None else b - a

    old_targets = {t.target_id: t for t in strict.coverage.targets}
    old_classes = {c.allele_id: c for c in strict.coverage.classes}
    return {
        "relative_to": "strict",
        "mean_coverage_delta": delta(
            strict.coverage.mean_coverage, tier.coverage.mean_coverage
        ),
        "mean_utility_delta": delta(
            strict.coverage.mean_utility, tier.coverage.mean_utility
        ),
        "added_site_pool_instances": [list(x) for x in sorted(new - old)],
        "removed_site_pool_instances": [list(x) for x in sorted(old - new)],
        "moved_sites": sorted({s for s, p in old - new} & {s for s, p in new - old}),
        "target_changes": [
            {
                "target_id": t.target_id,
                "strict_fraction": old_targets[t.target_id].fraction,
                "tier_fraction": t.fraction,
                "delta": delta(old_targets[t.target_id].fraction, t.fraction),
            }
            for t in tier.coverage.targets
        ],
        "class_changes": [
            {
                "allele_id": c.allele_id,
                "target_id": c.target_id,
                "strict_covered": old_classes[c.allele_id].covered_count,
                "tier_covered": c.covered_count,
                "delta": c.covered_count - old_classes[c.allele_id].covered_count,
            }
            for c in tier.coverage.classes
        ],
        "exposure": {
            str(pool): {
                "strict_violating_edges": v.get("violating_edge_count", 0),
                "incident_species": v.get("incident_species_count", 0),
            }
            for pool, v in tier.validation.get("pools", {}).items()
        },
        "interpretation": "Coverage predictions and cumulative exposure relative to -26; not experimental performance or safety.",
    }


def run_salvage(
    strict,
    profile,
    options=None,
    *,
    search_options=None,
    history=None,
    on_tier=None,
    cancelled=None,
    score_cache=None,
    observer=None,
):
    options = options or SalvageOptions()
    if strict.validation.get("stage_policy") != asdict(
        StagePolicy()
    ) or strict.metadata.get("stage_policy") != asdict(StagePolicy()):
        raise ValueError("salvage requires exact strict-stage identity")
    if not strict.validation.get("valid") or any(
        a.stage_id != "strict" for a in strict.assignments
    ):
        raise ValueError("salvage requires a validated strict result")
    if profile.to_dict() != strict.validation.get("profile"):
        raise ValueError("salvage cannot change the strict scientific profile")
    if (
        search_options is not None
        and search_options.coverage_target != strict.coverage.goal
    ):
        raise ValueError("salvage cannot change the strict coverage objective")
    if options.mode == "off":
        return SalvageRun(strict, (), options, "off")
    # A supplied success flag cannot substitute for fresh baseline validation.
    baseline_check = validate_allele_assignments(
        strict.catalog,
        strict.assignments,
        strict.ledger,
        profile,
        goal=strict.coverage.goal,
        expected_summary=strict.coverage,
    )
    if not baseline_check["valid"]:
        raise ValueError("strict baseline failed fresh validation")
    if search_options is None:
        raw = dict(strict.metadata["options"])
        from .allele_search import AlleleWorkLimits

        raw["work_limits"] = AlleleWorkLimits(**raw["work_limits"])
        search_options = AlleleSearchOptions(**raw)
    search_options = replace(search_options, time_limit=options.time_limit)
    history = history if history is not None else strict.history
    score_cache = {} if score_cache is None else score_cache
    previous = strict
    tiers = []
    stop_reason = "completed"
    for index, cutoff in enumerate(options.thresholds, 1):
        if cancelled is not None and cancelled():
            stop_reason = "cancelled"
            history.emit(
                stage_id="salvage",
                kind="salvage-run-stopped",
                entity_ids=(),
                changes={
                    "reason": "cancelled",
                    "next_stage": f"salvage-{index}",
                    "strict_preserved": True,
                },
            )
            if hasattr(history, "checkpoint"):
                history.checkpoint()
            break
        stage_id = f"salvage-{index}"
        policy = StagePolicy(
            stage_id,
            cutoff,
            -26,
            options.max_edges_per_pool,
            options.max_oligos_per_pool,
        )
        begin = monotonic()
        result = None
        comparison = {}
        try:
            history.emit(
                stage_id=stage_id,
                kind="salvage-stage-start",
                entity_ids=(),
                changes={
                    "policy": asdict(policy),
                    "baseline_stage": previous.metadata["stage_policy"]["stage_id"],
                    "limits_interpretation": "computational exposure; no assay safety guarantee",
                },
            )
            baseline = tuple(
                replace(a, stage_id=stage_id) for a in previous.assignments
            )
            result = search_allele_assignments(
                strict.catalog,
                profile,
                search_options,
                previous.ledger,
                baseline,
                policy=policy,
                history=history,
                cancelled=cancelled,
                score_cache=score_cache,
                observer=observer,
            )
            if not result.validation.get("valid"):
                raise ValueError("tier failed independent validation")
            comparison = compare_to_strict(strict, result)
            if on_tier is not None:
                on_tier(result)
            tiers.append(
                SalvageTier(
                    stage_id,
                    policy,
                    "validated",
                    result,
                    comparison,
                    elapsed_seconds=monotonic() - begin,
                )
            )
            previous = result
        except Exception as error:
            failed_ids = (
                tuple(a.configuration_id for a in result.assignments)
                if result is not None
                else ()
            )
            history.emit(
                stage_id=stage_id,
                kind="salvage-stage-failed",
                entity_ids=failed_ids,
                changes={"error": str(error), "strict_preserved": True},
            )
            prior = history.last_complete_stage
            touched = {
                entity
                for event in history.iter_events(
                    start=prior.events_count if prior else 0
                )
                for entity in event.entity_ids
            }
            history.complete_stage(
                stage_id=stage_id + "-failed",
                dispositions={e: "tier-publication-failed" for e in touched},
                catalog_digest=strict.catalog.semantic_digest,
                ledger_digest=(
                    result.ledger if result is not None else previous.ledger
                ).semantic_digest,
                prior_snapshot_id=prior.id if prior else None,
            )
            tiers.append(
                SalvageTier(
                    stage_id,
                    policy,
                    "failed",
                    result,
                    comparison,
                    str(error),
                    monotonic() - begin,
                )
            )
        if result is not None and result.metadata.get("stop_reason") == "cancelled":
            stop_reason = "cancelled"
            history.emit(
                stage_id="salvage",
                kind="salvage-run-stopped",
                entity_ids=(),
                changes={
                    "reason": "cancelled",
                    "last_stage": stage_id,
                    "strict_preserved": True,
                },
            )
            if hasattr(history, "checkpoint"):
                history.checkpoint()
            break
    return SalvageRun(strict, tuple(tiers), options, stop_reason)


def select_primary(run, name="strict"):
    if name == "strict":
        return run.strict
    for tier in run.tiers:
        if tier.stage_id == name:
            if tier.status != "validated" or tier.result is None:
                raise ValueError(
                    "requested primary tier failed validation or publication"
                )
            return tier.result
    raise ValueError("requested primary tier is unavailable")
