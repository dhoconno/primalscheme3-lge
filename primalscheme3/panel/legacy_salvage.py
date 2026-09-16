"""Bounded dimer salvage over the retained strict legacy candidate pool.

The legacy generator and selector stay authoritative.  This module only admits
additional already-generated pairs after the strict panel is complete, and it
keeps a compact audit of each bounded pass.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from time import monotonic
from typing import Any

from primalschemers import do_pool_interact

from primalscheme3.core.mismatches import detect_new_products
from primalscheme3.core.multiplex import PrimerPairCheck
from primalscheme3.panel.allele_validation import numerical_dimer_score


@dataclass(frozen=True)
class LegacySalvageOptions:
    mode: str = "off"
    thresholds: tuple[float, ...] = (-28.0, -30.0, -32.0)
    floor: float = -32.0
    max_edges_per_pool: int = 8
    max_incident_species_per_pool: int = 4
    min_reference_gain: int = 1
    max_candidate_evaluations: int = 10000

    def __post_init__(self):
        if self.mode not in {"off", "bounded"}:
            raise ValueError("legacy salvage mode must be off or bounded")
        if isinstance(self.floor, bool) or not math.isfinite(self.floor):
            raise ValueError("legacy salvage floor must be finite")
        if self.mode == "off":
            return
        if not self.thresholds or len(self.thresholds) > 16:
            raise ValueError("legacy salvage thresholds must contain 1 through 16 values")
        # Construction validates the self-contained defaults.  A caller may
        # use a different configured strict cutoff; that relationship is
        # checked by ``validate_for_strict_cutoff`` at run time.
        previous = float("inf")
        for cutoff in self.thresholds:
            if (
                isinstance(cutoff, bool)
                or not math.isfinite(cutoff)
                or cutoff >= previous
                or cutoff < self.floor
            ):
                raise ValueError(
                    "legacy salvage thresholds must be finite, decreasing, and stay above the floor"
                )
            previous = cutoff
        for name in (
            "max_edges_per_pool",
            "max_incident_species_per_pool",
            "min_reference_gain",
            "max_candidate_evaluations",
        ):
            value = getattr(self, name)
            if type(value) is not int or value < (1 if name in {"min_reference_gain", "max_candidate_evaluations"} else 0):
                raise ValueError(name + " must be an integer in the allowed range")
        object.__setattr__(self, "thresholds", tuple(float(x) for x in self.thresholds))
        object.__setattr__(self, "floor", float(self.floor))

    def validate_for_strict_cutoff(self, strict_cutoff: float) -> None:
        if isinstance(strict_cutoff, bool) or not math.isfinite(strict_cutoff):
            raise ValueError("legacy strict dimer cutoff must be finite")
        if self.mode == "off":
            return
        if self.floor > strict_cutoff:
            raise ValueError("legacy salvage floor must be <= the configured strict cutoff")
        previous = strict_cutoff
        for cutoff in self.thresholds:
            if cutoff >= previous or cutoff < self.floor:
                raise ValueError(
                    "legacy salvage thresholds must decrease below the configured strict cutoff and stay above the floor"
                )
            previous = cutoff

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["thresholds"] = list(self.thresholds)
        return result


@dataclass(frozen=True)
class LegacyCandidateEvaluation:
    candidate_id: str
    pool: int | None
    status: str
    reason: str | None
    marginal_trimmed_gain: int
    marginal_full_gain: int
    conflicts: tuple[dict[str, Any], ...] = ()
    cumulative_edges: int = 0
    cumulative_incident_species: int = 0


@dataclass(frozen=True)
class LegacySalvageStage:
    stage_id: str
    cutoff: float
    accepted: tuple[str, ...]
    evaluations: tuple[LegacyCandidateEvaluation, ...]
    rejections: dict[str, dict[str, Any]]
    coverage: dict[str, dict[str, Any]]
    elapsed_seconds: float


@dataclass(frozen=True)
class LegacySalvageRun:
    options: LegacySalvageOptions
    strict_candidate_ids: tuple[str, ...]
    strict_coverage: dict[str, dict[str, Any]]
    candidate_manifest: tuple[dict[str, Any], ...]
    stages: tuple[LegacySalvageStage, ...]
    stop_reason: str

    def to_dict(self) -> dict[str, Any]:
        manifest = []
        accepted_by_id = {
            candidate_id: (stage.stage_id, evaluation.pool)
            for stage in self.stages
            for candidate_id in stage.accepted
            for evaluation in stage.evaluations
            if evaluation.candidate_id == candidate_id
        }
        rejected_by_id = {
            candidate_id: stage.stage_id
            for stage in self.stages
            for candidate_id in stage.rejections
            if candidate_id != "__meta__"
        }
        for item in self.candidate_manifest:
            entry = dict(item)
            candidate_id = entry["candidateId"]
            if candidate_id in accepted_by_id:
                stage_id, pool = accepted_by_id[candidate_id]
                entry.update({"status": "accepted", "acceptedStage": stage_id, "pool": pool})
            elif candidate_id in rejected_by_id and entry.get("status") == "available":
                entry.update({"status": "rejected", "lastRejectedStage": rejected_by_id[candidate_id]})
            manifest.append(entry)
        return {
            "schemaVersion": "primalscheme3.legacy-dimer-salvage/v1",
            "options": self.options.to_dict(),
            "strictCandidateIds": list(self.strict_candidate_ids),
            "strictCoverage": self.strict_coverage,
            "candidateManifest": manifest,
            "stages": [
                {
                    "stageId": stage.stage_id,
                    "cutoff": stage.cutoff,
                    "accepted": list(stage.accepted),
                    "evaluations": [asdict(item) for item in stage.evaluations],
                    "rejections": stage.rejections,
                    "coverage": stage.coverage,
                    "elapsedSeconds": stage.elapsed_seconds,
                }
                for stage in self.stages
            ],
            "stopReason": self.stop_reason,
        }


def _canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def legacy_candidate_id(primer_pair, msa) -> str:
    """Stable ID for one post-generation legacy candidate."""
    payload = {
        "msaIndex": int(msa.msa_index),
        "targetOccurrence": str(getattr(msa, "_chrom_name", msa.msa_index)),
        "fullInterval": [min(primer_pair.fprimer.starts()), max(primer_pair.rprimer.ends())],
        "interiorInterval": [int(primer_pair.fprimer.end), int(primer_pair.rprimer.start)],
        "forwardOligos": sorted(str(item).upper() for item in primer_pair.fprimer.seqs()),
        "reverseOligos": sorted(str(item).upper() for item in primer_pair.rprimer.seqs()),
    }
    digest = hashlib.sha256(_canonical_json(payload)).hexdigest()
    return "legacy-candidate-" + digest


def native_dimer_score(first: str, second: str):
    """Return the exact numeric score used by the installed native kernel."""
    return numerical_dimer_score(first, second)


def _reference_length(msa) -> int:
    mapping = getattr(msa, "_mapping_array", None)
    if mapping is not None:
        observed = [int(value) for value in mapping if value is not None]
        return max(observed, default=-1) + 1
    array = getattr(msa, "array", ())
    return len(array[0]) if array else 0


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def _coverage(panel, msa_dict):
    result = {}
    for index, msa in sorted(msa_dict.items()):
        length = _reference_length(msa)
        trim, full = [], []
        for primer_pair in panel.all_primerpairs():
            if int(primer_pair.msa_index) != int(index):
                continue
            intervals = _valid_intervals(primer_pair, msa, reference_length=length)
            if intervals is None:
                raise ValueError(
                    f"primer pair interval is outside ungapped reference for MSA {index}"
                )
            trim.append(intervals[0])
            full.append(intervals[1])
        trim_union, full_union = _merge(trim), _merge(full)
        result[str(index)] = {
            "referenceLength": length,
            "trimmedIntervals": [list(item) for item in trim_union],
            "fullIntervals": [list(item) for item in full_union],
            "trimmedBases": sum(end - start for start, end in trim_union),
            "fullBases": sum(end - start for start, end in full_union),
            "trimmedFraction": (sum(end - start for start, end in trim_union) / length if length else 0.0),
            "fullFraction": (sum(end - start for start, end in full_union) / length if length else 0.0),
        }
    return result


def _valid_intervals(primer_pair, msa, *, reference_length: int | None = None):
    length = _reference_length(msa) if reference_length is None else reference_length
    trim = (int(primer_pair.fprimer.end), int(primer_pair.rprimer.start))
    full = (min(primer_pair.fprimer.starts()), max(primer_pair.rprimer.ends()))
    if not (0 <= trim[0] < trim[1] <= length):
        return None
    if not (0 <= full[0] < full[1] <= length):
        return None
    return trim, full


def _gain(panel, primer_pair, msa_dict):
    before = _coverage(panel, msa_dict)
    index = str(primer_pair.msa_index)
    msa = msa_dict[primer_pair.msa_index]
    intervals = _valid_intervals(primer_pair, msa)
    if intervals is None:
        return None
    trim, full = intervals
    trim_old = before[index]["trimmedBases"]
    full_old = before[index]["fullBases"]
    trim_new = sum(end - start for start, end in _merge((*[tuple(x) for x in before[index]["trimmedIntervals"]], trim)))
    full_new = sum(end - start for start, end in _merge((*[tuple(x) for x in before[index]["fullIntervals"]], full)))
    return trim_new - trim_old, full_new - full_old


def _gain_from_snapshot(primer_pair, msa_dict, before):
    msa = msa_dict[primer_pair.msa_index]
    intervals = _valid_intervals(
        primer_pair,
        msa,
        reference_length=before[str(primer_pair.msa_index)]["referenceLength"],
    )
    if intervals is None:
        return None
    trim, full = intervals
    index = str(primer_pair.msa_index)
    trim_union = _merge((*[tuple(x) for x in before[index]["trimmedIntervals"]], trim))
    full_union = _merge((*[tuple(x) for x in before[index]["fullIntervals"]], full))
    trim_new = sum(end - start for start, end in trim_union)
    full_new = sum(end - start for start, end in full_union)
    return trim_new - before[index]["trimmedBases"], full_new - before[index]["fullBases"]


def _pool_risk(panel, primer_pair, pool, cutoff, strict_cutoff, edge_keys, species):
    existing = sorted({str(seq).upper() for item in panel._pools[pool] for seq in item.all_seqs()})
    conflicts = []
    new_edges = set()
    new_species = set()
    for sequence in sorted({str(seq).upper() for seq in primer_pair.all_seqs()}):
        for other in existing:
            first, second = sequence.encode(), other.encode()
            if do_pool_interact([first], [second], cutoff):
                continue
            if do_pool_interact([first], [second], strict_cutoff):
                score = native_dimer_score(sequence, other)
                value = float(score.score)
                key = tuple(sorted(score.sequences))
                new_edges.add(key)
                new_species.update(key)
                conflicts.append({"sequences": list(score.sequences), "score": value, "reason": "relaxed-edge"})
    return conflicts, new_edges, new_species


def run_legacy_salvage(panel, msa_dict, options: LegacySalvageOptions | None = None, *, strict_cutoff: float | None = None, logger=None) -> LegacySalvageRun:
    options = options or LegacySalvageOptions()
    configured_strict_cutoff = -26.0 if strict_cutoff is None else float(strict_cutoff)
    options.validate_for_strict_cutoff(configured_strict_cutoff)
    reference_lengths = {index: _reference_length(msa) for index, msa in msa_dict.items()}
    strict_ids = tuple(sorted(legacy_candidate_id(item, msa_dict[item.msa_index]) for item in panel.all_primerpairs()))
    strict_coverage = _coverage(panel, msa_dict)
    if options.mode == "off":
        return LegacySalvageRun(options, strict_ids, strict_coverage, (), (), "off")

    candidates = []
    for _index, msa in sorted(msa_dict.items()):
        for item in getattr(msa, "primerpairs", ()):
            candidates.append((legacy_candidate_id(item, msa), item, msa))
    candidates.sort(key=lambda item: item[0])
    strict_objects = {id(item) for item in panel.all_primerpairs()}
    strict_id_set = set(strict_ids)
    retained_ids = set()
    retained_candidates = []
    for item in candidates:
        candidate_id = item[0]
        if id(item[1]) in strict_objects or candidate_id in strict_id_set:
            continue
        # Content-identical pairs can occur more than once in a generated
        # pool.  A stable candidate ID represents one salvage opportunity.
        if candidate_id in retained_ids:
            continue
        retained_ids.add(candidate_id)
        retained_candidates.append(item)
    candidates = retained_candidates

    def manifest_item(candidate_id, primer_pair, msa, *, status, pool=None):
        intervals = _valid_intervals(
            primer_pair,
            msa,
            reference_length=reference_lengths[msa.msa_index],
        )
        amplicon_number = getattr(primer_pair, "amplicon_number", None)
        amplicon_prefix = getattr(primer_pair, "amplicon_prefix", None)
        published_name = (
            f"{amplicon_prefix}_{amplicon_number}"
            if amplicon_prefix is not None and amplicon_number is not None
            else None
        )
        return {
            "candidateId": candidate_id,
            "msaIndex": int(msa.msa_index),
            "targetOccurrence": str(getattr(msa, "_chrom_name", msa.msa_index)),
            "fullInterval": list(intervals[1]) if intervals else None,
            "trimmedInterval": list(intervals[0]) if intervals else None,
            "forwardOligos": sorted(str(item).upper() for item in primer_pair.fprimer.seqs()),
            "reverseOligos": sorted(str(item).upper() for item in primer_pair.rprimer.seqs()),
            "status": status,
            "pool": pool,
            "ampliconNumber": amplicon_number,
            "ampliconPrefix": amplicon_prefix,
            "publishedAmpliconName": published_name,
            "scope": "post-generation-strict-legacy-pool",
        }
    strict_manifest = []
    for pool, pool_pairs in enumerate(panel._pools):
        for primer_pair in pool_pairs:
            strict_manifest.append(
                manifest_item(
                    legacy_candidate_id(primer_pair, msa_dict[primer_pair.msa_index]),
                    primer_pair,
                    msa_dict[primer_pair.msa_index],
                    status="strict-incumbent",
                    pool=pool,
                )
            )
    candidate_manifest = [manifest_item(*item, status="available") for item in candidates]
    manifest_entries = sorted((*strict_manifest, *candidate_manifest), key=lambda item: item["candidateId"])

    edge_keys = [set() for _ in range(panel.n_pools)]
    species = [set() for _ in range(panel.n_pools)]
    stages = []
    accepted_ids = set()
    stop_reason = "completed"
    work_limited = False
    for stage_index, cutoff in enumerate(options.thresholds, 1):
        started = monotonic()
        accepted = []
        evaluations = []
        rejections = {}
        remaining = [item for item in candidates if item[0] not in accepted_ids]
        evaluated = 0
        stage_work_limited = False
        rejected_this_stage = set()
        while remaining:
            if evaluated >= options.max_candidate_evaluations:
                work_limited = True
                stage_work_limited = True
                break
            before = _coverage(panel, msa_dict)
            gains = []
            for candidate_id, primer_pair, msa in remaining:
                if candidate_id in rejected_this_stage or candidate_id in accepted_ids:
                    continue
                gain = _gain_from_snapshot(primer_pair, msa_dict, before)
                if gain is None:
                    gains.append((0, 0, candidate_id, candidate_id, primer_pair, msa, None, None))
                    continue
                trim_gain, full_gain = gain
                gains.append((-(trim_gain / max(_reference_length(msa), 1)), -trim_gain, candidate_id, candidate_id, primer_pair, msa, trim_gain, full_gain))
            if not gains:
                break
            gains.sort()
            best_item = None
            for _, _, candidate_id, _, primer_pair, _msa, trim_gain, full_gain in gains:
                if evaluated >= options.max_candidate_evaluations:
                    work_limited = True
                    stage_work_limited = True
                    break
                rejected_this_stage.add(candidate_id)
                evaluated += 1
                if trim_gain is None:
                    rejections[candidate_id] = {"reason": "geometry"}
                    continue
                if trim_gain < options.min_reference_gain:
                    rejections[candidate_id] = {"reason": "reference-gain", "marginalTrimmedGain": trim_gain}
                    continue
                pool_reasons = {}
                admissible = []
                for pool in range(panel.n_pools):
                    check = panel.check_primerpair_can_be_added(
                        primer_pair, pool, dimer_score=cutoff
                    )
                    if check is not PrimerPairCheck.OK:
                        pool_reasons[str(pool)] = check.name.lower()
                        continue
                    conflicts, new_edges, new_species = _pool_risk(
                        panel, primer_pair, pool, cutoff, configured_strict_cutoff, edge_keys[pool], species[pool]
                    )
                    if len(edge_keys[pool] | new_edges) > options.max_edges_per_pool:
                        pool_reasons[str(pool)] = "edge-budget"
                        continue
                    if len(species[pool] | new_species) > options.max_incident_species_per_pool:
                        pool_reasons[str(pool)] = "species-budget"
                        continue
                    admissible.append(
                        (
                            len(new_edges - edge_keys[pool]),
                            len(new_species - species[pool]),
                            pool,
                            conflicts,
                            new_edges,
                            new_species,
                        )
                    )
                if admissible:
                    _, _, pool, conflicts, new_edges, new_species = min(admissible)
                    best_item = (
                        candidate_id,
                        primer_pair,
                        pool,
                        conflicts,
                        new_edges,
                        new_species,
                        trim_gain,
                        full_gain,
                    )
                if best_item is None:
                    rejections[candidate_id] = {
                        "reason": "pool-rejected",
                        "marginalTrimmedGain": trim_gain,
                        "pools": pool_reasons,
                    }
                else:
                    break
            if best_item is None:
                continue
            candidate_id, primer_pair, pool, conflicts, new_edges, new_species, trim_gain, full_gain = best_item
            panel._add_primerpair(primer_pair, pool, int(primer_pair.msa_index))
            edge_keys[pool].update(new_edges)
            species[pool].update(new_species)
            accepted.append(candidate_id)
            accepted_ids.add(candidate_id)
            rejected_this_stage.add(candidate_id)
            for entry in manifest_entries:
                if entry["candidateId"] == candidate_id:
                    entry.update(
                        {
                            "status": "accepted",
                            "pool": pool,
                            "ampliconNumber": getattr(primer_pair, "amplicon_number", None),
                            "ampliconPrefix": getattr(primer_pair, "amplicon_prefix", None),
                            "publishedAmpliconName": (
                                f"{getattr(primer_pair, 'amplicon_prefix', None)}_{getattr(primer_pair, 'amplicon_number', None)}"
                                if getattr(primer_pair, "amplicon_prefix", None) is not None
                                else None
                            ),
                        }
                    )
                    break
            evaluations.append(
                LegacyCandidateEvaluation(
                    candidate_id,
                    pool,
                    "accepted",
                    None,
                    trim_gain,
                    full_gain,
                    tuple(conflicts),
                    len(edge_keys[pool]),
                    len(species[pool]),
                )
            )
        stages.append(
            LegacySalvageStage(
                f"legacy-salvage-{stage_index}",
                cutoff,
                tuple(accepted),
                tuple(evaluations),
                {
                    **rejections,
                    **(
                        {
                            "__meta__": {
                                "stopReason": "work-limit",
                                "evaluated": evaluated,
                                "notExplored": max(0, len(remaining) - len(rejected_this_stage)),
                            }
                        }
                        if stage_work_limited
                        else {}
                    ),
                },
                _coverage(panel, msa_dict),
                monotonic() - started,
            )
        )
    if work_limited:
        stop_reason = "work-limit"
    return LegacySalvageRun(
        options,
        strict_ids,
        strict_coverage,
        tuple(manifest_entries),
        tuple(stages),
        stop_reason,
    )


def validate_legacy_salvage(
    panel,
    msa_dict,
    run: LegacySalvageRun,
    *,
    strict_cutoff: float = -26.0,
) -> dict[str, Any]:
    """Freshly audit strict assignments and replay every accepted addition.

    Validation intentionally works from the recorded strict pool assignments,
    then replays the accepted sequence at each recorded cutoff.  It therefore
    catches a changed final pool assignment or an active dimer introduced after
    the salvage run, without rebuilding candidate generation.
    """
    violations: list[str] = []
    run.options.validate_for_strict_cutoff(float(strict_cutoff))
    if tuple(stage.cutoff for stage in run.stages) != tuple(run.options.thresholds):
        violations.append("salvage stage cutoff sequence changed")
    all_pairs = list(panel.all_primerpairs())
    final_counts = Counter(
        legacy_candidate_id(item, msa_dict[item.msa_index]) for item in all_pairs
    )
    pair_ids = {
        legacy_candidate_id(item, msa_dict[item.msa_index]): item for item in all_pairs
    }
    accepted_ids = [candidate_id for stage in run.stages for candidate_id in stage.accepted]
    if not set(run.strict_candidate_ids).issubset(pair_ids):
        violations.append("strict candidate was removed")
    if not set(accepted_ids).issubset(pair_ids):
        violations.append("accepted salvage candidate is missing")
    expected_counts = Counter(run.strict_candidate_ids)
    expected_counts.update(accepted_ids)
    if final_counts != expected_counts:
        violations.append("final panel contains an unrecorded candidate")

    strict_ids = set(run.strict_candidate_ids)
    expected_pool_lists: dict[str, list[int]] = {}
    for item in run.candidate_manifest:
        if item.get("status") == "strict-incumbent":
            expected_pool_lists.setdefault(item["candidateId"], []).append(item.get("pool"))
    if set(expected_pool_lists) != strict_ids:
        violations.append("strict pool snapshot is incomplete")

    final_pools = getattr(panel, "_pools", ())
    pools: list[list[Any]] = [[] for _ in range(panel.n_pools)]
    actual_strict_by_pool: list[Counter[str]] = [Counter() for _ in range(panel.n_pools)]
    for pool, pool_pairs in enumerate(final_pools):
        if pool >= panel.n_pools:
            violations.append("final panel has an out-of-range pool")
            continue
        for item in pool_pairs:
            candidate_id = legacy_candidate_id(item, msa_dict[item.msa_index])
            if candidate_id in strict_ids:
                actual_strict_by_pool[pool][candidate_id] += 1
                pools[pool].append(item)
    expected_strict_by_pool: list[Counter[str]] = [Counter() for _ in range(panel.n_pools)]
    for candidate_id, pool_list in expected_pool_lists.items():
        for pool in pool_list:
            if pool is None or pool < 0 or pool >= panel.n_pools:
                violations.append(candidate_id + " has an invalid strict pool assignment")
            else:
                expected_strict_by_pool[pool][candidate_id] += 1
    if actual_strict_by_pool != expected_strict_by_pool:
        violations.append("strict pool assignments changed")
    strict_snapshot = _coverage(SimplePanelView(pools), msa_dict)
    if strict_snapshot != run.strict_coverage:
        violations.append("strict coverage snapshot changed")

    edge_keys = [set() for _ in range(panel.n_pools)]
    species = [set() for _ in range(panel.n_pools)]
    replay_matches: list[set] = [set() for _ in range(panel.n_pools)]
    matchdb = getattr(panel, "_matchDB", None)
    config = getattr(panel, "config", None)
    if matchdb is not None and config is not None:
        for pool, pool_pairs in enumerate(pools):
            for incumbent in pool_pairs:
                replay_matches[pool].update(
                    incumbent.find_matches(
                        matchdb,
                        fuzzy=config.mismatch_fuzzy,
                        remove_expected=True,
                        kmersize=config.mismatch_kmersize,
                    )
                )
    seen_accepted: set[str] = set()
    for stage in run.stages:
        if stage.cutoff < run.options.floor:
            violations.append(stage.stage_id + " is below the hard floor")
        for evaluation in stage.evaluations:
            candidate = pair_ids.get(evaluation.candidate_id)
            if candidate is None or evaluation.pool is None:
                violations.append(stage.stage_id + " has an unresolved accepted candidate")
                continue
            if evaluation.candidate_id in seen_accepted:
                violations.append(evaluation.candidate_id + " was accepted more than once")
            seen_accepted.add(evaluation.candidate_id)
            pool = evaluation.pool
            if pool < 0 or pool >= panel.n_pools or candidate not in final_pools[pool]:
                violations.append(evaluation.candidate_id + " has the wrong pool")
                continue
            for incumbent in pools[pool]:
                if int(incumbent.msa_index) != int(candidate.msa_index):
                    continue
                incumbent_full = (
                    min(incumbent.fprimer.starts()), max(incumbent.rprimer.ends())
                )
                candidate_full = (
                    min(candidate.fprimer.starts()), max(candidate.rprimer.ends())
                )
                if incumbent_full[0] < candidate_full[1] and candidate_full[0] < incumbent_full[1]:
                    violations.append(evaluation.candidate_id + " violates overlap admission")
                    break
            if matchdb is not None and config is not None:
                matches = candidate.find_matches(
                    matchdb,
                    fuzzy=config.mismatch_fuzzy,
                    remove_expected=False,
                    kmersize=config.mismatch_kmersize,
                )
                if detect_new_products(
                    matches,
                    replay_matches[pool],
                    config.mismatch_product_size,
                ):
                    violations.append(evaluation.candidate_id + " violates MatchDB admission")
            actual_gain = _gain_from_snapshot(candidate, msa_dict, _coverage(SimplePanelView(pools), msa_dict))
            if actual_gain is None or actual_gain[0] < run.options.min_reference_gain:
                violations.append(evaluation.candidate_id + " has insufficient reference gain")
            elif evaluation.marginal_trimmed_gain != actual_gain[0] or evaluation.marginal_full_gain != actual_gain[1]:
                violations.append(evaluation.candidate_id + " has inconsistent marginal gain")
            existing_sequences = [sequence for item in pools[pool] for sequence in item.all_seqs()]
            if do_pool_interact(
                candidate.all_seq_bytes(),
                [sequence.encode() for sequence in existing_sequences],
                stage.cutoff,
            ):
                violations.append(evaluation.candidate_id + " violates active cutoff")
            conflicts, new_edges, new_species = _pool_risk(
                SimplePanelView(pools),
                candidate,
                pool,
                stage.cutoff,
                float(strict_cutoff),
                edge_keys[pool],
                species[pool],
            )
            edge_keys[pool].update(new_edges)
            species[pool].update(new_species)
            if len(edge_keys[pool]) > run.options.max_edges_per_pool:
                violations.append(evaluation.candidate_id + " exceeds edge budget")
            if len(species[pool]) > run.options.max_incident_species_per_pool:
                violations.append(evaluation.candidate_id + " exceeds species budget")
            if (len(edge_keys[pool]), len(species[pool])) != (
                evaluation.cumulative_edges,
                evaluation.cumulative_incident_species,
            ):
                violations.append(evaluation.candidate_id + " has inconsistent exposure totals")
            pools[pool].append(candidate)
            if matchdb is not None and config is not None:
                replay_matches[pool].update(
                    candidate.find_matches(
                        matchdb,
                        fuzzy=config.mismatch_fuzzy,
                        remove_expected=True,
                        kmersize=config.mismatch_kmersize,
                    )
                )
        expected_coverage = stage.coverage
        actual_coverage = _coverage(SimplePanelView(pools), msa_dict)
        if actual_coverage != expected_coverage:
            violations.append(stage.stage_id + " has inconsistent coverage")

    if seen_accepted != set(accepted_ids):
        violations.append("accepted salvage manifest does not match evaluations")
    return {
        "schemaVersion": "primalscheme3.legacy-dimer-salvage-validation/v1",
        "valid": not violations,
        "violations": violations,
        "strictCandidateIds": list(run.strict_candidate_ids),
        "strictPoolAssignments": {
            candidate_id: sorted(expected_pool_lists.get(candidate_id, ()))
            for candidate_id in sorted(strict_ids)
        },
        "coverage": _coverage(panel, msa_dict),
    }


class SimplePanelView:
    def __init__(self, pools):
        self._pools = pools

    def all_primerpairs(self):
        return [item for pool in self._pools for item in pool]
