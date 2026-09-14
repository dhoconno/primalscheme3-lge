#!/usr/bin/env python3
"""Run the abstract-oracle coverage search benchmark with full provenance."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import shlex
import sys
import time
import traceback
import tracemalloc
from collections import Counter
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET_LENGTH = 300
NOMINAL_SIZE = 200
MINIMUM_SIZE = 150
MAXIMUM_SIZE = 280
POOL_COUNT = 1
COVERAGE_TARGET = 0.9
TEMPLATES = {
    "A": (20, 280),
    "B": (0, 150),
    "C": (150, 300),
}

sys.path.insert(0, str(ROOT))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def descriptor(path: Path, *, relative_to: Path | None = None) -> dict:
    return {
        "path": (
            path.relative_to(relative_to).as_posix()
            if relative_to is not None
            else str(path.resolve())
        ),
        "sha256": sha256(path),
        "size": path.stat().st_size,
    }


def scientific_source_inputs() -> list[Path]:
    return [
        ROOT / "primalscheme3/core/config.py",
        ROOT / "primalscheme3/panel/coverage_search.py",
        ROOT / "primalscheme3/panel/coverage_types.py",
        ROOT / "primalscheme3/panel/coverage_validation.py",
    ]


def build_fixture(
    *,
    target_count: int,
    seed: int,
    starts: int,
    repair_rounds: int,
    time_limit: float,
    search_limits: dict,
) -> dict:
    targets = []
    candidates = []
    conflicts = []
    for index in range(target_count):
        target_id = f"target-{index:03d}"
        targets.append({"id": target_id, "referenceLength": TARGET_LENGTH})
        ids = {}
        for label, (start, end) in TEMPLATES.items():
            candidate_id = f"{target_id}-{label}"
            ids[label] = candidate_id
            candidates.append(
                {
                    "id": candidate_id,
                    "targetID": target_id,
                    "fullInterval": [start, end],
                    "interiorInterval": [start + 1, end - 1],
                    "forwardOligos": [candidate_id + "-F"],
                    "reverseOligos": [candidate_id + "-R"],
                }
            )
        conflicts.extend([[ids["A"], ids["B"]], [ids["A"], ids["C"]]])
    return {
        "schemaVersion": "primalscheme3.abstract-search-fixture/v1",
        "scope": (
            "ABSTRACT ORACLE ONLY: no biological sequence, MSA, thermochemistry, "
            "or specificity evidence."
        ),
        "targetCount": target_count,
        "targetLength": TARGET_LENGTH,
        "candidateTemplate": {
            "greedy": list(TEMPLATES["A"]),
            "left": list(TEMPLATES["B"]),
            "right": list(TEMPLATES["C"]),
        },
        "abstractConstraints": {
            "nPools": POOL_COUNT,
            "ampliconSize": NOMINAL_SIZE,
            "ampliconSizeMin": MINIMUM_SIZE,
            "ampliconSizeMax": MAXIMUM_SIZE,
            "allCandidatesIntrinsicallyValid": True,
            "conflictRule": "strict full-interval overlap within the same target",
            "interTargetConflicts": False,
            "maxAmplicons": None,
            "maxAmpliconsPerTarget": None,
        },
        "searchLimits": dict(search_limits),
        "runs": {
            "greedy": {
                "metric": "full-span",
                "coverageTarget": COVERAGE_TARGET,
                "seed": seed,
                "starts": starts,
                "repairRounds": 0,
                "timeLimitSeconds": time_limit,
            },
            "repaired": {
                "metric": "full-span",
                "coverageTarget": COVERAGE_TARGET,
                "seed": seed,
                "starts": starts,
                "repairRounds": repair_rounds,
                "timeLimitSeconds": time_limit,
            },
        },
        "targets": targets,
        "candidates": candidates,
        "conflicts": conflicts,
    }


def build_catalog(fixture: dict):
    from primalscheme3.panel.coverage_types import Candidate, Catalog, Target

    targets = tuple(
        Target(
            id=item["id"],
            source_msa_index=index,
            occurrence=0,
            row_ids=(),
            rows=(),
            reference_sequence="",
            reference_length=item["referenceLength"],
            mapping=(),
            ref_to_alignment=(),
        )
        for index, item in enumerate(fixture["targets"])
    )
    candidates = tuple(
        Candidate(
            id=item["id"],
            target_id=item["targetID"],
            full_interval=tuple(item["fullInterval"]),
            interior_interval=tuple(item["interiorInterval"]),
            forward_oligos=tuple(item["forwardOligos"]),
            reverse_oligos=tuple(item["reverseOligos"]),
            joint_rows=(),
            unknown_rows=(),
            row_product_spans=(),
        )
        for item in fixture["candidates"]
    )
    semantic_bytes = json.dumps(fixture, sort_keys=True, separators=(",", ":")).encode()
    return Catalog(
        targets=targets,
        candidates=candidates,
        semantic_digest=hashlib.sha256(semantic_bytes).hexdigest(),
        resolved_config_json=json.dumps(
            fixture["abstractConstraints"], sort_keys=True, separators=(",", ":")
        ),
        source_mapping=(),
        _native_pair_items=(),
    )


class AbstractGraphOracle:
    """Fixture-defined unary and pair predicates; never a scientific oracle."""

    def __init__(self, conflicts):
        self.conflicts = {frozenset(pair) for pair in conflicts}
        self.counters = Counter(candidateCalls=0, pairCalls=0)

    def candidate_valid(self, candidate_id: str) -> bool:
        self.counters["candidateCalls"] += 1
        return True

    def conflict(self, first: str, second: str) -> bool:
        self.counters["pairCalls"] += 1
        return frozenset((first, second)) in self.conflicts


def literal_union(intervals, reference_length: int) -> dict:
    bits = bytearray(reference_length)
    for start, end in intervals:
        assert 0 <= start <= end <= reference_length
        bits[start:end] = b"\x01" * (end - start)
    bitstring = "".join("1" if bit else "0" for bit in bits)
    return {
        "bitset": bitstring,
        "bitsetLength": len(bitstring),
        "bitsetSha256": hashlib.sha256(bitstring.encode()).hexdigest(),
        "coveredBases": bitstring.count("1"),
        "fraction": bitstring.count("1") / reference_length,
    }


def independently_validate(catalog, assignments, conflicts, *, pool_count: int) -> dict:
    conflict_set = {frozenset(pair) for pair in conflicts}
    violations = []
    seen = set()
    selected = {assignment.candidate_id: assignment for assignment in assignments}
    for assignment in assignments:
        if assignment.candidate_id not in catalog.candidate_by_id:
            violations.append(
                {"reason": "unknown-candidate", "candidateID": assignment.candidate_id}
            )
            continue
        if assignment.candidate_id in seen:
            violations.append(
                {
                    "reason": "duplicate-candidate",
                    "candidateID": assignment.candidate_id,
                }
            )
        seen.add(assignment.candidate_id)
        if type(assignment.pool) is not int or not 0 <= assignment.pool < pool_count:
            violations.append(
                {
                    "reason": "invalid-pool",
                    "candidateID": assignment.candidate_id,
                    "pool": assignment.pool,
                }
            )
    selected_ids = sorted(selected)
    for index, first in enumerate(selected_ids):
        for second in selected_ids[index + 1 :]:
            if (
                selected[first].pool == selected[second].pool
                and frozenset((first, second)) in conflict_set
            ):
                violations.append(
                    {
                        "reason": "same-pool-conflict",
                        "first": first,
                        "second": second,
                        "pool": selected[first].pool,
                    }
                )
    per_target = {}
    for target in sorted(catalog.targets, key=lambda item: item.id):
        target_candidates = [
            catalog.candidate_by_id[assignment.candidate_id]
            for assignment in assignments
            if assignment.candidate_id in catalog.candidate_by_id
            and catalog.candidate_by_id[assignment.candidate_id].target_id == target.id
        ]
        for candidate in target_candidates:
            span = candidate.full_interval[1] - candidate.full_interval[0]
            if not MINIMUM_SIZE <= span <= MAXIMUM_SIZE:
                violations.append(
                    {
                        "reason": "span-out-of-bounds",
                        "candidateID": candidate.id,
                        "span": span,
                    }
                )
        per_target[target.id] = {
            "selectedCandidateIDs": sorted(
                candidate.id for candidate in target_candidates
            ),
            "fullReferenceUnion": literal_union(
                [candidate.full_interval for candidate in target_candidates],
                target.reference_length,
            ),
            "interiorReferenceUnion": literal_union(
                [candidate.interior_interval for candidate in target_candidates],
                target.reference_length,
            ),
        }
    return {
        "scope": "independent abstract interval/graph/pool validation",
        "valid": not violations,
        "violations": violations,
        "selectedCount": len(assignments),
        "uniqueSelectedCount": len(seen),
        "perTarget": per_target,
        "totalFullCoveredBases": sum(
            item["fullReferenceUnion"]["coveredBases"] for item in per_target.values()
        ),
        "totalInteriorCoveredBases": sum(
            item["interiorReferenceUnion"]["coveredBases"]
            for item in per_target.values()
        ),
    }


def peak_rss_record() -> dict:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "value": value,
        "unit": "bytes" if platform.system() == "Darwin" else "kilobytes",
        "kind": "process cumulative maximum resident set size",
    }


def run_search(catalog, profile, fixture, options_record: dict) -> dict:
    from primalscheme3.panel.coverage_search import SearchOptions, search_assignments

    options = SearchOptions(
        metric=options_record["metric"],
        coverage_target=options_record["coverageTarget"],
        seed=options_record["seed"],
        starts=options_record["starts"],
        repair_rounds=options_record["repairRounds"],
        time_limit=options_record["timeLimitSeconds"],
    )
    oracle = AbstractGraphOracle(fixture["conflicts"])
    tracemalloc.reset_peak()
    begin = time.perf_counter()
    search_result = search_assignments(catalog, profile, options, oracle)
    elapsed = time.perf_counter() - begin
    _, python_peak = tracemalloc.get_traced_memory()
    independent = independently_validate(
        catalog,
        search_result.assignments,
        fixture["conflicts"],
        pool_count=POOL_COUNT,
    )
    if search_result.metadata["validation"]["scope"] != "abstract-oracle":
        raise RuntimeError("search did not retain its abstract-oracle scope")
    if not search_result.metadata["validation"]["valid"] or not independent["valid"]:
        raise RuntimeError("search returned an invalid abstract assignment")
    if search_result.metadata["search_limits"] != fixture["searchLimits"]:
        raise RuntimeError("search work limits differ from the declared fixture")
    return {
        "label": "greedy" if options.repair_rounds == 0 else "repaired",
        "options": search_result.metadata["options"],
        "assignments": [
            {"candidateID": item.candidate_id, "pool": item.pool}
            for item in search_result.assignments
        ],
        "searchMetadata": search_result.metadata,
        "oracleCounters": dict(oracle.counters),
        "independentValidation": independent,
        "elapsedSeconds": elapsed,
        "peakPythonAllocatedBytes": python_peak,
        "processPeakRss": peak_rss_record(),
    }


def execute(args: argparse.Namespace, argv: list[str]) -> int:
    from primalscheme3.core.config import Config
    from primalscheme3.panel.coverage_provenance import (
        runtime_identity,
        source_identity,
    )
    from primalscheme3.panel.coverage_search import LIMITS
    from primalscheme3.panel.coverage_types import Assignment
    from primalscheme3.panel.coverage_validation import ConstraintProfile

    output = args.output.resolve()
    output.mkdir(parents=False)
    fixture_path = output / "fixture.json"
    result_path = output / "result.json"
    provenance_path = output / "provenance.json"
    stdout_path = output / "stdout.txt"
    stderr_path = output / "stderr.txt"
    clock = time.perf_counter()
    source_inputs = [Path(__file__).resolve(), *scientific_source_inputs()]
    input_start = [descriptor(path) for path in source_inputs]
    source_start = source_identity()
    runtime_start = runtime_identity()
    record = {
        "schemaVersion": "primalscheme3.abstract-search-diagnostic-provenance/v1",
        "workflow": "synthetic-coverage-search-benchmark",
        "workflowVersion": "1",
        "toolVersion": version("primalscheme3"),
        "scientificScope": (
            "ABSTRACT ORACLE ONLY: no biological sequences, MSAs, chemistry, "
            "specificity, or panel-v1 scientific validation."
        ),
        "command": {
            "argv": argv,
            "shell": shlex.join(argv),
            "workingDirectory": os.getcwd(),
            "stdoutPath": str(stdout_path),
            "stderrPath": str(stderr_path),
        },
        "resolvedOptions": {
            "output": str(output),
            "targetCount": args.target_count,
            "targetLength": TARGET_LENGTH,
            "nominalSize": NOMINAL_SIZE,
            "minimumSize": MINIMUM_SIZE,
            "maximumSize": MAXIMUM_SIZE,
            "poolCount": POOL_COUNT,
            "coverageTarget": COVERAGE_TARGET,
            "seed": args.seed,
            "starts": args.starts,
            "greedyRepairRounds": 0,
            "repairedRepairRounds": args.repair_rounds,
            "timeLimitSecondsPerRun": args.time_limit,
        },
        "startedAt": datetime.now(UTC).isoformat(),
        "outputDirectory": str(output),
        "source": source_start,
        "runtime": runtime_start,
        "inputsAtStart": input_start,
        "outputs": [],
    }
    status = 1
    stdout = ""
    stderr = ""
    tracemalloc.start()
    try:
        fixture = build_fixture(
            target_count=args.target_count,
            seed=args.seed,
            starts=args.starts,
            repair_rounds=args.repair_rounds,
            time_limit=args.time_limit,
            search_limits=dict(LIMITS),
        )
        fixture_path.write_text(json.dumps(fixture, indent=2, sort_keys=True) + "\n")
        catalog = build_catalog(fixture)
        profile = ConstraintProfile.from_config(
            Config(
                amplicon_size=NOMINAL_SIZE,
                amplicon_size_min=MINIMUM_SIZE,
                amplicon_size_max=MAXIMUM_SIZE,
                n_pools=POOL_COUNT,
                mismatch_product_size=2000,
            )
        )
        greedy = run_search(catalog, profile, fixture, fixture["runs"]["greedy"])
        repaired = run_search(catalog, profile, fixture, fixture["runs"]["repaired"])
        expected_assignments = tuple(
            Assignment(f"target-{index:03d}-{label}", 0)
            for index in range(args.target_count)
            for label in ("B", "C")
        )
        expected_validation = independently_validate(
            catalog,
            expected_assignments,
            fixture["conflicts"],
            pool_count=POOL_COUNT,
        )
        if not expected_validation["valid"]:
            raise RuntimeError("declared feasible alternative is invalid")
        repaired_ids = {item["candidateID"] for item in repaired["assignments"]}
        expected_ids = {assignment.candidate_id for assignment in expected_assignments}
        result = {
            "schemaVersion": "primalscheme3.abstract-search-diagnostic/v1",
            "scope": (
                "ABSTRACT ORACLE ONLY: interval unions and an explicit conflict graph; "
                "not biological evidence and not panel-v1 scientific validation."
            ),
            "sourceSearchCore": {
                "function": "primalscheme3.panel.coverage_search.search_assignments",
                "algorithm": "bounded-multistart-coverage/v1",
                "toolVersion": version("primalscheme3"),
            },
            "fixture": {**descriptor(fixture_path), "path": "fixture.json"},
            "catalogSemanticDigest": catalog.semantic_digest,
            "abstractProfile": fixture["abstractConstraints"],
            "defaultWorkLimits": dict(LIMITS),
            "knownFeasibleAlternative": {
                "description": (
                    "Select the disjoint left [0,150) and right [150,300) "
                    "candidate for every target in pool zero."
                ),
                "independentValidation": expected_validation,
            },
            "runs": [greedy, repaired],
            "comparison": {
                "greedyTotalFullCoveredBases": greedy["independentValidation"][
                    "totalFullCoveredBases"
                ],
                "repairedTotalFullCoveredBases": repaired["independentValidation"][
                    "totalFullCoveredBases"
                ],
                "expectedFeasibleTotalFullCoveredBases": args.target_count
                * TARGET_LENGTH,
                "repairedReachedExpectedFeasibleAlternative": repaired_ids
                == expected_ids,
            },
            "provenance": "provenance.json",
        }
        result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        stdout = (
            json.dumps(
                {
                    "comparison": result["comparison"],
                    "runs": [
                        {
                            "label": run["label"],
                            "selectedCount": len(run["assignments"]),
                            "elapsedSeconds": run["elapsedSeconds"],
                            "completedStarts": run["searchMetadata"][
                                "completed_starts"
                            ],
                            "completedRepairRounds": run["searchMetadata"][
                                "completed_repair_rounds"
                            ],
                            "repairsAccepted": run["searchMetadata"][
                                "repairs_accepted"
                            ],
                            "stopReason": run["searchMetadata"]["stop_reason"],
                            "work": run["searchMetadata"]["work"],
                            "peakPythonAllocatedBytes": run["peakPythonAllocatedBytes"],
                            "processPeakRss": run["processPeakRss"],
                        }
                        for run in result["runs"]
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        status = 0
    except BaseException:
        stderr = traceback.format_exc()
    finally:
        source_end = source_identity()
        runtime_end = runtime_identity()
        input_end = [descriptor(path) for path in source_inputs]
        record.update(
            {
                "endedAt": datetime.now(UTC).isoformat(),
                "wallSeconds": time.perf_counter() - clock,
                "sourceAtEnd": source_end,
                "sourceChangedDuringRun": source_start != source_end,
                "runtimeAtEnd": runtime_end,
                "runtimeChangedDuringRun": runtime_start != runtime_end,
                "inputsAtEnd": input_end,
                "inputsChangedDuringRun": input_start != input_end,
                "processPeakRss": peak_rss_record(),
            }
        )
        if status == 0 and any(
            record[key]
            for key in (
                "sourceChangedDuringRun",
                "runtimeChangedDuringRun",
                "inputsChangedDuringRun",
            )
        ):
            stderr += "source, runtime, or inputs changed during run\n"
            status = 1
        stdout_path.write_text(stdout)
        stderr_path.write_text(stderr)
        record["exitStatus"] = status
        record["status"] = "success" if status == 0 else "failure"
        record["stderr"] = stderr
        record["outputs"] = [
            descriptor(path, relative_to=output)
            for path in sorted(output.rglob("*"))
            if path.is_file() and path != provenance_path
        ]
        record["selfHashPolicy"] = (
            "provenance.json is excluded to avoid a self-hash cycle"
        )
        provenance_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        tracemalloc.stop()
    if stdout:
        sys.stdout.write(stdout)
    if stderr:
        sys.stderr.write(stderr)
    return status


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Run the real coverage search core against a declared abstract "
            "interval/conflict fixture; this is not a biological benchmark."
        )
    )
    result.add_argument("--output", required=True, type=Path)
    result.add_argument("--target-count", type=int, default=100)
    result.add_argument("--seed", type=int, default=0)
    result.add_argument("--starts", type=int, default=1)
    result.add_argument("--repair-rounds", type=int, default=1)
    result.add_argument("--time-limit", type=float, default=60.0)
    return result


def main() -> int:
    argument_parser = parser()
    args = argument_parser.parse_args()
    if args.output.exists():
        argument_parser.error(f"output already exists: {args.output}")
    if args.target_count <= 0:
        argument_parser.error("--target-count must be positive")
    if args.starts <= 0:
        argument_parser.error("--starts must be positive")
    if args.repair_rounds <= 0:
        argument_parser.error("--repair-rounds must be positive")
    if not args.time_limit > 0:
        argument_parser.error("--time-limit must be positive")
    argv = [sys.executable, *sys.argv]
    return execute(args, argv)


if __name__ == "__main__":
    raise SystemExit(main())
