#!/usr/bin/env python3
"""Run a predeclared cached allele matrix only after a separate source freeze.

No scientific kernels live here: native panel-create produces each fresh panel,
and native panel-audit independently checks it. Inputs are original frozen MSAs.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

TOOL = "benchmark-cached-allele-panel"
VERSION = "1"


def descriptor(path, root=None):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(path.relative_to(Path(root).resolve()) if root else path),
        "sha256": digest.hexdigest(),
        "size": path.stat().st_size,
    }


def write(path, value):
    Path(path).write_text(
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )


def read(path):
    with (
        gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()
    ) as stream:
        return json.load(stream)


def contained(root, name):
    root = Path(root).resolve()
    relative = Path(name)
    path = (root / relative).resolve()
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or not path.is_relative_to(root)
    ):
        raise ValueError("unsafe relative input path: " + str(name))
    return path


def snapshot_inputs(root):
    root = Path(root).resolve()
    label_path = root / "snapshot/source-row-label-map.json"
    manifest_path = root / "snapshot/hash-manifest.json"
    labels, manifest = read(label_path), read(manifest_path)
    if labels.get("schemaVersion") != "allele-coverage-source-label-map/v1":
        raise ValueError("unsupported snapshot label map")
    declared = {d["path"]: d for d in manifest["artifacts"]}

    def verify(path):
        actual = descriptor(path, root)
        expected = declared.get(actual["path"])
        if (
            expected is None
            or actual["sha256"] != expected["sha256"]
            or actual["size"] != expected.get("byteSize", expected.get("size"))
        ):
            raise ValueError("snapshot hash/size mismatch: " + actual["path"])
        return descriptor(path)

    verify(label_path)
    records = []
    identities = set()
    for item in labels["inputs"]:
        identity = item["sourceOccurrenceID"]
        if identity in identities:
            raise ValueError("duplicate input occurrence")
        identities.add(identity)
        choices = [
            name
            for name in item["snapshotArtifactPaths"]
            if name.startswith("snapshot/source-inputs/")
            and name.endswith("/alignment/primary.aligned.fasta")
        ]
        if len(choices) != 1:
            raise ValueError("each input requires exactly one original aligned FASTA")
        path = contained(root, choices[0])
        records.append(
            verify(path)
            | {
                "label": item["label"],
                "sourceOccurrenceID": identity,
                "sourceIndex": len(records),
            }
        )
    if not records:
        raise ValueError("snapshot contains no original inputs")
    return records, [
        descriptor(label_path),
        descriptor(manifest_path),
        *[{k: r[k] for k in ("path", "sha256", "size")} for r in records],
    ]


def cache_inputs(cache):
    cache = Path(cache).resolve()
    manifest = read(cache / "manifest.json")
    if manifest.get("schemaVersion") != "primalscheme3.discovery-cache/v1":
        raise ValueError("unsupported cache manifest")
    paths = [
        cache / "manifest.json",
        contained(cache, manifest["exportProvenancePath"]),
    ]
    for relative, expected in manifest["artifacts"].items():
        path = contained(cache, relative)
        actual = descriptor(path)
        if any(actual[k] != expected[k] for k in ("sha256", "size")):
            raise ValueError("cache artifact differs: " + relative)
        paths.append(path)
    catalog = read(contained(cache, manifest["catalogPath"]))
    settings = json.loads(catalog["resolved_config_json"])
    receipt = read(contained(cache, manifest["sourcePanelProvenance"]))
    diagnostics = {
        "generatedSites": len(catalog["sites"]),
        "generated17to18Sites": sum(
            17 <= len(s["sequence"]) <= 18 for s in catalog["sites"]
        ),
        "normalGeneratedSites": sum(
            "normal" in s["generated_profile_ids"] for s in catalog["sites"]
        ),
        "normalGeneratedBelow19": sum(
            "normal" in s["generated_profile_ids"] and len(s["sequence"]) < 19
            for s in catalog["sites"]
        ),
        "normalProfileProjectionExcludedSites": sum(
            "normal" not in s["generated_profile_ids"] for s in catalog["sites"]
        ),
    }
    return (
        manifest,
        settings,
        receipt,
        [descriptor(p) for p in sorted(set(paths))],
        diagnostics,
    )


def arm_manifest():
    base = {
        "candidate_profiles": "union",
        "variant_selection": "subsets",
        "phase_scheduling": "serial",
        "optimizer_seed": 0,
        "optimizer_starts": 1,
        "optimizer_repair_rounds": 2,
        "optimizer_time_limit": 120,
        "specificity_terminal_k": 17,
        "secondary_product_policy": "ordered-disjoint-intended-sites",
        "salvage": "off",
        "primary_tier": "strict",
    }
    arms = {}

    def add(
        name, changes=None, scope="same immutable original inputs and cache", repeat=1
    ):
        arms[name] = {
            "options": base | (changes or {}),
            "scope": scope,
            "repeat": repeat,
        }

    add(
        "normal-full",
        {"candidate_profiles": "normal", "variant_selection": "full-cloud"},
    )
    add("union-full", {"variant_selection": "full-cloud"})
    add("normal-subsets", {"candidate_profiles": "normal"})
    add("union-subsets-repair0", {"optimizer_repair_rounds": 0})
    add("union-subsets-120")
    add("union-subsets-600", {"optimizer_time_limit": 600})
    add(
        "union-subsets-salvage",
        {
            "salvage": "bounded",
            "salvage_max_stages": 3,
            "salvage_max_edges_per_pool": 8,
            "salvage_max_oligos_per_pool": 4,
            "salvage_time_limit": 60,
        },
    )
    for k in (17, 19):
        add(
            f"normal-k{k}",
            {"candidate_profiles": "normal", "specificity_terminal_k": k},
            scope="narrower matched subset: normal-only catalog, every generated primer >=19 bases; not union-wide k comparison",
        )
    add(
        "union-secondary-reject",
        {"secondary_product_policy": "reject-secondary-products/v1"},
    )
    add(
        "fixed-work-repeat",
        {
            "optimizer_time_limit": 3600,
            "work_frontier_candidates": 8,
            "work_construction_candidate_attempts": 16,
            "work_repair_candidate_probes_per_round": 8,
            "work_repair_neighborhoods_per_round": 8,
            "work_repair_trials_per_round": 8,
            "work_pool_lookahead_candidates": 2,
            "work_cleanup_moves_per_round": 8,
            "work_families_per_refresh": 4,
            "subset_expansion_limit": 16,
            "subset_beam_width": 4,
        },
        repeat=2,
    )
    return arms


def run_child(argv, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    process = None
    usage = None
    error = ""
    interrupted = False
    code = 127
    with (
        (directory / "stdout.txt").open("wb") as stdout,
        (directory / "stderr.txt").open("wb") as stderr,
    ):
        try:
            process = subprocess.Popen(
                argv, stdout=stdout, stderr=stderr, start_new_session=True
            )
            if hasattr(os, "wait4"):
                _, status, usage = os.wait4(process.pid, 0)
                code = os.waitstatus_to_exitcode(status)
                process.returncode = code
            else:
                code = process.wait()
        except BaseException as exc:
            error = f"{type(exc).__name__}: {exc}"
            interrupted = isinstance(exc, KeyboardInterrupt | SystemExit)
            if process is not None and process.returncode is None:
                os.killpg(process.pid, signal.SIGKILL)
                if hasattr(os, "wait4"):
                    _, status, usage = os.wait4(process.pid, 0)
                    code = os.waitstatus_to_exitcode(status)
                    process.returncode = code
                else:
                    code = process.wait()
    raw = usage.ru_maxrss if usage is not None else None
    unit = (
        "bytes"
        if sys.platform == "darwin"
        else "KiB"
        if sys.platform.startswith("linux")
        else "platform-native"
    )
    peak = {
        "value": raw,
        "unit": unit,
        "platform": platform.system(),
        "scope": "wait4 direct-child rusage; may include reaped descendants per OS; not a simultaneous process-tree aggregate"
        if usage is not None
        else "unavailable on this platform",
    }
    result = {
        "argv": list(argv),
        "shell": shlex.join(argv),
        "workingDirectory": os.getcwd(),
        "exitStatus": code,
        "wallSeconds": time.monotonic() - start,
        "stderr": (directory / "stderr.txt").read_text(errors="replace") + error,
        "interrupted": interrupted,
        "peakRSS": peak,
        "outputs": [
            descriptor(directory / name, directory)
            for name in ("stdout.txt", "stderr.txt")
        ],
    }
    write(directory / "command.json", result)
    return result


def probe(executable, directory):
    command = run_child([str(executable), "--capabilities-json"], directory)
    if command["exitStatus"] != 0:
        raise ValueError("native capability probe failed: " + command["stderr"])
    capabilities = read(Path(directory) / "stdout.txt")
    if "allele-coverage" not in capabilities.get("selectionAlgorithms", []):
        raise ValueError("native executable lacks allele-coverage")
    for key in ("source", "runtime", "toolVersion"):
        if key not in capabilities:
            raise ValueError("native identity incomplete: " + key)
    return capabilities


def wrapper_identity():
    return {
        "script": descriptor(__file__),
        "python": platform.python_version(),
        "executable": descriptor(sys.executable),
        "platform": platform.platform(),
        "condaPrefix": os.environ.get("CONDA_PREFIX"),
        "container": os.environ.get("container"),
    }


def frozen_identity(capabilities, executable, inputs, wrapper):
    return {
        "nativeSourceDigest": capabilities["source"]["sourceDigest"],
        "nativeRuntime": capabilities["runtime"],
        "nativeVersion": capabilities["toolVersion"],
        "executable": descriptor(executable),
        "inputs": [dict(item) for item in inputs],
        "wrapper": wrapper,
    }


def completed_fixed_work(metadata):
    options = metadata["options"]
    return (
        metadata.get("fixed_work_completed", True) is True
        and metadata.get("stop_reason") == "completed"
        and metadata.get("completed_starts") == options["starts"]
        and metadata.get("completed_repair_rounds")
        == options["starts"] * options["repair_rounds"]
        and metadata.get("completed_seed_modes") == metadata.get("effective_seed_modes")
    )


def command_options(options):
    result = []
    for name, value in sorted(options.items()):
        result.extend(["--" + name.replace("_", "-"), str(value)])
    return result


def summarize(panel):
    optimizer = read(panel / "panel-optimizer.json")
    coverage = read(panel / "row-coverage.json")
    validation = read(panel / "panel-validation.json")
    selected = validation.get("selected_sites", [])
    stages = []
    for stage in optimizer["stages"]:
        metadata = stage["optimizer"]
        stages.append(
            {
                "stage": stage["stage_id"],
                "coverage": stage["coverage"],
                "assignments": stage["assignments"],
                "work": metadata["work"],
                "proposalWork": metadata["proposal_work"],
                "searchSeconds": metadata["search_seconds"],
                "stopReason": metadata["stop_reason"],
                "completedStarts": metadata["completed_starts"],
                "completedRepairs": metadata["completed_repair_rounds"],
                "fixedWorkCompleted": completed_fixed_work(metadata),
                "workTruncated": metadata["work_truncated"],
            }
        )
    class_rows = coverage["classes"]
    pool_statistics = {}
    for pool, details in validation.get("pools", {}).items():
        scores = [edge["score"] for edge in details.get("dimer_edges", ())]
        minimum = min(scores) if scores else None
        pool_statistics[pool] = {
            "physicalEdges": len(scores),
            "minimumDimerScore": minimum,
            "marginAboveActiveCutoff": None
            if minimum is None
            else minimum - validation["stage_policy"]["active_cutoff"],
            "strictViolatingEdges": details.get("violating_edge_count"),
            "incidentSpecies": details.get("incident_species_count"),
        }
    return {
        "metric": optimizer["metric"],
        "primaryTier": optimizer["primaryTier"],
        "coverage": coverage,
        "dropoutClasses": sum(
            c["covered_count"] == 0 for c in class_rows if c["observed_count"] > 0
        ),
        "selectedSiteCount": len({s["site_id"] for s in selected}),
        "sequencePoolInstances": len({(s["sequence"], s["pool"]) for s in selected}),
        "assignmentCount": len(validation.get("assignments", [])),
        "allowedSecondaryProductCount": len(
            validation.get("allowed_secondary_products", [])
        ),
        "uncertaintyBlockCount": len(validation.get("uncertainty_blocks", [])),
        "stages": stages,
        "poolStatistics": pool_statistics,
        "salvage": read(panel / "salvage-comparison.json"),
    }


def fixed_signature(summary):
    return [
        {k: s[k] for k in ("stage", "coverage", "assignments", "work", "proposalWork")}
        for s in summary["stages"]
    ]


def check_native_receipt(path, before, *, options=None):
    receipt = read(path)
    if (
        receipt.get("status") != "success"
        or receipt.get("exitStatus") != 0
        or receipt.get("sourceChangedDuringRun") is not False
        or receipt.get("runtimeChangedDuringRun") is not False
    ):
        raise ValueError("native receipt not stable successful: " + str(path))
    if (
        receipt["source"]["sourceDigest"] != before["source"]["sourceDigest"]
        or receipt["runtime"] != before["runtime"]
    ):
        raise ValueError("native execution disagrees with frozen identity")
    if options is not None:
        # Check effective controls, not merely their presence in requested argv.
        config = read(path.parent / "config.json")
        resolved = json.loads(config["allele_options_json"])
        for key, value in options.items():
            if resolved.get(key) != value:
                raise ValueError("native resolved option mismatch: " + key)
        return resolved
    return receipt.get("resolvedOptions")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "cache", "native-executable", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--arm", action="append", choices=tuple(arm_manifest()))
    parser.add_argument("--source-freeze", type=Path)
    parser.add_argument("--freeze-only", action="store_true")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    snapshot = args.snapshot.resolve()
    cache = args.cache.resolve()
    executable = args.native_executable.absolute()
    if (
        output.exists()
        or output.is_relative_to(snapshot)
        or output.is_relative_to(cache)
    ):
        parser.error("output must be new and outside immutable inputs")
    output.mkdir(parents=True)
    started = time.monotonic()
    inputs = []
    results = []
    records = []
    catalog_diagnostics = {}
    repeat_report = None
    failure = ""
    before = None
    after = None
    wrapper = wrapper_identity()
    frozen = None
    try:
        inputs = [
            descriptor(path)
            for path in (
                snapshot / "snapshot/source-row-label-map.json",
                snapshot / "snapshot/hash-manifest.json",
                cache / "manifest.json",
            )
            if path.is_file()
        ]
        if not args.freeze_only and args.source_freeze is None:
            raise ValueError(
                "matrix requires a prior --source-freeze; use --freeze-only first"
            )
        records, snapshot_descriptors = snapshot_inputs(snapshot)
        manifest, settings, cache_receipt, cache_descriptors, catalog_diagnostics = (
            cache_inputs(cache)
        )
        inputs = snapshot_descriptors + cache_descriptors
        before = probe(executable, output / "native-before")
        identity = frozen_identity(before, executable, inputs, wrapper)
        chosen = args.arm or list(arm_manifest())
        if len(set(chosen)) != len(chosen):
            raise ValueError("duplicate matrix arm")
        arms = {name: arm_manifest()[name] for name in chosen}
        write(
            output / "arm-manifest.json",
            {"schemaVersion": "primalscheme3.cached-matrix/v1", "arms": arms},
        )
        if args.freeze_only:
            frozen = {
                "schemaVersion": "primalscheme3.cached-matrix-freeze/v1",
                "identity": identity,
            }
            write(output / "freeze.json", frozen)
        else:
            frozen = read(args.source_freeze)
            freeze_receipt = read(args.source_freeze.parent / "provenance.json")
            freeze_descriptor = descriptor(
                args.source_freeze, args.source_freeze.parent
            )
            if (
                freeze_receipt.get("status") != "success"
                or freeze_receipt.get("exitStatus") != 0
                or freeze_descriptor not in freeze_receipt.get("outputs", [])
            ):
                raise ValueError(
                    "source freeze lacks a successful byte-bound preparation receipt"
                )
            inputs.extend(
                [
                    descriptor(args.source_freeze),
                    descriptor(args.source_freeze.parent / "provenance.json"),
                ]
            )
            if (
                frozen.get("schemaVersion") != "primalscheme3.cached-matrix-freeze/v1"
                or frozen["identity"] != identity
            ):
                raise ValueError(
                    "source/runtime/input identity differs from prior freeze"
                )
            common = {
                "selection_algorithm": "allele-coverage",
                "mode": "equal",
                "amplicon_size": cache_receipt["resolvedOptions"]["amplicon_size"],
                "amplicon_size_min": settings["amplicon_size_min"],
                "amplicon_size_max": settings["amplicon_size_max"],
                "ncores": 1,
                "n_pools": 2,
                "min_base_freq": settings["minimum_frequency"]["normal"],
                "discovery_length_mode": settings["discovery_length_mode"],
                "coverage_target": 0.95,
                "reuse_discovery": str(cache),
            }
            if len(set(settings["minimum_frequency"].values())) != 1:
                raise ValueError("matrix requires shared profile minimum frequency")
            for name, arm in arms.items():
                if name in ("normal-k17", "normal-k19") and (
                    settings["profiles"]["normal"]["primer_size_min"] < 19
                    or catalog_diagnostics["normalGeneratedBelow19"]
                ):
                    raise ValueError(
                        "matched k17/k19 requires normal catalog primers >=19"
                    )
                for repeat in range(arm["repeat"]):
                    run_dir = output / "runs" / f"{name}-{repeat + 1}"
                    run_dir.mkdir(parents=True)
                    panel = run_dir / "panel"
                    audit = run_dir / "audit"
                    options = common | arm["options"]
                    command = [
                        str(executable),
                        "panel-create",
                        *command_options(options),
                        "--output",
                        str(panel),
                    ]
                    for record in records:
                        command.extend(["--msa", record["path"]])
                    current = probe(executable, run_dir / "native-before")
                    if (
                        frozen_identity(
                            current, executable, identity["inputs"], wrapper_identity()
                        )
                        != identity
                    ):
                        raise ValueError(
                            "native source/runtime changed before arm " + name
                        )
                    native = run_child(command, run_dir / "execution")
                    result = {
                        "arm": name,
                        "repeat": repeat + 1,
                        "scope": arm["scope"],
                        "requestedOptions": options,
                        "native": native,
                        "status": "failed",
                    }
                    results.append(result)
                    try:
                        if native["exitStatus"] != 0:
                            raise ValueError("native panel failed")
                        audited = run_child(
                            [
                                str(executable),
                                "panel-audit",
                                "--bundle",
                                str(panel),
                                "--output",
                                str(audit),
                            ],
                            run_dir / "audit-execution",
                        )
                        result["audit"] = audited
                        result["resolvedOptions"] = check_native_receipt(
                            panel / "panel-provenance.json",
                            before,
                            options={
                                k: v
                                for k, v in options.items()
                                if k not in ("selection_algorithm", "mode")
                            },
                        )
                        if (
                            audited["exitStatus"] != 0
                            or not read(audit / "validation.json")["valid"]
                        ):
                            raise ValueError("fresh native audit failed")
                        check_native_receipt(audit / "provenance.json", before)
                        result["summary"] = summarize(panel)
                        result["status"] = "success"
                    except Exception as exc:
                        result["error"] = str(exc)
                    finally:
                        result["boundReceipts"] = [
                            descriptor(p, output)
                            for p in (
                                panel / "panel-provenance.json",
                                panel / "config.json",
                                panel / "panel-optimizer.json",
                                audit / "provenance.json",
                                audit / "validation.json",
                            )
                            if p.is_file()
                        ]
                        write(run_dir / "receipt.json", result)
                    if native["interrupted"] or result.get("audit", {}).get(
                        "interrupted"
                    ):
                        raise KeyboardInterrupt("matrix interrupted")
            repeated = [r for r in results if r["arm"] == "fixed-work-repeat"]
            repeat_report = None
            if repeated:
                complete = len(repeated) == 2 and all(
                    r["status"] == "success"
                    and all(s["fixedWorkCompleted"] for s in r["summary"]["stages"])
                    for r in repeated
                )
                same = complete and fixed_signature(
                    repeated[0]["summary"]
                ) == fixed_signature(repeated[1]["summary"])
                repeat_report = {
                    "bothCompleted": complete,
                    "sameScientificResultAndWork": same,
                    "interpretation": "completed fixed-work comparison"
                    if complete
                    else "inconclusive: incomplete/interrupted work, not a fixed-work reproduction",
                }
                if not same:
                    failure = "fixed-work repeat did not complete identically"
            write(
                output / "summary.json",
                {
                    "arms": results,
                    "fixedWorkComparison": repeat_report,
                    "inputOrder": records,
                    "catalogDiagnostics": catalog_diagnostics,
                    "limitations": [
                        "No global-optimum claim.",
                        "Normal k17/k19 is a narrower matched catalog subset, not a union-wide comparison.",
                        "Peak RSS is OS wait4 direct-child usage, not simultaneous descendant sum.",
                    ],
                },
            )
            if any(r["status"] != "success" for r in results):
                failure = failure or "one or more arms failed"
    except (Exception, KeyboardInterrupt) as exc:
        failure = f"{type(exc).__name__}: {exc}"
    finally:
        try:
            if before is not None:
                after = probe(executable, output / "native-after")
                if (
                    after["source"]["sourceDigest"] != before["source"]["sourceDigest"]
                    or after["runtime"] != before["runtime"]
                    or descriptor(executable) != identity["executable"]
                ):
                    failure = failure or "native source/runtime changed during matrix"
            if wrapper_identity() != wrapper:
                failure = failure or "wrapper identity changed during matrix"
            for item in inputs:
                if descriptor(item["path"]) != item:
                    failure = failure or "immutable input changed: " + item["path"]
        except Exception as exc:
            failure = failure or f"final identity verification failed: {exc}"
        if not args.freeze_only:
            summary = (
                read(output / "summary.json")
                if (output / "summary.json").is_file()
                else {
                    "arms": results,
                    "fixedWorkComparison": repeat_report,
                    "inputOrder": records,
                    "catalogDiagnostics": catalog_diagnostics,
                }
            )
            summary.update(status="failed" if failure else "success", error=failure)
            write(output / "summary.json", summary)
        receipt = {
            "schemaVersion": "primalscheme3.cached-matrix-provenance/v1",
            "tool": TOOL,
            "toolVersion": VERSION,
            "argv": [sys.executable, str(Path(__file__).resolve()), *argv],
            "shell": shlex.join([sys.executable, str(Path(__file__).resolve()), *argv]),
            "workingDirectory": os.getcwd(),
            "requestedOptions": vars(args),
            "wrapperAtStart": wrapper,
            "wrapperAtEnd": wrapper_identity(),
            "nativeAtStart": before,
            "nativeAtEnd": after,
            "inputs": inputs,
            "outputs": [
                descriptor(p, output)
                for p in sorted(output.rglob("*"))
                if p.is_file() and p != output / "provenance.json"
            ],
            "wallSeconds": time.monotonic() - started,
            "exitStatus": 1 if failure else 0,
            "status": "failed" if failure else "success",
            "stderr": failure,
            "selfHashPolicy": "root provenance.json excluded",
        }
        receipt["requestedOptions"] = {
            k: str(v) if isinstance(v, Path) else v
            for k, v in receipt["requestedOptions"].items()
        }
        write(output / "provenance.json", receipt)
        if failure:
            print(failure, file=sys.stderr)
    return 1 if failure else 0


if __name__ == "__main__":
    raise SystemExit(main())
