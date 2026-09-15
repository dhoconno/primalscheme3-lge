#!/usr/bin/env python3
"""Freeze labelled MSA fixtures and run an explicit allele-coverage matrix.

The harness performs no scientific scoring.  It verifies the source bundle's
declared bytes, creates an immutable copy, and invokes matrix commands with
argument arrays while retaining a receipt for success and failure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
from datetime import UTC, datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

WORKFLOW = "allele-coverage-benchmark-fixture"
WORKFLOW_VERSION = "1"
CONTROL_DEFINITIONS = {
    "historical-independent-lge.2": {
        "class": "saved-lge-legacy-artifact",
        "grouping": "independent",
        "controlledComparison": False,
    },
    "historical-combined-lge.2": {
        "class": "saved-lge-legacy-artifact",
        "grouping": "combined",
        "controlledComparison": False,
    },
    "existing-lge.3": {
        "class": "saved-or-reproduced-lge.3-execution",
        "controlledComparison": False,
    },
    "upstream-original-3.3.0": {
        "class": "upstream-original-execution",
        "repository": "https://github.com/artic-network/primalscheme3.git",
        "version": "3.3.0",
        "expectedCommit": "60455e9",
        "controlledComparison": False,
    },
}


class VerificationError(ValueError):
    """A manifest verification failure with durable expected/observed evidence."""

    def __init__(self, message: str, evidence: dict[str, Any]):
        super().__init__(message)
        self.evidence = evidence


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def descriptor(path: Path, *, relative_to: Path | None = None) -> dict[str, Any]:
    resolved = path.resolve()
    shown = str(resolved if relative_to is None else resolved.relative_to(relative_to.resolve()))
    return {"path": shown, "sha256": sha256(resolved), "byteSize": resolved.stat().st_size}


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def safe_relative(value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe relative path in manifest: {value}")
    return path


def source_identity() -> dict[str, Any]:
    script = Path(__file__).resolve()
    identity: dict[str, Any] = {"script": descriptor(script)}
    try:
        root = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"], cwd=script.parent,
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root,
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--short"], cwd=root,
            capture_output=True, text=True, check=True,
        ).stdout.splitlines()
        identity.update({"repository": root, "commit": commit, "dirty": bool(status), "status": status})
    except (OSError, subprocess.CalledProcessError):
        identity["repository"] = None
    return identity


def runtime_identity() -> dict[str, Any]:
    return {
        "python": platform.python_version(),
        "pythonExecutable": str(Path(sys.executable).resolve()),
        "platform": platform.platform(),
        "condaPrefix": os.environ.get("CONDA_PREFIX"),
        "container": os.environ.get("container"),
    }


def load_and_verify_inputs(
    fixture_root: Path, matrix: dict[str, Any], attempted: list[dict[str, Any]]
) -> tuple[Path, list[dict[str, Any]], list[dict[str, Any]]]:
    analysis = (fixture_root / safe_relative(matrix["fixtureAnalysis"])).resolve()
    manifest_path = analysis / "manifest.json"
    attempted.append({**descriptor(manifest_path), "role": "fixture-manifest"})
    manifest = json.loads(manifest_path.read_text())
    artifacts = {item["relativePath"]: item for item in manifest.get("artifacts", [])}
    requested = set(matrix["inputLabels"])
    selected = [item for item in manifest.get("inputs", []) if item.get("label") in requested]
    found = {item["label"] for item in selected}
    if found != requested:
        raise ValueError(f"manifest labels not found: {sorted(requested - found)}")
    verified: list[dict[str, Any]] = [descriptor(manifest_path)]
    for source in selected:
        for relative_string in source.get("artifactPaths", []):
            declared = artifacts.get(relative_string)
            if declared is None:
                raise ValueError(f"input artifact absent from manifest: {relative_string}")
            relative = safe_relative(relative_string)
            path = analysis / relative
            actual = descriptor(path)
            evidence = {
                "path": str(path.resolve()),
                "manifestPath": relative_string,
                "expectedByteSize": declared["byteSize"],
                "observedByteSize": actual["byteSize"],
                "expectedSha256": declared["sha256"],
                "observedSha256": actual["sha256"],
            }
            attempted.append({**actual, "role": "attempted-scientific-input", "expected": declared})
            if actual["byteSize"] != declared["byteSize"]:
                raise VerificationError(f"size mismatch for {relative_string}", evidence)
            if actual["sha256"] != declared["sha256"]:
                raise VerificationError(f"checksum mismatch for {relative_string}", evidence)
            verified.append({**actual, "manifestPath": relative_string, "role": declared.get("role")})
    return analysis, selected, verified


def copy_snapshot(
    analysis: Path,
    selected: list[dict[str, Any]],
    verified: list[dict[str, Any]],
    output: Path,
) -> list[dict[str, Any]]:
    snapshot = output / "snapshot"
    snapshot.mkdir()
    copied: list[dict[str, Any]] = []
    for item in verified[1:]:
        relative = safe_relative(item["manifestPath"])
        destination = snapshot / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(analysis / relative, destination)
        current = descriptor(destination, relative_to=output)
        if current["sha256"] != item["sha256"] or current["byteSize"] != item["byteSize"]:
            raise RuntimeError(f"snapshot verification failed for {relative}")
        copied.append(current)
    label_map = {
        "schemaVersion": "allele-coverage-source-label-map/v1",
        "inputs": [
            {
                "label": item["label"],
                "sourceOccurrenceID": item["id"],
                "snapshotArtifactPaths": [
                    str(Path("snapshot") / safe_relative(path)) for path in item.get("artifactPaths", [])
                ],
            }
            for item in selected
        ],
    }
    label_path = snapshot / "source-row-label-map.json"
    write_json(label_path, label_map)
    copied.append(descriptor(label_path, relative_to=output))
    manifest_path = snapshot / "hash-manifest.json"
    write_json(
        manifest_path,
        {
            "schemaVersion": "allele-coverage-fixture-snapshot/v1",
            "createdAt": datetime.now(UTC).isoformat(),
            "artifacts": copied,
        },
    )
    copied.append(descriptor(manifest_path, relative_to=output))
    return copied


def freeze_saved_control(
    entry: dict[str, Any], fixture_root: Path, output: Path
) -> dict[str, Any]:
    control_dir = output / "controls" / entry["id"]
    artifact_dir = control_dir / "artifacts"
    artifact_dir.mkdir(parents=True)
    source_root = (fixture_root / safe_relative(entry["root"])).resolve()
    manifest_relative = safe_relative(entry["manifest"])
    manifest_path = source_root / manifest_relative
    manifest = json.loads(manifest_path.read_text())
    source_artifacts = [{**descriptor(manifest_path), "role": "control-manifest"}]
    frozen_artifacts: list[dict[str, Any]] = []
    for declared in manifest.get("artifacts", []):
        relative = safe_relative(declared["relativePath"])
        source = source_root / relative
        actual = descriptor(source)
        evidence = {
            "path": str(source.resolve()),
            "manifestPath": str(relative),
            "expectedByteSize": declared["byteSize"],
            "observedByteSize": actual["byteSize"],
            "expectedSha256": declared["sha256"],
            "observedSha256": actual["sha256"],
        }
        if actual["byteSize"] != declared["byteSize"]:
            raise VerificationError(f"size mismatch for control {entry['id']}: {relative}", evidence)
        if actual["sha256"] != declared["sha256"]:
            raise VerificationError(
                f"checksum mismatch for control {entry['id']}: {relative}", evidence
            )
        source_artifacts.append(actual)
        destination = artifact_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        frozen = descriptor(destination, relative_to=output)
        if frozen["sha256"] != declared["sha256"] or frozen["byteSize"] != declared["byteSize"]:
            raise RuntimeError(f"frozen control verification failed: {entry['id']}: {relative}")
        frozen_artifacts.append(frozen)
    manifest_destination = artifact_dir / manifest_relative
    manifest_destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(manifest_path, manifest_destination)
    frozen_artifacts.append(descriptor(manifest_destination, relative_to=output))
    receipt = {
        "schemaVersion": "allele-coverage-saved-control/v1",
        "controlIdentity": CONTROL_DEFINITIONS.get(entry["id"], entry.get("controlIdentity")),
        "status": "frozen",
        "sourceArtifacts": source_artifacts,
        "artifacts": frozen_artifacts,
    }
    write_json(control_dir / "receipt.json", receipt)
    return receipt


def write_pending_control(entry: dict[str, Any], output: Path) -> dict[str, Any]:
    control_dir = output / "controls" / entry["id"]
    control_dir.mkdir(parents=True)
    receipt = {
        "schemaVersion": "allele-coverage-pending-control/v1",
        "controlIdentity": CONTROL_DEFINITIONS.get(entry["id"], entry.get("controlIdentity")),
        "status": "pending-execution",
        "artifacts": [],
    }
    write_json(control_dir / "receipt.json", receipt)
    return receipt


def scientific_input_artifacts(
    inputs: list[dict[str, Any]], output: Path
) -> list[dict[str, Any]]:
    artifacts = []
    for item in inputs:
        for relative in item["snapshotArtifactPaths"]:
            if relative.endswith("primary.aligned.fasta"):
                path = output / relative
                artifacts.append(
                    {
                        **descriptor(path, relative_to=output),
                        "role": "scientific-input",
                        "label": item["label"],
                        "sourceOccurrenceID": item["sourceOccurrenceID"],
                    }
                )
    return artifacts


def expand_argv(
    entry: dict[str, Any],
    native_executable: Path,
    run_dir: Path,
    inputs: list[dict[str, Any]],
    output: Path,
) -> list[str]:
    executable = str(Path(entry.get("executable", native_executable)).absolute())
    values = {
        "native_executable": executable,
        "output": str((run_dir / "scientific-output").resolve()),
        "snapshot": str((output / "snapshot").resolve()),
    }
    argv = [executable]
    for value in entry.get("argv", []):
        if value == "{msa_args}":
            for artifact in scientific_input_artifacts(inputs, output):
                argv.extend(["--msa", str((output / artifact["path"]).resolve())])
        else:
            argv.append(str(value).format(**values))
    return argv


def validate_tool_identity(entry: dict[str, Any], probe: dict[str, Any]) -> dict[str, str]:
    expected = entry.get("expectedToolIdentity")
    required_expected = {"name", "version", "gitCommit"}
    if not isinstance(expected, dict) or not required_expected.issubset(expected):
        raise ValueError(
            f"{entry['id']} expectedToolIdentity requires name, version, and gitCommit"
        )
    tool = probe["tool"]
    tool_name = tool if isinstance(tool, str) else tool["name"]
    source = probe["source"]
    runtime = probe["runtime"]
    for key in ("gitCommit", "sourceDigest"):
        if not isinstance(source.get(key), str) or not source[key]:
            raise ValueError(f"{entry['id']} measured source requires nonempty {key}")
    if not isinstance(runtime.get("pythonExecutable"), str) or not runtime["pythonExecutable"]:
        raise ValueError(f"{entry['id']} measured runtime requires pythonExecutable")
    kernel = runtime.get("kernel")
    if not isinstance(kernel, dict) or any(not kernel.get(key) for key in ("system", "release", "version")):
        raise ValueError(f"{entry['id']} measured runtime requires complete kernel identity")
    dependencies = runtime.get("declaredRuntimeDependencies")
    if not isinstance(dependencies, list) or not dependencies or any(
        not isinstance(item, dict) or not item.get("distribution") or not item.get("version")
        for item in dependencies
    ):
        raise ValueError(f"{entry['id']} measured runtime requires dependency identities")
    observed = {
        "name": tool_name,
        "version": probe["toolVersion"],
        "gitCommit": probe["source"].get("gitCommit"),
    }
    for key, value in expected.items():
        actual = observed.get(key)
        if key == "gitCommit" and actual is not None:
            matches = actual.startswith(value)
        else:
            matches = actual == value
        if not matches:
            raise ValueError(
                f"{entry['id']} identity mismatch for {key}: expected {value!r}, "
                f"observed {actual!r}"
            )
    return {"name": observed["name"], "version": observed["version"]}


def option_values_match(argv_value: str, resolved_value: Any) -> bool:
    if isinstance(resolved_value, bool):
        return argv_value.lower() == str(resolved_value).lower()
    if isinstance(resolved_value, int | float):
        try:
            return Decimal(argv_value) == Decimal(str(resolved_value))
        except InvalidOperation:
            return False
    return argv_value == str(resolved_value)


def validate_option_contract(
    entry: dict[str, Any], argv: list[str], measured_tool_version: str
) -> None:
    contract = entry.get("optionContract")
    if not isinstance(contract, dict) or contract.get("schemaVersion") != (
        "allele-coverage-resolved-options/v1"
    ):
        raise ValueError(f"{entry['id']} requires the resolved-options/v1 contract")
    if contract.get("toolVersion") != measured_tool_version:
        raise ValueError(f"{entry['id']} option contract toolVersion mismatch")
    required = contract.get("requiredResolvedOptionKeys")
    option_map = contract.get("argvOptionMap")
    if not isinstance(required, list) or not required:
        raise ValueError(f"{entry['id']} option contract requires resolved keys")
    if not isinstance(option_map, dict) or not option_map:
        raise ValueError(f"{entry['id']} option contract requires argvOptionMap")
    resolved = entry["resolvedOptions"]
    missing = [key for key in required if key not in resolved]
    if missing:
        raise ValueError(f"{entry['id']} missing required resolved options: {missing}")
    for flag, key in option_map.items():
        if key not in required:
            raise ValueError(f"{entry['id']} argv mapping uses undeclared option {key}")
        positions = [index for index, value in enumerate(argv) if value == flag]
        if len(positions) != 1 or positions[0] + 1 >= len(argv):
            raise ValueError(f"{entry['id']} requires exactly one {flag} value")
        argv_value = argv[positions[0] + 1]
        if not option_values_match(argv_value, resolved[key]):
            raise ValueError(
                f"{entry['id']} argv mismatch for {key}: {argv_value!r} != {resolved[key]!r}"
            )


def validate_actual_resolved_options(
    entry: dict[str, Any], run_dir: Path, output: Path
) -> dict[str, Any] | None:
    declaration = entry.get("actualResolvedOptions")
    if declaration is None:
        return None
    path = run_dir / "scientific-output" / safe_relative(declaration["path"])
    payload: Any = json.loads(path.read_text())
    for key in declaration.get("jsonPath", []):
        payload = payload[key]
    if not isinstance(payload, dict):
        raise ValueError(f"{entry['id']} actual resolved options are not an object")
    required = entry["optionContract"]["requiredResolvedOptionKeys"]
    for key in required:
        if key not in payload or payload[key] != entry["resolvedOptions"][key]:
            raise ValueError(
                f"{entry['id']} actual resolved option mismatch for {key}: "
                f"{payload.get(key)!r} != {entry['resolvedOptions'][key]!r}"
            )
    return {**descriptor(path, relative_to=output), "jsonPath": declaration.get("jsonPath", [])}


def run_entry(
    entry: dict[str, Any],
    native_executable: Path,
    output: Path,
    inputs: list[dict[str, Any]],
    source_artifacts: list[dict[str, Any]],
) -> int:
    run_dir = output / "runs" / entry["id"]
    run_dir.mkdir(parents=True)
    stderr_path = run_dir / "stderr.txt"
    stdout_path = run_dir / "stdout.txt"
    argv = expand_argv(entry, native_executable, run_dir, inputs, output)
    probe_argv = [argv[0], *entry.get("identityProbeArgv", [])]
    consumed = scientific_input_artifacts(inputs, output)
    if not entry.get("optionsFullyResolved") or not entry.get("resolvedOptions"):
        raise ValueError(f"{entry['id']} requires complete resolvedOptions")
    started = time.perf_counter()
    completed: subprocess.CompletedProcess[str] | None = None
    probe: subprocess.CompletedProcess[str] | None = None
    probe_payload: dict[str, Any] = {}
    executable_artifact: dict[str, Any] | None = None
    tool_identity: dict[str, str] | None = None
    actual_options_artifact: dict[str, Any] | None = None
    status = "launch-error"
    error = ""
    try:
        executable_artifact = descriptor(Path(argv[0]))
        probe = subprocess.run(probe_argv, capture_output=True, text=True, check=False)
        if probe.returncode != 0:
            raise RuntimeError(
                f"identity probe exited {probe.returncode}: {probe.stderr.strip()}"
            )
        probe_payload = json.loads(probe.stdout)
        required = ("tool", "toolVersion", "source", "runtime")
        missing = [key for key in required if key not in probe_payload]
        if missing:
            raise ValueError(f"identity probe missing fields: {missing}")
        tool_identity = validate_tool_identity(entry, probe_payload)
        validate_option_contract(entry, argv, probe_payload["toolVersion"])
        completed = subprocess.run(argv, capture_output=True, text=True, check=False)
        status = "success" if completed.returncode == 0 else "failed"
        stdout_path.write_text(completed.stdout)
        stderr_path.write_text(completed.stderr)
        if completed.returncode == 0:
            actual_options_artifact = validate_actual_resolved_options(
                entry, run_dir, output
            )
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}\n"
        if completed is None:
            stdout_path.write_text("")
            stderr_path.write_text(error)
        else:
            status = "output-validation-error"
            stdout_path.write_text(completed.stdout)
            stderr_path.write_text(completed.stderr + error)
            declaration = entry.get("actualResolvedOptions")
            if declaration:
                actual_path = (
                    run_dir / "scientific-output" / safe_relative(declaration["path"])
                )
                if actual_path.is_file():
                    actual_options_artifact = {
                        **descriptor(actual_path, relative_to=output),
                        "jsonPath": declaration.get("jsonPath", []),
                    }
    elapsed = time.perf_counter() - started
    output_artifacts = [
        descriptor(stdout_path, relative_to=output),
        descriptor(stderr_path, relative_to=output),
    ]
    for path in sorted(run_dir.rglob("*")):
        if path.is_file() and path.name not in {"receipt.json", "stdout.txt", "stderr.txt"}:
            output_artifacts.append(descriptor(path, relative_to=output))
    receipt = {
        "schemaVersion": "allele-coverage-execution-receipt/v1",
        "workflow": WORKFLOW,
        "workflowVersion": WORKFLOW_VERSION,
        "controlIdentity": CONTROL_DEFINITIONS.get(entry["id"], entry.get("controlIdentity")),
        "status": status,
        "argv": argv,
        "command": shlex.join(argv),
        "identityProbe": {
            "argv": probe_argv,
            "exitStatus": None if probe is None else probe.returncode,
            "stdout": None if probe is None else probe.stdout,
            "stderr": None if probe is None else probe.stderr,
        },
        "toolIdentity": tool_identity,
        "resolvedOptions": entry["resolvedOptions"],
        "optionsResolution": "matrix-explicit-complete",
        "inputArtifacts": consumed,
        "sourceArtifacts": source_artifacts,
        "outputArtifacts": output_artifacts,
        "sourceIdentity": probe_payload.get("source"),
        "runtimeIdentity": probe_payload.get("runtime"),
        "executedToolArtifact": (
            None if executable_artifact is None else {**executable_artifact, "role": "executable"}
        ),
        "actualResolvedOptionsArtifact": actual_options_artifact,
        "harnessIdentity": {"source": source_identity(), "runtime": runtime_identity()},
        "exitStatus": None if completed is None else completed.returncode,
        "wallTimeSeconds": elapsed,
        "stderrPath": str(stderr_path.relative_to(output)),
        "stdoutPath": str(stdout_path.relative_to(output)),
    }
    write_json(run_dir / "receipt.json", receipt)
    if status == "output-validation-error":
        return 1
    return 1 if completed is None else completed.returncode


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", required=True, type=Path)
    parser.add_argument("--native-executable", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    return parser.parse_args(argv)


def execute(args: argparse.Namespace, argv: list[str]) -> int:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.mkdir(parents=True)
    stderr_path = output / "runner-stderr.txt"
    started = time.perf_counter()
    status = 1
    inputs: list[dict[str, Any]] = []
    outputs: list[dict[str, Any]] = []
    error = ""
    verification_failure: dict[str, Any] | None = None
    try:
        inputs.append({**descriptor(args.matrix.resolve()), "role": "matrix"})
        matrix = json.loads(args.matrix.read_text())
        if matrix.get("schemaVersion") != "allele-coverage-benchmark-matrix/v1":
            raise ValueError("unsupported benchmark matrix schemaVersion")
        analysis, selected, verified = load_and_verify_inputs(
            args.fixture_root.resolve(), matrix, inputs
        )
        outputs = copy_snapshot(analysis, selected, verified, output)
        label_map = json.loads((output / "snapshot" / "source-row-label-map.json").read_text())
        matrix_copy = output / "matrix.json"
        shutil.copyfile(args.matrix.resolve(), matrix_copy)
        if sha256(matrix_copy) != inputs[0]["sha256"]:
            raise RuntimeError("matrix copy verification failed")
        outputs.append(descriptor(matrix_copy, relative_to=output))
        write_json(output / "control-definitions.json", CONTROL_DEFINITIONS)
        outputs.append(descriptor(output / "control-definitions.json", relative_to=output))
        status = 0
        for entry in matrix.get("comparisons", []):
            if entry.get("kind") == "saved-artifact":
                freeze_saved_control(entry, args.fixture_root.resolve(), output)
                outputs.append(
                    descriptor(
                        output / "controls" / entry["id"] / "receipt.json",
                        relative_to=output,
                    )
                )
                continue
            if entry.get("kind") == "pending-execution":
                write_pending_control(entry, output)
                outputs.append(
                    descriptor(
                        output / "controls" / entry["id"] / "receipt.json",
                        relative_to=output,
                    )
                )
                continue
            if entry.get("kind") != "native-execution":
                raise ValueError(f"unknown comparison kind: {entry.get('kind')}")
            status = run_entry(entry, args.native_executable, output, label_map["inputs"], inputs)
            if status:
                break
    except VerificationError as exc:
        error = f"{type(exc).__name__}: {exc}\n"
        verification_failure = exc.evidence
        status = 1
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}\n"
        status = 1
    stderr_path.write_text(error)
    receipt = {
        "schemaVersion": "allele-coverage-runner-receipt/v1",
        "workflow": WORKFLOW,
        "workflowVersion": WORKFLOW_VERSION,
        "argv": argv,
        "command": shlex.join(argv),
        "resolvedOptions": {
            "fixtureRoot": str(args.fixture_root.resolve()),
            "nativeExecutable": str(args.native_executable.absolute()),
            "output": str(output),
            "matrix": str(args.matrix.resolve()),
        },
        "inputArtifacts": inputs,
        "outputArtifacts": outputs,
        "sourceIdentity": source_identity(),
        "runtimeIdentity": runtime_identity(),
        "exitStatus": status,
        "wallTimeSeconds": time.perf_counter() - started,
        "stderrPath": str(stderr_path.relative_to(output)),
        "verificationFailure": verification_failure,
    }
    write_json(output / "runner-receipt.json", receipt)
    if error:
        sys.stderr.write(error)
    return status


def main(argv: list[str] | None = None) -> int:
    raw = sys.argv[1:] if argv is None else argv
    args = parse_args(raw)
    full_argv = [str(Path(sys.argv[0]).resolve()), *raw]
    return execute(args, full_argv)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FileExistsError as exc:
        sys.stderr.write(f"{exc}\n")
        raise SystemExit(2) from None
