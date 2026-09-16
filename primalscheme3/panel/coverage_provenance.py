"""Runtime identity and final-byte provenance for native coverage outputs."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import logging
import os
import platform
import shlex
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Any

from primalscheme3.core.logger import close_owned_file_handlers

from .phase_scheduling import scheduling_policy

CAPABILITIES_SCHEMA = "primalscheme3.capabilities/v1"
PROVENANCE_SCHEMA = "primalscheme3.panel-provenance/v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _descriptor(path: Path, *, relative_to: Path | None = None) -> dict[str, Any]:
    return {
        "path": path.relative_to(relative_to).as_posix()
        if relative_to
        else str(path.resolve()),
        "sha256": _sha256(path),
        "size": path.stat().st_size,
    }


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def source_identity() -> dict[str, Any]:
    """Describe the editable source tree, rather than only the launcher."""

    root = _repo_root()
    paths = sorted((root / "primalscheme3").rglob("*.py"))
    if (root / "pyproject.toml").is_file():
        paths.append(root / "pyproject.toml")
    files = [_descriptor(path, relative_to=root) for path in sorted(paths)]
    semantic = json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        commit = None
    try:
        status_lines = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
    except (OSError, subprocess.SubprocessError):
        status_lines = []
    return {
        "kind": "editable-source-tree",
        "root": str(root),
        "gitCommit": commit,
        "gitDirty": bool(status_lines),
        "gitStatus": status_lines,
        "sourceDigest": hashlib.sha256(semantic).hexdigest(),
        "files": files,
        "build": _descriptor(root / "pyproject.toml", relative_to=root),
    }


def runtime_identity() -> dict[str, Any]:
    uname = platform.uname()
    native_kernels = []
    for package, distribution, modules in (
        (
            "primalschemers",
            "primalschemers",
            ("primalschemers", "primalschemers.primalschemers"),
        ),
        ("primer3", "primer3-py", ("primer3.thermoanalysis",)),
    ):
        kernel_files = []
        for module in modules:
            spec = importlib.util.find_spec(module)
            if spec is not None and spec.origin is not None:
                path = Path(spec.origin).resolve()
                if path.is_file():
                    kernel_files.append(_descriptor(path))
        native_kernels.append(
            {
                "package": package,
                "distribution": distribution,
                "version": version(distribution),
                "files": kernel_files,
            }
        )
    declared_dependencies = []
    for distribution in (
        "primer3-py",
        "networkx",
        "kaleido",
        "plotly",
        "tqdm",
        "typer",
        "dnaio",
        "matplotlib",
        "primalschemers",
        "numpy",
        "primalbedtools",
    ):
        declared_dependencies.append(
            {"distribution": distribution, "version": version(distribution)}
        )
    return {
        "pythonVersion": platform.python_version(),
        "pythonImplementation": platform.python_implementation(),
        "pythonExecutable": sys.executable,
        "pythonExecutableResolved": str(Path(sys.executable).resolve()),
        "pythonPrefix": sys.prefix,
        "pythonBasePrefix": sys.base_prefix,
        "platform": platform.platform(),
        "machine": uname.machine,
        "kernel": {
            "system": uname.system,
            "release": uname.release,
            "version": uname.version,
        },
        "nativeKernels": native_kernels,
        "declaredRuntimeDependencies": declared_dependencies,
    }


def capabilities_document() -> dict[str, Any]:
    return {
        "schemaVersion": CAPABILITIES_SCHEMA,
        "tool": "primalscheme3",
        "toolVersion": version("primalscheme3"),
        "selectionAlgorithms": ["legacy", "coverage", "allele-coverage"],
        "alleleCoverage": {
            "intendedProductPolicies": {
                "default": "exact-supported",
                "policies": {
                    "exact-supported": {"id": "exact-supported/v1"},
                    "concrete-designated-sites": {
                        "id": "concrete-designated-sites/v1",
                        "coverageCredit": 0,
                    },
                },
            },
            "secondaryProductPolicies": {
                "default": "ordered-disjoint-intended-sites",
                "policies": {
                    "ordered-disjoint-intended-sites": {
                        "id": "ordered-disjoint-intended-sites/v1"
                    },
                    "reject-secondary-products/v1": {
                        "id": "reject-secondary-products/v1"
                    },
                    "ordered-disjoint-concrete-designated-sites/v1": {
                        "id": "ordered-disjoint-concrete-designated-sites/v1",
                        "coverageCredit": 0,
                    },
                },
            },
            "phaseScheduling": {
                "default": "serial",
                "policies": {
                    name: scheduling_policy(name) for name in ("serial", "reserved")
                },
            },
            "sourceContractVersion": "3.3.0+lge.4",
            "algorithm": "bounded-allele-coverage/v1",
            "metric": "observed-allele-primer-trimmed/v1",
            "preset": "allele-balanced-v1",
            "catalogSchemaVersion": "primalscheme3.variant-catalog/v2",
            "configurationLedgerSchemaVersion": "primalscheme3.configuration-ledger/v2",
            "validationSchemaVersion": "primalscheme3.allele-panel-validation/v2",
            "profile": {
                "name": "allele-panel-v1",
                "specificityRevision": "selected-sites-v2",
            },
            "supportedScope": {
                "mode": "equal",
                "mapping": "first",
                "ampliconSizeMetric": "reference-span",
                "terminalGapPolicies": ["observed-only"],
                "linear": True,
                "freshDesign": True,
                "suppliedMsaSpecificity": True,
            },
        },
        "coverage": {
            "algorithm": "bounded-multistart-coverage/v1",
            "catalogSchemaVersion": "primalscheme3.coverage-catalog/v1",
            "optimizerSchemaVersion": "primalscheme3.panel-optimizer/v1",
            "validationSchemaVersion": "primalscheme3.panel-validation/v1",
            "provenanceSchemaVersion": PROVENANCE_SCHEMA,
            "profile": {"name": "panel-v1", "specificityRevision": "intended-sites-v1"},
            "supportedScope": {
                "mode": "equal",
                "mapping": "first",
                "ampliconSizeMetric": "reference-span",
                "terminalGapPolicies": ["legacy", "observed-only"],
                "linear": True,
                "freshDesign": True,
                "suppliedMsaSpecificity": True,
            },
        },
        "source": source_identity(),
        "runtime": runtime_identity(),
    }


def capture_execution_identity(input_paths: list[Path]) -> dict[str, Any]:
    """Snapshot executed code/runtime and historical input origins before work."""

    return {
        "source": source_identity(),
        "runtime": runtime_identity(),
        "inputs": [
            {
                "sourcePath": str(Path(path).resolve()),
                "sha256": _sha256(Path(path).resolve()),
                "size": Path(path).resolve().stat().st_size,
            }
            for path in input_paths
        ],
    }


def _flush_logger(logger: Any | None) -> None:
    handlers = list(logging.getLogger().handlers)
    if logger is not None:
        handlers.extend(logger.handlers)
    for handler in dict.fromkeys(handlers):
        try:
            handler.flush()
            console = getattr(handler, "console", None)
            file = getattr(console, "file", None)
            if file is not None:
                file.flush()
        except (AttributeError, OSError):
            continue
    if logger is not None:
        close_owned_file_handlers(logger)


def finalize_provenance(
    *,
    output_dir: Path,
    argv: list[str],
    resolved_options: dict[str, Any],
    inputs: list[dict[str, str]],
    started_at: float,
    ended_at: float,
    status: str,
    exit_status: int,
    stderr: str,
    scientific: dict[str, Any],
    logger: Any | None = None,
    execution_start: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Hash completed bytes and finally write the self-excluding record."""

    output_dir = Path(output_dir).resolve()
    _flush_logger(logger)
    if execution_start is None:
        execution_start = capture_execution_identity(
            [Path(item["sourcePath"]) for item in inputs]
        )
    end_source = source_identity()
    end_runtime = runtime_identity()
    source_changed = (
        execution_start["source"]["sourceDigest"] != end_source["sourceDigest"]
    )
    runtime_changed = execution_start["runtime"] != end_runtime
    input_records = []
    for index, item in enumerate(inputs):
        source = Path(item["sourcePath"]).resolve()
        stored = output_dir / item["storedPath"]
        start_input = execution_start["inputs"][index]
        if start_input["sourcePath"] != str(source):
            raise ValueError(
                "input provenance order disagrees with execution-start snapshot"
            )
        source_exists = source.is_file()
        source_end_sha = _sha256(source) if source_exists else None
        source_end_size = source.stat().st_size if source_exists else None
        record = {
            "sourcePath": str(source),
            "sourceAtStartSha256": start_input["sha256"],
            "sourceAtStartSize": start_input["size"],
            "sourceAtEndSha256": source_end_sha,
            "sourceAtEndSize": source_end_size,
            "sourceChangedDuringRun": (
                source_end_sha != start_input["sha256"]
                or source_end_size != start_input["size"]
            ),
            "storedPath": item["storedPath"],
            "sourceIndex": item.get("sourceIndex", index),
        }
        if stored.is_file():
            stored_descriptor = _descriptor(stored, relative_to=output_dir)
            record.update(
                {
                    "storedPath": stored_descriptor["path"],
                    "sha256": stored_descriptor["sha256"],
                    "size": stored_descriptor["size"],
                    "storedMissing": False,
                    "sourceAtStartMatchesStored": (
                        stored_descriptor["sha256"] == start_input["sha256"]
                        and stored_descriptor["size"] == start_input["size"]
                    ),
                }
            )
        else:
            record["storedMissing"] = True
        if status == "success" and (
            record.get("storedMissing") or not record.get("sourceAtStartMatchesStored")
        ):
            raise ValueError(
                "successful coverage provenance requires an exact execution-start input copy"
            )
        input_records.append(record)
    provenance_path = output_dir / "panel-provenance.json"
    outputs = [
        _descriptor(path, relative_to=output_dir)
        for path in sorted(output_dir.rglob("*"))
        if path.is_file() and path != provenance_path
    ]
    record = {
        "schemaVersion": PROVENANCE_SCHEMA,
        "tool": "primalscheme3 panel-create",
        "toolVersion": version("primalscheme3"),
        "command": {
            "argv": list(argv),
            "shell": shlex.join(argv),
            "workingDirectory": os.getcwd(),
        },
        "resolvedOptions": resolved_options,
        "source": execution_start["source"],
        "sourceAtEnd": end_source,
        "sourceChangedDuringRun": source_changed,
        "runtime": execution_start["runtime"],
        "runtimeAtEnd": end_runtime,
        "runtimeChangedDuringRun": runtime_changed,
        "inputs": input_records,
        "outputs": outputs,
        "scientific": scientific,
        "wallSeconds": max(0.0, ended_at - started_at),
        "status": status,
        "exitStatus": exit_status,
        "stderr": stderr,
        "selfHashPolicy": "panel-provenance.json is excluded to avoid a self-hash cycle",
    }
    if status == "success" and (source_changed or runtime_changed):
        raise ValueError(
            "executed source or runtime identity changed during coverage workflow"
        )
    provenance_path.write_text(
        json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
    )
    return record
