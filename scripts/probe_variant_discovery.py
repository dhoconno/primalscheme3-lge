"""Development-only bounded discovery probe, with reproducibility receipt."""

from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import json
import os
import resource
import shlex
import subprocess
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def descriptor(path):
    path = Path(path).resolve()
    return {
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "size": path.stat().st_size,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--msa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--anchor-start", type=int, default=100)
    parser.add_argument("--anchor-count", type=int, default=10)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument(
        "--history-backend", choices=("jsonl", "sqlite"), default="jsonl"
    )
    args = parser.parse_args()
    if args.anchor_start < 0 or args.anchor_count < 1 or args.cores < 1:
        parser.error("positive count/cores and nonnegative start required")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    sources = []
    history = None
    receipt = {
        "workflow": "development-variant-discovery-probe",
        "version": 1,
        "argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "command": shlex.join(
            [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]
        ),
        "cwd": str(Path.cwd()),
        "explicitOptions": {
            k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
        },
        "inputs": [],
        "sourceFiles": [],
        "exitStatus": 1,
    }
    try:
        with (
            (args.output / "stderr.txt").open("w") as stderr,
            contextlib.redirect_stderr(stderr),
        ):
            try:
                sources = [
                    descriptor(p)
                    for p in sorted((ROOT / "primalscheme3").rglob("*.py"))
                ]
                sources.append(descriptor(__file__))
                from primalscheme3.panel.coverage_provenance import runtime_identity

                runtime = runtime_identity()
                runtime["pythonFile"] = descriptor(sys.executable)
                runtime["condaPrefix"] = os.environ.get("CONDA_PREFIX")
                receipt.update(
                    inputs=[descriptor(args.msa)],
                    sourceFiles=sources,
                    runtime=runtime,
                    gitCommit=subprocess.check_output(
                        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
                    ).strip(),
                    gitStatus=subprocess.check_output(
                        ["git", "status", "--porcelain"], cwd=ROOT, text=True
                    ),
                )
                from primalscheme3.core.config import (
                    AmpliconSizeMetric,
                    Config,
                    TerminalGapPolicy,
                )
                from primalscheme3.core.mapping import create_mapping
                from primalscheme3.core.msa import parse_msa
                from primalscheme3.panel.coverage_discovery import (
                    build_variant_catalog,
                    variant_targets,
                )
                from primalscheme3.panel.coverage_history import (
                    CoverageHistory,
                    SQLiteCoverageHistory,
                )

                config = Config(
                    terminal_gap_policy=TerminalGapPolicy.OBSERVED_ONLY,
                    amplicon_size_metric=AmpliconSizeMetric.REFERENCE_SPAN,
                    amplicon_size=200,
                    amplicon_size_min=150,
                    amplicon_size_max=250,
                    ncores=args.cores,
                )
                receipt["resolvedOptions"] = config.items()
                array, _ = parse_msa(args.msa)
                mapping, _ = create_mapping(array)
                targets = variant_targets(
                    {
                        "probe": SimpleNamespace(
                            array=array, _mapping_array=mapping, msa_index=0
                        )
                    }
                )
                indexes = (
                    None
                    if args.full
                    else (
                        list(
                            range(
                                args.anchor_start, args.anchor_start + args.anchor_count
                            )
                        ),
                        list(
                            range(
                                args.anchor_start + 150,
                                args.anchor_start + 150 + args.anchor_count,
                            )
                        ),
                    )
                )
                receipt["resolvedOptions"].update(
                    discovery_length_mode="first-compatible",
                    candidate_profiles="union",
                    indexes=indexes,
                )
                history_class = (
                    SQLiteCoverageHistory
                    if args.history_backend == "sqlite"
                    else CoverageHistory
                )
                history = history_class(
                    args.output / "history", run_id="development-discovery-probe"
                )
                catalog = build_variant_catalog(
                    targets, config, indexes=indexes, history=history
                )
                with gzip.open(args.output / "catalog.json.gz", "wt") as handle:
                    json.dump(catalog.to_dict(), handle, sort_keys=True)
                summary = {
                    "targets": len(catalog.targets),
                    "sites": len(catalog.sites),
                    "selectableSites": len(catalog.selectable_site_ids),
                    "families": len(catalog.families),
                    "evidence": len(history.evidence),
                    "assessments": len(history.assessments),
                    "events": len(history.events),
                }
                (args.output / "summary.json").write_text(
                    json.dumps(summary, indent=2) + "\n"
                )
                receipt["exitStatus"] = 0
            except BaseException:
                traceback.print_exc()
    finally:
        if history is not None and hasattr(history, "close"):
            history.close()
        receipt["wallSeconds"] = time.monotonic() - start
        receipt["peakRSSRaw"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        receipt["peakRSSScope"] = "parent-process-only; excludes discovery worker RSS"
        receipt["peakRSSUnit"] = "bytes" if sys.platform == "darwin" else "KiB"
        try:
            receipt["sourceUnchanged"] = bool(sources) and all(
                descriptor(s["path"])["sha256"] == s["sha256"] for s in sources
            )
            receipt["inputUnchanged"] = (
                bool(receipt["inputs"]) and descriptor(args.msa) == receipt["inputs"][0]
            )
            runtime_files = (
                [
                    receipt["runtime"]["pythonFile"],
                    *(
                        f
                        for k in receipt["runtime"]["nativeKernels"]
                        for f in k["files"]
                    ),
                ]
                if "runtime" in receipt
                else []
            )
            receipt["runtimeUnchanged"] = bool(runtime_files) and all(
                descriptor(f["path"])["sha256"] == f["sha256"] for f in runtime_files
            )
            if not all(
                receipt[k]
                for k in ("sourceUnchanged", "inputUnchanged", "runtimeUnchanged")
            ):
                receipt["exitStatus"] = 1
                receipt["integrityError"] = (
                    "Source, input or runtime missing or changed during probe"
                )
            receipt["outputs"] = [
                descriptor(p) for p in sorted(args.output.rglob("*")) if p.is_file()
            ]
        except BaseException as error:
            receipt["exitStatus"] = 1
            receipt["receiptError"] = repr(error)
        (args.output / "provenance.json").write_text(
            json.dumps(receipt, indent=2, default=str) + "\n"
        )
    return receipt["exitStatus"]


if __name__ == "__main__":
    raise SystemExit(main())
