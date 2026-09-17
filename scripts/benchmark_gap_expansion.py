#!/usr/bin/env python3
"""Prepare (or explicitly run) the bounded legacy gap-expansion benchmark."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from scripts.benchmark_gap_completion import (
    EXPECTED_INPUTS,
    MSAS,
    PARENT,
    combined_geometry,
    geometry,
    stamp,
)

BASE = PARENT.parents[1]
PRIOR = BASE / "mhc-a1-a2-gap-completion-benchmark-pending-v1/gap-completed"
OUT = BASE / "mhc-a1-a2-gap-expansion-benchmark-pending-v1"
EXE = Path(
    "/Users/dho/Documents/lungfish-genome-explorer/.worktrees/"
    "allele-aware-primalscheme/.venv/bin/primalscheme3"
)
ANCHOR_BOUND = 2000
PAIR_BOUND = 1000
SCIENTIFIC_FILES = (
    "primer.bed",
    "amplicon.bed",
    "primertrim.amplicon.bed",
    "reference.fasta",
    "config.json",
)


COMMON = [
    "panel-create",
    "--msa",
    str(MSAS[0]),
    "--msa",
    str(MSAS[1]),
    "--selection-algorithm",
    "legacy",
    "--mode",
    "equal",
    "--amplicon-size",
    "200",
    "--amplicon-size-min",
    "150",
    "--amplicon-size-max",
    "250",
    "--n-pools",
    "2",
    "--ncores",
    "1",
    "--no-high-gc",
    "--min-base-freq",
    "0",
    "--mapping",
    "first",
    "--dimer-score",
    "-26",
    "--terminal-gap-policy",
    "legacy",
    "--use-matchdb",
    "--use-tm",
    "--offline-plots",
]


def parent_payload() -> dict[str, dict[str, object]]:
    return {name: stamp(PARENT / name) for name in SCIENTIFIC_FILES}


def prior_payload() -> dict[str, dict[str, object]]:
    return {name: stamp(PRIOR / name) for name in SCIENTIFIC_FILES}


def commands() -> list[str]:
    return [
        str(EXE),
        *COMMON,
        "--output",
        str(OUT / "expanded"),
        "--gap-completion-parent",
        str(PARENT),
        "--gap-expansion",
        "bounded",
        "--gap-expansion-max-anchors-per-msa",
        str(ANCHOR_BOUND),
        "--gap-expansion-max-pairs-per-msa",
        str(PAIR_BOUND),
    ]


def prepare() -> Path:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite {OUT}")
    if not PARENT.is_dir() or not PRIOR.is_dir():
        raise RuntimeError("missing preserved parent or prior gap-completion output")
    for path, (digest, size) in zip(MSAS, EXPECTED_INPUTS, strict=True):
        item = stamp(path)
        if item["sha256"] != digest or item["sizeBytes"] != size:
            raise RuntimeError(f"input receipt mismatch: {path}")
    # The prior result is an independent comparison artifact, never an input to
    # native expansion. Verify the same scientific payload is available before launch.
    parent_refs = geometry(PARENT)["referenceLengths"]
    prior_refs = geometry(PRIOR)["referenceLengths"]
    if parent_refs != prior_refs:
        raise RuntimeError("prior gap-completion references differ from parent")
    OUT.mkdir(parents=True)
    manifest = {
        "schemaVersion": "lge.mhc-a1-a2-gap-expansion-benchmark/v1",
        "status": "prepared; native launch requires explicit --run",
        "parent": str(PARENT.resolve()),
        "parentScientificPayloadBefore": parent_payload(),
        "priorGapCompletion": str(PRIOR.resolve()),
        "priorScientificPayload": prior_payload(),
        "inputs": [stamp(path) for path in MSAS],
        "argv": commands(),
        "command": shlex.join(commands()),
        "limits": {"wallTimeSeconds": 600, "rssBytes": 8 * 1024**3, "ncores": 1},
        "expansionBounds": {
            "maxAnchorsPerMsa": ANCHOR_BOUND,
            "maxGeometryEligibleSameObservedRowPairChecksPerMsa": PAIR_BOUND,
        },
        "preparedAt": datetime.now(timezone.utc).isoformat(),
    }
    path = OUT / "preparation-manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def _tree_rss(proc, psutil_module) -> int:
    process = psutil_module.Process(proc.pid)
    return sum(
        item.memory_info().rss
        for item in [process, *process.children(recursive=True)]
    )


def run() -> dict[str, object]:
    argv = commands()
    output = OUT / "expanded"
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    stdout_path = OUT / "expanded.stdout.log"
    stderr_path = OUT / "expanded.stderr.log"
    stdout_file = stdout_path.open("w")
    stderr_file = stderr_path.open("w")
    proc = subprocess.Popen(
        argv,
        cwd=Path(__file__).resolve().parents[1],
        stdout=stdout_file,
        stderr=stderr_file,
        start_new_session=True,
    )
    peak_rss = 0
    deadline = started + 600
    kill_reason = None
    try:
        import psutil
    except ImportError:
        psutil = None
    try:
        while proc.poll() is None:
            if psutil is not None:
                try:
                    peak_rss = max(peak_rss, _tree_rss(proc, psutil))
                except psutil.Error:
                    pass
            else:
                rows = subprocess.check_output(
                    ["ps", "-axo", "pid=,pgid=,rss="], text=True
                ).splitlines()
                pgid = str(proc.pid)
                rss = sum(
                    int(row.split()[2])
                    for row in rows
                    if len(row.split()) >= 3 and row.split()[1] == pgid
                ) * 1024
                peak_rss = max(peak_rss, rss)
            if time.monotonic() >= deadline:
                kill_reason = "timeout"
                os.killpg(proc.pid, 9)
                break
            if peak_rss > 8 * 1024**3:
                kill_reason = "rss"
                os.killpg(proc.pid, 9)
                break
            time.sleep(0.25)
    finally:
        proc.wait()
        stdout_file.close()
        stderr_file.close()
    elapsed = time.monotonic() - started
    stdout = stdout_path.read_text()
    stderr = stderr_path.read_text()
    receipt: dict[str, object] = {
        "schemaVersion": "lge.mhc-a1-a2-gap-expansion-benchmark-receipt/v1",
        "argv": argv,
        "command": shlex.join(argv),
        "startedAt": started_at,
        "wallTimeSeconds": elapsed,
        "exitStatus": proc.returncode,
        "rssBytes": peak_rss,
        "killReason": kill_reason,
        "stdoutSha256": hashlib.sha256(stdout.encode()).hexdigest(),
        "stderr": stderr,
        "runtime": {"python": sys.version, "platform": platform.platform()},
        "inputs": [stamp(path) for path in MSAS],
        "parentScientificPayloadBefore": parent_payload(),
        "priorScientificPayload": prior_payload(),
        "parentGeometry": geometry(PARENT),
        "priorFollowupGeometry": geometry(PRIOR),
        "priorCombinedGeometry": combined_geometry(PRIOR, PARENT),
    }
    if output.is_dir():
        receipt["followupGeometry"] = geometry(output)
        receipt["combinedGeometry"] = combined_geometry(output, PARENT)
        receipt["outputs"] = [
            stamp(path) for path in sorted(output.rglob("*")) if path.is_file()
        ]
        for name in ("gap-expansion.json", "gap-completion.json", "gap-completion-coverage.json"):
            path = output / name
            if path.is_file():
                receipt[name.replace(".json", "")] = json.loads(path.read_text())
    receipt["parentScientificPayloadAfter"] = parent_payload()
    receipt["parentUnchanged"] = receipt["parentScientificPayloadBefore"] == receipt["parentScientificPayloadAfter"]
    path = OUT / "expanded.benchmark-receipt.json"
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    manifest = prepare()
    print(json.dumps({"prepared": str(manifest), "command": shlex.join(commands())}, indent=2))
    if args.run:
        result = run()
        print(json.dumps(result, indent=2))
        if result["exitStatus"] != 0 or not result.get("parentUnchanged", False):
            raise SystemExit(1)
