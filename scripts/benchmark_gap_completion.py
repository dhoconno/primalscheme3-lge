#!/usr/bin/env python3
"""Prepare (or explicitly run) the two-pool legacy gap-completion benchmark."""
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

BASE = Path("/Users/dho/Desktop/sandbox/mhc-primal-scheme/allele-coverage-development")
PARENT = BASE / "mhc-a1-a2-legacy-salvage-benchmark-pending-v1/legacy-off"
OUT = BASE / "mhc-a1-a2-gap-completion-benchmark-pending-v1"
EXE = Path("/Users/dho/Documents/lungfish-genome-explorer/.worktrees/allele-aware-primalscheme/.venv/bin/primalscheme3")
MSAS = [
    BASE / "mhc-fixture-snapshot-v3/snapshot/source-inputs/C64AA855-086A-4BBB-8AB7-803F0FAC0212/source.lungfishmsa/alignment/primary.aligned.fasta",
    BASE / "mhc-fixture-snapshot-v3/snapshot/source-inputs/2E748544-24F0-4CA0-9FE6-14FB8FB3CBB2/source.lungfishmsa/alignment/primary.aligned.fasta",
]
EXPECTED_INPUTS = [
    ("668a009c6e5adf706dca189100f0e0d5e3df63d68ced96865cfb73ad2f67e13e", 17940),
    ("2dd4690dc72e97442afc14f240b6abc101963363dab737a470fd6b968f9b9ea7", 28233),
]
COMMON = [
    "panel-create", "--msa", str(MSAS[0]), "--msa", str(MSAS[1]),
    "--selection-algorithm", "legacy", "--mode", "equal", "--amplicon-size", "200",
    "--amplicon-size-min", "150", "--amplicon-size-max", "250", "--n-pools", "2",
    "--ncores", "1", "--no-high-gc", "--min-base-freq", "0", "--mapping", "first",
    "--dimer-score", "-26", "--terminal-gap-policy", "legacy", "--use-matchdb",
    "--use-tm", "--offline-plots",
]
SCIENTIFIC_FILES = ("primer.bed", "amplicon.bed", "primertrim.amplicon.bed", "reference.fasta", "config.json")


def stamp(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "sizeBytes": path.stat().st_size}


def parent_payload() -> dict[str, dict[str, object]]:
    return {name: stamp(PARENT / name) for name in SCIENTIFIC_FILES}


def fasta_lengths(path: Path) -> dict[str, int]:
    result = {}; current = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]; result[current] = 0
        elif current is not None:
            result[current] += len(line.strip().replace("-", ""))
    return result


def bed_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"): continue
        fields = line.split("\t")
        result.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return result


def union(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]: merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else: merged.append((start, end))
    return merged


def geometry(panel: Path) -> dict[str, object]:
    lengths = fasta_lengths(panel / "reference.fasta")
    result = {"referenceLengths": lengths}
    for filename, key in (("primertrim.amplicon.bed", "trimmed"), ("amplicon.bed", "full")):
        by_target = bed_intervals(panel / filename); result[key] = {}
        for target, length in lengths.items():
            merged = union(by_target.get(target, [])); bases = sum(e - s for s, e in merged)
            gaps = []; cursor = 0
            for start, end in merged:
                if cursor < start: gaps.append([cursor, start])
                cursor = max(cursor, end)
            if cursor < length: gaps.append([cursor, length])
            result[key][target] = {"intervals": [list(x) for x in merged], "unionBases": bases, "referenceLength": length, "percentReferenceCovered": 100 * bases / length if length else 0.0, "gaps": gaps}
    return result


def commands() -> list[str]:
    return [str(EXE), *COMMON, "--output", str(OUT / "gap-completed"), "--gap-completion-parent", str(PARENT)]


def prepare() -> Path:
    if OUT.exists(): raise RuntimeError(f"refusing to overwrite {OUT}")
    for path, (digest, size) in zip(MSAS, EXPECTED_INPUTS, strict=True):
        item = stamp(path)
        if item["sha256"] != digest or item["sizeBytes"] != size: raise RuntimeError(f"input receipt mismatch: {path}")
    if not PARENT.is_dir(): raise RuntimeError(f"missing preserved parent: {PARENT}")
    OUT.mkdir(parents=True)
    manifest = {
        "schemaVersion": "lge.mhc-a1-a2-gap-completion-benchmark/v1",
        "status": "prepared; native launch requires explicit --run",
        "parent": str(PARENT.resolve()), "parentScientificPayloadBefore": parent_payload(),
        "inputs": [stamp(path) for path in MSAS], "argv": commands(), "command": shlex.join(commands()),
        "limits": {"wallTimeSeconds": 600, "rssBytes": 8 * 1024**3, "ncores": 1},
        "parentUpperBounds": {"A1": {"unionBases": 2172, "referenceLength": 2923, "maxGain": 135}, "A2": {"unionBases": 2291, "referenceLength": 3082, "maxGain": 113}},
        "preparedAt": datetime.now(timezone.utc).isoformat(),
    }
    path = OUT / "preparation-manifest.json"; path.write_text(json.dumps(manifest, indent=2) + "\n"); return path


def run() -> dict[str, object]:
    argv = commands(); output = OUT / "gap-completed"; started = time.monotonic(); started_at = datetime.now(timezone.utc).isoformat()
    stdout_path, stderr_path = OUT / "gap-completed.stdout.log", OUT / "gap-completed.stderr.log"
    stdout_file, stderr_file = stdout_path.open("w"), stderr_path.open("w")
    proc = subprocess.Popen(argv, cwd=Path(__file__).resolve().parents[1], stdout=stdout_file, stderr=stderr_file, start_new_session=True)
    peak_rss = 0; deadline = started + 600
    try:
        import psutil
    except ImportError:
        psutil = None
    while proc.poll() is None:
        if psutil is not None:
            try:
                process = psutil.Process(proc.pid)
                peak_rss = max(peak_rss, sum([process.memory_info().rss, *(child.memory_info().rss for child in process.children(recursive=True))]))
            except psutil.Error: pass
        else:
            try: peak_rss = max(peak_rss, int(subprocess.check_output(["ps", "-o", "rss=", "-p", str(proc.pid)], text=True).strip() or 0) * 1024)
            except (ValueError, subprocess.SubprocessError): pass
        if time.monotonic() >= deadline or peak_rss > 8 * 1024**3:
            os.killpg(proc.pid, 9); break
        time.sleep(0.25)
    proc.wait(); stdout_file.close(); stderr_file.close(); elapsed = time.monotonic() - started
    stdout, stderr = stdout_path.read_text(), stderr_path.read_text()
    receipt = {"argv": argv, "command": shlex.join(argv), "startedAt": started_at, "wallTimeSeconds": elapsed, "exitStatus": proc.returncode, "rssBytes": peak_rss, "stderr": stderr, "runtime": {"python": sys.version, "platform": platform.platform()}, "inputs": [stamp(path) for path in MSAS], "parentScientificPayloadBefore": parent_payload()}
    if output.is_dir():
        receipt["geometry"] = geometry(output); receipt["outputs"] = [stamp(path) for path in sorted(output.rglob("*")) if path.is_file()]
        audits = [path for path in output.glob("*gap*") if path.suffix == ".json"]
        receipt["nativeAudits"] = {path.name: json.loads(path.read_text()) for path in audits}
    receipt["parentScientificPayloadAfter"] = parent_payload(); receipt["parentUnchanged"] = receipt["parentScientificPayloadBefore"] == receipt["parentScientificPayloadAfter"]
    path = OUT / "gap-completed.benchmark-receipt.json"; path.write_text(json.dumps(receipt, indent=2) + "\n"); return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--run", action="store_true"); args = parser.parse_args()
    manifest = prepare(); print(json.dumps({"prepared": str(manifest), "command": shlex.join(commands())}, indent=2))
    if args.run: print(json.dumps(run(), indent=2))
