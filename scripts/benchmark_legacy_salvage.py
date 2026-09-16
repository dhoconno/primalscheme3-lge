#!/usr/bin/env python3
"""Prepare (and, only with ``--run``, execute) the bounded legacy salvage gate.

The default action writes a manifest and performs no native design.  It is
intended for review before launching the two A1+A2 runs.
"""
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
BASELINE = BASE / "mhc-a1-a2-legacy-two-pool-01"
OUT = BASE / "mhc-a1-a2-legacy-salvage-benchmark-pending-v1"
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


def stamp(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"path": str(path.resolve()), "sha256": digest, "sizeBytes": path.stat().st_size}


def canonical_bed(path: Path) -> list[tuple[str, int, int, str, str, str]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        # Keep pool, strand and sequence; generated amplicon names are ignored.
        rows.append((fields[0], int(fields[1]), int(fields[2]), fields[4], fields[5], fields[6] if len(fields) > 6 else ""))
    return sorted(rows)


def strict_fingerprint(panel: Path) -> str:
    """Fingerprint coordinate/sequence rows while ignoring generated names."""
    rows = canonical_bed(panel / ("strict-primer.bed" if (panel / "strict-primer.bed").exists() else "primer.bed"))
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def geometry(panel: Path) -> dict[str, object]:
    """Use the established evaluator's interval convention on each output."""
    refs = {}
    current = None
    for line in (panel / "reference.fasta").read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            refs[current] = ""
        elif current is not None:
            refs[current] += line.strip().replace("-", "")
    result = {"referenceLengths": {name: len(seq) for name, seq in refs.items()}}
    for rel, key in (("amplicon.bed", "fullSpan"), ("primertrim.amplicon.bed", "primerTrimmedInterior")):
        rows = [x.split("\t") for x in (panel / rel).read_text().splitlines() if x and not x.startswith("#")]
        by_chrom = {}
        for row in rows:
            by_chrom.setdefault(row[0], []).append((int(row[1]), int(row[2])))
        result[key] = {}
        for chrom, intervals in sorted(by_chrom.items()):
            if chrom not in refs:
                raise ValueError(f"{rel}: output target {chrom!r} missing from reference.fasta")
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            bases = sum(e - s for s, e in merged)
            result[key][chrom] = {"intervals": [list(x) for x in merged], "unionBases": bases, "referenceLength": len(refs[chrom]), "percentReferenceCovered": 100 * bases / len(refs[chrom])}
    return result


def commands() -> dict[str, list[str]]:
    return {
        "legacy-off": [str(EXE), *COMMON, "--output", str(OUT / "legacy-off"), "--legacy-salvage", "off"],
        "legacy-bounded": [str(EXE), *COMMON, "--output", str(OUT / "legacy-bounded"), "--legacy-salvage", "bounded"],
    }


def prepare() -> Path:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite {OUT}")
    for path, (digest, size) in zip(MSAS, EXPECTED_INPUTS, strict=True):
        actual = stamp(path)
        if actual["sha256"] != digest or actual["sizeBytes"] != size:
            raise RuntimeError(f"input receipt mismatch: {path}")
    baseline_config = json.loads((BASELINE / "panel/config.json").read_text())
    expected = {
        "selection_algorithm": "legacy", "mode": "equal", "mapping": "first",
        "amplicon_size": 200, "amplicon_size_min": 150, "amplicon_size_max": 250,
        "n_pools": 2, "dimer_score": -26.0, "terminal_gap_policy": "legacy",
        "high_gc": False, "min_base_freq": 0.0, "use_matchdb": True,
        "mismatch_product_size": 0, "circular": False,
    }
    drift = {key: (baseline_config.get(key), value) for key, value in expected.items() if baseline_config.get(key) != value}
    if drift:
        raise RuntimeError(f"baseline config drift: {drift}")
    OUT.mkdir(parents=True)
    manifest = {
        "schemaVersion": "lge.mhc-a1-a2-legacy-salvage-benchmark/v1",
        "status": "prepared; native launch requires explicit --run",
        "baselineReceipt": str((BASELINE / "harness-receipt.json").resolve()),
        "baselinePanel": str((BASELINE / "panel").resolve()),
        "inputs": [stamp(x) for x in MSAS],
        "commands": {k: {"argv": v, "command": shlex.join(v), "output": v[v.index("--output") + 1]} for k, v in commands().items()},
        "limits": {"wallTimeSeconds": 600, "rssBytes": 8 * 1024**3, "ncores": 1},
        "resolvedSalvage": {"mode": "bounded", "thresholds": [-28, -30, -32], "floor": -32, "maxEdgesPerPool": 8, "maxIncidentSpeciesPerPool": 4, "minReferenceGain": 1, "maxCandidateEvaluations": 10000},
        "comparison": {"strictFingerprint": "canonical primer BED coordinates/sequences; generated names ignored", "geometry": "zero-based half-open ungapped first-row reference; per-MSA interval unions", "risk": "legacy-salvage.json stage accepted conflicts and cumulative counts"},
        "preparedAt": datetime.now(timezone.utc).isoformat(),
    }
    path = OUT / "preparation-manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def run_one(name: str, argv: list[str]) -> dict[str, object]:
    output = Path(argv[argv.index("--output") + 1])
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    proc = subprocess.Popen(argv, cwd=Path(__file__).resolve().parents[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    peak_rss = 0
    deadline = started + 600
    try:
        import psutil
    except ImportError:
        psutil = None
    while proc.poll() is None:
        if psutil is not None:
            try:
                process = psutil.Process(proc.pid)
                peak_rss = max(peak_rss, process.memory_info().rss)
                peak_rss = max(peak_rss, *(child.memory_info().rss for child in process.children(recursive=True)))
            except psutil.Error:
                pass
        if time.monotonic() >= deadline or peak_rss > 8 * 1024**3:
            os.killpg(proc.pid, 9)
            break
        time.sleep(0.25)
    stdout, stderr = proc.communicate()
    if proc.returncode is None:
        proc.wait()
    elapsed = time.monotonic() - started
    output.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.stdout.log").write_text(stdout)
    (OUT / f"{name}.stderr.log").write_text(stderr)
    receipt = {"name": name, "argv": argv, "command": shlex.join(argv), "startedAt": started_at, "wallTimeSeconds": elapsed, "exitStatus": proc.returncode, "stderr": stderr, "rssBytes": peak_rss, "runtime": {"python": sys.version, "platform": platform.platform()}, "inputs": [stamp(x) for x in MSAS], "outputs": [stamp(x) for x in sorted(output.rglob("*")) if x.is_file()]}
    if (output / "primer.bed").exists():
        receipt["strictFingerprintOrFinalFingerprint"] = strict_fingerprint(output)
        receipt["geometry"] = geometry(output)
    if (output / "legacy-salvage.json").exists():
        receipt["salvage"] = json.loads((output / "legacy-salvage.json").read_text())
    (OUT / f"{name}.benchmark-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="launch both bounded native runs after preparation")
    args = parser.parse_args()
    manifest = prepare()
    print(json.dumps({"prepared": str(manifest), "commands": {k: shlex.join(v) for k, v in commands().items()}}, indent=2))
    if args.run:
        for name, argv in commands().items():
            print(json.dumps(run_one(name, argv), indent=2), flush=True)
