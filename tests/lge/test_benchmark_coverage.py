"""Behavior tests for the retained abstract coverage benchmark utility."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "benchmark_coverage.py"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(output),
            "--target-count",
            "4",
            "--seed",
            "7",
            "--starts",
            "1",
            "--repair-rounds",
            "1",
            "--time-limit",
            "20",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_small_abstract_benchmark_uses_real_search_and_verifies_final_bytes(tmp_path):
    output = tmp_path / "abstract-benchmark"

    completed = _run(output)

    assert completed.returncode == 0, completed.stderr
    fixture = json.loads((output / "fixture.json").read_text())
    result = json.loads((output / "result.json").read_text())
    provenance = json.loads((output / "provenance.json").read_text())
    assert fixture["scope"].startswith("ABSTRACT ORACLE ONLY")
    assert fixture["targetCount"] == 4
    assert len(fixture["candidates"]) == 12
    assert len(fixture["conflicts"]) == 8
    assert fixture["runs"]["repaired"] == {
        "metric": "full-span",
        "coverageTarget": 0.9,
        "seed": 7,
        "starts": 1,
        "repairRounds": 1,
        "timeLimitSeconds": 20.0,
    }
    greedy, repaired = result["runs"]
    assert greedy["independentValidation"]["valid"] is True
    assert repaired["independentValidation"]["valid"] is True
    assert greedy["independentValidation"]["totalFullCoveredBases"] == 1_040
    assert repaired["independentValidation"]["totalFullCoveredBases"] == 1_200
    assert repaired["independentValidation"]["totalInteriorCoveredBases"] == 1_184
    assert result["comparison"]["repairedReachedExpectedFeasibleAlternative"] is True
    assert repaired["searchMetadata"]["validation"]["scope"] == "abstract-oracle"
    assert repaired["searchMetadata"]["stop_reason"] == "completed"
    assert repaired["searchMetadata"]["completed_starts"] == 1
    assert repaired["searchMetadata"]["repairs_accepted"] == 4
    assert result["defaultWorkLimits"] == repaired["searchMetadata"]["search_limits"]
    assert provenance["exitStatus"] == 0
    assert provenance["status"] == "success"
    assert provenance["scientificScope"].startswith("ABSTRACT ORACLE ONLY")
    assert provenance["inputsChangedDuringRun"] is False
    assert provenance["sourceChangedDuringRun"] is False
    assert provenance["runtimeChangedDuringRun"] is False
    assert provenance["runtime"]["nativeKernels"]
    assert any(
        file["path"].endswith(".so")
        for kernel in provenance["runtime"]["nativeKernels"]
        for file in kernel["files"]
    )
    expected = {record["path"]: record for record in provenance["outputs"]}
    actual = {
        str(path.relative_to(output)): path
        for path in output.rglob("*")
        if path.is_file() and path.name != "provenance.json"
    }
    assert set(expected) == set(actual)
    for relative, path in actual.items():
        assert expected[relative]["size"] == path.stat().st_size
        assert expected[relative]["sha256"] == _sha256(path)


def test_abstract_benchmark_refuses_existing_output_without_mutating_it(tmp_path):
    output = tmp_path / "abstract-benchmark"
    first = _run(output)
    assert first.returncode == 0, first.stderr
    before = {
        str(path.relative_to(output)): (path.stat().st_size, _sha256(path))
        for path in output.rglob("*")
        if path.is_file()
    }

    second = _run(output)

    assert second.returncode != 0
    assert "already exists" in second.stderr
    after = {
        str(path.relative_to(output)): (path.stat().st_size, _sha256(path))
        for path in output.rglob("*")
        if path.is_file()
    }
    assert after == before
