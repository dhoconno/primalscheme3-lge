"""Contract tests for the immutable allele-coverage benchmark runner."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "benchmark_allele_coverage.py"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(root: Path, *, corrupt: bool = False) -> Path:
    analysis = root / "analysis.lungfishprimeranalysis"
    artifacts = []
    inputs = []
    for occurrence, label in (("occurrence-a", "duplicate.lungfishmsa"), ("occurrence-b", "duplicate.lungfishmsa")):
        relative = f"source-inputs/{occurrence}/source.lungfishmsa/alignment/primary.aligned.fasta"
        path = analysis / relative
        path.parent.mkdir(parents=True)
        path.write_text(f">same-label\nACGT{occurrence[-1]}\n")
        artifacts.append(
            {
                "relativePath": relative,
                "role": "input",
                "format": "fasta",
                "byteSize": path.stat().st_size,
                "sha256": "0" * 64 if corrupt and occurrence == "occurrence-a" else _sha(path),
            }
        )
        inputs.append({"id": occurrence, "label": label, "artifactPaths": [relative]})
    manifest = {"schemaVersion": 1, "grouping": "combined", "inputs": inputs, "artifacts": artifacts}
    (analysis / "manifest.json").write_text(json.dumps(manifest))
    return analysis


def _matrix(root: Path, *, fail: bool = False) -> Path:
    marker = root / "executed"
    code = (
        "import pathlib,sys; pathlib.Path(sys.argv[1]).write_text('yes'); "
        + ("sys.stderr.write('controlled failure\\n');sys.exit(7)" if fail else "sys.exit(0)")
    )
    matrix = {
        "schemaVersion": "allele-coverage-benchmark-matrix/v1",
        "fixtureAnalysis": "analysis.lungfishprimeranalysis",
        "inputLabels": ["duplicate.lungfishmsa"],
        "comparisons": [
            {"id": "historical-independent-lge.2", "kind": "historical-artifact"},
            {"id": "historical-combined-lge.2", "kind": "historical-artifact"},
            {
                "id": "existing-lge.3",
                "kind": "native-execution",
                "argv": ["-c", code, str(marker)],
                "resolvedOptions": {"selectionAlgorithm": "coverage"},
            },
        ],
    }
    path = root / "matrix.json"
    path.write_text(json.dumps(matrix))
    return path


def _run(root: Path, output: Path, *, corrupt: bool = False, fail: bool = False):
    _fixture(root, corrupt=corrupt)
    matrix = _matrix(root, fail=fail)
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--fixture-root",
            str(root),
            "--native-executable",
            sys.executable,
            "--output",
            str(output),
            "--matrix",
            str(matrix),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_existing_output_is_refused_without_mutation(tmp_path):
    output = tmp_path / "result"
    output.mkdir()
    sentinel = output / "keep"
    sentinel.write_text("unchanged")

    completed = _run(tmp_path / "fixtures", output)

    assert completed.returncode != 0
    assert "already exists" in completed.stderr
    assert sentinel.read_text() == "unchanged"
    assert list(output.iterdir()) == [sentinel]


def test_mismatched_manifest_checksum_aborts_before_scientific_execution(tmp_path):
    root = tmp_path / "fixtures"
    output = tmp_path / "result"

    completed = _run(root, output, corrupt=True)

    assert completed.returncode != 0
    assert not (root / "executed").exists()
    receipt = json.loads((output / "runner-receipt.json").read_text())
    assert receipt["exitStatus"] != 0
    assert "checksum mismatch" in (output / receipt["stderrPath"]).read_text()


def test_duplicate_labels_keep_distinct_source_occurrence_ids(tmp_path):
    root = tmp_path / "fixtures"
    output = tmp_path / "result"

    completed = _run(root, output)

    assert completed.returncode == 0, completed.stderr
    labels = json.loads((output / "snapshot" / "source-row-label-map.json").read_text())
    assert [item["sourceOccurrenceID"] for item in labels["inputs"]] == [
        "occurrence-a",
        "occurrence-b",
    ]
    assert [item["label"] for item in labels["inputs"]] == [
        "duplicate.lungfishmsa",
        "duplicate.lungfishmsa",
    ]


def test_failed_subprocess_records_status_and_stderr(tmp_path):
    output = tmp_path / "result"

    completed = _run(tmp_path / "fixtures", output, fail=True)

    assert completed.returncode == 7
    receipt = json.loads((output / "runs" / "existing-lge.3" / "receipt.json").read_text())
    assert receipt["exitStatus"] == 7
    assert receipt["wallTimeSeconds"] >= 0
    assert receipt["argv"][0] == sys.executable
    assert receipt["resolvedOptions"] == {"selectionAlgorithm": "coverage"}
    stderr = output / receipt["stderrPath"]
    assert stderr.read_text() == "controlled failure\n"
    assert receipt["inputArtifacts"]
    assert receipt["sourceIdentity"]
    assert receipt["runtimeIdentity"]
