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


def _saved_control(root: Path, name: str) -> None:
    control = root / name
    result = control / "native" / "result.bed"
    result.parent.mkdir(parents=True)
    result.write_text(f"{name}\n")
    (control / "manifest.json").write_text(
        json.dumps(
            {
                "artifacts": [
                    {
                        "relativePath": "native/result.bed",
                        "byteSize": result.stat().st_size,
                        "sha256": _sha(result),
                    }
                ]
            }
        )
    )


def _matrix(root: Path, *, fail: bool = False, executable: str | None = None) -> Path:
    marker = root / "executed"
    code = (
        "import pathlib,sys; pathlib.Path(sys.argv[1]).write_text('yes'); "
        + ("sys.stderr.write('controlled failure\\n');sys.exit(7)" if fail else "sys.exit(0)")
    )
    probe = (
        "import json; print(json.dumps({"
        "'tool': 'fixture-tool', 'toolVersion':'1.2.3',"
        "'source': {'gitCommit':'abc123','gitDirty':False},"
        "'runtime': {'pythonVersion':'measured-runtime'}}))"
    )
    for name in ("historical-independent", "historical-combined"):
        _saved_control(root, name)
    matrix = {
        "schemaVersion": "allele-coverage-benchmark-matrix/v1",
        "fixtureAnalysis": "analysis.lungfishprimeranalysis",
        "inputLabels": ["duplicate.lungfishmsa"],
        "comparisons": [
            {
                "id": "historical-independent-lge.2",
                "kind": "saved-artifact",
                "root": "historical-independent",
                "manifest": "manifest.json",
            },
            {
                "id": "historical-combined-lge.2",
                "kind": "saved-artifact",
                "root": "historical-combined",
                "manifest": "manifest.json",
            },
            {"id": "upstream-original-3.3.0", "kind": "pending-execution"},
            {
                "id": "existing-lge.3",
                "kind": "native-execution",
                **({"executable": executable} if executable else {}),
                "identityProbeArgv": ["-c", probe],
                "expectedToolIdentity": {
                    "name": "fixture-tool",
                    "version": "1.2.3",
                    "gitCommit": "abc123",
                },
                "argv": ["-c", code, str(marker)],
                "optionsFullyResolved": True,
                "resolvedOptions": {
                    "selectionAlgorithm": "coverage",
                    "poolCount": 2,
                },
            },
        ],
    }
    path = root / "matrix.json"
    path.write_text(json.dumps(matrix))
    return path


def _run(
    root: Path,
    output: Path,
    *,
    corrupt: bool = False,
    fail: bool = False,
    executable: str | None = None,
):
    _fixture(root, corrupt=corrupt)
    matrix = _matrix(root, fail=fail, executable=executable)
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
    assert receipt["inputArtifacts"][0]["path"].endswith("matrix.json")
    mismatch = receipt["verificationFailure"]
    assert mismatch["expectedSha256"] == "0" * 64
    assert mismatch["observedSha256"] != mismatch["expectedSha256"]
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
    assert receipt["resolvedOptions"] == {
        "selectionAlgorithm": "coverage",
        "poolCount": 2,
    }
    stderr = output / receipt["stderrPath"]
    assert stderr.read_text() == "controlled failure\n"
    consumed = [item for item in receipt["inputArtifacts"] if item.get("role") == "scientific-input"]
    assert len(consumed) == 2
    assert all(item["path"].startswith("snapshot/") for item in consumed)
    assert receipt["toolIdentity"] == {"name": "fixture-tool", "version": "1.2.3"}
    assert receipt["sourceIdentity"] == {"gitCommit": "abc123", "gitDirty": False}
    assert receipt["runtimeIdentity"] == {"pythonVersion": "measured-runtime"}
    assert receipt["harnessIdentity"]["source"]


def test_saved_controls_are_verified_and_frozen_while_unrun_control_is_pending(tmp_path):
    output = tmp_path / "result"

    completed = _run(tmp_path / "fixtures", output)

    assert completed.returncode == 0, completed.stderr
    control = json.loads(
        (output / "controls" / "historical-independent-lge.2" / "receipt.json").read_text()
    )
    frozen = output / control["artifacts"][0]["path"]
    assert frozen.read_text() == "historical-independent\n"
    assert _sha(frozen) == control["artifacts"][0]["sha256"]
    pending = json.loads((output / "controls" / "upstream-original-3.3.0" / "receipt.json").read_text())
    assert pending["status"] == "pending-execution"
    assert pending["artifacts"] == []


def test_launch_error_still_writes_expanded_per_run_receipt(tmp_path):
    root = tmp_path / "fixtures"
    missing = root / "missing-primalscheme3"
    output = tmp_path / "result"

    completed = _run(root, output, executable=str(missing))

    assert completed.returncode != 0
    receipt = json.loads((output / "runs" / "existing-lge.3" / "receipt.json").read_text())
    assert receipt["argv"][0] == str(missing)
    assert receipt["exitStatus"] is None
    assert receipt["status"] == "launch-error"
    assert "No such file" in (output / receipt["stderrPath"]).read_text()


def test_msa_placeholder_expands_to_repeated_native_options(tmp_path):
    root = tmp_path / "fixtures"
    _fixture(root)
    matrix_path = _matrix(root)
    matrix = json.loads(matrix_path.read_text())
    native = next(item for item in matrix["comparisons"] if item["kind"] == "native-execution")
    native["argv"] = ["panel-create", "{msa_args}", "--output", "{output}"]
    matrix_path.write_text(json.dumps(matrix))
    output = tmp_path / "result"

    completed = subprocess.run(
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
            str(matrix_path),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    receipt = json.loads((output / "runs" / "existing-lge.3" / "receipt.json").read_text())
    pairs = list(zip(receipt["argv"], receipt["argv"][1:], strict=False))
    msa_paths = [value for option, value in pairs if option == "--msa"]
    assert len(msa_paths) == 2
    assert all("/snapshot/" in path for path in msa_paths)
    assert completed.returncode != 0  # Python rejects the fake panel-create command.
