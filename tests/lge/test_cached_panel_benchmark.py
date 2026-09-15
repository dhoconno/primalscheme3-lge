"""Controlled cached matrix contract; expensive scientific panels are not fixtures."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest


def api():
    path = (
        Path(__file__).resolve().parents[2] / "scripts/benchmark_cached_allele_panel.py"
    )
    spec = importlib.util.spec_from_file_location("cached_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_arms_have_effective_ablation_controls_and_matched_seed_scope():
    a = api()
    arms = a.arm_manifest()
    assert arms["normal-full"]["options"]["variant_selection"] == "full-cloud"
    assert arms["union-full"]["options"]["candidate_profiles"] == "union"
    assert arms["union-subsets-repair0"]["options"]["optimizer_repair_rounds"] == 0
    assert arms["union-subsets-600"]["options"]["optimizer_time_limit"] == 600
    assert arms["union-subsets-salvage"]["options"]["primary_tier"] == "strict"
    assert (
        arms["normal-k17"]["options"] | {"specificity_terminal_k": 19}
        == arms["normal-k19"]["options"]
    )
    assert (
        "normal-only" in arms["normal-k19"]["scope"]
        and "19" in arms["normal-k19"]["scope"]
    )
    assert arms["fixed-work-repeat"]["repeat"] == 2
    assert arms["fixed-work-repeat"]["options"]["optimizer_time_limit"] == 3600


def test_completed_fixed_work_rejects_interruption_and_compares_actual_work():
    a = api()
    metadata = {
        "stop_reason": "completed",
        "completed_starts": 1,
        "completed_repair_rounds": 2,
        "effective_seed_modes": ["full", "normal"],
        "completed_seed_modes": ["full", "normal"],
        "options": {"starts": 1, "repair_rounds": 2},
    }
    assert a.completed_fixed_work(metadata)
    assert not a.completed_fixed_work(metadata | {"stop_reason": "time-limit"})
    assert not a.completed_fixed_work(metadata | {"completed_repair_rounds": 1})


def make_snapshot(root, fasta):
    a = api()
    (root / "snapshot/source-inputs/id/source.lungfishmsa/alignment").mkdir(
        parents=True
    )
    source = (
        root
        / "snapshot/source-inputs/id/source.lungfishmsa/alignment/primary.aligned.fasta"
    )
    source.write_bytes(fasta.read_bytes())
    normalized = root / "snapshot/inputs/id.fasta"
    normalized.parent.mkdir()
    normalized.write_text(">WRONG\nCCCC\n")
    labels = {
        "schemaVersion": "allele-coverage-source-label-map/v1",
        "inputs": [
            {
                "label": "target",
                "sourceOccurrenceID": "id",
                "snapshotArtifactPaths": [
                    str(source.relative_to(root)),
                    str(normalized.relative_to(root)),
                ],
            }
        ],
    }
    label = root / "snapshot/source-row-label-map.json"
    label.write_text(json.dumps(labels))
    manifest = {
        "artifacts": [
            dict(a.descriptor(p, root), byteSize=p.stat().st_size)
            for p in (source, label)
        ]
    }
    (root / "snapshot/hash-manifest.json").write_text(json.dumps(manifest))
    return source


def test_original_manifest_order_and_bytes_are_required(tmp_path):
    a = api()
    raw = tmp_path / "raw.fa"
    raw.write_text(">original\nAAAA\n")
    source = make_snapshot(tmp_path / "fixture", raw)
    records, _ = a.snapshot_inputs(tmp_path / "fixture")
    assert records[0]["path"] == str(source.resolve())
    source.write_text("changed")
    with pytest.raises(ValueError, match="snapshot"):
        a.snapshot_inputs(tmp_path / "fixture")


def test_subprocess_receipt_measures_one_real_child_and_retains_failure(tmp_path):
    result = api().run_child(
        [
            sys.executable,
            "-c",
            'import sys;print("hello");sys.stderr.write("diagnostic");sys.exit(3)',
        ],
        tmp_path / "child",
    )
    assert result["exitStatus"] == 3
    assert result["stderr"] == "diagnostic"
    assert (
        result["peakRSS"]["scope"].startswith("wait4")
        or result["peakRSS"]["value"] is None
    )
    assert (
        json.loads((tmp_path / "child/command.json").read_text())["argv"][0]
        == sys.executable
    )


def test_matrix_requires_prior_freeze_and_retains_failure_receipt(tmp_path):
    a = api()
    out = tmp_path / "failed"
    result = a.main(
        [
            "--snapshot",
            str(tmp_path / "missing"),
            "--cache",
            str(tmp_path / "cache"),
            "--native-executable",
            sys.executable,
            "--output",
            str(out),
            "--arm",
            "normal-full",
        ]
    )
    assert result == 1
    assert json.loads((out / "provenance.json").read_text())["exitStatus"] == 1


def test_real_tiny_cache_freeze_panel_and_audit(tmp_path):
    import subprocess

    root = Path(__file__).resolve().parents[2]
    native = root / ".venv/bin/primalscheme3"
    raw = tmp_path / "original.fa"
    raw.write_text(">original\n" + "A" * 30 + "\n")
    original = tmp_path / "original-panel"
    cache = tmp_path / "cache"
    completed = subprocess.run(
        [
            str(native),
            "panel-create",
            "--selection-algorithm",
            "allele-coverage",
            "--mode",
            "equal",
            "--msa",
            str(raw),
            "--output",
            str(original),
            "--amplicon-size",
            "200",
            "--amplicon-size-min",
            "150",
            "--amplicon-size-max",
            "250",
            "--ncores",
            "1",
            "--optimizer-starts",
            "1",
            "--optimizer-repair-rounds",
            "0",
        ],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    completed = subprocess.run(
        [str(native), "panel-cache", "--bundle", str(original), "--output", str(cache)],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    source = make_snapshot(tmp_path / "fixture", raw)
    a = api()
    common = [
        "--snapshot",
        str(tmp_path / "fixture"),
        "--cache",
        str(cache),
        "--native-executable",
        str(native),
    ]
    freeze = tmp_path / "freeze"
    assert a.main(common + ["--output", str(freeze), "--freeze-only"]) == 0
    output = tmp_path / "matrix"
    assert (
        a.main(
            common
            + [
                "--output",
                str(output),
                "--source-freeze",
                str(freeze / "freeze.json"),
                "--arm",
                "normal-full",
            ]
        )
        == 0
    )
    report = json.loads((output / "summary.json").read_text())
    arm = report["arms"][0]
    assert arm["status"] == "success"
    assert arm["native"]["argv"].count(str(source.resolve())) == 1
    assert "WRONG" not in json.dumps(report)
    assert arm["resolvedOptions"]["variant_selection"] == "full-cloud"
    assert arm["audit"]["exitStatus"] == 0
    assert arm["boundReceipts"]
    assert (output / "runs/normal-full-1/panel/config.json").is_file()


@pytest.mark.parametrize("change", ["source", "input"])
def test_freeze_blocks_changed_identity_before_any_arm(tmp_path, monkeypatch, change):
    a = api()
    raw = tmp_path / "original.fa"
    raw.write_text(">original\nAAAA\n")
    source = make_snapshot(tmp_path / "fixture", raw)
    cache = tmp_path / "cache"
    cache.mkdir()
    manifest = cache / "manifest.json"
    manifest.write_text("{}")
    monkeypatch.setattr(
        a, "cache_inputs", lambda path: ({}, {}, {}, [a.descriptor(manifest)], {})
    )
    capability = {
        "source": {"sourceDigest": "initial"},
        "runtime": {"python": "test"},
        "toolVersion": "test",
    }
    monkeypatch.setattr(a, "probe", lambda *args: capability)
    common = [
        "--snapshot",
        str(tmp_path / "fixture"),
        "--cache",
        str(cache),
        "--native-executable",
        sys.executable,
    ]
    freeze = tmp_path / "freeze"
    assert a.main(common + ["--output", str(freeze), "--freeze-only"]) == 0
    if change == "source":
        capability["source"]["sourceDigest"] = "changed"
    else:
        manifest.write_text('{"changed":true}')
    monkeypatch.setattr(
        a, "run_child", lambda *args: pytest.fail("must reject before arm launch")
    )
    out = tmp_path / "run"
    assert (
        a.main(
            common
            + [
                "--output",
                str(out),
                "--source-freeze",
                str(freeze / "freeze.json"),
                "--arm",
                "normal-full",
            ]
        )
        == 1
    )
    assert json.loads((out / "provenance.json").read_text())["exitStatus"] == 1
    assert source.read_bytes() == raw.read_bytes()


def test_child_launch_failure_retains_exit_and_error(tmp_path):
    result = api().run_child([str(tmp_path / "absent-executable")], tmp_path / "launch")
    assert result["exitStatus"] == 127 and "FileNotFoundError" in result["stderr"]
    assert (tmp_path / "launch/command.json").is_file()
