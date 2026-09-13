"""Public CLI contract for the opt-in coverage panel selector."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from click import UsageError
from typer.testing import CliRunner

from primalscheme3.cli import app
from primalscheme3.core.config import Config
from primalscheme3.panel.coverage_validation import ConstraintProfile
from primalscheme3.panel.panel_classes import PanelRunModes
from primalscheme3.panel.panel_main import panelcreate


def _input(tmp_path: Path) -> Path:
    path = tmp_path / "input.fasta"
    path.write_text(">reference\nACGTACGT\n")
    return path


def test_capabilities_reports_actual_local_engine_and_runtime_identity():
    result = CliRunner().invoke(app, ["--capabilities-json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["schemaVersion"] == "primalscheme3.capabilities/v1"
    assert payload["toolVersion"] == "3.3.0+lge.3"
    assert payload["selectionAlgorithms"] == ["legacy", "coverage"]
    assert payload["coverage"]["profile"] == {
        "name": "panel-v1",
        "specificityRevision": "intended-sites-v1",
    }
    assert (
        payload["coverage"]["optimizerSchemaVersion"]
        == "primalscheme3.panel-optimizer/v1"
    )
    assert (
        payload["coverage"]["validationSchemaVersion"]
        == "primalscheme3.panel-validation/v1"
    )
    assert payload["source"]["sourceDigest"]
    assert payload["source"]["files"]
    assert payload["runtime"]["pythonVersion"]
    assert payload["runtime"]["platform"]
    primal = next(
        item
        for item in payload["runtime"]["nativeKernels"]
        if item["package"] == "primalschemers"
    )
    assert any(
        item["path"].endswith((".so", ".dylib", ".pyd")) for item in primal["files"]
    )
    assert isinstance(payload["source"]["gitDirty"], bool)


def test_legacy_cli_defaults_are_unchanged(tmp_path):
    msa = _input(tmp_path)
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app,
            ["panel-create", "--msa", str(msa), "--output", str(tmp_path / "out")],
        )
    assert result.exit_code == 0, result.output
    config = run.call_args.kwargs["config"]
    assert config.selection_algorithm == "legacy"
    assert config.mismatch_product_size == 0
    assert config.amplicon_size_metric.value == "legacy-pairing"


def test_legacy_cli_preserves_positive_count_limit_contract(tmp_path):
    msa = _input(tmp_path)
    result = CliRunner().invoke(
        app,
        [
            "panel-create",
            "--msa",
            str(msa),
            "--output",
            str(tmp_path / "out"),
            "--max-amplicons",
            "0",
        ],
    )
    assert result.exit_code == 2
    assert not (tmp_path / "out").exists()


def test_legacy_explicit_zero_cli_and_programmatic_product_size_contract(tmp_path):
    msa = _input(tmp_path)
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app,
            [
                "panel-create",
                "--msa",
                str(msa),
                "--output",
                str(tmp_path / "out"),
                "--mispriming-product-size",
                "0",
            ],
        )
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["config"].mismatch_product_size == 0
    assert Config(mismatch_product_size=321).mismatch_product_size == 321


def test_coverage_cli_resolves_contract_and_dynamic_product_default(tmp_path):
    msa = _input(tmp_path)
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app,
            [
                "panel-create",
                "--msa",
                str(msa),
                "--output",
                str(tmp_path / "out"),
                "--mode",
                "equal",
                "--selection-algorithm",
                "coverage",
                "--terminal-gap-policy",
                "observed-only",
                "--amplicon-size",
                "200",
                "--amplicon-size-min",
                "150",
                "--amplicon-size-max",
                "280",
                "--coverage-metric",
                "primer-trimmed",
                "--coverage-target",
                "0.8",
                "--optimizer-seed",
                "7",
                "--optimizer-starts",
                "3",
                "--optimizer-repair-rounds",
                "1",
                "--optimizer-time-limit",
                "12.5",
            ],
        )
    assert result.exit_code == 0, result.output
    config = run.call_args.kwargs["config"]
    assert config.selection_algorithm == "coverage"
    assert config.mismatch_product_size == 2000
    assert config.coverage_metric == "primer-trimmed"
    assert config.coverage_target == 0.8
    assert config.optimizer_seed == 7
    assert config.optimizer_starts == 3
    assert config.optimizer_repair_rounds == 1
    assert config.optimizer_time_limit == 12.5


def test_coverage_accepts_both_discovery_backends_and_zero_caps(tmp_path):
    msa = _input(tmp_path)
    for policy in ("legacy", "observed-only"):
        with patch("primalscheme3.cli.panelcreate") as run:
            result = CliRunner().invoke(
                app,
                [
                    "panel-create",
                    "--msa",
                    str(msa),
                    "--output",
                    str(tmp_path / policy),
                    "--mode",
                    "equal",
                    "--selection-algorithm",
                    "coverage",
                    "--terminal-gap-policy",
                    policy,
                    "--amplicon-size",
                    "200",
                    "--amplicon-size-min",
                    "150",
                    "--amplicon-size-max",
                    "280",
                    "--max-amplicons",
                    "0",
                    "--max-amplicons-msa",
                    "0",
                ],
            )
        assert result.exit_code == 0, result.output
        assert run.call_args.kwargs["config"].terminal_gap_policy.value == policy
        assert run.call_args.kwargs["max_amplicons"] == 0
        assert run.call_args.kwargs["max_amplicons_msa"] == 0


@pytest.mark.parametrize(
    "field",
    [
        "coverage_target",
        "optimizer_seed",
        "optimizer_starts",
        "optimizer_repair_rounds",
        "optimizer_time_limit",
    ],
)
def test_programmatic_optimizer_numeric_options_reject_booleans(field):
    with pytest.raises(ValueError):
        Config(**{field: True})


@pytest.mark.parametrize(
    "flags",
    [
        ["--selection-algorithm", "coverage", "--mode", "entropy"],
        ["--selection-algorithm", "coverage", "--mode", "region-only"],
        ["--selection-algorithm", "coverage", "--mode", "equal", "--no-use-matchdb"],
        [
            "--selection-algorithm",
            "coverage",
            "--mode",
            "equal",
            "--terminal-gap-policy",
            "observed-only",
        ],
        ["--selection-algorithm", "legacy", "--coverage-target", "0.8"],
        ["--selection-algorithm", "legacy", "--mispriming-product-size", "1"],
        ["--selection-algorithm", "coverage", "--coverage-target", "nan"],
        ["--selection-algorithm", "coverage", "--optimizer-time-limit", "inf"],
        ["--selection-algorithm", "coverage", "--optimizer-starts", "0"],
        ["--selection-algorithm", "coverage", "--optimizer-repair-rounds", "-1"],
    ],
)
def test_invalid_selector_cli_fails_before_workflow_or_output(tmp_path, flags):
    msa = _input(tmp_path)
    output = tmp_path / "out"
    with patch(
        "primalscheme3.cli.panelcreate", side_effect=AssertionError("workflow called")
    ):
        result = CliRunner().invoke(
            app,
            ["panel-create", "--msa", str(msa), "--output", str(output), *flags],
        )
    assert result.exit_code == 2, result.output + repr(result.exception)
    assert not output.exists()


def test_direct_coverage_rejects_scope_before_output_creation(tmp_path):
    output = tmp_path / "out"
    config = Config(
        selection_algorithm="coverage",
        terminal_gap_policy="observed-only",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
    )
    with pytest.raises(UsageError, match="equal"):
        panelcreate([], output, config, None, mode=PanelRunModes.ENTROPY)
    assert not output.exists()


@pytest.mark.parametrize("cap", [-1, True])
def test_direct_coverage_rejects_invalid_caps_before_output_creation(tmp_path, cap):
    output = tmp_path / "out"
    config = Config(
        selection_algorithm="coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
    )
    with pytest.raises(UsageError, match="nonnegative integer"):
        panelcreate(
            [], output, config, None, mode=PanelRunModes.EQUAL, max_amplicons=cap
        )
    assert not output.exists()


def test_coverage_discovery_failure_retains_failure_provenance(tmp_path):
    msa = _input(tmp_path)
    output = tmp_path / "out"
    config = Config(
        selection_algorithm="coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
    )
    with patch(
        "primalscheme3.panel.panel_main.PanelMSA",
        side_effect=RuntimeError("synthetic discovery failure"),
    ):
        with pytest.raises(RuntimeError, match="synthetic discovery failure"):
            panelcreate(
                [msa],
                output,
                config,
                None,
                mode=PanelRunModes.EQUAL,
                executed_argv=["primalscheme3", "panel-create"],
            )
    record = json.loads((output / "panel-provenance.json").read_text())
    assert record["status"] == "failure"
    assert record["exitStatus"] == 1
    assert record["stderr"] == "synthetic discovery failure"
    assert record["inputs"][0]["storedMissing"] is False
    assert record["inputs"][0]["sourceAtStartMatchesStored"] is True


def _tree_bytes(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file()
    }


def test_forced_discovery_failure_replaces_stale_success_with_current_final_bytes(
    tmp_path,
):
    source = _input(tmp_path)
    output = tmp_path / "out"
    (output / "work").mkdir(parents=True)
    (output / "old-result.txt").write_text("previous bytes")
    (output / "panel-provenance.json").write_text(
        '{"status":"success","exitStatus":0,"marker":"previous run"}'
    )
    config = Config(
        selection_algorithm="coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
    )
    argv = ["primalscheme3", "panel-create", "--force", "--msa", str(source)]
    with patch(
        "primalscheme3.panel.panel_main.PanelMSA",
        side_effect=RuntimeError("synthetic new-run discovery failure"),
    ):
        with pytest.raises(RuntimeError, match="synthetic new-run discovery failure"):
            panelcreate(
                [source],
                output,
                config,
                None,
                mode=PanelRunModes.EQUAL,
                force=True,
                executed_argv=argv,
            )
    record = json.loads((output / "panel-provenance.json").read_text())
    assert record["status"] == "failure"
    assert record["exitStatus"] == 1
    assert record["command"]["argv"] == argv
    assert record["inputs"][0]["sourceAtStartMatchesStored"] is True
    assert record["inputs"][0]["sourcePath"] == str(source.resolve())
    for descriptor in record["outputs"]:
        data = (output / descriptor["path"]).read_bytes()
        assert len(data) == descriptor["size"]
        assert hashlib.sha256(data).hexdigest() == descriptor["sha256"]
    assert {item["path"] for item in record["outputs"]} == {
        "old-result.txt",
        "work/0000-input.fasta",
        "work/file.log",
    }


def test_nonforce_and_preflight_rejections_preserve_existing_output_bytes(tmp_path):
    source = _input(tmp_path)
    config = Config(
        selection_algorithm="coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
    )
    for label, mode, force in (
        ("nonforce", PanelRunModes.EQUAL, False),
        ("preflight", PanelRunModes.ENTROPY, True),
    ):
        output = tmp_path / label
        (output / "work").mkdir(parents=True)
        (output / "work" / "existing.bin").write_bytes(b"unchanged\x00bytes")
        (output / "panel-provenance.json").write_text('{"marker":"existing"}')
        before = _tree_bytes(output)
        with pytest.raises(UsageError):
            panelcreate([source], output, config, None, mode=mode, force=force)
        assert _tree_bytes(output) == before


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("optimizer_starts", 1.7),
        ("optimizer_seed", 1.7),
        ("optimizer_repair_rounds", 1.7),
    ],
)
def test_optimizer_integers_reject_fractional_values_before_coercion(field, value):
    with pytest.raises(ValueError, match="integer"):
        Config(**{field: value})


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("n_pools", True),
        ("mismatch_product_size", True),
        ("amplicon_size", 200.5),
        ("amplicon_size_min", True),
        ("amplicon_size_max", 280.5),
        ("ncores", True),
        ("dimer_score", True),
        ("dimer_score", float("nan")),
        ("min_base_freq", True),
        ("min_base_freq", float("nan")),
        ("min_base_freq", 1.1),
    ],
)
def test_coverage_scientific_numerics_reject_invalid_raw_values(field, value):
    options = {
        "selection_algorithm": "coverage",
        "amplicon_size": 200,
        "amplicon_size_min": 150,
        "amplicon_size_max": 280,
        "amplicon_size_metric": "reference-span",
        "mismatch_product_size": 2000,
    }
    options[field] = value
    with pytest.raises(ValueError):
        Config(**options)


def test_legacy_scientific_config_coercion_and_positive_product_size_remain_compatible():
    config = Config(n_pools=1.7, mismatch_product_size=True)
    assert config.n_pools == 1
    assert config.mismatch_product_size == 1


def test_coverage_rejects_boolean_homopolymer_limit_before_coercion():
    with pytest.raises(ValueError, match="primer_homopolymer_max.*integer"):
        Config(
            selection_algorithm="coverage",
            amplicon_size=200,
            amplicon_size_min=150,
            amplicon_size_max=280,
            amplicon_size_metric="reference-span",
            mismatch_product_size=2000,
            primer_homopolymer_max=True,
        )


def test_coverage_rejects_fractional_primer_walk_before_coercion():
    with pytest.raises(ValueError, match="primer_max_walk.*integer"):
        Config(
            selection_algorithm="coverage",
            amplicon_size=200,
            amplicon_size_min=150,
            amplicon_size_max=280,
            amplicon_size_metric="reference-span",
            mismatch_product_size=2000,
            primer_max_walk=1.7,
        )


def test_coverage_rejects_unsupported_edit_distance_at_config_boundary():
    with pytest.raises(ValueError, match="single-mismatch"):
        Config(
            selection_algorithm="coverage",
            amplicon_size=200,
            amplicon_size_min=150,
            amplicon_size_max=280,
            amplicon_size_metric="reference-span",
            mismatch_product_size=2000,
            editdist_max=2,
        )


def test_coverage_preserves_fractional_hairpin_threshold_without_changing_legacy():
    config = Config(
        selection_algorithm="coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
        primer_hairpin_th_max=51.9,
    )
    assert config.primer_hairpin_th_max == 51.9
    assert type(config.primer_hairpin_th_max) is float
    assert (
        ConstraintProfile.from_config(config).to_dict()["thermochemistry"][
            "primer_hairpin_th_max"
        ]
        == 51.9
    )
    legacy = Config(
        primer_hairpin_th_max=51.9,
        primer_homopolymer_max=True,
        primer_max_walk=1.7,
        editdist_max=2,
    )
    assert legacy.primer_hairpin_th_max == 51
    assert legacy.primer_homopolymer_max == 1
    assert legacy.primer_max_walk == 1
    assert legacy.editdist_max == 2


@pytest.mark.parametrize("flag", ["--dimer-score", "--min-base-freq"])
def test_coverage_cli_rejects_nan_discovery_settings_before_workflow(tmp_path, flag):
    source = _input(tmp_path)
    output = tmp_path / "out"
    args = [
        "panel-create",
        "--msa",
        str(source),
        "--output",
        str(output),
        "--mode",
        "equal",
        "--selection-algorithm",
        "coverage",
        "--amplicon-size",
        "200",
        "--amplicon-size-min",
        "150",
        "--amplicon-size-max",
        "280",
        flag,
        "nan",
    ]
    with patch(
        "primalscheme3.cli.panelcreate", side_effect=AssertionError("workflow called")
    ):
        result = CliRunner().invoke(app, args)
    assert result.exit_code == 2, result.output + repr(result.exception)
    assert not output.exists()
