"""Public CLI contract for the opt-in coverage panel selector."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from click import UsageError
from typer.testing import CliRunner

from primalscheme3.cli import app
from primalscheme3.core.config import Config
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
                "panel-create", "--msa", str(msa), "--output", str(tmp_path / "out"),
                "--mispriming-product-size", "0",
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
