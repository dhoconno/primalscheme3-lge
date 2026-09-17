from pathlib import Path
import sys
from dataclasses import dataclass
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from primalscheme3.cli import app
from primalscheme3.panel.gap_completion_cli import (
    resolve_gap_completion_parent,
    resolve_gap_expansion_options,
)


def _msa(tmp_path: Path) -> Path:
    path = tmp_path / "input.fasta"
    path.write_text(">reference\nACGTACGT\n")
    return path


@pytest.mark.parametrize("n_pools", [1, 2, 3])
def test_gap_completion_forwards_existing_pool_count(tmp_path, n_pools):
    parent = tmp_path / "parent"
    parent.mkdir()
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app,
            [
                "panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / f"out-{n_pools}"),
                "--selection-algorithm", "legacy", "--mode", "equal", "--mapping", "first",
                "--n-pools", str(n_pools), "--gap-completion-parent", str(parent),
            ],
        )
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["config"].n_pools == n_pools
    assert run.call_args.kwargs["gap_completion_parent"] == parent.resolve()


def test_gap_completion_is_absent_from_ordinary_legacy_call(tmp_path):
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app, ["panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / "out")]
        )
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["gap_completion_parent"] is None
    assert run.call_args.kwargs["gap_expansion_options"] is None


def test_gap_expansion_forwards_bounded_options_and_parent(tmp_path):
    parent = tmp_path / "parent"
    parent.mkdir()

    @dataclass(frozen=True)
    class FakeOptions:
        mode: str
        max_anchors_per_msa: int
        max_pairs_per_msa: int

    fake_module = SimpleNamespace(GapExpansionOptions=FakeOptions)
    with patch.dict(sys.modules, {"primalscheme3.panel.gap_expansion": fake_module}):
        with patch("primalscheme3.cli.panelcreate") as run:
            result = CliRunner().invoke(
                app,
                [
                    "panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / "expanded"),
                    "--selection-algorithm", "legacy", "--mode", "equal", "--mapping", "first",
                    "--gap-completion-parent", str(parent), "--gap-expansion", "bounded",
                    "--gap-expansion-max-anchors-per-msa", "12", "--gap-expansion-max-pairs-per-msa", "7",
                ],
            )
    assert result.exit_code == 0, result.output
    options = run.call_args.kwargs["gap_expansion_options"]
    assert options == FakeOptions("bounded", 12, 7)


def test_gap_expansion_advanced_options_require_bounded_mode(tmp_path):
    parent = tmp_path / "parent"
    parent.mkdir()
    with patch("primalscheme3.cli.panelcreate"):
        result = CliRunner().invoke(
            app,
            [
                "panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / "bad"),
                "--selection-algorithm", "legacy", "--mode", "equal", "--mapping", "first",
                "--gap-completion-parent", str(parent), "--gap-expansion-max-pairs-per-msa", "7",
            ],
        )
    assert result.exit_code != 0
    assert "require --gap-expansion bounded" in result.output


@pytest.mark.parametrize(
    "extra, message",
    [
        (["--selection-algorithm", "coverage"], "legacy"),
        (["--mode", "entropy"], "--mode equal"),
        (["--mapping", "consensus"], "mapping"),
    ],
)
def test_gap_completion_rejects_unsupported_scope(tmp_path, extra, message):
    parent = tmp_path / "parent"
    parent.mkdir()
    values = {"selection_algorithm": "legacy", "mode": "equal", "mapping": "first"}
    if extra[0] == "--selection-algorithm":
        values["selection_algorithm"] = extra[1]
    elif extra[0] == "--mode":
        values["mode"] = extra[1]
    else:
        values["mapping"] = extra[1]
    with pytest.raises(ValueError, match=message):
        resolve_gap_completion_parent(
            parent=parent, terminal_gap_policy="legacy", circular=False,
            region_bedfile=None, input_bedfile=None, legacy_salvage=None, **values,
        )
