from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from primalscheme3.cli import app


def _msa(tmp_path: Path) -> Path:
    path = tmp_path / "input.fasta"
    path.write_text(">reference\nACGTACGT\n")
    return path


def test_legacy_salvage_defaults_do_not_reach_panelcreate(tmp_path):
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app, ["panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / "out")]
        )
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["legacy_salvage_options"] is None


def test_legacy_salvage_resolves_native_options(tmp_path):
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(
            app,
            [
                "panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / "out"),
                "--selection-algorithm", "legacy", "--mode", "equal", "--mapping", "first",
                "--legacy-salvage", "bounded", "--legacy-salvage-threshold", "-28",
                "--legacy-salvage-threshold", "-31", "--legacy-salvage-floor", "-32",
                "--legacy-salvage-max-edges-per-pool", "3",
                "--legacy-salvage-max-incident-species-per-pool", "2",
                "--legacy-salvage-min-reference-gain", "2",
                "--legacy-salvage-max-candidate-evaluations", "77",
            ],
        )
    assert result.exit_code == 0, result.output
    options = run.call_args.kwargs["legacy_salvage_options"]
    assert options.mode == "bounded"
    assert options.thresholds == (-28.0, -31.0)
    assert options.floor == -32.0
    assert options.max_edges_per_pool == 3
    assert options.max_incident_species_per_pool == 2
    assert options.min_reference_gain == 2
    assert options.max_candidate_evaluations == 77


@pytest.mark.parametrize(
    "extra, message",
    [
        (["--legacy-salvage", "bounded"], "--mode equal"),
        (["--legacy-salvage", "bounded", "--mode", "equal", "--mapping", "consensus"], "--mapping first"),
        (["--legacy-salvage", "bounded", "--mode", "equal", "--region-bedfile", "regions.bed"], "region"),
        (["--legacy-salvage", "off", "--legacy-salvage-floor", "-32"], "require --legacy-salvage bounded"),
    ],
)
def test_legacy_salvage_rejects_unsupported_or_misconfigured_requests(tmp_path, extra, message):
    args = ["panel-create", "--msa", str(_msa(tmp_path)), "--output", str(tmp_path / "out"), *extra]
    if "regions.bed" in args:
        Path(args[args.index("regions.bed")]).write_text("reference\t0\t8\tr\n")
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 2
    assert message.lower() in result.output.lower()
