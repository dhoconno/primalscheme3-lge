"""Native workflow publishes complete empty science with durable provenance."""

import json

from primalscheme3.core.config import AmpliconSizeMetric, Config
from primalscheme3.panel.panel_classes import PanelRunModes
from primalscheme3.panel.panel_main import panelcreate


def test_empty_native_allele_workflow_has_portable_outputs_and_provenance(tmp_path):
    msa = tmp_path / "short.fasta"
    msa.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "output"
    config = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=250,
        amplicon_size_metric=AmpliconSizeMetric.REFERENCE_SPAN,
        ncores=1,
        optimizer_starts=1,
        optimizer_repair_rounds=0,
    )
    argv = [
        "primalscheme3",
        "panel-create",
        "--selection-algorithm",
        "allele-coverage",
        "--mode",
        "equal",
        "--msa",
        str(msa),
        "--output",
        str(out),
        "--amplicon-size-min",
        "150",
        "--amplicon-size-max",
        "250",
    ]
    panelcreate([msa], out, config, None, mode=PanelRunModes.EQUAL, executed_argv=argv)
    provenance = json.loads((out / "panel-provenance.json").read_text())
    assert provenance["exitStatus"] == 0
    assert (out / "history/history.sqlite").is_file()
    optimizer = json.loads((out / "panel-optimizer.json").read_text())
    assert optimizer["primaryTier"] == "strict"
    assert optimizer["metric"] == "observed-allele-primer-trimmed/v1"
    assert optimizer["stages"][0]["coverage"]["mean_coverage"] == 0
    assert not [
        line
        for line in (out / "primer.bed").read_text().splitlines()
        if not line.startswith("#")
    ]
    from primalscheme3.panel.allele_publication import audit_allele_stage

    assert audit_allele_stage(out / "stages/strict")["valid"]


def test_cancelled_input_loading_writes_failure_provenance(tmp_path, monkeypatch):
    import pytest
    from primalscheme3.panel import panel_main

    msa = tmp_path / "short.fasta"
    msa.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "cancelled"
    config = Config(
        selection_algorithm="allele-coverage",
        amplicon_size_min=150,
        amplicon_size_max=250,
        amplicon_size=200,
        amplicon_size_metric=AmpliconSizeMetric.REFERENCE_SPAN,
        ncores=1,
    )

    def cancel(**kwargs):
        raise KeyboardInterrupt("cancel during input loading")

    monkeypatch.setattr(panel_main, "PanelMSA", cancel)
    with pytest.raises(KeyboardInterrupt):
        panelcreate(
            [msa],
            out,
            config,
            None,
            mode=PanelRunModes.EQUAL,
            executed_argv=["primalscheme3", "panel-create"],
        )
    result = json.loads((out / "panel-provenance.json").read_text())
    assert result["exitStatus"] != 0
    assert result["inputs"][0]["storedMissing"] is False
