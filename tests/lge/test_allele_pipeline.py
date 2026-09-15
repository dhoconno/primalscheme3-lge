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


def test_search_failure_retains_complete_discovery_catalog(tmp_path, monkeypatch):
    import gzip

    import pytest

    from primalscheme3.panel import allele_pipeline
    from primalscheme3.panel.coverage_types import VariantCatalog

    msa = tmp_path / "short.fasta"
    msa.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "failed-search"
    config = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=250,
        ncores=1,
    )

    def fail(*args, **kwargs):
        raise RuntimeError("injected search failure")

    monkeypatch.setattr(allele_pipeline, "search_allele_assignments", fail)
    with pytest.raises(RuntimeError, match="injected search failure"):
        panelcreate(
            [msa],
            out,
            config,
            None,
            mode=PanelRunModes.EQUAL,
            executed_argv=["primalscheme3", "panel-create"],
        )
    with gzip.open(out / "discovery-catalog.json.gz", "rt") as handle:
        catalog = VariantCatalog.from_dict(json.load(handle))
    assert len(catalog.targets) == 1
    assert (out / "history/history.sqlite").is_file()
    receipt = json.loads((out / "panel-provenance.json").read_text())
    assert receipt["exitStatus"] != 0
    assert any(p["path"] == "discovery-catalog.json.gz" for p in receipt["outputs"])


def test_pipeline_reuse_skips_discovery_and_preserves_origin_history(
    tmp_path, monkeypatch
):
    from primalscheme3.panel import allele_pipeline
    from primalscheme3.panel.allele_catalog_cache import export_panel_discovery_cache
    from primalscheme3.panel.allele_inspection import (
        audit_allele_bundle,
        query_allele_history,
    )

    msa = tmp_path / "short.fasta"
    msa.write_text(">short\n" + "A" * 30 + "\n")

    def config(**extra):
        return Config(
            selection_algorithm="allele-coverage",
            amplicon_size=200,
            amplicon_size_min=150,
            amplicon_size_max=250,
            ncores=1,
            optimizer_starts=1,
            optimizer_repair_rounds=0,
            **extra,
        )

    original, cache, reused = (
        tmp_path / name for name in ("original", "cache", "reused")
    )
    panelcreate(
        [msa],
        original,
        config(),
        None,
        mode=PanelRunModes.EQUAL,
        executed_argv=["primalscheme3", "panel-create"],
    )
    export_panel_discovery_cache(original, cache, argv=["primalscheme3", "panel-cache"])

    def never(*args, **kwargs):
        raise AssertionError("discovery must not run")

    monkeypatch.setattr(allele_pipeline, "build_variant_catalog", never)
    panelcreate(
        [msa],
        reused,
        config(reuse_discovery=str(cache)),
        None,
        mode=PanelRunModes.EQUAL,
        executed_argv=["primalscheme3", "panel-create"],
    )
    saved = json.loads((reused / "config.json").read_text())
    assert saved["discovery_core_count"] == 0
    assert saved["discovery_reused"] is True
    optimizer = json.loads((reused / "panel-optimizer.json").read_text())
    assert (
        optimizer["history"]["originDiscovery"]["path"]
        == "origin-discovery/manifest.json"
    )
    assert audit_allele_bundle(reused)["valid"]
    result = query_allele_history(reused, stage="discovery", limit=2)
    assert "origin_discovery" in result
    assert not result["origin_discovery"]["scope"]["current_decisions"]
