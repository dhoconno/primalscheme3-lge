"""Portable reuse preserves complete origin history and profile equivalence."""

import gzip
import json
from types import SimpleNamespace

import pytest

from primalscheme3.core.config import Config
from primalscheme3.core.mapping import create_mapping
from primalscheme3.core.msa import parse_msa
from primalscheme3.panel.coverage_discovery import (
    build_variant_catalog,
    discovery_profiles,
    variant_targets,
)
from primalscheme3.panel.coverage_history import CoverageHistory, SQLiteCoverageHistory
from primalscheme3.panel.coverage_types import canonical_json


def replace_manifest_receipts(cache, manifest):
    """Tamper with both local receipts to exercise the original-source binding."""
    import hashlib

    manifest_path = cache / "manifest.json"
    manifest_path.write_text(canonical_json(manifest))
    receipt_path = cache / "cache-export-provenance.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["outputs"] = dict(manifest["artifacts"]) | {
        "manifest.json": {
            "sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            "size": manifest_path.stat().st_size,
        }
    }
    receipt_path.write_text(canonical_json(receipt))


def fixture_panel(path):
    from time import monotonic

    from primalscheme3.panel.coverage_provenance import (
        capture_execution_identity,
        finalize_provenance,
    )

    path.mkdir()
    (path / "work").mkdir()
    (path / "stages/strict").mkdir(parents=True)
    raw = path / "work/0000-input.fasta"
    shared = "CAACGGCGGACTTTATTGTATCTCC"
    seq = shared + "ATATATATAT" + shared.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    raw.write_text(">one\n" + seq + "\n>two\n" + seq + "\n")
    array, _ = parse_msa(raw)
    mapping, _ = create_mapping(array)
    targets = variant_targets(
        {0: SimpleNamespace(array=array, _mapping_array=mapping, msa_index=0)}
    )
    config = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=60,
        amplicon_size_min=45,
        amplicon_size_max=75,
    )
    start = capture_execution_identity([raw])
    started = monotonic()
    history = SQLiteCoverageHistory(path / "history", run_id="fixture")
    catalog = build_variant_catalog(targets, config, history=history)
    history.close()
    with gzip.open(path / "stages/strict/catalog.json.gz", "wt") as f:
        f.write(canonical_json(catalog.to_dict()))
    with gzip.open(path / "configuration-ledger.json.gz", "wt") as f:
        f.write("{}")
    (path / "stages/strict/validation.json").write_text('{"valid":true}')
    from unittest.mock import patch

    # Other agents may edit unrelated source while this synthetic fixture runs.
    with patch(
        "primalscheme3.panel.coverage_provenance.source_identity",
        return_value=start["source"],
    ):
        finalize_provenance(
            output_dir=path,
            argv=["fixture"],
            resolved_options=config.to_dict(),
            inputs=[
                {
                    "sourcePath": str(raw),
                    "storedPath": "work/0000-input.fasta",
                    "sourceIndex": 0,
                }
            ],
            started_at=started,
            ended_at=monotonic(),
            status="success",
            exit_status=0,
            stderr="",
            scientific={
                "selectionAlgorithm": "allele-coverage",
                "validationValid": True,
            },
            execution_start=start,
        )
    return targets, config, catalog


def test_union_cache_projection_exactly_matches_fresh_profile_discovery(tmp_path):
    from primalscheme3.panel.allele_catalog_cache import (
        export_panel_discovery_cache,
        load_discovery_cache,
        materialize_cache_reuse,
    )

    targets, config, union = fixture_panel(tmp_path / "panel")
    cache = tmp_path / "cache"
    export_panel_discovery_cache(tmp_path / "panel", cache, argv=["export-cache"])
    for name in ("normal", "high-gc"):
        profiles = {name: discovery_profiles(config)[name]}
        from unittest.mock import patch

        with patch.object(
            SQLiteCoverageHistory,
            "_reload",
            side_effect=AssertionError("must not replay history"),
        ):
            reuse = load_discovery_cache(
                cache, targets=targets, config=config, profiles=profiles
            )
        fresh = build_variant_catalog(targets, config, profiles=profiles)
        assert fresh.sites and fresh.families
        assert reuse.catalog.to_dict() == fresh.to_dict()
        history = CoverageHistory(None, run_id="new-selection")
        destination = materialize_cache_reuse(reuse, tmp_path / name, history=history)
        assert (destination / "history/history.sqlite").read_bytes() == (
            cache / "history/history.sqlite"
        ).read_bytes()
        assert (
            len(history.events) == 1
            and history.events[0].kind == "reused-discovery-cache"
        )
        assert not history.assessments
    changed = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=60,
        amplicon_size_min=45,
        amplicon_size_max=75,
        optimizer_seed=7,
        specificity_terminal_k=19,
        mismatch_product_size=1500,
    )
    assert (
        load_discovery_cache(cache, targets=targets, config=changed).catalog.to_dict()
        == union.to_dict()
    )


def test_cache_rejects_scientific_mismatch_corruption_and_incomplete_source(tmp_path):
    from dataclasses import replace

    from primalscheme3.panel.allele_catalog_cache import (
        export_panel_discovery_cache,
        load_discovery_cache,
    )

    targets, config, _ = fixture_panel(tmp_path / "panel")
    cache = tmp_path / "cache"
    export_panel_discovery_cache(tmp_path / "panel", cache, argv=["export-cache"])
    wrong = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=60,
        amplicon_size_min=45,
        amplicon_size_max=75,
        min_base_freq=0.5,
    )
    with pytest.raises(ValueError, match="discovery settings"):
        load_discovery_cache(cache, targets=targets, config=wrong)
    with pytest.raises(ValueError, match="targets"):
        load_discovery_cache(
            cache, targets=(replace(targets[0], occurrence=3),), config=config
        )
    with (cache / "history/history.sqlite").open("ab") as f:
        f.write(b"corrupt")
    with pytest.raises(ValueError, match="artifact"):
        load_discovery_cache(cache, targets=targets, config=config)
    import hashlib

    manifest_path = cache / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    altered = (cache / "history/history.sqlite").read_bytes()
    manifest["artifacts"]["history/history.sqlite"] = {
        "sha256": hashlib.sha256(altered).hexdigest(),
        "size": len(altered),
    }
    replace_manifest_receipts(cache, manifest)
    with pytest.raises(ValueError, match="source receipt"):
        load_discovery_cache(cache, targets=targets, config=config)
    provenance = tmp_path / "panel/panel-provenance.json"
    data = json.loads(provenance.read_text())
    data["sourceChangedDuringRun"] = True
    provenance.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="stable successful"):
        export_panel_discovery_cache(
            tmp_path / "panel", tmp_path / "bad", argv=["export-cache"]
        )


def test_portable_cache_rejects_other_discovery_policies_and_keeps_local_origins(
    tmp_path,
):
    import shutil

    from primalscheme3.panel.allele_catalog_cache import (
        export_panel_discovery_cache,
        load_discovery_cache,
        materialize_cache_reuse,
    )

    targets, config, _ = fixture_panel(tmp_path / "panel")
    cache = tmp_path / "cache"
    export_panel_discovery_cache(tmp_path / "panel", cache, argv=["export-cache"])
    shutil.rmtree(tmp_path / "panel")
    relocated = tmp_path / "relocated"
    cache.rename(relocated)
    for changes in ({"discovery_length_mode": "all"}, {"amplicon_size_min": 46}):
        changed = Config(
            selection_algorithm="allele-coverage",
            amplicon_size=60,
            amplicon_size_min=45,
            amplicon_size_max=75,
            **{k: v for k, v in changes.items() if k != "amplicon_size_min"},
        )
        if "amplicon_size_min" in changes:
            changed.amplicon_size_min = 46
        with pytest.raises(ValueError, match="discovery settings"):
            load_discovery_cache(relocated, targets=targets, config=changed)
    reuse = load_discovery_cache(
        relocated,
        targets=targets,
        config=config,
        profiles={"normal": discovery_profiles(config)["normal"]},
    )
    history = CoverageHistory(None, run_id="fresh")
    local = materialize_cache_reuse(reuse, tmp_path / "new-run", history=history)
    assert {e.id for e in history.evidence} == {e.id for e in reuse.derived_evidence}
    assert all(
        (local / item["storedPath"]).is_file() for item in reuse.manifest["inputs"]
    )
    shutil.rmtree(relocated)
    assert (local / "history/history.sqlite").is_file()
    # A local cache remains reopenable after both source panel and cache vanish.
    assert (
        load_discovery_cache(local, targets=targets, config=config).catalog.targets
        == reuse.catalog.targets
    )


def test_cache_rejects_manifest_path_escape_and_kernel_identity_change(tmp_path):
    from primalscheme3.panel.allele_catalog_cache import (
        export_panel_discovery_cache,
        load_discovery_cache,
    )

    targets, config, _ = fixture_panel(tmp_path / "panel")
    cache = tmp_path / "cache"
    manifest = export_panel_discovery_cache(
        tmp_path / "panel", cache, argv=["export-cache"]
    )
    original = canonical_json(manifest)
    manifest["scientificIdentity"]["kernels"][0]["version"] = "wrong"
    replace_manifest_receipts(cache, manifest)
    with pytest.raises(ValueError, match="kernels"):
        load_discovery_cache(cache, targets=targets, config=config)
    manifest = json.loads(original)
    manifest["artifacts"]["../outside"] = {"sha256": "x", "size": 0}
    replace_manifest_receipts(cache, manifest)
    with pytest.raises(ValueError, match="escapes"):
        load_discovery_cache(cache, targets=targets, config=config)


def test_export_preflight_failures_have_receipts_without_mutating_source(tmp_path):
    from primalscheme3.panel.allele_catalog_cache import export_panel_discovery_cache

    panel = tmp_path / "invalid-panel"
    panel.mkdir()
    original = panel / "panel-provenance.json"
    original.write_text("not-json")
    output = tmp_path / "failed-cache"
    with pytest.raises(ValueError):
        export_panel_discovery_cache(
            panel, output, argv=["panel-cache", "--source", str(panel)]
        )
    receipt = json.loads((output / "cache-export-provenance.json").read_text())
    assert receipt["status"] == "failure" and receipt["exitStatus"] == 1
    assert receipt["inputs"][0]["path"] == str(original)
    assert receipt["inputs"][0]["size"] == len("not-json")
    assert receipt["stderr"] and receipt["shell"]
    assert original.read_text() == "not-json"
    with pytest.raises(ValueError, match="separate"):
        export_panel_discovery_cache(panel, panel / "nested", argv=["panel-cache"])
    assert not (panel / "nested").exists()
    before = (output / "cache-export-provenance.json").read_bytes()
    with pytest.raises(FileExistsError):
        export_panel_discovery_cache(panel, output, argv=["panel-cache"])
    assert (output / "cache-export-provenance.json").read_bytes() == before
