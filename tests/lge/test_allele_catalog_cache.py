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


def fixture_panel(path, format_version=2, history_detail="full"):
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
        discovery_history=history_detail,
    )
    start = capture_execution_identity([raw])
    started = monotonic()
    history = SQLiteCoverageHistory(path / "history", run_id="fixture", format_version=format_version)
    catalog = build_variant_catalog(
        targets, config, history=history, history_detail=history_detail
    )
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


@pytest.mark.parametrize("format_version", [1, 2])
def test_union_cache_projection_exactly_matches_fresh_profile_discovery(tmp_path, format_version):
    from primalscheme3.panel.allele_catalog_cache import (
        export_panel_discovery_cache,
        load_discovery_cache,
        materialize_cache_reuse,
    )

    targets, config, union = fixture_panel(tmp_path / "panel", format_version=format_version)
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
        history = SQLiteCoverageHistory(tmp_path / (name + "-history"), run_id="new-selection")
        assert history.format_version == 2
        destination = materialize_cache_reuse(reuse, tmp_path / name, history=history)
        assert (destination / "history/history.sqlite").read_bytes() == (
            cache / "history/history.sqlite"
        ).read_bytes()
        assert (
            len(history.events) == 1
            and history.events[0].kind == "reused-discovery-cache"
        )
        assert not history.assessments
        history.close()
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


def test_compact_cache_projects_membership_without_reading_sqlite_history(tmp_path, monkeypatch):
    from contextlib import contextmanager

    from primalscheme3.panel.allele_catalog_cache import (
        export_panel_discovery_cache,
        load_discovery_cache,
    )

    targets, config, _ = fixture_panel(tmp_path / "compact-panel", history_detail="compact")
    cache = tmp_path / "compact-cache"
    export_panel_discovery_cache(tmp_path / "compact-panel", cache, argv=["export-cache"])
    manifest = json.loads((cache / "manifest.json").read_text())
    assert manifest["discoveryHistoryDetail"] == "compact"
    assert manifest["discoveryHistoryScope"] == "compact-target-profile-summaries"
    profiles = {"normal": discovery_profiles(config)["normal"]}
    import primalscheme3.panel.allele_catalog_cache as cache_api

    calls = []
    original_connection = cache_api._connection

    @contextmanager
    def observe_connection(path):
        calls.append(path)
        with original_connection(path) as db:
            yield db

    with monkeypatch.context() as patcher:
        patcher.setattr(cache_api, "_connection", observe_connection)
        reuse = load_discovery_cache(
            cache,
            targets=targets,
            config=config,
            profiles=profiles,
            history_detail="compact",
        )
    assert len(calls) == 1  # snapshot integrity only; profile projection did not query SQLite
    fresh = build_variant_catalog(
        targets,
        config,
        profiles=profiles,
        history_detail="compact",
    )
    assert reuse.catalog.to_dict() == fresh.to_dict()
    with pytest.raises(ValueError, match="history detail"):
        load_discovery_cache(
            cache,
            targets=targets,
            config=config,
            profiles=profiles,
            history_detail="full",
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


def test_cache_invalidates_changed_config_dna_table_dependency(tmp_path, monkeypatch):
    from copy import deepcopy

    from primalscheme3.panel import allele_catalog_cache as cache_api

    targets, config, _ = fixture_panel(tmp_path / "panel")
    cache = tmp_path / "cache"
    cache_api.export_panel_discovery_cache(tmp_path / "panel", cache, argv=["cache"])
    changed = deepcopy(cache_api.source_identity())
    entry = next(
        f for f in changed["files"] if f["path"] == "primalscheme3/core/config.py"
    )
    entry["sha256"] = "changed-DNA-interpretation-table"
    monkeypatch.setattr(cache_api, "source_identity", lambda: changed)
    with pytest.raises(ValueError, match="implementation or kernels"):
        cache_api.load_discovery_cache(cache, targets=targets, config=config)


@pytest.mark.parametrize("identity", ["source", "runtime"])
def test_export_identity_drift_fails_with_exact_retained_output_receipt(
    tmp_path, monkeypatch, identity
):
    import hashlib
    from copy import deepcopy

    from primalscheme3.panel import allele_catalog_cache as cache_api

    targets, config, _ = fixture_panel(tmp_path / "panel")
    original = getattr(cache_api, identity + "_identity")()
    changed = deepcopy(original)
    changed["sourceDigest" if identity == "source" else "pythonVersion"] = (
        "changed-during-export"
    )
    calls = 0

    def drifting_identity():
        nonlocal calls
        calls += 1
        return original if calls == 1 else changed

    monkeypatch.setattr(cache_api, identity + "_identity", drifting_identity)
    cache = tmp_path / "cache"
    with pytest.raises(ValueError, match="identity changed"):
        cache_api.export_panel_discovery_cache(
            tmp_path / "panel", cache, argv=["cache"]
        )
    receipt = json.loads((cache / "cache-export-provenance.json").read_text())
    assert receipt["status"] == "failure" and receipt["exitStatus"] == 1
    assert receipt[identity + "ChangedDuringRun"] is True
    assert receipt[identity] == original and receipt[identity + "AtEnd"] == changed
    manifest = json.loads((cache / "manifest.json").read_text())
    assert manifest["export"]["source"] == receipt["source"]
    assert manifest["export"]["runtime"] == receipt["runtime"]
    assert manifest["export"]["exitStatus"] == 1
    for name, descriptor in receipt["outputs"].items():
        payload = (cache / name).read_bytes()
        assert descriptor == {
            "sha256": hashlib.sha256(payload).hexdigest(),
            "size": len(payload),
        }
    # Even a contradictory successful status cannot hide declared identity drift.
    receipt.update(status="success", exitStatus=0)
    (cache / "cache-export-provenance.json").write_text(canonical_json(receipt))
    with pytest.raises(ValueError, match="complete successfully"):
        cache_api.load_discovery_cache(cache, targets=targets, config=config)


@pytest.mark.parametrize("fallback", [False, True])
def test_cache_copy_mutations_and_relocation_are_independent(
    tmp_path, monkeypatch, fallback
):
    import hashlib
    import shutil

    from primalscheme3.panel import allele_catalog_cache as cache_api
    from primalscheme3.panel import immutable_copy

    if fallback:
        monkeypatch.setattr(immutable_copy, "_clonefile_api", lambda: None)
    source = tmp_path / "panel"
    targets, config, _ = fixture_panel(source)
    cache = tmp_path / "cache"
    cache_api.export_panel_discovery_cache(source, cache, argv=["export-cache"])
    source_db, cached_db = (
        source / "history/history.sqlite",
        cache / "history/history.sqlite",
    )
    expected = hashlib.sha256(cached_db.read_bytes()).hexdigest()
    assert source_db.stat().st_ino != cached_db.stat().st_ino
    with source_db.open("r+b") as handle:
        handle.write(b"modified-source")
    assert hashlib.sha256(cached_db.read_bytes()).hexdigest() == expected
    shutil.rmtree(source)
    reuse = cache_api.load_discovery_cache(cache, targets=targets, config=config)
    local = cache_api.materialize_cache_reuse(
        reuse, tmp_path / "first", history=CoverageHistory(None, run_id="first")
    )
    local_db = local / "history/history.sqlite"
    assert local_db.stat().st_ino != cached_db.stat().st_ino
    with local_db.open("r+b") as handle:
        handle.write(b"modified-destination")
    assert hashlib.sha256(cached_db.read_bytes()).hexdigest() == expected
    shutil.rmtree(local)
    second = cache_api.materialize_cache_reuse(
        reuse, tmp_path / "second", history=CoverageHistory(None, run_id="second")
    )
    shutil.rmtree(cache)
    moved = tmp_path / "moved-origin"
    second.rename(moved)
    assert (
        cache_api.load_discovery_cache(moved, targets=targets, config=config).catalog
        == reuse.catalog
    )


def test_copied_cache_bytes_still_verified_after_clone_or_fallback(
    tmp_path, monkeypatch
):
    from primalscheme3.panel import allele_catalog_cache as cache_api

    targets, config, _ = fixture_panel(tmp_path / "panel")
    cache_api.export_panel_discovery_cache(
        tmp_path / "panel", tmp_path / "cache", argv=["export-cache"]
    )
    reuse = cache_api.load_discovery_cache(
        tmp_path / "cache", targets=targets, config=config
    )
    original = cache_api.copy_immutable_file

    def altered(source, destination):
        result = original(source, destination)
        if destination.name == "history.sqlite":
            with destination.open("r+b") as handle:
                handle.write(b"changed-after-copy")
        return result

    monkeypatch.setattr(cache_api, "copy_immutable_file", altered)
    with pytest.raises(ValueError, match="checksum mismatch"):
        cache_api.materialize_cache_reuse(
            reuse, tmp_path / "result", history=CoverageHistory(None, run_id="tamper")
        )
