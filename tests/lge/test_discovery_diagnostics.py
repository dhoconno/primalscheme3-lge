"""Focused contracts for bounded native discovery replay."""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict
from types import SimpleNamespace

import pytest

from primalscheme3.core.config import Config
from primalscheme3.panel.coverage_types import (
    CandidateFamily,
    OligoSite,
    Target,
    VariantCatalog,
)


def _real_source_bundle(tmp_path, module):
    """Build the smallest descriptor-bound source panel for real preflight."""
    from primalscheme3.core.mapping import create_mapping
    from primalscheme3.core.msa import parse_msa
    from primalscheme3.panel.coverage_discovery import (
        build_variant_catalog,
        variant_targets,
    )
    from primalscheme3.panel.coverage_history import CoverageHistory

    bundle = tmp_path / "real-source"
    stage = bundle / "stages/strict"
    stage.mkdir(parents=True)
    raw = bundle / "input.fa"
    sequence = "CAACGGCGGACTTTATTGTATCTCC" * 4
    raw.write_text(">reference\n" + sequence + "\n")
    array, _ = parse_msa(raw)
    mapping, _ = create_mapping(array)
    target = variant_targets(
        {0: SimpleNamespace(array=array, _mapping_array=mapping, msa_index=0)}
    )[0]
    config = Config(
        selection_algorithm="allele-coverage",
        amplicon_size=100,
        amplicon_size_min=80,
        amplicon_size_max=120,
    )
    source_catalog = build_variant_catalog(
        (target,),
        config,
        indexes=([24], [76]),
        history=CoverageHistory(None, run_id="source"),
        history_detail="full",
    )
    with gzip.open(stage / "catalog.json.gz", "wt") as stream:
        json.dump(
            source_catalog.to_dict(), stream, sort_keys=True, separators=(",", ":")
        )
    with gzip.open(stage / "authoritative-targets.json.gz", "wt") as stream:
        json.dump([asdict(target)], stream, sort_keys=True, separators=(",", ":"))
    (bundle / "config.json").write_text(json.dumps(config.to_dict(), sort_keys=True))
    optimizer = {
        "schemaVersion": "primalscheme3.panel-optimizer/v2",
        "primaryTier": "strict",
        "stages": [{"stage_id": "strict", "path": "stages/strict"}],
        "options": {
            "candidate_profiles": "union",
            "discovery_length_mode": "first-compatible",
        },
    }
    (bundle / "panel-optimizer.json").write_text(json.dumps(optimizer, sort_keys=True))
    descriptors = [
        module._descriptor(path, bundle)
        for path in sorted(bundle.rglob("*"))
        if path.is_file()
    ]
    input_descriptor = next(item for item in descriptors if item["path"] == "input.fa")
    source = module.source_identity()
    runtime = module.runtime_identity()
    provenance = {
        "schemaVersion": "primalscheme3.panel-provenance/v1",
        "status": "success",
        "exitStatus": 0,
        "source": source,
        "sourceAtEnd": source,
        "runtime": runtime,
        "runtimeAtEnd": runtime,
        "sourceChangedDuringRun": False,
        "runtimeChangedDuringRun": False,
        "inputs": [
            {
                **input_descriptor,
                "sourcePath": str(raw.resolve()),
                "sourceAtStartSha256": input_descriptor["sha256"],
                "sourceAtStartSize": input_descriptor["size"],
                "storedPath": "input.fa",
                "sourceIndex": 0,
            }
        ],
        "outputs": descriptors,
    }
    (bundle / "panel-provenance.json").write_text(
        json.dumps(provenance, sort_keys=True)
    )
    return bundle, source_catalog, target, config


def _fixture():
    target = Target(
        "target-test",
        0,
        0,
        ("row-test",),
        (tuple("ACGT" * 20),),
        "ACGT" * 20,
        80,
        tuple(range(80)),
        tuple(range(80)) + (80,),
    )
    forward = OligoSite(
        target.id,
        "ACGTACGTACGTACGTACG",
        "+",
        20,
        (0, 20),
        accepting_profile_ids=("normal",),
        generated_profile_ids=("normal",),
    )
    reverse = OligoSite(
        target.id,
        "TACGTACGTACGTACGTAC",
        "-",
        60,
        (60, 80),
        accepting_profile_ids=("normal",),
        generated_profile_ids=("normal",),
    )
    family = CandidateFamily(
        target.id,
        (forward.alignment_anchor, reverse.alignment_anchor),
        (forward.id,),
        (reverse.id,),
        (("normal", "normal"),),
    )
    catalog = VariantCatalog(
        (target,),
        (),
        (forward, reverse),
        (family,),
        json.dumps(
            {
                "profiles": {"normal": {}},
                "amplicon_size_min": 20,
                "amplicon_size_max": 100,
            }
        ),
    )
    return catalog, target, family, forward


@pytest.fixture
def harness(tmp_path, monkeypatch):
    import primalscheme3.panel.discovery_diagnostics as module

    catalog, target, family, site = _fixture()
    source = tmp_path / "source"
    source.mkdir()
    source_input = source / "stored.fa"
    source_input.write_text(">x\nACGT\n")
    source_inputs = [{"path": source_input, "sourceIndex": 0}]
    provenance = {"schemaVersion": "primalscheme3.panel-provenance/v1"}
    optimizer = {
        "schemaVersion": "primalscheme3.panel-optimizer/v2",
        "primaryTier": "strict",
    }
    monkeypatch.setattr(
        module,
        "_source_preflight",
        lambda _: (provenance, optimizer, catalog, (target,), source_inputs),
    )
    monkeypatch.setattr(
        module,
        "_config_from_bundle",
        lambda *_: SimpleNamespace(discovery_length_mode="first-compatible"),
    )
    monkeypatch.setattr(
        module, "_profile_binding", lambda *_: {"profiles": ["normal"], "resolved": {}}
    )
    monkeypatch.setattr(
        module, "discovery_profiles", lambda _: {"normal": SimpleNamespace()}
    )
    monkeypatch.setattr(module, "source_identity", lambda: {"sourceDigest": "test"})
    monkeypatch.setattr(
        module, "runtime_identity", lambda: {"nativeKernels": [{"package": "test"}]}
    )
    return module, tmp_path, catalog, family, site, source_input


def test_family_replay_is_bounded_and_preserves_membership(harness, monkeypatch):
    module, tmp_path, source_catalog, family, _, source_input = harness
    calls = {}

    def replay(targets, config, **kwargs):
        calls.update(kwargs)
        return source_catalog

    monkeypatch.setattr(module, "_invoke_discovery", replay)
    output = tmp_path / "diagnostic"
    report = module.diagnose_discovery(
        bundle=tmp_path / "source",
        output=output,
        family_id=family.id,
        argv=["primalscheme3", "panel-discovery-diagnose"],
    )
    assert report["valid"]
    assert report["scope"]["allOriginalTargetRows"] is True
    assert report["scope"]["bounded"] is True
    assert calls["indexes"] == ([20], [60])
    assert (output / "history").is_dir()
    receipt = json.loads((output / "provenance.json").read_text())
    assert receipt["tool"] == "primalscheme3 panel-discovery-diagnose"
    assert receipt["status"] == "success"
    assert (output / "source/0000-stored.fa").read_bytes() == source_input.read_bytes()


def test_site_replay_uses_only_requested_strand_anchor(harness, monkeypatch):
    module, tmp_path, source_catalog, _, site, _ = harness
    calls = {}
    monkeypatch.setattr(
        module,
        "_invoke_discovery",
        lambda targets, config, **kwargs: (calls.update(kwargs) or source_catalog),
    )
    report = module.diagnose_discovery(
        bundle=tmp_path / "source", output=tmp_path / "site", site_id=site.id
    )
    assert report["valid"]
    assert calls["indexes"] == ([20], [])


@pytest.mark.parametrize(
    "kwargs", [{"family_id": "missing"}, {"family_id": "x", "site_id": "y"}, {}]
)
def test_invalid_entity_request_writes_failure_receipt_and_preserves_source(
    harness, kwargs
):
    module, tmp_path, _, _, _, source_input = harness
    before = source_input.read_bytes()
    output = tmp_path / "failed"
    report = module.diagnose_discovery(
        bundle=tmp_path / "source", output=output, **kwargs
    )
    assert report["valid"] is False
    receipt = json.loads((output / "provenance.json").read_text())
    assert receipt["status"] == "failure"
    assert receipt["exitStatus"] != 0
    assert receipt["stderr"]
    assert source_input.read_bytes() == before


def test_membership_drift_and_builder_failure_are_failed_provenance(
    harness, monkeypatch
):
    module, tmp_path, source_catalog, family, _, _ = harness
    changed = VariantCatalog(
        source_catalog.targets,
        source_catalog.observations,
        source_catalog.sites,
        (),
        source_catalog.resolved_config_json,
    )
    monkeypatch.setattr(module, "_invoke_discovery", lambda *args, **kwargs: changed)
    output = tmp_path / "drift"
    report = module.diagnose_discovery(
        bundle=tmp_path / "source", output=output, family_id=family.id
    )
    assert not report["valid"]
    assert json.loads((output / "provenance.json").read_text())["status"] == "failure"

    output = tmp_path / "error"
    monkeypatch.setattr(
        module,
        "_invoke_discovery",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("replay failed")),
    )
    report = module.diagnose_discovery(
        bundle=tmp_path / "source", output=output, family_id=family.id
    )
    assert not report["valid"]
    receipt = json.loads((output / "provenance.json").read_text())
    assert "replay failed" in receipt["stderr"]


def test_real_preflight_replay_and_config_binding(tmp_path):
    import primalscheme3.panel.discovery_diagnostics as module

    bundle, source_catalog, _, config = _real_source_bundle(tmp_path, module)
    loaded_provenance, _, loaded_catalog, _, _ = module._source_preflight(bundle)
    assert loaded_provenance["status"] == "success"
    assert loaded_catalog.semantic_digest == source_catalog.semantic_digest
    replay_config = module._config_from_bundle(bundle, {"options": {}})
    assert replay_config.amplicon_size_min == config.amplicon_size_min
    assert replay_config.amplicon_size_max == config.amplicon_size_max
    family = source_catalog.families[0]
    report = module.diagnose_discovery(
        bundle=bundle, output=tmp_path / "real-replay", family_id=family.id
    )
    assert report["valid"], report
    assert report["scope"]["allOriginalTargetRows"] is True
    receipt = json.loads((tmp_path / "real-replay/provenance.json").read_text())
    assert receipt["status"] == "success"
    assert receipt["sourceChangedDuringRun"] is False
    assert receipt["runtimeChangedDuringRun"] is False
    assert all(item["sourceChangedDuringRun"] is False for item in receipt["inputs"])


@pytest.mark.parametrize("tamper", ["input", "catalog", "runtime", "source"])
def test_real_preflight_rejects_tampered_or_missing_bindings(tmp_path, tamper):
    import primalscheme3.panel.discovery_diagnostics as module

    bundle, _, _, _ = _real_source_bundle(tmp_path, module)
    if tamper == "input":
        (bundle / "input.fa").write_text(">reference\nACGT\n")
    elif tamper == "catalog":
        with gzip.open(bundle / "stages/strict/catalog.json.gz", "wt") as stream:
            json.dump({"tampered": True}, stream)
    else:
        path = bundle / "panel-provenance.json"
        provenance = json.loads(path.read_text())
        provenance.pop("runtime" if tamper == "runtime" else "source")
        path.write_text(json.dumps(provenance))
    with pytest.raises(ValueError):
        module._source_preflight(bundle)


def test_replay_requires_new_history_detail_keyword(harness, monkeypatch):
    module, _, catalog, _, _, _ = harness
    calls = {}

    def builder(targets, config, **kwargs):
        calls.update(kwargs)
        return catalog

    monkeypatch.setattr(module, "build_variant_catalog", builder)
    module._invoke_discovery(
        (catalog.targets[0],),
        SimpleNamespace(discovery_length_mode="first-compatible"),
        profiles={"normal": SimpleNamespace()},
        indexes=([20], []),
        history=SimpleNamespace(),
    )
    assert calls["history_detail"] == "full"
