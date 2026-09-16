"""Focused contracts for bounded native discovery replay."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from primalscheme3.panel.coverage_types import (
    CandidateFamily,
    OligoSite,
    Target,
    VariantCatalog,
)


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
