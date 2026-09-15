"""Read-only bounded history inspection and fresh portable bundle audits."""

import json
import shutil
import sqlite3

import pytest
from typer.testing import CliRunner

from primalscheme3.panel.coverage_history import SQLiteCoverageHistory


def api():
    from primalscheme3.panel import allele_inspection

    return allele_inspection


def history_bundle(path):
    path.mkdir()
    (path / "panel-optimizer.json").write_text(
        json.dumps({"history": {"path": "history/history.sqlite"}})
    )
    with SQLiteCoverageHistory(path / "history", run_id="test") as h:
        evidence = h.record_evidence(
            entity_ids=("site",),
            measurement="native",
            dependency_key={},
            values={"score": -28},
            status="evaluated",
        )
        parent = None
        for stage, outcome in [("strict", "fail"), ("salvage-1", "pass")]:
            assessment = h.assess(
                stage_id=stage,
                entity_ids=("site",),
                pool=0,
                context_digest=stage,
                profile_id="normal",
                kernel_versions={},
                thresholds={},
                check_name="dimer",
                outcome=outcome,
                reason=outcome,
                evidence_ids=(evidence.id,),
            )
            parent = h.emit(
                stage_id=stage,
                kind="assessed",
                entity_ids=("site",),
                assessment_ids=(assessment.id,),
                parent_event_ids=(parent.id,) if parent else (),
            )
        for i in range(2000):
            h.emit(stage_id="unrelated", kind="noise", entity_ids=(f"noise-{i}",))
    return path


def test_indexed_history_filters_linked_evidence_and_bounded_lineage(
    tmp_path, monkeypatch
):
    bundle = history_bundle(tmp_path / "bundle")
    monkeypatch.setattr(
        SQLiteCoverageHistory,
        "__init__",
        lambda *a, **k: pytest.fail("full reload forbidden"),
    )
    report = api().query_allele_history(
        bundle,
        entity="site",
        stage="salvage-1",
        pool=1,
        profile="normal",
        lineage=True,
        limit=1,
    )
    assert report["assessments"][0]["outcome"] == "pass"
    assert len(report["evidence"]) == 1
    assert {e["stage_id"] for e in report["events"]} == {"strict", "salvage-1"}
    assert report["query"]["pool"] == 1
    assert report["scope"]["full_history_verified"] is False
    page = api().query_allele_history(bundle, stage="unrelated", limit=3, offset=2)
    assert len(page["events"]) == 3
    assert page["pagination"]["events"]["truncated"]
    assert page["events"][0]["sequence_number"] == 4
    with pytest.raises(ValueError, match="filter"):
        api().query_allele_history(bundle)


def test_query_failure_writes_receipt_and_does_not_mutate_bundle(tmp_path):
    bundle = history_bundle(tmp_path / "bundle")
    database = bundle / "history/history.sqlite"
    before = database.read_bytes()
    out = tmp_path / "receipt"
    result = api().run_inspection(
        "panel-history",
        bundle,
        out,
        argv=["primalscheme3", "panel-history"],
        entity="site",
        pool=0,
    )
    assert not result["valid"]
    receipt = json.loads((out / "provenance.json").read_text())
    assert receipt["tool"] == "primalscheme3 panel-history"
    assert receipt["exitStatus"] != 0 and receipt["stderr"]
    assert receipt["inputs"][0]["sha256"]
    assert database.read_bytes() == before


@pytest.fixture
def bundle(tmp_path):
    from primalscheme3.core.config import Config
    from primalscheme3.panel.panel_classes import PanelRunModes
    from primalscheme3.panel.panel_main import panelcreate

    msa = tmp_path / "input.fa"
    msa.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "bundle"
    panelcreate(
        [msa],
        out,
        Config(
            selection_algorithm="allele-coverage",
            amplicon_size=200,
            amplicon_size_min=150,
            amplicon_size_max=250,
            ncores=1,
            optimizer_starts=1,
            optimizer_repair_rounds=0,
        ),
        None,
        mode=PanelRunModes.EQUAL,
        executed_argv=["primalscheme3", "panel-create"],
    )
    return out


def test_fresh_audit_survives_relocation_and_missing_original(bundle, tmp_path):
    moved = tmp_path / "moved"
    shutil.move(bundle, moved)
    (tmp_path / "input.fa").unlink()
    result = api().audit_allele_bundle(moved)
    assert result["valid"], result
    assert result["stages"]["strict"]["valid"]
    assert result["raw_inputs_reparsed"]
    with pytest.raises(ValueError, match="tier"):
        api().audit_allele_bundle(moved, tier="salvage-99")


def test_raw_input_change_detected_even_if_output_descriptors_updated(bundle):
    from primalscheme3.panel.allele_publication import artifact_descriptor

    prov = json.loads((bundle / "panel-provenance.json").read_text())
    record = prov["inputs"][0]
    raw = bundle / record["storedPath"]
    raw.write_text(">short\n" + "C" * 30 + "\n")
    desc = artifact_descriptor(raw, bundle)
    record.update(sha256=desc["sha256"], size=desc["size"])
    prov["outputs"] = [
        desc if d["path"] == desc["path"] else d for d in prov["outputs"]
    ]
    (bundle / "panel-provenance.json").write_text(json.dumps(prov))
    result = api().audit_allele_bundle(bundle)
    assert not result["valid"]
    assert any("input" in str(v) for v in result["violations"])


@pytest.mark.parametrize(
    "tamper", ["runtimeChangedDuringRun", "sourceChangedDuringRun"]
)
def test_audit_fails_changed_execution_identity(bundle, tamper):
    p = bundle / "panel-provenance.json"
    data = json.loads(p.read_text())
    data[tamper] = True
    p.write_text(json.dumps(data))
    assert not api().audit_allele_bundle(bundle)["valid"]


def test_cli_audit_failure_retains_provenance(bundle, tmp_path):
    from primalscheme3.cli import app

    (bundle / "primer.bed").write_text("corrupt\n")
    out = tmp_path / "audit"
    result = CliRunner().invoke(
        app, ["panel-audit", "--bundle", str(bundle), "--output", str(out)]
    )
    assert result.exit_code != 0
    assert (out / "validation.json").is_file()
    assert json.loads((out / "provenance.json").read_text())["exitStatus"] != 0


def test_region_queries_include_unselected_catalog_sites_and_grouped_conditions(
    tmp_path,
):
    from dataclasses import replace

    from test_allele_validation import dna, fixture

    from primalscheme3.panel.allele_publication import _write
    from primalscheme3.panel.coverage_discovery import HISTORY_ORIGIN_POLICY

    cat, configs, ledger = fixture(seq=dna(seed=0))
    site = cat.sites[0]
    rejected = replace(site, sequence="C" * 20, accepting_profile_ids=())
    cat = replace(cat, sites=cat.sites + (rejected,))
    ledger = replace(ledger, catalog_digest=cat.semantic_digest)
    bundle = tmp_path / "catalog-bundle"
    (bundle / "stages/strict").mkdir(parents=True)
    _write(bundle / "stages/strict/catalog.json.gz", cat.to_dict())
    _write(bundle / "configuration-ledger.json.gz", ledger.to_dict())
    (bundle / "panel-optimizer.json").write_text(
        json.dumps(
            {
                "history": {"path": "history/history.sqlite"},
                "publication": {"targetToReference": {site.target_id: "ref"}},
            }
        )
    )
    with SQLiteCoverageHistory(bundle / "history", run_id="r") as h:
        ev = h.record_evidence(
            entity_ids=(rejected.id,),
            measurement="row-enumeration",
            dependency_key={},
            values={
                "origin_policy": HISTORY_ORIGIN_POLICY,
                "common": {
                    "checks": {"gc": {"passed": False}, "length": {"passed": True}}
                },
            },
            status="evaluated",
        )
        h.assess(
            stage_id="discovery",
            entity_ids=(rejected.id,),
            pool=None,
            context_digest="c",
            profile_id="normal",
            kernel_versions={},
            thresholds={},
            check_name="profile-and-reference",
            outcome="fail",
            reason="gc",
            evidence_ids=(ev.id,),
        )
    result = api().query_allele_history(
        bundle, target="ref", region=site.reference_footprint, stage="discovery"
    )
    assert rejected.id in result["assessments"][0]["entity_ids"]
    assert {r["check_name"] for r in result["discovery_conditions"]} == {"gc", "length"}
    assert result["context"]["matched_entities"] >= 2


def test_selected_corrupt_link_fails_closed(tmp_path):
    bundle = history_bundle(tmp_path / "bundle")
    db = sqlite3.connect(bundle / "history/history.sqlite")
    db.execute("DELETE FROM record_links WHERE kind='evidence'")
    db.commit()
    db.close()
    with pytest.raises(ValueError, match="link"):
        api().query_allele_history(bundle, entity="site", stage="strict")


def test_profile_filters_events_by_linked_assessment(tmp_path):
    bundle = history_bundle(tmp_path / "bundle")
    result = api().query_allele_history(bundle, entity="site", profile="nonexistent")
    assert result["assessments"] == []
    assert result["events"] == []


def test_raw_reparse_detects_forged_input_provenance(bundle):
    from primalscheme3.panel.allele_publication import artifact_descriptor

    p = bundle / "panel-provenance.json"
    prov = json.loads(p.read_text())
    record = prov["inputs"][0]
    raw = bundle / record["storedPath"]
    raw.write_text(">short\n" + "C" * 30 + "\n")
    desc = artifact_descriptor(raw, bundle)
    record.update(
        sha256=desc["sha256"],
        size=desc["size"],
        sourceAtStartSha256=desc["sha256"],
        sourceAtStartSize=desc["size"],
    )
    prov["outputs"] = [
        desc if d["path"] == desc["path"] else d for d in prov["outputs"]
    ]
    p.write_text(json.dumps(prov))
    result = api().audit_allele_bundle(bundle)
    assert any(v["reason"] == "raw-input-target-mismatch" for v in result["violations"])


def test_audit_rejects_unsafe_descriptor_without_reading_outside(bundle):
    p = bundle / "panel-provenance.json"
    prov = json.loads(p.read_text())
    prov["outputs"][0]["path"] = "../outside"
    p.write_text(json.dumps(prov))
    with pytest.raises(ValueError, match="unsafe"):
        api().audit_allele_bundle(bundle)


def test_nonempty_bundle_omitted_variant_fails_even_after_hash_refresh(tmp_path):
    from dataclasses import replace

    from test_allele_validation import dna, fixture, profile

    from primalscheme3.panel.allele_coverage import canonical_observations
    from primalscheme3.panel.allele_pipeline import _authoritative_targets
    from primalscheme3.panel.allele_publication import (
        _write,
        artifact_descriptor,
        publish_allele_stage,
    )
    from primalscheme3.panel.coverage_types import (
        AlleleAssignment,
        CandidateFamily,
        ConfigurationLedger,
    )
    from primalscheme3.panel.coverage_variants import make_configuration

    bundle = tmp_path / "selected"
    (bundle / "work").mkdir(parents=True)
    raw = bundle / "work/0.fa"
    raw.write_text(">first\n" + dna(seed=0) + "\n")
    target = _authoritative_targets(
        bundle, [{"storedPath": "work/0.fa", "sourceIndex": 0}]
    )[0]
    cat, _, _ = fixture(seq=dna(seed=0))
    sites = tuple(replace(s, target_id=target.id) for s in cat.sites)
    fam = CandidateFamily(
        target.id,
        cat.families[0].anchor_pair,
        tuple(s.id for s in sites if s.strand == "+"),
        tuple(s.id for s in sites if s.strand == "-"),
    )
    cat = replace(
        cat,
        targets=(target,),
        observations=canonical_observations(target),
        sites=sites,
        families=(fam,),
    )
    config = make_configuration(cat, fam.id, fam.forward_site_ids, fam.reverse_site_ids)
    ledger = ConfigurationLedger(cat.semantic_digest, (config,))
    publish_allele_stage(
        bundle / "stages/strict",
        catalog=cat,
        ledger=ledger,
        assignments=(AlleleAssignment(config.id, 0, "strict"),),
        profile=profile(),
        authoritative_targets=(target,),
        references={target.id: "first"},
    )
    stage = bundle / "stages/strict"
    for name in (
        "primer.bed",
        "primers.fasta",
        "order-sheet.tsv",
        "amplicon.bed",
        "primertrim.amplicon.bed",
        "reference.fasta",
    ):
        shutil.copyfile(stage / name, bundle / name)
    with SQLiteCoverageHistory(bundle / "history", run_id="run"):
        pass
    _write(bundle / "configuration-ledger.json.gz", ledger.to_dict())
    manifest = json.loads((stage / "stage.json").read_text())
    optimizer = {
        "metric": "observed-allele-primer-trimmed/v1",
        "profile": profile().to_dict(),
        "primaryTier": "strict",
        "history": {"path": "history/history.sqlite"},
        "stages": [
            {
                "stage_id": "strict",
                "path": "stages/strict",
                "optimizer": {"stage_policy": manifest["stage_policy"]},
            }
        ],
    }
    _write(bundle / "panel-optimizer.json", optimizer)
    _write(bundle / "config.json", {"msa_data": {"0": {"msa_path": "work/0.fa"}}})
    descriptor = artifact_descriptor(raw, bundle)
    prov = {
        "status": "success",
        "exitStatus": 0,
        "sourceChangedDuringRun": False,
        "runtimeChangedDuringRun": False,
        "source": {"sourceDigest": "source"},
        "sourceAtEnd": {"sourceDigest": "source"},
        "runtime": {"pythonVersion": "test"},
        "runtimeAtEnd": {"pythonVersion": "test"},
        "inputs": [
            {
                "storedPath": "work/0.fa",
                "sourceIndex": 0,
                "sha256": descriptor["sha256"],
                "size": descriptor["size"],
                "sourceAtStartSha256": descriptor["sha256"],
                "sourceAtStartSize": descriptor["size"],
            }
        ],
    }

    def refresh():
        prov["outputs"] = [
            artifact_descriptor(p, bundle)
            for p in sorted(bundle.rglob("*"))
            if p.is_file() and p.name != "panel-provenance.json"
        ]
        _write(bundle / "panel-provenance.json", prov)

    refresh()
    result = api().audit_allele_bundle(bundle)
    assert result["valid"], result
    bed = stage / "primer.bed"
    lines = bed.read_text().splitlines()
    bed.write_text("\n".join(lines[:-1]) + "\n")
    manifest["artifacts"]["primer.bed"] = artifact_descriptor(bed, stage)
    _write(stage / "stage.json", manifest)
    shutil.copyfile(bed, bundle / "primer.bed")
    refresh()
    result = api().audit_allele_bundle(bundle)
    assert not result["valid"]
    assert not result["stages"]["strict"]["valid"]


def test_cli_region_failure_receipt_contains_requested_options(tmp_path):
    from primalscheme3.cli import app

    bundle = history_bundle(tmp_path / "bundle")
    out = tmp_path / "query"
    result = CliRunner().invoke(
        app,
        [
            "panel-history",
            "--bundle",
            str(bundle),
            "--entity",
            "site",
            "--region",
            "wrong",
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 1
    receipt = json.loads((out / "provenance.json").read_text())
    assert receipt["resolvedOptions"]["region"] == "wrong"
    assert receipt["tool"] == "primalscheme3 panel-history"


def test_unrelated_corrupt_payload_is_not_decoded_by_bounded_query(tmp_path):
    bundle = history_bundle(tmp_path / "bundle")
    db = sqlite3.connect(bundle / "history/history.sqlite")
    db.execute(
        "UPDATE records SET payload=? WHERE stage_id='unrelated'", (b"not zlib",)
    )
    db.commit()
    db.close()
    result = api().query_allele_history(bundle, entity="site", stage="strict", limit=1)
    assert result["valid"] and len(result["assessments"]) == 1
    assert result["scope"]["full_history_verified"] is False


def test_query_rejects_forged_index_context(tmp_path):
    bundle = history_bundle(tmp_path / "bundle")
    db = sqlite3.connect(bundle / "history/history.sqlite")
    db.execute("UPDATE records SET profile_id='forged' WHERE stream='assessments'")
    db.commit()
    db.close()
    with pytest.raises(ValueError, match="index"):
        api().query_allele_history(bundle, entity="site", profile="forged")


def test_origin_history_is_separate_and_independently_paginated(tmp_path):
    from test_allele_validation import dna, fixture

    from primalscheme3.panel.allele_publication import _write, artifact_descriptor

    bundle = history_bundle(tmp_path / "bundle")
    cat, _, ledger = fixture(seq=dna(seed=0))
    origin = bundle / "origin-discovery"
    (origin / "history").mkdir(parents=True)
    shutil.copyfile(
        bundle / "history/history.sqlite", origin / "history/history.sqlite"
    )
    _write(origin / "catalog.json.gz", cat.to_dict())
    _write(origin / "configuration-ledger.json.gz", ledger.to_dict())
    manifest = {
        "schemaVersion": "primalscheme3.discovery-cache/v1",
        "catalogPath": "catalog.json.gz",
        "historyPath": "history/history.sqlite",
        "originLedgerPath": "configuration-ledger.json.gz",
        "catalogSemanticDigest": cat.semantic_digest,
        "historyRole": "immutable-origin-history-including-prior-selection; not-current-decisions",
        "artifacts": {},
    }
    for name in (
        "catalog.json.gz",
        "history/history.sqlite",
        "configuration-ledger.json.gz",
    ):
        d = artifact_descriptor(origin / name, origin)
        manifest["artifacts"][name] = {k: v for k, v in d.items() if k != "path"}
    _write(origin / "manifest.json", manifest)
    optimizer = json.loads((bundle / "panel-optimizer.json").read_text())
    optimizer["history"]["originDiscovery"] = {"path": "origin-discovery/manifest.json"}
    _write(bundle / "panel-optimizer.json", optimizer)
    result = api().query_allele_history(bundle, entity="site", stage="strict", limit=1)
    assert result["origin_discovery"]["assessments"][0]["outcome"] == "fail"
    assert result["origin_discovery"]["scope"]["current_decisions"] is False
    assert result["assessments"][0]["outcome"] == "fail"


def test_region_inside_amplicon_finds_family_without_overlapping_primers(tmp_path):
    from test_allele_validation import dna, fixture

    from primalscheme3.panel.allele_publication import _write

    cat, _, ledger = fixture(seq=dna(seed=0))
    bundle = tmp_path / "region"
    (bundle / "stages/strict").mkdir(parents=True)
    _write(bundle / "stages/strict/catalog.json.gz", cat.to_dict())
    _write(bundle / "configuration-ledger.json.gz", ledger.to_dict())
    _write(
        bundle / "panel-optimizer.json", {"history": {"path": "history/history.sqlite"}}
    )
    with SQLiteCoverageHistory(bundle / "history", run_id="r") as h:
        h.emit(stage_id="discovery", kind="family", entity_ids=(cat.families[0].id,))
    report = api().query_allele_history(bundle, region=(100, 110), stage="discovery")
    assert len(report["events"]) == 1


def test_missing_execution_identity_is_not_a_successful_audit(bundle):
    path = bundle / "panel-provenance.json"
    data = json.loads(path.read_text())
    for key in ("source", "sourceAtEnd", "runtime", "runtimeAtEnd"):
        data.pop(key)
    path.write_text(json.dumps(data))
    assert not api().audit_allele_bundle(bundle)["valid"]


def test_early_unsafe_history_path_retains_manifest_descriptors(tmp_path):
    bundle = history_bundle(tmp_path / "bundle")
    optimizer = bundle / "panel-optimizer.json"
    data = json.loads(optimizer.read_text())
    data["history"]["path"] = "../unsafe"
    optimizer.write_text(json.dumps(data))
    output = tmp_path / "failure"
    report = api().run_inspection(
        "panel-history",
        bundle,
        output,
        argv=["primalscheme3", "panel-history"],
        entity="site",
    )
    assert not report["valid"]
    provenance = json.loads((output / "provenance.json").read_text())
    assert any(d["path"] == "panel-optimizer.json" for d in provenance["inputs"])
