"""Synthetic region diagnosis; no real assay sequence fixtures."""

import importlib.util
import json
from pathlib import Path

import pytest
from test_allele_validation import dna, fixture, profile

from primalscheme3.panel.coverage_types import AlleleAssignment


def api():
    path = Path(__file__).parents[2] / "scripts/diagnose_allele_region.py"
    assert path.is_file(), "standalone region diagnostic is missing"
    spec = importlib.util.spec_from_file_location("region_diagnostic", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def legacy(cat, c):
    f, r = (cat.site_by_id[ids[0]] for ids in (c.forward_site_ids, c.reverse_site_ids))
    return {
        "catalogIndex": 4,
        "spanBED": list(c.full_interval),
        "forward": {
            "regionBED": list(f.reference_footprint),
            "sequences": [f.sequence],
        },
        "reverse": {
            "regionBED": list(r.reference_footprint),
            "sequences": [r.sequence],
        },
        "pools": {"1": {"acceptedAtCurrentSettings": True, "worstDimerScore": 0}},
    }


def test_same_class_region_coverage_and_exact_variant_losses():
    seq = dna(seed=0)
    # Reverse footprint mismatch: a forward-only hit must earn zero coverage.
    other = seq[:180] + "A" * 20 + seq[200:]
    cat, cs, ledger = fixture(seq=seq, rows=(seq, seq, other))
    assignments = (AlleleAssignment(cs[0].id, 0, "strict"),)
    report = api().region_coverage(cat, ledger, assignments, cat.targets[0], (0, 200))
    assert sorted(x["covered_bases"] for x in report["classes"]) == [0, 160]
    assert sorted(x["multiplicity"] for x in report["classes"]) == [1, 2]
    assert all(x["observed_bases"] == 200 for x in report["classes"])
    failed = next(x for x in report["classes"] if x["covered_bases"] == 0)
    assert any(x["status"] == "mismatch" for x in failed["site_support"])


def test_legacy_matching_requires_exact_sequence_anchor_and_existing_family():
    cat, cs, _ = fixture(seq=dna(seed=0))
    item = legacy(cat, cs[0])
    matched = api().match_alternative(cat, cat.targets[0], item)
    assert matched["configuration"].id == cs[0].id
    item["forward"]["regionBED"][1] += 1
    bad = api().match_alternative(cat, cat.targets[0], item)
    assert bad["status"] == "unavailable"
    assert bad["issues"][0]["reason"] == "no-exact-site"


def test_insertion_uses_current_panel_and_preserves_raw_dimer_and_assessments(tmp_path):
    from primalscheme3.panel.allele_validation import StagePolicy
    from primalscheme3.panel.coverage_history import CoverageHistory

    cat, cs, ledger = fixture(seq=dna(seed=0))
    chosen = (AlleleAssignment(cs[0].id, 0, "strict"),)
    history = CoverageHistory(None, run_id="diagnostic-test")
    result = api().evaluate_alternative(
        cat, ledger, chosen, cs[0], profile(), StagePolicy(), 0.95, history
    )
    assert result["pools"][0]["status"] == "already-selected"
    # Other pool is an actual insertion, never a historical verdict.
    assert result["pools"][1]["validation"]["assignments"] == [
        chosen[0].to_dict(),
        AlleleAssignment(cs[0].id, 1, "strict").to_dict(),
    ]
    assert result["pools"][1]["validation"]["pools"][1]["dimer_edges"]
    assert len(history.assessments) > 0
    assert result["unary"]["profile"]["specificity_revision"] == "selected-sites-v2"


def test_failure_receipt_and_fresh_output_guard(tmp_path):
    module = api()
    bundle = tmp_path / "source"
    bundle.mkdir()
    output = tmp_path / "failed"
    result = module.run(
        bundle,
        output,
        target_index=0,
        region=(0, 10),
        argv=["diagnose", "--bundle", str(bundle)],
    )
    assert not result["valid"]
    receipt = json.loads((output / "provenance.json").read_text())
    assert receipt["status"] == "failed"
    assert receipt["command"]["argv"][0] == "diagnose"
    assert receipt["outputs"][0]["sha256"]
    with pytest.raises(ValueError):
        module.run(bundle, output, target_index=0, region=(0, 10), argv=[])
    with pytest.raises(ValueError):
        module.run(bundle, bundle / "nested", target_index=0, region=(0, 10), argv=[])
    assert not (bundle / "nested").exists()


def test_completed_native_bundle_diagnosis_is_portable_and_bounded(
    tmp_path, monkeypatch
):
    import hashlib

    from primalscheme3.core.config import Config
    from primalscheme3.panel.coverage_history import SQLiteCoverageHistory
    from primalscheme3.panel.panel_main import PanelRunModes, panelcreate

    raw = tmp_path / "input.fasta"
    raw.write_text(">one\n" + "AT" * 40 + "\n>two\n" + "AT" * 40 + "\n")
    bundle = tmp_path / "panel"
    panelcreate(
        [raw],
        bundle,
        Config(
            selection_algorithm="allele-coverage",
            amplicon_size=60,
            amplicon_size_min=40,
            amplicon_size_max=80,
            optimizer_starts=1,
            optimizer_repair_rounds=0,
            optimizer_time_limit=0.1,
        ),
        None,
        mode=PanelRunModes.EQUAL,
        executed_argv=["primalscheme3", "panel-create"],
    )
    module = api()
    original = SQLiteCoverageHistory.__init__

    def guard(self, directory, *args, **kwargs):
        assert not Path(directory).resolve().is_relative_to(bundle.resolve()), (
            "full source history reload forbidden"
        )
        return original(self, directory, *args, **kwargs)

    monkeypatch.setattr(SQLiteCoverageHistory, "__init__", guard)
    output = tmp_path / "diagnosis"
    report = module.run(
        bundle,
        output,
        target_index=0,
        region=(10, 60),
        max_candidates=0,
        argv=["diagnose", "--target-index", "0", "--region", "10:60"],
    )
    assert report["valid"], report
    assert report["selected"]["classes"][0]["observed_bases"] == 50
    assert report["selected"]["classes"][0]["covered_bases"] == 0
    assert report["selected"]["classes"][0]["multiplicity"] == 2
    assert report["work"]["examined_unique_configurations"] == 0
    receipt = json.loads((output / "provenance.json").read_text())
    assert receipt["status"] == "success"
    for item in receipt["outputs"]:
        path = output / item["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"]


def test_budget_ranks_only_matching_configurations_and_deduplicates():
    rows = [
        {
            "status": "matched",
            "configuration_id": str(i),
            "primer_trimmed_reference_overlap": i,
            "source_index": 0,
            "ordinal": i,
            "legacy_candidate": {"catalogIndex": i},
        }
        for i in range(20)
    ]
    rows.append(dict(rows[-1], source_index=1))
    rows.append({"status": "unavailable", "primer_trimmed_reference_overlap": 1000})
    assert api().select_alternatives(rows, 16) == [str(i) for i in range(19, 3, -1)]


def test_projection_includes_internal_insertions_and_excludes_unknowns():
    from allele_fixtures import target

    from primalscheme3.panel.allele_coverage import canonical_observations
    from primalscheme3.panel.coverage_types import ConfigurationLedger, VariantCatalog

    t = target("AA--CCGG", "AATNCCGG")
    cat = VariantCatalog((t,), canonical_observations(t), (), ())
    report = api().region_coverage(
        cat, ConfigurationLedger(cat.semantic_digest, ()), (), t, (1, 4)
    )
    assert report["alignment_region"] == [1, 6]
    assert sorted(c["observed_bases"] for c in report["classes"]) == [3, 4]
    assert any(c["unknown_alignment_columns"] == [3] for c in report["classes"])


def test_unary_and_insertion_do_not_trust_historical_acceptance():
    from primalscheme3.panel.allele_validation import StagePolicy
    from primalscheme3.panel.coverage_history import CoverageHistory

    cat, cs, ledger = fixture(seq=dna(seed=0), spans=((0, 200), (30, 230)))
    item = legacy(cat, cs[1])
    matched = api().match_alternative(cat, cat.targets[0], item)
    history = CoverageHistory(None, run_id="fresh")
    result = api().evaluate_alternative(
        cat,
        ledger,
        (AlleleAssignment(cs[0].id, 0, "strict"),),
        matched["configuration"],
        profile(),
        StagePolicy(),
        0.95,
        history,
    )
    assert not result["pools"][0]["validation"]["valid"]
    assert any(
        v["reason"] == "same-target-overlap"
        for v in result["pools"][0]["validation"]["violations"]
    )


def test_diagnostic_persists_matching_budgeted_and_unavailable_origins(
    tmp_path, monkeypatch
):
    import gzip

    from primalscheme3.panel.allele_publication import publish_allele_stage
    from primalscheme3.panel.coverage_history import SQLiteCoverageHistory

    cat, cs, ledger = fixture(seq=dna(seed=0), spans=((0, 200), (30, 230)))
    bundle = tmp_path / "source"
    stage = bundle / "stages/strict"
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=(AlleleAssignment(cs[0].id, 0, "strict"),),
        profile=profile(),
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "synthetic"},
    )
    with gzip.open(bundle / "configuration-ledger.json.gz", "wt") as stream:
        json.dump(ledger.to_dict(), stream)
    (bundle / "panel-optimizer.json").write_text(
        json.dumps(
            {
                "primaryTier": "strict",
                "stages": [{"stage_id": "strict", "path": "stages/strict"}],
                "history": {"configurationLedger": "configuration-ledger.json.gz"},
            }
        )
    )
    unavailable = legacy(cat, cs[0])
    unavailable["forward"]["sequences"] = ["A" * 20]
    old = tmp_path / "legacy.json"
    old.write_text(
        json.dumps(
            {"candidates": [legacy(cat, cs[0]), legacy(cat, cs[1]), unavailable]}
        )
    )
    module = api()
    # Isolate report assembly with a real freshly audited published stage;
    # the full bundle/raw-input audit is exercised by the native end-to-end test.
    monkeypatch.setattr(module, "audit_allele_bundle", lambda b: {"valid": True})
    monkeypatch.setattr(
        module,
        "query_allele_history",
        lambda *a, **k: {"bounded": True, "limit": k["limit"]},
    )
    output = tmp_path / "diagnostic"
    output.mkdir()
    report = module.diagnose(
        bundle,
        output,
        target_index=0,
        region=(20, 180),
        legacy_alternatives=(old,),
        max_candidates=1,
    )
    assert report["work"]["total_legacy_records"] == 3
    assert report["work"]["examined_unique_configurations"] == 1
    assert report["work"]["omitted_unique_configurations"] == 1
    assert report["work"]["unavailable_legacy_records"] == 1
    assert report["alternatives"][1]["insertion"]["status"] == "not-evaluated-budget"
    with gzip.open(
        output / report["alternatives"][0]["insertion"]["path"], "rt"
    ) as handle:
        measured = json.load(handle)
    assert measured["configuration"]["id"] == cs[0].id
    assert measured["unary"]["pools"]["0"]["dimer_edges"]
    with SQLiteCoverageHistory(
        output / "history", run_id="diagnose-allele-region"
    ) as history:
        assert history.last_complete_stage.stage_id == "region-diagnostic"
        assert len(history.assessments) > 0


@pytest.mark.parametrize("changed", ["input", "source", "runtime"])
def test_identity_drift_fails_closed_with_retained_output_receipt(
    tmp_path, monkeypatch, changed
):
    module = api()
    bundle = tmp_path / "source"
    bundle.mkdir()
    (bundle / "panel-provenance.json").write_text(
        json.dumps({"status": "success", "exitStatus": 0})
    )
    legacy_path = tmp_path / "legacy.json"
    legacy_path.write_text("{}")
    calls = {"source": 0, "runtime": 0}

    def source():
        calls["source"] += 1
        return {
            "native": {
                "sourceDigest": str(calls["source"] if changed == "source" else 1)
            },
            "script": {"sha256": "same"},
        }

    def runtime():
        calls["runtime"] += 1
        return {"version": calls["runtime"] if changed == "runtime" else 1}

    monkeypatch.setattr(module, "_identity", source)
    monkeypatch.setattr(module, "runtime_identity", runtime)

    def diagnose(*args, **kwargs):
        if changed == "input":
            legacy_path.write_text('{"changed": true}')
        return {"valid": True}

    monkeypatch.setattr(module, "diagnose", diagnose)
    output = tmp_path / "result"
    report = module.run(
        bundle,
        output,
        target_index=0,
        region=(1, 2),
        legacy_alternatives=(legacy_path,),
        argv=["diagnose"],
    )
    assert not report["valid"]
    receipt = json.loads((output / "provenance.json").read_text())
    assert receipt["status"] == "failed" and receipt["exitStatus"] == 1
    assert (
        receipt["changedInputPaths"]
        or receipt["sourceChangedDuringRun"]
        or receipt["runtimeChangedDuringRun"]
    )
    assert receipt["outputs"][0]["sha256"]


def test_legacy_records_index_families_once_and_keep_ambiguity(tmp_path):
    cat, cs, _ = fixture(seq=dna(seed=0))

    # A proxy makes repeated whole-family enumeration observable without
    # changing the immutable catalog used by make_configuration.
    class Families:
        scans = 0

        def __iter__(self):
            self.scans += 1
            assert self.scans == 1, "repeated whole-catalog family scan"
            return iter(cat.families)

    families = Families()

    class Proxy:
        def __getattr__(self, name):
            return families if name == "families" else getattr(cat, name)

    old = tmp_path / "old.json"
    old.write_text(json.dumps({"candidates": [legacy(cat, cs[0])] * 8}))
    module = api()
    records = module._legacy_records((old,), Proxy(), cat.targets[0], (0, 200))
    assert len(records) == 8 and all(r["status"] == "matched" for r in records)
    assert families.scans == 1
    f = cat.families[0]
    ambiguous = module.match_alternative(
        cat,
        cat.targets[0],
        legacy(cat, cs[0]),
        family_index={(cat.targets[0].id, f.anchor_pair): [f, f]},
    )
    assert ambiguous["issues"][0]["reason"] == "no-unique-existing-family"


def test_omitted_conflicting_third_right_reports_per_class_support_only():
    from dataclasses import replace

    from primalscheme3.core.seq_functions import reverse_complement
    from primalscheme3.panel.allele_validation import numerical_dimer_score
    from primalscheme3.panel.coverage_types import ConfigurationLedger, OligoSite
    from primalscheme3.panel.coverage_variants import make_configuration

    seq = dna(seed=0)
    second = seq[180:200][:-1] + ("A" if seq[199] != "A" else "C")
    bad_sequence = "GCGCGCGCGCGCGCGCGCGC"
    cat, _, _ = fixture(
        seq=seq,
        rows=(
            seq,
            seq[:180] + second + seq[200:],
            seq[:180] + reverse_complement(bad_sequence) + seq[200:],
        ),
    )
    target = cat.targets[0]
    good2 = OligoSite(
        target.id,
        reverse_complement(second),
        "-",
        180,
        (180, 200),
        accepting_profile_ids=("normal",),
    )
    bad = OligoSite(
        target.id, bad_sequence, "-", 180, (180, 200), accepting_profile_ids=("normal",)
    )
    family = replace(
        cat.families[0],
        reverse_site_ids=cat.families[0].reverse_site_ids + (good2.id, bad.id),
    )
    cat = replace(cat, sites=cat.sites + (good2, bad), families=(family,))
    selected = make_configuration(
        cat,
        family.id,
        family.forward_site_ids,
        tuple(s for s in family.reverse_site_ids if s != bad.id),
    )
    ledger = ConfigurationLedger(cat.semantic_digest, (selected,))
    report = api().region_coverage(
        cat, ledger, (AlleleAssignment(selected.id, 0, "strict"),), target, (20, 180)
    )
    assert numerical_dimer_score(bad_sequence, bad_sequence).score <= -26
    omission = report["selected_families"][0]
    assert omission["family_id"] == family.id
    assert omission["omitted_eligible_site_ids"] == [bad.id]
    support = omission["omitted_eligible_sites"][0]["class_support"]
    confirmed = {s["allele_id"] for s in support if s["status"] == "confirmed"}
    assert len(confirmed) == 1
    uncovered = next(r for r in report["classes"] if r["allele_id"] in confirmed)
    assert uncovered["covered_bases"] == 0
    assert "not counterfactual" in omission["scope"]
