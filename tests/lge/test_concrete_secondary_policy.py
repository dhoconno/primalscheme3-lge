"""Four actual selected footprints; screening-only secondary certificates."""

from dataclasses import replace

import pytest
from test_allele_validation import api, dna, fixture, profile
from test_intended_site_policy import changed

NEW = "ordered-disjoint-concrete-designated-sites/v1"


def instance(position=199, cell=None):
    sequence = dna(seed=0)
    other = (
        changed(sequence, position)
        if cell is None
        else sequence[:position] + cell + sequence[position + 1 :]
    )
    return fixture(((0, 200), (300, 500)), seq=sequence, rows=(sequence, other))


def checker(cat, secondary=NEW, intended="concrete-designated-sites", cache=True):
    return api().SelectedSiteSpecificityChecker(
        cat,
        profile(secondary_product_policy=secondary, intended_product_policy=intended),
        cache=cache,
    )


@pytest.mark.parametrize("position", [0, 185, 199, 305, 499])
def test_complete_nonexact_secondary_is_separate_and_zero_credit(position):
    from primalscheme3.panel.allele_coverage import configuration_coverage

    cat, (a, b), _ = instance(position)
    before = [configuration_coverage(cat, c) for c in (a, b)]
    old = checker(cat, "ordered-disjoint-intended-sites").pair(a, b)
    assert old.rejected_products
    result = checker(cat).pair(a, b)
    assert not result.rejected_products
    assert len(result.allowed_secondary_products) == 2
    exact = [
        w
        for w in result.allowed_secondary_products
        if w["reason"] == "ordered-disjoint-intended-sites"
    ]
    new = [w for w in result.allowed_secondary_products if w["reason"] == NEW]
    assert len(exact) == len(new) == 1
    assert new[0]["coverage_credit"] == 0 and not new[0]["uncertain"]
    assert new[0]["classification"] == "nonexact-ordered-concrete-secondary-product"
    cert = new[0]["certificate"]
    assert (
        cert["left_configuration_id"] == a.id and cert["right_configuration_id"] == b.id
    )
    assert len(new[0]["site_ids"]) == 4
    for name in ("left_forward", "left_reverse", "right_forward", "right_reverse"):
        p = cert[name]
        hit = p["terminal_hit"]
        assert len(p["terminal_mismatch_positions"]) == hit["mismatches"] <= 1
        assert any(
            h.signature
            == (
                cat.targets[0].id,
                cert["row_id"],
                p["oligo"],
                hit["orientation"],
                hit["start"],
                hit["end"],
            )
            for h in checker(cat).hits(p["oligo"])
        )
    assert before == [configuration_coverage(cat, c) for c in (a, b)]
    assert checker(cat, cache=False).pair(a, b) == result
    assert checker(cat).pair(b, a) == result


@pytest.mark.parametrize("kind", ["two-terminal", "N", "IUPAC", "missing"])
def test_internal_partner_not_merely_concrete_geometry(kind):
    seq = dna(seed=0)
    if kind == "two-terminal":
        row = changed(changed(seq, 181), 182)
    else:
        row = seq[:199] + {"N": "N", "IUPAC": "R", "missing": ""}[kind] + seq[200:]
    if kind == "missing":
        row = tuple(seq[:199]) + ("",) + tuple(seq[200:])
    cat, (a, b), _ = fixture(((0, 200), (300, 500)), seq=seq, rows=(seq, row))
    result = checker(cat).pair(a, b)
    assert any(
        w["row_id"] == "row-1"
        and w["start"] == 0
        and w["end"] == (499 if kind == "missing" else 500)
        for w in result.rejected_products
    )
    assert not any(w.get("reason") == NEW for w in result.allowed_secondary_products)


def test_policy_identity_mapping_and_capability():
    from test_allele_options import kwargs

    from primalscheme3.core.config import Config
    from primalscheme3.panel.coverage_provenance import capabilities_document

    p = profile(secondary_product_policy=NEW)
    assert p.cache_key != profile().cache_key
    config = Config(**kwargs(secondary_product_policy=NEW))
    assert (
        api().AlleleConstraintProfile.from_config(config).secondary_product_policy
        == NEW
    )
    assert Config(**config.to_dict()).secondary_product_policy == NEW
    cap = capabilities_document()["alleleCoverage"]["secondaryProductPolicies"]
    assert cap["default"] == "ordered-disjoint-intended-sites"
    assert cap["policies"][NEW] == {"id": NEW, "coverageCredit": 0}


def test_fresh_publication_rejects_secondary_certificate_tamper(tmp_path):
    import json

    from primalscheme3.panel.allele_publication import (
        artifact_descriptor,
        audit_allele_stage,
        publish_allele_stage,
    )
    from primalscheme3.panel.coverage_types import AlleleAssignment

    cat, cs, ledger = instance()
    p = profile(
        secondary_product_policy=NEW,
        intended_product_policy="concrete-designated-sites",
    )
    assignments = tuple(AlleleAssignment(c.id, 0, "strict") for c in cs)
    stage = tmp_path / "stage"
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=assignments,
        profile=p,
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "test"},
    )
    report = audit_allele_stage(stage)
    assert report["valid"] and report["allowed_secondary_product_count"] == 2
    saved = json.loads((stage / "validation.json").read_text())
    new = next(w for w in saved["allowed_secondary_products"] if w["reason"] == NEW)
    new["certificate"]["left_reverse"]["expected_footprint"] = [181, 201]
    (stage / "validation.json").write_text(json.dumps(saved))
    manifest = json.loads((stage / "stage.json").read_text())
    manifest["artifacts"]["validation.json"] = artifact_descriptor(
        stage / "validation.json", stage
    )
    (stage / "stage.json").write_text(json.dumps(manifest))
    audited = audit_allele_stage(stage)
    assert not audited["valid"]
    assert any(
        v["reason"] == "saved-secondary-products-mismatch"
        for v in audited["violations"]
    )


@pytest.mark.parametrize("intended", ["exact-supported", "concrete-designated-sites"])
@pytest.mark.parametrize(
    "secondary",
    ["ordered-disjoint-intended-sites", "reject-secondary-products/v1", NEW],
)
def test_independent_policy_axes_preserve_unary_gate(intended, secondary):
    cat, (a, b), ledger = instance(199)
    p = profile(intended_product_policy=intended, secondary_product_policy=secondary)
    oracle = api().AlleleCompatibilityOracle(cat, ledger, p)
    unary = oracle.candidate_diagnostics(a.id)
    assert ("specificity" in unary["reasons"]) == (intended == "exact-supported")
    pair = oracle.pair_diagnostics(a.id, b.id)
    if intended == "exact-supported":
        assert "invalid-configuration" in pair["reasons"]
        assert pair["checks"]["specificity"] == "not-evaluated"
    else:
        assert ("specificity" in pair["reasons"]) == (secondary != NEW)


def test_overlap_third_target_wrong_occurrence_and_missing_partner():
    from allele_fixtures import target

    from primalscheme3.panel.allele_coverage import canonical_observations

    cat, (a, b), _ = instance(199)
    check = checker(cat, cache=False)
    old = checker(cat, "ordered-disjoint-intended-sites").pair(a, b)
    witness = next(w for w in old.rejected_products if w["row_id"] == "row-1")
    va, vb = check._view(a), check._view(b)
    assert check._concrete_secondary_certificate((va, vb), witness)
    assert (
        check._concrete_secondary_certificate(
            (va, replace(vb, forward_site_ids=())), witness
        )
        is None
    )
    assert check._concrete_secondary_certificate((va, vb, va), witness) is None
    shifted = witness | {"plus": witness["plus"] | {"start": 1}}
    assert check._concrete_secondary_certificate((va, vb), shifted) is None
    opposite = witness | {"minus": witness["minus"] | {"orientation": "+"}}
    assert check._concrete_secondary_certificate((va, vb), opposite) is None
    other = target(dna(seed=0), id="other")
    cross = replace(
        cat,
        targets=cat.targets + (other,),
        observations=cat.observations + canonical_observations(other),
    )
    assert (
        checker(cross, cache=False)._concrete_secondary_certificate(
            (va, vb), witness | {"target_id": other.id}
        )
        is None
    )
    seq = dna(seed=0)
    overlapping, (a, b), _ = fixture(
        ((0, 200), (190, 390)), seq=seq, rows=(seq, changed(seq, 199))
    )
    result = checker(overlapping).pair(a, b)
    assert result.rejected_products and not any(
        w["reason"] == NEW for w in result.allowed_secondary_products
    )


def test_half_open_inner_boundary_and_D():
    seq = dna(seed=0)
    cat, (a, b), _ = fixture(
        ((0, 200), (200, 400)), seq=seq, rows=(seq, changed(seq, 0))
    )
    result = checker(cat).pair(a, b)
    new = next(w for w in result.allowed_secondary_products if w["reason"] == NEW)
    assert (
        new["certificate"]["left_reverse"]["expected_footprint"][1]
        == new["certificate"]["right_forward"]["expected_footprint"][0]
    )
    check = api().SelectedSiteSpecificityChecker(
        cat,
        profile(
            399,
            secondary_product_policy=NEW,
            intended_product_policy="concrete-designated-sites",
        ),
    )
    assert not check.pair(a, b).allowed_secondary_products
    check = api().SelectedSiteSpecificityChecker(
        cat,
        profile(
            400,
            secondary_product_policy=NEW,
            intended_product_policy="concrete-designated-sites",
        ),
    )
    assert len(check.pair(a, b).allowed_secondary_products) == 2


def test_same_anchor_different_lengths_are_occurrences_not_envelopes():
    from primalscheme3.panel.coverage_types import OligoSite
    from primalscheme3.panel.coverage_variants import make_configuration

    cat, (a, b), _ = instance(199)
    seq = dna(seed=0)
    short = OligoSite(
        a.target_id, seq[1:20], "+", 20, (1, 20), accepting_profile_ids=("normal",)
    )
    af = cat.family_by_id[a.family_id]
    af = replace(af, forward_site_ids=af.forward_site_ids + (short.id,))
    cat = replace(
        cat,
        sites=cat.sites + (short,),
        families=tuple(af if f.id == a.family_id else f for f in cat.families),
    )
    a = make_configuration(cat, af.id, af.forward_site_ids, af.reverse_site_ids)
    result = checker(cat).pair(a, b)
    new = [w for w in result.allowed_secondary_products if w["reason"] == NEW]
    assert {
        tuple(w["certificate"]["left_forward"]["expected_footprint"]) for w in new
    } == {(0, 20), (1, 20)}
    assert all(
        w["plus"]["start"] == w["certificate"]["left_forward"]["expected_footprint"][0]
        for w in new
    )


def test_mixed_actual_offtarget_still_rejects_pair():
    seq = dna(seed=0)
    seq = seq[:600] + seq[:20] + seq[620:]
    cat, (a, b), _ = fixture(
        ((0, 200), (300, 500)), seq=seq, rows=(seq, changed(seq, 199))
    )
    result = checker(cat).pair(a, b)
    # The late F copies point inward to no R; introduce a shifted R downstream.
    seq = seq[:700] + seq[480:500] + seq[720:]
    cat, (a, b), _ = fixture(
        ((0, 200), (300, 500)), seq=seq, rows=(seq, changed(seq, 199))
    )
    result = checker(cat).pair(a, b)
    assert any(w["reason"] == NEW for w in result.allowed_secondary_products)
    assert any(w["end"] == 720 for w in result.rejected_products)


def test_authoritative_rows_history_and_coverage_are_independent(tmp_path):
    from primalscheme3.panel.coverage_history import SQLiteCoverageHistory
    from primalscheme3.panel.coverage_types import AlleleAssignment, ConfigurationLedger

    cat, cs, ledger = instance()
    assignments = tuple(AlleleAssignment(c.id, 0, "strict") for c in cs)
    p = profile(
        secondary_product_policy=NEW,
        intended_product_policy="concrete-designated-sites",
    )
    with SQLiteCoverageHistory(tmp_path / "history", run_id="concrete-secondary") as h:
        report = api().validate_allele_assignments(
            cat, assignments, ledger, p, authoritative_targets=cat.targets, history=h
        )
        measured = [
            e
            for e in h.evidence
            if any(
                w.get("reason") == NEW
                for w in e.values.get("allowed_secondary_products", ())
            )
        ]
        assert measured and any(len(e.entity_ids) == 2 for e in measured)
        assert all(e.dependency_key["profile"] == p.cache_key for e in measured)
    old = api().validate_allele_assignments(
        cat,
        assignments,
        ledger,
        replace(p, secondary_product_policy="ordered-disjoint-intended-sites"),
    )
    assert report["valid"] and not old["valid"]
    assert old["coverage"] == report["coverage"]
    seq = dna(seed=0)
    t = replace(
        cat.targets[0], rows=(tuple(seq), tuple(seq[:199]) + ("N",) + tuple(seq[200:]))
    )
    invalid = api().validate_allele_assignments(
        cat, assignments, ledger, p, authoritative_targets=(t,)
    )
    assert not invalid["valid"]
    assert not any(
        w.get("reason") == NEW for w in invalid["allowed_secondary_products"]
    )
    oracle = api().AlleleCompatibilityOracle(
        cat, ConfigurationLedger(cat.semantic_digest, ()), p
    )
    for c in cs:
        oracle.register(c)
    fresh = api().AlleleCompatibilityOracle(cat, ledger, p, cache=False)
    assert oracle.pair_diagnostics(cs[0].id, cs[1].id) == fresh.pair_diagnostics(
        cs[0].id, cs[1].id
    )


def test_no_discovery_identity_change_and_legacy_roundtrip():
    import json

    from test_allele_options import kwargs

    from primalscheme3.core.config import Config
    from primalscheme3.panel.allele_catalog_cache import discovery_settings
    from primalscheme3.panel.coverage_discovery import discovery_profiles

    old = Config(**kwargs())
    new = Config(**kwargs(secondary_product_policy=NEW))
    assert discovery_settings(old, discovery_profiles(old)) == discovery_settings(
        new, discovery_profiles(new)
    )
    serialized = old.to_dict()
    serialized.pop("secondary_product_policy")
    options = json.loads(serialized["allele_options_json"])
    options.pop("secondary_product_policy")
    serialized["allele_options_json"] = json.dumps(options)
    assert (
        Config(**serialized).secondary_product_policy
        == "ordered-disjoint-intended-sites"
    )
    with pytest.raises(ValueError):
        Config(secondary_product_policy=NEW)
    with pytest.raises(ValueError):
        Config(**kwargs(secondary_product_policy="unknown"))


def test_old_stage_missing_policy_and_secondary_count_keeps_exact_rule(tmp_path):
    import json

    from primalscheme3.panel.allele_publication import (
        artifact_descriptor,
        audit_allele_stage,
        publish_allele_stage,
    )
    from primalscheme3.panel.coverage_types import AlleleAssignment

    seq = dna(seed=0)
    cat, cs, ledger = fixture(((0, 200), (300, 500)), seq=seq)
    stage = tmp_path / "stage"
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=tuple(AlleleAssignment(c.id, 0, "strict") for c in cs),
        profile=profile(),
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "test"},
    )
    manifest = json.loads((stage / "stage.json").read_text())
    manifest["constraints"].pop("secondary_product_policy")
    saved = json.loads((stage / "validation.json").read_text())
    saved["profile"].pop("secondary_product_policy")
    saved.pop("allowed_secondary_product_count")
    (stage / "validation.json").write_text(json.dumps(saved))
    manifest["artifacts"]["validation.json"] = artifact_descriptor(
        stage / "validation.json", stage
    )
    (stage / "stage.json").write_text(json.dumps(manifest))
    report = audit_allele_stage(stage)
    assert (
        report["valid"]
        and report["profile"]["secondary_product_policy"]
        == "ordered-disjoint-intended-sites"
    )
    assert report["allowed_secondary_product_count"] == 1


def test_real_cli_capability_receipts_and_whole_bundle_audit(tmp_path):
    import json

    from typer.testing import CliRunner

    from primalscheme3.cli import app

    raw = tmp_path / "original.fa"
    raw.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "panel"
    argv = [
        "panel-create",
        "--msa",
        str(raw),
        "--output",
        str(out),
        "--selection-algorithm",
        "allele-coverage",
        "--amplicon-size",
        "200",
        "--amplicon-size-min",
        "150",
        "--amplicon-size-max",
        "250",
        "--secondary-product-policy",
        NEW,
        "--intended-product-policy",
        "concrete-designated-sites",
        "--offline-plots",
    ]
    result = CliRunner().invoke(app, argv)
    assert result.exit_code == 0, result.output + repr(result.exception)
    receipt = json.loads((out / "panel-provenance.json").read_text())
    assert receipt["resolvedOptions"]["secondary_product_policy"] == NEW
    audit = tmp_path / "audit"
    result = CliRunner().invoke(
        app, ["panel-audit", "--bundle", str(out), "--output", str(audit)]
    )
    assert result.exit_code == 0, result.output + repr(result.exception)
    report = json.loads((audit / "validation.json").read_text())
    assert report["valid"]
    assert report["stages"]["strict"]["profile"]["secondary_product_policy"] == NEW

    from primalscheme3.panel.allele_inspection import audit_allele_bundle
    from primalscheme3.panel.allele_publication import artifact_descriptor

    optimizer_path = out / "panel-optimizer.json"
    optimizer = json.loads(optimizer_path.read_text())
    optimizer["options"]["secondary_product_policy"] = "ordered-disjoint-intended-sites"
    optimizer_path.write_text(json.dumps(optimizer))
    for item in receipt["outputs"]:
        item.update(artifact_descriptor(out / item["path"], out))
    (out / "panel-provenance.json").write_text(json.dumps(receipt))
    tampered = audit_allele_bundle(out)
    assert not tampered["valid"]
    assert any(
        v["reason"] == "secondary-product-policy-mismatch"
        for v in tampered["violations"]
    )


def test_gaps_preserve_all_four_projected_occurrences():
    from primalscheme3.panel.allele_coverage import canonical_observations
    from primalscheme3.panel.coverage_variants import make_configuration

    cat, _, _ = instance()
    target = replace(
        cat.targets[0],
        rows=tuple(row[:5] + ("-",) + row[5:] for row in cat.targets[0].rows),
    )
    target = api()._rebuild_target(target)
    sites = {
        s.id: replace(s, alignment_anchor=s.alignment_anchor + 1) for s in cat.sites
    }
    families = tuple(
        replace(
            f,
            anchor_pair=tuple(x + 1 for x in f.anchor_pair),
            forward_site_ids=tuple(sites[s].id for s in f.forward_site_ids),
            reverse_site_ids=tuple(sites[s].id for s in f.reverse_site_ids),
        )
        for f in cat.families
    )
    cat = replace(
        cat,
        targets=(target,),
        observations=canonical_observations(target),
        sites=tuple(sites.values()),
        families=families,
    )
    a, b = tuple(
        make_configuration(cat, f.id, f.forward_site_ids, f.reverse_site_ids)
        for f in families
    )
    result = checker(cat).pair(a, b)
    cert = next(
        w["certificate"]
        for w in result.allowed_secondary_products
        if w["reason"] == NEW
    )
    assert not result.rejected_products
    assert [
        cert[name]["expected_footprint"]
        for name in ("left_forward", "left_reverse", "right_forward", "right_reverse")
    ] == [[0, 20], [180, 200], [300, 320], [480, 500]]
