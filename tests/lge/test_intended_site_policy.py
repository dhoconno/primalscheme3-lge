"""Opt-in screening certificates never broaden exact allele coverage."""

import json
from dataclasses import replace

import pytest
from allele_fixtures import target
from test_allele_validation import api, dna, fixture, profile

from primalscheme3.core.seq_functions import reverse_complement as rc
from primalscheme3.panel.allele_coverage import allele_summary, canonical_observations
from primalscheme3.panel.coverage_types import (
    AlleleAssignment,
    ConfigurationLedger,
    OligoSite,
)
from primalscheme3.panel.coverage_variants import make_configuration

POLICY = "concrete-designated-sites"


def changed(sequence, position):
    return (
        sequence[:position]
        + ("A" if sequence[position] != "A" else "C")
        + sequence[position + 1 :]
    )


def checker(cat, **kw):
    return api().SelectedSiteSpecificityChecker(
        cat, profile(intended_product_policy=POLICY, **kw)
    )


@pytest.mark.parametrize(
    "position,terminal", [(0, False), (5, True), (199, False), (185, True)]
)
def test_same_anchor_near_match_is_screening_only(position, terminal):
    seq = dna(seed=0)
    cat, (c,), ledger = fixture(seq=seq, rows=(seq, changed(seq, position)))
    old = api().SelectedSiteSpecificityChecker(cat, profile()).intrinsic(c)
    assert old.rejected_products
    result = checker(cat).intrinsic(c)
    assert not result.rejected_products
    assert not result.allowed_secondary_products
    assert len(result.allowed_intended_products) == 1
    w = result.allowed_intended_products[0]
    assert w["classification"] == "nonexact-designated-intended-product"
    assert w["coverage_credit"] == 0 and w["uncertain"] is False
    assert w["policy_id"] == POLICY + "/v1"
    cert = w["certificate"]
    assert cert["configuration_id"] == c.id
    assert cert["forward"]["expected_footprint"] == [0, 20]
    assert cert["reverse"]["expected_footprint"] == [180, 200]
    end = cert["forward" if position < 20 else "reverse"]
    assert len(end["mismatch_positions"]) == 1
    assert bool(end["terminal_mismatch_positions"]) == terminal
    assert bool(end["outside_terminal_mismatch_positions"]) != terminal
    assignments = (AlleleAssignment(c.id, 0, "strict"),)
    fresh = api().validate_allele_assignments(
        cat, assignments, ledger, profile(intended_product_policy=POLICY)
    )
    assert fresh["valid"], fresh["violations"]
    assert (
        fresh["coverage"]
        == allele_summary(cat, assignments, ledger, goal=0.95).to_dict()
    )
    assert sorted(v["covered_count"] for v in fresh["coverage"]["classes"]) == [0, 160]
    assert fresh["allowed_intended_product_count"] == 1
    assert (
        api()
        .AlleleConstraintProfile(
            **{
                k: v
                for k, v in profile().__dict__.items()
                if k != "intended_product_policy"
            }
        )
        .intended_product_policy
        == "exact-supported"
    )


def test_another_retained_exact_variant_can_certify_same_span_without_new_coverage():
    seq = dna(seed=0)
    other = changed(seq, 0)
    cat, (c,), _ = fixture(seq=seq, rows=(seq, other))
    f = OligoSite(
        c.target_id, other[:20], "+", 20, (0, 20), accepting_profile_ids=("normal",)
    )
    family = replace(
        cat.family_by_id[c.family_id], forward_site_ids=c.forward_site_ids + (f.id,)
    )
    cat = replace(cat, sites=cat.sites + (f,), families=(family,))
    full = make_configuration(
        cat, family.id, family.forward_site_ids, family.reverse_site_ids
    )
    result = checker(cat).intrinsic(full)
    assert not result.rejected_products
    assert len(result.allowed_intended_products) == 2
    assert len(checker(cat).intended(full)[0]) == 2
    assert all(w["coverage_credit"] == 0 for w in result.allowed_intended_products)


@pytest.mark.parametrize("cell", ["N", "R", ""])
def test_ambiguous_or_missing_full_footprint_blocks_even_with_concrete_terminal(cell):
    seq = dna()
    row = (cell,) + tuple(seq[1:])
    cat, (c,), _ = fixture(seq=seq, rows=(seq, row))
    result = checker(cat).intrinsic(c)
    assert result.rejected_products and any(
        w["uncertain"] for w in result.rejected_products
    )
    assert not result.allowed_intended_products


def test_shifted_exact_locus_and_wrong_orientation_remain_blocked():
    seq = dna()
    seq = seq[:300] + seq[:20] + seq[320:480] + seq[180:200] + seq[500:]
    cat, (c,), _ = fixture(seq=seq)
    result = checker(cat).intrinsic(c)
    assert any(w["start"] == 300 and w["end"] == 500 for w in result.rejected_products)
    seq = dna()
    seq = seq[:300] + rc(seq[:20]) + seq[320:]
    cat, (c,), _ = fixture(seq=seq)
    result = checker(cat).intrinsic(c)
    assert any(w["start"] == 0 and w["end"] == 320 for w in result.rejected_products)
    assert not result.allowed_intended_products


def test_no_other_target_or_third_configuration_rescue_and_secondary_unchanged():
    seq = dna()
    cat, (a, b, c), _ = fixture(
        ((0, 200), (300, 500), (150, 350)), seq=seq, rows=(seq, changed(seq, 5))
    )
    result = checker(cat).pair(a, b)
    assert (
        result.rejected_products
    )  # F(A) near-match + R(B) lacks EXACT secondary certificate
    assert len(result.allowed_secondary_products) == 1
    assert all(
        w["reason"] == "ordered-disjoint-intended-sites"
        for w in result.allowed_secondary_products
    )
    assert any(
        w["start"] == 0 and w["end"] == 500 and w["row_id"] == "row-1"
        for w in result.rejected_products
    )
    assert (
        checker(cat, secondary_product_policy="reject-secondary-products/v1")
        .pair(a, b)
        .rejected_products
    )
    reduced = replace(
        cat, families=tuple(f for f in cat.families if f.id != c.family_id)
    )
    assert checker(reduced).pair(a, b) == result
    other = target(changed(seq, 5), id="other")
    both = replace(
        cat,
        targets=cat.targets + (other,),
        observations=cat.observations + canonical_observations(other),
    )
    rejected = checker(both).intrinsic(a).rejected_products
    assert any(w["target_id"] == other.id for w in rejected)


def test_profile_identity_and_invalid_policy():
    assert profile().to_dict()["intended_product_policy"] == "exact-supported"
    assert profile().cache_key != profile(intended_product_policy=POLICY).cache_key
    with pytest.raises(ValueError, match="intended"):
        profile(intended_product_policy="anything")


def test_cli_option_capability_roundtrip_and_legacy_rejection(tmp_path):
    from test_allele_options import invoke, kwargs

    from primalscheme3.core.config import Config
    from primalscheme3.panel.allele_options import AlleleOptions
    from primalscheme3.panel.coverage_provenance import capabilities_document

    result, run = invoke(tmp_path, ("--intended-product-policy", POLICY))
    assert result.exit_code == 0, result.output
    config = run.call_args.kwargs["config"]
    assert AlleleOptions.from_config(config).intended_product_policy == POLICY
    assert Config(**config.to_dict()).intended_product_policy == POLICY
    assert (
        api().AlleleConstraintProfile.from_config(config).intended_product_policy
        == POLICY
    )
    cap = capabilities_document()["alleleCoverage"]["intendedProductPolicies"]
    assert cap["default"] == "exact-supported"
    assert cap["policies"][POLICY] == {"id": POLICY + "/v1", "coverageCredit": 0}
    with pytest.raises(ValueError):
        Config(intended_product_policy=POLICY)
    with pytest.raises(ValueError):
        Config(**kwargs(intended_product_policy="unknown"))
    old = Config(**kwargs()).to_dict()
    old.pop("intended_product_policy")
    saved = json.loads(old["allele_options_json"])
    saved.pop("intended_product_policy")
    old["allele_options_json"] = json.dumps(saved)
    assert Config(**old).intended_product_policy == "exact-supported"


def test_publication_fresh_certificate_audit_and_tamper(tmp_path):
    from primalscheme3.panel.allele_publication import (
        artifact_descriptor,
        audit_allele_stage,
        publish_allele_stage,
    )

    seq = dna(seed=0)
    cat, (c,), ledger = fixture(seq=seq, rows=(seq, changed(seq, 0)))
    assignments = (AlleleAssignment(c.id, 0, "strict"),)
    stage = tmp_path / "stage"
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=assignments,
        profile=profile(intended_product_policy=POLICY),
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "test"},
    )
    report = audit_allele_stage(stage)
    assert report["valid"] and report["allowed_intended_product_count"] == 1
    saved = json.loads((stage / "validation.json").read_text())
    saved["allowed_intended_products"][0]["certificate"]["forward"][
        "expected_footprint"
    ] = [1, 21]
    (stage / "validation.json").write_text(json.dumps(saved))
    manifest = json.loads((stage / "stage.json").read_text())
    manifest["artifacts"]["validation.json"] = artifact_descriptor(
        stage / "validation.json", stage
    )
    (stage / "stage.json").write_text(json.dumps(manifest))
    report = audit_allele_stage(stage)
    assert not report["valid"]
    assert any(
        v["reason"] == "saved-intended-products-mismatch" for v in report["violations"]
    )


def test_selected_side_membership_is_required_even_for_shared_physical_species():
    seq = dna()
    cat, (c,), _ = fixture(seq=seq, rows=(seq, changed(seq, 0)))
    check = checker(cat)
    view = check._view(c)
    old = api().SelectedSiteSpecificityChecker(cat, profile()).intrinsic(c)
    witness = old.rejected_products[0]
    # The complete-list association is insufficient if the selected F site is
    # only on the wrong side of a malformed input view.
    wrong = replace(view, forward_site_ids=(), reverse_site_ids=view.selected_site_ids)
    assert check._designated_certificate((wrong,), witness) is None


def test_cached_uncached_registration_history_and_authoritative_rows(tmp_path):
    from primalscheme3.panel.allele_inspection import query_allele_history
    from primalscheme3.panel.allele_publication import _write
    from primalscheme3.panel.coverage_history import SQLiteCoverageHistory

    seq = dna(seed=0)
    cat, (c,), ledger = fixture(seq=seq, rows=(seq, changed(seq, 0)))
    p = profile(intended_product_policy=POLICY)
    cached = api().AlleleCompatibilityOracle(
        cat, ConfigurationLedger(cat.semantic_digest, ()), p
    )
    cached.register(c)
    fresh = api().AlleleCompatibilityOracle(cat, ledger, p, cache=False)
    assert cached.candidate_diagnostics(c.id) == fresh.candidate_diagnostics(c.id)
    assert (
        cached.cache_identity
        != api().AlleleCompatibilityOracle(cat, ledger, profile()).cache_identity
    )
    bundle = tmp_path / "history-bundle"
    bundle.mkdir()
    _write(bundle / "catalog.json.gz", cat.to_dict())
    _write(bundle / "configuration-ledger.json.gz", ledger.to_dict())
    _write(
        bundle / "panel-optimizer.json",
        {
            "catalogPath": "catalog.json.gz",
            "history": {"path": "history/history.sqlite"},
        },
    )
    assignments = (AlleleAssignment(c.id, 0, "strict"),)
    with SQLiteCoverageHistory(bundle / "history", run_id="intended-policy") as h:
        result = api().validate_allele_assignments(
            cat, assignments, ledger, p, history=h, authoritative_targets=cat.targets
        )
        assert result["valid"] and result["allowed_intended_product_count"] == 1
    queried = query_allele_history(bundle, entity=c.id, stage="strict", limit=100)
    evidence = [
        r for r in queried["evidence"] if r["values"].get("allowed_intended_products")
    ]
    assert evidence
    certificate = evidence[0]["values"]["allowed_intended_products"][0]["certificate"]
    assert certificate["forward_site_id"] in cat.site_by_id
    assert certificate["row_id"] in cat.targets[0].row_ids
    assert any(
        a["thresholds"]["intended_product_policy"] == POLICY
        for a in queried["assessments"]
    )
    # Fresh validation never trusts cached certificates: authoritative N outside
    # terminal17 must replace the old concrete row's tolerance with a block.
    altered = replace(cat.targets[0], rows=(tuple(seq), ("N",) + tuple(seq[1:])))
    rejected = api().validate_allele_assignments(
        cat, assignments, ledger, p, authoritative_targets=(altered,)
    )
    assert not rejected["valid"] and rejected["uncertainty_blocks"]
    assert rejected["allowed_intended_products"] == []


def test_gapped_projection_missing_elsewhere_and_product_bound():
    seq = dna()
    other = changed(seq, 0)
    cat, (c,), _ = fixture(seq=seq, rows=(seq, other))
    t = replace(
        cat.targets[0],
        rows=tuple(row[:5] + ("-",) + row[5:] for row in cat.targets[0].rows),
    )
    # Rebuild reference mapping independently, preserving declared occurrence.
    t = api()._rebuild_target(t)
    sites = tuple(
        replace(s, alignment_anchor=s.alignment_anchor + 1) for s in cat.sites
    )
    family = replace(
        cat.families[0],
        anchor_pair=(21, 181),
        forward_site_ids=tuple(s.id for s in sites if s.strand == "+"),
        reverse_site_ids=tuple(s.id for s in sites if s.strand == "-"),
    )
    cat = replace(
        cat,
        targets=(t,),
        observations=canonical_observations(t),
        sites=sites,
        families=(family,),
    )
    c = make_configuration(
        cat, family.id, family.forward_site_ids, family.reverse_site_ids
    )
    result = checker(cat).intrinsic(c)
    assert len(result.allowed_intended_products) == 1 and not result.rejected_products
    assert result.allowed_intended_products[0]["certificate"]["forward"][
        "expected_footprint"
    ] == [0, 20]
    assert (
        api()
        .SelectedSiteSpecificityChecker(
            cat, profile(199, intended_product_policy=POLICY)
        )
        .intrinsic(c)
        .allowed_intended_products
        == ()
    )
    assert (
        len(
            api()
            .SelectedSiteSpecificityChecker(
                cat, profile(200, intended_product_policy=POLICY)
            )
            .intrinsic(c)
            .allowed_intended_products
        )
        == 1
    )
    # Missing sequence beyond the two concrete footprints is not an uncertainty
    # shortcut that rejects an otherwise in-place concrete product.
    t = replace(t, rows=(t.rows[0], t.rows[1][:600] + ("",) + t.rows[1][601:]))
    cat = replace(cat, targets=(t,), observations=canonical_observations(t))
    assert checker(cat).intrinsic(c).allowed_intended_products


def test_new_policy_does_not_change_discovery_scientific_settings():
    from test_allele_options import kwargs

    from primalscheme3.core.config import Config
    from primalscheme3.panel.allele_catalog_cache import discovery_settings
    from primalscheme3.panel.coverage_discovery import discovery_profiles

    old = Config(**kwargs())
    new = Config(**kwargs(intended_product_policy=POLICY))
    assert discovery_settings(old, discovery_profiles(old)) == discovery_settings(
        new, discovery_profiles(new)
    )


def test_raw_fasta_authority_roundtrip_keeps_coverage_exact(tmp_path):
    from primalscheme3.panel.allele_pipeline import _authoritative_targets
    from primalscheme3.panel.allele_publication import (
        audit_allele_stage,
        publish_allele_stage,
    )
    from primalscheme3.panel.coverage_types import CandidateFamily

    seq = dna(seed=0)
    other = changed(seq, 0)
    raw = tmp_path / "original.fa"
    raw.write_text(">exact\n" + seq + "\n>near\n" + other + "\n")
    inputs = [{"storedPath": "original.fa", "sourceIndex": 0}]
    original = _authoritative_targets(tmp_path, inputs)
    cat, _, _ = fixture(seq=seq, rows=(seq, other))
    sites = tuple(replace(s, target_id=original[0].id) for s in cat.sites)
    family = CandidateFamily(
        original[0].id,
        (20, 180),
        tuple(s.id for s in sites if s.strand == "+"),
        tuple(s.id for s in sites if s.strand == "-"),
    )
    cat = replace(
        cat,
        targets=original,
        observations=canonical_observations(original[0]),
        sites=sites,
        families=(family,),
    )
    c = make_configuration(
        cat, family.id, family.forward_site_ids, family.reverse_site_ids
    )
    ledger = ConfigurationLedger(cat.semantic_digest, (c,))
    assignments = (AlleleAssignment(c.id, 0, "strict"),)
    stage = tmp_path / "stage"
    p = profile(intended_product_policy=POLICY)
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=assignments,
        profile=p,
        authoritative_targets=_authoritative_targets(tmp_path, inputs),
        references={original[0].id: "exact"},
    )
    report = audit_allele_stage(stage)
    assert report["valid"] and report["allowed_intended_product_count"] == 1
    assert sorted(r["covered_bases"] for r in report["literal_coverage"]) == [0, 160]
    # A changed raw non-reference row cannot inherit the prior concrete certificate.
    raw.write_text(">exact\n" + seq + "\n>near\nN" + other[1:] + "\n")
    invalid = api().validate_allele_assignments(
        cat,
        assignments,
        ledger,
        p,
        authoritative_targets=_authoritative_targets(tmp_path, inputs),
        selected_sites=report["selected_sites"],
    )
    assert not invalid["valid"]
    assert invalid["allowed_intended_products"] == []


def test_real_cli_policy_receipts_and_whole_bundle_audit(tmp_path):
    from typer.testing import CliRunner

    from primalscheme3.cli import app

    path = tmp_path / "short.fa"
    path.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "panel"
    run = CliRunner().invoke(
        app,
        [
            "panel-create",
            "--msa",
            str(path),
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
            "--intended-product-policy",
            POLICY,
            "--offline-plots",
        ],
    )
    assert run.exit_code == 0, run.output + repr(run.exception)
    optimizer = json.loads((out / "panel-optimizer.json").read_text())
    assert optimizer["profile"]["intended_product_policy"] == POLICY
    assert optimizer["options"]["intended_product_policy"] == POLICY
    receipt = json.loads((out / "panel-provenance.json").read_text())
    assert (
        receipt["exitStatus"] == 0
        and receipt["resolvedOptions"]["intended_product_policy"] == POLICY
    )
    audit = tmp_path / "audit"
    run = CliRunner().invoke(
        app, ["panel-audit", "--bundle", str(out), "--output", str(audit)]
    )
    assert run.exit_code == 0, run.output + repr(run.exception)
    report = json.loads((audit / "validation.json").read_text())
    assert report["valid"] and report["raw_inputs_reparsed"]
    assert report["stages"]["strict"]["allowed_intended_product_count"] == 0
    assert report["stages"]["strict"]["allowed_intended_products"] == []
    assert json.loads((audit / "provenance.json").read_text())["exitStatus"] == 0
