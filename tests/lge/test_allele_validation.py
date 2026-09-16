"""Selected-site specificity and independent validation scientific regressions."""

import json
import random
from copy import deepcopy
from dataclasses import replace

import pytest
from allele_fixtures import target
from primalschemers import do_pool_interact

from primalscheme3.core.config import Config
from primalscheme3.core.seq_functions import reverse_complement as rc
from primalscheme3.panel.allele_coverage import canonical_observations
from primalscheme3.panel.coverage_types import (
    AlleleAssignment,
    CandidateFamily,
    ConfigurationLedger,
    OligoSite,
    VariantCatalog,
)
from primalscheme3.panel.coverage_variants import make_configuration


def api():
    from primalscheme3.panel import allele_validation

    return allele_validation


def dna(n=700, seed=91):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def chemistry():
    from primalscheme3.core.variant_thermo import _FIELDS

    config = Config()
    values = {k: getattr(config, k) for k in _FIELDS}
    values.update(
        primer_size_min=17,
        primer_size_max=36,
        primer_gc_min=0,
        primer_gc_max=100,
        primer_tm_min=-100,
        primer_tm_max=200,
        primer_hairpin_th_max=200,
        primer_homopolymer_max=40,
    )
    return values


def fixture(spans=((0, 200),), seq=None, rows=None):
    seq = seq or dna()
    t = target(*(rows or (seq,)))
    sites, families = [], []
    for start, end in spans:
        f = OligoSite(
            t.id,
            seq[start : start + 20],
            "+",
            start + 20,
            (start, start + 20),
            accepting_profile_ids=("normal",),
        )
        r = OligoSite(
            t.id,
            rc(seq[end - 20 : end]),
            "-",
            end - 20,
            (end - 20, end),
            accepting_profile_ids=("normal",),
        )
        sites.extend((f, r))
        families.append(CandidateFamily(t.id, (start + 20, end - 20), (f.id,), (r.id,)))
    sites = tuple({s.id: s for s in sites}.values())
    cat = VariantCatalog(
        (t,),
        canonical_observations(t),
        sites,
        tuple(families),
        json.dumps(
            {
                "profiles": {"normal": chemistry()},
                "amplicon_size_min": 40,
                "amplicon_size_max": 1000,
            }
        ),
    )
    cs = tuple(
        make_configuration(cat, f.id, f.forward_site_ids, f.reverse_site_ids)
        for f in families
    )
    return cat, cs, ConfigurationLedger(cat.semantic_digest, cs)


def profile(bound=2000, **kw):
    return api().AlleleConstraintProfile(
        amplicon_size_min=40,
        amplicon_size_max=1000,
        n_pools=2,
        mismatch_product_size=bound,
        **kw,
    )


def test_specificity_pair_local_certificates_and_secondary_policy():
    a = api()
    cat, (x, y, z), _ = fixture(((0, 200), (300, 500), (150, 350)))
    check = a.SelectedSiteSpecificityChecker(cat, profile())
    result = check.pair(x, y)
    assert not result.rejected_products
    assert [
        (p["start"], p["end"], p["coverage_credit"])
        for p in result.allowed_secondary_products
    ] == [(0, 500, 0)]
    assert check.pair(x, y) == check.pair(y, x)
    assert check.pair(x, z).rejected_products
    small = replace(
        cat,
        families=tuple(f for f in cat.families if f.id in (x.family_id, y.family_id)),
    )
    assert a.SelectedSiteSpecificityChecker(small, profile()).pair(x, y) == result
    strict = a.SelectedSiteSpecificityChecker(
        cat, profile(secondary_product_policy="reject-secondary-products/v1")
    )
    assert strict.pair(x, y).rejected_products


def test_inclusive_bound_ambiguity_and_unavailable_full_footprints():
    a = api()
    seq = dna()
    seq = seq[:300] + seq[:20] + seq[320:480] + seq[180:200] + seq[500:]
    cat, (c,), _ = fixture(seq=seq)
    assert (
        a.SelectedSiteSpecificityChecker(cat, profile(200))
        .intrinsic(c)
        .rejected_products
    )
    assert (
        not a.SelectedSiteSpecificityChecker(cat, profile(199))
        .intrinsic(c)
        .rejected_products
    )
    original = dna()
    uncertain = original[:60] + "N" + original[1:20] + original[80:]
    cat, (c,), _ = fixture(seq=original, rows=(original, uncertain))
    result = a.SelectedSiteSpecificityChecker(cat, profile()).intrinsic(c)
    assert any(p["uncertain"] and p["start"] == 60 for p in result.rejected_products)
    row = ("",) + tuple(original[1:])
    cat, (c,), _ = fixture(seq=original, rows=(original, row))
    result = a.SelectedSiteSpecificityChecker(cat, profile()).intrinsic(c)
    assert any(
        p["plus"]["start"] == -1 and p["uncertain"] for p in result.rejected_products
    )


def test_no_seed_invented_beyond_terminal_padding_and_positive_bound():
    a = api()
    original = dna()
    row = ("",) * 100 + tuple(original[100:])
    cat, (c,), _ = fixture(seq=original, rows=(original, row))
    hits = a.SelectedSiteSpecificityChecker(cat, profile()).hits(original[:20])
    assert not any(h.row_id == "row-1" for h in hits)
    assert (
        profile().to_dict()["uncertainty_scope"]
        == "supplied-row-seeds-and-seeded-full-footprints/v1"
    )
    for bound in (0, -1):
        with pytest.raises(ValueError):
            profile(bound)


def test_numeric_native_dimer_both_directions_and_equality():
    a = api()
    sequences = (
        "CAGTTGTCTTGTTCGCGAGTCGT",
        "ACGCGCGCGCTATATATACC",
        "GCGCGCGCGCGCGCGCGCGC",
    )
    for s in sequences:
        for t in sequences:
            score = a.numerical_dimer_score(s, t)
            assert len(score.orientation_scores) == 2
            for threshold in (-26, -28, score.score, score.score - 0.01):
                assert (score.score <= threshold) == bool(
                    do_pool_interact([s.encode()], [t.encode()], threshold)
                )
    assert a.StagePolicy().rejects(-26)
    assert not a.StagePolicy().rejects(-25.99)
    with pytest.raises(ValueError):
        a.StagePolicy(strict_cutoff=-25)


def test_pool_exposure_counts_physical_edges_not_site_copies():
    a = api()
    cat, (c,), ledger = fixture()
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=8,
        max_incident_species=4,
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile(), policy)
    result = oracle.pool_diagnostics((c.id, c.id))
    assert (
        len(result["dimer_edges"]) == 3
    )  # two self edges and one unordered cross edge
    assert len({tuple(e["sequences"]) for e in result["dimer_edges"]}) == 3


def test_fresh_validation_rebuilds_support_coverage_and_exact_selected_output():
    a = api()
    cat, (c,), ledger = fixture()
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=8,
        max_incident_species=4,
    )
    assignments = (AlleleAssignment(c.id, 0, "salvage-1"),)
    result = a.validate_allele_assignments(
        cat, assignments, ledger, profile(), policy=policy
    )
    assert result["valid"], result["violations"]
    assert result["coverage"]["classes"][0]["covered_count"] == 160
    assert result["coverage"]["classes"][0]["observed_count"] == 700
    assert a.validate_allele_assignments(cat, (), ledger, profile())["valid"]
    assert not a.validate_allele_assignments(
        cat,
        assignments,
        ledger,
        profile(),
        policy=policy,
        selected_sites=result["selected_sites"][:-1],
    )["valid"]
    wrong_summary = deepcopy(result["coverage"])
    wrong_summary["classes"][0]["observed_count"] = 699
    assert not a.validate_allele_assignments(
        cat,
        assignments,
        ledger,
        profile(),
        policy=policy,
        expected_summary=wrong_summary,
    )["valid"]
    stale = replace(c, supported_products=())
    stale_ledger = ConfigurationLedger(cat.semantic_digest, (stale,))
    assert not a.validate_allele_assignments(
        cat, assignments, stale_ledger, profile(), policy=policy
    )["valid"]
    outside = replace(c, full_interval=(-1, 200))
    assert not a.validate_allele_assignments(
        cat,
        assignments,
        ConfigurationLedger(cat.semantic_digest, (outside,)),
        profile(),
        policy=policy,
    )["valid"]


def test_complete_profile_membership_is_remeasured_not_rectangular_union():
    a = api()
    cat, (c,), _ = fixture()
    resolved = json.loads(cat.resolved_config_json)
    resolved["profiles"]["normal"]["primer_gc_max"] = 100
    resolved["profiles"]["high-gc"] = dict(chemistry(), primer_gc_min=100)
    cat = replace(cat, resolved_config_json=json.dumps(resolved))
    ledger = ConfigurationLedger(cat.semantic_digest, (c,))
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=8,
        max_incident_species=4,
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile(), policy)
    assert oracle.candidate_valid(c.id), oracle.candidate_diagnostics(c.id)
    badsites = tuple(replace(s, accepting_profile_ids=("high-gc",)) for s in cat.sites)
    badcat = replace(cat, sites=badsites)
    result = a.validate_allele_assignments(
        badcat,
        (AlleleAssignment(c.id, 0, "salvage-1"),),
        ConfigurationLedger(badcat.semantic_digest, (c,)),
        profile(),
        policy=policy,
    )
    assert not result["valid"]
    assert any(v["reason"] == "false-profile-membership" for v in result["violations"])


def test_fresh_validation_does_not_trust_oracle_cache_or_observations():
    a = api()
    cat, (x, y), ledger = fixture(((0, 200), (150, 350)))
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=100,
        max_incident_species=100,
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile(), policy)
    assert oracle.conflict(x.id, y.id)
    oracle._pairs[tuple(sorted((x.id, y.id)))] = a._DiagnosticCacheEntry.capture({"reasons": []})
    assert not oracle.conflict(x.id, y.id)
    result = a.validate_allele_assignments(
        cat,
        (
            AlleleAssignment(x.id, 0, "salvage-1"),
            AlleleAssignment(y.id, 0, "salvage-1"),
        ),
        ledger,
        profile(),
        policy=policy,
    )
    assert not result["valid"] and any(
        v["reason"] == "specificity" for v in result["violations"]
    )
    badcat = replace(cat, observations=())
    result = a.validate_allele_assignments(
        badcat, (), ConfigurationLedger(badcat.semantic_digest, ()), profile()
    )
    assert not result["valid"] and any(
        v["reason"] == "observation-mismatch" for v in result["violations"]
    )


def test_actual_pool_aggregate_exceeds_individually_passing_exposure():
    a = api()
    cat, cs, ledger = fixture(((0, 200), (250, 450), (480, 680)), seq=dna(seed=13))
    generous = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=100,
        max_incident_species=100,
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile(), generous)
    singles = [oracle.pool_diagnostics((c.id,)) for c in cs]
    all_edges = oracle.pool_diagnostics(tuple(c.id for c in cs))
    bound = max(x["violating_edge_count"] for x in singles)
    assert all_edges["violating_edge_count"] > bound
    restricted = a.AlleleCompatibilityOracle(
        cat, ledger, profile(), replace(generous, max_violating_edges=bound)
    )
    assert all(restricted.pool_valid((c.id,)) for c in cs)
    assert not restricted.pool_valid(tuple(c.id for c in cs))


def test_exact_variant_removal_and_shared_sequence_sites():
    a = api()
    seq = dna()
    seq = seq[:300] + seq[:20] + seq[320:]
    cat, (x, y), _ = fixture(((0, 200), (300, 500)), seq=seq)
    checker = a.SelectedSiteSpecificityChecker(cat, profile(200))
    sx, _ = checker.intended(x)
    sy, _ = checker.intended(y)
    assert len(sx) == len(sy) == 1
    assert next(iter(sx))[0][2] == next(iter(sy))[0][2]
    assert next(iter(sx))[0][4] != next(iter(sy))[0][4]
    assert not checker.pair(x, y).rejected_products
    # A second reverse variant supports only a distinct allele at the same site.
    original = dna()
    changed = original[:180] + "A" * 20 + original[200:]
    cat, (c,), _ = fixture(seq=original, rows=(original, changed))
    alt = OligoSite(
        c.target_id, "T" * 20, "-", 180, (180, 200), accepting_profile_ids=("normal",)
    )
    family = replace(
        cat.family_by_id[c.family_id], reverse_site_ids=c.reverse_site_ids + (alt.id,)
    )
    cat = replace(cat, sites=cat.sites + (alt,), families=(family,))
    full = make_configuration(
        cat, family.id, c.forward_site_ids, family.reverse_site_ids
    )
    sub = make_configuration(cat, family.id, c.forward_site_ids, c.reverse_site_ids)
    checker = a.SelectedSiteSpecificityChecker(cat, profile())
    assert len(checker.intended(full)[0]) == 2
    assert len(checker.intended(sub)[0]) == 1


def test_no_cross_row_product_and_unsupported_short_seed_rejects():
    a = api()
    original = dna()
    cat, (c,), _ = fixture(seq=original)
    rows = (
        tuple(original[:20]) + ("",) * 680,
        ("",) * 180 + tuple(original[180:200]) + ("",) * 500,
    )
    t = replace(cat.targets[0], rows=rows, row_ids=("f-only", "r-only"))
    small = replace(cat, targets=(t,), observations=canonical_observations(t))
    checker = a.SelectedSiteSpecificityChecker(small, profile())
    assert not checker.intrinsic(c).rejected_products
    p = profile(mismatch_kmersize=21)
    result = a.AlleleCompatibilityOracle(
        cat, ConfigurationLedger(cat.semantic_digest, (c,)), p
    ).candidate_diagnostics(c.id)
    assert "oligo-shorter-than-terminal-kmer" in result["reasons"]
    assert result["checks"]["specificity"] == "not-evaluated"


def test_native_complete_profile_and_effective_tm_tolerance():
    a = api()
    from primalscheme3.core.variant_thermo import _FIELDS

    config = Config()
    values = {k: getattr(config, k) for k in _FIELDS}
    measurements = a.chemistry_measurements("CAGTTGTCTTGTTCGCGAGTCGT", values)
    assert 62.5 < measurements["tm"]["value"] < 64.5
    assert all(check["outcome"] == "pass" for check in measurements.values())


def test_history_records_short_circuit_and_authoritative_output_scope():
    from primalscheme3.panel.coverage_history import CoverageHistory

    a = api()
    cat, (c,), ledger = fixture()
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=8,
        max_incident_species=4,
    )
    assignments = (AlleleAssignment(c.id, 0, "salvage-1"),)
    h = CoverageHistory(run_id="validation")
    report = a.validate_allele_assignments(
        cat,
        assignments,
        ledger,
        profile(mismatch_kmersize=21),
        policy=policy,
        history=h,
    )
    assert not report["valid"]
    assert any(
        x.check_name == "specificity" and x.outcome == "not-evaluated"
        for x in h.assessments
    )
    report = a.validate_allele_assignments(
        cat, assignments, ledger, profile(), policy=policy
    )
    full = a.validate_allele_assignments(
        cat,
        assignments,
        ledger,
        profile(),
        policy=policy,
        authoritative_targets=cat.targets,
        selected_sites=report["selected_sites"],
        expected_summary=report["coverage"],
    )
    assert (
        full["valid"]
        and full["validation_scope"]
        == "external-authoritative-inputs-and-selected-artifacts"
    )


def test_third_complete_certificate_cannot_rescue_pair_and_all_orientation_edges():
    a = api()
    cat, (x, y, z), ledger = fixture(((0, 200), (150, 350), (150, 200)))
    checker = a.SelectedSiteSpecificityChecker(cat, profile())
    with_c = checker.pair(x, y)
    small = replace(cat, families=tuple(f for f in cat.families if f.id != z.family_id))
    assert (
        with_c.rejected_products
        and a.SelectedSiteSpecificityChecker(small, profile()).pair(x, y) == with_c
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile())
    edges = oracle.pool_diagnostics((x.id, y.id))["dimer_edges"]
    assert len(edges) == 10
    strandpairs = set()
    for edge in edges:
        owned = [
            cat.site_by_id[sid].strand
            for owners in edge["owners"].values()
            for _, sid in owners
        ]
        strandpairs.add(tuple(sorted(set(owned))))
    assert {("+",), ("-",), ("+", "-")} <= strandpairs
    assert all(edge["id"] and edge["site_ids"] for edge in edges)
    assert all(w["id"] and w["site_ids"] for w in with_c.rejected_products)


def test_raw_score_cache_is_reusable_across_stage_cutoffs():
    a = api()
    cat, (c,), ledger = fixture()
    shared = {}
    strict = a.AlleleCompatibilityOracle(cat, ledger, profile(), score_cache=shared)
    strict_result = strict.pool_diagnostics((c.id,))
    salvage = a.AlleleCompatibilityOracle(
        cat,
        ledger,
        profile(),
        a.StagePolicy(
            stage_id="salvage-1",
            active_cutoff=-1000,
            max_violating_edges=8,
            max_incident_species=4,
        ),
        score_cache=shared,
    )
    assert shared
    assert (
        strict_result["dimer_edges"] == salvage.pool_diagnostics((c.id,))["dimer_edges"]
    )
    assert strict.cache_identity != salvage.cache_identity


def test_real_normal_only_chemistry_pass_is_preserved():
    a = api()
    from primalscheme3.core.variant_thermo import _FIELDS

    c = Config()
    normal = {k: getattr(c, k) for k in _FIELDS}
    high = dict(
        normal,
        primer_size_min=17,
        primer_size_max=30,
        primer_gc_min=40,
        primer_gc_max=65,
    )
    oligo = "TAAATCAGAAATGGAACAAAGCACCCT"
    assert all(
        x["outcome"] == "pass" for x in a.chemistry_measurements(oligo, normal).values()
    )
    assert a.chemistry_measurements(oligo, high)["gc"]["outcome"] == "fail"


def test_fresh_validation_measures_each_selected_site_once_before_pair_checks(
    monkeypatch,
):
    a = api()
    cat, (x, y), ledger = fixture(((0, 200), (300, 500)))
    measured = []
    real = a.chemistry_measurements

    def counted(sequence, settings):
        measured.append(sequence)
        return real(sequence, settings)

    monkeypatch.setattr(a, "chemistry_measurements", counted)
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=100,
        max_incident_species=100,
    )
    report = a.validate_allele_assignments(
        cat,
        (
            AlleleAssignment(x.id, 0, "salvage-1"),
            AlleleAssignment(y.id, 0, "salvage-1"),
        ),
        ledger,
        profile(),
        policy=policy,
    )
    assert report["valid"]
    assert len(measured) == 4


def test_pool_guard_reuses_valid_base_and_matches_full_diagnostics(monkeypatch):
    a = api()
    cat, (x, y), ledger = fixture(((0, 200), (300, 500)))
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=100,
        max_incident_species=100,
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile(), policy)
    measured = []
    real = oracle._score

    def counted(left, right):
        measured.append((left, right))
        return real(left, right)

    monkeypatch.setattr(oracle, "_score", counted)
    assert oracle.pool_valid((x.id,))
    assert len(measured) == 3
    measured.clear()
    assert oracle.pool_valid((x.id, y.id))
    assert len(measured) == 7
    measured.clear()
    assert oracle.pool_valid((y.id, x.id))
    assert not measured
    for ids in ((), (x.id,), (y.id,), (x.id, y.id)):
        assert oracle.pool_valid(ids) == (not oracle.pool_diagnostics(ids)["reasons"])
    assert profile().to_dict()["name"] == "allele-panel-v1"


def test_compact_exposure_counts_never_claim_invalid_pool_is_valid():
    a = api()
    cat, (c,), ledger = fixture(seq=dna(seed=13))
    policy = a.StagePolicy(
        stage_id="salvage-1",
        active_cutoff=-1000,
        max_violating_edges=100,
        max_incident_species=100,
    )
    oracle = a.AlleleCompatibilityOracle(cat, ledger, profile(), policy)
    details = oracle.pool_diagnostics((c.id,))
    assert oracle.pool_exposure_counts((c.id,)) == (
        details["violating_edge_count"],
        details["incident_species_count"],
    )
    assert oracle.pool_exposure_counts(()) == (0, 0)
    strict = a.AlleleCompatibilityOracle(cat, ledger, profile())
    assert not strict.pool_valid((c.id,))
    with pytest.raises(ValueError, match="invalid pool"):
        strict.pool_exposure_counts((c.id,))
