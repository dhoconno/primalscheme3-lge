"""Literal synthetic specificity geometries; chemistry isolated where indicated."""

import random
from dataclasses import replace

import pytest

from primalscheme3.core.config import Config
from primalscheme3.core.seq_functions import reverse_complement as rc
from primalscheme3.core.thermo import THERMO_RESULT
from primalscheme3.panel.coverage_specificity import SpecificityChecker
from primalscheme3.panel.coverage_types import Assignment, Candidate, Catalog, Target
from primalscheme3.panel.coverage_validation import (
    CompatibilityOracle,
    ConstraintProfile,
    validate_assignments,
)


def dna(length=1000, seed=91):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def target(sequence, ident="t", rows=None):
    rows = rows or [sequence]
    mapping = []
    coordinates = []
    position = 0
    for i, base in enumerate(rows[0]):
        mapping.append(None if base in {"-", ""} else position)
        if base not in {"-", ""}:
            coordinates.append(i)
            position += 1
    coordinates.append(coordinates[-1] + 1)
    return Target(
        ident,
        0,
        0,
        tuple(f"{ident}-r{i}" for i in range(len(rows))),
        tuple(tuple(r) for r in rows),
        sequence,
        len(sequence),
        tuple(mapping),
        tuple(coordinates),
    )


def candidate(t, start, end, ident="a", size=20):
    # Deliberately stale cached support: validator must reconstruct it.
    return Candidate(
        ident,
        t.id,
        (start, end),
        (start + size, end - size),
        (t.reference_sequence[start : start + size],),
        (rc(t.reference_sequence[end - size : end]),),
        (),
        (),
        (),
    )


def catalog(targets, candidates):
    return Catalog(tuple(targets), tuple(candidates), "synthetic", "{}", (), ())


def profile(bound=2000, **kwargs):
    config = Config(
        amplicon_size=200,
        amplicon_size_min=40,
        amplicon_size_max=1000,
        mismatch_product_size=bound,
        mismatch_fuzzy=False,
    )
    config.mismatch_kmersize = 20
    return ConstraintProfile.from_config(config, **kwargs)


@pytest.fixture
def geometry_only(monkeypatch):
    # Geometry fixtures use arbitrary DNA, not thermochemically selected primers.
    monkeypatch.setattr(
        "primalscheme3.panel.coverage_validation.thermo_check",
        lambda *_: THERMO_RESULT.PASS,
    )
    monkeypatch.setattr(
        "primalscheme3.panel.coverage_validation.do_pool_interact", lambda *_: False
    )


def test_declared_and_disjoint_products_and_reverse_orientation():
    t = target(dna())
    a, b = candidate(t, 0, 200), candidate(t, 300, 500, "b")
    checker = SpecificityChecker(catalog([t], [a, b]), profile())
    assert not checker.intrinsic(a).rejected_products
    assert not checker.intrinsic(b).rejected_products
    result = checker.pair(a, b)
    assert not result.rejected_products
    assert result == checker.pair(b, a)
    assert [
        (p["start"], p["end"], p["length"]) for p in result.allowed_secondary_products
    ] == [(0, 500, 500)]
    hits = checker.hits(a.reverse_oligos[0])
    assert any(h.orientation == "-" and h.start == 180 and h.end == 200 for h in hits)


def test_inward_overlap_rejected_in_both_orders():
    t = target(dna())
    a, b = candidate(t, 0, 200), candidate(t, 150, 350, "b")
    checker = SpecificityChecker(catalog([t], [a, b]), profile())
    result = checker.pair(a, b)
    assert result == checker.pair(b, a)
    assert any(p["start"] == 150 and p["end"] == 200 for p in result.rejected_products)


def test_extra_within_candidate_off_anchor_product_is_rejected():
    seq = dna()
    seq = seq[:60] + seq[:20] + seq[80:]
    t = target(seq)
    a = candidate(t, 0, 200)
    result = SpecificityChecker(catalog([t], [a]), profile()).intrinsic(a)
    assert any(p["start"] == 60 and p["end"] == 200 for p in result.rejected_products)


def test_cross_msa_and_duplicate_occurrence_hits_are_not_discarded():
    t, other = target(dna()), target(dna(), "other")
    a = candidate(t, 0, 200)
    result = SpecificityChecker(catalog([t, other], [a]), profile()).intrinsic(a)
    assert any(p["target_id"] == "other" for p in result.rejected_products)


def test_different_rows_never_join_hits_and_unsupported_candidate_fails(geometry_only):
    seq = dna()
    f_only = seq[:180] + "N" * 20 + seq[200:]
    r_only = "N" * 20 + seq[20:]
    t = target(seq, rows=[f_only, r_only])
    a = candidate(t, 0, 200)
    cat = catalog([t], [a])
    assert not SpecificityChecker(cat, profile()).intrinsic(a).rejected_products
    oracle = CompatibilityOracle(cat, profile())
    assert not oracle.candidate_valid(a.id)
    assert "no_joint_support" in oracle.candidate_reasons(a.id)


@pytest.mark.parametrize("bound", [0, -1])
def test_nonpositive_bound_rejected(bound):
    with pytest.raises(ValueError, match="positive"):
        profile(bound)


def test_product_bound_inclusive_and_crossed_three_prime_ends_excluded():
    seq = dna()
    seq = seq[:300] + seq[:20] + seq[320:480] + seq[180:200] + seq[500:]
    t = target(seq)
    a = candidate(t, 0, 200)
    cat = catalog([t], [a])
    assert SpecificityChecker(cat, profile(200)).intrinsic(a).rejected_products
    assert not SpecificityChecker(cat, profile(199)).intrinsic(a).rejected_products
    adjacent = candidate(t, 0, 40, "adjacent")
    crossed = replace(
        adjacent,
        full_interval=(0, 30),
        interior_interval=(20, 10),
        reverse_oligos=(rc(seq[10:30]),),
    )
    assert (
        not SpecificityChecker(catalog([t], [crossed]), profile())
        .intrinsic(crossed)
        .rejected_products
    )
    assert (
        not SpecificityChecker(catalog([t], [adjacent]), profile())
        .intrinsic(adjacent)
        .rejected_products
    )


def test_shared_oligo_declared_signature_union_and_no_third_candidate_rescue():
    seq = dna()
    seq = seq[:300] + seq[:20] + seq[320:]
    t = target(seq)
    a, b = candidate(t, 0, 200), candidate(t, 300, 500, "b")
    checker = SpecificityChecker(catalog([t], [a, b]), profile(200))
    assert not checker.intrinsic(a).rejected_products
    assert not checker.intrinsic(b).rejected_products
    assert not checker.pair(a, b).rejected_products
    assert checker.pair(a, b) == checker.pair(b, a)
    # C declares the inward off-anchor cross product but cannot rescue P(A,B).
    t2 = target(dna())
    a2, b2, c2 = (
        candidate(t2, 0, 200),
        candidate(t2, 150, 350, "b"),
        candidate(t2, 150, 200, "c"),
    )
    small = SpecificityChecker(catalog([t2], [a2, b2]), profile()).pair(a2, b2)
    large = SpecificityChecker(catalog([t2], [a2, b2, c2]), profile()).pair(a2, b2)
    assert small.rejected_products and small == large


def test_fuzzy_and_ambiguous_hits_cannot_grant_intended_support():
    seq = dna()
    mutated = ("A" if seq[0] != "A" else "C") + seq[1:20]
    row = seq[:60] + mutated + seq[80:]
    t = target(seq, rows=[row])
    a = candidate(t, 0, 200)
    cat = catalog([t], [a])
    assert not SpecificityChecker(cat, profile()).intrinsic(a).rejected_products
    fuzzy = replace(profile(), mismatch_fuzzy=True)
    assert SpecificityChecker(cat, fuzzy).intrinsic(a).rejected_products
    ambiguous = seq[:5] + ("R" if seq[5] in "AG" else "Y") + seq[6:]
    ta = target(seq, rows=[ambiguous])
    result = SpecificityChecker(catalog([ta], [a]), profile()).intrinsic(a)
    assert result.rejected_products
    assert result.support["joint_rows"] == []
    assert result.support["unknown_rows"] == ["t-r0"]


def test_gapped_row_products_use_ungapped_coordinates():
    seq = dna()
    t = target(seq, rows=[seq[:10] + "---" + seq[10:]])
    a, b = candidate(t, 0, 200), candidate(t, 300, 500, "b")
    result = SpecificityChecker(catalog([t], [a, b]), profile()).pair(a, b)
    assert not result.rejected_products
    assert [(p["start"], p["end"]) for p in result.allowed_secondary_products] == [
        (0, 500)
    ]


def test_validator_rebuilds_coverage_after_removal_and_checks_pools(geometry_only):
    t = target(dna(100))
    a, b = candidate(t, 0, 80), candidate(t, 10, 90, "b")
    cat = catalog([t], [a, b])
    assignments = [Assignment("a", 0), Assignment("b", 1)]
    report = validate_assignments(cat, assignments, profile(), "primer-trimmed", 0.9)
    assert report["valid"]
    remaining = validate_assignments(
        cat, assignments[1:], profile(), "primer-trimmed", 0.9
    )
    assert remaining["per_target"]["t"]["interior"]["covered_bases"] == 40
    assert remaining["per_target"]["t"]["interior"]["intervals"] == [[30, 70]]
    assert remaining["per_target"]["t"]["full"]["covered_bases"] == 80
    assert not validate_assignments(
        cat, [Assignment("b", 2)], profile(), "full-span", 0.9
    )["valid"]
    assert not validate_assignments(
        cat, [Assignment("a", 0), Assignment("b", 0)], profile(), "full-span", 0.9
    )["valid"]
    assert remaining["support_diagnostics"]["b"]["joint_rows"] == ["t-r0"]


def test_caps_duplicates_bounds_and_fresh_validation(geometry_only):
    t = target(dna())
    a, b = candidate(t, 0, 200), candidate(t, 300, 500, "b")
    cat = catalog([t], [a, b])
    p = profile(max_amplicons=1, max_amplicons_msa=1)
    oracle = CompatibilityOracle(cat, p)
    assert oracle.candidate_valid("a")
    assert not oracle.conflict("a", "b")
    assert not oracle.conflict("b", "a")
    assert oracle.counters["pair_evaluations"] == 1
    assert not validate_assignments(
        cat, [Assignment("a", 0), Assignment("b", 1)], p, "full-span", 0.9
    )["valid"]
    assert not validate_assignments(
        cat, [Assignment("a", 0), Assignment("a", 1)], p, "full-span", 0.9
    )["valid"]
    assert not validate_assignments(
        cat, [Assignment("missing", 0)], p, "full-span", 0.9
    )["valid"]
    bad = replace(a, full_interval=(-1, 200))
    assert not CompatibilityOracle(catalog([t], [bad]), profile()).candidate_valid("a")
    short = replace(
        a, forward_oligos=(a.forward_oligos[0][-18:],), full_interval=(2, 200)
    )
    assert "oligo_too_short" in CompatibilityOracle(
        catalog([t], [short]), profile()
    ).candidate_reasons("a")


def test_profile_serializes_resolved_thresholds_and_freezes_source_config():
    config = Config(mismatch_product_size=2000)
    p = ConstraintProfile.from_config(config, max_amplicons=7, max_amplicons_msa=2)
    serialized = p.to_dict()
    assert serialized["name"] == "panel-v1"
    assert serialized["specificity_revision"] == "intended-sites-v1"
    assert (
        serialized["allowed_secondary_product_policy"]
        == "ordered-disjoint-intended-sites"
    )
    assert serialized["thermochemistry"]["primer_tm_max"] == 62.5
    assert serialized["thermochemistry"]["effective_tm_upper_offset"] == 2
    config.primer_tm_max = 100
    assert p.to_dict() == serialized
    with pytest.raises(ValueError):
        ConstraintProfile.from_config(
            Config(use_matchdb=False, mismatch_product_size=2000)
        )


def test_real_chemistry_rejects_unselected_arbitrary_oligos():
    t = target("A" * 100)
    a = candidate(t, 0, 80)
    oracle = CompatibilityOracle(catalog([t], [a]), profile())
    assert not oracle.candidate_valid("a")
    assert any(reason.startswith("thermo:") for reason in oracle.candidate_reasons("a"))
    assert oracle.counters["thermo_calls"] > 0
    assert oracle.counters["interaction_calls"] > 0


def test_pair_off_anchor_shared_sequence_is_not_intended_and_is_symmetric():
    seq = dna()
    seq = seq[:350] + seq[:20] + seq[370:]
    t = target(seq)
    a, b = candidate(t, 0, 200), candidate(t, 300, 500, "b")
    checker = SpecificityChecker(catalog([t], [a, b]), profile(200))
    assert not checker.intrinsic(a).rejected_products
    assert not checker.intrinsic(b).rejected_products
    result = checker.pair(a, b)
    assert result == checker.pair(b, a)
    assert [(p["start"], p["end"]) for p in result.rejected_products] == [(350, 500)]


def test_terminal_only_anchor_match_and_unavailable_projection_not_exempt():
    seq = dna()
    changed = ("A" if seq[0] != "A" else "C") + seq[1:]
    t = target(seq, rows=[changed])
    a = candidate(t, 0, 200)
    p = replace(profile(), mismatch_kmersize=19)
    result = SpecificityChecker(catalog([t], [a]), p).intrinsic(a)
    assert result.rejected_products and not result.support["joint_rows"]
    shifted = target(seq, rows=[seq[1:]])
    hits = SpecificityChecker(catalog([shifted], [a]), p).hits(a.forward_oligos[0])
    assert any(
        h.start == -1 and h.classification == "uncertain-footprint" for h in hits
    )


def test_same_oligo_opposite_orientation_extra_product_is_checked():
    seq = dna()
    seq = seq[:80] + rc(seq[:20]) + seq[100:]
    t = target(seq)
    a = candidate(t, 0, 200)
    result = SpecificityChecker(catalog([t], [a]), profile()).intrinsic(a)
    assert any(
        p["start"] == 0
        and p["end"] == 100
        and p["plus"]["oligo"] == p["minus"]["oligo"]
        for p in result.rejected_products
    )


def test_final_validator_never_trusts_cached_oracle_verdicts(geometry_only):
    t = target(dna())
    a, b = candidate(t, 0, 200), candidate(t, 150, 350, "b")
    cat = catalog([t], [a, b])
    oracle = CompatibilityOracle(cat, profile())
    assert oracle.conflict("a", "b")
    oracle._pair_results[("a", "b")] = {"reasons": []}
    oracle._evaluator.specificity._hits.clear()
    assert not oracle.conflict("a", "b")
    report = validate_assignments(
        cat, [Assignment("a", 0), Assignment("b", 0)], profile(), "full-span", 0.9
    )
    assert not report["valid"]
    assert any(v["reason"] == "specificity" for v in report["violations"])
    assert report["counters"]["pair_evaluations"] == 1


def test_actual_tm_upper_offset_and_native_interaction_are_preserved():
    # This concrete 23-mer has Tm 63.9175 under the existing default chemistry;
    # native thermo_check accepts up to primer_tm_max + 2, here 64.5.
    oligo = "CAGTTGTCTTGTTCGCGAGTCGT"
    seq = dna()
    seq = oligo + seq[len(oligo) : 177] + rc(oligo) + seq[200:]
    t = target(seq)
    a = candidate(t, 0, 200, size=len(oligo))
    oracle = CompatibilityOracle(catalog([t], [a]), profile())
    assert "thermo:HIGH_TM" not in oracle.candidate_reasons("a")
    assert oracle.counters["thermo_calls"] == 1
    assert oracle.counters["interaction_calls"] == 1
    assert oracle.candidate_valid("a")


def test_lazy_hit_index_does_not_repeat_row_window_scans():
    t = target(dna())
    a = candidate(t, 0, 200)
    checker = SpecificityChecker(catalog([t], [a]), profile())
    checker.hits(a.forward_oligos[0])
    scanned = checker.counters["row_windows_indexed"]
    assert scanned == 981
    checker.hits(a.reverse_oligos[0])
    assert checker.counters["row_windows_indexed"] == scanned
