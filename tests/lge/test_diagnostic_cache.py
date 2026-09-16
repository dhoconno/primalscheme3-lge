"""Lossless private diagnostic storage; no scientific or logical cache changes."""

import zlib
from copy import deepcopy
from dataclasses import FrozenInstanceError, dataclass, replace

import pytest
from test_allele_validation import dna, fixture, profile
from test_intended_site_policy import changed

from primalscheme3.panel import allele_validation as validation


@dataclass(frozen=True)
class RawEntry:
    """Reference representation: original full deepcopy retained and returned."""

    reasons: tuple
    result: dict

    @classmethod
    def capture(cls, result):
        return cls(tuple(result["reasons"]), deepcopy(result))

    def unpack(self):
        return deepcopy(self.result)


def test_private_entry_preserves_shapes_and_return_isolation():
    result = {
        "reasons": ["specificity"],
        "nested": [
            {"pair": ("F", "R"), "set": frozenset({"A"}), "number": -26.00000000001}
        ],
    }
    result["alias"] = result["nested"]
    entry = validation._DiagnosticCacheEntry.capture(result)
    assert entry.reasons == ("specificity",) and isinstance(entry.payload, bytes)
    with pytest.raises(FrozenInstanceError):
        entry.payload = b""
    assert entry.unpack() == result
    result["nested"][0]["pair"] = ("changed",)
    one, two = entry.unpack(), entry.unpack()
    assert one["alias"] is one["nested"] and two["alias"] is two["nested"]
    one["nested"][0]["pair"] = ("other",)
    one["reasons"].clear()
    assert two["nested"][0]["pair"] == ("F", "R")
    assert entry.unpack() == two and entry.reasons == ("specificity",)
    with pytest.raises(zlib.error):
        replace(entry, payload=b"not-compressed").unpack()


@pytest.mark.parametrize("intended", ["exact-supported", "concrete-designated-sites"])
@pytest.mark.parametrize(
    "secondary",
    [
        "ordered-disjoint-intended-sites",
        "reject-secondary-products/v1",
        "ordered-disjoint-concrete-designated-sites/v1",
    ],
)
@pytest.mark.parametrize("salvage", [False, True])
def test_native_diagnostics_are_exact_and_isolated(
    monkeypatch, intended, secondary, salvage
):
    seq = dna(seed=0)
    cat, (a, b), ledger = fixture(
        ((0, 200), (300, 500)), seq=seq, rows=(seq, changed(seq, 199))
    )
    p = profile(intended_product_policy=intended, secondary_product_policy=secondary)
    policy = (
        validation.StagePolicy(
            stage_id="salvage-1",
            active_cutoff=-28,
            max_violating_edges=8,
            max_incident_species=4,
        )
        if salvage
        else validation.StagePolicy()
    )
    oracle = validation.AlleleCompatibilityOracle(cat, ledger, p, policy)
    uncached = validation.AlleleCompatibilityOracle(cat, ledger, p, policy, cache=False)
    unary = oracle.candidate_diagnostics(a.id)
    pair = oracle.pair_diagnostics(a.id, b.id)
    assert unary == uncached.candidate_diagnostics(a.id)
    assert pair == uncached.pair_diagnostics(a.id, b.id)
    counters = dict(oracle.specificity.counters)
    saved_unary, saved_pair = deepcopy(unary), deepcopy(pair)
    unary["chemistry"].clear()
    unary["support"].clear()
    unary["reasons"].append("mutated")
    for witness in pair["allowed_secondary_products"]:
        if "certificate" in witness:
            witness["certificate"]["left_forward"]["expected_footprint"][0] = -1
    pair["checks"].clear()
    pair["allowed_secondary_products"].clear()
    pair["rejected_products"].append({"mutated": True})
    assert oracle.candidate_diagnostics(a.id) == saved_unary
    assert oracle.pair_diagnostics(b.id, a.id) == saved_pair
    assert oracle.specificity.counters == counters

    def forbidden(_):
        raise AssertionError("boolean cache hit decoded diagnostics")

    monkeypatch.setattr(validation._DiagnosticCacheEntry, "unpack", forbidden)
    assert oracle.candidate_valid(a.id) == (not saved_unary["reasons"])
    assert oracle.conflict(b.id, a.id) == bool(saved_pair["reasons"])
    assert oracle.specificity.counters == counters
    assert not uncached._candidates and not uncached._pairs


def test_registration_and_cache_false_never_capture(monkeypatch):
    cat, (a, b), ledger = fixture(((0, 200), (300, 500)), seq=dna(seed=0))
    oracle = validation.AlleleCompatibilityOracle(
        cat, replace(ledger, configurations=(a,)), profile()
    )
    assert oracle.candidate_diagnostics(a.id) == oracle.candidate_diagnostics(a.id)
    oracle.register(b)
    assert oracle.pair_diagnostics(a.id, b.id) == oracle.pair_diagnostics(b.id, a.id)
    with pytest.raises(ValueError, match="collision"):
        oracle.register(replace(a, full_interval=(0, 250)))
    unknown = oracle.candidate_diagnostics("unknown")
    assert "unknown-configuration" in unknown["reasons"]
    assert oracle.candidate_diagnostics("unknown") == unknown

    def forbidden(_):
        raise AssertionError("cache=False captured diagnostics")

    monkeypatch.setattr(validation._DiagnosticCacheEntry, "capture", forbidden)
    fresh = validation.AlleleCompatibilityOracle(cat, ledger, profile(), cache=False)
    fresh.candidate_diagnostics(a.id)
    fresh.pair_diagnostics(a.id, b.id)


class Clock:
    def __init__(self):
        self.calls = 0

    def __call__(self):
        self.calls += 1
        return self.calls / 10000


@pytest.mark.parametrize("scheduling", ["serial", "reserved"])
@pytest.mark.parametrize("variant_selection", ["subsets", "full-cloud"])
@pytest.mark.parametrize("salvage", [False, True])
def test_native_fixed_work_search_history_parity(
    monkeypatch, scheduling, variant_selection, salvage
):
    from primalscheme3.core.seq_functions import reverse_complement
    from primalscheme3.panel import allele_search as search
    from primalscheme3.panel.coverage_history import CoverageHistory
    from primalscheme3.panel.coverage_types import ConfigurationLedger, OligoSite
    from primalscheme3.panel.coverage_variants import make_configuration

    seq = dna(seed=0)
    row = changed(seq, 199)
    cat, configs, _ = fixture(
        ((0, 200), (300, 500), (500, 700)), seq=seq, rows=(seq, row)
    )
    family = cat.family_by_id[configs[0].family_id]
    other = OligoSite(
        family.target_id,
        reverse_complement(row[180:200]),
        "-",
        180,
        (180, 200),
        accepting_profile_ids=("normal",),
    )
    family = replace(family, reverse_site_ids=family.reverse_site_ids + (other.id,))
    cat = replace(
        cat,
        sites=cat.sites + (other,),
        families=tuple(family if f.id == family.id else f for f in cat.families),
    )
    configs = tuple(
        make_configuration(cat, f.id, f.forward_site_ids, f.reverse_site_ids)
        for f in cat.families
    )
    ledger = ConfigurationLedger(cat.semantic_digest, configs)
    p = profile(
        intended_product_policy="concrete-designated-sites",
        secondary_product_policy="ordered-disjoint-concrete-designated-sites/v1",
    )
    policy = (
        validation.StagePolicy(
            stage_id="salvage-1",
            active_cutoff=-28,
            max_violating_edges=8,
            max_incident_species=4,
        )
        if salvage
        else validation.StagePolicy()
    )
    options = search.AlleleSearchOptions(
        starts=2,
        repair_rounds=1,
        time_limit=1000,
        phase_scheduling=scheduling,
        variant_selection=variant_selection,
    )
    compressed = validation._DiagnosticCacheEntry
    results = []
    clocks = []
    for entry in (RawEntry, compressed):
        monkeypatch.setattr(validation, "_DiagnosticCacheEntry", entry)
        clock = Clock()
        clocks.append(clock)
        results.append(
            search.search_allele_assignments(
                cat,
                p,
                options,
                ledger,
                policy=policy,
                clock=clock,
                history=CoverageHistory(run_id="diagnostic-parity"),
            )
        )
    a, b = results
    assert a.assignments == b.assignments and a.assignments
    assert (
        a.objective == b.objective
        and a.coverage == b.coverage
        and a.validation == b.validation
    )
    assert a.metadata == b.metadata and a.ledger == b.ledger
    assert clocks[0].calls == clocks[1].calls
    for stream in ("evidence", "assessments", "events", "snapshots"):
        assert getattr(a.history, stream) == getattr(b.history, stream)
    assert b.metadata["work"]["pair_evaluations"] > 0
    if variant_selection == "subsets":
        assert b.metadata["work"]["repair_trials"] > 0
        assert any(
            len(c.reverse_site_ids)
            < len(cat.family_by_id[c.family_id].reverse_site_ids)
            for c in b.ledger.configurations
        )


@pytest.mark.parametrize("cancelled", [False, True])
def test_clock_and_cancellation_are_unchanged(monkeypatch, cancelled):
    from primalscheme3.panel import allele_search as search
    from primalscheme3.panel.coverage_history import CoverageHistory

    cat, _, ledger = fixture(((0, 200), (300, 500)), seq=dna(seed=0))
    compressed = validation._DiagnosticCacheEntry
    results = []
    clocks = []
    for entry in (RawEntry, compressed):
        monkeypatch.setattr(validation, "_DiagnosticCacheEntry", entry)
        clock = Clock()
        clocks.append(clock)
        results.append(
            search.search_allele_assignments(
                cat,
                profile(),
                search.AlleleSearchOptions(time_limit=0.002),
                ledger,
                clock=clock,
                cancelled=(lambda: True) if cancelled else None,
                history=CoverageHistory(run_id="clock-parity"),
            )
        )
    a, b = results
    assert a.assignments == b.assignments and a.metadata == b.metadata
    assert clocks[0].calls == clocks[1].calls
    assert b.metadata["stop_reason"] == ("cancelled" if cancelled else "time-limit")
