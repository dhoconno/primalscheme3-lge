"""Focused compact/full discovery history contracts."""

from primalscheme3.core.config import Config
from primalscheme3.panel.coverage_discovery import build_variant_catalog
from primalscheme3.panel.coverage_history import CoverageHistory
from primalscheme3.panel.coverage_types import Target


def _target(rows):
    width = max(map(len, rows))
    cells = tuple(tuple(row.ljust(width, "-")) for row in rows)
    mapping = tuple(i if cell != "-" else None for i, cell in enumerate(cells[0]))
    reference = "".join(cell for cell in cells[0] if cell != "-")
    reverse = tuple(i for i, value in enumerate(mapping) if value is not None)
    return Target(
        "compact-target",
        0,
        0,
        tuple(str(i) for i in range(len(cells))),
        cells,
        reference,
        len(reference),
        mapping,
        reverse + ((reverse[-1] + 1) if reverse else 0,),
    )


def _projection(catalog):
    return {
        "targets": tuple((t.id, t.rows, t.mapping) for t in catalog.targets),
        "observations": tuple((o.id, o.row_ids, o.cells) for o in catalog.observations),
        "sites": tuple(
            (
                s.id,
                s.target_id,
                s.sequence,
                s.strand,
                s.alignment_anchor,
                s.reference_footprint,
                s.mapping_failure,
                s.accepting_profile_ids,
                s.generated_profile_ids,
            )
            for s in catalog.sites
        ),
        "families": tuple(
            (
                f.id,
                f.target_id,
                f.anchor_pair,
                f.forward_site_ids,
                f.reverse_site_ids,
                f.discovery_profile_combinations,
            )
            for f in catalog.families
        ),
    }


def _fixture():
    shared = "CAACGGCGGACTTTATTGTATCTCC"
    normal = "CATGCTCATTAGGTATATCTTTCAATAAGTTGCA"
    high = "AATCGGCCGACTGTACGCGA"
    rows = [normal + "ATATATATAT" + high, shared + "ATATATATAT" + high]
    target = _target(rows)
    config = Config(amplicon_size=55, amplicon_size_min=50, amplicon_size_max=60)
    end = len(normal)
    return target, config, ([end], [end + 10])


def test_compact_and_full_have_identical_scientific_projection_and_bounded_history():
    target, config, indexes = _fixture()
    compact_history = CoverageHistory(None, run_id="compact")
    full_history = CoverageHistory(None, run_id="full")
    compact = build_variant_catalog(
        (target,), config, indexes=indexes, history=compact_history, history_detail="compact"
    )
    full = build_variant_catalog(
        (target,), config, indexes=indexes, history=full_history, history_detail="full"
    )

    assert compact.sites and compact.families
    assert _projection(compact) == _projection(full)
    assert compact.resolved_config_json != full.resolved_config_json
    assert not compact_history.evidence
    assert not compact_history.assessments
    assert all(e.kind not in {"generated", "assessed", "generated-family", "family-geometry"}
               for e in compact_history.events)
    assert any(e.kind == "compact-discovery-summary" for e in compact_history.events)
    assert compact_history.snapshots[-1].completeness == "complete"
    assert set(compact_history.snapshots[-1].dispositions) == {target.id}


def test_compact_retains_rejected_mapped_sites_with_empty_evidence_ids():
    shared = "CAACGGCGGACTTTATTGTATCTCC"
    target = _target([shared, "A" * len(shared)])
    history = CoverageHistory(None, run_id="compact")
    catalog = build_variant_catalog(
        (target,), Config(min_base_freq=0.5), indexes=([len(shared)], []),
        length_mode="all", history=history, history_detail="compact"
    )
    assert all(not site.intrinsic_evidence_ids for site in catalog.sites)
    assert any(site.mapping_failure is None and not site.accepting_profile_ids for site in catalog.sites)


def test_compact_preserves_chemistry_membership_for_reference_rejected_sites():
    shared = "CAACGGCGGACTTTATTGTATCTCC"
    target = _target([shared + "-", shared + "-"])
    config = Config()
    compact = build_variant_catalog(
        (target,), config, indexes=([len(target.rows[0])], []),
        length_mode="all", history_detail="compact"
    )
    full = build_variant_catalog(
        (target,), config, indexes=([len(target.rows[0])], []),
        length_mode="all", history_detail="full"
    )
    compact_rejected = [s for s in compact.sites if s.mapping_failure]
    full_rejected = [s for s in full.sites if s.mapping_failure]
    assert compact_rejected and full_rejected
    assert any(s.accepting_profile_ids for s in compact_rejected)
    assert {
        (s.id, s.accepting_profile_ids, s.generated_profile_ids)
        for s in compact_rejected
    } == {
        (s.id, s.accepting_profile_ids, s.generated_profile_ids)
        for s in full_rejected
    }


def test_compact_summary_counts_raw_attempts_and_failed_reasons():
    target, config, indexes = _fixture()
    history = CoverageHistory(None, run_id="compact")
    build_variant_catalog((target,), config, indexes=indexes, history=history, history_detail="compact")
    summary = next(e for e in history.events if e.kind == "compact-discovery-summary")
    assert summary.changes["raw_attempt_count"] >= summary.changes["concrete_attempt_count"]
    assert summary.changes["profile_id"] in {"normal", "high-gc"}
    assert "failed_attempts_by_reason" in summary.changes


def test_compact_family_members_are_stable_when_raw_records_arrive_reversed(monkeypatch):
    from primalscheme3.panel import coverage_discovery

    target, config, indexes = _fixture()
    expected = build_variant_catalog(
        (target,), config, indexes=indexes, history_detail="compact"
    )
    original = coverage_discovery.discover

    def reversed_discover(*args, **kwargs):
        records, workers = original(*args, **kwargs)
        return list(reversed(records)), workers

    monkeypatch.setattr(coverage_discovery, "discover", reversed_discover)
    actual = build_variant_catalog(
        (target,), config, indexes=indexes, history_detail="compact"
    )
    assert _projection(actual) == _projection(expected)


def test_compact_skips_diagnostic_record_construction(monkeypatch):
    from primalscheme3.panel import coverage_discovery

    target, config, indexes = _fixture()
    history = CoverageHistory(None, run_id="compact")

    def forbidden(*args, **kwargs):
        raise AssertionError("compact discovery constructed diagnostic history")

    monkeypatch.setattr(history, "record_evidence", forbidden)
    monkeypatch.setattr(history, "assess", forbidden)
    monkeypatch.setattr(coverage_discovery, "_group_variant_records", forbidden)
    catalog = build_variant_catalog(
        (target,), config, indexes=indexes, history=history, history_detail="compact"
    )
    assert catalog.sites and catalog.families
