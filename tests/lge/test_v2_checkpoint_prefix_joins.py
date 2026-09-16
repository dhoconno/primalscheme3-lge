"""Checkpoint closure has the same verdicts with physical v2 integer links."""

import sqlite3

import pytest

from primalscheme3.panel.coverage_history import SQLiteCoverageHistory, StageSnapshot

KINDS = (
    ("evidence", "assessments", "assessments_count", "evidence_count"),
    ("assessment", "events", "events_count", "assessments_count"),
    ("parent", "events", "events_count", "events_count"),
)


def _history_with_all_link_kinds(path, version):
    store = SQLiteCoverageHistory(path, run_id="prefix", format_version=version)
    evidence = [
        store.record_evidence(
            entity_ids=(f"site-{i}",),
            measurement="tm",
            dependency_key={"row": i},
            values={"tm": 60 + i},
            status="measured",
        )
        for i in range(2)
    ]
    assessments = [
        store.assess(
            stage_id="strict",
            entity_ids=(f"site-{i}",),
            pool=0,
            context_digest="context",
            profile_id="normal",
            kernel_versions={},
            thresholds={},
            check_name="tm",
            outcome="pass",
            reason="ok",
            evidence_ids=(evidence[i].id,),
        )
        for i in range(2)
    ]
    first = store.emit(
        stage_id="strict",
        kind="generated",
        entity_ids=("site-0",),
        assessment_ids=(assessments[0].id,),
    )
    store.emit(
        stage_id="strict",
        kind="generated",
        entity_ids=("site-1",),
        assessment_ids=(assessments[1].id,),
        parent_event_ids=(first.id,),
    )
    return store


def _snapshot(store, evidence_count, assessments_count, events_count):
    return StageSnapshot(
        stage_id="strict",
        completeness="incomplete",
        prior_snapshot_id=None,
        dispositions={},
        catalog_digest="catalog",
        ledger_digest="ledger",
        evidence_digest=store._prefix_digest("evidence", evidence_count),
        assessments_digest=store._prefix_digest("assessments", assessments_count),
        events_digest=store._prefix_digest("events", events_count),
        evidence_count=evidence_count,
        assessments_count=assessments_count,
        events_count=events_count,
    )


def _original_closure_hits(store, snapshot):
    """Original ID/view SQL, independently applied to each of the three kinds."""
    return {
        kind: store._db.execute(
            """SELECT 1 FROM record_links l JOIN records s ON s.id=l.source_id
                JOIN records t ON t.id=l.target_id WHERE l.kind=? AND s.stream=?
                AND s.ordinal<? AND t.ordinal>=? LIMIT 1""",
            (
                kind,
                source_stream,
                getattr(snapshot, source_count),
                getattr(snapshot, target_count),
            ),
        ).fetchone()
        is not None
        for kind, source_stream, source_count, target_count in KINDS
    }


@pytest.mark.parametrize("version", [1, 2])
def test_checkpoint_closure_matches_original_queries_for_full_and_truncated_prefixes(
    tmp_path, version
):
    with _history_with_all_link_kinds(tmp_path / str(version), version) as store:
        cases = (
            ((2, 2, 2), {"evidence": False, "assessment": False, "parent": False}),
            ((2, 1, 1), {"evidence": False, "assessment": False, "parent": False}),
            ((1, 2, 2), {"evidence": True, "assessment": False, "parent": False}),
            ((2, 1, 2), {"evidence": False, "assessment": True, "parent": False}),
        )
        for counts, expected_hits in cases:
            snapshot = _snapshot(store, *counts)
            assert _original_closure_hits(store, snapshot) == expected_hits
            statements = []
            store._db.set_trace_callback(statements.append)
            try:
                if any(expected_hits.values()):
                    with pytest.raises(ValueError, match="checkpoint references uncommitted prefix"):
                        store._validate_checkpoint(snapshot)
                else:
                    store._validate_checkpoint(snapshot)
            finally:
                store._db.set_trace_callback(None)
            closure_sql = [sql for sql in statements if "record_links" in sql]
            assert closure_sql
            if version == 2:
                assert all("record_links_int" in sql for sql in closure_sql)
                assert all("s.position=l.source_key" in sql for sql in closure_sql)
                assert all("t.position=l.target_key" in sql for sql in closure_sql)
            else:
                assert all("record_links l" in sql for sql in closure_sql)
                assert all("s.id=l.source_id" in sql for sql in closure_sql)


@pytest.mark.parametrize("version", [1, 2])
def test_parent_ordinal_corruption_still_rejected_at_full_equal_count_prefix(
    tmp_path, version
):
    with _history_with_all_link_kinds(tmp_path / str(version), version) as store:
        snapshot = _snapshot(store, 2, 2, 2)
        assert not any(_original_closure_hits(store, snapshot).values())
        store.checkpoint()
        store._db.execute(
            "UPDATE records SET ordinal=3 WHERE stream='events' AND ordinal=0"
        )
        assert _original_closure_hits(store, snapshot) == {
            "evidence": False,
            "assessment": False,
            "parent": True,
        }
        with pytest.raises(ValueError, match="checkpoint references uncommitted prefix"):
            store._validate_checkpoint(snapshot)
        store.rollback()
        assert len(store.events) == 2


def test_v2_full_stage_export_and_reopen_match_v1(tmp_path):
    exports = []
    snapshots = []
    for version in (1, 2):
        path = tmp_path / str(version)
        with _history_with_all_link_kinds(path, version) as store:
            snapshot = store.complete_stage(
                stage_id="strict",
                dispositions={"site-0": "selected", "site-1": "selected"},
                catalog_digest="catalog",
                ledger_digest="ledger",
            )
            snapshots.append(snapshot)
            exports.append(store.export(tmp_path / f"export-{version}"))
        with SQLiteCoverageHistory(path, run_id="prefix") as loaded:
            assert loaded.last_complete_stage == snapshot
            assert loaded.query(entity_id="site-1", include_lineage=True)["events"]
    assert snapshots[0] == snapshots[1]
    assert exports[0] == exports[1]


@pytest.mark.parametrize("version", [1, 2])
def test_valid_truncated_checkpoint_survives_cold_reload(tmp_path, version):
    path = tmp_path / str(version)
    with _history_with_all_link_kinds(path, version) as store:
        valid = _snapshot(store, 2, 1, 1)
        store._append("snapshots", valid)
    with SQLiteCoverageHistory(path, run_id="prefix") as loaded:
        assert tuple(loaded.snapshots) == (valid,)
        assert len(loaded.evidence) == 2
        assert len(loaded.assessments) == 2
        assert len(loaded.events) == 2


@pytest.mark.parametrize("version", [1, 2])
@pytest.mark.parametrize("corruption", ["unique-id", "dangling-link"])
def test_cold_reload_still_rejects_identity_and_foreign_key_corruption(
    tmp_path, version, corruption
):
    path = tmp_path / f"{version}-{corruption}"
    with _history_with_all_link_kinds(path, version) as store:
        store.complete_stage(
            stage_id="strict",
            dispositions={"site-0": "selected", "site-1": "selected"},
            catalog_digest="catalog",
            ledger_digest="ledger",
        )
    with sqlite3.connect(path / "history.sqlite") as db:
        if corruption == "unique-id":
            db.execute(
                "UPDATE records SET id='forged-record-id' "
                "WHERE stream='evidence' AND ordinal=0"
            )
        elif version == 1:
            db.execute(
                "UPDATE record_links SET target_id='missing-record-id' "
                "WHERE kind='evidence'"
            )
        else:
            db.execute(
                "UPDATE record_links_int SET target_key=999999 "
                "WHERE kind='evidence'"
            )
    with pytest.raises(ValueError, match="corrupt history|record index"):
        SQLiteCoverageHistory(path, run_id="prefix")
