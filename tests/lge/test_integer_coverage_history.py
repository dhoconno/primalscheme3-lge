"""Physical history versions preserve public scientific bytes and durability."""

import gc
import sqlite3
import weakref

import pytest

from primalscheme3.panel import coverage_history as history
from tests.lge.test_sqlite_coverage_history import populate


@pytest.mark.parametrize("version", [1, 2])
def test_existing_format_reopens_and_appends_without_migration(tmp_path, version):
    path = tmp_path / str(version)
    with history.SQLiteCoverageHistory(path, run_id="run", format_version=version) as h:
        populate(h)
        expected = h.query(entity_id="site", include_lineage=True)
        metadata = dict(h._db.execute("SELECT key,value FROM metadata"))
        assert metadata["schema"] == f"primalscheme3.sqlite-history/v{version}"
    with history.SQLiteCoverageHistory(path, run_id="run") as h:
        assert h.format_version == version
        assert h.query(entity_id="site", include_lineage=True) == expected
        h.emit(stage_id="later", kind="generated", entity_ids=("new",))
        assert h._db.execute(
            "SELECT type FROM sqlite_master WHERE name='record_links'"
        ).fetchone()[0] == ("table" if version == 1 else "view")
    with history.SQLiteCoverageHistory(path, run_id="run") as h:
        assert h.format_version == version and len(h.events) == 3


def test_default_integer_format_matches_v1_compressed_bytes_and_all_indexes(tmp_path):
    stores = []
    for version in (1, 2):
        h = history.SQLiteCoverageHistory(
            tmp_path / str(version), run_id="run", format_version=version
        )
        populate(h)
        stores.append(h)
    try:
        for sql in (
            "SELECT position,stream,ordinal,id,payload,stage_id,pool,profile_id FROM records ORDER BY position",
            "SELECT * FROM entity_records ORDER BY entity_id,record_id",
            "SELECT * FROM record_links ORDER BY source_id,target_id,kind",
        ):
            assert list(stores[0]._db.execute(sql)) == list(stores[1]._db.execute(sql))
        assert stores[1].format_version == 2
        assert stores[1]._db.execute("PRAGMA foreign_key_check").fetchall() == []
    finally:
        for h in stores:
            h.close()
    with history.SQLiteCoverageHistory(tmp_path / "default", run_id="run") as h:
        assert h.format_version == 2


def test_entity_key_cache_is_instance_bounded_and_invalidated_on_rollback_close(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(history, "SQLITE_ENTITY_CACHE_SIZE", 2)
    h = history.SQLiteCoverageHistory(tmp_path, run_id="run")
    h.emit(stage_id="s", kind="generated", entity_ids=("committed",))
    h.checkpoint()
    for name in ("rolled-back", "two", "three"):
        h.emit(stage_id="s", kind="generated", entity_ids=(name,))
    assert len(h._entity_key_cache) == 2
    h.rollback()
    assert not h._entity_key_cache and len(h.events) == 1
    h.emit(stage_id="s", kind="generated", entity_ids=("reuses-rowid",))
    h.emit(stage_id="s", kind="generated", entity_ids=("rolled-back",))
    assert h.latest_event("reuses-rowid").entity_ids == ("reuses-rowid",)
    assert h.latest_event("rolled-back").entity_ids == ("rolled-back",)
    h.close()
    assert not h._entity_key_cache
    ref = weakref.ref(h)
    del h
    gc.collect()
    assert ref() is None
    with history.SQLiteCoverageHistory(tmp_path, run_id="run") as loaded:
        assert len(loaded.events) == 3
        assert loaded.latest_event("rolled-back").sequence_number == 2


def test_missing_integer_relation_insert_raises_and_rolls_back_pending_batch(
    tmp_path, monkeypatch
):
    h = history.SQLiteCoverageHistory(tmp_path, run_id="run")
    parent = h.emit(stage_id="s", kind="generated", entity_ids=("parent",))
    h.checkpoint()
    # Inject a storage race/fault after ordinary reference validation: no silent
    # INSERT...SELECT omission is allowed even if upstream validation passed.
    original = h._validate_disk_record

    def fault(stream, record, ordinal, **kwargs):
        original(stream, record, ordinal, **kwargs)
        monkeypatch.setattr(h, "_validate_disk_record", original)
        h._db.execute(
            "DELETE FROM entity_records_int WHERE record_key=(SELECT position FROM records WHERE id=?)",
            (parent.id,),
        )
        h._db.execute("DELETE FROM records WHERE id=?", (parent.id,))

    monkeypatch.setattr(h, "_validate_disk_record", fault)
    with pytest.raises(ValueError, match="missing.*relation|relation.*missing"):
        h.emit(
            stage_id="s",
            kind="generated",
            entity_ids=("child",),
            parent_event_ids=(parent.id,),
        )
    assert not h._entity_key_cache
    monkeypatch.setattr(h, "_validate_disk_record", original)
    assert len(h.events) == 1 and h.latest_event("parent") == parent
    h.close()


@pytest.mark.parametrize("version", [1, 2])
def test_explicit_rollback_retains_complete_stage_and_committed_incomplete_tail(
    tmp_path, version
):
    with history.SQLiteCoverageHistory(
        tmp_path, run_id="run", format_version=version
    ) as h:
        populate(h)
        before = tuple(h.events)
        h.emit(stage_id="tail", kind="generated", entity_ids=("pending",))
        h.rollback()
        assert tuple(h.events) == before
        h.emit(stage_id="tail", kind="generated", entity_ids=("durable-tail",))
        h.checkpoint()
    with history.SQLiteCoverageHistory(tmp_path, run_id="run") as h:
        assert len(h.events) == 3 and h.events[-1].entity_ids == ("durable-tail",)
        assert h.last_complete_stage.stage_id == "strict"


def test_failed_explicit_checkpoint_invalidates_cached_rowids_and_pending_counts(
    tmp_path,
):
    h = history.SQLiteCoverageHistory(tmp_path, run_id="run")
    h.emit(stage_id="s", kind="generated", entity_ids=("committed",))
    h.checkpoint()
    h.emit(stage_id="s", kind="generated", entity_ids=("pending",))
    db = h._db

    class FailingCommit:
        def __getattr__(self, name):
            return getattr(db, name)

        def commit(self):
            raise sqlite3.OperationalError("injected commit failure")

    h._db = FailingCommit()
    with pytest.raises(sqlite3.OperationalError, match="commit failure"):
        h.checkpoint()
    assert not h._entity_key_cache and len(h.events) == 1
    h._db = db
    h.emit(stage_id="s", kind="generated", entity_ids=("replacement",))
    h.close()
    with history.SQLiteCoverageHistory(tmp_path, run_id="run") as loaded:
        assert loaded.latest_event("pending") is None
        assert loaded.latest_event("replacement").sequence_number == 1


def test_misdeclared_physical_format_fails_without_schema_mutation(tmp_path):
    import hashlib

    with history.SQLiteCoverageHistory(tmp_path, run_id="run", format_version=1) as h:
        populate(h)
    path = tmp_path / "history.sqlite"
    with sqlite3.connect(path) as db:
        db.execute(
            "UPDATE metadata SET value='primalscheme3.sqlite-history/v2' WHERE key='schema'"
        )
    before = hashlib.sha256(path.read_bytes()).digest()
    with pytest.raises(ValueError, match="schema|format"):
        history.SQLiteCoverageHistory(tmp_path, run_id="run")
    assert hashlib.sha256(path.read_bytes()).digest() == before


def test_missing_entity_relation_source_does_not_commit_partial_record(
    tmp_path, monkeypatch
):
    h = history.SQLiteCoverageHistory(tmp_path, run_id="run")
    original = h._insert_relations

    def fault(stream, ordinal, record):
        h._db.execute("DELETE FROM records WHERE id=?", (record.id,))
        original(stream, ordinal, record)

    monkeypatch.setattr(h, "_insert_relations", fault)
    with pytest.raises(ValueError, match="missing entity relation source"):
        h.emit(stage_id="s", kind="generated", entity_ids=("site",))
    assert len(h.events) == 0 and not h._entity_key_cache
    assert h._db.execute("SELECT COUNT(*) FROM entities").fetchone()[0] == 0
    h.close()


@pytest.mark.parametrize("version", [1, 2])
def test_killed_writer_retains_checkpoint_and_committed_tail_only(tmp_path, version):
    import subprocess
    import sys

    script = """
import sys, time
from primalscheme3.panel.coverage_history import SQLiteCoverageHistory
h = SQLiteCoverageHistory(sys.argv[1], run_id='killed', format_version=int(sys.argv[2]), batch_size=2)
h.emit(stage_id='strict', kind='generated', entity_ids=('base',))
h.complete_stage(stage_id='strict', dispositions={'base':'selected'}, catalog_digest='c', ledger_digest='l')
for name in ('tail-one', 'tail-two', 'pending'):
    h.emit(stage_id='tail', kind='attempted', entity_ids=(name,))
print('ready', flush=True)
time.sleep(30)
"""
    child = subprocess.Popen(
        [sys.executable, "-c", script, str(tmp_path), str(version)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        assert child.stdout.readline().strip() == "ready"
        child.kill()
        child.wait(timeout=10)
        assert child.returncode != 0
    finally:
        if child.poll() is None:
            child.kill()
            child.wait(timeout=10)
    with history.SQLiteCoverageHistory(tmp_path, run_id="killed") as loaded:
        assert loaded.format_version == version
        assert len(loaded.events) == 3 and len(loaded.snapshots) == 1
        assert loaded.last_complete_stage.stage_id == "strict"
        assert loaded.latest_event("pending") is None
        loaded.emit(stage_id="tail", kind="continued", entity_ids=("pending",))
    with history.SQLiteCoverageHistory(tmp_path, run_id="killed") as loaded:
        assert loaded.latest_event("pending").sequence_number == 3


@pytest.mark.parametrize("failure", [KeyboardInterrupt, SystemExit])
@pytest.mark.parametrize("version", [1, 2])
def test_interrupted_relation_write_never_commits_partial_record(
    tmp_path, monkeypatch, failure, version
):
    h = history.SQLiteCoverageHistory(tmp_path, run_id="run", format_version=version)
    h.emit(stage_id="s", kind="generated", entity_ids=("committed",))
    h.checkpoint()

    def interrupt(*args):
        raise failure("interrupted insertion")

    monkeypatch.setattr(h, "_insert_relations", interrupt)
    try:
        with pytest.raises(failure):
            h.emit(stage_id="s", kind="generated", entity_ids=("partial",))
        assert h._db.execute("SELECT COUNT(*) FROM records").fetchone()[0] == 1
        assert not h._entity_key_cache
    finally:
        h.rollback()
        h.close()
    with history.SQLiteCoverageHistory(tmp_path, run_id="run") as loaded:
        assert len(loaded.events) == 1
        assert loaded.latest_event("partial") is None
