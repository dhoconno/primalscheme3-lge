"""Disk history preserves scientific stream identity without retaining records."""
import json
import sqlite3
import tracemalloc
import zlib
import pytest
from primalscheme3.panel import coverage_history as history


def populate(h):
    e = h.record_evidence(entity_ids=('site',), measurement='tm', dependency_key={'v': 1},
                          values={'tm': 61}, status='measured')
    a = h.assess(stage_id='strict', entity_ids=('site',), pool=0, context_digest='context',
                 profile_id='normal', kernel_versions={}, thresholds={'tm': 60}, check_name='tm',
                 outcome='pass', reason='ok', evidence_ids=(e.id,))
    event = h.emit(stage_id='strict', kind='generated', entity_ids=('site',), assessment_ids=(a.id,))
    complete = h.complete_stage(stage_id='strict', dispositions={'site': 'selected'},
                                catalog_digest='catalog', ledger_digest='ledger')
    tail = h.emit(stage_id='salvage', kind='reconsidered', entity_ids=('child',),
                  parent_event_ids=(event.id,))
    failed = h.complete_stage(stage_id='salvage', dispositions={}, catalog_digest='catalog',
                              ledger_digest='ledger', prior_snapshot_id=complete.id, completeness='failed')
    return e, a, event, complete, tail, failed


def test_sqlite_matches_memory_streams_hashes_queries_and_reload(tmp_path):
    memory = history.CoverageHistory(None, run_id='run')
    disk = history.SQLiteCoverageHistory(tmp_path, run_id='run', batch_size=3)
    assert populate(memory) == populate(disk)
    for stream in ('evidence', 'assessments', 'events', 'snapshots'):
        assert tuple(getattr(memory, stream)) == tuple(getattr(disk, stream))
    assert disk.query(entity_id='child', include_lineage=True) == memory.query(entity_id='child', include_lineage=True)
    assert disk.query(entity_id='site', pool=0) == memory.query(entity_id='site', pool=0)
    expected = disk.last_complete_stage
    assert disk.latest_event('child').kind == 'reconsidered'
    disk.close()
    loaded = history.SQLiteCoverageHistory(tmp_path, run_id='run')
    assert loaded.last_complete_stage == expected
    assert len(loaded.events) == 2
    assert loaded.latest_event('child').sequence_number == 1
    assert loaded.export(tmp_path/'export1') == loaded.export(tmp_path/'export2')
    loaded.close()


def test_sqlite_reference_and_fresh_disposition_validation(tmp_path):
    h = history.SQLiteCoverageHistory(tmp_path, run_id='run')
    with pytest.raises(ValueError, match='reference'):
        h.emit(stage_id='strict', kind='bad', entity_ids=('site',), parent_event_ids=('missing',))
    h.emit(stage_id='strict', kind='generated', entity_ids=('site',))
    first = h.complete_stage(stage_id='strict', dispositions={'site': 'selected'}, catalog_digest='c', ledger_digest='l')
    h.emit(stage_id='salvage', kind='moved', entity_ids=('site',))
    with pytest.raises(ValueError, match='disposition'):
        h.complete_stage(stage_id='salvage', dispositions={}, catalog_digest='c', ledger_digest='l', prior_snapshot_id=first.id)
    h.close()


@pytest.mark.parametrize('corruption', ['payload', 'deleted_prefix', 'entity_index', 'reference_index'])
def test_sqlite_corruption_is_rejected_on_reload(tmp_path, corruption):
    h = history.SQLiteCoverageHistory(tmp_path, run_id='run')
    populate(h)
    h.close()
    db = sqlite3.connect(tmp_path/'history.sqlite')
    if corruption == 'payload':
        db.execute("UPDATE records SET payload=? WHERE stream='events' AND ordinal=0", (zlib.compress(b'{}'),))
    elif corruption == 'deleted_prefix':
        db.execute("DELETE FROM records WHERE stream='evidence'")
    elif corruption == 'entity_index':
        db.execute("DELETE FROM entity_records WHERE entity_id='site'")
    else:
        db.execute("DELETE FROM record_links WHERE kind='evidence'")
    db.commit(); db.close()
    with pytest.raises(ValueError, match='corrupt|checkpoint|reference|index|prefix'):
        history.SQLiteCoverageHistory(tmp_path, run_id='run')


def test_sqlite_checkpoint_tail_and_portable_backup(tmp_path):
    h = history.SQLiteCoverageHistory(tmp_path/'source', run_id='run', batch_size=100)
    strict = h.emit(stage_id='strict', kind='generated', entity_ids=('site',))
    snapshot = h.complete_stage(stage_id='strict', dispositions={'site': 'selected'}, catalog_digest='c', ledger_digest='l')
    h.emit(stage_id='salvage', kind='attempted', entity_ids=('child',), parent_event_ids=(strict.id,))
    h.checkpoint()
    h.export_database(tmp_path/'portable'/'history.sqlite')
    moved = history.SQLiteCoverageHistory(tmp_path/'portable', run_id='run')
    assert moved.last_complete_stage == snapshot
    assert len(moved.events) == 2 and len(moved.snapshots) == 1
    moved.close(); h.close()


def test_sqlite_record_storage_memory_does_not_grow_with_payload_volume(tmp_path):
    h = history.SQLiteCoverageHistory(tmp_path, run_id='run', batch_size=100)
    tracemalloc.start()
    for i in range(1200):
        h.record_evidence(entity_ids=('site',), measurement='large', dependency_key={'i': i},
                          values={'payload': 'ACGT' * 4000 + str(i)}, status='measured')
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert len(h.evidence) == 1200
    assert current < 3_000_000 and peak < 5_000_000
    assert not isinstance(h._records['evidence'], list)
    h.close()


def test_crash_retains_committed_batch_and_discards_only_pending_tail(tmp_path):
    import subprocess
    import sys
    script = '''
import os, sys
from primalscheme3.panel.coverage_history import SQLiteCoverageHistory
h=SQLiteCoverageHistory(sys.argv[1], run_id='crash', batch_size=2)
for i in range(3):
    h.emit(stage_id='discovery',kind='generated',entity_ids=(str(i),))
os._exit(0)
'''
    subprocess.run([sys.executable, '-c', script, str(tmp_path)], check=True)
    with history.SQLiteCoverageHistory(tmp_path, run_id='crash') as h:
        assert len(h.events) == 2
        assert h.events[-1].entity_ids == ('1',)
        assert h.last_complete_stage is None
        h.emit(stage_id='discovery', kind='continued', entity_ids=('2',))
    with history.SQLiteCoverageHistory(tmp_path, run_id='crash') as h:
        assert len(h.events) == 3
        assert h.events[-1].sequence_number == 2


def test_valid_snapshot_identity_with_wrong_prefix_hash_is_rejected(tmp_path):
    from dataclasses import replace
    from primalscheme3.panel.coverage_types import canonical_json
    with history.SQLiteCoverageHistory(tmp_path, run_id='run') as h:
        h.emit(stage_id='strict', kind='generated', entity_ids=('site',))
        original = h.complete_stage(stage_id='strict', dispositions={'site': 'selected'}, catalog_digest='c', ledger_digest='l')
    corrupt = replace(original, events_digest='history-events-wrong')
    with sqlite3.connect(tmp_path/'history.sqlite') as db:
        db.execute("UPDATE records SET id=?,payload=? WHERE stream='snapshots'",
                   (corrupt.id, zlib.compress(canonical_json(corrupt).encode())))
    with pytest.raises(ValueError, match='prefix hash mismatch'):
        history.SQLiteCoverageHistory(tmp_path, run_id='run')


def test_sqlite_discovery_catalog_and_exported_records_match_memory(tmp_path):
    import gzip
    from primalscheme3.core.config import Config
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.coverage_types import Target
    sequence = 'CAACGGCGGACTTTATTGTATCTCC'
    target = Target('target', 0, 0, ('row',), (tuple(sequence),), sequence, len(sequence),
                    tuple(range(len(sequence))), tuple(range(len(sequence)+1)))
    memory = history.CoverageHistory(None, run_id='discovery')
    with history.SQLiteCoverageHistory(tmp_path/'db', run_id='discovery', batch_size=7) as disk:
        first = build_variant_catalog((target,), Config(), indexes=([len(sequence)], []), history=memory)
        second = build_variant_catalog((target,), Config(), indexes=([len(sequence)], []), history=disk)
        assert first.semantic_digest == second.semantic_digest
        assert first == second
        disk.export(tmp_path/'disk-export')
        memory.export(tmp_path/'memory-export')
        for name in ('evidence', 'assessments', 'events', 'snapshots'):
            assert gzip.decompress((tmp_path/'disk-export'/f'{name}.jsonl.gz').read_bytes()) == gzip.decompress(
                (tmp_path/'memory-export'/f'{name}.jsonl.gz').read_bytes())


def test_event_range_iterator_preserves_chronological_offsets(tmp_path):
    for h in (history.CoverageHistory(None, run_id='range'), history.SQLiteCoverageHistory(tmp_path, run_id='range')):
        events = [h.emit(stage_id='strict', kind='generated', entity_ids=(str(i),)) for i in range(7)]
        assert list(h.iter_events(start=3)) == events[3:]
        assert list(h.iter_events(start=2, stop=4)) == events[2:4]
        assert list(h.iter_events(start=7)) == []
        if hasattr(h, 'close'):
            h.close()


def test_hot_caches_match_disabled_and_skip_duplicate_decode(tmp_path, monkeypatch):
    with history.SQLiteCoverageHistory(tmp_path/'hot', run_id='run') as hot, history.SQLiteCoverageHistory(tmp_path/'cold', run_id='run', hot_cache_bytes=0) as cold:
        assert populate(hot) == populate(cold)
        evidence = hot.evidence[0]
        def unexpected(*args):
            raise AssertionError('hot duplicate must not decode the stored payload')
        monkeypatch.setattr(hot, '_decode', unexpected)
        assert hot._append('evidence', evidence) == evidence
        monkeypatch.undo()
        for stream in hot._streams:
            assert tuple(getattr(hot, stream)) == tuple(getattr(cold, stream))
            assert hot._prefix_digest(stream, hot._counts[stream]) == cold._prefix_digest(stream, cold._counts[stream])


def test_hot_cache_budget_eviction_oversize_and_close(tmp_path):
    h = history.SQLiteCoverageHistory(tmp_path, run_id='run', hot_cache_bytes=4096)
    records = [h.record_evidence(entity_ids=('site',), measurement='x', dependency_key={'i': i}, values={'x': 'a'*100}, status='measured') for i in range(20)]
    assert h._payload_cache.bytes_used <= 3072
    assert h._reference_cache.bytes_used <= 1024
    assert h._payload_cache.get(records[0].id) is None
    large = h.record_evidence(entity_ids=('site',), measurement='large', dependency_key={}, values={'x': 'a'*5000}, status='measured')
    assert h._payload_cache.get(large.id) is None
    assert h._append('evidence', records[0]) == records[0]
    h.close()
    assert h._payload_cache.bytes_used == h._reference_cache.bytes_used == 0


def test_hot_cache_collision_wrong_stream_and_external_corruption(tmp_path):
    from dataclasses import replace
    h = history.SQLiteCoverageHistory(tmp_path, run_id='run')
    e = h.record_evidence(entity_ids=('site',), measurement='x', dependency_key={}, values={'x': 1}, status='measured')
    forged = replace(e, values={'x': 2})
    object.__setattr__(forged, 'id', e.id)
    with pytest.raises(ValueError, match='collision'):
        h._append('evidence', forged)
    with pytest.raises(ValueError, match='reference'):
        h.emit(stage_id='strict', kind='bad', entity_ids=('site',), parent_event_ids=(e.id,))
    h.checkpoint()
    with sqlite3.connect(h.database_path) as db:
        db.execute('UPDATE records SET payload=? WHERE id=?', (zlib.compress(b'{}'), e.id))
    with pytest.raises(ValueError, match='corrupt'):
        h._append('evidence', e)
    h.close()
    with pytest.raises(ValueError, match='corrupt'):
        history.SQLiteCoverageHistory(tmp_path, run_id='run')


def test_hot_cache_rollback_and_reload_are_cold(tmp_path):
    h = history.SQLiteCoverageHistory(tmp_path, run_id='run')
    e = h.record_evidence(entity_ids=('site',), measurement='x', dependency_key={}, values={'x': 1}, status='measured')
    h.rollback()
    assert h._payload_cache.bytes_used == h._reference_cache.bytes_used == 0
    assert len(h.evidence) == 0
    assert h._append('evidence', e) == e
    h.checkpoint()
    h._reload()
    assert h._payload_cache.bytes_used == h._reference_cache.bytes_used == 0
    assert len(h.evidence) == 1
    h.close()


def test_hot_references_skip_sql_but_cold_reload_revalidates(tmp_path):
    with history.SQLiteCoverageHistory(tmp_path, run_id='run') as h:
        parent = h.emit(stage_id='strict', kind='generated', entity_ids=('site',))
        statements = []
        h._db.set_trace_callback(statements.append)
        h.emit(stage_id='strict', kind='child', entity_ids=('site',), parent_event_ids=(parent.id,))
        assert not any('SELECT stream,position' in sql for sql in statements)
        h.checkpoint()
        statements.clear()
        h._reload()
        assert any('SELECT stream,position' in sql for sql in statements)
        assert h._reference_cache.bytes_used == 0
        assert h._db.execute('PRAGMA synchronous').fetchone()[0] == 2
        assert h._db.execute('PRAGMA cache_size').fetchone()[0] == -8192
        assert h.batch_size == 1000


def test_byte_lru_refreshes_recency_and_bypasses_large_entries():
    cache = history._ByteLRU(1100)
    cache.put('a', ('v', 1))
    cache.put('b', ('v', 2))
    assert cache.get('a') == ('v', 1)
    cache.put('c', ('v', 3))
    assert cache.get('b') is None
    assert cache.get('a') == ('v', 1)
    cache.put('large', ('a'*2000, 1))
    assert cache.get('large') is None
    assert cache.bytes_used <= cache.budget


def test_direct_sql_and_rollback_invalidate_hot_payloads(tmp_path):
    with history.SQLiteCoverageHistory(tmp_path, run_id='run') as h:
        e = h.record_evidence(entity_ids=('site',), measurement='x', dependency_key={}, values={}, status='measured')
        h._db.rollback()
        assert h._append('evidence', e) == e
        assert len(h.evidence) == 1
        h._db.execute('UPDATE records SET payload=? WHERE id=?', (zlib.compress(b'{}'), e.id))
        with pytest.raises(ValueError, match='corrupt'):
            h._append('evidence', e)
        assert h._payload_cache.bytes_used == h._reference_cache.bytes_used == 0
        h.rollback()


@pytest.mark.parametrize('budget', [-1, True, 1.5])
def test_hot_cache_budget_requires_nonnegative_integer(tmp_path, budget):
    with pytest.raises(ValueError, match='hot_cache_bytes'):
        history.SQLiteCoverageHistory(tmp_path, run_id='run', hot_cache_bytes=budget)


def test_checkpoint_does_not_mask_direct_sql_changes(tmp_path):
    with history.SQLiteCoverageHistory(tmp_path, run_id='run') as h:
        e = h.record_evidence(entity_ids=('site',), measurement='x', dependency_key={}, values={}, status='measured')
        h._db.execute('UPDATE records SET payload=? WHERE id=?', (zlib.compress(b'{}'), e.id))
        h.checkpoint()
        with pytest.raises(ValueError, match='corrupt'):
            h._append('evidence', e)
