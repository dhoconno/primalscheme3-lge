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


def test_sqlite_page_cache_default_preserves_durability_settings(tmp_path):
    with history.SQLiteCoverageHistory(tmp_path, run_id='run') as h:
        assert history.SQLITE_PAGE_CACHE_KIB == 65536
        assert h._db.execute('PRAGMA cache_size').fetchone()[0] == -65536
        assert h._db.execute('PRAGMA synchronous').fetchone()[0] == 2
        assert h._db.execute('PRAGMA journal_mode').fetchone()[0] == 'delete'
        assert h._db.execute('PRAGMA foreign_keys').fetchone()[0] == 1
        assert h.batch_size == 1000
