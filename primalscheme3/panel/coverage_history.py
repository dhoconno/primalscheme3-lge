"""Append-only scientific evidence, contextual verdicts and decision lineage.

A single coordinator owns a history writer. Workers return evidence batches;
callers merge them deterministically before assigning event sequence numbers.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import sqlite3
import zlib
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path
from typing import Any

from primalscheme3.panel.coverage_types import (
    ScientificRecord,
    canonical_json,
    semantic_id,
)

OUTCOMES = frozenset(('pass', 'fail', 'unknown', 'not-evaluated'))


@dataclass(frozen=True)
class IntrinsicEvidence(ScientificRecord):
    entity_ids: tuple[str, ...]
    measurement: str
    dependency_key: Mapping[str, Any]
    values: Mapping[str, Any]
    status: str
    id: str = field(init=False)
    _identity_fields = ('entity_ids', 'measurement', 'dependency_key', 'values', 'status')


@dataclass(frozen=True)
class AssessmentRecord(ScientificRecord):
    run_id: str
    stage_id: str
    entity_ids: tuple[str, ...]
    pool: int | None
    context_digest: str
    profile_id: str
    kernel_versions: Mapping[str, Any]
    thresholds: Mapping[str, Any]
    check_name: str
    outcome: str
    reason: str
    witness_ids: tuple[str, ...] = ()
    evidence_ids: tuple[str, ...] = ()
    id: str = field(init=False)
    _identity_fields = ('run_id', 'stage_id', 'entity_ids', 'pool', 'context_digest', 'profile_id',
                        'kernel_versions', 'thresholds', 'check_name', 'outcome', 'reason',
                        'witness_ids', 'evidence_ids')

    def __post_init__(self):
        if self.outcome not in OUTCOMES:
            raise ValueError('assessment requires an explicit pass/fail/unknown/not-evaluated outcome')
        if self.pool is not None and self.pool < 0:
            raise ValueError('negative pool')
        super().__post_init__()


@dataclass(frozen=True)
class DecisionEvent(ScientificRecord):
    run_id: str
    sequence_number: int
    stage_id: str
    kind: str
    entity_ids: tuple[str, ...]
    parent_event_ids: tuple[str, ...] = ()
    assessment_ids: tuple[str, ...] = ()
    changes: Mapping[str, Any] = field(default_factory=dict)
    id: str = field(init=False)
    _identity_fields = ('run_id', 'sequence_number', 'stage_id', 'kind', 'entity_ids',
                        'parent_event_ids', 'assessment_ids', 'changes')


@dataclass(frozen=True)
class StageSnapshot(ScientificRecord):
    stage_id: str
    completeness: str
    prior_snapshot_id: str | None
    dispositions: Mapping[str, str]
    catalog_digest: str
    ledger_digest: str
    evidence_digest: str
    assessments_digest: str
    events_digest: str
    evidence_count: int = 0
    assessments_count: int = 0
    events_count: int = 0
    id: str = field(init=False)
    _identity_fields = ('stage_id', 'completeness', 'prior_snapshot_id', 'dispositions', 'catalog_digest',
                        'ledger_digest', 'evidence_digest', 'assessments_digest', 'events_digest',
                        'evidence_count', 'assessments_count', 'events_count')

    def __post_init__(self):
        if self.completeness not in ('complete', 'incomplete', 'failed'):
            raise ValueError('invalid snapshot completeness')
        super().__post_init__()


class CoverageHistory:
    """Persist each record before exposing it. None directory gives an in-memory writer.

    Loading fails closed on corruption/truncation; source files are never repaired
    or overwritten. Complete preceding snapshots remain available in their stream.
    """
    _streams = {'evidence': IntrinsicEvidence, 'assessments': AssessmentRecord,
                'events': DecisionEvent, 'snapshots': StageSnapshot}

    def __init__(self, directory: Path | None = None, *, run_id: str):
        self.directory = Path(directory) if directory is not None else None
        self.run_id = run_id
        self._records = {name: [] for name in self._streams}
        self._indexes = {name: {} for name in self._streams}
        self._latest_event_by_entity: dict[str, DecisionEvent] = {}
        if self.directory is not None:
            self.directory.mkdir(parents=True, exist_ok=True)
            for name, cls in self._streams.items():
                path = self.directory / (name + '.jsonl')
                if not path.exists():
                    continue
                for number, line in enumerate(path.read_text().splitlines(), 1):
                    try:
                        record = cls.from_dict(json.loads(line))
                        self._validate(name, record)
                        self._remember(name, record)
                    except (ValueError, TypeError, KeyError) as exc:
                        raise ValueError(f'incomplete or corrupt history: {path.name}:{number}: {exc}') from exc

    def _validate(self, name, record):
        if record.id in self._indexes[name]:
            raise ValueError('duplicate persisted record')
        if hasattr(record, 'run_id') and record.run_id != self.run_id:
            raise ValueError('history belongs to another run')
        if name == 'assessments' and not set(record.evidence_ids) <= self._indexes['evidence'].keys():
            raise ValueError('unknown evidence reference')
        if name == 'events':
            if record.sequence_number != len(self._records['events']):
                raise ValueError('non-contiguous decision sequence')
            if not set(record.parent_event_ids) <= self._indexes['events'].keys():
                raise ValueError('unknown parent event')
            if not set(record.assessment_ids) <= self._indexes['assessments'].keys():
                raise ValueError('unknown assessment')
        if name == 'snapshots' and record.prior_snapshot_id is not None:
            prior = self._indexes['snapshots'].get(record.prior_snapshot_id)
            if prior is None or prior.completeness != 'complete':
                raise ValueError('snapshot inheritance requires a complete prior snapshot')
        if name == 'snapshots':
            self._validate_checkpoint(record)

    def _validate_checkpoint(self, snapshot):
        prefixes = {}
        for name in ('evidence', 'assessments', 'events'):
            count = getattr(snapshot, name + '_count')
            if type(count) is not int or not 0 <= count <= len(self._records[name]):
                raise ValueError('checkpoint committed prefix is missing: ' + name)
            prefix = tuple(self._records[name][:count])
            if semantic_id('history-' + name, prefix) != getattr(snapshot, name + '_digest'):
                raise ValueError('checkpoint committed prefix hash mismatch: ' + name)
            prefixes[name] = prefix
        evidence_ids = {r.id for r in prefixes['evidence']}
        assessment_ids = {r.id for r in prefixes['assessments']}
        if any(not set(a.evidence_ids) <= evidence_ids for a in prefixes['assessments']):
            raise ValueError('checkpoint assessment references uncommitted evidence')
        if any(not set(e.assessment_ids) <= assessment_ids for e in prefixes['events']):
            raise ValueError('checkpoint event references uncommitted assessment')
        prior = self._indexes['snapshots'].get(snapshot.prior_snapshot_id)
        inherited = self._resolved_dispositions(prior) if prior is not None else {}
        if prior is not None and any(getattr(snapshot, name + '_count') < getattr(prior, name + '_count')
                                     for name in ('evidence', 'assessments', 'events')):
            raise ValueError('checkpoint prefixes precede inherited snapshot')
        if snapshot.completeness == 'complete':
            seen = {entity for event in prefixes['events'] for entity in event.entity_ids}
            touched = {entity for event in prefixes['events'][prior.events_count if prior else 0:]
                       for entity in event.entity_ids}
            if not seen <= (set(snapshot.dispositions) | set(inherited)) or not touched <= set(snapshot.dispositions):
                raise ValueError('complete snapshot requires fresh dispositions for touched entities')

    def _remember(self, name, record):
        self._records[name].append(record)
        self._indexes[name][record.id] = record
        if name == 'events':
            for entity_id in record.entity_ids:
                self._latest_event_by_entity[entity_id] = record

    def _append(self, name, record):
        existing = self._indexes[name].get(record.id)
        if existing is not None:
            if existing != record:
                raise ValueError('record identity collision')
            return existing
        self._validate(name, record)
        if self.directory is not None:
            with (self.directory / (name + '.jsonl')).open('a') as handle:
                handle.write(canonical_json(record) + '\n')
                handle.flush()
                os.fsync(handle.fileno())
        self._remember(name, record)
        return record

    @property
    def evidence(self): return tuple(self._records['evidence'])

    @property
    def assessments(self): return tuple(self._records['assessments'])

    @property
    def events(self): return tuple(self._records['events'])

    @property
    def snapshots(self): return tuple(self._records['snapshots'])

    @property
    def last_complete_stage(self) -> StageSnapshot | None:
        return next((s for s in reversed(self.snapshots) if s.completeness == 'complete'), None)

    def latest_event(self, entity_id: str) -> DecisionEvent | None:
        """O(1) latest causal event lookup, maintained on emit and reload."""
        return self._latest_event_by_entity.get(entity_id)

    def iter_events(self, start=0, stop=None):
        """Iterate a chronological event range without copying a stream slice."""
        return islice(self._records['events'], start, stop)

    def record_evidence(self, *, entity_ids, measurement, dependency_key, values, status):
        return self._append('evidence', IntrinsicEvidence(entity_ids, measurement, dependency_key, values, status))

    def assess(self, *, stage_id, entity_ids, pool, context_digest, profile_id, kernel_versions,
               thresholds, check_name, outcome, reason, witness_ids=(), evidence_ids=()):
        return self._append('assessments', AssessmentRecord(self.run_id, stage_id, entity_ids, pool,
            context_digest, profile_id, kernel_versions, thresholds, check_name, outcome, reason,
            witness_ids, evidence_ids))

    def emit(self, *, stage_id, kind, entity_ids, parent_event_ids=(), assessment_ids=(), changes=None):
        return self._append('events', DecisionEvent(self.run_id, len(self._records['events']), stage_id, kind,
            entity_ids, parent_event_ids, assessment_ids, {} if changes is None else changes))

    def _resolved_dispositions(self, snapshot):
        result = {}
        if snapshot.prior_snapshot_id is not None:
            result.update(self._resolved_dispositions(self._indexes['snapshots'][snapshot.prior_snapshot_id]))
        result.update(snapshot.dispositions)
        return result

    def complete_stage(self, *, stage_id, dispositions, catalog_digest, ledger_digest,
                       prior_snapshot_id=None, completeness='complete'):
        inherited = {}
        if prior_snapshot_id is not None:
            prior = self._indexes['snapshots'].get(prior_snapshot_id)
            if prior is None or prior.completeness != 'complete':
                raise ValueError('snapshot inheritance requires complete prior snapshot')
            inherited = self._resolved_dispositions(prior)
        seen = {entity for event in self.events for entity in event.entity_ids}
        if completeness == 'complete' and not seen <= (set(dispositions) | set(inherited)):
            raise ValueError('complete snapshot requires disposition for every generated/explored entity')
        return self._append('snapshots', StageSnapshot(stage_id, completeness, prior_snapshot_id,
            dispositions, catalog_digest, ledger_digest,
            semantic_id('history-evidence', self.evidence), semantic_id('history-assessments', self.assessments),
            semantic_id('history-events', self.events), len(self._records['evidence']),
            len(self._records['assessments']), len(self._records['events'])))

    def query(self, *, entity_id=None, pool=None, stage_id=None, profile_id=None, include_lineage=False):
        assessments = tuple(a for a in self.assessments
            if (entity_id is None or entity_id in a.entity_ids)
            and (pool is None or pool == a.pool) and (stage_id is None or stage_id == a.stage_id)
            and (profile_id is None or profile_id == a.profile_id))
        events = [e for e in self.events if (entity_id is None or entity_id in e.entity_ids)
                  and (stage_id is None or stage_id == e.stage_id)]
        if include_lineage:
            ids = {e.id for e in events}
            while True:
                expanded = ids | {p for e in self.events if e.id in ids for p in e.parent_event_ids}
                expanded |= {e.id for e in self.events if set(e.parent_event_ids) & ids}
                if expanded == ids:
                    break
                ids = expanded
            events = [e for e in self.events if e.id in ids]
        evidence_ids = {i for a in assessments for i in a.evidence_ids}
        evidence = tuple(e for e in self.evidence if e.id in evidence_ids or entity_id in e.entity_ids)
        return {'evidence': evidence, 'assessments': assessments, 'events': tuple(events),
                'snapshots': tuple(s for s in self.snapshots if stage_id is None or s.stage_id == stage_id)}

    def export(self, directory: Path) -> dict[str, str]:
        """Deterministically compressed streams; returned hashes cover final bytes.

        This low-level serializer does not replace invocation provenance or atomic
        scientific publication in the pipeline. Snapshot hashes exclude snapshots.
        """
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        hashes = {}
        for name, records in self._records.items():
            data = ''.join(canonical_json(r)+'\n' for r in records).encode()
            compressed = gzip.compress(data, mtime=0)
            filename = name + '.jsonl.gz'
            (directory / filename).write_bytes(compressed)
            hashes[filename] = hashlib.sha256(compressed).hexdigest()
        index = {s.id: {'stage_id': s.stage_id, 'completeness': s.completeness,
                        'dispositions': self._resolved_dispositions(s)} for s in self.snapshots}
        data = canonical_json(index).encode()
        (directory/'snapshot-index.json').write_bytes(data)
        hashes['snapshot-index.json'] = hashlib.sha256(data).hexdigest()
        return hashes


# The JSONL implementation above remains the compatibility default. Production
# runs may explicitly select this disk-backed implementation of the same API.


# Bounded 100/500-anchor probes favored 64 MiB over 8 MiB; 256 MiB did
# not improve the larger probe. Keep FULL/DELETE and the commit cadence intact.
SQLITE_PAGE_CACHE_KIB = 64 * 1024
SQLITE_ENTITY_CACHE_SIZE = 32768
SQLITE_HISTORY_SCHEMAS = frozenset(('primalscheme3.sqlite-history/v1', 'primalscheme3.sqlite-history/v2'))


def sqlite_history_format(connection, metadata=None):
    """Validate the explicit physical format without writing or migrating it."""
    if metadata is None:
        metadata = dict(connection.execute('SELECT key,value FROM metadata'))
    schema = metadata.get('schema')
    if schema not in SQLITE_HISTORY_SCHEMAS:
        raise ValueError('unsupported history schema')
    version = int(schema.rsplit('v', 1)[1])
    objects = dict(connection.execute("SELECT name,type FROM sqlite_master WHERE type IN ('table','view')"))
    required = {'metadata': 'table', 'records': 'table',
                'entity_records': 'table' if version == 1 else 'view',
                'record_links': 'table' if version == 1 else 'view'}
    if version == 2:
        required.update(entities='table', entity_records_int='table', record_links_int='table')
    if any(objects.get(name) != kind for name, kind in required.items()):
        raise ValueError('corrupt history physical schema')
    return version


class _DiskSequence(Sequence):
    def __init__(self, owner, stream):
        self.owner, self.stream = owner, stream

    def __len__(self):
        return self.owner._counts[self.stream]

    def __iter__(self):
        cursor = self.owner._db.execute('SELECT payload FROM records WHERE stream=? ORDER BY ordinal', (self.stream,))
        for (payload,) in cursor:
            yield self.owner._decode(self.stream, payload)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(self[i] for i in range(*index.indices(len(self))))
        if index < 0:
            index += len(self)
        row = self.owner._db.execute('SELECT payload FROM records WHERE stream=? AND ordinal=?', (self.stream, index)).fetchone()
        if row is None:
            raise IndexError(index)
        return self.owner._decode(self.stream, row[0])


class _DiskIndex:
    def __init__(self, owner, stream):
        self.owner, self.stream = owner, stream

    def get(self, identity, default=None):
        row = self.owner._db.execute('SELECT payload FROM records WHERE stream=? AND id=?', (self.stream, identity)).fetchone()
        return self.owner._decode(self.stream, row[0]) if row else default

    def __getitem__(self, identity):
        result = self.get(identity)
        if result is None:
            raise KeyError(identity)
        return result


def _history_hash_start(stream):
    # Exactly semantic_id('history-'+stream, records), with its array streamed.
    prefix = '{"canonicalVersion":1,"kind":' + json.dumps('history-' + stream) + ',"semanticFields":['
    return hashlib.sha256(prefix.encode())


class SQLiteCoverageHistory(CoverageHistory):
    """Compressed, indexed full history with bounded record-storage memory.

    Stream properties are lazy read-only sequences. Every batch, checkpoint(),
    complete_stage(), and close() commits with SQLite synchronous=FULL. A crash
    can lose only the current uncommitted batch; completed stages are durable.
    New histories use integer-relation v2; format_version=1 creates legacy v1.
    Existing histories always retain their detected format, without migration.
    rollback() discards the entire pending batch and restores the durable prefix.
    Reload verifies payload identities, references, indexes and snapshot hashes
    without materializing full record streams. One writer owns the connection.
    """
    def __init__(self, directory: Path, *, run_id: str, batch_size: int = 1000, format_version: int = 2):
        if type(batch_size) is not int or batch_size < 1:
            raise ValueError('batch_size must be a positive integer')
        if type(format_version) is not int or format_version not in (1, 2):
            raise ValueError('format_version must be 1 or 2 for a new history')
        self._entity_key_cache = OrderedDict()
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.run_id, self.batch_size = run_id, batch_size
        self.database_path = self.directory / 'history.sqlite'
        self._db = sqlite3.connect(self.database_path)
        self._pending = 0
        self._closed = False
        self._counts = {name: 0 for name in self._streams}
        self._hashes = {name: _history_hash_start(name) for name in self._streams}
        self._records = {name: _DiskSequence(self, name) for name in self._streams}
        self._indexes = {name: _DiskIndex(self, name) for name in self._streams}
        try:
            self._db.execute('PRAGMA foreign_keys=ON')
            self._db.execute('PRAGMA synchronous=FULL')
            self._db.execute('PRAGMA journal_mode=DELETE')
            self._db.execute(f'PRAGMA cache_size=-{SQLITE_PAGE_CACHE_KIB}')
            self._db.execute('PRAGMA temp_store=FILE')
            # Detect an existing format before DDL. In particular v2 exposes
            # read views where v1 has tables; opening either never migrates it.
            tables = {row[0] for row in self._db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if not tables:
                self.format_version = format_version
                self._create_schema()
                self._db.executemany('INSERT INTO metadata VALUES (?,?)',
                    [('schema', f'primalscheme3.sqlite-history/v{format_version}'), ('run_id', run_id),
                     ('committed_counts', canonical_json(self._counts))])
                self._db.commit()
            else:
                if 'metadata' not in tables:
                    raise ValueError('corrupt history: missing metadata')
                metadata = dict(self._db.execute('SELECT key,value FROM metadata'))
                if metadata.get('schema') not in SQLITE_HISTORY_SCHEMAS or metadata.get('run_id') != run_id:
                    raise ValueError('history belongs to another run or unsupported schema')
                self.format_version = sqlite_history_format(self._db, metadata)
            self._reload()
        except Exception:
            self._db.close()
            self._closed = True
            raise

    def _create_schema(self):
        self._db.executescript("""
            CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE records(
                position INTEGER PRIMARY KEY, stream TEXT NOT NULL, ordinal INTEGER NOT NULL,
                id TEXT NOT NULL UNIQUE, payload BLOB NOT NULL, stage_id TEXT, pool INTEGER, profile_id TEXT,
                UNIQUE(stream, ordinal));
            CREATE INDEX records_filter ON records(stream, stage_id, pool, profile_id, ordinal);
        """)
        if self.format_version == 1:
            self._db.executescript("""
                CREATE TABLE entity_records(
                    entity_id TEXT NOT NULL, stream TEXT NOT NULL, ordinal INTEGER NOT NULL,
                    record_id TEXT NOT NULL REFERENCES records(id), PRIMARY KEY(entity_id, record_id));
                CREATE INDEX entity_records_lookup ON entity_records(entity_id, stream, ordinal);
                CREATE INDEX entity_records_record ON entity_records(record_id);
                CREATE TABLE record_links(
                    source_id TEXT NOT NULL REFERENCES records(id), target_id TEXT NOT NULL REFERENCES records(id),
                    kind TEXT NOT NULL, PRIMARY KEY(source_id, target_id, kind));
                CREATE INDEX record_links_target ON record_links(target_id, kind, source_id);
            """)
        else:
            self._db.executescript("""
                CREATE TABLE entities(entity_key INTEGER PRIMARY KEY, canonical_id TEXT NOT NULL UNIQUE);
                CREATE TABLE entity_records_int(
                    entity_key INTEGER NOT NULL REFERENCES entities(entity_key), stream TEXT NOT NULL,
                    ordinal INTEGER NOT NULL, record_key INTEGER NOT NULL REFERENCES records(position),
                    PRIMARY KEY(entity_key, record_key));
                CREATE INDEX entity_int_lookup ON entity_records_int(entity_key, stream, ordinal);
                CREATE INDEX entity_int_record ON entity_records_int(record_key);
                CREATE TABLE record_links_int(
                    source_key INTEGER NOT NULL REFERENCES records(position),
                    target_key INTEGER NOT NULL REFERENCES records(position), kind TEXT NOT NULL,
                    PRIMARY KEY(source_key, target_key, kind));
                CREATE INDEX links_int_target ON record_links_int(target_key, kind, source_key);
                CREATE VIEW entity_records AS
                    SELECT e.canonical_id entity_id,x.stream,x.ordinal,r.id record_id
                    FROM entity_records_int x JOIN entities e ON e.entity_key=x.entity_key
                    JOIN records r ON r.position=x.record_key;
                CREATE VIEW record_links AS
                    SELECT s.id source_id,t.id target_id,x.kind
                    FROM record_links_int x JOIN records s ON s.position=x.source_key
                    JOIN records t ON t.position=x.target_key;
            """)

    def _entity_key(self, identity):
        cached = self._entity_key_cache.get(identity)
        if cached is not None:
            self._entity_key_cache.move_to_end(identity)
            return cached
        self._db.execute('INSERT OR IGNORE INTO entities(canonical_id) VALUES (?)', (identity,))
        key = self._db.execute('SELECT entity_key FROM entities WHERE canonical_id=?', (identity,)).fetchone()[0]
        self._entity_key_cache[identity] = key
        if len(self._entity_key_cache) > SQLITE_ENTITY_CACHE_SIZE:
            self._entity_key_cache.popitem(last=False)
        return key

    def _insert_relations(self, stream, ordinal, record):
        entities = set(getattr(record, 'entity_ids', ()))
        links = self._links(record)
        if self.format_version == 1:
            self._db.executemany('INSERT INTO entity_records VALUES (?,?,?,?)',
                                 ((e, stream, ordinal, record.id) for e in entities))
            self._db.executemany('INSERT INTO record_links VALUES (?,?,?)',
                                 ((record.id, target, kind) for target, kind, _ in links))
            return
        for entity in entities:
            written = self._db.execute(
                'INSERT INTO entity_records_int SELECT ?,?,?,position FROM records WHERE id=?',
                (self._entity_key(entity), stream, ordinal, record.id))
            if written.rowcount != 1:
                raise ValueError('missing entity relation source record')
        for target, kind, _ in links:
            written = self._db.execute("""INSERT INTO record_links_int
                SELECT s.position,t.position,? FROM records s,records t WHERE s.id=? AND t.id=?""",
                (kind, record.id, target))
            if written.rowcount != 1:
                raise ValueError('missing history relation source or target record')

    def _decode(self, stream, payload):
        try:
            raw = zlib.decompress(payload)
            record = self._streams[stream].from_dict(json.loads(raw))
            if canonical_json(record).encode() != raw:
                raise ValueError('noncanonical payload')
            return record
        except (ValueError, TypeError, KeyError, zlib.error) as exc:
            raise ValueError('corrupt history payload: ' + stream) from exc

    @staticmethod
    def _links(record):
        links = []
        for attribute, kind, stream in (('evidence_ids', 'evidence', 'evidence'),
                                       ('assessment_ids', 'assessment', 'assessments'),
                                       ('parent_event_ids', 'parent', 'events')):
            links.extend((identity, kind, stream) for identity in getattr(record, attribute, ()))
        if getattr(record, 'prior_snapshot_id', None) is not None:
            links.append((record.prior_snapshot_id, 'prior', 'snapshots'))
        return set(links)

    def _validate_disk_record(self, stream, record, ordinal, *, position=None):
        if hasattr(record, 'run_id') and record.run_id != self.run_id:
            raise ValueError('history belongs to another run')
        if stream == 'events' and record.sequence_number != ordinal:
            raise ValueError('non-contiguous decision sequence')
        for identity, kind, expected in self._links(record):
            row = self._db.execute('SELECT stream,position FROM records WHERE id=?', (identity,)).fetchone()
            if row is None or row[0] != expected or (position is not None and row[1] >= position):
                raise ValueError('unknown or forward ' + kind + ' reference')
        if stream == 'snapshots':
            self._validate_checkpoint(record)

    def _reload(self):
        if self._db.execute('PRAGMA integrity_check').fetchone()[0] != 'ok' or self._db.execute('PRAGMA foreign_key_check').fetchone():
            raise ValueError('corrupt history database/reference')
        if self._db.execute('SELECT 1 FROM records WHERE stream NOT IN (?,?,?,?) LIMIT 1', tuple(self._streams)).fetchone():
            raise ValueError('corrupt history stream')
        for stream in self._streams:
            for position, ordinal, identity, payload, stage, pool, profile in self._db.execute(
                    'SELECT position,ordinal,id,payload,stage_id,pool,profile_id FROM records WHERE stream=? ORDER BY ordinal', (stream,)):
                if ordinal != self._counts[stream]:
                    raise ValueError('corrupt history ordinal prefix')
                record = self._decode(stream, payload)
                if (identity, stage, pool, profile) != (record.id, getattr(record, 'stage_id', None),
                                                      getattr(record, 'pool', None), getattr(record, 'profile_id', None)):
                    raise ValueError('corrupt history record index')
                actual = set(self._db.execute('SELECT entity_id,stream,ordinal FROM entity_records WHERE record_id=?', (identity,)))
                expected = {(entity, stream, ordinal) for entity in getattr(record, 'entity_ids', ())}
                if actual != expected:
                    raise ValueError('corrupt history entity index')
                links = set(self._db.execute('SELECT target_id,kind FROM record_links WHERE source_id=?', (identity,)))
                if links != {(target, kind) for target, kind, _ in self._links(record)}:
                    raise ValueError('corrupt history reference index')
                self._validate_disk_record(stream, record, ordinal, position=position)
                self._update_hash(stream, zlib.decompress(payload))
                self._counts[stream] += 1
        committed = self._db.execute("SELECT value FROM metadata WHERE key='committed_counts'").fetchone()
        if committed is None or json.loads(committed[0]) != self._counts:
            raise ValueError('corrupt history committed prefix counts')

    def _update_hash(self, stream, raw):
        if self._counts[stream]:
            self._hashes[stream].update(b',')
        self._hashes[stream].update(raw)

    def _prefix_digest(self, stream, count):
        if type(count) is not int or not 0 <= count <= self._counts[stream]:
            raise ValueError('checkpoint committed prefix is missing: ' + stream)
        if count == self._counts[stream]:
            digest = self._hashes[stream].copy()
        else:
            digest = _history_hash_start(stream)
            for ordinal, payload in self._db.execute('SELECT ordinal,payload FROM records WHERE stream=? AND ordinal<? ORDER BY ordinal', (stream, count)):
                if ordinal:
                    digest.update(b',')
                digest.update(zlib.decompress(payload))
        digest.update(b']}')
        return 'history-' + stream + '-' + digest.hexdigest()

    def _validate_checkpoint(self, snapshot):
        for stream in ('evidence', 'assessments', 'events'):
            count = getattr(snapshot, stream + '_count')
            if self._prefix_digest(stream, count) != getattr(snapshot, stream + '_digest'):
                raise ValueError('checkpoint committed prefix hash mismatch: ' + stream)
        prior = self._indexes['snapshots'].get(snapshot.prior_snapshot_id) if snapshot.prior_snapshot_id else None
        if snapshot.prior_snapshot_id and (prior is None or prior.completeness != 'complete'):
            raise ValueError('snapshot inheritance requires complete prior snapshot')
        if prior and any(getattr(snapshot, s + '_count') < getattr(prior, s + '_count') for s in ('evidence', 'assessments', 'events')):
            raise ValueError('checkpoint prefixes precede inherited snapshot')
        for kind, source_stream, source_count, target_count in (
                ('evidence', 'assessments', snapshot.assessments_count, snapshot.evidence_count),
                ('assessment', 'events', snapshot.events_count, snapshot.assessments_count),
                ('parent', 'events', snapshot.events_count, snapshot.events_count)):
            if self._db.execute('''SELECT 1 FROM record_links l JOIN records s ON s.id=l.source_id
                    JOIN records t ON t.id=l.target_id WHERE l.kind=? AND s.stream=?
                    AND s.ordinal<? AND t.ordinal>=? LIMIT 1''', (kind, source_stream, source_count, target_count)).fetchone():
                raise ValueError('checkpoint references uncommitted prefix')
        if snapshot.completeness == 'complete':
            inherited = self._resolved_dispositions(prior) if prior else {}
            for (entity,) in self._db.execute("SELECT DISTINCT entity_id FROM entity_records WHERE stream='events' AND ordinal<?", (snapshot.events_count,)):
                if entity not in snapshot.dispositions and entity not in inherited:
                    raise ValueError('complete snapshot requires disposition for every generated/explored entity')
            for (entity,) in self._db.execute("SELECT DISTINCT entity_id FROM entity_records WHERE stream='events' AND ordinal>=? AND ordinal<?", (prior.events_count if prior else 0, snapshot.events_count)):
                if entity not in snapshot.dispositions:
                    raise ValueError('complete snapshot requires fresh dispositions for touched entities')

    def _append(self, stream, record):
        existing = self._indexes[stream].get(record.id)
        if existing is not None:
            if existing != record:
                raise ValueError('record identity collision')
            return existing
        ordinal = self._counts[stream]
        self._validate_disk_record(stream, record, ordinal)
        raw = canonical_json(record).encode()
        try:
            self._db.execute('INSERT INTO records(stream,ordinal,id,payload,stage_id,pool,profile_id) VALUES (?,?,?,?,?,?,?)',
                             (stream, ordinal, record.id, zlib.compress(raw, level=1), getattr(record, 'stage_id', None),
                              getattr(record, 'pool', None), getattr(record, 'profile_id', None)))
            self._insert_relations(stream, ordinal, record)
            self._update_hash(stream, raw)
            self._counts[stream] += 1
            self._pending += 1
            if self._pending >= self.batch_size or stream == 'snapshots':
                self.checkpoint()
        except BaseException:
            # Cancellation also must not leave a partial record for close().
            # A partially written record must never survive into a later commit.
            # Restore hashes/counts along with the entire pending SQL batch.
            self.rollback()
            raise
        return record

    @property
    def evidence(self): return self._records['evidence']

    @property
    def assessments(self): return self._records['assessments']

    @property
    def events(self): return self._records['events']

    @property
    def snapshots(self): return self._records['snapshots']

    def latest_event(self, entity_id):
        row = self._db.execute("SELECT r.payload FROM entity_records e JOIN records r ON r.id=e.record_id WHERE e.entity_id=? AND e.stream='events' ORDER BY e.ordinal DESC LIMIT 1", (entity_id,)).fetchone()
        return self._decode('events', row[0]) if row else None

    def iter_events(self, start=0, stop=None):
        if type(start) is not int or start < 0 or (stop is not None and (type(stop) is not int or stop < start)):
            raise ValueError('event range requires nonnegative ordered offsets')
        cursor = self._db.execute("SELECT payload FROM records WHERE stream='events' AND ordinal>=? AND ordinal<? ORDER BY ordinal",
                                  (start, self._counts['events'] if stop is None else stop))
        return (self._decode('events', payload) for (payload,) in cursor)

    def complete_stage(self, *, stage_id, dispositions, catalog_digest, ledger_digest,
                       prior_snapshot_id=None, completeness='complete'):
        return self._append('snapshots', StageSnapshot(stage_id, completeness, prior_snapshot_id, dispositions,
            catalog_digest, ledger_digest, *(self._prefix_digest(s, self._counts[s]) for s in ('evidence', 'assessments', 'events')),
            *(self._counts[s] for s in ('evidence', 'assessments', 'events'))))

    def checkpoint(self):
        """Durably commit all current records, including an incomplete stage tail."""
        try:
            self._db.execute("UPDATE metadata SET value=? WHERE key='committed_counts'", (canonical_json(self._counts),))
            self._db.commit()
            self._pending = 0
        except BaseException:
            self.rollback()
            raise

    def rollback(self):
        """Discard the pending batch and restore the last durable prefix.

        Entity row IDs may be reused after rollback, so no cached key survives.
        """
        self._entity_key_cache.clear()
        self._db.rollback()
        self._pending = 0
        self._counts = {name: 0 for name in self._streams}
        self._hashes = {name: _history_hash_start(name) for name in self._streams}
        self._reload()

    def close(self):
        if not self._closed:
            try:
                self.checkpoint()
            finally:
                self._entity_key_cache.clear()
                self._db.close()
                self._closed = True

    def __enter__(self): return self

    def __exit__(self, *exc): self.close()

    def query(self, *, entity_id=None, pool=None, stage_id=None, profile_id=None, include_lineage=False):
        def filtered(stream, *, assessment=False):
            terms, args = ['r.stream=?'], [stream]
            if entity_id is not None:
                terms.append('EXISTS(SELECT 1 FROM entity_records e WHERE e.record_id=r.id AND e.entity_id=?)')
                args.append(entity_id)
            for column, value in (('stage_id', stage_id), ('pool', pool if assessment else None),
                                  ('profile_id', profile_id if assessment else None)):
                if value is not None:
                    terms.append('r.' + column + '=?')
                    args.append(value)
            return 'SELECT r.id FROM records r WHERE ' + ' AND '.join(terms), args
        aq, ap = filtered('assessments', assessment=True)
        eq, ep = filtered('events')
        if include_lineage:
            eq = '''WITH RECURSIVE walk(id) AS (''' + eq + ''' UNION
                SELECT CASE WHEN l.source_id=w.id THEN l.target_id ELSE l.source_id END
                FROM record_links l JOIN walk w ON l.source_id=w.id OR l.target_id=w.id
                WHERE l.kind='parent') SELECT id FROM walk'''
        assessments = tuple(self._decode('assessments', p) for (p,) in self._db.execute(
            'SELECT payload FROM records WHERE id IN (' + aq + ') ORDER BY ordinal', ap))
        events = tuple(self._decode('events', p) for (p,) in self._db.execute(
            'SELECT payload FROM records WHERE id IN (' + eq + ') ORDER BY ordinal', ep))
        evidence = tuple(self._decode('evidence', p) for (p,) in self._db.execute('''SELECT r.payload FROM records r
            WHERE r.stream='evidence' AND (r.id IN (SELECT target_id FROM record_links WHERE kind='evidence'
            AND source_id IN (''' + aq + ''')) OR EXISTS(SELECT 1 FROM entity_records e WHERE e.record_id=r.id
            AND e.entity_id=?)) ORDER BY r.ordinal''', [*ap, entity_id]))
        snapshots = tuple(s for s in self.snapshots if stage_id is None or s.stage_id == stage_id)
        return dict(evidence=evidence, assessments=assessments, events=events, snapshots=snapshots)

    def export_database(self, path: Path):
        """Portable, self-contained SQLite backup; no WAL sidecar required."""
        path = Path(path)
        if path.resolve() == self.database_path.resolve():
            raise ValueError('backup path must differ from active database')
        path.parent.mkdir(parents=True, exist_ok=True)
        self.checkpoint()
        destination = sqlite3.connect(path)
        try:
            self._db.backup(destination)
        finally:
            destination.close()
        return _file_sha256(path)

    def export(self, directory: Path):
        """Stream deterministic gzip JSONL exports without materializing streams."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.checkpoint()
        hashes = {}
        for stream in self._streams:
            path = directory / (stream + '.jsonl.gz')
            with path.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as output:
                for (payload,) in self._db.execute('SELECT payload FROM records WHERE stream=? ORDER BY ordinal', (stream,)):
                    output.write(zlib.decompress(payload) + b'\n')
            hashes[path.name] = _file_sha256(path)
        path = directory / 'snapshot-index.json'
        with path.open('w') as output:
            output.write('{')
            for i, (identity,) in enumerate(self._db.execute("SELECT id FROM records WHERE stream='snapshots' ORDER BY id")):
                snapshot = self._indexes['snapshots'][identity]
                if i:
                    output.write(',')
                output.write(canonical_json(identity) + ':' + canonical_json(dict(stage_id=snapshot.stage_id,
                    completeness=snapshot.completeness, dispositions=self._resolved_dispositions(snapshot))))
            output.write('}')
        hashes[path.name] = _file_sha256(path)
        return hashes


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()
