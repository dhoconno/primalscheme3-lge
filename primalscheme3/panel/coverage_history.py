"""Append-only scientific evidence, contextual verdicts and decision lineage.

A single coordinator owns a history writer. Workers return evidence batches;
callers merge them deterministically before assigning event sequence numbers.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import gzip
import hashlib
import json
import os
from pathlib import Path
from typing import Mapping, Any
from primalscheme3.panel.coverage_types import ScientificRecord, canonical_json, semantic_id

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
