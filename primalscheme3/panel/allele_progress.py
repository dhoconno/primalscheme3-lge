"""Observational search progress, separate from canonical scientific history."""

import json
from datetime import UTC, datetime

from .coverage_types import semantic_id

SCHEMA = "primalscheme3.allele-search-progress/v1"


class ProgressWriter:
    """Exclusive, append-only JSONL sink; flushed records survive workflow failure."""

    def __init__(self, path):
        self.stream = open(path, "x", encoding="utf-8")

    def __call__(self, record):
        record = dict(
            record,
            schemaVersion=SCHEMA,
            timestamp=datetime.now(UTC).isoformat(),
        )
        self.stream.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
        self.stream.flush()

    def close(self):
        self.stream.close()


class SearchProgress:
    """Uses existing search clock samples; never evaluates compatibility predicates."""

    def __init__(self, search, proposals, stage, observer, begin):
        self.search, self.proposals = search, proposals
        self.stage, self.observer, self.begin = stage, observer, begin
        self.now = self.last_emit = begin
        self.phase = None
        self.latest_improvement = None
        self.capture(search.state(search.best))
        self.emit("search-started", assignments=True)

    def capture(self, state):
        self.assignments = [
            dict(configuration_id=a.candidate_id, pool=a.pool) for a in self.search.best
        ]
        self.incumbent = semantic_id(
            "incumbent",
            tuple(sorted((a.candidate_id, a.pool) for a in self.search.best)),
        )
        self.coverage = {t.id: state._fraction(t.id) for t in state.catalog.targets}
        values = [v for v in self.coverage.values() if v is not None]
        self.mean = sum(values) / len(values) if values else None

    def emit(self, kind, *, assignments=False, **extra):
        record = dict(
            kind=kind,
            stage=self.stage,
            elapsed_seconds=max(0, self.now - self.begin),
            active_phase=self.phase,
            incumbent_id=self.incumbent,
            objective=list(self.search.best_key[:-1]),
            objective_version="soft-observed-allele-coverage/v1",
            per_target_coverage=dict(self.coverage),
            mean_coverage=self.mean,
            work=dict(self.search.work),
            proposal_work=dict(self.proposals.work),
            family_cursors=dict(self.proposals.cursors),
            unique_expanded_families=len(self.proposals.searched_families),
            latest_improvement_seconds=self.latest_improvement,
            validation_scope="provisional-oracle-checked",
            **extra,
        )
        if assignments:
            record["assignments"] = list(self.assignments)
        self.observer(record)
        self.last_emit = self.now

    def tick(self, now):
        self.now = now
        if now - self.last_emit >= 30:
            self.emit("heartbeat")

    def improved(self, state):
        self.capture(state)
        self.latest_improvement = max(0, self.now - self.begin)
        self.emit("incumbent-improved", assignments=True)

    def phase_event(self, name, kind, now, **extra):
        self.now, self.phase = now, name
        self.emit(kind, **extra)
