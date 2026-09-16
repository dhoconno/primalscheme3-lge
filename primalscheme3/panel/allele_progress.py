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

    def __init__(self, search, proposals, stage, observer, begin, *, policy):
        self.search, self.proposals = search, proposals
        self.stage, self.observer, self.begin = stage, observer, begin
        self.policy = dict(policy)
        self.catalog_digest = proposals.catalog.semantic_digest
        self.interrupted_phase = None
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
            interrupted_phase=self.interrupted_phase,
            catalog_semantic_digest=self.catalog_digest,
            stage_policy=dict(self.policy),
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
            # Persist only the selected design definitions. Exact support is
            # reconstructed from the already retained catalog, not copied here.
            record["configuration_definitions"] = [
                {
                    key: getattr(self.proposals.records[a.candidate_id], key)
                    for key in (
                        "id",
                        "target_id",
                        "family_id",
                        "forward_site_ids",
                        "reverse_site_ids",
                        "full_interval",
                        "anchor_pair",
                    )
                }
                for a in self.search.best
            ]
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

    def emit_preserving_error(self, error, kind, **extra):
        if error is None:
            self.emit(kind, **extra)
            return
        try:
            self.emit(kind, **extra)
        except BaseException as observation_error:
            error.add_note(f"Progress {kind} failed: {observation_error}")

    def phase_event(self, name, kind, now, *, primary_error=None, **extra):
        self.now, self.phase = now, name
        if extra.get("outcome") in ("failed", "cancelled", "time-limit"):
            self.interrupted_phase = name
        try:
            self.emit_preserving_error(primary_error, kind, **extra)
        finally:
            if kind == "phase-finished":
                self.phase = None
