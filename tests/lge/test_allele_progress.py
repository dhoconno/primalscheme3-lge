import json

import pytest
from test_allele_search import GraphOracle, instance

from primalscheme3.panel.allele_search import (
    AlleleSearchOptions,
    search_allele_configurations,
)
from primalscheme3.panel.coverage_history import CoverageHistory
from primalscheme3.panel.coverage_types import ConfigurationLedger


def execute(observer=None, **kwargs):
    cat, configs, profile = instance({"a": (2, 10), "b": (12, 22)})
    history = CoverageHistory(run_id="equivalence")
    result = search_allele_configurations(
        cat,
        profile,
        AlleleSearchOptions(starts=1, repair_rounds=0, time_limit=1000),
        ConfigurationLedger(cat.semantic_digest, tuple(configs.values())),
        GraphOracle(configs.values()),
        history=history,
        observer=observer,
        **kwargs,
    )
    return result


def test_observation_preserves_assignments_objective_work_and_canonical_history():
    events = []
    plain, watched = execute(), execute(events.append)
    assert plain.assignments == watched.assignments
    assert plain.objective == watched.objective
    assert plain.metadata["work"] == watched.metadata["work"]
    for stream in ("evidence", "assessments", "events", "snapshots"):
        assert getattr(plain.history, stream) == getattr(watched.history, stream)
    assert events[0]["kind"] == "search-started"
    assert events[-1]["kind"] == "search-finished"
    assert any(e["kind"] == "incumbent-improved" for e in events)
    assert any(e["kind"] == "phase-started" for e in events)
    assert all(e["validation_scope"] == "provisional-oracle-checked" for e in events)
    assert events[-1]["objective"] == list(watched.objective[:-1])


def test_jsonl_writer_flushes_and_preserves_partial_failure(tmp_path):
    from primalscheme3.panel.allele_progress import ProgressWriter

    path = tmp_path / "progress.jsonl"
    writer = ProgressWriter(path)
    writer({"kind": "partial"})
    assert json.loads(path.read_text())["kind"] == "partial"
    writer.close()
    with pytest.raises(FileExistsError):
        ProgressWriter(path)


@pytest.mark.parametrize("ending", ["failure", "cancelled"])
def test_cooperative_ticks_throttled_and_terminal_events(monkeypatch, ending):
    from test_phase_scheduling import Clock

    from primalscheme3.panel.coverage_search import _Cancelled, _Search

    clock, records = Clock(), []

    def fill(search, *args):
        for _ in range(65):
            clock.advance()
            search.tick()
        if ending == "failure":
            raise RuntimeError("injected")
        raise _Cancelled

    monkeypatch.setattr(_Search, "fill", fill)
    if ending == "failure":
        with pytest.raises(RuntimeError, match="injected"):
            execute(records.append, clock=clock)
        assert records[-1]["kind"] == "search-failed"
    else:
        result = execute(records.append, clock=clock)
        assert result.metadata["stop_reason"] == "cancelled"
        assert records[-1]["outcome"] == "cancelled"
    assert [r["elapsed_seconds"] for r in records if r["kind"] == "heartbeat"] == [
        30,
        60,
    ]
    assert records[-2]["kind"] == "phase-finished"


def test_empty_assessable_coverage_uses_null():
    cat, configs, profile = instance({}, rows=("N" * 24,))
    records = []
    search_allele_configurations(
        cat,
        profile,
        AlleleSearchOptions(starts=1),
        ConfigurationLedger(cat.semantic_digest, tuple(configs.values())),
        GraphOracle(configs.values()),
        observer=records.append,
    )
    assert records[0]["mean_coverage"] is None
    assert list(records[0]["per_target_coverage"].values()) == [None]
    json.dumps(records, allow_nan=False)


def test_observer_does_not_add_search_clock_calls():
    class Clock:
        calls = 0

        def __call__(self):
            self.calls += 1
            return self.calls / 1000

    a, b = Clock(), Clock()
    plain, watched = execute(clock=a), execute(lambda record: None, clock=b)
    assert a.calls == b.calls
    assert plain.assignments == watched.assignments
    assert plain.history.events == watched.history.events
