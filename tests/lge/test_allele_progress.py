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


def test_completed_record_has_no_active_phase():
    records = []
    execute(records.append)
    assert records[-1]["active_phase"] is None


def test_terminal_observer_failure_preserves_scientific_error(monkeypatch):
    from primalscheme3.panel.coverage_search import _Search

    def fail(*args):
        raise ValueError("primary scientific failure")

    def observer(record):
        if record["kind"] in ("phase-finished", "search-failed"):
            raise OSError("secondary sink failure")

    monkeypatch.setattr(_Search, "fill", fail)
    with pytest.raises(ValueError, match="primary scientific failure") as exc:
        execute(observer)
    assert any("secondary sink failure" in note for note in exc.value.__notes__)


def test_observer_latency_is_charged_to_wall_budget():
    from test_phase_scheduling import Clock

    clock = Clock()

    def slow(record):
        clock.advance(1001)

    watched = execute(slow, clock=clock)
    assert watched.metadata["stop_reason"] == "time-limit"
    assert not watched.assignments


def test_generated_incumbent_failure_retains_self_resolving_definitions(tmp_path):
    import gzip

    from primalscheme3.panel.allele_progress import ProgressWriter
    from primalscheme3.panel.allele_validation import StagePolicy
    from primalscheme3.panel.coverage_types import (
        SelectionConfiguration,
        VariantCatalog,
    )
    from primalscheme3.panel.coverage_variants import make_configuration

    cat, _, profile = instance({"a": (2, 10), "b": (12, 22)})
    catalog_path = tmp_path / "discovery-catalog.json.gz"
    with gzip.open(catalog_path, "wt") as stream:
        json.dump(cat.to_dict(), stream)
    writer = ProgressWriter(tmp_path / "search-progress.jsonl")

    def observer(record):
        writer(record)
        if record["kind"] == "incumbent-improved":
            raise RuntimeError("failed after generated selection")

    with pytest.raises(RuntimeError, match="failed after generated selection"):
        search_allele_configurations(
            cat,
            profile,
            AlleleSearchOptions(starts=1),
            ConfigurationLedger(cat.semantic_digest, ()),
            GraphOracle(()),
            observer=observer,
        )
    writer.close()
    with gzip.open(catalog_path, "rt") as stream:
        retained = VariantCatalog.from_dict(json.load(stream))
    records = [
        json.loads(line)
        for line in (tmp_path / "search-progress.jsonl").read_text().splitlines()
    ]
    assert records[0]["stage_policy"]["stage_id"] == StagePolicy().stage_id
    improved = next(r for r in records if r["kind"] == "incumbent-improved")
    assert improved["catalog_semantic_digest"] == retained.semantic_digest
    definitions = {
        d["id"]: SelectionConfiguration.from_dict(d)
        for d in improved["configuration_definitions"]
    }
    for assignment in improved["assignments"]:
        config = definitions[assignment["configuration_id"]]
        rebuilt = make_configuration(
            retained, config.family_id, config.forward_site_ids, config.reverse_site_ids
        )
        assert rebuilt.id == config.id
        assert rebuilt.full_interval == config.full_interval
        assert rebuilt.anchor_pair == config.anchor_pair
    assert not (tmp_path / "configuration-ledger.json.gz").exists()


def test_salvage_progress_retains_explicit_stage_policy():
    from dataclasses import asdict

    from primalscheme3.panel.allele_validation import StagePolicy

    policy = StagePolicy("salvage-28", -28, -26, 8, 4)
    records = []
    execute(records.append, policy=policy)
    assert all(record["stage_policy"] == asdict(policy) for record in records)
