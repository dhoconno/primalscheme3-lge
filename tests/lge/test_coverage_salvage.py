"""Tier orchestration preserves strict results and uses cumulative pool guards."""

import pytest
from test_allele_validation import dna, fixture, profile

from primalscheme3.panel.allele_search import (
    AlleleSearchOptions,
    search_allele_assignments,
)


def strict_result(seed=91):
    cat, configs, ledger = fixture(seq=dna(seed=seed))
    return search_allele_assignments(
        cat,
        profile(),
        AlleleSearchOptions(starts=1, repair_rounds=0, time_limit=10, seed_modes=()),
        ledger,
    )


def test_disabled_does_not_mutate_history_or_strict_result():
    from primalscheme3.panel.coverage_salvage import SalvageOptions, run_salvage

    strict = strict_result()
    before = len(strict.history.events)
    result = run_salvage(strict, profile(), SalvageOptions())
    assert result.strict is strict
    assert result.tiers == ()
    assert len(strict.history.events) == before


@pytest.mark.parametrize(
    "kw",
    [
        {"thresholds": (-25,)},
        {"thresholds": (-28, -27)},
        {"thresholds": (-28, -28)},
        {"thresholds": (float("nan"),)},
        {"max_stages": 4},
        {"max_edges_per_pool": True},
        {"time_limit": 0},
    ],
)
def test_invalid_policy_rejected(kw):
    from primalscheme3.panel.coverage_salvage import SalvageOptions

    with pytest.raises(ValueError):
        SalvageOptions(mode="bounded", **kw)


def test_rescue_and_strict_rejections_both_remain_queryable():
    from primalscheme3.panel.coverage_salvage import (
        SalvageOptions,
        run_salvage,
        select_primary,
    )

    strict = strict_result()
    assert strict.assignments == ()
    run = run_salvage(
        strict,
        profile(),
        SalvageOptions(
            mode="bounded",
            thresholds=(-1000,),
            max_edges_per_pool=8,
            max_oligos_per_pool=4,
            time_limit=10,
        ),
    )
    assert len(run.tiers) == 1
    rescued = run.tiers[0].result
    assert rescued.validation["valid"] and rescued.assignments
    assert all(a.stage_id == "salvage-1" for a in rescued.assignments)
    assert strict.assignments == ()
    assert run.tiers[0].comparison["mean_coverage_delta"] > 0
    assert select_primary(run, "strict") is strict
    assert select_primary(run, "salvage-1") is rescued
    with pytest.raises(ValueError):
        select_primary(run, "salvage-2")
    stages = {a.stage_id for a in strict.history.assessments}
    assert {"strict", "salvage-1"} <= stages


def test_cumulative_exposure_cap_can_prevent_rescue_despite_permissive_score():
    from primalscheme3.panel.coverage_salvage import SalvageOptions, run_salvage

    strict = strict_result()
    run = run_salvage(
        strict,
        profile(),
        SalvageOptions(
            mode="bounded",
            thresholds=(-1000,),
            max_edges_per_pool=0,
            max_oligos_per_pool=0,
            time_limit=10,
        ),
    )
    assert not run.tiers[0].result.assignments
    assert run.tiers[0].result.validation["valid"]
    assert run.tiers[0].comparison["mean_coverage_delta"] == 0


def test_existing_strict_coverage_retained_and_tier_callback_failure_cannot_replace_it():
    from primalscheme3.panel.coverage_salvage import (
        SalvageOptions,
        run_salvage,
        select_primary,
    )

    strict = strict_result(seed=0)
    saved = tuple(strict.assignments)

    def fail(result):
        raise OSError("simulated publication failure")

    run = run_salvage(
        strict,
        profile(),
        SalvageOptions(mode="bounded", thresholds=(-28,), time_limit=10),
        on_tier=fail,
    )
    assert run.tiers[0].status == "failed"
    assert (
        run.tiers[0].result is not None
    )  # Retain explored ledger for diagnostic history.
    assert "publication failure" in run.tiers[0].error
    assert select_primary(run, "strict") is strict
    with pytest.raises(ValueError):
        select_primary(run, "salvage-1")
    assert strict.assignments == saved


def test_salvage_cannot_change_strict_specificity_profile():
    from primalscheme3.panel.coverage_salvage import SalvageOptions, run_salvage

    strict = strict_result()
    with pytest.raises(ValueError, match="profile"):
        run_salvage(strict, profile(199), SalvageOptions(mode="bounded"))


def test_failure_disposition_does_not_survive_as_selected_in_cancelled_next_tier():
    from primalscheme3.panel.coverage_salvage import SalvageOptions, run_salvage

    strict = strict_result()
    cancel = False
    emit = strict.history.emit

    def wrapped(**kw):
        nonlocal cancel
        if kw["stage_id"] == "salvage-2" and kw["kind"] == "salvage-stage-start":
            cancel = True
        return emit(**kw)

    strict.history.emit = wrapped

    def publish(result):
        if result.metadata["stage_policy"]["stage_id"] == "salvage-1":
            raise OSError("failed first publication")

    run = run_salvage(
        strict,
        profile(),
        SalvageOptions(mode="bounded", thresholds=(-1000, -1001)),
        on_tier=publish,
        cancelled=lambda: cancel,
    )
    assert run.tiers[0].status == "failed"
    assert (
        run.tiers[0].result is not None
    )  # Retain explored ledger for diagnostic history.
    assert not run.tiers[1].result.assignments
    statuses = strict.history._resolved_dispositions(strict.history.last_complete_stage)
    assert "selected-salvage" not in statuses.values()


def test_cancellation_before_first_tier_has_durable_reason():
    from primalscheme3.panel.coverage_salvage import SalvageOptions, run_salvage

    strict = strict_result()
    run = run_salvage(
        strict, profile(), SalvageOptions(mode="bounded"), cancelled=lambda: True
    )
    assert run.stop_reason == "cancelled"
    assert any(e.kind == "salvage-run-stopped" for e in strict.history.events)
