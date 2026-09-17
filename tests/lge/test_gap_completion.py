from __future__ import annotations

from types import SimpleNamespace

from primalschemers import FKmer, RKmer

from primalscheme3.core.classes import PrimerPair
from primalscheme3.panel.gap_completion import (
    FollowupSelection,
    gap_candidate_id,
    select_followup_candidates,
    trimmed_coverage,
)


def pair(start: int, end: int, forward: str = "A" * 20, reverse: str = "C" * 20):
    return PrimerPair(
        FKmer([forward.encode()], start + len(forward)),
        RKmer([reverse.encode()], end - len(reverse)),
        0,
    )


def target(length: int = 200):
    return SimpleNamespace(
        msa_index=0,
        _chrom_name="target",
        _mapping_array=list(range(length)),
    )


def test_trimmed_coverage_uses_ungapped_half_open_union():
    msa = target(120)
    incumbent = pair(0, 80)
    overlapping = pair(40, 120)
    summary = trimmed_coverage({0: msa}, [incumbent, overlapping])
    assert summary["0"]["trimmedIntervals"] == [[20, 100]]
    assert summary["0"]["trimmedBases"] == 80
    assert summary["0"]["referenceLength"] == 120


def test_selector_accepts_one_base_gain_against_parent_union():
    msa = target(100)
    parent = pair(0, 60)
    candidate = pair(40, 81)
    result = select_followup_candidates(
        {0: msa},
        [parent],
        [candidate],
        pool_count=2,
        admission=lambda _candidate, _pool: None,
    )
    assert isinstance(result, FollowupSelection)
    assert result.accepted[0].candidate_id == gap_candidate_id(candidate, msa)
    assert result.accepted[0].trimmed_gain == 1
    assert result.followup_coverage["0"]["trimmedBases"] == 1
    assert result.combined_coverage["0"]["trimmedBases"] == 21
    assert result.candidate_union_coverage["0"]["trimmedBases"] == 21


def test_selector_uses_deterministic_candidate_and_pool_ties():
    msa = target(160)
    first = pair(0, 80, "A" * 20, "C" * 20)
    second = pair(80, 160, "G" * 20, "T" * 20)
    calls = []

    def admission(candidate, pool):
        calls.append((gap_candidate_id(candidate, msa), pool))
        return None

    result = select_followup_candidates(
        {0: msa}, [], [second, first], pool_count=2, admission=admission
    )
    assert [item.candidate_id for item in result.accepted] == sorted(
        item.candidate_id for item in result.accepted
    )
    assert {item.pool for item in result.accepted} == {0}
    assert all(pool in {0, 1} for _, pool in calls)


def test_selector_records_invalid_geometry_without_sorting_error():
    msa = target(100)
    invalid = pair(90, 130)
    valid = pair(0, 80)
    result = select_followup_candidates(
        {0: msa}, [], [invalid, valid], pool_count=1, admission=lambda *_: None
    )
    assert result.statuses[gap_candidate_id(invalid, msa)]["reason"] == "geometry"
