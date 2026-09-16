"""Cooperative, initial-cycle reservations; never extend the global deadline."""

import math
from copy import deepcopy

from .coverage_search import _Cancelled, _PhaseLimit, _TimeLimit

POLICIES = {
    "serial": {"id": "serial/v1"},
    "reserved": {
        "id": "initial-phase-reservations/v1",
        "initial_weights": {"seeds": 0.2, "construction": 0.4, "repair": 0.4},
        "repair_weights": {"preparation": 0.2, "cleanup": 0.2, "exchange": 0.6},
    },
}


def scheduling_policy(name):
    return deepcopy(POLICIES[name])


class Reservation:
    """Unused time flows to remaining weights; inactive phases have no weight."""

    def __init__(self, deadline, weight):
        self.deadline = deadline
        self.remaining = weight

    def cutoff(self, now, weight):
        grant = max(0.0, self.deadline - now) * weight / self.remaining
        self.remaining = max(0.0, self.remaining - weight)
        return min(self.deadline, now + grant)


class PhaseScheduler:
    def __init__(self, search, proposals, history, stage_id):
        self.search, self.proposals = search, proposals
        self.history, self.stage_id = history, stage_id
        self.epoch = search.clock()
        self.progress = []
        self.active_at_stop = None

    def run(self, name, action, reservation=None, weight=None):
        search = self.search
        begin = search.clock()
        prior = search.phase_deadline
        cutoff = prior
        if reservation is not None:
            cutoff = min(prior, reservation.cutoff(begin, weight))
        search.phase_deadline = cutoff
        before = dict(search.work)
        proposal_before = dict(self.proposals.work)
        item = {
            "phase": name,
            "outcome": "running",
            "started_seconds": max(0.0, begin - self.epoch),
            "local_budget_seconds": None
            if math.isinf(cutoff)
            else max(0.0, cutoff - begin),
        }
        self.progress.append(item)
        self.history.emit(
            stage_id=self.stage_id,
            kind="phase-started",
            entity_ids=(),
            changes={"phase": name},
        )
        try:
            search.tick()
            item["outcome"] = action() or "completed"
            # Indivisible calls may return after a cutoff without an internal tick.
            search.tick()
        except _PhaseLimit:
            item["outcome"] = "phase-time-limit"
        except (_TimeLimit, _Cancelled) as error:
            item["outcome"] = (
                "time-limit" if isinstance(error, _TimeLimit) else "cancelled"
            )
            if self.active_at_stop is None:
                self.active_at_stop = name
            raise
        except BaseException:
            item["outcome"] = "failed"
            raise
        finally:
            item["elapsed_seconds"] = max(0.0, search.clock() - begin)
            item["local_overshoot_seconds"] = max(
                0.0, begin + item["elapsed_seconds"] - cutoff
            )
            item["work_delta"] = {
                key: value - before.get(key, 0) for key, value in search.work.items()
            }
            item["proposal_work_delta"] = {
                key: value - proposal_before.get(key, 0)
                for key, value in self.proposals.work.items()
            }
            item["family_cursors"] = dict(self.proposals.cursors)
            self.history.emit(
                stage_id=self.stage_id,
                kind="phase-finished",
                entity_ids=(),
                changes={
                    key: item[key]
                    for key in (
                        "phase",
                        "outcome",
                        "work_delta",
                        "proposal_work_delta",
                        "family_cursors",
                    )
                },
            )
            search.phase_deadline = prior
        return item["outcome"]
