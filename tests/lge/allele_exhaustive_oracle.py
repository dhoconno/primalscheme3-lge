"""Tiny literal allele-set oracle; no production coverage/objective/state helpers."""

from itertools import product
from statistics import mean


def literal_objective(
    catalog, configurations, assignments, profile, goal=0.95, requested_size=None
):
    concrete = {}
    targets = {}
    for observation in catalog.observations:
        row = [cell for cell in observation.cells if cell not in ("", "-")]
        concrete[observation.id] = {i for i, cell in enumerate(row) if cell in "ACGT"}
        targets.setdefault(observation.target_id, []).append(observation.id)
    covered = {key: set() for key in concrete}
    species = [set() for _ in range(profile.n_pools)]
    deviation = 0
    desired = (
        requested_size
        if requested_size is not None
        else (profile.amplicon_size_min + profile.amplicon_size_max) / 2
    )
    for cid, pool in assignments:
        config = configurations[cid]
        for support in config.supported_products:
            covered[support.allele_id].update(
                range(support.forward_footprint[1], support.reverse_footprint[0])
            )
        species[pool].update(
            catalog.site_by_id[sid].sequence
            for sid in config.forward_site_ids + config.reverse_site_ids
        )
        deviation += abs(config.full_interval[1] - config.full_interval[0] - desired)
    fractions = []
    for tid in sorted(targets):
        assessable = [
            len(concrete[aid] & covered[aid]) / len(concrete[aid])
            for aid in targets[tid]
            if concrete[aid]
        ]
        if assessable:
            fractions.append(mean(assessable))
    utility = (
        mean(c - max(0, goal - c) ** 2 / goal for c in fractions) if fractions else 0
    )
    groups = sorted(
        tuple(sorted(cid for cid, p in assignments if p == pool))
        for pool in range(profile.n_pools)
        if any(p == pool for _, p in assignments)
    )
    canonical = tuple(
        sorted((cid, pool) for pool, group in enumerate(groups) for cid in group)
    )
    numeric = (
        -utility,
        -mean(fractions) if fractions else 0,
        -min(fractions) if fractions else 0,
        0,
        0,
        sum(map(len, species)),
        len(assignments),
        max(map(len, species), default=0) - min(map(len, species), default=0),
        deviation,
    )
    return tuple(round(x, 12) for x in numeric) + (canonical,)


def exhaustive(
    catalog,
    configurations,
    profile,
    *,
    edges=(),
    invalid=(),
    goal=0.95,
    requested_size=None,
):
    keys = tuple(sorted(configurations))
    forbidden = {frozenset(edge) for edge in edges}
    best = None
    for pools in product(range(-1, profile.n_pools), repeat=len(keys)):
        placed = tuple((cid, pool) for cid, pool in zip(keys, pools) if pool >= 0)
        if any(cid in invalid for cid, _ in placed):
            continue
        if profile.max_amplicons is not None and len(placed) > profile.max_amplicons:
            continue
        if profile.max_amplicons_msa is not None and any(
            sum(configurations[cid].target_id == t.id for cid, _ in placed)
            > profile.max_amplicons_msa
            for t in catalog.targets
        ):
            continue
        if any(
            p == q and frozenset((a, b)) in forbidden
            for i, (a, p) in enumerate(placed)
            for b, q in placed[i + 1 :]
        ):
            continue
        objective = literal_objective(
            catalog, configurations, placed, profile, goal, requested_size
        )
        if best is None or objective < best[0]:
            best = objective, placed
    return best
