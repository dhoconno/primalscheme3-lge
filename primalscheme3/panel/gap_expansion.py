"""Bounded, gap-focused legacy candidate expansion.

This module deliberately keeps the native legacy chemistry and pair geometry.
It only recovers row-confirmed single-oligo alternatives at anchors where the
ordinary cloud discovery stopped at its first compatible length or rejected a
whole cloud because one sibling failed chemistry.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from typing import Any

from primalschemers import FKmer, RKmer, do_pool_interact

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import AmpliconSizeMetric
from primalscheme3.core.digestion import variant_digest_index
from primalscheme3.panel.allele_validation import numerical_dimer_score


@dataclass(frozen=True)
class GapExpansionOptions:
    """Validated bounds for the opt-in expansion stage.

    ``max_anchors_per_msa`` counts one direction plus one alignment anchor.
    ``max_pairs_per_msa`` counts generated same-row pair attempts, including
    strict chemistry or geometry failures.
    """

    mode: str = "off"
    max_anchors_per_msa: int = 2000
    max_pairs_per_msa: int = 1000

    def __post_init__(self) -> None:
        if self.mode not in {"off", "bounded"}:
            raise ValueError("gap expansion mode must be off or bounded")
        for name in ("max_anchors_per_msa", "max_pairs_per_msa"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class GapExpansionResult:
    candidates: tuple[PrimerPair, ...]
    evidence: tuple[dict[str, Any], ...]
    report: dict[str, Any]


@dataclass(frozen=True)
class _Oligo:
    direction: str
    anchor: int
    row_index: int
    sequence: str
    footprint: tuple[int, ...]


def _gaps(coverage: dict[str, dict[str, Any]]) -> list[tuple[int, int, int]]:
    result = []
    for target, summary in sorted(coverage.items(), key=lambda item: int(item[0])):
        cursor = 0
        for start, end in summary["trimmedIntervals"]:
            if cursor < start:
                result.append((int(target), cursor, int(start)))
            cursor = max(cursor, int(end))
        length = int(summary["referenceLength"])
        if cursor < length:
            result.append((int(target), cursor, length))
    return result


def _mapped_span(msa, footprint: tuple[int, ...]) -> tuple[int, int] | None:
    mapping = getattr(msa, "_mapping_array", None)
    if mapping is None or not footprint:
        return None
    mapped = []
    for column in footprint:
        if not 0 <= int(column) < len(mapping) or mapping[int(column)] is None:
            return None
        mapped.append(int(mapping[int(column)]))
    if mapped != list(range(mapped[0], mapped[-1] + 1)):
        return None
    return mapped[0], mapped[-1] + 1


def _spatial_order(values: list[int]) -> list[int]:
    """Visit a sorted anchor list by midpoint subdivision, not prefix order."""
    ordered: list[int] = []
    pending = deque([(0, len(values))])
    while pending:
        start, stop = pending.popleft()
        if start >= stop:
            continue
        middle = (start + stop - 1) // 2
        ordered.append(values[middle])
        pending.append((start, middle))
        pending.append((middle + 1, stop))
    return ordered


def _make_kmer(msa, oligo: _Oligo):
    span = _mapped_span(msa, oligo.footprint)
    if span is None or span[1] - span[0] != len(oligo.sequence):
        return None
    if oligo.direction == "f":
        kmer = FKmer([oligo.sequence.encode()], oligo.anchor)
        kmer.remap(span[1])
        if kmer.region() != span:
            return None
    else:
        kmer = RKmer([oligo.sequence.encode()], oligo.anchor)
        kmer.remap(span[0])
        if kmer.region() != span:
            return None
    return kmer


def _anchor_tasks(msa, gaps: list[tuple[int, int, int]], config):
    mapping = getattr(msa, "_mapping_array", None)
    if mapping is None:
        return []
    by_gap: dict[int, dict[str, list[int]]] = {}
    for gap_id, (_msa_index, start, end) in enumerate(gaps):
        window_start = max(0, start - int(config.amplicon_size_max))
        window_end = end + int(config.amplicon_size_max)
        forward = [
            index
            for index in range(1, len(mapping) + 1)
            if mapping[index - 1] is not None
            and window_start <= int(mapping[index - 1]) + 1 <= window_end
        ]
        reverse = [
            index
            for index in range(len(mapping))
            if mapping[index] is not None
            and window_start <= int(mapping[index]) <= window_end
        ]
        by_gap[gap_id] = {
            "f": deque(_spatial_order(sorted(set(forward)))),
            "r": deque(_spatial_order(sorted(set(reverse)))),
        }
    # Round-robin gap windows and directions so a cap cannot spend all work on
    # one long terminal gap or on forward discovery alone.
    tasks = []
    seen: set[tuple[str, int]] = set()
    while True:
        added = False
        for gap_id in sorted(by_gap):
            for direction in ("f", "r"):
                values = by_gap[gap_id][direction]
                if not values:
                    continue
                anchor = values.popleft()
                key = (direction, anchor)
                if key in seen:
                    continue
                seen.add(key)
                tasks.append((gap_id, direction, anchor))
                added = True
        if not added and not any(by_gap[gap_id][direction] for gap_id in by_gap for direction in ("f", "r")):
            break
    return tasks


def _pair_geometry_passes(fkmer, rkmer, config) -> bool:
    fstart = min(fkmer.starts())
    if config.amplicon_size_metric == AmpliconSizeMetric.REFERENCE_SPAN:
        span = rkmer.region()[1] - fstart
        return config.amplicon_size_min <= span <= config.amplicon_size_max
    distance = rkmer.start - fstart
    return config.amplicon_size_min <= distance <= config.amplicon_size_max


def expand_gap_candidates(msa, primary_coverage, config, options: GapExpansionOptions) -> GapExpansionResult:
    """Discover bounded row-confirmed candidates for one MSA occurrence."""
    if options.mode == "off":
        return GapExpansionResult((), (), {"mode": "off", "status": "disabled"})
    gaps = [gap for gap in _gaps(primary_coverage) if gap[0] == int(msa.msa_index)]
    if not gaps:
        return GapExpansionResult((), (), {"mode": options.mode, "status": "no-gaps"})
    tasks = _anchor_tasks(msa, gaps, config)
    attempted_anchor_tasks = tasks[: options.max_anchors_per_msa]
    statuses = Counter()
    oligos: dict[str, dict[int, list[_Oligo]]] = {"f": defaultdict(list), "r": defaultdict(list)}
    for _gap_id, direction, anchor in attempted_anchor_tasks:
        records = variant_digest_index(
            msa.array,
            config,
            direction,
            anchor,
            expansion_limit=1,
            length_mode="all",
        )
        for record in records:
            sequence = record.get("sequence")
            if sequence is None:
                statuses[record.get("reason", "no-sequence")] += 1
                continue
            if record.get("support") != "confirmed":
                statuses["ambiguous-support-excluded"] += 1
                continue
            if not record.get("accepted"):
                failed_checks = [
                    name
                    for name, check in record.get("checks", {}).items()
                    if check.get("outcome") == "fail"
                ]
                if failed_checks:
                    statuses.update(f"chemistry-{name}" for name in failed_checks)
                else:
                    statuses["chemistry-rejected"] += 1
                continue
            footprint = tuple(int(column) for column in record["alignment_footprint"])
            if _mapped_span(msa, footprint) is None or len(footprint) != len(sequence):
                statuses["unrepresentable-footprint"] += 1
                continue
            self_score = record.get("checks", {}).get("self-dimer", {}).get("value")
            if self_score is None:
                self_score = numerical_dimer_score(sequence, sequence).score
            if float(self_score) <= config.dimer_score:
                statuses["strict-self-dimer"] += 1
                continue
            oligos[direction][anchor].append(
                _Oligo(direction, anchor, int(record["row_index"]), sequence, footprint)
            )
            statuses["oligos-accepted"] += 1
    for direction in oligos:
        for anchor in oligos[direction]:
            unique = {(item.row_index, item.sequence, item.footprint): item for item in oligos[direction][anchor]}
            oligos[direction][anchor] = [unique[key] for key in sorted(unique)]

    # Pair attempts are generated per gap with a quota, then interleaved. This
    # keeps pair-cap coverage diverse without materializing a global graph.
    per_gap_quota = max(1, math.ceil(options.max_pairs_per_msa / len(gaps)))
    pair_tasks: dict[int, deque[tuple[_Oligo, _Oligo]]] = defaultdict(deque)
    pair_search_truncated: dict[str, bool] = {}
    pair_task_seen: set[tuple[str, str, int, int]] = set()
    for gap_id, (_msa_index, start, end) in enumerate(gaps):
        made = 0
        truncated = False
        families: list[deque[tuple[_Oligo, _Oligo]]] = []
        for f_anchor in _spatial_order(sorted(oligos["f"])):
            for r_anchor in _spatial_order(sorted(oligos["r"])):
                if len(families) >= per_gap_quota:
                    truncated = True
                    break
                f_ref_end = int(msa._mapping_array[f_anchor - 1]) + 1
                r_ref_start = int(msa._mapping_array[r_anchor])
                if f_ref_end > end or r_ref_start < start:
                    continue
                f_by_row = defaultdict(list)
                r_by_row = defaultdict(list)
                for item in oligos["f"][f_anchor]:
                    f_by_row[item.row_index].append(item)
                for item in oligos["r"][r_anchor]:
                    r_by_row[item.row_index].append(item)
                family: deque[tuple[_Oligo, _Oligo]] = deque()
                for row_index in sorted(set(f_by_row) & set(r_by_row)):
                    for f_item in f_by_row[row_index]:
                        for r_item in r_by_row[row_index]:
                            fkmer = _make_kmer(msa, f_item)
                            rkmer = _make_kmer(msa, r_item)
                            if fkmer is None or rkmer is None:
                                statuses["unrepresentable-footprint"] += 1
                                continue
                            if not _pair_geometry_passes(fkmer, rkmer, config):
                                statuses["amplicon-size-rejected"] += 1
                                continue
                            pair_key = (
                                f_item.sequence,
                                r_item.sequence,
                                f_item.anchor,
                                r_item.anchor,
                            )
                            if pair_key in pair_task_seen:
                                continue
                            pair_task_seen.add(pair_key)
                            family.append((f_item, r_item))
                if family:
                    families.append(family)
            if len(families) >= per_gap_quota:
                break
        # One member from each anchor-pair family per round prevents a single
        # early family with many row/length variants consuming the gap quota.
        while families and made < per_gap_quota:
            next_families: list[deque[tuple[_Oligo, _Oligo]]] = []
            for family in families:
                if not family or made >= per_gap_quota:
                    continue
                pair_tasks[gap_id].append(family.popleft())
                made += 1
                if family:
                    next_families.append(family)
            families = next_families
        if families:
            truncated = True
        pair_search_truncated[str(gap_id)] = truncated
    candidates: list[PrimerPair] = []
    evidence: list[dict[str, Any]] = []
    attempted_pairs = 0
    seen_pairs: set[tuple[str, str, int, int]] = set()
    while any(pair_tasks.values()) and attempted_pairs < options.max_pairs_per_msa:
        for gap_id in sorted(pair_tasks):
            if attempted_pairs >= options.max_pairs_per_msa:
                break
            if not pair_tasks[gap_id]:
                continue
            f_item, r_item = pair_tasks[gap_id].popleft()
            key = (f_item.sequence, r_item.sequence, f_item.anchor, r_item.anchor)
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            attempted_pairs += 1
            fkmer = _make_kmer(msa, f_item)
            rkmer = _make_kmer(msa, r_item)
            if fkmer is None or rkmer is None:
                statuses["unrepresentable-footprint"] += 1
                continue
            if do_pool_interact(fkmer.seqs_bytes(), rkmer.seqs_bytes(), config.dimer_score):
                statuses["strict-pair-dimer"] += 1
                continue
            pair = PrimerPair(fkmer, rkmer, int(msa.msa_index))
            pair.chrom_name = msa._chrom_name
            pair.amplicon_prefix = msa._uuid
            candidates.append(pair)
            row_ids = list(getattr(msa, "_seq_dict", {}).keys())
            evidence.append(
                {
                    "rowIndex": f_item.row_index,
                    "rowId": str(row_ids[f_item.row_index]) if f_item.row_index < len(row_ids) else None,
                    "forwardFootprint": list(f_item.footprint),
                    "reverseFootprint": list(r_item.footprint),
                    "origin": "expanded-lengths-and-confirmed-single-oligos",
                }
            )
            statuses["pairs-accepted"] += 1
    report = {
        "mode": options.mode,
        "maxAnchorsPerMsa": options.max_anchors_per_msa,
        "maxPairsPerMsa": options.max_pairs_per_msa,
        "anchorWorkUnit": "direction-plus-alignment-anchor",
        "pairWorkUnit": "same-observed-row-forward-reverse-attempt",
        "gapCount": len(gaps),
        "anchorTasksAvailable": len(tasks),
        "anchorTasksEvaluated": len(attempted_anchor_tasks),
        "anchorTasksNotExplored": max(0, len(tasks) - len(attempted_anchor_tasks)),
        "pairAttempts": attempted_pairs,
        "pairChecksNotExplored": sum(len(queue) for queue in pair_tasks.values()),
        "pairSearchTruncatedByGap": pair_search_truncated,
        "statusCounts": dict(sorted(statuses.items())),
        "candidateCount": len(candidates),
    }
    return GapExpansionResult(tuple(candidates), tuple(evidence), report)
