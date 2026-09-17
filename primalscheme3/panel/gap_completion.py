"""Independent follow-up pool selection for uncovered legacy reference gaps."""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from time import monotonic
from typing import Any

import dnaio

from primalscheme3.core.bedfiles import read_bedlines_to_bedprimerpairs
from primalscheme3.core.config import MappingType, TerminalGapPolicy
from primalscheme3.core.mapping import generate_reference
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex


@dataclass(frozen=True)
class FollowupCandidate:
    candidate_id: str
    pool: int
    trimmed_gain: int
    full_gain: int
    status: str = "accepted"
    reason: str | None = None


@dataclass(frozen=True)
class FollowupSelection:
    accepted: tuple[FollowupCandidate, ...]
    statuses: dict[str, dict[str, Any]]
    primary_coverage: dict[str, dict[str, Any]]
    followup_coverage: dict[str, dict[str, Any]]
    combined_coverage: dict[str, dict[str, Any]]
    candidate_union_coverage: dict[str, dict[str, Any]]


def _reference_length(msa) -> int:
    mapping = getattr(msa, "_mapping_array", None)
    if mapping is not None:
        observed = [int(value) for value in mapping if value is not None]
        return max(observed, default=-1) + 1
    array = getattr(msa, "array", ())
    return len(array[0]) if len(array) else 0


def _merge(intervals: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def _intervals(
    pair, msa, *, reference_length: int | None = None
) -> tuple[tuple[int, int], tuple[int, int]] | None:
    length = _reference_length(msa) if reference_length is None else reference_length
    trimmed = (int(pair.fprimer.end), int(pair.rprimer.start))
    full = (min(pair.fprimer.starts()), max(pair.rprimer.ends()))
    if not (0 <= trimmed[0] < trimmed[1] <= length):
        return None
    if not (0 <= full[0] < full[1] <= length):
        return None
    return trimmed, full


def _coverage(msa_dict, pairs, *, reference_lengths=None) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for index, msa in sorted(msa_dict.items()):
        trimmed: list[tuple[int, int]] = []
        full: list[tuple[int, int]] = []
        for pair in pairs:
            if int(pair.msa_index) != int(index):
                continue
            intervals = _intervals(
                pair, msa, reference_length=(reference_lengths or {}).get(index)
            )
            if intervals is None:
                raise ValueError(f"primer pair interval is outside MSA {index}")
            trimmed.append(intervals[0])
            full.append(intervals[1])
        trimmed_union = _merge(trimmed)
        full_union = _merge(full)
        length = (
            reference_lengths[index]
            if reference_lengths is not None and index in reference_lengths
            else _reference_length(msa)
        )
        trimmed_bases = sum(end - start for start, end in trimmed_union)
        full_bases = sum(end - start for start, end in full_union)
        result[str(index)] = {
            "referenceLength": length,
            "trimmedIntervals": [list(item) for item in trimmed_union],
            "fullIntervals": [list(item) for item in full_union],
            "trimmedBases": trimmed_bases,
            "fullBases": full_bases,
            "trimmedFraction": trimmed_bases / length if length else 0.0,
            "fullFraction": full_bases / length if length else 0.0,
        }
    return result


def trimmed_coverage(msa_dict, pairs):
    """Return coverage from primer-trimmed and full ungapped intervals."""
    return _coverage(msa_dict, pairs)


def _canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def gap_candidate_id(pair, msa, *, reference_length: int | None = None) -> str:
    intervals = _intervals(pair, msa, reference_length=reference_length)
    payload = {
        "msaIndex": int(msa.msa_index),
        "targetOccurrence": str(getattr(msa, "_chrom_name", msa.msa_index)),
        "trimmedInterval": list(intervals[0]) if intervals else None,
        "fullInterval": list(intervals[1]) if intervals else None,
        "forwardOligos": sorted(str(seq).upper() for seq in pair.fprimer.seqs()),
        "reverseOligos": sorted(str(seq).upper() for seq in pair.rprimer.seqs()),
    }
    return "gap-candidate-" + hashlib.sha256(_canonical_json(payload)).hexdigest()


def _gain(before, pair, msa, index: str) -> tuple[int, int] | None:
    intervals = _intervals(
        pair,
        msa,
        reference_length=before[index]["referenceLength"],
    )
    if intervals is None:
        return None
    trimmed, full = intervals
    trim_union = _merge([tuple(item) for item in before[index]["trimmedIntervals"]] + [trimmed])
    full_union = _merge([tuple(item) for item in before[index]["fullIntervals"]] + [full])
    return (
        sum(end - start for start, end in trim_union) - before[index]["trimmedBases"],
        sum(end - start for start, end in full_union) - before[index]["fullBases"],
    )


def select_followup_candidates(
    msa_dict,
    primary_pairs,
    candidates,
    *,
    pool_count: int,
    admission: Callable[[Any, int], str | None],
    on_accept: Callable[[Any, int], None] | None = None,
) -> FollowupSelection:
    """Select deterministic positive-gain additions into independent pools.

    ``admission`` is supplied by the native runner and sees follow-up pools
    only.  Parent pairs are used solely to seed the combined coverage union.
    """
    if type(pool_count) is not int or pool_count < 1:
        raise ValueError("gap-completion pool_count must be a positive integer")
    candidates = list(candidates)
    reference_lengths = {
        index: _reference_length(msa) for index, msa in msa_dict.items()
    }
    primary_coverage = _coverage(
        msa_dict, primary_pairs, reference_lengths=reference_lengths
    )
    candidate_union_coverage = _coverage(
        msa_dict,
        [
            *primary_pairs,
            *[
                pair
                for pair in candidates
                if _intervals(
                    pair,
                    msa_dict[pair.msa_index],
                    reference_length=reference_lengths[pair.msa_index],
                )
                is not None
            ],
        ],
    )
    followup_pairs: list[Any] = []
    statuses: dict[str, dict[str, Any]] = {}
    accepted: list[FollowupCandidate] = []
    remaining = []
    seen_ids: set[str] = set()
    for pair in candidates:
        msa = msa_dict[pair.msa_index]
        candidate_id = gap_candidate_id(
            pair, msa, reference_length=reference_lengths[pair.msa_index]
        )
        if candidate_id in seen_ids:
            continue
        seen_ids.add(candidate_id)
        remaining.append((candidate_id, pair, msa))

    while remaining:
        combined_pairs = [*primary_pairs, *followup_pairs]
        before = _coverage(
            msa_dict, combined_pairs, reference_lengths=reference_lengths
        )
        ranked = []
        for candidate_id, pair, msa in remaining:
            gain = _gain(before, pair, msa, str(pair.msa_index))
            if gain is None:
                ranked.append((0, "", candidate_id, pair, msa, None))
            else:
                ranked.append((-gain[0], candidate_id, candidate_id, pair, msa, gain))
        ranked.sort(key=lambda item: item[:3])
        chosen = None
        for _, _, candidate_id, pair, _msa, gain in ranked:
            if gain is None:
                statuses[candidate_id] = {"status": "rejected", "reason": "geometry"}
                continue
            trimmed_gain, full_gain = gain
            if trimmed_gain <= 0:
                statuses[candidate_id] = {
                    "status": "rejected",
                    "reason": "nonpositive-reference-gain",
                    "trimmedGain": trimmed_gain,
                    "fullGain": full_gain,
                }
                continue
            pool_reasons: dict[str, str] = {}
            for pool in range(pool_count):
                reason = admission(pair, pool)
                if reason is None:
                    chosen = (candidate_id, pair, pool, trimmed_gain, full_gain)
                    break
                pool_reasons[str(pool)] = reason
            if chosen is None:
                statuses[candidate_id] = {
                    "status": "rejected",
                    "reason": "admission",
                    "trimmedGain": trimmed_gain,
                    "fullGain": full_gain,
                    "pools": pool_reasons,
                }
            else:
                break
        if chosen is None:
            break
        candidate_id, pair, pool, trimmed_gain, full_gain = chosen
        followup_pairs.append(pair)
        if on_accept is not None:
            on_accept(pair, pool)
        accepted.append(FollowupCandidate(candidate_id, pool, trimmed_gain, full_gain))
        statuses[candidate_id] = {
            "status": "accepted",
            "pool": pool,
            "trimmedGain": trimmed_gain,
            "fullGain": full_gain,
        }
        remaining = [item for item in remaining if item[0] != candidate_id]

    followup_coverage = _coverage(
        msa_dict, followup_pairs, reference_lengths=reference_lengths
    )
    combined_coverage = _coverage(
        msa_dict, [*primary_pairs, *followup_pairs], reference_lengths=reference_lengths
    )
    return FollowupSelection(
        tuple(accepted),
        statuses,
        primary_coverage,
        followup_coverage,
        combined_coverage,
        candidate_union_coverage,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path, *, relative_to: Path | None = None) -> dict[str, Any]:
    return {
        "path": (
            path.relative_to(relative_to).as_posix() if relative_to is not None else str(path)
        ),
        "sha256": _sha256(path),
        "size": path.stat().st_size,
    }


def _remaining_gaps(coverage: dict[str, dict[str, Any]]) -> dict[str, list[list[int]]]:
    gaps: dict[str, list[list[int]]] = {}
    for index, summary in coverage.items():
        length = int(summary["referenceLength"])
        cursor = 0
        target_gaps: list[list[int]] = []
        for start, end in summary["trimmedIntervals"]:
            if cursor < start:
                target_gaps.append([cursor, start])
            cursor = max(cursor, end)
        if cursor < length:
            target_gaps.append([cursor, length])
        gaps[index] = target_gaps
    return gaps


def _parent_reference(parent_dir: Path) -> dict[str, str]:
    with dnaio.open(parent_dir / "reference.fasta") as records:
        return {
            str(getattr(record, "id", getattr(record, "name", ""))): record.sequence.upper()
            for record in records
        }


def _validate_parent(
    parent_dir: Path,
    msa_dict: dict[int, MSA],
    config,
    logger,
) -> tuple[list[Any], dict[str, Any]]:
    required = ("primer.bed", "primertrim.amplicon.bed", "reference.fasta")
    missing = [name for name in required if not (parent_dir / name).is_file()]
    if missing:
        raise ValueError("parent scheme is missing " + ", ".join(missing))
    parent_pairs, _headers = read_bedlines_to_bedprimerpairs(parent_dir / "primer.bed")
    msa_by_chrom = {msa._chrom_name: index for index, msa in msa_dict.items()}
    for pair in parent_pairs:
        if pair.chrom_name not in msa_by_chrom:
            raise ValueError(
                f"parent primer target {pair.chrom_name!r} is absent from supplied MSAs"
            )
        pair.msa_index = msa_by_chrom[pair.chrom_name]
    references = _parent_reference(parent_dir)
    supplied_targets = set(msa_by_chrom)
    if set(references) != supplied_targets:
        raise ValueError("parent reference targets do not match supplied MSA targets")
    reference_checks = {}
    for index, msa in msa_dict.items():
        expected = generate_reference(msa.array).upper()
        observed = references.get(msa._chrom_name)
        if observed != expected:
            raise ValueError(
                f"parent reference does not match supplied first-row reference for {msa._chrom_name}"
            )
        reference_checks[str(index)] = {
            "targetOccurrence": msa._chrom_name,
            "length": len(expected),
            "sha256": hashlib.sha256(expected.encode()).hexdigest(),
            "verified": True,
        }
    # Published legacy amplicon BEDs use the historical inclusive right edge
    # (one less than ``create_amplicon_str``'s half-open coordinate). Validate
    # target, name, pool, and primer-derived boundaries explicitly.
    observed_trimmed = []
    for line in (parent_dir / "primertrim.amplicon.bed").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise ValueError("parent trimmed amplicon BED has a malformed row")
        observed_trimmed.append(
            (fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]))
        )
    expected_trimmed = []
    for pair in sorted(parent_pairs, key=lambda pair: (pair.chrom_name, pair.amplicon_number)):
        expected_trimmed.append(
            (
                pair.chrom_name,
                int(pair.fprimer.region()[1]),
                int(pair.rprimer.region()[0]),
                f"{pair.amplicon_prefix}_{pair.amplicon_number}",
                int(pair.pool) + 1,
            )
        )
    if len(observed_trimmed) != len(expected_trimmed):
        raise ValueError("parent trimmed amplicon BED is inconsistent with parent primer BED")
    accepted_trimmed = {
        row
        for row in expected_trimmed
        for row in (row, (row[0], row[1], row[2] - 1, row[3], row[4]))
    }
    if any(row not in accepted_trimmed for row in observed_trimmed):
        raise ValueError("parent trimmed amplicon BED is inconsistent with parent primer BED")
    expected_keys = [(row[0], row[3], row[4]) for row in expected_trimmed]
    observed_keys = [(row[0], row[3], row[4]) for row in observed_trimmed]
    if sorted(observed_keys) != sorted(expected_keys):
        raise ValueError("parent trimmed amplicon BED is inconsistent with parent primer BED")
    for pair in parent_pairs:
        if _intervals(pair, msa_dict[pair.msa_index]) is None:
            raise ValueError("parent primer interval is outside supplied reference")
    receipt = parent_dir / "panel-provenance.json"
    verification = {
        "referenceAndPrimerBed": True,
        "verificationLevel": "reference-and-primer-bed",
        "referenceChecks": reference_checks,
        "parentReceiptPresent": receipt.is_file(),
    }
    if receipt.is_file():
        try:
            receipt_data = json.loads(receipt.read_text())
        except json.JSONDecodeError as error:
            raise ValueError("parent panel-provenance.json is invalid JSON") from error
        verification["parentReceiptStatus"] = receipt_data.get("status")
        verification["parentReceiptValidated"] = False
    logger.info("Validated parent scheme against supplied first-row references")
    return parent_pairs, verification


def _validate_followup_output(
    output: Path,
    msa_dict: dict[int, MSA],
    local_paths: list[Path],
    config,
    expected: FollowupSelection,
) -> dict[str, Any]:
    """Replay the published follow-up independently from selector state."""
    pairs, _headers = read_bedlines_to_bedprimerpairs(output / "primer.bed")
    by_chrom = {msa._chrom_name: index for index, msa in msa_dict.items()}
    for pair in pairs:
        if pair.chrom_name not in by_chrom:
            raise ValueError("follow-up primer BED contains an unknown MSA target")
        pair.msa_index = by_chrom[pair.chrom_name]
    replay_db = MatchDB(
        output / "work" / "fresh-validation-mismatch",
        [str(path) for path in local_paths] if config.use_matchdb else [],
        config,
    )
    replay = Multiplex(config, replay_db, msa_dict)
    for pair in sorted(
        pairs,
        key=lambda item: (int(item.pool), item.chrom_name, int(item.amplicon_number)),
    ):
        result = replay.check_primerpair_can_be_added(pair, int(pair.pool))
        if result.name != "OK":
            raise ValueError(
                f"fresh follow-up validation rejected {pair.amplicon_prefix}_{pair.amplicon_number}: {result.name}"
            )
        replay.add_primer_pair_to_pool(pair, int(pair.pool), int(pair.msa_index))
    actual_coverage = trimmed_coverage(msa_dict, pairs)
    expected_ids = {item.candidate_id for item in expected.accepted}
    actual_ids = {
        gap_candidate_id(pair, msa_dict[pair.msa_index]) for pair in pairs
    }
    if actual_ids != expected_ids:
        raise ValueError("fresh follow-up validation does not match selected candidate IDs")
    if actual_coverage != expected.followup_coverage:
        raise ValueError("fresh follow-up validation coverage differs from the selection snapshot")
    return {
        "status": "passed",
        "publishedPairCount": len(pairs),
        "publishedCandidateIds": sorted(actual_ids),
        "coverage": actual_coverage,
    }


def _assert_unchanged(path: Path, baseline: dict[str, Any]) -> None:
    current = _file_record(path)
    if current["sha256"] != baseline["sha256"] or current["size"] != baseline["size"]:
        raise ValueError(f"input changed during gap-completion run: {path}")


def run_gap_completion(
    *,
    msa: list[Path],
    output_dir: Path,
    parent_dir: Path,
    config,
    pm,
    force: bool = False,
    offline_plots: bool = True,
    executed_argv: list[str] | None = None,
    execution_start: dict[str, Any] | None = None,
    workflow_started_at: float | None = None,
    invocation_state=None,
) -> None:
    """Create independent follow-up pools for a completed primary scheme."""
    started_at = workflow_started_at if workflow_started_at is not None else monotonic()
    parent_dir = Path(parent_dir).resolve()
    output = Path(output_dir).resolve()
    if not parent_dir.is_dir():
        raise ValueError("gap-completion parent must be an existing directory")
    if output == parent_dir or output.is_relative_to(parent_dir) or parent_dir.is_relative_to(output):
        raise ValueError("gap-completion output and parent paths must be disjoint")
    # The follow-up bundle owns its complete output tree.  Reusing a stale
    # directory would make old reports and BED rows indistinguishable from this
    # run, so even ``force`` cannot make an existing path eligible.
    if output.exists():
        raise ValueError(f"{output} already exists; choose a new output directory")
    if config.selection_algorithm != "legacy":
        raise ValueError("gap completion requires the legacy selection algorithm")
    if config.mapping != MappingType.FIRST:
        raise ValueError("gap completion requires --mapping first")
    if config.circular or config.input_bedfile is not None:
        raise ValueError("gap completion does not support circular or imported primers")
    if config.terminal_gap_policy != TerminalGapPolicy.LEGACY:
        raise ValueError("gap completion requires the legacy MSA generation policy")
    for source in msa:
        source = Path(source).resolve()
        if source == output or source.is_relative_to(output):
            raise ValueError("gap-completion output cannot contain an input MSA")
    parent_files = [
        "primer.bed",
        "amplicon.bed",
        "primertrim.amplicon.bed",
        "reference.fasta",
        "config.json",
        "panel-provenance.json",
    ]
    parent_baseline = {
        name: _file_record(parent_dir / name)
        for name in parent_files
        if (parent_dir / name).is_file()
    }
    input_baseline = [_file_record(Path(source).resolve()) for source in msa]
    output.mkdir(parents=True)
    (output / "work").mkdir(exist_ok=True)
    from primalscheme3.core.logger import setup_rich_logger

    logger = setup_rich_logger(str(output / "work" / "file.log"))
    if invocation_state is not None:
        invocation_state.output_owned = True
        invocation_state.logger = logger
    if pm is None:
        from primalscheme3.core.progress_tracker import ProgressManager

        pm = ProgressManager()

    input_records = []
    local_paths = []
    for index, source in enumerate(msa):
        local_name = f"input-{index:04d}-{source.name}"
        local_path = output / "work" / local_name
        shutil.copyfile(source, local_path)
        local_paths.append(local_path)
        input_records.append(
            {
                "sourcePath": str(source.resolve()),
                "storedPath": f"work/{local_name}",
                "sourceIndex": index,
            }
        )

    msa_dict: dict[int, MSA] = {}
    for index, local_path in enumerate(local_paths):
        msa_dict[index] = MSA(
            name=msa[index].stem,
            path=local_path,
            msa_index=index,
            logger=logger,
            progress_manager=pm,
            config=config,
        )
    parent_pairs, parent_verification = _validate_parent(parent_dir, msa_dict, config, logger)
    parent_snapshot = output / "parent"
    parent_snapshot.mkdir(exist_ok=True)
    for name in parent_files:
        source = parent_dir / name
        if source.is_file():
            shutil.copy2(source, parent_snapshot / name)
    parent_manifest = {
        name: _file_record(parent_dir / name)
        for name in parent_files
        if (parent_dir / name).is_file()
    }
    if parent_manifest != parent_baseline:
        raise ValueError("parent scheme changed while it was being copied")
    for name, baseline in parent_baseline.items():
        copied = _file_record(parent_snapshot / name, relative_to=output)
        if copied["sha256"] != baseline["sha256"] or copied["size"] != baseline["size"]:
            raise ValueError(f"copied parent file does not match source: {name}")

    candidates = []
    for msa_obj in msa_dict.values():
        msa_obj.digest_rs(config, None)
        msa_obj.generate_primerpairs(
            amplicon_size_min=config.amplicon_size_min,
            amplicon_size_max=config.amplicon_size_max,
            dimerscore=config.dimer_score,
            amplicon_size_metric=config.amplicon_size_metric,
        )
        candidates.extend(msa_obj.primerpairs)
    match_db = MatchDB(
        output / "work" / "mismatch",
        [str(path) for path in local_paths] if config.use_matchdb else [],
        config,
    )
    followup = Multiplex(config, match_db, msa_dict)

    def admission(pair, pool):
        result = followup.check_primerpair_can_be_added(pair, pool)
        return None if result.name == "OK" else result.name.lower()

    def on_accept(pair, pool):
        followup.add_primer_pair_to_pool(pair, pool, int(pair.msa_index))

    selection = select_followup_candidates(
        msa_dict,
        parent_pairs,
        candidates,
        pool_count=config.n_pools,
        admission=admission,
        on_accept=on_accept,
    )
    followup_bed = followup.to_bed()
    (output / "primer.bed").write_text(followup_bed)
    (output / "amplicon.bed").write_text(followup.to_amplicons(trim_primers=False))
    (output / "primertrim.amplicon.bed").write_text(
        followup.to_amplicons(trim_primers=True)
    )
    with dnaio.FastaWriter(output / "reference.fasta", line_length=60) as writer:
        for msa_obj in msa_dict.values():
            writer.write(
                dnaio.SequenceRecord(
                    name=msa_obj._chrom_name, sequence=generate_reference(msa_obj.array)
                )
            )
    fresh_validation = _validate_followup_output(
        output, msa_dict, local_paths, config, selection
    )
    for source, baseline in zip(msa, input_baseline, strict=True):
        _assert_unchanged(Path(source).resolve(), baseline)
    for name, baseline in parent_baseline.items():
        _assert_unchanged(parent_dir / name, baseline)
    coverage_report = {
        "schemaVersion": "primalscheme3.gap-completion-coverage/v1",
        "primary": selection.primary_coverage,
        "followup": selection.followup_coverage,
        "combined": selection.combined_coverage,
        "candidateUnionCeiling": selection.candidate_union_coverage,
        "remainingGaps": _remaining_gaps(selection.combined_coverage),
        "candidateUnionRemainingGaps": _remaining_gaps(selection.candidate_union_coverage),
        "freshValidation": fresh_validation,
        "sourceDrift": {"parentChangedDuringRun": False, "inputsChangedDuringRun": False},
    }
    (output / "gap-completion-coverage.json").write_text(
        json.dumps(coverage_report, sort_keys=True, separators=(",", ":")) + "\n"
    )
    accepted_by_id = {item.candidate_id: item for item in selection.accepted}
    manifest = []
    manifest_ids: set[str] = set()
    for pair in candidates:
        candidate_id = gap_candidate_id(pair, msa_dict[pair.msa_index])
        if candidate_id in manifest_ids:
            continue
        manifest_ids.add(candidate_id)
        intervals = _intervals(pair, msa_dict[pair.msa_index])
        status = selection.statuses.get(candidate_id, {"status": "not-evaluated"})
        record = {
            "candidateId": candidate_id,
            "msaIndex": int(pair.msa_index),
            "targetOccurrence": msa_dict[pair.msa_index]._chrom_name,
            "trimmedInterval": list(intervals[0]) if intervals else None,
            "fullInterval": list(intervals[1]) if intervals else None,
            "forwardOligos": sorted(str(seq).upper() for seq in pair.fprimer.seqs()),
            "reverseOligos": sorted(str(seq).upper() for seq in pair.rprimer.seqs()),
            **status,
        }
        accepted = accepted_by_id.get(candidate_id)
        if accepted is not None:
            record.update(
                {
                    "pool": accepted.pool,
                    "poolId": f"followup-pool-{accepted.pool + 1}",
                    "publishedAmpliconName": f"{pair.amplicon_prefix}_{pair.amplicon_number}",
                }
            )
        manifest.append(record)
    followup_id = "followup-scheme-" + hashlib.sha256(followup_bed.encode()).hexdigest()[:16]
    parent_id = "primary-scheme-" + _sha256(parent_dir / "primer.bed")[:16]
    report = {
        "schemaVersion": "primalscheme3.gap-completion/v1",
        "algorithm": "bounded-legacy-gap-completion/v1",
        "parent": {
            "schemeId": parent_id,
            "path": str(parent_dir),
            "files": parent_manifest,
            "verification": parent_verification,
            "preRunFiles": parent_baseline,
            "unchangedAfterRun": True,
        },
        "inputSnapshot": input_baseline,
        "followup": {
            "schemeId": followup_id,
            "poolCount": config.n_pools,
            "poolIds": [f"followup-pool-{pool + 1}" for pool in range(config.n_pools)],
        },
        "candidateManifest": manifest,
        "coverage": coverage_report,
    }
    (output / "gap-completion.json").write_text(
        json.dumps(report, sort_keys=True, separators=(",", ":")) + "\n"
    )
    config_dict = config.to_dict()
    config_dict.update(
        {
            "gap_completion_parent": str(parent_dir),
            "gap_completion_pool_count": int(config.n_pools),
            "gap_completion_report": "gap-completion.json",
            "gap_completion_coverage": "gap-completion-coverage.json",
        }
    )
    (output / "config.json").write_text(json.dumps(config_dict, sort_keys=True))
    from primalscheme3.panel.coverage_provenance import finalize_provenance

    finalize_provenance(
        output_dir=output,
        argv=executed_argv or list(sys.argv),
        resolved_options=config_dict,
        inputs=input_records,
        started_at=started_at,
        ended_at=monotonic(),
        status="success",
        exit_status=0,
        stderr="",
        scientific={
            "selectionAlgorithm": "legacy-gap-completion",
            "algorithm": "bounded-legacy-gap-completion/v1",
            "parent": report["parent"],
            "coverage": coverage_report,
        },
        logger=logger,
        execution_start=execution_start,
    )
    if invocation_state is not None:
        invocation_state.provenance_finalized = True
