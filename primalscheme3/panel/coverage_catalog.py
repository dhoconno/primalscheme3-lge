"""Build and persist deterministic coverage-panel candidate catalogues."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from primalscheme3.core.config import ALL_DNA_WITH_N, AMBIGUOUS_DNA_COMPLEMENT
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.panel.coverage_types import Candidate, Catalog, Target


CATALOG_SCHEMA = "primalscheme3.coverage-catalog/v1"
_CONCRETE_BASES = frozenset("ACGT")
_AnchoredMatch = tuple[str, int, int]
_RowSupport = tuple[tuple[_AnchoredMatch, ...], bool]
_SupportCache = dict[tuple[bool, int, tuple[str, ...]], tuple[_RowSupport, ...]]


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value)).hexdigest()


def _json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            str(key): _json_value(item)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "value"):
        return _json_value(value.value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _rows(msa: Any) -> tuple[tuple[str, ...], ...]:
    return tuple(
        tuple(str(cell).upper() for cell in row.tolist()) for row in msa.array
    )


def _target_records(msa_dict: dict[Any, Any]) -> tuple[tuple[Target, Any], ...]:
    entries = sorted(
        msa_dict.items(), key=lambda item: (int(item[1].msa_index), str(item[0]))
    )
    content_counts: Counter[str] = Counter()
    records: list[tuple[Target, Any]] = []
    for _, msa in entries:
        rows = _rows(msa)
        content_key = _digest({"rows": rows})
        occurrence = content_counts[content_key]
        content_counts[content_key] += 1
        target_id = f"target-{content_key}-{occurrence}"
        row_ids = tuple(
            "row-"
            + _digest(
                {
                    "target_id": target_id,
                    "row_index": row_index,
                    "row": row,
                }
            )
            for row_index, row in enumerate(rows)
        )
        mapping = tuple(
            None if value is None else int(value) for value in msa._mapping_array
        )
        reference_sequence = "".join(
            cell for cell in rows[0] if cell not in {"", "-"}
        )
        reference_length = len(reference_sequence)
        by_reference = {
            int(reference_index): alignment_index
            for alignment_index, reference_index in enumerate(mapping)
            if reference_index is not None
        }
        if reference_length:
            ref_to_alignment = tuple(
                by_reference[index] for index in range(reference_length)
            ) + (by_reference[reference_length - 1] + 1,)
        else:
            ref_to_alignment = (0,)
        records.append(
            (
                Target(
                    id=target_id,
                    source_msa_index=int(msa.msa_index),
                    occurrence=occurrence,
                    row_ids=row_ids,
                    rows=rows,
                    reference_sequence=reference_sequence,
                    reference_length=reference_length,
                    mapping=mapping,
                    ref_to_alignment=ref_to_alignment,
                ),
                msa,
            )
        )
    return tuple(sorted(records, key=lambda record: record[0].id))


def _canonical_oligos(kmer: Any) -> tuple[str, ...]:
    return tuple(sorted({str(sequence).upper() for sequence in kmer.seqs()}))


def _base_match(primer_base: str, row_base: str) -> str:
    if primer_base == "N" or row_base == "N" or row_base == "":
        return "unknown"
    primer_set = set(ALL_DNA_WITH_N.get(primer_base, ""))
    row_set = set(ALL_DNA_WITH_N.get(row_base, ""))
    if not primer_set or not row_set or primer_set.isdisjoint(row_set):
        return "mismatch"
    if primer_base not in _CONCRETE_BASES or row_base not in _CONCRETE_BASES:
        return "unknown"
    return "exact"


def _sequence_match(primer: str, observed: tuple[str, ...]) -> str:
    if len(primer) != len(observed) or "-" in observed:
        return "mismatch"
    statuses = tuple(
        _base_match(primer_base, row_base)
        for primer_base, row_base in zip(primer, observed, strict=True)
    )
    if "mismatch" in statuses:
        return "mismatch"
    if "unknown" in statuses:
        return "unknown"
    return "exact"


def _reverse_observed(observed: tuple[str, ...]) -> tuple[str, ...]:
    if all(base not in {"", "-"} for base in observed):
        return tuple(reverse_complement("".join(observed)))
    return tuple(
        "" if base == "" else AMBIGUOUS_DNA_COMPLEMENT[base]
        for base in reversed(observed)
    )


def _walk_anchor(
    row: tuple[str, ...], anchor_alignment: int, length: int, *, reverse: bool
) -> tuple[tuple[str, ...], int, int]:
    selected: list[tuple[int, str]] = []
    cursor = anchor_alignment if reverse else anchor_alignment - 1
    step = 1 if reverse else -1
    while 0 <= cursor < len(row) and len(selected) < length:
        cell = row[cursor]
        if cell != "-":
            selected.append((cursor, cell))
        cursor += step
    if len(selected) < length:
        selected.extend((cursor, "") for _ in range(length - len(selected)))
    if not reverse:
        selected.reverse()
    real_indexes = [index for index, _ in selected if 0 <= index < len(row)]
    start_alignment = min(real_indexes, default=max(0, anchor_alignment))
    end_alignment = max(real_indexes, default=anchor_alignment - 1) + 1
    return tuple(cell for _, cell in selected), start_alignment, end_alignment


def _anchored_matches(
    target: Target,
    row: tuple[str, ...],
    oligos: Iterable[str],
    anchor: int,
    *,
    reverse: bool,
) -> _RowSupport:
    """Return exact ``(oligo, row_start, row_end)`` hits and uncertainty."""

    matches: list[tuple[str, int, int]] = []
    uncertain = False
    for oligo in oligos:
        if anchor < 0 or anchor > target.reference_length:
            uncertain = True
            continue
        if reverse:
            if anchor == target.reference_length:
                anchor_alignment = target.ref_to_alignment[-1]
            else:
                anchor_alignment = target.ref_to_alignment[anchor]
        else:
            if anchor == 0:
                anchor_alignment = target.ref_to_alignment[0]
            else:
                anchor_alignment = target.ref_to_alignment[anchor - 1] + 1
        observed, start_alignment, end_alignment = _walk_anchor(
            row, anchor_alignment, len(oligo), reverse=reverse
        )
        if reverse:
            observed = _reverse_observed(observed)
        status = _sequence_match(oligo, observed)
        if status == "exact":
            row_start = sum(
                base not in {"", "-"} for base in row[:start_alignment]
            )
            row_end = sum(base not in {"", "-"} for base in row[:end_alignment])
            matches.append((oligo, row_start, row_end))
        elif status == "unknown":
            uncertain = True
    return tuple(matches), uncertain


def _cloud_support(
    target: Target,
    oligos: tuple[str, ...],
    anchor: int,
    *,
    reverse: bool,
    cache: _SupportCache,
) -> tuple[_RowSupport, ...]:
    key = (reverse, anchor, oligos)
    if key not in cache:
        cache[key] = tuple(
            _anchored_matches(
                target, row, oligos, anchor, reverse=reverse
            )
            for row in target.rows
        )
    return cache[key]


def _candidate(
    target: Target,
    pair: Any,
    support_cache: _SupportCache,
) -> Candidate:
    forward_oligos = _canonical_oligos(pair.fprimer)
    reverse_oligos = _canonical_oligos(pair.rprimer)
    forward_end = int(pair.fprimer.end)
    reverse_start = int(pair.rprimer.start)
    full_interval = (
        forward_end - max(map(len, forward_oligos)),
        reverse_start + max(map(len, reverse_oligos)),
    )
    interior_interval = (forward_end, reverse_start)
    identity = {
        "target_id": target.id,
        "full_interval": full_interval,
        "interior_interval": interior_interval,
        "forward_oligos": forward_oligos,
        "reverse_oligos": reverse_oligos,
    }
    joint_rows: list[str] = []
    unknown_rows: list[str] = []
    product_spans: list[tuple[str, int, int]] = []
    forward_support = _cloud_support(
        target,
        forward_oligos,
        forward_end,
        reverse=False,
        cache=support_cache,
    )
    reverse_support = _cloud_support(
        target,
        reverse_oligos,
        reverse_start,
        reverse=True,
        cache=support_cache,
    )
    for row_id, forward_row, reverse_row in zip(
        target.row_ids, forward_support, reverse_support, strict=True
    ):
        forward_matches, forward_unknown = forward_row
        reverse_matches, reverse_unknown = reverse_row
        if forward_matches and reverse_matches:
            joint_rows.append(row_id)
            row_spans = {
                (forward_match[1], reverse_match[2])
                for forward_match in forward_matches
                for reverse_match in reverse_matches
            }
            product_spans.extend(
                (row_id, start, end) for start, end in sorted(row_spans)
            )
        elif forward_unknown or reverse_unknown:
            unknown_rows.append(row_id)
    return Candidate(
        id="candidate-" + _digest(identity),
        target_id=target.id,
        full_interval=full_interval,
        interior_interval=interior_interval,
        forward_oligos=forward_oligos,
        reverse_oligos=reverse_oligos,
        joint_rows=tuple(joint_rows),
        unknown_rows=tuple(unknown_rows),
        row_product_spans=tuple(product_spans),
    )


def _target_json(target: Target, *, include_source: bool) -> dict[str, Any]:
    record = {
        "id": target.id,
        "occurrence": target.occurrence,
        "row_ids": target.row_ids,
        "rows": target.rows,
        "reference_sequence": target.reference_sequence,
        "reference_length": target.reference_length,
        "mapping": target.mapping,
        "ref_to_alignment": target.ref_to_alignment,
    }
    if include_source:
        record["source_msa_index"] = target.source_msa_index
    return record


def _candidate_json(candidate: Candidate, oligo_ids: dict[str, str]) -> dict[str, Any]:
    return {
        "id": candidate.id,
        "target_id": candidate.target_id,
        "full_interval": candidate.full_interval,
        "interior_interval": candidate.interior_interval,
        "forward_oligo_ids": tuple(
            oligo_ids[sequence] for sequence in candidate.forward_oligos
        ),
        "reverse_oligo_ids": tuple(
            oligo_ids[sequence] for sequence in candidate.reverse_oligos
        ),
        "joint_rows": candidate.joint_rows,
        "unknown_rows": candidate.unknown_rows,
        "row_product_spans": candidate.row_product_spans,
    }


def _semantic_payload(
    targets: tuple[Target, ...],
    candidates: tuple[Candidate, ...],
) -> dict[str, Any]:
    oligos = sorted(
        {
            sequence
            for candidate in candidates
            for sequence in candidate.forward_oligos + candidate.reverse_oligos
        }
    )
    oligo_ids = {sequence: "oligo-" + _digest(sequence) for sequence in oligos}
    return {
        "schema": CATALOG_SCHEMA,
        "targets": tuple(
            _target_json(target, include_source=False) for target in targets
        ),
        "oligos": tuple(
            {"id": oligo_ids[sequence], "sequence": sequence} for sequence in oligos
        ),
        "candidates": tuple(
            _candidate_json(candidate, oligo_ids) for candidate in candidates
        ),
    }


def build_catalog(msa_dict: dict[Any, Any], config: Any) -> Catalog:
    """Canonicalize native MSA primer pairs without mutating native objects."""

    target_records = _target_records(msa_dict)
    targets = tuple(target for target, _ in target_records)
    native_pairs: dict[str, object] = {}
    candidates: dict[str, Candidate] = {}
    for target, msa in target_records:
        support_cache: _SupportCache = {}
        for pair in msa.primerpairs:
            candidate = _candidate(target, pair, support_cache)
            candidates.setdefault(candidate.id, candidate)
            native_pairs.setdefault(candidate.id, pair)
    ordered_candidates = tuple(sorted(candidates.values(), key=lambda item: item.id))
    source_mapping = tuple(
        (target.source_msa_index, target.id) for target in targets
    )
    resolved_config = _json_value(config.to_dict())
    resolved_config_json = _canonical_json(resolved_config).decode("utf-8")
    semantic = _semantic_payload(targets, ordered_candidates)
    return Catalog(
        targets=targets,
        candidates=ordered_candidates,
        semantic_digest=_digest(semantic),
        resolved_config_json=resolved_config_json,
        source_mapping=source_mapping,
        _native_pair_items=tuple(
            (candidate.id, native_pairs[candidate.id])
            for candidate in ordered_candidates
        ),
    )


def write_catalog(catalog: Catalog, path: str | Path) -> dict[str, str | int]:
    """Write deterministic gzip JSON and return provenance for the actual bytes."""

    output_path = Path(path)
    payload = _semantic_payload(catalog.targets, catalog.candidates)
    payload["targets"] = tuple(
        _target_json(target, include_source=True) for target in catalog.targets
    )
    payload["metadata"] = {
        "resolved_config": json.loads(catalog.resolved_config_json),
        "source_mapping": catalog.source_mapping,
    }
    payload["semantic_digest"] = catalog.semantic_digest
    buffer = io.BytesIO()
    with gzip.GzipFile(
        filename="", mode="wb", fileobj=buffer, mtime=0
    ) as compressed:
        compressed.write(_canonical_json(payload))
    data = buffer.getvalue()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(data)
    return {
        "path": str(output_path),
        "semantic_digest": catalog.semantic_digest,
        "file_sha256": hashlib.sha256(data).hexdigest(),
        "file_size": len(data),
    }
