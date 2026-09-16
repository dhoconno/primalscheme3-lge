"""Display-only source labels for immutable allele and row identifiers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import dnaio

from primalscheme3.core.msa import parse_chrom_name, parse_msa

SCHEMA = "primalscheme3.allele-label-map/v1"


def _contained(root: Path, relative: str) -> Path:
    root = root.resolve()
    path = (root / relative).resolve()
    if path == root or root not in path.parents:
        raise ValueError("stored input path escapes panel bundle")
    return path


def _descriptor(path: Path, relative: str) -> dict[str, Any]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"stored_path": relative, "sha256": digest.hexdigest(), "size": size}


def _records(path: Path) -> list[dict[str, str | None]]:
    records = []
    with dnaio.open(path) as stream:
        for record in stream:
            records.append(
                {
                    "fasta_description": record.name,
                    "fasta_record_id": record.id,
                    "fasta_comment": record.comment,
                    "normalized_record_id": parse_chrom_name(record.id),
                }
            )
    return records


def _row_digest(cells) -> str:
    encoded = json.dumps(list(cells), separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_allele_label_map(catalog, bundle_root, input_records):
    """Bind saved FASTA headers to catalog rows without changing scientific IDs."""

    root = Path(bundle_root).resolve()
    indexes = [
        item.get("sourceIndex", index) for index, item in enumerate(input_records)
    ]
    if (
        any(type(index) is not int or index < 0 for index in indexes)
        or len(set(indexes)) != len(indexes)
        or set(indexes) != {target.source_msa_index for target in catalog.targets}
    ):
        raise ValueError("input source occurrences differ from catalog targets")
    target_by_source = {target.source_msa_index: target for target in catalog.targets}
    if len(target_by_source) != len(catalog.targets):
        raise ValueError("catalog has ambiguous source occurrences")

    inputs = []
    targets = []
    for position, input_record in enumerate(input_records):
        source_index = input_record.get("sourceIndex", position)
        relative = input_record.get("storedPath")
        if not isinstance(relative, str) or not relative:
            raise ValueError("input stored path is unavailable")
        path = _contained(root, relative)
        if not path.is_file():
            raise ValueError("stored input is unavailable: " + relative)
        inputs.append({"source_msa_index": source_index, **_descriptor(path, relative)})

        # parse_msa is the authoritative native parser. It retains source row
        # order and rejects duplicate normalized record IDs exactly as design did.
        array, parsed_index = parse_msa(path)
        labels = _records(path)
        target = target_by_source[source_index]
        cells = tuple(tuple(str(cell) for cell in row) for row in array)
        if cells != target.rows or len(labels) != len(target.row_ids):
            raise ValueError("saved input row order/count/content differs from catalog")
        normalized_ids = tuple(item["normalized_record_id"] for item in labels)
        if normalized_ids != tuple(parsed_index):
            raise ValueError(
                "saved input labels differ from authoritative parser order"
            )

        rows = []
        by_row_id = {}
        for row_index, (row_id, row_cells, label) in enumerate(
            zip(target.row_ids, target.rows, labels, strict=True)
        ):
            row = {
                "row_index": row_index,
                "row_id": row_id,
                "row_content_sha256": _row_digest(row_cells),
                **label,
            }
            rows.append(row)
            by_row_id[row_id] = row

        classes = []
        for observation in catalog.observations:
            if observation.target_id != target.id:
                continue
            row_ids = [
                row_id for row_id in target.row_ids if row_id in observation.row_ids
            ]
            if len(row_ids) != observation.multiplicity or set(row_ids) != set(
                observation.row_ids
            ):
                raise ValueError("observed class row aliases differ from target rows")
            classes.append(
                {
                    "allele_id": observation.id,
                    "multiplicity": observation.multiplicity,
                    "row_ids": row_ids,
                    "fasta_record_ids": [
                        by_row_id[row_id]["fasta_record_id"] for row_id in row_ids
                    ],
                    "fasta_descriptions": [
                        by_row_id[row_id]["fasta_description"] for row_id in row_ids
                    ],
                }
            )
        targets.append(
            {
                "target_id": target.id,
                "source_msa_index": source_index,
                "occurrence": target.occurrence,
                "rows": rows,
                "classes": classes,
            }
        )

    return {
        "schemaVersion": SCHEMA,
        "scientificIdentityRole": "display-only-excluded",
        "scope": (
            "Source labels and immutable ID aliases for display only; excluded from "
            "catalog, cache, coverage, objective, configuration, and history identities."
        ),
        "inputs": inputs,
        "targets": sorted(targets, key=lambda item: item["source_msa_index"]),
    }
