#!/usr/bin/env python3
"""Score a historical native primer BED with the observed-allele metric.

This command deliberately evaluates binding and primer-trimmed coverage only.
It does not apply current discovery eligibility, chemistry, dimer, specificity,
or panel-validity filters to historical primer records.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shlex
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from statistics import mean
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from primalscheme3.panel.allele_coverage import (  # noqa: E402
    canonical_observations,
    coverage_utility,
)
from primalscheme3.panel.coverage_types import (  # noqa: E402
    ClassCoverage,
    CoverageSummary,
    Target,
    TargetCoverage,
)

TOOL = "score-legacy-allele-coverage"
TOOL_VERSION = "1"
SCHEMA = "primalscheme3.legacy-allele-coverage/v1"
PROVENANCE_SCHEMA = "primalscheme3.scientific-command-provenance/v1"
METRIC = "observed-allele-primer-trimmed/v1"
GOAL = 0.95
SCOPE = (
    "BINDING-MODEL COVERAGE ONLY: historical primers receive exact same-row "
    "primer-trimmed coverage credit without current eligibility, chemistry, "
    "dimer, specificity, or panel-validity filtering. Invalid historical "
    "chemistry can therefore earn coverage in this report."
)
_BASES = frozenset("ACGTRYSWKMBDHVN")
_CONCRETE = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_PRIMER_NAME = re.compile(r"^(?P<family>.+)_(?P<side>LEFT|RIGHT)_(?P<variant>[1-9][0-9]*)$")


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _descriptor(path: Path, *, role: str | None = None, relative_to: Path | None = None) -> dict[str, Any]:
    resolved = path.resolve()
    record: dict[str, Any] = {
        "path": resolved.relative_to(relative_to.resolve()).as_posix() if relative_to else str(resolved),
        "sha256": _sha256(resolved),
        "size": resolved.stat().st_size,
    }
    if role is not None:
        record["role"] = role
    return record


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _native_name(header: str) -> str:
    return header.split()[0].replace("/", "_").replace("-", "_")


def _parse_fasta(path: Path) -> tuple[tuple[str, ...], tuple[tuple[str, ...], ...]]:
    headers: list[str] = []
    sequences: list[str] = []
    current: list[str] | None = None
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip().split()[0] if line[1:].strip() else ""
            if not header:
                raise ValueError(f"{path}: empty FASTA header at line {line_number}")
            headers.append(header)
            sequences.append("")
            current = []
            continue
        if current is None:
            raise ValueError(f"{path}: sequence before first FASTA header")
        sequences[-1] += line.upper()
    if not headers:
        raise ValueError(f"{path}: empty FASTA")
    if len(set(map(_native_name, headers))) != len(headers):
        raise ValueError(f"{path}: duplicate native-normalized FASTA header")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ValueError(f"{path}: MSA rows must have one positive aligned length")
    for row_index, sequence in enumerate(sequences):
        invalid = set(sequence) - _BASES - {"-"}
        if invalid:
            raise ValueError(f"{path}: invalid FASTA bases in row {row_index}: {sorted(invalid)}")
    cells = []
    for sequence in sequences:
        first = next((i for i, base in enumerate(sequence) if base != "-"), len(sequence))
        last = next((i for i in range(len(sequence) - 1, -1, -1) if sequence[i] != "-"), -1)
        cells.append(
            tuple("" if base == "-" and (i < first or i > last) else base for i, base in enumerate(sequence))
        )
    return tuple(headers), tuple(cells)


def _target(identity: str, index: int, occurrence: int, headers: tuple[str, ...], rows: tuple[tuple[str, ...], ...]) -> Target:
    reference_columns = tuple(i for i, cell in enumerate(rows[0]) if cell not in ("", "-"))
    mapping_by_column = {column: position for position, column in enumerate(reference_columns)}
    mapping = tuple(mapping_by_column.get(i) for i in range(len(rows[0])))
    reference = "".join(rows[0][i] for i in reference_columns)
    return Target(
        identity,
        index,
        occurrence,
        headers,
        rows,
        reference,
        len(reference),
        mapping,
        reference_columns + ((reference_columns[-1] + 1) if reference_columns else 0,),
    )


def _load_msas(paths: list[Path]) -> list[dict[str, Any]]:
    content_counts: Counter[str] = Counter()
    records = []
    for index, path in enumerate(paths):
        headers, rows = _parse_fasta(path)
        digest = hashlib.sha256(
            json.dumps(rows, separators=(",", ":")).encode()
        ).hexdigest()
        occurrence = content_counts[digest]
        content_counts[digest] += 1
        identity = f"target-{digest}-{occurrence}"
        records.append(
            {
                "path": path.resolve(),
                "headers": headers,
                "target": _target(identity, index, occurrence, headers, rows),
            }
        )
    return records


def _load_row_map(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError(f"{path}: invalid row-map JSON: {error}") from error
    if not isinstance(value, dict) or not isinstance(value.get("inputID"), str) or not isinstance(value.get("rows"), list):
        raise ValueError(f"{path}: invalid row-map structure")
    return value


def _headers_match_map(headers: tuple[str, ...], rows: list[dict[str, Any]]) -> bool:
    if len(headers) != len(rows):
        return False
    for index, (header, row) in enumerate(zip(headers, rows, strict=True)):
        if not isinstance(row, dict) or row.get("rowIndex") != index:
            return False
        normalized = row.get("normalizedHeader")
        original = row.get("originalHeader")
        if not isinstance(normalized, str) or not isinstance(original, str):
            return False
        if _native_name(header) not in {_native_name(normalized), _native_name(original)}:
            return False
    return True


def _pair_row_maps(msas: list[dict[str, Any]], paths: list[Path]) -> dict[int, dict[str, Any]]:
    paired: dict[int, dict[str, Any]] = {}
    for path in paths:
        value = _load_row_map(path)
        input_id = value["inputID"]
        by_identity = [
            index
            for index, msa in enumerate(msas)
            if input_id == msa["path"].stem or input_id in msa["path"].parts
        ]
        by_headers = [
            index for index, msa in enumerate(msas) if _headers_match_map(msa["headers"], value["rows"])
        ]
        candidates = sorted(set(by_identity or by_headers))
        if len(candidates) != 1:
            raise ValueError(f"{path}: row map must pair to exactly one --msa by sourceOccurrenceID, stem, or headers")
        index = candidates[0]
        if index in paired:
            raise ValueError(f"{path}: multiple row maps pair to one --msa")
        msa = msas[index]
        if len(value["rows"]) != len(msa["headers"]):
            raise ValueError(f"{path}: row count does not match paired MSA")
        if not _headers_match_map(msa["headers"], value["rows"]):
            raise ValueError(f"{path}: row headers or indexes do not match paired MSA")
        paired[index] = {"path": path.resolve(), "value": value}
    return paired


def _load_labels(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError(f"{path}: invalid labels JSON: {error}") from error
    if value.get("schemaVersion") != "allele-coverage-source-label-map/v1" or not isinstance(value.get("inputs"), list):
        raise ValueError(f"{path}: unsupported labels schema")
    result = {}
    for item in value["inputs"]:
        identity, label = item.get("sourceOccurrenceID"), item.get("label")
        if not isinstance(identity, str) or not identity or not isinstance(label, str) or not label:
            raise ValueError(f"{path}: invalid labels entry")
        if identity in result:
            raise ValueError(f"{path}: duplicate sourceOccurrenceID {identity}")
        result[identity] = label
    return result


def _parse_bed(path: Path) -> list[dict[str, Any]]:
    records = []
    names = set()
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 8:
            raise ValueError(f"{path}:{line_number}: native primer BED requires exactly 8 columns")
        chrom, start_text, end_text, name, pool_text, strand, sequence, attributes = columns
        match = _PRIMER_NAME.fullmatch(name)
        if match is None:
            raise ValueError(f"{path}:{line_number}: primer name must end in _LEFT_N or _RIGHT_N")
        try:
            start, end, pool = int(start_text), int(end_text), int(pool_text)
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: BED start, end, and pool must be integers") from error
        sequence = sequence.upper()
        side = match.group("side")
        expected_strand = "+" if side == "LEFT" else "-"
        if not chrom or start < 0 or end <= start or pool < 1:
            raise ValueError(f"{path}:{line_number}: invalid BED coordinates or pool")
        if strand != expected_strand:
            raise ValueError(f"{path}:{line_number}: {side} primer requires strand {expected_strand}")
        if len(sequence) != end - start or not sequence or set(sequence) - _CONCRETE:
            raise ValueError(f"{path}:{line_number}: sequence must be concrete ACGT and match BED span")
        if name in names:
            raise ValueError(f"{path}:{line_number}: duplicate primer name {name}")
        names.add(name)
        records.append(
            {
                "chromosome": chrom,
                "start": start,
                "end": end,
                "name": name,
                "pool": pool,
                "strand": strand,
                "sequence": sequence,
                "attributes": attributes,
                "familyID": match.group("family"),
                "side": side,
                "variant": int(match.group("variant")),
            }
        )
    return records


def _resolve_chromosomes(
    msas: list[dict[str, Any]],
    row_maps: dict[int, dict[str, Any]],
    bed_records: list[dict[str, Any]],
    labels: dict[str, str],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    bed_chromosomes = {record["chromosome"] for record in bed_records}
    mappings = []
    chromosome_to_index: dict[str, int] = {}
    for index, msa in enumerate(msas):
        direct = _native_name(msa["headers"][0])
        row_map = row_maps.get(index)
        source_id = None
        mapped = None
        if row_map is not None:
            source_id = row_map["value"]["inputID"]
            mapped = _native_name(row_map["value"]["rows"][0]["originalHeader"])
        else:
            source_ids = sorted(
                identity
                for identity in labels
                if identity == msa["path"].stem or identity in msa["path"].parts
            )
            if len(source_ids) > 1:
                raise ValueError(
                    f"{msa['path']}: multiple reporting labels match source path"
                )
            source_id = source_ids[0] if source_ids else None
        matches = {candidate for candidate in (direct, mapped) if candidate in bed_chromosomes}
        if len(matches) > 1:
            raise ValueError(
                f"{msa['path']}: direct and explicit row-map chromosomes both occur in BED; mapping is ambiguous"
            )
        if matches:
            chromosome = next(iter(matches))
            evidence = "direct-first-header" if chromosome == direct else "explicit-row-map"
        elif mapped is not None:
            chromosome = mapped
            evidence = "explicit-row-map"
        else:
            chromosome = direct
            evidence = "direct-first-header"
        if chromosome in chromosome_to_index:
            raise ValueError(f"chromosome {chromosome} maps to multiple supplied MSAs")
        chromosome_to_index[chromosome] = index
        mappings.append(
            {
                "msaPath": str(msa["path"]),
                "nativeInputIndex": index,
                "firstRowHeader": msa["headers"][0],
                "directNativeChromosome": direct,
                "chromosome": chromosome,
                "mappingEvidence": evidence,
                "sourceOccurrenceID": source_id,
                "rowMapPath": str(row_map["path"]) if row_map else None,
                "rowMapFirstOriginalHeader": row_map["value"]["rows"][0]["originalHeader"] if row_map else None,
                "label": labels.get(source_id) if source_id else None,
            }
        )
    missing = sorted(bed_chromosomes - set(chromosome_to_index))
    if missing:
        raise ValueError(
            "BED chromosomes do not map to supplied first rows: "
            + ", ".join(missing)
            + "; provide explicit --row-map files for normalized snapshots"
        )
    return mappings, chromosome_to_index


def _families(
    records: list[dict[str, Any]], msas: list[dict[str, Any]], chromosome_to_index: dict[str, int]
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        target = msas[chromosome_to_index[record["chromosome"]]]["target"]
        if record["end"] > target.reference_length:
            raise ValueError(f"{record['name']}: BED coordinates exceed mapped first-row reference")
        record["targetID"] = target.id
        record["alignmentAnchor"] = (
            target.ref_to_alignment[record["end"] - 1] + 1
            if record["side"] == "LEFT"
            else target.ref_to_alignment[record["start"]]
        )
        grouped[(record["chromosome"], record["familyID"])].append(record)
    families = []
    for (chromosome, family_id), primers in sorted(grouped.items()):
        forward = sorted((p for p in primers if p["side"] == "LEFT"), key=lambda p: p["name"])
        reverse = sorted((p for p in primers if p["side"] == "RIGHT"), key=lambda p: p["name"])
        if not forward or not reverse:
            raise ValueError(f"{chromosome}:{family_id}: amplicon family requires LEFT and RIGHT primers")
        if len({p["pool"] for p in primers}) != 1:
            raise ValueError(f"{chromosome}:{family_id}: amplicon family spans multiple pools")
        if len({p["end"] for p in forward}) != 1 or len({p["start"] for p in reverse}) != 1:
            raise ValueError(f"{chromosome}:{family_id}: variants do not share native LEFT/RIGHT anchors")
        if forward[0]["end"] > reverse[0]["start"]:
            raise ValueError(f"{chromosome}:{family_id}: primer anchors are reversed or overlap")
        families.append(
            {
                "chromosome": chromosome,
                "familyID": family_id,
                "targetID": primers[0]["targetID"],
                "pool": primers[0]["pool"],
                "forward": forward,
                "reverse": reverse,
            }
        )
    return families


def _binding(primer: dict[str, Any], cells: tuple[str, ...], alignment_to_row: tuple[int | None, ...]) -> tuple[str, tuple[int, int] | None]:
    forward = primer["strand"] == "+"
    cursor = primer["alignmentAnchor"] - 1 if forward else primer["alignmentAnchor"]
    step = -1 if forward else 1
    selected = []
    while len(selected) < len(primer["sequence"]):
        if cursor < 0 or cursor >= len(cells) or cells[cursor] == "":
            return "unavailable", None
        if cells[cursor] != "-":
            selected.append(cursor)
        cursor += step
    selected.sort()
    genomic = tuple(cells[i] for i in selected)
    expected = (
        tuple(primer["sequence"])
        if forward
        else tuple(primer["sequence"].translate(_COMPLEMENT)[::-1])
    )
    statuses = [
        "unknown" if observed not in _CONCRETE else "exact" if observed == wanted else "mismatch"
        for observed, wanted in zip(genomic, expected, strict=True)
    ]
    if "mismatch" in statuses:
        return "mismatch", None
    if "unknown" in statuses:
        return "unknown", None
    positions = [alignment_to_row[i] for i in selected]
    if any(position is None for position in positions):
        raise AssertionError("selected concrete footprint lacks row coordinates")
    return "confirmed", (positions[0], positions[-1] + 1)


def _intervals(positions: set[int]) -> tuple[tuple[int, int], ...]:
    result: list[tuple[int, int]] = []
    for position in sorted(positions):
        if result and result[-1][1] == position:
            result[-1] = (result[-1][0], position + 1)
        else:
            result.append((position, position + 1))
    return tuple(result)


def _score(
    msas: list[dict[str, Any]],
    mappings: list[dict[str, Any]],
    bed_records: list[dict[str, Any]],
    families: list[dict[str, Any]],
) -> dict[str, Any]:
    by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for family in families:
        by_target[family["targetID"]].append(family)
    class_records: list[ClassCoverage] = []
    target_records: list[TargetCoverage] = []
    target_reports = []
    aggregate_support = Counter(confirmed=0, mismatch=0, unknown=0, unavailable=0)
    amp_reports = {
        (family["chromosome"], family["familyID"]): {
            "chromosome": family["chromosome"],
            "familyID": family["familyID"],
            "targetID": family["targetID"],
            "pool": family["pool"],
            "forwardVariantCount": len(family["forward"]),
            "reverseVariantCount": len(family["reverse"]),
            "possibleVariantCombinations": len(family["forward"]) * len(family["reverse"]),
            "supportedAlleleClasses": 0,
            "supportedProductCombinations": 0,
        }
        for family in families
    }
    for msa, mapping in zip(msas, mappings, strict=True):
        target = msa["target"]
        observations = canonical_observations(target)
        per_class = []
        for observation in observations:
            observed = set(observation.observed_positions)
            covered: set[int] = set()
            support_counts = Counter(confirmed=0, mismatch=0, unknown=0, unavailable=0)
            supported_families = []
            for family in by_target[target.id]:
                forward = [(primer, _binding(primer, observation.cells, observation.alignment_to_row)) for primer in family["forward"]]
                reverse = [(primer, _binding(primer, observation.cells, observation.alignment_to_row)) for primer in family["reverse"]]
                for _, (status, _) in forward + reverse:
                    support_counts[status] += 1
                    aggregate_support[status] += 1
                products = 0
                for _, (forward_status, forward_footprint) in forward:
                    for _, (reverse_status, reverse_footprint) in reverse:
                        if (
                            forward_status == reverse_status == "confirmed"
                            and forward_footprint[1] <= reverse_footprint[0]
                        ):
                            products += 1
                            covered.update(range(forward_footprint[1], reverse_footprint[0]))
                if products:
                    supported_families.append(family["familyID"])
                    amp = amp_reports[(family["chromosome"], family["familyID"])]
                    amp["supportedAlleleClasses"] += 1
                    amp["supportedProductCombinations"] += products
            recovered = observed & covered
            coverage = ClassCoverage(
                target.id,
                observation.id,
                len(observed),
                len(recovered),
                len(observation.row_to_alignment) - len(observed),
                observation.cells.count(""),
                observation.cells.count("-"),
                len(recovered) / len(observed) if observed else None,
                _intervals(recovered),
            )
            class_records.append(coverage)
            per_class.append(
                coverage.to_dict()
                | {
                    "rowIDs": list(observation.row_ids),
                    "multiplicity": observation.multiplicity,
                    "bindingSupportStatusCounts": dict(support_counts),
                    "supportingAmpliconFamilyIDs": sorted(supported_families),
                }
            )
        fractions = [record["fraction"] for record in per_class if record["fraction"] is not None]
        target_coverage = TargetCoverage(
            target.id,
            mean(fractions) if fractions else None,
            "assessable" if fractions else "unassessable-target",
            len(fractions),
            len(per_class) - len(fractions),
        )
        target_records.append(target_coverage)
        target_families = by_target[target.id]
        target_primers = [primer for family in target_families for primer in family["forward"] + family["reverse"]]
        pool_counts = []
        for pool in sorted({family["pool"] for family in target_families}):
            pool_counts.append(
                {
                    "pool": pool,
                    "ampliconFamilies": sum(family["pool"] == pool for family in target_families),
                    "primerRecords": sum(primer["pool"] == pool for primer in target_primers),
                }
            )
        target_reports.append(
            {
                "targetID": target.id,
                "nativeInputIndex": mapping["nativeInputIndex"],
                "sourceOccurrenceID": mapping["sourceOccurrenceID"],
                "label": mapping["label"],
                "msaPath": mapping["msaPath"],
                "chromosome": mapping["chromosome"],
                "inputRowCount": len(target.rows),
                "distinctObservedAlleleClassCount": len(observations),
                "ampliconFamilyCount": len(target_families),
                "primerRecordCount": len(target_primers),
                "pools": pool_counts,
                "coverage": target_coverage.to_dict(),
                "classes": per_class,
            }
        )
    target_fractions = tuple(record.fraction for record in target_records if record.fraction is not None)
    summary = CoverageSummary(
        tuple(class_records),
        tuple(target_records),
        coverage_utility(target_fractions, goal=GOAL),
        mean(target_fractions) if target_fractions else None,
        min(target_fractions) if target_fractions else None,
        sum(fraction > GOAL for fraction in target_fractions),
        sum(fraction >= GOAL for fraction in target_fractions),
        GOAL,
        "ok" if target_fractions else "no-assessable-targets",
    )
    pools = []
    for pool in sorted({record["pool"] for record in bed_records}):
        pools.append(
            {
                "pool": pool,
                "ampliconFamilies": sum(family["pool"] == pool for family in families),
                "primerRecords": sum(record["pool"] == pool for record in bed_records),
            }
        )
    return {
        "schemaVersion": SCHEMA,
        "metricID": METRIC,
        "valid": True,
        "invalidReasons": [],
        "scientificScope": SCOPE,
        "limitations": [
            "Historical primers are not filtered by current discovery eligibility or chemistry.",
            "Dimer, specificity, secondary-product, and panel-validity checks are outside this score.",
            "Only full exact strand-correct primer matches on the same observed row establish support.",
            "Distinct observed aligned allele classes have equal weight; duplicate multiplicities are reporting-only.",
        ],
        "coverage": summary.to_dict(),
        "counts": {
            "targets": len(msas),
            "inputRows": sum(len(msa["target"].rows) for msa in msas),
            "distinctObservedAlleleClasses": len(class_records),
            "primerRecords": len(bed_records),
            "ampliconFamilies": len(families),
        },
        "bindingSupportStatusCounts": dict(aggregate_support),
        "pools": pools,
        "targets": target_reports,
        "amplicons": list(amp_reports.values()),
    }


def _runtime_identity() -> dict[str, Any]:
    try:
        package_version = version("primalscheme3")
    except PackageNotFoundError:
        package_version = None
    container_markers = [str(path) for path in (Path("/.dockerenv"), Path("/run/.containerenv")) if path.exists()]
    return {
        "pythonVersion": platform.python_version(),
        "pythonImplementation": platform.python_implementation(),
        "pythonExecutable": sys.executable,
        "pythonPrefix": sys.prefix,
        "platform": platform.platform(),
        "kernel": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "primalscheme3Version": package_version,
        "conda": {
            "prefix": os.environ.get("CONDA_PREFIX"),
            "defaultEnvironment": os.environ.get("CONDA_DEFAULT_ENV"),
        },
        "container": {"detected": bool(container_markers), "markers": container_markers},
    }


def _source_identity() -> dict[str, Any]:
    script = Path(__file__).resolve()
    scientific_sources = (
        script,
        ROOT / "primalscheme3" / "panel" / "allele_coverage.py",
        ROOT / "primalscheme3" / "panel" / "coverage_types.py",
    )
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        commit = None
    return {
        "script": _descriptor(script),
        "scientificSources": [_descriptor(path) for path in scientific_sources],
        "gitCommit": commit,
    }


def _provenance(
    *,
    started_at: str,
    started: float,
    argv: list[str],
    inputs: list[dict[str, Any]],
    outputs: list[dict[str, Any]],
    status: str,
    exit_status: int,
    stderr: str,
    source_at_start: dict[str, Any],
    output_directory: Path,
    working_directory: Path,
    scientific_report_valid: bool,
) -> dict[str, Any]:
    current_inputs = []
    for record in inputs:
        path = Path(record["path"])
        current_inputs.append(_descriptor(path, role=record["role"]) if path.is_file() else None)
    source_at_end = _source_identity()
    return {
        "schemaVersion": PROVENANCE_SCHEMA,
        "tool": {"name": TOOL, "version": TOOL_VERSION},
        "workflow": {"name": "historical-native-bed-observed-allele-scoring", "version": "1"},
        "metricID": METRIC,
        "scientificScope": SCOPE,
        "command": {
            "argv": argv,
            "shell": shlex.join(argv),
            "workingDirectory": str(working_directory.resolve()),
        },
        "resolvedOptions": {
            "goal": GOAL,
            "mapping": "first",
            "terminalGapPolicy": "observed-only",
            "duplicateWeighting": "equal-distinct-observed-allele-classes",
            "chemistryEligibilityFiltering": False,
            "specificityFiltering": False,
            "requireExactSameRowForwardReverseSupport": True,
        },
        "runtime": _runtime_identity(),
        "source": source_at_start,
        "sourceChangedDuringRun": source_at_start["scientificSources"]
        != source_at_end["scientificSources"],
        "inputs": inputs,
        "inputsChangedDuringRun": any(before != after for before, after in zip(inputs, current_inputs, strict=True)),
        "outputs": outputs,
        "outputDirectory": str(output_directory.resolve()),
        "selfHashPolicy": "provenance.json is excluded to avoid a self-hash cycle",
        "startedAt": started_at,
        "finishedAt": _now(),
        "wallTimeSeconds": time.monotonic() - started,
        "status": status,
        "exitStatus": exit_status,
        "stderr": stderr,
        "scientificReportValid": scientific_report_valid,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--msa", action="append", required=True, type=Path, help="Aligned FASTA in native input order; repeat for each target")
    parser.add_argument("--row-map", action="append", default=[], type=Path, help="Explicit row map for a normalized MSA; repeat as needed")
    parser.add_argument("--labels", type=Path, help="Optional source-row-label-map.json used only for report labels")
    parser.add_argument("--bed", required=True, type=Path, help="Historical native 8-column primer.bed")
    parser.add_argument("--output", required=True, type=Path, help="Fresh output directory")
    return parser


def main() -> int:
    args = _parser().parse_args()
    output = args.output.resolve()
    if output.exists():
        print(f"output directory already exists: {output}", file=sys.stderr)
        return 2
    output.mkdir(parents=True)
    started = time.monotonic()
    started_at = _now()
    working_directory = Path.cwd()
    argv = [sys.executable, *sys.argv]
    source_at_start = _source_identity()
    input_specs = [(Path(path), "msa") for path in args.msa]
    input_specs.extend((Path(path), "row-map") for path in args.row_map)
    if args.labels is not None:
        input_specs.append((args.labels, "labels"))
    input_specs.append((args.bed, "bed"))
    inputs: list[dict[str, Any]] = []
    try:
        inputs = [_descriptor(path, role=role) for path, role in input_specs]
        msas = _load_msas([Path(path) for path in args.msa])
        row_maps = _pair_row_maps(msas, [Path(path) for path in args.row_map])
        labels = _load_labels(args.labels)
        bed_records = _parse_bed(args.bed)
        mappings, chromosome_to_index = _resolve_chromosomes(msas, row_maps, bed_records, labels)
        families = _families(bed_records, msas, chromosome_to_index)
        report = _score(msas, mappings, bed_records, families)
        coverage_path = output / "coverage.json"
        mapping_path = output / "resolved-chromosome-maps.json"
        _write_json(coverage_path, report)
        _write_json(
            mapping_path,
            {
                "schemaVersion": "primalscheme3.resolved-chromosome-maps/v1",
                "mapping": "first",
                "mappings": mappings,
            },
        )
        outputs = [
            _descriptor(coverage_path, relative_to=output),
            _descriptor(mapping_path, relative_to=output),
        ]
        provenance = _provenance(
            started_at=started_at,
            started=started,
            argv=argv,
            inputs=inputs,
            outputs=outputs,
            status="success",
            exit_status=0,
            stderr="",
            source_at_start=source_at_start,
            output_directory=output,
            working_directory=working_directory,
            scientific_report_valid=True,
        )
        invalid_reasons = []
        if provenance["inputsChangedDuringRun"]:
            invalid_reasons.append("scientific inputs changed during execution")
        if provenance["sourceChangedDuringRun"]:
            invalid_reasons.append("scientific source changed during execution")
        if invalid_reasons:
            report["valid"] = False
            report["invalidReasons"] = invalid_reasons
            _write_json(coverage_path, report)
            outputs = [
                _descriptor(coverage_path, relative_to=output),
                _descriptor(mapping_path, relative_to=output),
            ]
            stderr = "; ".join(invalid_reasons)
            failed = _provenance(
                started_at=started_at,
                started=started,
                argv=argv,
                inputs=inputs,
                outputs=outputs,
                status="error",
                exit_status=1,
                stderr=stderr,
                source_at_start=source_at_start,
                output_directory=output,
                working_directory=working_directory,
                scientific_report_valid=False,
            )
            for field in ("inputsChangedDuringRun", "sourceChangedDuringRun"):
                failed[field] = failed[field] or provenance[field]
            _write_json(output / "provenance.json", failed)
            print(stderr, file=sys.stderr)
            return 1
        _write_json(output / "provenance.json", provenance)
        return 0
    except Exception as error:
        stderr = str(error)
        print(stderr, file=sys.stderr)
        outputs = [
            _descriptor(path, relative_to=output)
            for path in sorted(output.rglob("*"))
            if path.is_file() and path.name != "provenance.json"
        ]
        provenance = _provenance(
            started_at=started_at,
            started=started,
            argv=argv,
            inputs=inputs,
            outputs=outputs,
            status="error",
            exit_status=1,
            stderr=stderr,
            source_at_start=source_at_start,
            output_directory=output,
            working_directory=working_directory,
            scientific_report_valid=False,
        )
        _write_json(output / "provenance.json", provenance)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
