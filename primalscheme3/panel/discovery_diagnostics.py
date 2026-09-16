"""Bounded, on-demand discovery diagnostics for completed allele panels.

The diagnostic is deliberately separate from panel creation.  It reloads the
immutable target/catalog bytes from a completed bundle, validates the saved
source/runtime binding, and then runs the normal discovery implementation over
one requested anchor (or anchor pair).  The resulting history is therefore a
reconstructed anchored trace; it is never presented as the historical reason
an optimizer did or did not select an entity.
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path
from time import monotonic
from typing import Any

from .allele_publication import _read, _targets, artifact_descriptor
from .coverage_discovery import build_variant_catalog, discovery_profiles
from .coverage_history import SQLiteCoverageHistory
from .coverage_provenance import (
    capture_execution_identity,
    finalize_provenance,
    runtime_identity,
    source_identity,
)
from .coverage_types import CandidateFamily, OligoSite, Target, VariantCatalog

SCHEMA = "primalscheme3.discovery-diagnostic/v1"
_OUTPUT_RECEIPT = "provenance.json"
_REPORT = "discovery-diagnostic.json"
_CATALOG = "discovery-catalog.json.gz"
_HISTORY = "history"


def _contained(root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("unsafe bundle-relative path: " + str(value))
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError("path escapes bundle: " + str(value))
    return resolved


def _json_write(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, sort_keys=True, indent=2, default=str) + "\n")


def _gzip_write(path: Path, value: Any) -> None:
    with gzip.open(path, "wt") as stream:
        json.dump(value, stream, sort_keys=True, separators=(",", ":"))


def _descriptor(path: Path, root: Path) -> dict[str, Any]:
    return artifact_descriptor(path, root)


def _input_records(bundle: Path, provenance: dict[str, Any]) -> list[dict[str, Any]]:
    records = provenance.get("inputs")
    if not isinstance(records, list) or not records:
        raise ValueError("completed source panel has no input provenance")
    resolved: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise ValueError("source input provenance record is not an object")
        stored = record.get("storedPath")
        if not isinstance(stored, str) or not stored:
            raise ValueError("source input provenance has no stored path")
        path = _contained(bundle, stored)
        if not path.is_file():
            raise ValueError("missing stored source input: " + stored)
        actual = _descriptor(path, bundle)
        expected_sha = record.get("sha256", record.get("sourceAtStartSha256"))
        expected_size = record.get("size", record.get("sourceAtStartSize"))
        if expected_sha != actual["sha256"] or expected_size != actual["size"]:
            raise ValueError("stored source input integrity mismatch: " + stored)
        source_index = record.get("sourceIndex", index)
        if type(source_index) is not int or source_index < 0:
            raise ValueError("source input index is invalid")
        resolved.append(
            {
                "path": path,
                "sourceIndex": source_index,
                "sourcePath": record.get("sourcePath", str(path)),
                "sha256": actual["sha256"],
                "size": actual["size"],
            }
        )
    if len({record["sourceIndex"] for record in resolved}) != len(resolved):
        raise ValueError("source input indexes are not unique")
    return sorted(resolved, key=lambda item: item["sourceIndex"])


def _source_preflight(
    bundle: Path,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    VariantCatalog,
    tuple[Target, ...],
    list[dict[str, Any]],
]:
    provenance_path = bundle / "panel-provenance.json"
    optimizer_path = bundle / "panel-optimizer.json"
    if not provenance_path.is_file() or not optimizer_path.is_file():
        raise ValueError("bundle is missing completed panel provenance or optimizer")
    provenance = _read(provenance_path)
    if (
        provenance.get("status") != "success"
        or provenance.get("exitStatus") != 0
        or provenance.get("sourceChangedDuringRun") is not False
        or provenance.get("runtimeChangedDuringRun") is not False
    ):
        raise ValueError("source panel provenance is not a completed stable run")
    for key in ("source", "sourceAtEnd", "runtime", "runtimeAtEnd"):
        if not isinstance(provenance.get(key), dict):
            raise ValueError("source panel provenance lacks " + key)
    if provenance["source"] != provenance["sourceAtEnd"]:
        raise ValueError("source panel source identity changed during its run")
    if provenance["runtime"] != provenance["runtimeAtEnd"]:
        raise ValueError("source panel runtime identity changed during its run")
    current_source = source_identity()
    current_runtime = runtime_identity()
    if provenance["source"] != current_source:
        raise ValueError("source runtime drift: editable source identity differs")
    if provenance["runtime"] != current_runtime:
        raise ValueError("source runtime drift: kernel/runtime identity differs")
    if not provenance["runtime"].get("nativeKernels"):
        raise ValueError("source runtime provenance has no native kernel manifest")
    inputs = _input_records(bundle, provenance)

    optimizer = _read(optimizer_path)
    if optimizer.get("schemaVersion") != "primalscheme3.panel-optimizer/v2":
        raise ValueError("unsupported panel optimizer schema")
    primary = optimizer.get("primaryTier")
    stages = optimizer.get("stages", [])
    stage = next((item for item in stages if item.get("stage_id") == primary), None)
    if stage is None:
        raise ValueError("completed panel has no primary stage")
    stage_path = stage.get("path")
    if not isinstance(stage_path, str) or not stage_path:
        raise ValueError("primary stage has no safe artifact path")
    stage_dir = _contained(bundle, stage_path)
    catalog_path = stage_dir / "catalog.json.gz"
    targets_path = stage_dir / "authoritative-targets.json.gz"
    if not catalog_path.is_file() or not targets_path.is_file():
        raise ValueError("primary stage lacks authoritative catalog/target bytes")
    output_descriptors = {
        item.get("path"): item
        for item in provenance.get("outputs", [])
        if isinstance(item, dict) and isinstance(item.get("path"), str)
    }
    for path in (
        bundle / "panel-optimizer.json",
        bundle / "config.json",
        catalog_path,
        targets_path,
    ):
        relative = path.relative_to(bundle).as_posix()
        expected = output_descriptors.get(relative)
        if expected is None or _descriptor(path, bundle) != expected:
            raise ValueError("source output integrity mismatch: " + relative)
    catalog = VariantCatalog.from_dict(_read(catalog_path))
    saved_targets = _targets(_read(targets_path))
    if tuple(sorted(saved_targets, key=lambda item: item.id)) != tuple(
        sorted(catalog.targets, key=lambda item: item.id)
    ):
        raise ValueError("stage target/catalog membership mismatch")
    # Reparse the stored inputs and require byte-identical scientific targets.
    from .allele_pipeline import _authoritative_targets

    parsed_targets = _authoritative_targets(
        bundle,
        [
            {
                "storedPath": str(item["path"].relative_to(bundle)),
                "sourceIndex": item["sourceIndex"],
            }
            for item in inputs
        ],
    )
    if tuple(sorted(parsed_targets, key=lambda item: item.id)) != tuple(
        sorted(saved_targets, key=lambda item: item.id)
    ):
        raise ValueError("stored/raw target correspondence mismatch")
    return provenance, optimizer, catalog, tuple(saved_targets), inputs


def _config_from_bundle(bundle: Path, optimizer: dict[str, Any]) -> Any:
    """Reconstruct only persisted Config/AlleleOptions fields.

    Arbitrary command strings and launcher state are intentionally excluded.
    The Config constructor validates the resulting resolved option set.
    """
    from primalscheme3.core.config import Config

    from .allele_options import AlleleOptions

    config_data = (
        _read(bundle / "config.json") if (bundle / "config.json").is_file() else {}
    )
    allowed = (
        set(dir(Config))
        | {field.name for field in __import__("dataclasses").fields(AlleleOptions)}
        | {"allele_options_json", "allele_requested_options_json"}
    )
    kwargs = {
        key: value
        for key, value in config_data.items()
        if key in allowed and not key.startswith("_")
    }
    kwargs["selection_algorithm"] = "allele-coverage"
    kwargs["terminal_gap_policy"] = "observed-only"
    kwargs["mapping"] = "first"
    kwargs["amplicon_size_metric"] = "reference-span"
    # The source options snapshot is the authority for discovery length/profile
    # selection; these fields are safe, explicit controls.
    options = optimizer.get("options", {})
    for key in ("candidate_profiles", "discovery_length_mode", "discovery_history"):
        if key in options:
            kwargs[key] = options[key]
    return Config(**kwargs)


def _profile_binding(catalog: VariantCatalog, config: Any) -> dict[str, Any]:
    source = json.loads(catalog.resolved_config_json)
    profiles = discovery_profiles(config)
    expected = source.get("profiles", {})
    actual = {
        name: {key: getattr(profile, key) for key in next(iter(expected.values()), {})}
        for name, profile in profiles.items()
        if name in expected
    }
    if expected and actual != expected:
        raise ValueError(
            "resolved profile/chemistry binding differs from source catalog"
        )
    if (
        source.get("amplicon_size_min") != config.amplicon_size_min
        or source.get("amplicon_size_max") != config.amplicon_size_max
    ):
        raise ValueError("resolved amplicon bounds differ from source catalog")
    return {"profiles": sorted(expected), "resolved": expected}


def _site_signature(site: OligoSite) -> tuple[Any, ...]:
    return (
        site.id,
        site.target_id,
        site.sequence,
        site.strand,
        site.alignment_anchor,
        site.reference_footprint,
        tuple(site.accepting_profile_ids),
        tuple(site.generated_profile_ids),
    )


def _family_signature(
    family: CandidateFamily, catalog: VariantCatalog
) -> dict[str, Any]:
    return {
        "id": family.id,
        "target_id": family.target_id,
        "anchor_pair": list(family.anchor_pair),
        "forward_sites": [
            _site_signature(catalog.site_by_id[sid]) for sid in family.forward_site_ids
        ],
        "reverse_sites": [
            _site_signature(catalog.site_by_id[sid]) for sid in family.reverse_site_ids
        ],
        "profiles": [list(item) for item in family.discovery_profile_combinations],
    }


def _site_report(site: OligoSite) -> dict[str, Any]:
    return {
        "id": site.id,
        "target_id": site.target_id,
        "sequence": site.sequence,
        "strand": site.strand,
        "alignment_anchor": site.alignment_anchor,
        "reference_footprint": list(site.reference_footprint)
        if site.reference_footprint is not None
        else None,
        "accepting_profile_ids": list(site.accepting_profile_ids),
        "generated_profile_ids": list(site.generated_profile_ids),
    }


def _entity_binding(
    source: VariantCatalog, family_id: str | None, site_id: str | None
) -> tuple[dict[str, Any], tuple[list[int], list[int]], tuple[str, ...]]:
    if (family_id is None) == (site_id is None):
        raise ValueError("exactly one of family_id or site_id is required")
    if family_id is not None:
        family = source.family_by_id.get(family_id)
        if family is None:
            raise ValueError("unknown family_id: " + family_id)
        return (
            {
                "kind": "family",
                "id": family.id,
                "family": _family_signature(family, source),
            },
            ([family.anchor_pair[0]], [family.anchor_pair[1]]),
            (family.target_id,),
        )
    site = source.site_by_id.get(site_id)
    if site is None:
        raise ValueError("unknown site_id: " + site_id)
    indexes = (
        ([site.alignment_anchor], [])
        if site.strand == "+"
        else ([], [site.alignment_anchor])
    )
    return (
        {"kind": "site", "id": site.id, "site": _site_report(site)},
        (indexes[0], indexes[1]),
        (site.target_id,),
    )


def _invoke_discovery(targets, config, *, profiles, indexes, history):
    kwargs = {
        "profiles": profiles,
        "indexes": indexes,
        "history": history,
        "length_mode": getattr(config, "discovery_length_mode", None),
        "history_detail": "full",
    }
    return build_variant_catalog(targets, config, **kwargs)


def _write_receipt_path(output: Path) -> None:
    generated = output / "panel-provenance.json"
    receipt = output / _OUTPUT_RECEIPT
    if generated.is_file():
        data = json.loads(generated.read_text())
        data["tool"] = "primalscheme3 panel-discovery-diagnose"
        data["workflow"] = "panel-discovery-diagnose"
        data["selfHashPolicy"] = (
            "provenance.json is excluded to avoid a self-hash cycle"
        )
        generated.unlink()
        receipt.write_text(
            json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n"
        )


def _receipt_is_stable_success(record: dict[str, Any]) -> bool:
    return bool(
        record.get("status") == "success"
        and record.get("exitStatus") == 0
        and record.get("sourceChangedDuringRun") is False
        and record.get("runtimeChangedDuringRun") is False
        and all(
            item.get("sourceChangedDuringRun") is False
            for item in record.get("inputs", [])
        )
    )


def _failure_receipt(
    output: Path,
    *,
    argv: list[str],
    error: BaseException,
    inputs: list[dict[str, Any]],
    started: float,
    execution_start: dict[str, Any] | None,
) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    normalized_inputs = []
    for index, item in enumerate(inputs):
        source = (
            Path(item["path"]).resolve()
            if "path" in item
            else Path(item.get("sourcePath", "")).resolve()
        )
        stored = f"source/{index:04d}-{source.name}"
        if source.is_file():
            destination = output / stored
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            normalized_inputs.append(
                {
                    "sourcePath": str(source),
                    "storedPath": stored,
                    "sourceIndex": item.get("sourceIndex", index),
                }
            )
        else:
            normalized_inputs.append(
                {"sourcePath": str(source), "storedPath": stored, "sourceIndex": index}
            )
    if execution_start is None:
        try:
            execution_start = capture_execution_identity(
                [Path(item["sourcePath"]) for item in normalized_inputs]
            )
        except Exception:
            execution_start = None
    try:
        record = finalize_provenance(
            output_dir=output,
            argv=argv,
            resolved_options={},
            inputs=normalized_inputs,
            started_at=started,
            ended_at=monotonic(),
            status="failure",
            exit_status=1,
            stderr=str(error),
            scientific={"schemaVersion": SCHEMA, "valid": False},
            execution_start=execution_start,
        )
    except Exception as provenance_error:
        record = {
            "schemaVersion": SCHEMA,
            "tool": "primalscheme3 panel-discovery-diagnose",
            "command": {"argv": list(argv)},
            "inputs": [],
            "outputs": [],
            "status": "failure",
            "exitStatus": 1,
            "stderr": f"{error}; provenance: {provenance_error}",
        }
        (output / "panel-provenance.json").write_text(
            json.dumps(record, sort_keys=True) + "\n"
        )
    _write_receipt_path(output)
    return record


def diagnose_discovery(*, bundle, output, family_id=None, site_id=None, argv=None):
    """Replay one source family or site at its original anchor(s).

    Returns a report dictionary.  Invalid/unsupported requests produce a
    failed report and failure provenance in a new output directory. Existing
    output directories and paths inside the source bundle are rejected before
    any filesystem mutation.
    """
    bundle = Path(bundle).resolve()
    output = Path(output).resolve()
    if argv is None:
        argv = [
            "primalscheme3",
            "panel-discovery-diagnose",
            "--bundle",
            str(bundle),
            "--output",
            str(output),
        ]
        if family_id is not None:
            argv.extend(("--family-id", family_id))
        elif site_id is not None:
            argv.extend(("--site-id", site_id))
    else:
        argv = list(argv)
    if output.exists():
        raise ValueError("diagnostic output already exists: " + str(output))
    if output == bundle or output.is_relative_to(bundle):
        raise ValueError("diagnostic output must be separate from source bundle")
    started = monotonic()
    source_inputs: list[dict[str, Any]] = []
    execution_start = None
    report: dict[str, Any] = {"schemaVersion": SCHEMA, "valid": False}
    try:
        provenance, optimizer, source_catalog, targets, source_inputs = (
            _source_preflight(bundle)
        )
        execution_start = capture_execution_identity(
            [item["path"] for item in source_inputs]
        )
        entity, indexes, target_ids = _entity_binding(
            source_catalog, family_id, site_id
        )
        config = _config_from_bundle(bundle, optimizer)
        profile_binding = _profile_binding(source_catalog, config)
        target = tuple(item for item in targets if item.id in target_ids)
        all_profiles = discovery_profiles(config)
        profiles = {
            name: all_profiles[name]
            for name in profile_binding["profiles"]
            if name in all_profiles
        }
        if set(profiles) != set(profile_binding["profiles"]):
            raise ValueError("source catalog requests an unavailable discovery profile")
        output.mkdir(parents=True)
        history = SQLiteCoverageHistory(
            output / _HISTORY, run_id="discovery-diagnostic"
        )
        try:
            replay_catalog = _invoke_discovery(
                target, config, profiles=profiles, indexes=indexes, history=history
            )
            history.checkpoint()
            history_counts = {
                name: len(getattr(history, name)) for name in history._streams
            }
        finally:
            history.close()
        # Compare scientific membership only. Detailed evidence IDs and event
        # chronology are intentionally excluded from parity.
        if entity["kind"] == "family":
            expected = source_catalog.family_by_id[entity["id"]]
            replayed = replay_catalog.family_by_id.get(expected.id)
            membership_match = replayed is not None and _family_signature(
                replayed, replay_catalog
            ) == _family_signature(expected, source_catalog)
            source_entity = _family_signature(expected, source_catalog)
            replay_entity = (
                _family_signature(replayed, replay_catalog)
                if replayed is not None
                else None
            )
        else:
            expected = source_catalog.site_by_id[entity["id"]]
            replayed = replay_catalog.site_by_id.get(expected.id)
            membership_match = replayed is not None and _site_signature(
                replayed
            ) == _site_signature(expected)
            source_entity = _site_report(expected)
            replay_entity = _site_report(replayed) if replayed is not None else None
        report = {
            "schemaVersion": SCHEMA,
            "valid": bool(membership_match),
            "status": "success" if membership_match else "failure",
            "source": {
                "bundle": str(bundle),
                "catalogSemanticDigest": source_catalog.semantic_digest,
                "provenanceSchema": provenance.get("schemaVersion"),
                "primaryTier": optimizer.get("primaryTier"),
                "inputs": [
                    {
                        "storedPath": str(item["path"].relative_to(bundle)),
                        "sourcePath": item.get("sourcePath", str(item["path"])),
                        "sha256": item.get("sha256"),
                        "size": item.get("size"),
                    }
                    for item in source_inputs
                ],
            },
            "entity": entity,
            "scope": {
                "kind": "reconstructed-anchored-discovery",
                "label": "reconstructed anchored discovery; not the original historical trace",
                "historicalTrace": False,
                "bounded": True,
                "allOriginalTargetRows": True,
                "targetRowCount": sum(len(item.rows) for item in target),
                "targetIds": list(target_ids),
                "anchors": {"forward": list(indexes[0]), "reverse": list(indexes[1])},
                "profiles": profile_binding["profiles"],
            },
            "scientific": {
                "membershipMatch": bool(membership_match),
                "source": source_entity,
                "replayed": replay_entity,
                "comparison": "scientific site/profile membership; evidence IDs and chronology ignored",
            },
            "history": {
                "path": _HISTORY,
                "detail": "full",
                "records": history_counts,
            },
            "error": None
            if membership_match
            else "replayed requested entity scientific membership differs from source",
        }
        _gzip_write(output / _CATALOG, replay_catalog.to_dict())
        _json_write(output / _REPORT, report)
        normalized_inputs = []
        for index, item in enumerate(source_inputs):
            stored = f"source/{index:04d}-{item['path'].name}"
            destination = output / stored
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item["path"], destination)
            normalized_inputs.append(
                {
                    "sourcePath": str(item["path"]),
                    "storedPath": stored,
                    "sourceIndex": item["sourceIndex"],
                }
            )
        receipt = finalize_provenance(
            output_dir=output,
            argv=argv,
            resolved_options={
                "family_id": family_id,
                "site_id": site_id,
                "history_detail": "full",
                "sourceResolvedOptions": optimizer.get("options", {}),
                "replayConfig": config.to_dict() if hasattr(config, "to_dict") else {},
            },
            inputs=normalized_inputs,
            started_at=started,
            ended_at=monotonic(),
            status="success" if membership_match else "failure",
            exit_status=0 if membership_match else 1,
            stderr="" if membership_match else report["error"],
            scientific={
                "schemaVersion": SCHEMA,
                "scope": report["scope"],
                "membershipMatch": bool(membership_match),
            },
            execution_start=execution_start,
        )
        if not _receipt_is_stable_success(receipt) and membership_match:
            report["valid"] = False
            report["status"] = "failure"
            report["scientific"]["membershipMatch"] = False
            report["error"] = "source/runtime/input stability changed during replay"
            _json_write(output / _REPORT, report)
            stored_receipt = json.loads((output / "panel-provenance.json").read_text())
            stored_receipt.update(
                {
                    "status": "failure",
                    "exitStatus": 1,
                    "stderr": report["error"],
                }
            )
            (output / "panel-provenance.json").write_text(
                json.dumps(stored_receipt, sort_keys=True, separators=(",", ":")) + "\n"
            )
        _write_receipt_path(output)
        return report
    except BaseException as error:
        if not source_inputs:
            try:
                provenance = _read(bundle / "panel-provenance.json")
                source_inputs = _input_records(bundle, provenance)
            except Exception:
                source_inputs = []
        report = {
            "schemaVersion": SCHEMA,
            "valid": False,
            "status": "failure",
            "error": str(error),
            "scope": {
                "kind": "reconstructed-anchored-discovery",
                "label": "reconstructed anchored discovery; not the original historical trace",
                "historicalTrace": False,
            },
        }
        output.mkdir(parents=True, exist_ok=True)
        _json_write(output / _REPORT, report)
        _failure_receipt(
            output,
            argv=argv,
            error=error,
            inputs=source_inputs,
            started=started,
            execution_start=execution_start,
        )
        return report


__all__ = ["diagnose_discovery", "SCHEMA"]
