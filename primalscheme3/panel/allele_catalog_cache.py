"""Portable immutable discovery origins, independent of current selection verdicts."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shlex
import shutil
import sqlite3
import tempfile
import zlib
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from importlib.metadata import version
from pathlib import Path
from time import monotonic
from types import SimpleNamespace

from primalscheme3.core.digestion import VARIANT_FREQUENCY_POLICY
from primalscheme3.core.variant_thermo import _FIELDS

from .coverage_discovery import (
    HISTORY_ASSESSMENT_POLICY,
    HISTORY_ORIGIN_POLICY,
    discovery_profiles,
    variant_targets,
)
from .coverage_history import IntrinsicEvidence, StageSnapshot
from .coverage_provenance import runtime_identity, source_identity
from .coverage_types import CandidateFamily, VariantCatalog, canonical_json

SCHEMA = "primalscheme3.discovery-cache/v1"
_DISCOVERY_FILES = tuple(
    "primalscheme3/" + name
    for name in (
        "core/config.py",
        "core/digestion.py",
        "core/variant_thermo.py",
        "core/thermo.py",
        "core/seq_functions.py",
        "core/parallel_discovery.py",
        "core/msa.py",
        "core/mapping.py",
        "panel/coverage_discovery.py",
        "panel/allele_coverage.py",
        "panel/coverage_types.py",
    )
)


def discovery_settings(config, profiles, *, indexes=None, length_mode=None):
    length_mode = length_mode or getattr(
        config, "discovery_length_mode", "first-compatible"
    )
    return {
        "profiles": {
            name: {k: getattr(p, k) for k in _FIELDS}
            for name, p in sorted(profiles.items())
        },
        "minimum_frequency": {
            name: p.min_base_freq for name, p in sorted(profiles.items())
        },
        "frequency_policy": VARIANT_FREQUENCY_POLICY,
        "history_origin_policy": HISTORY_ORIGIN_POLICY,
        "history_assessment_policy": HISTORY_ASSESSMENT_POLICY,
        "maximum_alignment_walk": {
            name: p.primer_max_walk for name, p in sorted(profiles.items())
        },
        "enumeration": "row-anchor-profile-length/v2",
        "discovery_length_mode": length_mode,
        "ambiguous_expansion_limit": 256,
        "amplicon_size_min": config.amplicon_size_min,
        "amplicon_size_max": config.amplicon_size_max,
        "indexes": indexes,
    }


@dataclass(frozen=True)
class CacheReuse:
    catalog: VariantCatalog
    manifest: dict
    source_directory: Path
    derived_evidence: tuple[IntrinsicEvidence, ...] = ()


def _descriptor(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"sha256": digest.hexdigest(), "size": Path(path).stat().st_size}


def _local(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if path == root or root not in path.parents:
        raise ValueError("artifact path escapes cache")
    return path


def _verify(root, relative, expected):
    path = _local(root, relative)
    if not path.is_file() or _descriptor(path) != {
        k: expected[k] for k in ("sha256", "size")
    }:
        raise ValueError("artifact checksum mismatch: " + relative)
    return path


def _read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        return json.load(stream)


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(canonical_json(value) + "\n")


def _scientific_identity(source, runtime):
    files = {
        item["path"]: {k: item[k] for k in ("sha256", "size")}
        for item in source["files"]
    }
    return {
        "files": {name: files[name] for name in _DISCOVERY_FILES},
        "kernels": [
            {k: item[k] for k in ("package", "distribution", "version")}
            | {"files": [{k: f[k] for k in ("sha256", "size")} for f in item["files"]]}
            for item in runtime["nativeKernels"]
        ],
    }


@contextmanager
def _connection(path):
    # Immutable read-only handle: no migration, journal, record replay or writes.
    connection = sqlite3.connect(
        Path(path).resolve().as_uri() + "?mode=ro&immutable=1", uri=True
    )
    try:
        yield connection
    finally:
        connection.close()


def _discovery_snapshot(path, catalog):
    if any(Path(str(path) + suffix).exists() for suffix in ("-wal", "-journal")):
        raise ValueError("origin history must be closed without journal sidecars")
    with _connection(path) as db:
        metadata = dict(db.execute("SELECT key,value FROM metadata"))
        if metadata.get("schema") != "primalscheme3.sqlite-history/v1":
            raise ValueError("unsupported origin history schema")
        found = None
        for (payload,) in db.execute(
            "SELECT payload FROM records WHERE stream='snapshots' AND stage_id='discovery' ORDER BY ordinal"
        ):
            snapshot = StageSnapshot.from_dict(json.loads(zlib.decompress(payload)))
            if (
                snapshot.catalog_digest == catalog.semantic_digest
                and snapshot.completeness == "complete"
            ):
                found = snapshot
        if found is None:
            raise ValueError("complete discovery snapshot missing")
        counts = json.loads(metadata["committed_counts"])
        if any(
            counts[name] < getattr(found, name + "_count")
            for name in ("evidence", "assessments", "events")
        ):
            raise ValueError("discovery snapshot exceeds committed history")
        return found


def _targets_from_inputs(root, inputs):
    from primalscheme3.core.mapping import create_mapping
    from primalscheme3.core.msa import parse_msa

    sources = {}
    for index, item in enumerate(inputs):
        array, _ = parse_msa(_local(root, item["storedPath"]))
        mapping, _ = create_mapping(array)
        sources[index] = SimpleNamespace(
            array=array, _mapping_array=mapping, msa_index=item["sourceIndex"]
        )
    return variant_targets(sources)


def _target_signature(targets):
    # Explicit source order is not the catalog's identity-sorted presentation.
    return [asdict(t) for t in sorted(targets, key=lambda t: t.source_msa_index)]


def _export_panel_discovery_cache(
    panel_dir, cache_dir, *, argv, executed_source, executed_runtime
):
    """Export a closed successful native panel, never a partial/in-flight run."""
    started = monotonic()
    panel_dir, cache_dir = Path(panel_dir).resolve(), Path(cache_dir).resolve()
    source_receipt_bytes = (panel_dir / "panel-provenance.json").read_bytes()
    provenance = json.loads(source_receipt_bytes)
    if (
        provenance.get("status") != "success"
        or provenance.get("exitStatus") != 0
        or provenance.get("sourceChangedDuringRun") is not False
        or provenance.get("runtimeChangedDuringRun") is not False
        or provenance.get("scientific", {}).get("validationValid") is not True
        or provenance.get("scientific", {}).get("selectionAlgorithm")
        != "allele-coverage"
    ):
        raise ValueError("cache requires stable successful allele panel provenance")
    if cache_dir.exists():
        raise FileExistsError(cache_dir)
    outputs = {item["path"]: item for item in provenance["outputs"]}
    catalog_path = "stages/strict/catalog.json.gz"
    required = (
        catalog_path,
        "history/history.sqlite",
        "configuration-ledger.json.gz",
        "stages/strict/validation.json",
    )
    for name in required:
        _verify(panel_dir, name, outputs[name])
    if _read(panel_dir / "stages/strict/validation.json").get("valid") is not True:
        raise ValueError("strict publication was not validated")
    catalog = VariantCatalog.from_dict(_read(panel_dir / catalog_path))
    snapshot = _discovery_snapshot(panel_dir / "history/history.sqlite", catalog)
    for item in provenance["inputs"]:
        if (
            item.get("storedMissing")
            or item.get("sourceAtStartMatchesStored") is not True
        ):
            raise ValueError("source input copy was not authoritative")
        _verify(panel_dir, item["storedPath"], item)
    if _target_signature(
        _targets_from_inputs(panel_dir, provenance["inputs"])
    ) != _target_signature(catalog.targets):
        raise ValueError("source catalog targets disagree with raw inputs")
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    pending = Path(
        tempfile.mkdtemp(
            prefix="." + cache_dir.name + "-pending-", dir=cache_dir.parent
        )
    )
    try:
        (pending / "history").mkdir()
        (pending / "inputs").mkdir()
        shutil.copyfile(panel_dir / catalog_path, pending / "catalog.json.gz")
        shutil.copyfile(
            panel_dir / "history/history.sqlite", pending / "history/history.sqlite"
        )
        shutil.copyfile(
            panel_dir / "configuration-ledger.json.gz",
            pending / "configuration-ledger.json.gz",
        )
        (pending / "provenance.json").write_bytes(source_receipt_bytes)
        inputs = []
        for index, item in enumerate(provenance["inputs"]):
            relative = f"inputs/{index:04d}-" + Path(item["storedPath"]).name
            shutil.copyfile(_local(panel_dir, item["storedPath"]), pending / relative)
            inputs.append(
                {
                    "storedPath": relative,
                    "sourceIndex": item["sourceIndex"],
                    **_descriptor(pending / relative),
                }
            )
        manifest = {
            "schemaVersion": SCHEMA,
            "catalogPath": "catalog.json.gz",
            "historyPath": "history/history.sqlite",
            "originLedgerPath": "configuration-ledger.json.gz",
            "historyRole": "immutable-origin-history-including-prior-selection; not-current-decisions",
            "sourcePanelProvenance": "provenance.json",
            "exportProvenancePath": "cache-export-provenance.json",
            "inputs": inputs,
            "catalogSemanticDigest": catalog.semantic_digest,
            "discoverySnapshot": snapshot.to_dict(),
            "scientificIdentity": _scientific_identity(
                provenance["source"], provenance["runtime"]
            ),
            "artifacts": {
                p.relative_to(pending).as_posix(): _descriptor(p)
                for p in pending.rglob("*")
                if p.is_file()
            },
            "export": {
                "tool": "primalscheme3 export-panel-discovery-cache",
                "argv": list(argv),
                "toolVersion": version("primalscheme3"),
                "workingDirectory": os.getcwd(),
                "stderr": "",
                "resolvedOptions": {
                    "panelDirectory": str(panel_dir),
                    "cacheDirectory": str(cache_dir),
                },
                "sourceProvenanceSha256": _descriptor(
                    panel_dir / "panel-provenance.json"
                )["sha256"],
                "wallSeconds": monotonic() - started,
                "exitStatus": 0,
                "source": executed_source,
                "runtime": executed_runtime,
            },
        }
        # Verify copies against the original immutable receipt, including a race
        # between the initial verification and copying a source payload.
        for source, dest in (
            (catalog_path, "catalog.json.gz"),
            ("history/history.sqlite", "history/history.sqlite"),
            ("configuration-ledger.json.gz", "configuration-ledger.json.gz"),
        ):
            _verify(pending, dest, outputs[source])
        for copied, original in zip(inputs, provenance["inputs"], strict=True):
            _verify(pending, copied["storedPath"], original)
        _write(pending / "manifest.json", manifest)
        pending.rename(cache_dir)
        return manifest
    except BaseException:
        # Keep staged payloads at the requested final directory so the outer
        # failure receipt describes stable, locally inspectable files.
        if pending.exists() and not cache_dir.exists():
            pending.rename(cache_dir)
        raise


def export_panel_discovery_cache(panel_dir, cache_dir, *, argv):
    """Export with success/failure provenance, including preflight failures."""
    panel_dir, cache_dir = Path(panel_dir).resolve(), Path(cache_dir).resolve()
    if cache_dir.exists():
        raise FileExistsError(cache_dir)
    if (
        panel_dir == cache_dir
        or panel_dir in cache_dir.parents
        or cache_dir in panel_dir.parents
    ):
        raise ValueError("cache destination must be separate from source panel")
    started = monotonic()
    executed_source, executed_runtime = source_identity(), runtime_identity()
    source_receipt = panel_dir / "panel-provenance.json"
    inputs = (
        [{"path": str(source_receipt), **_descriptor(source_receipt)}]
        if source_receipt.is_file()
        else []
    )
    result, failure = None, None
    try:
        result = _export_panel_discovery_cache(
            panel_dir,
            cache_dir,
            argv=argv,
            executed_source=executed_source,
            executed_runtime=executed_runtime,
        )
        return result
    except BaseException as error:
        failure = error
        raise
    finally:
        end_source, end_runtime = source_identity(), runtime_identity()
        source_changed = executed_source["sourceDigest"] != end_source["sourceDigest"]
        runtime_changed = executed_runtime != end_runtime
        drift_failure = failure is None and (source_changed or runtime_changed)
        if drift_failure:
            failure = ValueError(
                "source or runtime identity changed during cache export"
            )
        stability = {
            "sourceAtEnd": end_source,
            "runtimeAtEnd": end_runtime,
            "sourceChangedDuringRun": source_changed,
            "runtimeChangedDuringRun": runtime_changed,
        }
        cache_dir.mkdir(parents=True, exist_ok=True)
        receipt = cache_dir / "cache-export-provenance.json"
        if result is not None:
            result["export"].update(
                stability,
                status="success" if failure is None else "failure",
                exitStatus=0 if failure is None else 1,
                stderr="" if failure is None else str(failure),
            )
            _write(cache_dir / "manifest.json", result)
        if result is not None:
            # These copied payload descriptors were verified against the source
            # receipt by the exporter. Do not reread a multi-GB DB once more.
            outputs = dict(result["artifacts"]) | {
                "manifest.json": _descriptor(cache_dir / "manifest.json")
            }
            consumed = _read(cache_dir / "provenance.json")
            by_path = {item["path"]: item for item in consumed["outputs"]}
            for relative in (
                "stages/strict/catalog.json.gz",
                "history/history.sqlite",
                "configuration-ledger.json.gz",
                "stages/strict/validation.json",
            ):
                inputs.append(
                    {
                        "path": str(panel_dir / relative),
                        **{k: by_path[relative][k] for k in ("sha256", "size")},
                    }
                )
            inputs.extend(
                {
                    "path": str(_local(panel_dir, item["storedPath"])),
                    **{k: item[k] for k in ("sha256", "size")},
                }
                for item in consumed["inputs"]
            )
        else:
            outputs = {
                p.relative_to(cache_dir).as_posix(): _descriptor(p)
                for p in cache_dir.rglob("*")
                if p.is_file() and p != receipt
            }
        _write(
            receipt,
            {
                "schemaVersion": "primalscheme3.discovery-cache-export-provenance/v1",
                "tool": "primalscheme3 export-panel-discovery-cache",
                "toolVersion": version("primalscheme3"),
                "argv": list(argv),
                "shell": shlex.join(argv),
                "workingDirectory": os.getcwd(),
                "resolvedOptions": {
                    "panelDirectory": str(panel_dir),
                    "cacheDirectory": str(cache_dir),
                },
                "source": executed_source,
                "runtime": executed_runtime,
                **stability,
                "inputs": inputs,
                "outputs": outputs,
                "outputDirectory": str(cache_dir),
                "status": "success" if failure is None else "failure",
                "exitStatus": 0 if failure is None else 1,
                "wallSeconds": monotonic() - started,
                "stderr": "" if failure is None else str(failure),
                "selfHashPolicy": "cache-export-provenance.json excluded; manifest references it without a checksum to avoid a cycle",
            },
        )

        if drift_failure:
            raise failure


def _project(catalog, names, database):
    resolved = json.loads(catalog.resolved_config_json)
    if names == set(resolved["profiles"]):
        return catalog, ()
    sites = []
    with _connection(database) as db:

        @lru_cache(maxsize=4096)
        def belongs(identity):
            # Indexed assessment/evidence links supply profile membership without
            # decoding or replaying potentially millions of row-origin payloads.
            if (
                db.execute(
                    "SELECT 1 FROM records WHERE stream='evidence' AND id=?",
                    (identity,),
                ).fetchone()
                is None
            ):
                raise ValueError("catalog evidence missing from origin history")
            return (
                db.execute(
                    "SELECT 1 FROM record_links l JOIN records r ON r.id=l.source_id "
                    "WHERE l.kind='evidence' AND l.target_id=? AND r.stream='assessments' "
                    "AND r.stage_id='discovery' AND r.profile_id IN ("
                    + ",".join("?" for _ in names)
                    + ") LIMIT 1",
                    (identity, *sorted(names)),
                ).fetchone()
                is not None
            )

        for site in catalog.sites:
            if not names.intersection(site.generated_profile_ids):
                continue
            evidence_ids = [
                identity
                for identity in site.intrinsic_evidence_ids
                if belongs(identity)
            ]
            sites.append(
                replace(
                    site,
                    generated_profile_ids=tuple(
                        names.intersection(site.generated_profile_ids)
                    ),
                    accepting_profile_ids=tuple(
                        names.intersection(site.accepting_profile_ids)
                    ),
                    intrinsic_evidence_ids=tuple(evidence_ids),
                )
            )
    families, derived_evidence = [], []
    for target in catalog.targets:
        forward, reverse = defaultdict(list), defaultdict(list)
        for site in sites:
            if (
                site.target_id == target.id
                and site.accepting_profile_ids
                and site.reference_footprint is not None
            ):
                (forward if site.strand == "+" else reverse)[
                    site.alignment_anchor
                ].append(site)
        for fa, fs in sorted(forward.items()):
            for ra, rs in sorted(reverse.items()):
                if not (
                    fa <= ra
                    and any(
                        f.reference_footprint[1] <= r.reference_footprint[0]
                        and resolved["amplicon_size_min"]
                        <= r.reference_footprint[1] - f.reference_footprint[0]
                        <= resolved["amplicon_size_max"]
                        for f in fs
                        for r in rs
                    )
                ):
                    continue
                family = CandidateFamily(
                    target.id,
                    (fa, ra),
                    tuple(s.id for s in fs),
                    tuple(s.id for s in rs),
                    tuple(
                        (fp, rp)
                        for fp in sorted(
                            {p for s in fs for p in s.accepting_profile_ids}
                        )
                        for rp in sorted(
                            {p for s in rs for p in s.accepting_profile_ids}
                        )
                    ),
                )
                evidence = IntrinsicEvidence(
                    (family.id,),
                    "family-geometry-existence",
                    {
                        "sites": sorted((s.id, s.reference_footprint) for s in fs + rs),
                        "minimum": resolved["amplicon_size_min"],
                        "maximum": resolved["amplicon_size_max"],
                    },
                    {"some_selection_feasible": True},
                    "evaluated",
                )
                families.append(replace(family, geometry_evidence_ids=(evidence.id,)))
                if (
                    evidence.id
                    not in catalog.family_by_id[family.id].geometry_evidence_ids
                ):
                    derived_evidence.append(evidence)
    for key in ("profiles", "minimum_frequency", "maximum_alignment_walk"):
        resolved[key] = {
            name: value for name, value in resolved[key].items() if name in names
        }
    return replace(
        catalog,
        sites=tuple(sites),
        families=tuple(families),
        resolved_config_json=canonical_json(resolved),
    ), tuple(derived_evidence)


def load_discovery_cache(cache_dir, *, targets, config, profiles=None, indexes=None):
    cache_dir = Path(cache_dir).resolve()
    manifest = _read(cache_dir / "manifest.json")
    if manifest.get("schemaVersion") != SCHEMA:
        raise ValueError("unsupported discovery cache schema")
    for name, descriptor in manifest["artifacts"].items():
        _verify(cache_dir, name, descriptor)
    export_receipt = _read(_local(cache_dir, manifest["exportProvenancePath"]))
    if (
        export_receipt.get("status") != "success"
        or export_receipt.get("exitStatus") != 0
        or export_receipt.get("sourceChangedDuringRun") is not False
        or export_receipt.get("runtimeChangedDuringRun") is not False
    ):
        raise ValueError("cache export did not complete successfully")
    _verify(cache_dir, "manifest.json", export_receipt["outputs"]["manifest.json"])
    if any(
        export_receipt["outputs"].get(name) != descriptor
        for name, descriptor in manifest["artifacts"].items()
    ):
        raise ValueError("cache artifacts disagree with export receipt")

    for key in (
        "catalogPath",
        "historyPath",
        "originLedgerPath",
        "sourcePanelProvenance",
    ):
        if manifest[key] not in manifest["artifacts"]:
            raise ValueError("required cache artifact not fingerprinted")
    if manifest["scientificIdentity"] != _scientific_identity(
        source_identity(), runtime_identity()
    ):
        raise ValueError("discovery implementation or kernels changed")
    source_provenance = _read(_local(cache_dir, manifest["sourcePanelProvenance"]))
    if (
        source_provenance.get("status") != "success"
        or source_provenance.get("exitStatus") != 0
        or source_provenance.get("sourceChangedDuringRun") is not False
        or source_provenance.get("runtimeChangedDuringRun") is not False
        or source_provenance.get("scientific", {}).get("validationValid") is not True
        or _scientific_identity(
            source_provenance["source"], source_provenance["runtime"]
        )
        != manifest["scientificIdentity"]
    ):
        raise ValueError("cache lacks stable successful source provenance")
    source_outputs = {item["path"]: item for item in source_provenance["outputs"]}
    for cached, original in (
        (manifest["catalogPath"], "stages/strict/catalog.json.gz"),
        (manifest["historyPath"], "history/history.sqlite"),
        (manifest["originLedgerPath"], "configuration-ledger.json.gz"),
    ):
        if manifest["artifacts"][cached] != {
            k: source_outputs[original][k] for k in ("sha256", "size")
        }:
            raise ValueError("cache artifact differs from source receipt")
    if len(manifest["inputs"]) != len(source_provenance["inputs"]):
        raise ValueError("cache input occurrence count differs")
    for item, original in zip(
        manifest["inputs"], source_provenance["inputs"], strict=True
    ):
        if (
            item["storedPath"] not in manifest["artifacts"]
            or item["sourceIndex"] != original["sourceIndex"]
        ):
            raise ValueError("cache input fingerprints/order differ")
        _verify(cache_dir, item["storedPath"], item)
        _verify(cache_dir, item["storedPath"], original)
    catalog = VariantCatalog.from_dict(
        _read(_local(cache_dir, manifest["catalogPath"]))
    )
    if catalog.semantic_digest != manifest["catalogSemanticDigest"]:
        raise ValueError("catalog semantic digest mismatch")
    if isinstance(targets, dict):
        targets = variant_targets(targets)
    if _target_signature(targets) != _target_signature(
        catalog.targets
    ) or _target_signature(
        _targets_from_inputs(cache_dir, manifest["inputs"])
    ) != _target_signature(catalog.targets):
        raise ValueError("fresh targets differ from cache inputs/order/rows")
    if profiles is None:
        profiles = discovery_profiles(config)
        choice = getattr(config, "candidate_profiles", "union")
        if choice != "union":
            profiles = {choice: profiles[choice]}
    requested = discovery_settings(config, profiles, indexes=indexes)
    available = json.loads(catalog.resolved_config_json)
    names = set(profiles)
    if not names or not names <= set(available["profiles"]):
        raise ValueError("requested discovery profiles not available")
    for key in ("profiles", "minimum_frequency", "maximum_alignment_walk"):
        available[key] = {
            name: value for name, value in available[key].items() if name in names
        }
    if canonical_json(available) != canonical_json(requested):
        raise ValueError("discovery settings differ from cache")
    database = _local(cache_dir, manifest["historyPath"])
    if (
        _discovery_snapshot(database, catalog).to_dict()
        != manifest["discoverySnapshot"]
    ):
        raise ValueError("discovery snapshot mismatch")
    projected, derived = _project(catalog, names, database)
    return CacheReuse(projected, manifest, cache_dir, derived)


def materialize_cache_reuse(reuse, output_dir, *, history):
    destination = Path(output_dir) / "origin-discovery"
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Copied bytes are verified against the loaded receipt, catching source
    # changes after load without an additional whole-history pre-copy read.
    if _read(reuse.source_directory / "manifest.json") != reuse.manifest:
        raise ValueError("cache manifest changed after validation")
    destination.mkdir()
    for name in (
        *reuse.manifest["artifacts"],
        "manifest.json",
        reuse.manifest["exportProvenancePath"],
    ):
        target = _local(destination, name)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(_local(reuse.source_directory, name), target)
    for name, descriptor in reuse.manifest["artifacts"].items():
        _verify(destination, name, descriptor)
    receipt = _read(_local(destination, reuse.manifest["exportProvenancePath"]))
    if (
        receipt.get("status") != "success"
        or receipt.get("exitStatus") != 0
        or receipt.get("sourceChangedDuringRun") is not False
        or receipt.get("runtimeChangedDuringRun") is not False
    ):
        raise ValueError("copied cache export receipt is invalid")
    _verify(destination, "manifest.json", receipt["outputs"]["manifest.json"])
    for evidence in reuse.derived_evidence:
        recorded = history.record_evidence(
            entity_ids=evidence.entity_ids,
            measurement=evidence.measurement,
            dependency_key=evidence.dependency_key,
            values=evidence.values,
            status=evidence.status,
        )
        if recorded.id != evidence.id:
            raise ValueError("derived geometry evidence identity changed")
    history.emit(
        stage_id="discovery-reuse",
        kind="reused-discovery-cache",
        entity_ids=tuple(t.id for t in reuse.catalog.targets)
        + tuple(e.entity_ids[0] for e in reuse.derived_evidence),
        changes={
            "origin_manifest": "origin-discovery/manifest.json",
            "origin_catalog_digest": reuse.manifest["catalogSemanticDigest"],
            "catalog_digest": reuse.catalog.semantic_digest,
            "profiles": sorted(
                json.loads(reuse.catalog.resolved_config_json)["profiles"]
            ),
            "origin_history_role": reuse.manifest["historyRole"],
        },
    )
    return destination
