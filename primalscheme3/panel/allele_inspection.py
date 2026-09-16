"""Read-only bounded history queries and independent portable bundle audits."""

from __future__ import annotations

import json
import os
import shlex
import sqlite3
import time
import zlib
from collections import deque
from importlib.metadata import version
from pathlib import Path

from .allele_coverage import canonical_observations
from .allele_labels import SCHEMA as ALLELE_LABEL_SCHEMA
from .allele_labels import build_allele_label_map
from .allele_publication import _read, _targets, artifact_descriptor, audit_allele_stage
from .coverage_history import (
    CoverageHistory,
    SQLiteCoverageHistory,
    sqlite_history_format,
)
from .coverage_provenance import runtime_identity, source_identity
from .coverage_types import ConfigurationLedger, VariantCatalog, canonical_json


def _contained(root, value):
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("unsafe bundle-relative path: " + str(value))
    result = (root / path).resolve()
    if not result.is_relative_to(root.resolve()):
        raise ValueError("path escapes bundle: " + str(value))
    return result


def _write(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2, default=str) + "\n")


def _catalog_path(bundle, optimizer):
    return _contained(
        bundle, optimizer.get("catalogPath", "stages/strict/catalog.json.gz")
    )


def _context_entities(bundle, optimizer, target, region):
    catalog = VariantCatalog.from_dict(_read(_catalog_path(bundle, optimizer)))
    references = optimizer.get("publication", {}).get("targetToReference", {})
    targets = {
        t.id
        for t in catalog.targets
        if target is None or target in (t.id, references.get(t.id))
    }
    if target is not None and not targets:
        raise ValueError("unmatched target/reference: " + target)
    sites, unmapped = set(), []
    for site in catalog.sites:
        if site.target_id not in targets:
            continue
        if site.reference_footprint is None:
            unmapped.append(site.id)
            if region is None:
                sites.add(site.id)
        elif region is None or (
            site.reference_footprint[0] < region[1]
            and region[0] < site.reference_footprint[1]
        ):
            sites.add(site.id)
    families = set()
    for family in catalog.families:
        if family.target_id not in targets:
            continue
        footprints = [
            catalog.site_by_id[sid].reference_footprint
            for sid in family.forward_site_ids + family.reverse_site_ids
        ]
        mapped = [span for span in footprints if span is not None]
        if region is None or (
            mapped
            and min(p[0] for p in mapped) < region[1]
            and region[0] < max(p[1] for p in mapped)
        ):
            families.add(family.id)
    ledger_path = optimizer.get("history", {}).get(
        "configurationLedger", "configuration-ledger.json.gz"
    )
    ledger = ConfigurationLedger.from_dict(_read(_contained(bundle, ledger_path)))
    if ledger.catalog_digest != catalog.semantic_digest:
        raise ValueError("configuration ledger catalog identity mismatch")
    configs = {
        c.id
        for c in ledger.configurations
        if c.target_id in targets
        and (
            region is None
            or (c.full_interval[0] < region[1] and region[0] < c.full_interval[1])
        )
    }
    entities = sites | families | configs | (targets if region is None else set())
    return entities, {
        "matched_targets": sorted(targets),
        "matched_entities": len(entities),
        "unmapped_site_count": len(unmapped),
        "unmapped_site_ids": unmapped[:100],
        "unmapped_sites_truncated": len(unmapped) > 100,
        "region_coordinates": "zero-based half-open ungapped first reference",
    }


def _query_history(
    bundle,
    *,
    entity=None,
    target=None,
    region=None,
    pool=None,
    stage=None,
    profile=None,
    lineage=False,
    limit=100,
    offset=0,
    optimizer=None,
    scientific_ids=None,
):
    """Read indexed pages; dependency/parent closure is separately bounded.

    limit/offset apply per stream. Dependencies can cross the requested filters;
    they retain their original context and are never relabeled as matching rows.
    """
    bundle = Path(bundle).resolve()
    if isinstance(region, str):
        try:
            region = tuple(int(x) for x in region.split(":"))
        except ValueError as exc:
            raise ValueError("region must be zero-based half-open START:END") from exc
    if not (entity or target or region is not None or stage):
        raise ValueError("an entity, stage, target or region filter is required")
    if (
        type(limit) is not int
        or not 1 <= limit <= 1000
        or type(offset) is not int
        or offset < 0
    ):
        raise ValueError("limit must be 1..1000 and offset nonnegative")
    if pool is not None and (type(pool) is not int or pool < 1):
        raise ValueError("pool is one-based and must be positive")
    if region is not None and (
        len(region) != 2
        or any(type(x) is not int for x in region)
        or not 0 <= region[0] < region[1]
    ):
        raise ValueError("region must be zero-based half-open START:END")
    optimizer = optimizer or _read(_contained(bundle, "panel-optimizer.json"))
    entities, context = None, {}
    if target is not None or region is not None:
        entities, context = _context_entities(bundle, optimizer, target, region)
    if entity is not None:
        entities = {entity} if entities is None else entities.intersection((entity,))
    path = _contained(bundle, optimizer["history"]["path"])
    output = {name: [] for name in CoverageHistory._streams}
    pages, records, direct = {}, {}, set()
    snapshot_roots = []
    max_linked = min(10000, max(100, 10 * limit))
    db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    try:
        db.execute("PRAGMA query_only=ON")
        metadata = dict(db.execute("SELECT key,value FROM metadata"))
        sqlite_history_format(db, metadata)

        def decode(row):
            position, stream, ordinal, identity, payload = row
            raw = zlib.decompress(payload)
            record = CoverageHistory._streams[stream].from_dict(json.loads(raw))
            if record.id != identity or canonical_json(record).encode() != raw:
                raise ValueError("corrupt history semantic identity")
            if hasattr(record, "run_id") and record.run_id != metadata["run_id"]:
                raise ValueError("history run identity mismatch")
            if stream == "events" and record.sequence_number != ordinal:
                raise ValueError("history sequence mismatch")
            stored_context = db.execute(
                "SELECT stage_id,pool,profile_id FROM records WHERE id=?", (record.id,)
            ).fetchone()
            if stored_context != tuple(
                getattr(record, k, None) for k in ("stage_id", "pool", "profile_id")
            ):
                raise ValueError("corrupt history filter index")
            entities = set(
                db.execute(
                    "SELECT entity_id,stream,ordinal FROM entity_records WHERE record_id=?",
                    (record.id,),
                )
            )
            if entities != {
                (e, stream, ordinal) for e in getattr(record, "entity_ids", ())
            }:
                raise ValueError("corrupt history entity index")
            for entity_id in getattr(record, "entity_ids", ()):
                if scientific_ids is not None:
                    for kind, identities in scientific_ids.items():
                        if (
                            entity_id.startswith(kind + "-")
                            and entity_id not in identities
                        ):
                            raise ValueError(
                                "origin history references unknown " + kind
                            )
            links = SQLiteCoverageHistory._links(record)
            actual_links = set(
                db.execute(
                    "SELECT target_id,kind FROM record_links WHERE source_id=?",
                    (record.id,),
                )
            )
            if actual_links != {(identity, kind) for identity, kind, _ in links}:
                raise ValueError("corrupt indexed history links")
            for identity, kind, expected in links:
                linked = db.execute(
                    "SELECT stream,position FROM records WHERE id=?", (identity,)
                ).fetchone()
                if linked is None or linked[0] != expected or linked[1] >= position:
                    raise ValueError("missing/forward history " + kind + " reference")
                if (
                    db.execute(
                        "SELECT 1 FROM record_links WHERE source_id=? AND target_id=? AND kind=?",
                        (record.id, identity, kind),
                    ).fetchone()
                    is None
                ):
                    raise ValueError("missing indexed history link")
            return record, links

        for stream in output:
            terms, args = ["r.stream=?"], [stream]
            if entities is not None and stream != "snapshots":
                # json_each avoids SQL variable limits; entity_records' primary
                # index narrows candidates without decoding unrelated payloads.
                terms.append(
                    "r.id IN (SELECT record_id FROM entity_records WHERE entity_id IN (SELECT value FROM json_each(?)))"
                )
                args.append(json.dumps(sorted(entities)))
            if stage is not None:
                terms.append("r.stage_id=?")
                args.append(stage)
            if stream == "events" and (pool is not None or profile is not None):
                linked_terms = ["l.source_id=r.id", "l.kind='assessment'"]
                for column, value in [
                    ("pool", None if pool is None else pool - 1),
                    ("profile_id", profile),
                ]:
                    if value is not None:
                        linked_terms.append("a." + column + "=?")
                        args.append(value)
                terms.append(
                    "EXISTS(SELECT 1 FROM record_links l JOIN records a ON a.id=l.target_id WHERE "
                    + " AND ".join(linked_terms)
                    + ")"
                )
            if stream == "assessments":
                for column, value in [
                    ("pool", None if pool is None else pool - 1),
                    ("profile_id", profile),
                ]:
                    if value is not None:
                        terms.append("r." + column + "=?")
                        args.append(value)
            sql = (
                "SELECT r.position,r.stream,r.ordinal,r.id,r.payload FROM records r WHERE "
                + " AND ".join(terms)
                + (
                    " ORDER BY r.ordinal DESC"
                    if stream == "snapshots"
                    else " ORDER BY r.ordinal"
                )
                + " LIMIT ? OFFSET ?"
            )
            rows = db.execute(sql, (*args, limit + 1, offset)).fetchall()
            pages[stream] = {
                "limit": limit,
                "offset": offset,
                "returned": min(limit, len(rows)),
                "truncated": len(rows) > limit,
                "next_offset": offset + limit if len(rows) > limit else None,
            }
            if stream == "snapshots":
                pages[stream].update(
                    contextual=True, entity_filter_applied=False, order="newest-first"
                )
            for row in rows[:limit]:
                record, links = decode(row)
                records[record.id] = (row[1], record, links)
                direct.add(record.id)
                if stream == "snapshots":
                    snapshot_roots.append(record.id)
        pending = deque(records)
        truncated = False
        while pending:
            source = pending.popleft()
            for identity, kind, expected in sorted(records[source][2]):
                if kind in ("parent", "prior") and not lineage:
                    continue
                if identity in records:
                    continue
                if len(records) - len(direct) >= max_linked:
                    truncated = True
                    continue
                row = db.execute(
                    "SELECT position,stream,ordinal,id,payload FROM records WHERE id=?",
                    (identity,),
                ).fetchone()
                record, links = decode(row)
                records[identity] = (expected, record, links)
                pending.append(identity)
        # Snapshots have disposition keys rather than entity-index entries.
        # Preserve complete canonical records, and expose a separate filtered
        # view; never assign the original semantic ID to a changed payload.
        disposition_views = []
        for root_id in snapshot_roots:
            current_id = root_id
            dispositions, origins = {}, {}
            while current_id in records:
                snapshot = records[current_id][1]
                for identity, disposition in snapshot.dispositions.items():
                    if (
                        entities is None or identity in entities
                    ) and identity not in dispositions:
                        dispositions[identity] = disposition
                        origins[identity] = snapshot.id
                current_id = snapshot.prior_snapshot_id
                if current_id is None:
                    break
            disposition_views.append(
                {
                    "snapshot_id": root_id,
                    "dispositions": dispositions,
                    "source_snapshot_ids": origins,
                    "inheritance_complete": current_id is None,
                    "unresolved_prior_snapshot_id": current_id,
                }
            )
        for stream, record, _ in records.values():
            output[stream].append(record.to_dict())
        for stream in output:
            output[stream].sort(key=lambda r: (r.get("sequence_number", 0), r["id"]))
    finally:
        db.close()
    from .coverage_discovery import discovery_condition_rows

    output.update(
        valid=True,
        query=dict(
            entity=entity,
            target=target,
            region=region,
            pool=pool,
            stage=stage,
            profile=profile,
            lineage=lineage,
            limit=limit,
            offset=offset,
        ),
        pagination=pages,
        context=context,
        linked={
            "limit": max_linked,
            "returned": len(records) - len(direct),
            "truncated": truncated,
        },
        scope={
            "full_history_verified": False,
            "verification": "selected canonical records and their reference/index integrity",
            "linked_records_may_cross_filters": True,
            "lineage_direction": "ancestors",
            "snapshot_prefix_hashes_verified": False,
        },
        matched_record_ids=sorted(direct),
        snapshot_dispositions=disposition_views,
    )
    output["discovery_conditions"] = list(discovery_condition_rows(output))
    return output


def _origin_source(bundle, optimizer):
    pointer = optimizer.get("history", {}).get("originDiscovery")
    if not pointer:
        return None
    manifest_path = _contained(bundle, pointer["path"])
    manifest = _read(manifest_path)
    if manifest.get("schemaVersion") != "primalscheme3.discovery-cache/v1":
        raise ValueError("unsupported origin discovery schema")
    root = manifest_path.parent
    keys = ("catalogPath", "historyPath", "originLedgerPath")
    paths = {key: _contained(root, manifest[key]) for key in keys}
    for key, path in paths.items():
        expected = manifest.get("artifacts", {}).get(manifest[key])
        actual = artifact_descriptor(path, root)
        if expected is None or any(
            expected.get(field) != actual[field] for field in ("sha256", "size")
        ):
            raise ValueError("origin artifact integrity: " + key)
    return manifest, paths, manifest_path


def query_allele_history(
    bundle,
    *,
    entity=None,
    target=None,
    region=None,
    pool=None,
    stage=None,
    profile=None,
    lineage=False,
    limit=100,
    offset=0,
):
    """Query current and, when present, immutable origin histories separately."""
    bundle = Path(bundle).resolve()
    options = dict(
        entity=entity,
        target=target,
        region=region,
        pool=pool,
        stage=stage,
        profile=profile,
        lineage=lineage,
        limit=limit,
        offset=offset,
    )
    result = _query_history(bundle, **options)
    optimizer = _read(_contained(bundle, "panel-optimizer.json"))
    origin = _origin_source(bundle, optimizer)
    if origin:
        manifest, paths, _ = origin
        catalog = VariantCatalog.from_dict(_read(paths["catalogPath"]))
        ledger = ConfigurationLedger.from_dict(_read(paths["originLedgerPath"]))
        if (
            catalog.semantic_digest != manifest["catalogSemanticDigest"]
            or ledger.catalog_digest != catalog.semantic_digest
        ):
            raise ValueError("origin catalog/ledger identity mismatch")
        origin_optimizer = {
            "history": {
                "path": str(paths["historyPath"].relative_to(bundle)),
                "configurationLedger": str(
                    paths["originLedgerPath"].relative_to(bundle)
                ),
            },
            "catalogPath": str(paths["catalogPath"].relative_to(bundle)),
            "publication": optimizer.get("publication", {}),
        }
        result["origin_discovery"] = _query_history(
            bundle,
            optimizer=origin_optimizer,
            scientific_ids={
                "OligoSite": set(catalog.site_by_id),
                "CandidateFamily": set(catalog.family_by_id),
                "SelectionConfiguration": set(ledger.configuration_by_id),
            },
            **options,
        )
        result["origin_discovery"]["scope"].update(
            current_decisions=False, history_role=manifest["historyRole"]
        )
    return result


def audit_allele_bundle(bundle, *, tier=None):
    """Audit stored originals and all successful stages; never use cached validity."""
    from .allele_pipeline import _authoritative_targets

    bundle = Path(bundle).resolve()
    provenance = _read(_contained(bundle, "panel-provenance.json"))
    optimizer = _read(_contained(bundle, "panel-optimizer.json"))
    violations, reports = [], {}
    config = _read(_contained(bundle, "config.json"))
    stages = optimizer.get("stages", [])
    if len({s["stage_id"] for s in stages}) != len(stages):
        violations.append({"reason": "duplicate-stage-identities"})
    selected = [s for s in stages if tier is None or s["stage_id"] == tier]
    if not selected:
        raise ValueError("requested tier does not exist")
    primary = optimizer.get("primaryTier")
    primary_stage = next((s for s in stages if s["stage_id"] == primary), None)
    if primary_stage is not None and primary_stage not in selected:
        selected.append(primary_stage)  # Primary copies always require fresh science.
    if provenance.get("exitStatus") != 0 or provenance.get("status") != "success":
        violations.append({"reason": "source-workflow-not-successful"})
    if any(
        not provenance.get(key)
        for key in ("source", "sourceAtEnd", "runtime", "runtimeAtEnd")
    ):
        violations.append({"reason": "missing-execution-identity"})
    for flag in ("sourceChangedDuringRun", "runtimeChangedDuringRun"):
        if provenance.get(flag) is not False:
            violations.append({"reason": flag})
    if provenance.get("source", {}).get("sourceDigest") != provenance.get(
        "sourceAtEnd", {}
    ).get("sourceDigest") or provenance.get("runtime") != provenance.get(
        "runtimeAtEnd"
    ):
        violations.append({"reason": "execution-identity-mismatch"})
    descriptors = provenance.get("outputs", [])
    paths = set()
    for item in descriptors:
        name = item["path"]
        if name in paths:
            violations.append({"reason": "duplicate-output-descriptor", "path": name})
        paths.add(name)
        path = _contained(bundle, name)
        if not path.is_file() or artifact_descriptor(path, bundle) != item:
            violations.append({"reason": "output-integrity", "path": name})
    required = {
        "panel-optimizer.json",
        "config.json",
        "configuration-ledger.json.gz",
        optimizer["history"]["path"],
    }
    if not required.issubset(paths):
        violations.append(
            {
                "reason": "missing-required-output-descriptors",
                "paths": sorted(required - paths),
            }
        )
    inputs = provenance.get("inputs", [])
    if not inputs:
        violations.append({"reason": "missing-original-inputs"})
    indexes = [item.get("sourceIndex", i) for i, item in enumerate(inputs)]
    if any(type(i) is not int or i < 0 for i in indexes) or len(set(indexes)) != len(
        indexes
    ):
        violations.append({"reason": "input-source-order-invalid"})
    msa_data = config.get("msa_data", {})
    if {str(i) for i in indexes} != set(msa_data):
        violations.append({"reason": "config-input-source-order-mismatch"})
    for i, item in enumerate(inputs):
        if (
            msa_data.get(str(item.get("sourceIndex", i)), {}).get("msa_path")
            != item["storedPath"]
        ):
            violations.append({"reason": "config-input-path-mismatch"})
        path = _contained(bundle, item["storedPath"])
        if not path.is_file():
            violations.append(
                {"reason": "missing-stored-input", "path": item["storedPath"]}
            )
            continue
        desc = artifact_descriptor(path, bundle)
        if (
            desc["sha256"] != item.get("sourceAtStartSha256")
            or desc["size"] != item.get("sourceAtStartSize")
            or desc["sha256"] != item.get("sha256")
            or desc["size"] != item.get("size")
        ):
            violations.append(
                {"reason": "original-input-integrity", "path": item["storedPath"]}
            )
    try:
        targets = _authoritative_targets(bundle, inputs)
        raw_reparsed = True
    except Exception as exc:
        targets = ()
        raw_reparsed = False
        violations.append({"reason": "raw-input-reparse", "error": str(exc)})
    label_reference = optimizer.get("publication", {}).get("alleleLabelMap")
    label_report = {"advertised": False, "valid": None, "path": None}
    if label_reference is not None:
        label_report = {"advertised": True, "valid": False, "path": None}
        try:
            if not isinstance(label_reference, dict):
                raise ValueError("reference must be an object")
            label_path = label_reference.get("path")
            if (
                label_reference.get("schemaVersion") != ALLELE_LABEL_SCHEMA
                or not isinstance(label_path, str)
                or not label_path
            ):
                raise ValueError("reference schema/path is invalid")
            label_report["path"] = label_path
            if label_path not in paths:
                raise ValueError("artifact is absent from panel output descriptors")
            stored_labels = _read(_contained(bundle, label_path))
            if stored_labels.get("schemaVersion") != ALLELE_LABEL_SCHEMA:
                raise ValueError("stored schema is invalid")
            observations = tuple(
                observation
                for target in targets
                for observation in canonical_observations(target)
            )
            expected_labels = build_allele_label_map(
                targets, observations, bundle, inputs
            )
            if stored_labels != expected_labels:
                raise ValueError("saved labels differ from freshly parsed raw inputs")
            label_report.update(
                valid=True,
                targets=len(stored_labels["targets"]),
                rows=sum(len(target["rows"]) for target in stored_labels["targets"]),
                classes=sum(
                    len(target["classes"]) for target in stored_labels["targets"]
                ),
            )
        except Exception as exc:
            label_report["error"] = str(exc)
            violations.append(
                {"reason": "allele-label-map-validation", "error": str(exc)}
            )
    # Missing historical fields mean the original exact-supported policy.
    optimizer_profile = dict(optimizer.get("profile", {}))
    optimizer_profile.setdefault("intended_product_policy", "exact-supported")
    optimizer_profile.setdefault(
        "secondary_product_policy", "ordered-disjoint-intended-sites"
    )
    effective_intended = optimizer_profile["intended_product_policy"]
    option_sources = [config, optimizer.get("options", {})]
    if config.get("allele_options_json"):
        option_sources.append(json.loads(config["allele_options_json"]))
    if "resolvedOptions" in provenance:
        option_sources.append(provenance["resolvedOptions"])
    scientific = provenance.get("scientific", {})
    for key in ("profile", "resolvedAlleleOptions"):
        if key in scientific:
            option_sources.append(scientific[key])
    if effective_intended not in (
        "exact-supported",
        "concrete-designated-sites",
    ) or any(
        value.get("intended_product_policy", "exact-supported") != effective_intended
        for value in option_sources
    ):
        violations.append({"reason": "intended-product-policy-mismatch"})
    effective_secondary = optimizer_profile["secondary_product_policy"]
    if effective_secondary not in (
        "ordered-disjoint-intended-sites",
        "reject-secondary-products/v1",
        "ordered-disjoint-concrete-designated-sites/v1",
    ) or any(
        value.get("secondary_product_policy", "ordered-disjoint-intended-sites")
        != effective_secondary
        for value in option_sources
    ):
        violations.append({"reason": "secondary-product-policy-mismatch"})
    expected_metric = "observed-allele-primer-trimmed/v1"
    if (
        optimizer.get("metric") != expected_metric
        or optimizer.get("profile", {}).get("name") != "allele-panel-v1"
    ):
        violations.append({"reason": "optimizer-scientific-identity"})
    for stage in selected:
        name = stage["stage_id"]
        directory = _contained(bundle, stage["path"])
        try:
            manifest = _read(directory / "stage.json")
            if directory != _contained(bundle, "stages/" + name):
                raise ValueError("stage directory identity mismatch")
            expected_policy = stage.get("optimizer", {}).get("stage_policy")
            if expected_policy != manifest["stage_policy"]:
                raise ValueError("stage policy disagrees with optimizer")
            if (
                manifest["stage_id"] != name
                or manifest["stage_policy"]["stage_id"] != name
            ):
                raise ValueError("tier identity mismatch")
            catalog = VariantCatalog.from_dict(_read(directory / "catalog.json.gz"))
            saved_targets = _targets(_read(directory / "authoritative-targets.json.gz"))

            # Input order is scientifically meaningful inside each Target's
            # occurrence identity; container order differs because catalogs sort
            # IDs. Compare whole records in one canonical order, not a set that
            # could collapse duplicate source occurrences.
            def canonical_targets(values):
                return tuple(sorted(values, key=lambda t: t.id))

            if canonical_targets(targets) != canonical_targets(
                saved_targets
            ) or canonical_targets(targets) != canonical_targets(catalog.targets):
                violations.append(
                    {"reason": "raw-input-target-mismatch", "stage": name}
                )
            report = audit_allele_stage(directory)
            reports[name] = report
            if not report["valid"]:
                violations.append(
                    {"reason": "stage-scientific-validation", "stage": name}
                )
            if (
                report["profile"] != optimizer_profile
                or report["coverage"]["metric_id"] != expected_metric
            ):
                violations.append(
                    {"reason": "stage-profile-metric-mismatch", "stage": name}
                )
            if "stages/" + name + "/stage.json" not in paths:
                violations.append({"reason": "missing-stage-descriptor", "stage": name})
        except Exception as exc:
            reports[name] = {"valid": False, "error": str(exc)}
            violations.append(
                {"reason": "stage-audit-error", "stage": name, "error": str(exc)}
            )
    primary = optimizer.get("primaryTier")
    primary_stage = next((s for s in stages if s["stage_id"] == primary), None)
    if primary_stage is None:
        violations.append({"reason": "missing-primary-stage"})
    else:
        directory = _contained(bundle, primary_stage["path"])
        for name in (
            "primer.bed",
            "primers.fasta",
            "order-sheet.tsv",
            "amplicon.bed",
            "primertrim.amplicon.bed",
            "reference.fasta",
        ):
            path = _contained(bundle, name)
            if (
                not path.is_file()
                or not (directory / name).is_file()
                or path.read_bytes() != (directory / name).read_bytes()
            ):
                violations.append({"reason": "primary-copy-mismatch", "path": name})
            if name not in paths:
                violations.append(
                    {"reason": "missing-primary-descriptor", "path": name}
                )
    return {
        "valid": not violations,
        "violations": violations,
        "stages": reports,
        "raw_inputs_reparsed": raw_reparsed,
        "primary_tier": primary,
        "allele_label_map": label_report,
        "scope": "stored-original-inputs-and-fresh-selected-stage-kernels",
    }


def run_inspection(command, bundle, output, *, argv, **options):
    """Write independent results and a failure receipt; source bundle is read-only."""
    bundle, output = Path(bundle).resolve(), Path(output).resolve()
    if output == bundle or output.is_relative_to(bundle):
        raise ValueError("inspection output must be outside the source bundle")
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    source = source_identity()
    runtime = runtime_identity()
    inputs = []
    stderr = ""
    try:
        # Descriptors include exact source bytes even on a failed request. Only
        # history query inputs are hashed; audit verifies the full output bundle.
        candidates = [
            _contained(bundle, "panel-provenance.json"),
            _contained(bundle, "panel-optimizer.json"),
        ]
        inputs = [artifact_descriptor(p, bundle) for p in candidates if p.is_file()]
        if (bundle / "panel-optimizer.json").is_file():
            optimizer = _read(bundle / "panel-optimizer.json")
            if command == "panel-history":
                candidates.append(_contained(bundle, optimizer["history"]["path"]))
                origin = _origin_source(bundle, optimizer)
                if origin:
                    _, paths, manifest_path = origin
                    candidates.extend([manifest_path, *paths.values()])
                if (
                    options.get("target") is not None
                    or options.get("region") is not None
                ):
                    candidates.extend(
                        [
                            _catalog_path(bundle, optimizer),
                            _contained(
                                bundle,
                                optimizer["history"].get(
                                    "configurationLedger",
                                    "configuration-ledger.json.gz",
                                ),
                            ),
                        ]
                    )
            else:
                candidates.extend(
                    _contained(bundle, str(p.relative_to(bundle)))
                    for p in bundle.rglob("*")
                    if p.is_file()
                )
        inputs = [
            artifact_descriptor(p, bundle)
            for p in sorted(set(candidates))
            if p.is_file()
        ]
        if command == "panel-history":
            result = query_allele_history(bundle, **options)
        elif command == "panel-audit":
            result = audit_allele_bundle(bundle, **options)
        else:
            raise ValueError("unknown inspection command")
    except (Exception, KeyboardInterrupt) as exc:
        stderr = f"{type(exc).__name__}: {exc}"
        result = {"valid": False, "error": stderr}
    changed_inputs = []
    for item in inputs:
        path = _contained(bundle, item["path"])
        if not path.is_file() or artifact_descriptor(path, bundle) != item:
            changed_inputs.append(item["path"])
    if changed_inputs:
        result["valid"] = False
        result.setdefault("violations", []).append(
            {"reason": "inspection-input-changed", "paths": changed_inputs}
        )
    source_end = source_identity()
    runtime_end = runtime_identity()
    changed = (
        source["sourceDigest"] != source_end["sourceDigest"] or runtime != runtime_end
    )
    if changed:
        result["valid"] = False
        result.setdefault("violations", []).append(
            {"reason": "inspection-execution-identity-changed"}
        )
    filename = "query.json" if command == "panel-history" else "validation.json"
    _write(output / filename, result)
    receipt = {
        "schemaVersion": "primalscheme3.panel-inspection-provenance/v1",
        "tool": "primalscheme3 " + command,
        "toolVersion": version("primalscheme3"),
        "command": {
            "argv": list(argv),
            "shell": shlex.join(argv),
            "workingDirectory": os.getcwd(),
        },
        "requestedOptions": options,
        "resolvedOptions": result.get(
            "query", {"tier": None} | options if command == "panel-audit" else options
        ),
        "inputBundle": str(bundle),
        "inputs": inputs,
        "outputs": [artifact_descriptor(output / filename, output)],
        "source": source,
        "sourceAtEnd": source_end,
        "runtime": runtime,
        "runtimeAtEnd": runtime_end,
        "sourceChangedDuringRun": source["sourceDigest"] != source_end["sourceDigest"],
        "runtimeChangedDuringRun": runtime != runtime_end,
        "inputChangedDuringRun": bool(changed_inputs),
        "wallSeconds": time.time() - started,
        "status": "success" if result["valid"] else "failed",
        "exitStatus": 0 if result["valid"] else 1,
        "stderr": stderr
        or (json.dumps(result.get("violations", [])) if not result["valid"] else ""),
        "selfHashPolicy": "provenance.json excluded to avoid a self-hash cycle",
    }
    _write(output / "provenance.json", receipt)
    return result
