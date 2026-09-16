"""Native allele design workflow; strict publication precedes optional salvage."""

from __future__ import annotations

import gzip
import json
import shutil
from dataclasses import asdict
from pathlib import Path
from time import monotonic
from types import SimpleNamespace

from .allele_labels import SCHEMA as ALLELE_LABEL_SCHEMA
from .allele_labels import build_allele_label_map
from .allele_options import AlleleOptions
from .allele_progress import SCHEMA, ProgressWriter
from .allele_publication import publish_allele_stage
from .allele_search import search_allele_assignments
from .allele_validation import AlleleConstraintProfile, StagePolicy
from .coverage_discovery import (
    build_variant_catalog,
    discovery_profiles,
    variant_targets,
)
from .coverage_history import SQLiteCoverageHistory
from .coverage_pipeline import _export_references
from .coverage_provenance import finalize_provenance
from .coverage_salvage import run_salvage, select_primary
from .coverage_types import ConfigurationLedger


def _write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, default=str) + "\n")


def _authoritative_targets(output_dir, input_records):
    from primalscheme3.core.mapping import create_mapping
    from primalscheme3.core.msa import parse_msa

    sources = {}
    for index, record in enumerate(input_records):
        array, _ = parse_msa(output_dir / record["storedPath"])
        mapping, _ = create_mapping(array)
        sources[index] = SimpleNamespace(
            array=array,
            _mapping_array=mapping,
            msa_index=record.get("sourceIndex", index),
        )
    return variant_targets(sources)


def _gap_report(catalog, result):
    reports = []
    eligible = set(catalog.selectable_site_ids)
    for target in catalog.targets:
        sites = [s for s in catalog.sites if s.target_id == target.id]
        families = [f for f in catalog.families if f.target_id == target.id]
        configs = [c for c in result.ledger.configurations if c.target_id == target.id]
        coverage = next(t for t in result.coverage.targets if t.target_id == target.id)
        causes = []
        if coverage.fraction is None:
            causes.append("unavailable-observations")
        if not any(s.id in eligible for s in sites):
            causes.append("no-individually-eligible-binding-sites")
        elif not families:
            causes.append("no-valid-family-geometry")
        if families and not configs:
            causes.append("no-configurations-evaluated-within-work-bounds")
        if coverage.fraction is not None and coverage.fraction < result.coverage.goal:
            causes.append(
                "inspect-stage-assessments-for-specificity-dimer-overlap-support-and-work-limits"
            )
        reports.append(
            {
                "target_id": target.id,
                "coverage_fraction": coverage.fraction,
                "goal": result.coverage.goal,
                "generated_sites": len(sites),
                "eligible_sites": sum(s.id in eligible for s in sites),
                "families": len(families),
                "explored_configurations": len(configs),
                "diagnostic_categories": causes,
            }
        )
    return {
        "schemaVersion": "primalscheme3.allele-gap-diagnostics/v1",
        "targets": reports,
        "scope": "Catalog and evaluated-configuration diagnostics, not a feasible-union coverage bound or exhaustive causal attribution.",
        "search_work": result.metadata.get("proposal_work", {}),
    }


def run_allele_pipeline(
    *,
    msa_dict,
    msa_data,
    input_records,
    output_dir,
    config,
    config_dict,
    max_amplicons,
    max_amplicons_msa,
    logger,
    argv,
    started_at,
    execution_start,
    invocation_state=None,
):
    output_dir = Path(output_dir)
    history = None
    progress = None
    timings = {}
    origin_discovery = None
    scientific = {"selectionAlgorithm": "allele-coverage"}
    options = AlleleOptions.from_config(config)
    try:
        history = SQLiteCoverageHistory(output_dir / "history", run_id="allele-panel")
        start = monotonic()
        profiles = discovery_profiles(config)
        if options.candidate_profiles != "union":
            profiles = {
                options.candidate_profiles: profiles[options.candidate_profiles]
            }
        if options.reuse_discovery:
            from .allele_catalog_cache import (
                load_discovery_cache,
                materialize_cache_reuse,
            )

            logger.info(
                "Validating and reusing discovery catalog and immutable origin history"
            )
            reuse = load_discovery_cache(
                options.reuse_discovery,
                targets=variant_targets(msa_dict),
                config=config,
                profiles=profiles,
            )
            materialize_cache_reuse(reuse, output_dir, history=history)
            catalog = reuse.catalog
            origin_discovery = {"path": "origin-discovery/manifest.json"}
            config.discovery_workers_by_msa = {
                str(t.source_msa_index): 0 for t in catalog.targets
            }
            config.discovery_workers_by_target_profile = {
                t.id: {p: 0 for p in profiles} for t in catalog.targets
            }
            config_dict["discovery_reused"] = True
        else:
            logger.info(
                "Discovering concrete primer variants with complete profile membership"
            )
            catalog = build_variant_catalog(
                msa_dict,
                config,
                profiles=profiles,
                history=history,
                length_mode=options.discovery_length_mode,
            )
            config_dict["discovery_reused"] = False
        timings["discovery_seconds"] = monotonic() - start
        # Discovery is expensive and scientifically useful even when a later
        # selector or publication fails. Keep its complete catalog independently.
        with gzip.open(output_dir / "discovery-catalog.json.gz", "wt") as handle:
            json.dump(catalog.to_dict(), handle, sort_keys=True, separators=(",", ":"))
        _write(
            output_dir / "allele-label-map.json",
            build_allele_label_map(catalog, output_dir, input_records),
        )
        history.checkpoint()
        for name in (
            "discovery_workers_by_msa",
            "discovery_workers_by_target_profile",
            "discovery_core_count",
            "discovery_backend",
        ):
            if hasattr(config, name):
                config_dict[name] = getattr(config, name)
        if origin_discovery is not None:
            config_dict["discovery_core_count"] = 0
        authoritative = _authoritative_targets(output_dir, input_records)
        references = _export_references(catalog, msa_dict)
        profile = AlleleConstraintProfile.from_config(
            config,
            max_amplicons=max_amplicons,
            max_amplicons_msa=max_amplicons_msa,
            secondary_product_policy=options.secondary_product_policy,
        )
        search_options = options.search_options()
        scores = {}
        start = monotonic()
        progress = ProgressWriter(output_dir / "search-progress.jsonl")
        scientific["searchProgress"] = {
            "path": "search-progress.jsonl",
            "schemaVersion": SCHEMA,
            "validationScope": "provisional-oracle-checked",
        }
        strict = search_allele_assignments(
            catalog,
            profile,
            search_options,
            history=history,
            score_cache=scores,
            observer=progress,
        )
        timings["strict_search_and_validation_seconds"] = monotonic() - start
        publications = {}
        stage_results = {"strict": strict}

        def publish(result):
            stage = result.metadata["stage_policy"]["stage_id"]
            begin = monotonic()
            publications[stage] = publish_allele_stage(
                output_dir / "stages" / stage,
                catalog=catalog,
                ledger=result.ledger,
                assignments=result.assignments,
                profile=profile,
                policy=StagePolicy(**result.metadata["stage_policy"]),
                authoritative_targets=authoritative,
                references=references,
                goal=options.coverage_target,
                history=history,
            )
            timings[stage + "_publication_and_audit_seconds"] = monotonic() - begin
            stage_results[stage] = result

        publish(
            strict
        )  # Durable independently audited strict artifact before any relaxation.
        start = monotonic()
        salvage = run_salvage(
            strict,
            profile,
            options.salvage_options(),
            search_options=search_options,
            history=history,
            on_tier=publish,
            score_cache=scores,
            observer=progress,
        )
        timings["salvage_total_seconds"] = monotonic() - start
        # Failed publication still retains every explored configuration for history queries.
        explored = dict(strict.ledger.configuration_by_id)
        for tier in salvage.tiers:
            if tier.result is not None:
                explored.update(tier.result.ledger.configuration_by_id)
        complete_ledger = ConfigurationLedger(
            catalog.semantic_digest, tuple(explored.values())
        )
        with gzip.open(output_dir / "configuration-ledger.json.gz", "wt") as handle:
            json.dump(
                complete_ledger.to_dict(), handle, sort_keys=True, separators=(",", ":")
            )
        _write(
            output_dir / "salvage-status.json",
            {
                "stop_reason": salvage.stop_reason,
                "tiers": [
                    {
                        "stage_id": tier.stage_id,
                        "status": tier.status,
                        "error": tier.error,
                        "explored_ledger_digest": tier.result.ledger.semantic_digest
                        if tier.result is not None
                        else None,
                    }
                    for tier in salvage.tiers
                ],
            },
        )
        primary = select_primary(salvage, options.primary_tier)
        primary_stage = primary.metadata["stage_policy"]["stage_id"]
        primary_dir = output_dir / "stages" / primary_stage
        primary_files = (
            "primer.bed",
            "amplicon.bed",
            "primertrim.amplicon.bed",
            "primers.fasta",
            "order-sheet.tsv",
            "reference.fasta",
        )
        for name in primary_files:
            shutil.copyfile(primary_dir / name, output_dir / name)
        _write(
            output_dir / "panel-validation.json",
            publications[primary_stage]["validation"],
        )
        _write(output_dir / "row-coverage.json", primary.coverage.to_dict())
        _write(output_dir / "gap-diagnostics.json", _gap_report(catalog, primary))
        _write(
            output_dir / "salvage-comparison.json",
            {
                "options": asdict(salvage.options),
                "strictPreserved": True,
                "stopReason": salvage.stop_reason,
                "tiers": [
                    {
                        "stage_id": t.stage_id,
                        "status": t.status,
                        "policy": asdict(t.policy),
                        "comparison": t.comparison,
                        "error": t.error,
                        "elapsed_seconds": t.elapsed_seconds,
                    }
                    for t in salvage.tiers
                ],
            },
        )
        previous = history.last_complete_stage
        ids = tuple(a.configuration_id for a in primary.assignments)
        history.emit(
            stage_id="publication",
            kind="published-primary",
            entity_ids=ids,
            changes={
                "primary_tier": primary_stage,
                "relative_stage_path": str(Path("stages") / primary_stage),
            },
        )
        touched = {
            entity
            for event in history.iter_events(
                start=previous.events_count if previous else 0
            )
            for entity in event.entity_ids
        }
        history.complete_stage(
            stage_id="publication",
            dispositions={
                e: "published-primary"
                if e in ids
                else "recorded-publication-diagnostic"
                for e in touched
            },
            catalog_digest=catalog.semantic_digest,
            ledger_digest=complete_ledger.semantic_digest,
            prior_snapshot_id=previous.id if previous else None,
        )
        history_counts = {
            name: len(getattr(history, name))
            for name in ("evidence", "assessments", "events", "snapshots")
        }
        history.close()
        history = None
        stages = [
            {
                "stage_id": stage,
                "path": "stages/" + stage,
                "coverage": result.coverage.to_dict(),
                "assignments": [a.to_dict() for a in result.assignments],
                "optimizer": result.metadata,
                "validationValid": publications[stage]["validation"]["valid"],
            }
            for stage, result in stage_results.items()
        ]
        optimizer = {
            "schemaVersion": "primalscheme3.panel-optimizer/v2",
            "algorithm": "bounded-allele-coverage/v1",
            "metric": "observed-allele-primer-trimmed/v1",
            "primaryTier": primary_stage,
            "stages": stages,
            "profile": profile.to_dict(),
            "options": options.to_dict(),
            "timings": timings,
            "history": {
                "path": "history/history.sqlite",
                "counts": history_counts,
                "configurationLedger": "configuration-ledger.json.gz",
            },
            "publication": {
                "targetToReference": references,
                "alleleLabelMap": {
                    "path": "allele-label-map.json",
                    "schemaVersion": ALLELE_LABEL_SCHEMA,
                },
                "configurationToAmplicon": publications[primary_stage]["manifest"][
                    "configuration_names"
                ],
                "artifacts": {name: name for name in primary_files},
                "projection": "ungapped-first-reference-coordinates",
            },
        }
        if origin_discovery is not None:
            optimizer["history"]["originDiscovery"] = origin_discovery
        _write(output_dir / "panel-optimizer.json", optimizer)
        config_dict.update(
            output=".",
            msa_data=msa_data,
            panel_optimizer={
                "schemaVersion": "primalscheme3.panel-native-integration/v2",
                "selectionAlgorithm": "allele-coverage",
                "toolVersion": config.version,
                "profile": profile.to_dict(),
                "options": options.to_dict(),
                "primaryTier": primary_stage,
                "optimizer": {
                    "path": "panel-optimizer.json",
                    "schemaVersion": optimizer["schemaVersion"],
                },
                "validation": {
                    "path": "panel-validation.json",
                    "schemaVersion": primary.validation["schemaVersion"],
                    "valid": True,
                },
                "provenance": {
                    "path": "panel-provenance.json",
                    "schemaVersion": "primalscheme3.panel-provenance/v1",
                },
                "publication": optimizer["publication"],
            },
        )
        _write(output_dir / "config.json", config_dict)
        scientific.update(
            algorithm=optimizer["algorithm"],
            metric=optimizer["metric"],
            primaryTier=primary_stage,
            profile=profile.to_dict(),
            resolvedAlleleOptions=options.to_dict(),
            catalogSemanticDigest=catalog.semantic_digest,
            ledgerSemanticDigest=primary.ledger.semantic_digest,
            validationValid=True,
            timings=timings,
            history=optimizer["history"],
        )
        logger.info(
            "Completed allele-aware panel with independently validated primary tier %s",
            primary_stage,
        )
        if progress is not None:
            progress.close()
        finalize_provenance(
            output_dir=output_dir,
            argv=argv,
            resolved_options=config_dict,
            inputs=input_records,
            started_at=started_at,
            ended_at=monotonic(),
            status="success",
            exit_status=0,
            stderr="",
            scientific=scientific,
            logger=logger,
            execution_start=execution_start,
        )
        if invocation_state is not None:
            invocation_state.provenance_finalized = True
        return optimizer
    except BaseException as error:
        for name, resource in (("progress", progress), ("history", history)):
            if resource is not None:
                try:
                    resource.close()
                except BaseException as cleanup_error:
                    error.add_note(f"{name} close failed: {cleanup_error}")
        if getattr(error, "__notes__", None):
            scientific["secondaryErrors"] = list(error.__notes__)
        logger.exception("Allele-aware panel failed")
        finalize_provenance(
            output_dir=output_dir,
            argv=argv,
            resolved_options=config_dict,
            inputs=input_records,
            started_at=started_at,
            ended_at=monotonic(),
            status="failure",
            exit_status=1,
            stderr=str(error),
            scientific=scientific,
            logger=logger,
            execution_start=execution_start,
        )
        if invocation_state is not None:
            invocation_state.provenance_finalized = True
        raise
