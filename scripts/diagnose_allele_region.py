#!/usr/bin/env python3
"""Diagnose a reference BED region of a completed native allele-aware panel.

Read-only inputs. Alternative insertions are separate experiments against the
saved panel, never a combined feasible panel or coverage ceiling.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shlex
import sys
from collections import defaultdict
from dataclasses import asdict
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from statistics import mean
from time import monotonic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from primalscheme3.panel.allele_coverage import (  # noqa: E402
    _intervals,
    binding_support,
    configuration_coverage,
)
from primalscheme3.panel.allele_inspection import (  # noqa: E402
    _contained,
    audit_allele_bundle,
    query_allele_history,
)
from primalscheme3.panel.allele_publication import (  # noqa: E402
    _read,
    _targets,
    artifact_descriptor,
)
from primalscheme3.panel.allele_validation import (  # noqa: E402
    AlleleConstraintProfile,
    StagePolicy,
    validate_allele_assignments,
)
from primalscheme3.panel.coverage_history import SQLiteCoverageHistory  # noqa: E402
from primalscheme3.panel.coverage_provenance import (  # noqa: E402
    runtime_identity,
    source_identity,
)
from primalscheme3.panel.coverage_types import (  # noqa: E402
    AlleleAssignment,
    ConfigurationLedger,
    SelectionConfiguration,
    VariantCatalog,
)
from primalscheme3.panel.coverage_variants import make_configuration  # noqa: E402

TOOL = "diagnose-allele-region"
SCOPE = (
    "Separate insertions into the unchanged saved panel, not optimization, a "
    "jointly feasible alternative panel, or a feasible coverage ceiling. Supplied "
    "MSA rows only; no whole-genome specificity claim. Historical acceptance and "
    "dimer values are context only; historical specificity D=0 does not satisfy "
    "current selected-sites-v2 specificity."
)


def _write(path, value):
    text = json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def region_coverage(catalog, ledger, assignments, target, region):
    """Project BED boundaries to alignment, then count concrete row positions."""
    start, end = region
    if not 0 <= start < end <= target.reference_length:
        raise ValueError("region must lie within the ungapped first reference")
    left, right = target.ref_to_alignment[start], target.ref_to_alignment[end]
    configs = {
        a.configuration_id: ledger.configuration_by_id[a.configuration_id]
        for a in assignments
        if ledger.configuration_by_id[a.configuration_id].target_id == target.id
        and ledger.configuration_by_id[a.configuration_id].full_interval[0] < end
        and start < ledger.configuration_by_id[a.configuration_id].full_interval[1]
    }
    coverages = {cid: configuration_coverage(catalog, c) for cid, c in configs.items()}
    site_ids = sorted(
        {s for c in configs.values() for s in c.forward_site_ids + c.reverse_site_ids}
    )
    classes = []
    for allele in catalog.observations:
        if allele.target_id != target.id:
            continue
        observed = {
            p
            for p in allele.observed_positions
            if left <= allele.row_to_alignment[p] < right
        }
        covered = (
            set().union(*(set(v.get(allele.id, ())) for v in coverages.values()))
            & observed
        )
        contributions = {
            cid: _intervals(set(v.get(allele.id, ())) & observed)
            for cid, v in coverages.items()
        }
        classes.append(
            {
                "allele_id": allele.id,
                "row_ids": list(allele.row_ids),
                "multiplicity": allele.multiplicity,
                "observed_bases": len(observed),
                "covered_bases": len(covered),
                "fraction": len(covered) / len(observed) if observed else None,
                "observed_row_intervals": _intervals(observed),
                "covered_row_intervals": _intervals(covered),
                "uncovered_row_intervals": _intervals(observed - covered),
                "unknown_alignment_columns": [
                    i
                    for i in range(left, right)
                    if allele.cells[i] not in ("", "-", "A", "C", "G", "T")
                ],
                "unavailable_alignment_columns": [
                    i for i in range(left, right) if allele.cells[i] == ""
                ],
                "configuration_contributions": contributions,
                "site_support": [
                    binding_support(catalog.site_by_id[s], allele).to_dict()
                    for s in site_ids
                ],
            }
        )
    fractions = [r["fraction"] for r in classes if r["fraction"] is not None]
    return {
        "target_id": target.id,
        "reference_region": list(region),
        "alignment_region": [left, right],
        "projection": "reference inverse-mapping half-open boundaries; includes row insertions inside alignment span",
        "weighting": "equal distinct observed classes; row aliases do not multiply weight",
        "mean_class_fraction": mean(fractions) if fractions else None,
        "classes": classes,
        "configurations": [c.to_dict() for c in configs.values()],
        "assignments": [
            a.to_dict() for a in assignments if a.configuration_id in configs
        ],
        "sites": [catalog.site_by_id[s].to_dict() for s in site_ids],
    }


def _site_index(catalog, target):
    index = defaultdict(list)
    for site in catalog.sites:
        if site.target_id == target.id and site.reference_footprint is not None:
            anchor = site.reference_footprint[1 if site.strand == "+" else 0]
            index[(site.sequence, site.strand, anchor)].append(site)
    return index


def match_alternative(catalog, target, item, *, index=None):
    """All original variants must match; never silently delete failed variants."""
    index = _site_index(catalog, target) if index is None else index
    sides, issues = [], []
    for name, strand in (("forward", "+"), ("reverse", "-")):
        cloud = item.get(name, {})
        span, sequences = cloud.get("regionBED"), cloud.get("sequences")
        if (
            not isinstance(span, list)
            or len(span) != 2
            or any(type(v) is not int for v in span)
            or not 0 <= span[0] < span[1] <= target.reference_length
            or not isinstance(sequences, list)
            or not sequences
            or any(
                not isinstance(s, str) or not s or set(s) - set("ACGT")
                for s in sequences
            )
        ):
            issues.append({"side": name, "reason": "invalid-legacy-cloud"})
            sides.append(())
            continue
        anchor = span[1 if strand == "+" else 0]
        matched = []
        for seq in sorted(set(sequences)):
            sites = index.get((seq, strand, anchor), ())
            reason = (
                "no-exact-site"
                if not sites
                else "ambiguous-exact-site"
                if len(sites) != 1
                else "site-not-eligible"
                if not sites[0].accepting_profile_ids or sites[0].mapping_failure
                else None
            )
            if reason:
                issues.append(
                    {
                        "side": name,
                        "sequence": seq,
                        "reference_anchor": anchor,
                        "reason": reason,
                        "site_ids": [s.id for s in sites],
                    }
                )
            else:
                matched.append(sites[0])
        sides.append(tuple(matched))
    if issues:
        return {"status": "unavailable", "issues": issues}
    fs, rs = sides
    anchors = ({s.alignment_anchor for s in fs}, {s.alignment_anchor for s in rs})
    families = [
        f
        for f in catalog.families
        if f.target_id == target.id
        and anchors == ({f.anchor_pair[0]}, {f.anchor_pair[1]})
    ]
    if len(families) != 1:
        return {
            "status": "unavailable",
            "issues": [{"reason": "no-unique-existing-family"}],
        }
    try:
        c = make_configuration(
            catalog, families[0].id, tuple(s.id for s in fs), tuple(s.id for s in rs)
        )
    except ValueError as error:
        return {
            "status": "unavailable",
            "issues": [{"reason": "configuration-ineligible", "detail": str(error)}],
        }
    return {"status": "matched", "configuration": c, "issues": []}


def evaluate_alternative(
    catalog,
    ledger,
    assignments,
    candidate,
    profile,
    policy,
    goal,
    history,
    *,
    authoritative_targets=None,
):
    configs = dict(ledger.configuration_by_id)
    configs[candidate.id] = candidate
    ledger = ConfigurationLedger(catalog.semantic_digest, tuple(configs.values()))
    kwargs = dict(
        policy=policy,
        goal=goal,
        history=history,
        authoritative_targets=authoritative_targets,
    )
    unary = validate_allele_assignments(
        catalog,
        (AlleleAssignment(candidate.id, 0, policy.stage_id),),
        ledger,
        profile,
        **kwargs,
    )
    pools = []
    for pool in range(profile.n_pools):
        inserted = AlleleAssignment(candidate.id, pool, policy.stage_id)
        if any(
            a.configuration_id == candidate.id and a.pool == pool for a in assignments
        ):
            pools.append({"pool_zero_based": pool, "status": "already-selected"})
            continue
        result = validate_allele_assignments(
            catalog, (*assignments, inserted), ledger, profile, **kwargs
        )
        pools.append(
            {"pool_zero_based": pool, "status": "evaluated", "validation": result}
        )
    return {
        "configuration": candidate.to_dict(),
        "unary": unary,
        "pools": pools,
        "scope": SCOPE,
    }


def _identity():
    return {
        "native": source_identity(),
        "script": artifact_descriptor(Path(__file__), ROOT),
    }


def _source_changed(before, after):
    return (
        before["native"]["sourceDigest"] != after["native"]["sourceDigest"]
        or before["script"] != after["script"]
    )


def _legacy_records(paths, catalog, target, region):
    index = _site_index(catalog, target)
    records = []
    for source_index, path in enumerate(paths):
        data = _read(path)
        candidates = data.get("candidates")
        if not isinstance(candidates, list):
            raise ValueError("legacy alternatives require a candidates list")
        for ordinal, item in enumerate(candidates):
            if not isinstance(item, dict):
                raise ValueError("legacy candidate must be an object")
            match = match_alternative(catalog, target, item, index=index)
            c = match.pop("configuration", None)
            record = {
                "source_path": str(path),
                "source_index": source_index,
                "ordinal": ordinal,
                "legacy_candidate": item,
                **match,
            }
            if c is not None:
                f_end = max(
                    catalog.site_by_id[s].reference_footprint[1]
                    for s in c.forward_site_ids
                )
                r_start = min(
                    catalog.site_by_id[s].reference_footprint[0]
                    for s in c.reverse_site_ids
                )
                record.update(
                    configuration_id=c.id,
                    configuration=c.to_dict(),
                    primer_trimmed_reference_overlap=max(
                        0, min(region[1], r_start) - max(region[0], f_end)
                    ),
                )
            records.append(record)
    return records


def select_alternatives(records, budget):
    matched = [r for r in records if r["status"] == "matched"]

    def key(record):
        index = record["legacy_candidate"].get("catalogIndex")
        return (
            -record["primer_trimmed_reference_overlap"],
            record["source_index"],
            index if type(index) is int else record["ordinal"],
            record["ordinal"],
        )

    ranked = sorted(matched, key=key)
    return list(dict.fromkeys(r["configuration_id"] for r in ranked))[:budget]


def diagnose(
    bundle,
    output,
    *,
    target_index=None,
    target_id=None,
    region,
    legacy_alternatives=(),
    max_candidates=16,
    history_limit=100,
):
    if type(max_candidates) is not int or not 0 <= max_candidates <= 1000:
        raise ValueError("max_candidates must be between 0 and 1000")
    if type(history_limit) is not int or not 1 <= history_limit <= 1000:
        raise ValueError("history_limit must be between 1 and 1000")
    if (target_index is None) == (target_id is None):
        raise ValueError("exactly one target index or ID is required")
    audit = audit_allele_bundle(bundle)
    _write(output / "input-bundle-audit.json", audit)
    if not audit["valid"]:
        raise ValueError("input bundle failed fresh scientific audit")
    optimizer = _read(bundle / "panel-optimizer.json")
    stage = next(
        s for s in optimizer["stages"] if s["stage_id"] == optimizer["primaryTier"]
    )
    stage_path = _contained(bundle, stage["path"])
    manifest = _read(stage_path / "stage.json")
    catalog = VariantCatalog.from_dict(_read(stage_path / "catalog.json.gz"))
    ledger = ConfigurationLedger.from_dict(
        _read(_contained(bundle, optimizer["history"]["configurationLedger"]))
    )
    if ledger.catalog_digest != catalog.semantic_digest:
        raise ValueError("complete ledger catalog mismatch")
    targets = [
        t
        for t in catalog.targets
        if t.id == target_id
        or (target_id is None and t.source_msa_index == target_index)
    ]
    if len(targets) != 1:
        raise ValueError("target must resolve to one exact occurrence")
    target = targets[0]
    assignments = tuple(
        AlleleAssignment.from_dict(v) for v in _read(stage_path / "assignments.json")
    )
    selected = region_coverage(catalog, ledger, assignments, target, region)
    nearby = [
        c
        for c in ledger.configurations
        if c.target_id == target.id
        and c.full_interval[0] < region[1]
        and region[0] < c.full_interval[1]
    ]
    history_query = query_allele_history(
        bundle, target=target.id, region=region, limit=history_limit
    )
    _write(output / "nearby-history.json.gz", history_query)
    records = _legacy_records(legacy_alternatives, catalog, target, region)
    matched = [r for r in records if r["status"] == "matched"]
    chosen = select_alternatives(records, max_candidates)
    profile = AlleleConstraintProfile(**manifest["constraints"])
    policy = StagePolicy(**manifest["stage_policy"])
    targets = _targets(_read(stage_path / "authoritative-targets.json.gz"))
    diagnostic_configs = dict(ledger.configuration_by_id)
    diagnostic_configs.update(
        {
            r["configuration_id"]: SelectionConfiguration.from_dict(r["configuration"])
            for r in matched
        }
    )
    diagnostic_ledger = ConfigurationLedger(
        catalog.semantic_digest, tuple(diagnostic_configs.values())
    )
    _write(output / "diagnostic-ledger.json.gz", diagnostic_ledger.to_dict())
    results = {}
    (output / "alternatives").mkdir()
    with SQLiteCoverageHistory(output / "history", run_id=TOOL) as history:
        for number, cid in enumerate(chosen):
            candidate = next(c for c in matched if c["configuration_id"] == cid)
            config = SelectionConfiguration.from_dict(candidate["configuration"])
            result = evaluate_alternative(
                catalog,
                ledger,
                assignments,
                config,
                profile,
                policy,
                manifest["goal"],
                history,
                authoritative_targets=targets,
            )
            # Coverage here is binding-model attribution, including rejected insertions.
            # The adjacent fresh validation is authoritative for feasibility.
            result["candidate_region_binding_coverage"] = region_coverage(
                catalog,
                diagnostic_ledger,
                (AlleleAssignment(cid, 0, policy.stage_id),),
                target,
                region,
            )
            result["saved_panel_plus_candidate_region_binding_coverage"] = (
                region_coverage(
                    catalog,
                    diagnostic_ledger,
                    (*assignments, AlleleAssignment(cid, 0, policy.stage_id)),
                    target,
                    region,
                )
            )
            result["region_coverage_scope"] = (
                "Binding-model attribution only; a failed insertion has no feasible-panel coverage claim."
            )
            filename = f"alternatives/{number:04d}.json.gz"
            _write(output / filename, result)
            results[cid] = {
                "path": filename,
                "unary_valid": result["unary"]["valid"],
                "pool_results": [
                    {
                        "pool_zero_based": p["pool_zero_based"],
                        "status": p["status"],
                        "valid": p.get("validation", {}).get("valid"),
                    }
                    for p in result["pools"]
                ],
            }
            history.emit(
                stage_id="region-diagnostic",
                kind="alternative-evaluated",
                entity_ids=(cid,),
                changes={"result_path": filename, "saved_panel_stage": policy.stage_id},
            )
        history.complete_stage(
            stage_id="region-diagnostic",
            dispositions={c: "evaluated-independent-insertion" for c in chosen},
            catalog_digest=catalog.semantic_digest,
            ledger_digest=diagnostic_ledger.semantic_digest,
        )
    for record in matched:
        record["insertion"] = results.get(
            record["configuration_id"], {"status": "not-evaluated-budget"}
        )
    return {
        "schemaVersion": "primalscheme3.allele-region-diagnostic/v1",
        "valid": True,
        "scope": SCOPE,
        "selected_stage": policy.stage_id,
        "stage_policy": asdict(policy),
        "profile": profile.to_dict(),
        "selected": selected,
        "nearby_explored": {
            "count": len(nearby),
            "configuration_ids": [c.id for c in nearby],
            "history_path": "nearby-history.json.gz",
            "history_limit": history_limit,
            "scope": "Exact complete-ledger interval count; reasons are bounded contextual history results, not exhaustive causal attribution.",
        },
        "legacy_sources": [
            {
                "path": str(path),
                "metadata": {k: v for k, v in _read(path).items() if k != "candidates"},
            }
            for path in legacy_alternatives
        ],
        "alternatives": records,
        "work": {
            "total_legacy_records": len(records),
            "matched_legacy_records": len(matched),
            "unavailable_legacy_records": len(records) - len(matched),
            "eligible_unique_configurations": len(
                {r["configuration_id"] for r in matched}
            ),
            "examined_unique_configurations": len(chosen),
            "omitted_unique_configurations": len(
                {r["configuration_id"] for r in matched}
            )
            - len(chosen),
            "max_candidates": max_candidates,
            "ranking": "descending primer-trimmed reference overlap; source argument order; numeric catalog index (ordinal fallback); original ordinal",
            "exhaustive": len({r["configuration_id"] for r in matched}) == len(chosen),
        },
    }


def run(bundle, output, *, argv, **options):
    bundle, output = Path(bundle).resolve(), Path(output).resolve()
    alternatives = tuple(
        Path(p).resolve() for p in options.get("legacy_alternatives", ())
    )
    if (
        output.exists()
        or output == bundle
        or output.is_relative_to(bundle)
        or any(p.is_relative_to(output) for p in alternatives)
    ):
        raise ValueError("output must be fresh and outside all source inputs")
    output.mkdir(parents=True)
    started = monotonic()
    started_at = datetime.now(UTC).isoformat()
    source, runtime = _identity(), runtime_identity()
    inputs, stderr = [], ""
    try:
        # Reject active SQLite before reading scientific input or hashing a live DB.
        if any(bundle.rglob("*-wal")) or any(bundle.rglob("*-journal")):
            raise ValueError("source bundle contains an active SQLite transaction")
        for path in [bundle / "panel-provenance.json", *alternatives]:
            if path.is_file():
                inputs.append(
                    artifact_descriptor(path, path.parent) | {"path": str(path)}
                )
        receipt = _read(bundle / "panel-provenance.json")
        if receipt.get("status") != "success" or receipt.get("exitStatus") != 0:
            raise ValueError("source must be a completed successful native panel")
        paths = sorted(
            {p.resolve() for p in bundle.rglob("*") if p.is_file()} | set(alternatives)
        )
        inputs = []
        for path in paths:
            inputs.append(artifact_descriptor(path, path.parent) | {"path": str(path)})
        if isinstance(options.get("region"), str):
            options["region"] = tuple(int(v) for v in options["region"].split(":"))
        result = diagnose(
            bundle, output, **(options | {"legacy_alternatives": alternatives})
        )
    except (Exception, KeyboardInterrupt) as error:
        stderr = f"{type(error).__name__}: {error}"
        result = {"valid": False, "error": stderr, "scope": SCOPE}
    source_end, runtime_end = _identity(), runtime_identity()
    changed = []
    for item in inputs:
        path = Path(item["path"])
        if (
            not path.is_file()
            or (artifact_descriptor(path, path.parent) | {"path": str(path)}) != item
        ):
            changed.append(item["path"])
    if changed or _source_changed(source, source_end) or runtime != runtime_end:
        result["valid"] = False
        stderr = stderr or "input or execution identity changed during diagnosis"
        result["error"] = stderr
    _write(output / "diagnostic.json", result)
    outputs = [
        artifact_descriptor(p, output) for p in sorted(output.rglob("*")) if p.is_file()
    ]
    _write(
        output / "provenance.json",
        {
            "schemaVersion": "primalscheme3.scientific-command-provenance/v1",
            "tool": TOOL,
            "toolVersion": "1",
            "nativeVersion": version("primalscheme3"),
            "command": {
                "argv": list(argv),
                "shell": shlex.join(argv),
                "workingDirectory": os.getcwd(),
            },
            "resolvedOptions": {"max_candidates": 16, "history_limit": 100}
            | {
                k: [str(x) for x in v] if k == "legacy_alternatives" else v
                for k, v in options.items()
            },
            "inputBundle": str(bundle),
            "outputDirectory": str(output),
            "inputs": inputs,
            "outputs": outputs,
            "source": source,
            "sourceAtEnd": source_end,
            "runtime": runtime,
            "runtimeAtEnd": runtime_end,
            "sourceChangedDuringRun": _source_changed(source, source_end),
            "runtimeChangedDuringRun": runtime != runtime_end,
            "changedInputPaths": changed,
            "wallSeconds": monotonic() - started,
            "startedAt": started_at,
            "finishedAt": datetime.now(UTC).isoformat(),
            "status": "success" if result["valid"] else "failed",
            "exitStatus": 0 if result["valid"] else 1,
            "stderr": stderr,
            "selfHashPolicy": "provenance.json excluded to avoid a self-hash cycle",
        },
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument(
        "--target-index", type=int, help="Zero-based input occurrence index"
    )
    target.add_argument("--target-id")
    parser.add_argument(
        "--region", required=True, help="Zero-based half-open reference BED start:end"
    )
    parser.add_argument("--legacy-alternatives", action="append", default=[], type=Path)
    parser.add_argument("--max-candidates", type=int, default=16)
    parser.add_argument("--history-limit", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    args = vars(parser.parse_args())
    try:
        result = run(**args, argv=[sys.executable, *sys.argv])
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2
    if not result["valid"]:
        print(result["error"], file=sys.stderr)
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
