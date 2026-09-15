"""Portable selected-stage artifacts, checked from their saved bytes before publish."""

from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
from dataclasses import asdict
from pathlib import Path

from .allele_validation import (
    AlleleConstraintProfile,
    StagePolicy,
    validate_allele_assignments,
)
from .coverage_types import (
    AlleleAssignment,
    ConfigurationLedger,
    Target,
    VariantCatalog,
)


def artifact_descriptor(path, root):
    path, root = Path(path), Path(root)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": digest.hexdigest(),
        "size": path.stat().st_size,
    }


def _write(path, value):
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            json.dump(value, handle, sort_keys=True, separators=(",", ":"))
    else:
        path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def _read(path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return json.load(handle)
    return json.loads(path.read_text())


def _targets(values):
    return tuple(
        Target(
            **{
                **v,
                "row_ids": tuple(v["row_ids"]),
                "rows": tuple(tuple(r) for r in v["rows"]),
                "mapping": tuple(v["mapping"]),
                "ref_to_alignment": tuple(v["ref_to_alignment"]),
            }
        )
        for v in values
    )


def _selected(catalog, ledger, assignments):
    records = []
    for assignment in sorted(assignments, key=lambda a: (a.configuration_id, a.pool)):
        config = ledger.configuration_by_id[assignment.configuration_id]
        for sid in config.forward_site_ids + config.reverse_site_ids:
            site = catalog.site_by_id[sid]
            records.append(
                {
                    "configuration_id": config.id,
                    "pool": assignment.pool,
                    "site_id": site.id,
                    "sequence": site.sequence,
                    "strand": site.strand,
                    "alignment_anchor": site.alignment_anchor,
                    "reference_footprint": list(site.reference_footprint),
                }
            )
    return records


def _render(catalog, ledger, assignments, references, stage_id):
    bed = ["# artic-bed-version v3.0", "# tier=" + stage_id]
    fasta = []
    orders = ["name\tsequence\tpool\ttier\tsite_id\tconfiguration_id"]
    amplicons = []
    trimmed = []
    names = {}
    for index, a in enumerate(
        sorted(
            assignments,
            key=lambda a: (
                ledger.configuration_by_id[a.configuration_id].target_id,
                ledger.configuration_by_id[a.configuration_id].full_interval,
                a.configuration_id,
                a.pool,
            ),
        ),
        1,
    ):
        c = ledger.configuration_by_id[a.configuration_id]
        ref = references[c.target_id]
        name = f"allele_{index}"
        names[c.id] = name
        for side, ids in (("LEFT", c.forward_site_ids), ("RIGHT", c.reverse_site_ids)):
            for variant, sid in enumerate(sorted(ids), 1):
                s = catalog.site_by_id[sid]
                start, end = s.reference_footprint
                n = f"{name}_{side}_{variant}"
                attrs = f"site_id={sid};configuration_id={c.id};anchor={s.alignment_anchor};tier={stage_id}"
                bed.append(
                    "\t".join(
                        map(
                            str,
                            (
                                ref,
                                start,
                                end,
                                n,
                                a.pool + 1,
                                s.strand,
                                s.sequence,
                                attrs,
                            ),
                        )
                    )
                )
                fasta.extend((">" + n, s.sequence))
                orders.append(
                    "\t".join(
                        map(str, (n, s.sequence, a.pool + 1, stage_id, sid, c.id))
                    )
                )
        amplicons.append("\t".join(map(str, (ref, *c.full_interval, name, a.pool + 1))))
        f_end = max(
            catalog.site_by_id[s].reference_footprint[1] for s in c.forward_site_ids
        )
        r_start = min(
            catalog.site_by_id[s].reference_footprint[0] for s in c.reverse_site_ids
        )
        trimmed.append("\t".join(map(str, (ref, f_end, r_start, name, a.pool + 1))))
    texts = {
        "primer.bed": "\n".join(bed) + "\n",
        "primers.fasta": "\n".join(fasta) + ("\n" if fasta else ""),
        "order-sheet.tsv": "\n".join(orders) + "\n",
        "amplicon.bed": "\n".join(amplicons) + ("\n" if amplicons else ""),
        "primertrim.amplicon.bed": "\n".join(trimmed) + ("\n" if trimmed else ""),
        "reference.fasta": "".join(
            ">" + references[t.id] + "\n" + t.reference_sequence + "\n"
            for t in catalog.targets
        ),
    }
    return texts, names


def _bed_records(path, manifest, catalog, ledger):
    records = []
    violations = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) != 8:
            raise ValueError("invalid selected BED columns")
        attrs = dict(x.split("=", 1) for x in c[7].split(";"))
        cid = attrs["configuration_id"]
        configuration = ledger.configuration_by_id.get(cid)
        if (
            configuration is None
            or c[0] != manifest["references"].get(configuration.target_id)
            or attrs.get("tier") != manifest["stage_id"]
        ):
            violations.append({"reason": "selected-bed-reference-or-tier-mismatch"})
        records.append(
            {
                "configuration_id": cid,
                "pool": int(c[4]) - 1,
                "site_id": attrs["site_id"],
                "sequence": c[6],
                "strand": c[5],
                "alignment_anchor": int(attrs["anchor"]),
                "reference_footprint": [int(c[1]), int(c[2])],
            }
        )
    return records, violations


def _literal_coverage(targets, records, catalog):
    """Second implementation: walk saved row cells, literal full oligo equality.

    Does not call binding_support, configuration_products or allele_summary and
    does not read stored products/coverage caches. No specificity hits give credit.
    """
    complement = str.maketrans("ACGT", "TGCA")
    by_config = {}
    for record in records:
        by_config.setdefault(record["configuration_id"], []).append(record)
    result = []
    for target in targets:
        rows = {tuple(row): [] for row in target.rows}
        for row_id, row in zip(target.row_ids, target.rows, strict=True):
            rows[tuple(row)].append(row_id)
        obs_ids = {
            o.cells: o.id for o in catalog.observations if o.target_id == target.id
        }
        for cells, row_ids in rows.items():
            row_positions = {}
            observed = set()
            position = 0
            for i, base in enumerate(cells):
                if base not in ("", "-"):
                    row_positions[i] = position
                    if base in "ACGT":
                        observed.add(position)
                    position += 1
            covered = set()

            def binding(
                record, target=target, cells=cells, row_positions=row_positions
            ):
                sid = record["site_id"]
                if (
                    sid not in catalog.site_by_id
                    or catalog.site_by_id[sid].target_id != target.id
                ):
                    return None
                forward = record["strand"] == "+"
                index = (
                    record["alignment_anchor"] - 1
                    if forward
                    else record["alignment_anchor"]
                )
                selected = []
                while 0 <= index < len(cells) and len(selected) < len(
                    record["sequence"]
                ):
                    if cells[index] == "":
                        return None
                    if cells[index] != "-":
                        selected.append(index)
                    index += -1 if forward else 1
                if len(selected) != len(record["sequence"]):
                    return None
                selected.sort()
                genomic = "".join(cells[i] for i in selected)
                expected = (
                    record["sequence"]
                    if forward
                    else record["sequence"].translate(complement)[::-1]
                )
                if genomic != expected or any(b not in "ACGT" for b in genomic):
                    return None
                return row_positions[selected[0]], row_positions[selected[-1]] + 1

            for group in by_config.values():
                forward = [binding(r) for r in group if r["strand"] == "+"]
                reverse = [binding(r) for r in group if r["strand"] == "-"]
                for f in forward:
                    for r in reverse:
                        if f is not None and r is not None and f[1] <= r[0]:
                            covered.update(range(f[1], r[0]))
            result.append(
                {
                    "target_id": target.id,
                    "allele_id": obs_ids.get(cells),
                    "row_ids": row_ids,
                    "covered_bases": len(covered & observed),
                    "observed_bases": len(observed),
                }
            )
    return result


def audit_allele_stage(directory, *, history=None):
    directory = Path(directory)
    manifest = _read(directory / "stage.json")
    if manifest.get("schemaVersion") != "primalscheme3.allele-stage/v2":
        raise ValueError("unknown stage schema")
    required = {
        "catalog.json.gz",
        "ledger.json.gz",
        "authoritative-targets.json.gz",
        "assignments.json",
        "coverage.json",
        "primer.bed",
        "primers.fasta",
        "order-sheet.tsv",
        "amplicon.bed",
        "primertrim.amplicon.bed",
        "reference.fasta",
    }
    if not required.issubset(manifest.get("artifacts", {})):
        raise ValueError("artifact integrity: missing required descriptors")
    for name, expected in manifest["artifacts"].items():
        path = directory / name
        if path.resolve().parent != directory.resolve() or expected.get("path") != name:
            raise ValueError("artifact integrity: unsafe path")
        if not path.is_file() or artifact_descriptor(path, directory) != expected:
            raise ValueError("artifact integrity: " + name)
    catalog = VariantCatalog.from_dict(_read(directory / "catalog.json.gz"))
    ledger = ConfigurationLedger.from_dict(_read(directory / "ledger.json.gz"))
    targets = _targets(_read(directory / "authoritative-targets.json.gz"))
    assignments = tuple(
        AlleleAssignment.from_dict(v) for v in _read(directory / "assignments.json")
    )
    profile = AlleleConstraintProfile(**manifest["constraints"])
    policy = StagePolicy(**manifest["stage_policy"])
    if manifest["stage_id"] != policy.stage_id:
        raise ValueError("stage identity mismatch")
    if (
        manifest["catalog_semantic_digest"] != catalog.semantic_digest
        or manifest["ledger_semantic_digest"] != ledger.semantic_digest
    ):
        raise ValueError("scientific semantic digest mismatch")
    records, violations = _bed_records(
        directory / "primer.bed", manifest, catalog, ledger
    )
    report = validate_allele_assignments(
        catalog,
        assignments,
        ledger,
        profile,
        policy=policy,
        goal=manifest["goal"],
        expected_summary=_read(directory / "coverage.json"),
        selected_sites=records,
        authoritative_targets=targets,
        history=history,
    )
    texts, _ = _render(
        catalog, ledger, assignments, manifest["references"], manifest["stage_id"]
    )
    for name, expected in texts.items():
        if name != "primer.bed" and (directory / name).read_text() != expected:
            violations.append(
                {"reason": "derived-selected-output-mismatch", "artifact": name}
            )
    literal = _literal_coverage(targets, records, catalog)
    actual = {
        r["allele_id"]: (r["covered_bases"], r["observed_bases"]) for r in literal
    }
    for record in report["coverage"]["classes"]:
        if actual.get(record["allele_id"]) != (
            record["covered_count"],
            record["observed_count"],
        ):
            violations.append(
                {
                    "reason": "literal-coverage-mismatch",
                    "allele_id": record["allele_id"],
                }
            )
    report["violations"].extend(violations)
    report["valid"] = not report["violations"]
    report["literal_coverage"] = literal
    return report


def publish_allele_stage(
    directory,
    *,
    catalog,
    ledger,
    assignments,
    profile,
    authoritative_targets,
    references,
    policy=None,
    goal=0.95,
    history=None,
):
    policy = StagePolicy() if policy is None else policy
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError(directory)
    if set(references) != set(catalog.target_by_id) or len(
        set(references.values())
    ) != len(references):
        raise ValueError("unique reference name required for every target")
    if any(
        not name or any(c.isspace() for c in name) or ";" in name
        for name in references.values()
    ):
        raise ValueError("invalid export reference name")
    selected = _selected(catalog, ledger, assignments)
    validation = validate_allele_assignments(
        catalog,
        assignments,
        ledger,
        profile,
        policy=policy,
        goal=goal,
        authoritative_targets=authoritative_targets,
        selected_sites=selected,
        history=history,
    )
    if not validation["valid"]:
        raise ValueError(
            "independent stage validation failed: "
            + json.dumps(validation["violations"])
        )
    directory.parent.mkdir(parents=True, exist_ok=True)
    pending = Path(
        tempfile.mkdtemp(
            prefix="." + directory.name + "-pending-", dir=directory.parent
        )
    )
    # Failed staging directories are deliberately retained for failure provenance.
    _write(pending / "catalog.json.gz", catalog.to_dict())
    _write(pending / "ledger.json.gz", ledger.to_dict())
    _write(
        pending / "authoritative-targets.json.gz",
        [asdict(t) for t in authoritative_targets],
    )
    _write(pending / "assignments.json", [a.to_dict() for a in assignments])
    _write(pending / "coverage.json", validation["coverage"])
    texts, names = _render(catalog, ledger, assignments, references, policy.stage_id)
    for name, text in texts.items():
        (pending / name).write_text(text)
    manifest = {
        "schemaVersion": "primalscheme3.allele-stage/v2",
        "stage_id": policy.stage_id,
        "stage_policy": asdict(policy),
        "constraints": asdict(profile),
        "goal": goal,
        "references": references,
        "configuration_names": names,
        "catalog_semantic_digest": catalog.semantic_digest,
        "ledger_semantic_digest": ledger.semantic_digest,
        "artifacts": {
            p.name: artifact_descriptor(p, pending)
            for p in pending.iterdir()
            if p.is_file()
        },
    }
    _write(pending / "stage.json", manifest)
    report = audit_allele_stage(pending, history=history)
    if not report["valid"]:
        raise ValueError(
            "saved-byte stage validation failed: " + json.dumps(report["violations"])
        )
    _write(pending / "validation.json", report)
    manifest["artifacts"]["validation.json"] = artifact_descriptor(
        pending / "validation.json", pending
    )
    _write(pending / "stage.json", manifest)
    pending.rename(directory)
    return {
        "path": str(directory),
        "stage_id": policy.stage_id,
        "validation": report,
        "manifest": manifest,
    }
