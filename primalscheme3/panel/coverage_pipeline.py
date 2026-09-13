"""Coverage-only search, detached native publication, and artifact integration."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from time import monotonic
from types import SimpleNamespace
from typing import Any

import dnaio
import numpy as np
from click import UsageError
from primalbedtools.scheme import Scheme
from primalschemers import FKmer, RKmer

from primalscheme3.core.bedfiles import create_amplicon_str, create_bedfile_str
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.create_report_data import generate_all_plotdata
from primalscheme3.core.create_reports import generate_all_plots_html
from primalscheme3.core.primer_visual import (
    plot_primer_thermo_profile_html,
    primer_mismatch_heatmap,
)
from primalscheme3.panel.coverage_catalog import build_catalog, write_catalog
from primalscheme3.panel.coverage_provenance import finalize_provenance
from primalscheme3.panel.coverage_search import SearchOptions, optimize_catalog
from primalscheme3.panel.coverage_types import Assignment, Catalog
from primalscheme3.panel.coverage_validation import ConstraintProfile


@dataclass(frozen=True)
class CoveragePublication:
    target_to_reference: dict[str, str]
    candidate_to_amplicon: dict[str, str]
    artifact_paths: dict[str, str]
    projection: str = "ungapped-first-reference-coordinates"


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _write_json(path: Path, value: Any) -> None:
    path.write_bytes(_json_bytes(value))


def _counts_for(native: Any, sequences: tuple[str, ...]) -> list[float] | None:
    counts = list(native.counts())
    if not counts:
        return None
    count_by_sequence = {
        str(sequence).upper(): float(count)
        for sequence, count in zip(native.seqs(), counts, strict=True)
    }
    try:
        return [count_by_sequence[sequence] for sequence in sequences]
    except KeyError as error:
        raise ValueError(
            "native primer cloud disagrees with frozen candidate"
        ) from error


def _clone_fkmer(native: Any) -> FKmer:
    sequences = tuple(str(item).upper() for item in native.seqs())
    counts = list(native.counts()) or None
    return FKmer([item.encode() for item in sequences], int(native.end), counts)


def _clone_rkmer(native: Any) -> RKmer:
    sequences = tuple(str(item).upper() for item in native.seqs())
    counts = list(native.counts()) or None
    return RKmer([item.encode() for item in sequences], int(native.start), counts)


def _export_references(catalog: Catalog, msa_dict: dict[int, Any]) -> dict[str, str]:
    bases = {
        target.id: str(msa_dict[target.source_msa_index]._chrom_name)
        for target in catalog.targets
    }
    frequencies = {
        base: list(bases.values()).count(base) for base in set(bases.values())
    }
    result: dict[str, str] = {}
    occupied: set[str] = set()
    for target in sorted(catalog.targets, key=lambda item: item.id):
        base = bases[target.id]
        proposed = base if frequencies[base] == 1 else f"{base[:160]}__{target.id}"
        if proposed in occupied:
            proposed = f"{proposed}__occurrence-{target.occurrence}"
        if proposed in occupied:
            raise ValueError("coverage export reference IDs are not unique")
        result[target.id] = proposed
        occupied.add(proposed)
    return result


def _catalog_ranks(catalog: Catalog) -> dict[str, int]:
    result = {}
    for target in catalog.targets:
        candidates = sorted(
            (item for item in catalog.candidates if item.target_id == target.id),
            key=lambda item: (item.full_interval, item.interior_interval, item.id),
        )
        result.update(
            {candidate.id: index for index, candidate in enumerate(candidates, 1)}
        )
    return result


def _detached_pairs(
    catalog: Catalog,
    assignments: tuple[Assignment, ...],
    references: dict[str, str],
) -> tuple[list[PrimerPair], dict[str, str]]:
    ranks = _catalog_ranks(catalog)
    records = []
    names = {}
    for assignment in assignments:
        candidate = catalog.candidate_by_id[assignment.candidate_id]
        target = catalog.target_by_id[candidate.target_id]
        native = catalog.native_pair_by_candidate_id[candidate.id]
        forward = FKmer(
            [item.encode() for item in candidate.forward_oligos],
            candidate.interior_interval[0],
            _counts_for(native.fprimer, candidate.forward_oligos),
        )
        reverse = RKmer(
            [item.encode() for item in candidate.reverse_oligos],
            candidate.interior_interval[1],
            _counts_for(native.rprimer, candidate.reverse_oligos),
        )
        pair = PrimerPair(
            forward,
            reverse,
            target.source_msa_index,
            amplicon_number=ranks[candidate.id],
            pool=assignment.pool,
        )
        pair.chrom_name = references[target.id]
        pair.amplicon_prefix = f"coverage-{candidate.id}"
        if (min(forward.starts()), max(reverse.ends())) != candidate.full_interval:
            raise ValueError("detached full interval disagrees with frozen candidate")
        if (forward.end, reverse.start) != candidate.interior_interval:
            raise ValueError(
                "detached interior interval disagrees with frozen candidate"
            )
        if (
            tuple(forward.seqs()) != candidate.forward_oligos
            or tuple(reverse.seqs()) != candidate.reverse_oligos
        ):
            raise ValueError("detached primer cloud disagrees with frozen candidate")
        names[candidate.id] = f"{pair.amplicon_prefix}_{pair.amplicon_number}"
        records.append((target.id, ranks[candidate.id], candidate.id, pair))
    records.sort(key=lambda item: item[:3])
    return [record[-1] for record in records], names


def _projected_view(target: Any, source: Any, reference: str) -> Any:
    keep = [
        index
        for index, coordinate in enumerate(target.mapping)
        if coordinate is not None
    ]
    array = np.array(target.rows, dtype="U1")[:, keep]
    return SimpleNamespace(
        array=array,
        _mapping_array=None,
        fkmers=[_clone_fkmer(item) for item in source.fkmers],
        rkmers=[_clone_rkmer(item) for item in source.rkmers],
        regions=None,
        msa_index=target.source_msa_index,
        _chrom_name=reference,
    )


def publish_coverage_outputs(
    *,
    catalog: Catalog,
    assignments: tuple[Assignment, ...],
    validation: dict[str, Any],
    msa_dict: dict[int, Any],
    output_dir: Path,
    config: Config,
    offline_plots: bool,
) -> CoveragePublication:
    """Publish frozen assignments through detached native-format objects."""

    if validation.get(
        "schemaVersion"
    ) != "primalscheme3.panel-validation/v1" or not validation.get("valid"):
        raise ValueError("independent panel-v1 validation must pass before publication")
    expected = tuple(
        (item["candidate_id"], item["pool"])
        for item in validation.get("assignments", [])
    )
    actual = tuple((item.candidate_id, item.pool) for item in assignments)
    if expected != actual:
        raise ValueError("validation assignments disagree with frozen search result")

    output_dir = Path(output_dir)
    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    references = _export_references(catalog, msa_dict)
    pairs, names = _detached_pairs(catalog, assignments, references)

    primer = create_bedfile_str(
        ["# artic-bed-version v3.0", "# pc=PrimerCountInMSA"], pairs
    )
    # Native parser/serializer is an extra format check before bytes are public.
    primer = Scheme.from_str(primer).to_str()
    (output_dir / "primer.bed").write_text(primer)
    (output_dir / "amplicon.bed").write_text(create_amplicon_str(pairs, False))
    (output_dir / "primertrim.amplicon.bed").write_text(
        create_amplicon_str(pairs, True)
    )
    with dnaio.FastaWriter(output_dir / "reference.fasta", line_length=60) as stream:
        for target in sorted(catalog.targets, key=lambda item: item.id):
            stream.write(
                dnaio.SequenceRecord(references[target.id], target.reference_sequence)
            )

    (work / "primer_thermo.html").write_text(
        plot_primer_thermo_profile_html(
            output_dir / "primer.bed", config, offline_plots=offline_plots
        )
    )
    views = [
        _projected_view(
            target, msa_dict[target.source_msa_index], references[target.id]
        )
        for target in sorted(catalog.targets, key=lambda item: item.id)
    ]
    plot_data = generate_all_plotdata(views, work, last_pp_added=pairs)
    (output_dir / "plot.html").write_text(
        generate_all_plots_html(plot_data, output_dir, offline_plots=offline_plots)
    )
    mismatch = []
    for index, target in enumerate(sorted(catalog.targets, key=lambda item: item.id)):
        source = msa_dict[target.source_msa_index]
        original = dict(source._seq_dict)
        first_key = next(iter(original), references[target.id])
        original.pop(first_key, None)
        seqdict = {references[target.id]: "".join(target.rows[0]), **original}
        try:
            mismatch.append(
                primer_mismatch_heatmap(
                    array=np.array(target.rows, dtype="U1"),
                    seqdict=seqdict,
                    bedfile=output_dir / "primer.bed",
                    offline_plots=offline_plots and index == 0,
                    mapping=MappingType.FIRST,
                )
            )
        except UsageError:
            continue
    (output_dir / "primer.html").write_text("\n".join(mismatch))
    return CoveragePublication(
        target_to_reference=references,
        candidate_to_amplicon=names,
        artifact_paths={
            "primerBed": "primer.bed",
            "ampliconBed": "amplicon.bed",
            "primerTrimmedAmpliconBed": "primertrim.amplicon.bed",
            "referenceFasta": "reference.fasta",
            "plot": "plot.html",
            "mismatchPlot": "primer.html",
            "plotData": "work/plotdata.json.gz",
            "thermoPlot": "work/primer_thermo.html",
        },
    )


def run_coverage_pipeline(
    *,
    msa_dict: dict[int, Any],
    msa_data: dict[int, Any],
    input_records: list[dict[str, str]],
    output_dir: Path,
    config: Config,
    config_dict: dict[str, Any],
    max_amplicons: int | None,
    max_amplicons_msa: int | None,
    offline_plots: bool,
    logger: Any,
    argv: list[str],
    started_at: float,
    execution_start: dict[str, Any] | None,
) -> CoveragePublication:
    """Build, optimize, publish, and finally hash a coverage panel run."""

    try:
        resolved = dict(config.to_dict())
        resolved["output"] = "."
        catalog = build_catalog(msa_dict, config)
        catalog = replace(
            catalog,
            resolved_config_json=json.dumps(
                resolved, sort_keys=True, separators=(",", ":")
            ),
        )
        catalog_descriptor = write_catalog(
            catalog, output_dir / "candidate-catalog.json.gz"
        )
        profile = ConstraintProfile.from_config(
            config, max_amplicons=max_amplicons, max_amplicons_msa=max_amplicons_msa
        )
        options = SearchOptions(
            metric=config.coverage_metric,
            coverage_target=config.coverage_target,
            seed=config.optimizer_seed,
            starts=config.optimizer_starts,
            repair_rounds=config.optimizer_repair_rounds,
            time_limit=config.optimizer_time_limit,
        )
        result = optimize_catalog(catalog, profile, options)
        validation = result.metadata["validation"]
        if not validation.get("valid"):
            raise ValueError("independent panel-v1 validation rejected coverage result")
        publication = publish_coverage_outputs(
            catalog=catalog,
            assignments=result.assignments,
            validation=validation,
            msa_dict=msa_dict,
            output_dir=output_dir,
            config=config,
            offline_plots=offline_plots,
        )
        _write_json(output_dir / "panel-validation.json", validation)
        optimizer = dict(result.metadata)
        optimizer["assignments"] = [
            {"candidate_id": item.candidate_id, "pool": item.pool}
            for item in result.assignments
        ]
        optimizer["catalog"] = {
            "path": "candidate-catalog.json.gz",
            "schemaVersion": "primalscheme3.coverage-catalog/v1",
            "semanticDigest": catalog.semantic_digest,
            "fileSha256": catalog_descriptor["file_sha256"],
            "fileSize": catalog_descriptor["file_size"],
        }
        optimizer["publication"] = {
            "targetToReference": publication.target_to_reference,
            "candidateToAmplicon": publication.candidate_to_amplicon,
            "artifacts": publication.artifact_paths,
            "projection": publication.projection,
        }
        _write_json(output_dir / "panel-optimizer.json", optimizer)
        config_dict["output"] = "."
        config_dict["msa_data"] = msa_data
        config_dict["panel_optimizer"] = {
            "schemaVersion": "primalscheme3.panel-native-integration/v1",
            "selectionAlgorithm": "coverage",
            "toolVersion": config.version,
            "profile": profile.to_dict(),
            "options": options.__dict__,
            "catalog": optimizer["catalog"],
            "optimizer": {
                "path": "panel-optimizer.json",
                "schemaVersion": optimizer["schemaVersion"],
            },
            "validation": {
                "path": "panel-validation.json",
                "schemaVersion": validation["schemaVersion"],
                "valid": validation["valid"],
            },
            "provenance": {
                "path": "panel-provenance.json",
                "schemaVersion": "primalscheme3.panel-provenance/v1",
            },
            "publication": optimizer["publication"],
        }
        _write_json(output_dir / "config.json", config_dict)
        logger.info("Completed Successfully")
        scientific = {
            "algorithm": optimizer["algorithm"],
            "optimizerSchemaVersion": optimizer["schemaVersion"],
            "profile": profile.to_dict(),
            "metric": options.metric,
            "coverageTarget": options.coverage_target,
            "seed": options.seed,
            "budgets": {
                "starts": options.starts,
                "repairRounds": options.repair_rounds,
                "timeLimit": options.time_limit,
            },
            "completedWork": {
                "completedStarts": optimizer["completed_starts"],
                "completedRepairRounds": optimizer["completed_repair_rounds"],
                "work": optimizer["work"],
                "stopReason": optimizer["stop_reason"],
            },
            "catalogSemanticDigest": catalog.semantic_digest,
            "catalogFileSha256": catalog_descriptor["file_sha256"],
            "baselineObjective": optimizer["baseline_objective"],
            "finalObjective": optimizer["final_objective"],
            "validationValid": True,
        }
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
        return publication
    except Exception as error:
        logger.exception("Coverage panel failed")
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
            scientific={"selectionAlgorithm": "coverage"},
            logger=logger,
            execution_start=execution_start,
        )
        raise
