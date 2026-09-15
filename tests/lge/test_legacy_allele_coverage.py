"""Historical BED scoring on the production observed-allele metric."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

from primalscheme3.panel.allele_coverage import allele_summary, canonical_observations
from primalscheme3.panel.coverage_types import (
    AlleleAssignment,
    CandidateFamily,
    ConfigurationLedger,
    OligoSite,
    SelectionConfiguration,
    Target,
    VariantCatalog,
)

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "score_legacy_allele_coverage.py"


def _script_module():
    spec = importlib.util.spec_from_file_location("legacy_allele_coverage_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
    return path


def _bed(path: Path, rows: list[tuple[object, ...]]) -> Path:
    path.write_text(
        "# artic-bed-version v3.0\n"
        + "".join("\t".join(map(str, row)) + "\n" for row in rows)
    )
    return path


def _run(
    output: Path,
    msas: list[Path],
    bed: Path,
    *,
    row_maps: list[Path] | None = None,
    labels: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    argv = [sys.executable, str(SCRIPT)]
    for msa in msas:
        argv.extend(("--msa", str(msa)))
    for row_map in row_maps or []:
        argv.extend(("--row-map", str(row_map)))
    if labels is not None:
        argv.extend(("--labels", str(labels)))
    argv.extend(("--bed", str(bed), "--output", str(output)))
    return subprocess.run(
        argv, cwd=ROOT, text=True, capture_output=True, check=False
    )


def _production_summary() -> object:
    rows = tuple(
        tuple(row) for row in ("AA---CCNGTT", "AA---CCNGTT", "AAT-TCCNGTT")
    )
    ref_columns = tuple(i for i, base in enumerate(rows[0]) if base != "-")
    mapping = tuple(
        ref_columns.index(i) if i in ref_columns else None for i in range(len(rows[0]))
    )
    target = Target(
        "target-test",
        0,
        0,
        ("ref", "duplicate", "insertion"),
        rows,
        "AACCNGTT",
        8,
        mapping,
        ref_columns + (ref_columns[-1] + 1,),
    )
    forward = OligoSite(target.id, "AA", "+", 2, (0, 2))
    reverse = OligoSite(target.id, "AA", "-", 9, (6, 8))
    family = CandidateFamily(target.id, (2, 9), (forward.id,), (reverse.id,))
    config = SelectionConfiguration(
        target.id,
        family.id,
        (forward.id,),
        (reverse.id,),
        (0, 8),
        family.anchor_pair,
    )
    catalog = VariantCatalog(
        (target,), canonical_observations(target), (forward, reverse), (family,)
    )
    ledger = ConfigurationLedger(catalog.semantic_digest, (config,))
    return allele_summary(
        catalog, (AlleleAssignment(config.id, 0, "strict"),), ledger, goal=0.95
    )


def _production_empty_summary() -> object:
    rows = (tuple("AAAATTTT"),)
    target = Target(
        "target-empty",
        0,
        0,
        ("ref",),
        rows,
        "AAAATTTT",
        8,
        tuple(range(8)),
        tuple(range(9)),
    )
    catalog = VariantCatalog((target,), canonical_observations(target), (), ())
    return allele_summary(
        catalog, (), ConfigurationLedger(catalog.semantic_digest, ()), goal=0.95
    )


def test_scores_all_variants_with_production_metric_and_complete_provenance(tmp_path):
    msa = _fasta(
        tmp_path / "target.fasta",
        [
            ("ref-chrom", "AA---CCNGTT"),
            ("duplicate", "AA---CCNGTT"),
            ("insertion", "AAT-TCCNGTT"),
        ],
    )
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("ref_chrom", 0, 2, "amp_1_LEFT_1", 1, "+", "AA", "pc=2"),
            ("ref_chrom", 0, 2, "amp_1_LEFT_2", 1, "+", "CC", "pc=1"),
            ("ref_chrom", 6, 8, "amp_1_RIGHT_1", 1, "-", "AA", "pc=2"),
            ("ref_chrom", 6, 8, "amp_1_RIGHT_2", 1, "-", "CC", "pc=1"),
        ],
    )
    output = tmp_path / "score"

    completed = _run(output, [msa], bed)

    assert completed.returncode == 0, completed.stderr
    report = json.loads((output / "coverage.json").read_text())
    production = _production_summary()
    assert report["scientificScope"].startswith("BINDING-MODEL COVERAGE ONLY")
    assert report["coverage"]["metric_id"] == production.metric_id
    assert report["coverage"]["mean_coverage"] == production.mean_coverage
    assert sorted(c["fraction"] for c in report["coverage"]["classes"]) == sorted(
        c.fraction for c in production.classes
    )
    assert sorted(c["multiplicity"] for c in report["targets"][0]["classes"]) == [1, 2]
    assert report["counts"] == {
        "ampliconFamilies": 1,
        "distinctObservedAlleleClasses": 2,
        "inputRows": 3,
        "primerRecords": 4,
        "targets": 1,
    }
    assert report["pools"] == [
        {"ampliconFamilies": 1, "pool": 1, "primerRecords": 4}
    ]
    assert report["amplicons"][0]["possibleVariantCombinations"] == 4
    assert report["amplicons"][0]["supportedAlleleClasses"] == 2
    assert report["bindingSupportStatusCounts"] == {
        "confirmed": 4,
        "mismatch": 4,
        "unavailable": 0,
        "unknown": 0,
    }

    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["tool"] == {
        "name": "score-legacy-allele-coverage",
        "version": "1",
    }
    assert provenance["status"] == "success"
    assert provenance["exitStatus"] == 0
    assert provenance["stderr"] == ""
    assert provenance["wallTimeSeconds"] >= 0
    assert provenance["command"]["argv"] == [sys.executable, str(SCRIPT), "--msa", str(msa), "--bed", str(bed), "--output", str(output)]
    assert provenance["resolvedOptions"]["chemistryEligibilityFiltering"] is False
    assert provenance["resolvedOptions"]["specificityFiltering"] is False
    assert provenance["runtime"]["pythonExecutable"]
    assert "conda" in provenance["runtime"] and "container" in provenance["runtime"]
    assert {record["role"] for record in provenance["inputs"]} == {"bed", "msa"}
    outputs = {record["path"]: record for record in provenance["outputs"]}
    assert set(outputs) == {"coverage.json", "resolved-chromosome-maps.json"}
    for relative, descriptor in outputs.items():
        path = output / relative
        assert descriptor["size"] == path.stat().st_size
        assert descriptor["sha256"] == _sha256(path)


def test_explicit_row_map_and_labels_resolve_normalized_snapshot(tmp_path):
    msa = _fasta(
        tmp_path / "occurrence-a.fasta",
        [("input_occurrence_a_row_0", "AAAATTTT"), ("input_occurrence_a_row_1", "AAAATTTT")],
    )
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("original_ref", 0, 2, "legacy_1_LEFT_1", 2, "+", "AA", "pc=2"),
            ("original_ref", 6, 8, "legacy_1_RIGHT_1", 2, "-", "AA", "pc=2"),
        ],
    )
    row_map = tmp_path / "occurrence-a-row-map.json"
    row_map.write_text(
        json.dumps(
            {
                "schemaVersion": 1,
                "inputID": "occurrence-a",
                "rows": [
                    {
                        "rowIndex": 0,
                        "normalizedHeader": "input_occurrence_a_row_0",
                        "originalHeader": "original-ref",
                    },
                    {
                        "rowIndex": 1,
                        "normalizedHeader": "input_occurrence_a_row_1",
                        "originalHeader": "other-ref",
                    },
                ],
            }
        )
    )
    labels = tmp_path / "source-row-label-map.json"
    labels.write_text(
        json.dumps(
            {
                "schemaVersion": "allele-coverage-source-label-map/v1",
                "inputs": [
                    {"sourceOccurrenceID": "occurrence-a", "label": "MHC target A"}
                ],
            }
        )
    )
    output = tmp_path / "score"

    completed = _run(output, [msa], bed, row_maps=[row_map], labels=labels)

    assert completed.returncode == 0, completed.stderr
    mapping = json.loads((output / "resolved-chromosome-maps.json").read_text())["mappings"][0]
    assert mapping["chromosome"] == "original_ref"
    assert mapping["mappingEvidence"] == "explicit-row-map"
    assert mapping["sourceOccurrenceID"] == "occurrence-a"
    report = json.loads((output / "coverage.json").read_text())
    assert report["targets"][0]["label"] == "MHC target A"
    provenance = json.loads((output / "provenance.json").read_text())
    assert {record["role"] for record in provenance["inputs"]} == {
        "bed",
        "labels",
        "msa",
        "row-map",
    }


def test_explicit_row_map_resolves_normalized_bed_against_original_header_msa(tmp_path):
    msa = _fasta(
        tmp_path / "occurrence-a.fasta",
        [("original-ref", "AAAATTTT"), ("other-ref", "AAAATTTT")],
    )
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("input_occurrence_a_row_0", 0, 2, "legacy_1_LEFT_1", 2, "+", "AA", "pc=2"),
            ("input_occurrence_a_row_0", 6, 8, "legacy_1_RIGHT_1", 2, "-", "AA", "pc=2"),
        ],
    )
    row_map = tmp_path / "occurrence-a-row-map.json"
    row_map.write_text(
        json.dumps(
            {
                "schemaVersion": 1,
                "inputID": "occurrence-a",
                "rows": [
                    {
                        "rowIndex": 0,
                        "normalizedHeader": "input_occurrence_a_row_0",
                        "originalHeader": "original-ref",
                    },
                    {
                        "rowIndex": 1,
                        "normalizedHeader": "input_occurrence_a_row_1",
                        "originalHeader": "other-ref",
                    },
                ],
            }
        )
    )

    completed = _run(tmp_path / "score", [msa], bed, row_maps=[row_map])

    assert completed.returncode == 0, completed.stderr
    mapping = json.loads(
        (tmp_path / "score" / "resolved-chromosome-maps.json").read_text()
    )["mappings"][0]
    assert mapping["chromosome"] == "input_occurrence_a_row_0"
    assert mapping["mappingEvidence"] == "explicit-row-map"
    assert mapping["rowMapFirstOriginalHeader"] == "original-ref"
    assert mapping["rowMapFirstNormalizedHeader"] == "input_occurrence_a_row_0"


def test_explicit_row_map_rejects_bed_with_both_header_aliases(tmp_path):
    msa = _fasta(tmp_path / "occurrence-a.fasta", [("original-ref", "AAAATTTT")])
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("original_ref", 0, 2, "legacy_1_LEFT_1", 2, "+", "AA", "pc=2"),
            ("original_ref", 6, 8, "legacy_1_RIGHT_1", 2, "-", "AA", "pc=2"),
            ("input_occurrence_a_row_0", 0, 2, "legacy_2_LEFT_1", 2, "+", "AA", "pc=2"),
            ("input_occurrence_a_row_0", 6, 8, "legacy_2_RIGHT_1", 2, "-", "AA", "pc=2"),
        ],
    )
    row_map = tmp_path / "occurrence-a-row-map.json"
    row_map.write_text(
        json.dumps(
            {
                "schemaVersion": 1,
                "inputID": "occurrence-a",
                "rows": [
                    {
                        "rowIndex": 0,
                        "normalizedHeader": "input_occurrence_a_row_0",
                        "originalHeader": "original-ref",
                    }
                ],
            }
        )
    )

    completed = _run(tmp_path / "score", [msa], bed, row_maps=[row_map])

    assert completed.returncode != 0
    assert "row-map chromosome aliases occur in BED" in completed.stderr


def test_reports_unknown_and_unavailable_binding_without_coverage_credit(tmp_path):
    msa = _fasta(
        tmp_path / "target.fasta",
        [("ref", "AA--CCNGTT"), ("terminal-missing", "--TTCCNGTT")],
    )
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("ref", 0, 2, "amp_1_LEFT_1", 1, "+", "AA", "pc=1"),
            ("ref", 6, 8, "amp_1_RIGHT_1", 1, "-", "AA", "pc=2"),
            ("ref", 3, 5, "amp_2_LEFT_1", 1, "+", "CA", "pc=0"),
            ("ref", 6, 8, "amp_2_RIGHT_1", 1, "-", "AA", "pc=2"),
        ],
    )
    output = tmp_path / "score"

    completed = _run(output, [msa], bed)

    assert completed.returncode == 0, completed.stderr
    report = json.loads((output / "coverage.json").read_text())
    assert report["bindingSupportStatusCounts"] == {
        "confirmed": 5,
        "mismatch": 0,
        "unavailable": 1,
        "unknown": 2,
    }
    terminal = next(
        item
        for item in report["targets"][0]["classes"]
        if "terminal-missing" in item["rowIDs"]
    )
    assert terminal["bindingSupportStatusCounts"]["unavailable"] == 1
    assert terminal["fraction"] == 0


@pytest.mark.parametrize("failure", ["unmappable", "malformed-row-map"])
def test_mapping_failures_write_failure_provenance(tmp_path, failure):
    msa = _fasta(
        tmp_path / "occurrence-a.fasta",
        [("synthetic_row_0", "AAAATTTT"), ("synthetic_row_1", "AAAATTTT")],
    )
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("original_ref", 0, 2, "legacy_1_LEFT_1", 1, "+", "AA", "pc=2"),
            ("original_ref", 6, 8, "legacy_1_RIGHT_1", 1, "-", "AA", "pc=2"),
        ],
    )
    row_maps = []
    if failure == "malformed-row-map":
        row_map = tmp_path / "occurrence-a-row-map.json"
        row_map.write_text(
            json.dumps(
                {
                    "schemaVersion": 1,
                    "inputID": "occurrence-a",
                    "rows": [
                        {
                            "rowIndex": 0,
                            "normalizedHeader": "synthetic_row_0",
                            "originalHeader": "original-ref",
                        }
                    ],
                }
            )
        )
        row_maps = [row_map]
    output = tmp_path / "score"

    completed = _run(output, [msa], bed, row_maps=row_maps)

    assert completed.returncode != 0
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["status"] == "error"
    assert provenance["exitStatus"] == 1
    assert provenance["stderr"] in completed.stderr
    assert provenance["outputs"] == []
    assert all(record["sha256"] == _sha256(Path(record["path"])) for record in provenance["inputs"])
    if failure == "unmappable":
        assert "explicit --row-map" in completed.stderr
    else:
        assert "row count" in completed.stderr


def test_refuses_existing_output_without_mutation(tmp_path):
    msa = _fasta(tmp_path / "target.fasta", [("ref", "AAAATTTT")])
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("ref", 0, 2, "amp_1_LEFT_1", 1, "+", "AA", "pc=1"),
            ("ref", 6, 8, "amp_1_RIGHT_1", 1, "-", "AA", "pc=1"),
        ],
    )
    output = tmp_path / "score"
    output.mkdir()
    marker = output / "marker"
    marker.write_text("keep")

    completed = _run(output, [msa], bed)

    assert completed.returncode != 0
    assert "already exists" in completed.stderr
    assert list(output.iterdir()) == [marker]
    assert marker.read_text() == "keep"


def test_rejects_malformed_amplicon_family_with_failure_provenance(tmp_path):
    msa = _fasta(tmp_path / "target.fasta", [("ref", "AAAATTTT")])
    bed = _bed(
        tmp_path / "primer.bed",
        [("ref", 0, 2, "amp_1_LEFT_1", 1, "+", "AA", "pc=1")],
    )
    output = tmp_path / "score"

    completed = _run(output, [msa], bed)

    assert completed.returncode != 0
    assert "requires LEFT and RIGHT primers" in completed.stderr
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["status"] == "error"
    assert provenance["inputs"][0]["sha256"] == _sha256(msa)
    assert provenance["inputs"][1]["sha256"] == _sha256(bed)


def test_header_only_empty_bed_reports_production_zero_coverage(tmp_path):
    msa = _fasta(tmp_path / "target.fasta", [("ref", "AAAATTTT")])
    bed = tmp_path / "primer.bed"
    bed.write_text("# artic-bed-version v3.0\n# pc=PrimerCountInMSA\n")
    output = tmp_path / "score"

    completed = _run(output, [msa], bed)

    assert completed.returncode == 0, completed.stderr
    report = json.loads((output / "coverage.json").read_text())
    production = _production_empty_summary()
    assert report["valid"] is True
    assert report["coverage"]["mean_coverage"] == production.mean_coverage == 0
    assert report["coverage"]["classes"][0]["covered_count"] == 0
    assert report["counts"]["primerRecords"] == 0
    assert report["counts"]["ampliconFamilies"] == 0
    assert report["pools"] == []
    assert report["amplicons"] == []
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["status"] == "success"
    assert provenance["scientificReportValid"] is True


@pytest.mark.parametrize(
    ("changed_field", "reason"),
    [
        ("inputsChangedDuringRun", "scientific inputs changed during execution"),
        ("sourceChangedDuringRun", "scientific source changed during execution"),
    ],
)
def test_changed_inputs_or_source_invalidate_report_and_fail_closed(
    tmp_path, monkeypatch, capsys, changed_field, reason
):
    msa = _fasta(tmp_path / "target.fasta", [("ref", "AAAATTTT")])
    bed = _bed(
        tmp_path / "primer.bed",
        [
            ("ref", 0, 2, "amp_1_LEFT_1", 1, "+", "AA", "pc=1"),
            ("ref", 6, 8, "amp_1_RIGHT_1", 1, "-", "AA", "pc=1"),
        ],
    )
    output = tmp_path / "score"
    module = _script_module()
    actual_provenance = module._provenance

    def injected_change(**kwargs):
        receipt = actual_provenance(**kwargs)
        receipt[changed_field] = True
        return receipt

    monkeypatch.setattr(module, "_provenance", injected_change)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [str(SCRIPT), "--msa", str(msa), "--bed", str(bed), "--output", str(output)],
    )

    assert module.main() == 1

    assert reason in capsys.readouterr().err
    report = json.loads((output / "coverage.json").read_text())
    assert report["valid"] is False
    assert reason in report["invalidReasons"]
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["status"] == "error"
    assert provenance["exitStatus"] == 1
    assert provenance["scientificReportValid"] is False
    assert provenance["command"]["workingDirectory"] == str(tmp_path)
