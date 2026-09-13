"""Detached native publication and complete final-byte provenance."""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import shutil
import time
from types import SimpleNamespace

import numpy as np
from primalbedtools.bedfiles import BedLineParser
from primalschemers import FKmer, RKmer

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config
from primalscheme3.core.mapping import create_mapping, ref_index_to_msa
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.panel.coverage_catalog import build_catalog, write_catalog
from primalscheme3.panel.coverage_pipeline import (
    publish_coverage_outputs,
    run_coverage_pipeline,
)
from primalscheme3.panel.coverage_provenance import (
    capture_execution_identity,
    finalize_provenance,
)
from primalscheme3.panel.coverage_types import Assignment

REFERENCE = "ACGTTGCAACGTCAGTACGATCGTAGCTAGCATCGATGCTAGCTACGATCGTACGATGCATGCA"


def _msa(index, name="ref", *, gapped=False):
    row = REFERENCE[:22] + ("--" if gapped else "") + REFERENCE[22:]
    array = np.array([list(row)], dtype="U1")
    mapping, array = create_mapping(array, 0)
    forward = [REFERENCE[4:12], REFERENCE[6:12]]
    reverse = [
        reverse_complement(REFERENCE[48:56]),
        reverse_complement(REFERENCE[48:54]),
    ]
    pair = PrimerPair(
        FKmer([x.encode() for x in forward], 12, [0.6, 0.4]),
        RKmer([x.encode() for x in reverse], 48, [0.7, 0.3]),
        index,
    )
    return SimpleNamespace(
        array=array,
        _mapping_array=mapping,
        _ref_to_msa=ref_index_to_msa(mapping),
        _seq_dict={name: row},
        primerpairs=[pair],
        fkmers=[pair.fprimer],
        rkmers=[pair.rprimer],
        msa_index=index,
        name=name,
        _chrom_name=name,
        regions=None,
    )


def _validation(catalog, assignments):
    return {
        "schemaVersion": "primalscheme3.panel-validation/v1",
        "valid": True,
        "violations": [],
        "assignments": [
            {"candidate_id": item.candidate_id, "pool": item.pool}
            for item in assignments
        ],
    }


def _snapshot(msa):
    return (
        tuple((k.end, tuple(k.seqs()), tuple(k.counts())) for k in msa.fkmers),
        tuple((k.start, tuple(k.seqs()), tuple(k.counts())) for k in msa.rkmers),
    )


def test_detached_publication_preserves_clouds_gaps_and_shared_objects(tmp_path):
    msa = _msa(0, gapped=True)
    catalog = build_catalog({0: msa}, Config())
    candidate = catalog.candidates[0]
    assignments = (Assignment(candidate.id, 1),)
    before = _snapshot(msa)
    output = tmp_path / "native"
    (output / "work").mkdir(parents=True)

    publication = publish_coverage_outputs(
        catalog=catalog,
        assignments=assignments,
        validation=_validation(catalog, assignments),
        msa_dict={0: msa},
        output_dir=output,
        config=Config(),
        offline_plots=False,
    )

    assert _snapshot(msa) == before
    headers, bedlines = BedLineParser.from_file(output / "primer.bed")
    assert "# artic-bed-version v3.0" in headers
    assert len(bedlines) == 4
    assert [line.sequence for line in bedlines[:2]] == list(candidate.forward_oligos)
    assert [line.sequence for line in bedlines[2:]] == list(candidate.reverse_oligos)
    assert [float(line.attributes["pc"]) for line in bedlines] == [0.4, 0.6, 0.3, 0.7]
    assert {line.ipool for line in bedlines} == {1}
    assert publication.candidate_to_amplicon[candidate.id]
    assert (
        (output / "reference.fasta").read_text().replace("\n", "").endswith(REFERENCE)
    )
    with gzip.open(output / "work" / "plotdata.json.gz", "rt") as stream:
        plotdata = json.load(stream)
    record = next(iter(plotdata.values()))
    assert record["dims"][1] == len(REFERENCE)
    assert next(iter(record["amplicons"].values()))["s"] == candidate.full_interval[0]
    assert (output / "primer.html").exists()


def test_duplicate_reference_names_are_unique_and_empty_selection_is_publishable(
    tmp_path,
):
    first, second = _msa(0, "same"), _msa(1, "same")
    catalog = build_catalog({0: first, 1: second}, Config())
    output = tmp_path / "native"
    (output / "work").mkdir(parents=True)
    publication = publish_coverage_outputs(
        catalog=catalog,
        assignments=(),
        validation=_validation(catalog, ()),
        msa_dict={0: first, 1: second},
        output_dir=output,
        config=Config(),
        offline_plots=False,
    )
    assert len(set(publication.target_to_reference.values())) == 2
    assert (output / "primer.bed").read_text().startswith("# artic-bed-version v3.0")
    assert (output / "amplicon.bed").read_bytes() == b"\n"
    assert (output / "primertrim.amplicon.bed").read_bytes() == b"\n"
    with gzip.open(output / "work" / "plotdata.json.gz", "rt") as stream:
        plotdata = json.load(stream)
    assert len(plotdata) == 2
    assert all(not row["amplicons"] for row in plotdata.values())
    assert all(
        row["uncovered"] == {"0": len(REFERENCE) - 1} for row in plotdata.values()
    )


def test_provenance_describes_final_bytes_and_survives_move(tmp_path):
    output = tmp_path / "native"
    (output / "work").mkdir(parents=True)
    source = tmp_path / "source.fasta"
    source.write_bytes(b">x\nACGT\n")
    copied = output / "work" / "0000-source.fasta"
    shutil.copyfile(source, copied)
    (output / "config.json").write_text('{"panel_optimizer":{}}')
    (output / "panel-optimizer.json").write_text(
        '{"schemaVersion":"primalscheme3.panel-optimizer/v1"}'
    )
    descriptor = write_catalog(
        build_catalog({}, Config()), output / "candidate-catalog.json.gz"
    )

    provenance = finalize_provenance(
        output_dir=output,
        argv=["primalscheme3", "panel-create", "--msa", str(source)],
        resolved_options={"selection_algorithm": "coverage"},
        inputs=[
            {
                "sourcePath": str(source.resolve()),
                "storedPath": "work/0000-source.fasta",
            }
        ],
        started_at=1.0,
        ended_at=3.5,
        status="success",
        exit_status=0,
        stderr="",
        scientific={"catalogSemanticDigest": descriptor["semantic_digest"]},
    )
    assert (
        provenance["selfHashPolicy"]
        == "panel-provenance.json is excluded to avoid a self-hash cycle"
    )
    assert provenance["wallSeconds"] == 2.5
    assert provenance["command"]["argv"][3] == str(source)
    assert provenance["inputs"][0]["storedPath"] == "work/0000-source.fasta"
    assert (
        provenance["inputs"][0]["sha256"]
        == hashlib.sha256(source.read_bytes()).hexdigest()
    )
    assert (
        provenance["inputs"][0]["sourceAtStartSha256"]
        == provenance["inputs"][0]["sha256"]
    )
    assert provenance["sourceChangedDuringRun"] is False
    assert provenance["runtimeChangedDuringRun"] is False
    for record in provenance["outputs"]:
        path = output / record["path"]
        assert path.stat().st_size == record["size"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]

    moved = tmp_path / "moved"
    shutil.move(output, moved)
    persisted = json.loads((moved / "panel-provenance.json").read_text())
    for record in persisted["outputs"]:
        assert (
            hashlib.sha256((moved / record["path"]).read_bytes()).hexdigest()
            == record["sha256"]
        )


def test_durable_copy_remains_authoritative_if_historical_origin_changes(tmp_path):
    output = tmp_path / "native"
    (output / "work").mkdir(parents=True)
    source = tmp_path / "source.fasta"
    source.write_bytes(b">before\nACGT\n")
    execution_start = capture_execution_identity([source])
    shutil.copyfile(source, output / "work" / "0000-source.fasta")
    source.write_bytes(b">after\nTGCA\n")
    record = finalize_provenance(
        output_dir=output,
        argv=["primalscheme3"],
        resolved_options={},
        inputs=[{"sourcePath": str(source), "storedPath": "work/0000-source.fasta"}],
        started_at=1,
        ended_at=2,
        status="success",
        exit_status=0,
        stderr="",
        scientific={},
        execution_start=execution_start,
    )
    assert record["inputs"][0]["sourceAtStartMatchesStored"] is True
    assert record["inputs"][0]["sourceChangedDuringRun"] is True
    assert record["inputs"][0]["sha256"] == execution_start["inputs"][0]["sha256"]


def test_full_empty_pipeline_publishes_linked_versioned_contract(tmp_path):
    source = tmp_path / "input.fasta"
    source.write_text(f">ref\n{REFERENCE}\n")
    execution_start = capture_execution_identity([source])
    output = tmp_path / "native"
    (output / "work").mkdir(parents=True)
    stored = output / "work" / "0000-input.fasta"
    shutil.copyfile(source, stored)
    msa = _msa(0)
    msa.primerpairs = []
    msa.fkmers = []
    msa.rkmers = []
    config = Config(
        selection_algorithm="coverage",
        amplicon_size=200,
        amplicon_size_min=150,
        amplicon_size_max=280,
        amplicon_size_metric="reference-span",
        mismatch_product_size=2000,
        optimizer_starts=1,
        optimizer_repair_rounds=0,
        optimizer_time_limit=1,
    )
    run_coverage_pipeline(
        msa_dict={0: msa},
        msa_data={0: {"msa_path": "work/0000-input.fasta"}},
        input_records=[
            {
                "sourcePath": str(source),
                "storedPath": "work/0000-input.fasta",
                "sourceIndex": 0,
            }
        ],
        output_dir=output,
        config=config,
        config_dict={**config.to_dict(), "mode": "equal"},
        max_amplicons=0,
        max_amplicons_msa=0,
        offline_plots=False,
        logger=logging.getLogger("coverage-pipeline-test"),
        argv=["primalscheme3", "panel-create"],
        started_at=time.monotonic(),
        execution_start=execution_start,
    )
    config_json = json.loads((output / "config.json").read_text())
    assert config_json["offline_plots"] is False
    native = config_json["panel_optimizer"]
    assert native["toolVersion"] == "3.3.0+lge.3"
    assert native["options"]["coverage_target"] == 0.9
    assert native["catalog"]["schemaVersion"] == "primalscheme3.coverage-catalog/v1"
    assert native["validation"]["valid"] is True
    assert native["provenance"]["path"] == "panel-provenance.json"
    optimizer = json.loads((output / native["optimizer"]["path"]).read_text())
    assert optimizer["assignments"] == []
    assert optimizer["final_objective"][4] == 0
    provenance = json.loads((output / native["provenance"]["path"]).read_text())
    assert provenance["status"] == "success"
    assert provenance["resolvedOptions"]["offline_plots"] is False
    for descriptor in provenance["outputs"]:
        payload = (output / descriptor["path"]).read_bytes()
        assert len(payload) == descriptor["size"]
        assert hashlib.sha256(payload).hexdigest() == descriptor["sha256"]
