"""Saved selected oligos must survive relocation and independent byte auditing."""

import gzip
import json
import shutil
from dataclasses import replace

import pytest
from test_allele_validation import dna, fixture, profile

from primalscheme3.panel.coverage_types import AlleleAssignment, ConfigurationLedger


def selected():
    catalog, configs, ledger = fixture(seq=dna(seed=0))
    return catalog, ledger, (AlleleAssignment(configs[0].id, 0, "strict"),)


def publish(path, assignments=None):
    from primalscheme3.panel.allele_publication import publish_allele_stage

    cat, ledger, chosen = selected()
    return publish_allele_stage(
        path,
        catalog=cat,
        ledger=ledger,
        assignments=chosen if assignments is None else assignments,
        profile=profile(),
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "example"},
    )


def test_round_trip_relocation_and_literal_trimmed_coverage(tmp_path):
    from primalscheme3.panel.allele_publication import audit_allele_stage

    stage = tmp_path / "original"
    publish(stage)
    moved = tmp_path / "moved"
    shutil.move(stage, moved)
    report = audit_allele_stage(moved)
    assert report["valid"], report
    assert (
        report["validation_scope"]
        == "external-authoritative-inputs-and-selected-artifacts"
    )
    bed = [
        line.split("\t")
        for line in (moved / "primer.bed").read_text().splitlines()
        if not line.startswith("#")
    ]
    assert len(bed) == 2
    assert [(int(line[1]), int(line[2]), line[5]) for line in bed] == [
        (0, 20, "+"),
        (180, 200, "-"),
    ]
    assert (moved / "primertrim.amplicon.bed").read_text().split("\t")[1:3] == [
        "20",
        "180",
    ]
    assert report["literal_coverage"][0]["covered_bases"] == 160
    assert report["literal_coverage"][0]["observed_bases"] == 700


@pytest.mark.parametrize(
    "artifact", ["primer.bed", "primers.fasta", "order-sheet.tsv", "coverage.json"]
)
def test_changed_final_bytes_are_blocking(tmp_path, artifact):
    from primalscheme3.panel.allele_publication import audit_allele_stage

    stage = tmp_path / "stage"
    publish(stage)
    with (stage / artifact).open("a") as f:
        f.write("tampered\n")
    with pytest.raises(ValueError, match="integrity"):
        audit_allele_stage(stage)


def test_forged_bed_even_with_updated_hash_is_scientifically_rejected(tmp_path):
    from primalscheme3.panel.allele_publication import (
        artifact_descriptor,
        audit_allele_stage,
    )

    stage = tmp_path / "stage"
    publish(stage)
    bed = stage / "primer.bed"
    lines = bed.read_text().splitlines()
    first = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    cells = lines[first].split("\t")
    cells[6] = "A" * len(cells[6])
    lines[first] = "\t".join(cells)
    bed.write_text("\n".join(lines) + "\n")
    manifest = json.loads((stage / "stage.json").read_text())
    manifest["artifacts"]["primer.bed"] = artifact_descriptor(bed, stage)
    (stage / "stage.json").write_text(json.dumps(manifest))
    report = audit_allele_stage(stage)
    assert not report["valid"]
    assert any(v["reason"] == "selected-output-mismatch" for v in report["violations"])


def test_valid_empty_below_goal_is_publishable(tmp_path):
    from primalscheme3.panel.allele_publication import audit_allele_stage

    stage = tmp_path / "empty"
    publish(stage, assignments=())
    assert audit_allele_stage(stage)["valid"]
    assert not [
        line
        for line in (stage / "primer.bed").read_text().splitlines()
        if not line.startswith("#")
    ]


def test_invalid_stage_never_publishes_target_directory(tmp_path):
    from primalscheme3.panel.allele_publication import publish_allele_stage

    cat, configs, ledger = fixture()  # this fixture has a real strict dimer violation
    path = tmp_path / "invalid"
    with pytest.raises(ValueError, match="validation"):
        publish_allele_stage(
            path,
            catalog=cat,
            ledger=ledger,
            assignments=(AlleleAssignment(configs[0].id, 0, "strict"),),
            profile=profile(),
            authoritative_targets=cat.targets,
            references={cat.targets[0].id: "example"},
        )
    assert not path.exists()


def test_pruned_parent_variant_retained_only_in_catalog(tmp_path):
    from primalscheme3.panel.allele_publication import (
        audit_allele_stage,
        publish_allele_stage,
    )
    from primalscheme3.panel.coverage_types import CandidateFamily, OligoSite
    from primalscheme3.panel.coverage_variants import make_configuration

    cat, ledger, _ = selected()
    base = ledger.configurations[0]
    old = cat.site_by_id[base.reverse_site_ids[0]]
    removed = OligoSite(
        old.target_id,
        "A" * 20,
        "-",
        old.alignment_anchor,
        old.reference_footprint,
        accepting_profile_ids=("normal",),
    )
    family = CandidateFamily(
        base.target_id,
        base.anchor_pair,
        base.forward_site_ids,
        (*base.reverse_site_ids, removed.id),
    )
    cat = replace(cat, sites=(*cat.sites, removed), families=(family,))
    subset = make_configuration(
        cat, family.id, base.forward_site_ids, base.reverse_site_ids
    )
    ledger = ConfigurationLedger(cat.semantic_digest, (subset,))
    stage = tmp_path / "pruned"
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=(AlleleAssignment(subset.id, 0, "strict"),),
        profile=profile(),
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "example"},
    )
    for name in ("primer.bed", "primers.fasta", "order-sheet.tsv"):
        assert removed.sequence not in (stage / name).read_text()
    with gzip.open(stage / "catalog.json.gz", "rt") as f:
        assert removed.id in {r["id"] for r in json.load(f)["sites"]}
    assert audit_allele_stage(stage)["valid"]


def test_literal_audit_handles_internal_gap_unknown_and_duplicate_rows(tmp_path):
    from primalscheme3.panel.allele_publication import (
        audit_allele_stage,
        publish_allele_stage,
    )

    seq = dna(seed=0)
    other = seq[:80] + "-" + seq[81:100] + "N" + seq[101:]
    cat, configs, ledger = fixture(seq=seq, rows=(seq, other, seq))
    stage = tmp_path / "gapped"
    publish_allele_stage(
        stage,
        catalog=cat,
        ledger=ledger,
        assignments=(AlleleAssignment(configs[0].id, 0, "strict"),),
        profile=profile(),
        authoritative_targets=cat.targets,
        references={cat.targets[0].id: "example"},
    )
    report = audit_allele_stage(stage)
    assert report["valid"]
    assert sorted(
        (r["observed_bases"], r["covered_bases"]) for r in report["literal_coverage"]
    ) == [(698, 158), (700, 160)]
    assert len(report["literal_coverage"]) == 2


def test_bed_name_must_match_fasta_and_order_identity_even_with_updated_hash(tmp_path):
    from primalscheme3.panel.allele_publication import (
        audit_allele_stage,
        artifact_descriptor,
    )

    stage = tmp_path / "named"
    publish(stage)
    bed = stage / "primer.bed"
    lines = bed.read_text().splitlines()
    i = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    cells = lines[i].split("\t")
    cells[3] = "wrong_configuration_RIGHT_99"
    lines[i] = "\t".join(cells)
    bed.write_text("\n".join(lines) + "\n")
    manifest = json.loads((stage / "stage.json").read_text())
    manifest["artifacts"]["primer.bed"] = artifact_descriptor(bed, stage)
    (stage / "stage.json").write_text(json.dumps(manifest))
    result = audit_allele_stage(stage)
    assert not result["valid"]
    assert any(v.get("artifact") == "primer.bed" for v in result["violations"])
