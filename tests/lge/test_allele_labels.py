"""Display-only allele labels remain bound to saved raw input order."""

from __future__ import annotations

import pytest

from primalscheme3.panel.allele_coverage import canonical_observations
from primalscheme3.panel.allele_pipeline import _authoritative_targets
from primalscheme3.panel.coverage_types import VariantCatalog


def _catalog(bundle, records):
    targets = _authoritative_targets(bundle, records)
    return VariantCatalog(
        targets,
        tuple(
            observation
            for target in targets
            for observation in canonical_observations(target)
        ),
        (),
        (),
    )


def test_label_map_preserves_multi_input_order_descriptions_and_duplicate_classes(
    tmp_path,
):
    from primalscheme3.panel.allele_labels import build_allele_label_map

    (tmp_path / "work").mkdir()
    first = tmp_path / "work/0000-first.fa"
    second = tmp_path / "work/0001-second.fa"
    repeated = "ACGT" * 8
    first.write_text(
        ">alpha first display label\n"
        + repeated
        + "\n>beta second display label\n"
        + repeated
        + "\n"
    )
    second.write_text(">gamma third display label\n" + "TGCA" * 8 + "\n")
    records = [
        {"storedPath": "work/0000-first.fa", "sourceIndex": 0},
        {"storedPath": "work/0001-second.fa", "sourceIndex": 1},
    ]
    catalog = _catalog(tmp_path, records)

    result = build_allele_label_map(
        catalog.targets, catalog.observations, tmp_path, records
    )

    assert result["schemaVersion"] == "primalscheme3.allele-label-map/v1"
    assert result["scientificIdentityRole"] == "display-only-excluded"
    assert [item["source_msa_index"] for item in result["targets"]] == [0, 1]
    first_target = result["targets"][0]
    assert [row["row_index"] for row in first_target["rows"]] == [0, 1]
    assert [row["fasta_record_id"] for row in first_target["rows"]] == [
        "alpha",
        "beta",
    ]
    assert [row["fasta_description"] for row in first_target["rows"]] == [
        "alpha first display label",
        "beta second display label",
    ]
    assert [row["fasta_comment"] for row in first_target["rows"]] == [
        "first display label",
        "second display label",
    ]
    assert [row["normalized_record_id"] for row in first_target["rows"]] == [
        "alpha",
        "beta",
    ]
    assert len(first_target["classes"]) == 1
    assert first_target["classes"][0]["multiplicity"] == 2
    assert first_target["classes"][0]["row_ids"] == [
        row["row_id"] for row in first_target["rows"]
    ]
    assert first_target["classes"][0]["fasta_descriptions"] == [
        "alpha first display label",
        "beta second display label",
    ]


def test_reused_catalog_can_publish_current_headers_without_identity_change(tmp_path):
    from primalscheme3.panel.allele_labels import build_allele_label_map

    (tmp_path / "work").mkdir()
    raw = tmp_path / "work/input.fa"
    sequence = "ACGT" * 8
    raw.write_text(">old-name old display\n" + sequence + "\n")
    records = [{"storedPath": "work/input.fa", "sourceIndex": 0}]
    catalog = _catalog(tmp_path, records)
    old = build_allele_label_map(
        catalog.targets, catalog.observations, tmp_path, records
    )

    raw.write_text(">new-name renamed display\n" + sequence + "\n")
    renamed_catalog = _catalog(tmp_path, records)
    new = build_allele_label_map(
        catalog.targets, catalog.observations, tmp_path, records
    )

    assert renamed_catalog.semantic_digest == catalog.semantic_digest
    assert renamed_catalog.targets == catalog.targets
    assert old["targets"][0]["rows"][0]["fasta_description"] == "old-name old display"
    assert (
        new["targets"][0]["rows"][0]["fasta_description"] == "new-name renamed display"
    )
    assert (
        old["targets"][0]["rows"][0]["row_id"] == new["targets"][0]["rows"][0]["row_id"]
    )


@pytest.mark.parametrize("change", ["count", "content", "source-index"])
def test_label_map_rejects_raw_input_detached_from_catalog(tmp_path, change):
    from primalscheme3.panel.allele_labels import build_allele_label_map

    (tmp_path / "work").mkdir()
    raw = tmp_path / "work/input.fa"
    raw.write_text(">one\n" + "ACGT" * 8 + "\n")
    records = [{"storedPath": "work/input.fa", "sourceIndex": 0}]
    catalog = _catalog(tmp_path, records)
    if change == "count":
        raw.write_text(raw.read_text() + ">two\n" + "ACGT" * 8 + "\n")
    elif change == "content":
        raw.write_text(">one\n" + "TGCA" * 8 + "\n")
    else:
        records[0]["sourceIndex"] = 7

    with pytest.raises(ValueError, match="input|row|source"):
        build_allele_label_map(catalog.targets, catalog.observations, tmp_path, records)


def test_label_map_rejects_duplicate_native_record_ids_without_parser_change(tmp_path):
    (tmp_path / "work").mkdir()
    raw = tmp_path / "work/input.fa"
    raw.write_text(">same first\nACGT\n>same second\nACGT\n")
    # The protected native parser already rejects duplicate normalized record IDs;
    # the display bridge does not broaden that scientific input contract.
    with pytest.raises(Exception, match="Duplicate ID"):
        _authoritative_targets(
            tmp_path,
            [{"storedPath": "work/input.fa", "sourceIndex": 0}],
        )
