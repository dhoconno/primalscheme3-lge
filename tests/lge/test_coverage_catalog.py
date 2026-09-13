"""Synthetic contracts for the immutable coverage candidate catalogue."""

import gzip
import json
from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import numpy as np
import pytest
from primalschemers import FKmer, RKmer

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config
from primalscheme3.core.mapping import create_mapping, ref_index_to_msa
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.panel.coverage_catalog import build_catalog, write_catalog
from primalscheme3.panel.coverage_types import Assignment


REFERENCE = "AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC"


def synthetic_msa(
    source_index,
    rows,
    pairs,
    *,
    labels=None,
    name="synthetic",
):
    array = np.array([list(row) for row in rows], dtype="U1")
    mapping, array = create_mapping(array, 0)
    labels = labels or [f"source-{index}" for index in range(len(rows))]
    return SimpleNamespace(
        array=array,
        _mapping_array=mapping,
        _ref_to_msa=ref_index_to_msa(mapping),
        _seq_dict=dict(zip(labels, rows, strict=True)),
        primerpairs=pairs,
        msa_index=source_index,
        name=name,
        _chrom_name=name,
    )


def native_pair(source_index=0, *, forward=None, reverse=None):
    forward = forward or [REFERENCE[10:18]]
    reverse = reverse or [reverse_complement(REFERENCE[42:50])]
    return PrimerPair(
        FKmer([sequence.encode() for sequence in forward], 18),
        RKmer([sequence.encode() for sequence in reverse], 42),
        source_index,
    )


def basic_catalog(*, source_index=0, rows=None, pair=None, labels=None):
    rows = rows or [REFERENCE]
    pair = pair or native_pair(source_index)
    msa = synthetic_msa(source_index, rows, [pair], labels=labels)
    return build_catalog({source_index: msa}, Config())


def test_geometry_native_orientation_and_recreated_objects_have_stable_ids():
    first = basic_catalog(pair=native_pair())
    second = basic_catalog(pair=native_pair())

    candidate = first.candidates[0]
    assert candidate.full_interval == (10, 50)
    assert candidate.interior_interval == (18, 42)
    assert candidate.forward_oligos == (REFERENCE[10:18],)
    assert candidate.reverse_oligos == (
        reverse_complement(REFERENCE[42:50]),
    )
    assert candidate.reverse_oligos != (REFERENCE[42:50],)
    assert candidate.id == second.candidates[0].id
    assert first.semantic_digest == second.semantic_digest
    assert first.native_pair_by_candidate_id[candidate.id] is not second.native_pair_by_candidate_id[candidate.id]


def test_target_and_catalogue_ids_ignore_input_order_and_source_labels():
    pair_a = native_pair(3)
    pair_b = PrimerPair(
        FKmer([REFERENCE[20:28].encode()], 28),
        RKmer([reverse_complement(REFERENCE[50:58]).encode()], 50),
        8,
    )
    msa_a = synthetic_msa(3, [REFERENCE], [pair_a], labels=["alpha"], name="first")
    msa_b = synthetic_msa(8, [REFERENCE.replace("A", "T", 1)], [pair_b], labels=["beta"], name="second")
    forward = build_catalog({3: msa_a, 8: msa_b}, Config())

    renamed_a = synthetic_msa(
        101,
        [REFERENCE],
        [native_pair(101)],
        labels=["renamed"],
        name="elsewhere",
    )
    renamed_b_pair = PrimerPair(
        FKmer([REFERENCE[20:28].encode()], 28),
        RKmer([reverse_complement(REFERENCE[50:58]).encode()], 50),
        55,
    )
    renamed_b = synthetic_msa(
        55,
        [REFERENCE.replace("A", "T", 1)],
        [renamed_b_pair],
        labels=["other"],
        name="different",
    )
    reversed_input = build_catalog(
        {55: renamed_b, 101: renamed_a},
        Config(output="a-different-runtime-path", ncores=7),
    )

    assert tuple(target.id for target in forward.targets) == tuple(
        target.id for target in reversed_input.targets
    )
    assert forward.semantic_digest == reversed_input.semantic_digest
    assert tuple(forward.target_by_id) == tuple(reversed_input.target_by_id)


def test_duplicate_target_inputs_remain_distinct_deterministic_occurrences():
    first = synthetic_msa(2, [REFERENCE], [native_pair(2)])
    second = synthetic_msa(9, [REFERENCE], [native_pair(9)])

    catalog = build_catalog({9: second, 2: first}, Config())
    rebuilt = build_catalog({2: first, 9: second}, Config())

    assert len(catalog.targets) == 2
    assert len({target.id for target in catalog.targets}) == 2
    assert [target.source_msa_index for target in catalog.targets] == [2, 9]
    assert [target.id for target in catalog.targets] == [
        target.id for target in rebuilt.targets
    ]
    assert len(catalog.candidates) == 2
    assert len({candidate.target_id for candidate in catalog.candidates}) == 2


def test_gapped_first_row_uses_ungapped_reference_denominator_and_mapping():
    gapped_reference = REFERENCE[:15] + "--" + REFERENCE[15:]
    pair = native_pair(4)
    msa = synthetic_msa(4, [gapped_reference], [pair])

    target = build_catalog({4: msa}, Config()).targets[0]

    assert target.reference_sequence == REFERENCE
    assert target.reference_length == 60
    assert len(target.rows[0]) == 62
    assert target.mapping[15:18] == (None, None, 15)
    assert target.ref_to_alignment[14:17] == (14, 17, 18)
    assert target.ref_to_alignment[-1] == 62


def test_support_walks_row_bases_across_reference_gaps_and_indel_alleles():
    gapped_reference = list(REFERENCE[:13]) + ["-"] + list(REFERENCE[13:])
    deletion_row = (
        list(REFERENCE[:12]) + ["-", "-"] + list(REFERENCE[13:])
    )
    deletion_forward = REFERENCE[9:12] + REFERENCE[13:18]
    pair = native_pair(forward=[REFERENCE[10:18], deletion_forward])
    msa = synthetic_msa(0, [gapped_reference, deletion_row], [pair])

    catalog = build_catalog({0: msa}, Config())
    row_ids = catalog.targets[0].row_ids
    candidate = catalog.candidates[0]

    assert candidate.joint_rows == row_ids
    assert candidate.row_product_spans == (
        (row_ids[0], 10, 50),
        (row_ids[1], 9, 49),
    )


def test_joint_support_requires_both_exact_anchors_and_unknown_is_diagnostic():
    exact = REFERENCE
    forward_only = REFERENCE[:42] + "A" * 8 + REFERENCE[50:]
    reverse_only = REFERENCE[:10] + "T" * 8 + REFERENCE[18:]
    terminal_missing = [""] * 11 + list(REFERENCE[11:])
    ambiguous = REFERENCE[:12] + "R" + REFERENCE[13:]
    n_unknown = REFERENCE[:45] + "N" + REFERENCE[46:]
    reverse_terminal_missing = list(REFERENCE[:49]) + [""] * 11
    rows = [
        exact,
        forward_only,
        reverse_only,
        terminal_missing,
        ambiguous,
        n_unknown,
        reverse_terminal_missing,
    ]
    msa = synthetic_msa(0, rows, [native_pair()], labels=list("abcdefg"))

    candidate = build_catalog({0: msa}, Config()).candidates[0]
    row_ids = build_catalog({0: msa}, Config()).targets[0].row_ids

    assert candidate.joint_rows == (row_ids[0],)
    assert candidate.unknown_rows == (
        row_ids[3],
        row_ids[4],
        row_ids[5],
        row_ids[6],
    )
    assert candidate.row_product_spans == ((row_ids[0], 10, 50),)


def test_disjoint_forward_only_and_reverse_only_rows_have_no_joint_support():
    forward_only = REFERENCE[:42] + "A" * 8 + REFERENCE[50:]
    reverse_only = REFERENCE[:10] + "T" * 8 + REFERENCE[18:]
    catalog = basic_catalog(rows=[forward_only, reverse_only])

    assert catalog.candidates[0].joint_rows == ()
    assert catalog.candidates[0].row_product_spans == ()


def test_row_product_spans_retain_each_matching_cloud_length_combination():
    pair = native_pair(
        forward=[REFERENCE[10:18], REFERENCE[12:18]],
        reverse=[
            reverse_complement(REFERENCE[42:50]),
            reverse_complement(REFERENCE[42:48]),
        ],
    )
    catalog = basic_catalog(pair=pair)
    row_id = catalog.targets[0].row_ids[0]

    assert catalog.candidates[0].row_product_spans == (
        (row_id, 10, 48),
        (row_id, 10, 50),
        (row_id, 12, 48),
        (row_id, 12, 50),
    )


def test_records_are_frozen_and_id_maps_are_read_only():
    catalog = basic_catalog()

    with pytest.raises(FrozenInstanceError):
        catalog.targets[0].reference_length = 1
    with pytest.raises(FrozenInstanceError):
        catalog.candidates[0].id = "changed"
    with pytest.raises(TypeError):
        catalog.target_by_id["new"] = catalog.targets[0]
    assert catalog.target_by_id is catalog.target_by_id
    assert catalog.candidate_by_id is catalog.candidate_by_id
    assert (
        catalog.native_pair_by_candidate_id
        is catalog.native_pair_by_candidate_id
    )
    assignment = Assignment(candidate_id=catalog.candidates[0].id, pool=1)
    with pytest.raises(FrozenInstanceError):
        assignment.pool = 2


def test_writer_is_byte_deterministic_and_reports_actual_file_provenance(tmp_path):
    catalog = basic_catalog(labels=["header-is-not-identity"])
    first_path = tmp_path / "first.json.gz"
    second_path = tmp_path / "second.json.gz"

    first_metadata = write_catalog(catalog, first_path)
    second_metadata = write_catalog(catalog, second_path)

    assert first_path.read_bytes() == second_path.read_bytes()
    assert first_metadata["semantic_digest"] == catalog.semantic_digest
    assert first_metadata["file_sha256"] == second_metadata["file_sha256"]
    assert first_metadata["file_size"] == first_path.stat().st_size
    assert first_metadata["path"] == str(first_path)

    payload = json.loads(gzip.decompress(first_path.read_bytes()))
    assert payload["schema"] == "primalscheme3.coverage-catalog/v1"
    assert payload["semantic_digest"] == catalog.semantic_digest
    assert payload["targets"][0]["rows"] == [list(REFERENCE)]
    assert payload["targets"][0]["reference_sequence"] == REFERENCE
    assert payload["oligos"] == [
        {"id": payload["oligos"][0]["id"], "sequence": REFERENCE[10:18]},
        {
            "id": payload["oligos"][1]["id"],
            "sequence": reverse_complement(REFERENCE[42:50]),
        },
    ]
    assert payload["candidates"][0]["full_interval"] == [10, 50]
    assert payload["candidates"][0]["row_product_spans"][0][1:] == [10, 50]
