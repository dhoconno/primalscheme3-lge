from pathlib import Path

from scripts.benchmark_legacy_salvage import geometry, strict_fingerprint, stamp


def test_benchmark_fingerprint_uses_strict_bed_pool_and_ignores_generated_name(tmp_path):
    panel = tmp_path / "panel"
    panel.mkdir()
    rows = "chrA\t10\t14\tUUID_alpha\t1\t+\tACGT\nchrA\t30\t34\tUUID_alpha\t1\t-\tTGCA\n"
    (panel / "strict-primer.bed").write_text(rows)
    (panel / "primer.bed").write_text(rows.replace("UUID_alpha", "UUID_other"))
    first = strict_fingerprint(panel)
    (panel / "strict-primer.bed").write_text(rows.replace("\t1\t", "\t2\t").replace("UUID_alpha", "UUID_other"))
    assert first != strict_fingerprint(panel)


def test_benchmark_geometry_reports_per_target_unions_and_denominators(tmp_path):
    panel = tmp_path / "panel"
    panel.mkdir()
    (panel / "reference.fasta").write_text(">chrA\nAAAAAA\n>chrB\nAAAAAAAA\n")
    (panel / "amplicon.bed").write_text("chrA\t0\t4\tp\t1\nchrA\t3\t6\tq\t2\nchrB\t1\t5\tr\t1\n")
    (panel / "primertrim.amplicon.bed").write_text("chrA\t1\t3\tp\t1\nchrA\t2\t5\tq\t2\nchrB\t2\t4\tr\t1\n")
    report = geometry(panel)
    assert report["referenceLengths"] == {"chrA": 6, "chrB": 8}
    assert report["fullSpan"]["chrA"]["unionBases"] == 6
    assert report["primerTrimmedInterior"]["chrA"]["percentReferenceCovered"] == 100 * 4 / 6


def test_duplicate_basename_inputs_are_distinct_stamped_records(tmp_path):
    left = tmp_path / "a" / "primary.aligned.fasta"
    right = tmp_path / "b" / "primary.aligned.fasta"
    left.parent.mkdir(); right.parent.mkdir()
    left.write_text(">a\nAAAA\n"); right.write_text(">b\nCCCC\n")
    records = [stamp(left), stamp(right)]
    assert [Path(item["path"]).parent.name for item in records] == ["a", "b"]
    assert records[0]["sha256"] != records[1]["sha256"]
