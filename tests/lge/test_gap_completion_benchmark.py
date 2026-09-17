from scripts.benchmark_gap_completion import geometry
import pytest


def test_gap_completion_geometry_reports_trimmed_full_and_complements(tmp_path):
    panel = tmp_path / "panel"
    panel.mkdir()
    (panel / "reference.fasta").write_text(">a\nAAAAAAAAAA\n")
    (panel / "primertrim.amplicon.bed").write_text("a\t2\t5\tp\t1\na\t4\t8\tq\t2\n")
    (panel / "amplicon.bed").write_text("a\t1\t9\tp\t1\n")
    report = geometry(panel)
    assert report["trimmed"]["a"]["intervals"] == [[2, 8]]
    assert report["trimmed"]["a"]["unionBases"] == 6
    assert report["trimmed"]["a"]["gaps"] == [[0, 2], [8, 10]]
    assert report["full"]["a"]["percentReferenceCovered"] == 80.0


@pytest.mark.parametrize("bed", ["a\t-1\t3\tp\t1\n", "a\t3\t3\tp\t1\n", "a\t3\t11\tp\t1\n", "b\t1\t3\tp\t1\n"])
def test_gap_completion_geometry_rejects_malformed_intervals(tmp_path, bed):
    panel = tmp_path / "panel"
    panel.mkdir()
    (panel / "reference.fasta").write_text(">a\nAAAAAAAAAA\n")
    (panel / "primertrim.amplicon.bed").write_text(bed)
    (panel / "amplicon.bed").write_text("a\t1\t2\tp\t1\n")
    with pytest.raises(ValueError):
        geometry(panel)


def test_gap_completion_geometry_rejects_duplicate_reference_identifiers(tmp_path):
    panel = tmp_path / "panel"
    panel.mkdir()
    (panel / "reference.fasta").write_text(">a\nAAAA\n>a\nCCCC\n")
    (panel / "primertrim.amplicon.bed").write_text("a\t0\t2\tp\t1\n")
    (panel / "amplicon.bed").write_text("a\t0\t3\tp\t1\n")
    with pytest.raises(ValueError, match="duplicate"):
        geometry(panel)
