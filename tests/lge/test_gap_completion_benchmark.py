from scripts.benchmark_gap_completion import geometry


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
