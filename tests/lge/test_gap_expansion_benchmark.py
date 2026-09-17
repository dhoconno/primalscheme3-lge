from scripts.benchmark_gap_expansion import PAIR_BOUND, ANCHOR_BOUND, commands


def test_gap_expansion_benchmark_uses_explicit_bounded_search():
    argv = commands()
    assert "--gap-expansion" in argv
    assert argv[argv.index("--gap-expansion") + 1] == "bounded"
    assert argv[argv.index("--gap-expansion-max-anchors-per-msa") + 1] == str(ANCHOR_BOUND)
    assert argv[argv.index("--gap-expansion-max-pairs-per-msa") + 1] == str(PAIR_BOUND)
    assert "--gap-completion-parent" in argv
