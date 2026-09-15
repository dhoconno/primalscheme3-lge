from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from primalscheme3.cli import app, panel_create
from primalscheme3.core.config import Config


def options():
    from primalscheme3.panel.allele_options import AlleleOptions

    return AlleleOptions


def kwargs(**changes):
    return (
        dict(
            selection_algorithm="allele-coverage",
            amplicon_size=200,
            amplicon_size_min=150,
            amplicon_size_max=280,
        )
        | changes
    )


def test_programmatic_preset_resolves_and_roundtrips_exact_types():
    config = Config(**kwargs())
    assert config.coverage_target == 0.95
    assert config.coverage_metric == "observed-allele-primer-trimmed"
    assert config.mismatch_kmersize == 17 and config.mismatch_product_size == 2000
    assert config.terminal_gap_policy.value == "observed-only"
    assert config.amplicon_size_metric.value == "reference-span"
    resolved = options().from_config(config)
    assert resolved.preset == "allele-balanced-v1"
    assert resolved.variant_selection == "subsets"
    assert resolved.candidate_profiles == "union"
    assert resolved.salvage == "off" and resolved.primary_tier == "strict"
    assert resolved.salvage_thresholds == (-28.0, -30.0, -32.0)
    assert Config(**config.to_dict()).to_dict() == config.to_dict()
    assert resolved.to_dict()["requested_options"]["amplicon_size_min"] == 150
    assert "coverage_target" not in resolved.to_dict()["requested_options"]


def test_programmatic_overrides_and_search_mapping():
    config = Config(
        **kwargs(
            coverage_target=0.8,
            specificity_terminal_k=19,
            optimizer_seed=7,
            subset_expansion_limit=8,
            work_families_per_refresh=2,
            variant_selection="subsets",
            salvage="bounded",
            salvage_thresholds=(-29, -33),
            primary_tier="salvage-2",
        )
    )
    resolved = options().from_config(config)
    search = resolved.search_options()
    assert (
        search.seed == 7
        and search.coverage_target == 0.8
        and search.subset_expansion_limit == 8
    )
    assert search.work_limits.families_per_refresh == 2
    assert search.variant_selection == "subsets"
    assert (
        options()
        .from_config(Config(**kwargs(variant_selection="full-cloud")))
        .search_options()
        .variant_selection
        == "full-cloud"
    )
    assert config.mismatch_kmersize == 19
    assert [p.active_cutoff for p in resolved.stage_policies()] == [-26, -29, -33]


@pytest.mark.parametrize(
    "change",
    [
        {"high_gc": True},
        {"coverage_metric": "full-span"},
        {"dimer_score": -28},
        {"terminal_gap_policy": "legacy"},
        {"mode": "entropy"},
        {"mapping": "consensus"},
        {"circular": True},
        {"backtrack": True},
        {"downsample": True},
        {"use_annealing": True},
        {"use_matchdb": False},
        {"mismatch_product_size": 0},
        {"specificity_terminal_k": 0},
        {"allele_weighting": "row-frequency"},
        {"coverage_target": 0},
        {"coverage_target": float("nan")},
        {"subset_beam_width": True},
        {"exchange_width": 3},
        {"work_repair_trials_per_round": 0},
        {"primary_tier": "salvage-1"},
        {"salvage_thresholds": (-28,)},
        {"salvage": "bounded", "salvage_thresholds": (-28, -27)},
        {"salvage": "bounded", "salvage_thresholds": (-28, float("-inf"))},
        {"salvage": "bounded", "salvage_thresholds": (-28, -30, -32, -34)},
        {"salvage": "bounded", "primary_tier": "salvage-4"},
        {"primer_tm_min": 58},
        {"max_amplicons_region_group": 1},
    ],
)
def test_new_mode_rejects_incompatible_or_ignored_options(change):
    with pytest.raises(ValueError):
        Config(**kwargs(**change))


def invoke(tmp_path, extra=()):
    path = tmp_path / "input.fasta"
    path.write_text(">a\nACGT\n")
    args = [
        "panel-create",
        "--msa",
        str(path),
        "--output",
        str(tmp_path / "out"),
        "--selection-algorithm",
        "allele-coverage",
        "--amplicon-size",
        "200",
        "--amplicon-size-min",
        "150",
        "--amplicon-size-max",
        "280",
        *extra,
    ]
    with patch("primalscheme3.cli.panelcreate") as run:
        result = CliRunner().invoke(app, args)
    return result, run


def test_cli_default_sources_do_not_override_new_preset(tmp_path):
    result, run = invoke(tmp_path)
    assert result.exit_code == 0, result.output + repr(result.exception)
    config = run.call_args.kwargs["config"]
    resolved = options().from_config(config)
    assert (
        config.coverage_target == 0.95
        and config.terminal_gap_policy.value == "observed-only"
    )
    assert run.call_args.kwargs["mode"].value == "equal"
    assert "coverage_target" not in resolved.to_dict()["requested_options"]
    result, run = invoke(
        tmp_path, ("--coverage-target", ".90", "--specificity-terminal-k", "19")
    )
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["config"].coverage_target == 0.90
    assert (
        options()
        .from_config(run.call_args.kwargs["config"])
        .to_dict()["requested_options"]["coverage_target"]
        == 0.90
    )


@pytest.mark.parametrize(
    "flags",
    [
        ("--high-gc",),
        ("--terminal-gap-policy", "legacy"),
        ("--coverage-metric", "full-span"),
        ("--mispriming-product-size", "0"),
        ("--dimer-score", "-28"),
        ("--mode", "region-only"),
        ("--salvage-threshold", "-28"),
        ("--primary-tier", "salvage-1"),
    ],
)
def test_cli_preflight_rejects_before_workflow(tmp_path, flags):
    result, run = invoke(tmp_path, flags)
    assert result.exit_code == 2, result.output + repr(result.exception)
    run.assert_not_called()
    assert not (tmp_path / "out").exists()


def test_new_flags_rejected_for_legacy_and_direct_function_defaults(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">a\nACGT\n")
    result = CliRunner().invoke(
        app,
        [
            "panel-create",
            "--msa",
            str(path),
            "--output",
            str(tmp_path / "out"),
            "--candidate-profiles",
            "union",
        ],
    )
    assert result.exit_code == 2
    with patch("primalscheme3.cli.panelcreate") as run:
        panel_create(msa=[path], output=tmp_path / "out", **kwargs())
    assert run.call_args.kwargs["config"].coverage_target == 0.95


def test_allele_capabilities_are_additive_and_reserve_source_contract():
    from primalscheme3.panel.coverage_provenance import capabilities_document

    caps = capabilities_document()
    assert caps["selectionAlgorithms"] == ["legacy", "coverage", "allele-coverage"]
    assert caps["coverage"]["profile"]["name"] == "panel-v1"
    assert caps["alleleCoverage"]["profile"] == {
        "name": "allele-panel-v1",
        "specificityRevision": "selected-sites-v2",
    }
    assert caps["alleleCoverage"]["sourceContractVersion"] == "3.3.0+lge.4"


def test_ignored_legacy_controls_and_tampered_saved_chemistry_reject():
    for field, value in [
        ("downsample_target", 0.5),
        ("in_memory_db", False),
        ("min_overlap", 12),
        ("typo_option", 1),
    ]:
        with pytest.raises(ValueError):
            Config(**kwargs(**{field: value}))
    with pytest.raises(ValueError):
        Config(candidate_profiles="union")
    config = Config(**kwargs())
    saved = config.to_dict()
    saved["primer_tm_min"] = 58.0
    with pytest.raises(ValueError):
        Config(**saved)


def test_default_ladder_can_be_bounded_to_two_requested_stages():
    resolved = options().from_config(
        Config(**kwargs(salvage="bounded", salvage_max_stages=2))
    )
    assert resolved.salvage_thresholds == (-28.0, -30.0)
    assert resolved.salvage_options().max_stages == 2
    with pytest.raises(ValueError):
        Config(
            **kwargs(
                salvage="bounded",
                salvage_max_stages=2,
                salvage_thresholds=(-28, -30, -32),
            )
        )


def test_new_cli_help_groups_and_salvage_flags_resolve(tmp_path):
    result, run = invoke(
        tmp_path,
        (
            "--salvage",
            "bounded",
            "--salvage-threshold",
            "-29",
            "--salvage-threshold",
            "-33",
            "--primary-tier",
            "salvage-2",
            "--work-cleanup-moves-per-round",
            "3",
        ),
    )
    assert result.exit_code == 0, result.output
    resolved = options().from_config(run.call_args.kwargs["config"])
    assert resolved.salvage_thresholds == (-29.0, -33.0)
    assert resolved.search_options().work_limits.cleanup_moves_per_round == 3
    helptext = CliRunner().invoke(app, ["panel-create", "--help"]).output
    assert (
        "Allele science" in helptext
        and "Allele compute" in helptext
        and "Allele salvage" in helptext
    )


@pytest.mark.parametrize(
    "flags",
    [
        ("--force",),
        ("--downsample-target", ".5"),
        ("--variant-selection", "full-cloud", "--subset-beam-width", "1"),
    ],
)
def test_cli_rejects_inactive_or_immutable_output_controls(tmp_path, flags):
    result, run = invoke(tmp_path, flags)
    assert result.exit_code == 2, result.output
    run.assert_not_called()


def test_reloaded_options_track_changed_resolved_values_as_new_requests():
    config = Config(**kwargs())
    saved = config.to_dict() | {"coverage_target": 0.8}
    rebuilt = Config(**saved)
    assert (
        options().from_config(rebuilt).to_dict()["requested_options"]["coverage_target"]
        == 0.8
    )


def test_preset_requires_explicit_bounds_and_rejects_unknown_preset():
    with pytest.raises(ValueError, match="bounds"):
        Config(selection_algorithm="allele-coverage")
    with pytest.raises(ValueError, match="preset"):
        Config(**kwargs(preset="other"))
    config = Config(
        selection_algorithm="allele-coverage", amplicon_size=200, amplicon_size_min=150
    )
    assert config.amplicon_size_max == 220


@pytest.mark.parametrize(
    "name",
    [
        "_primer_size_default_min",
        "_primer_size_default_max",
        "_primer_size_hgc_min",
        "_primer_size_hgc_max",
        "_primer_gc_default_min",
        "_primer_gc_default_max",
        "_primer_gc_hgc_min",
        "_primer_gc_hgc_max",
    ],
)
@pytest.mark.parametrize(
    "source", ["direct", "saved", "saved-request", "internal-request"]
)
def test_private_preset_overrides_reject_without_silent_provenance_loss(name, source):
    import json

    values = kwargs() if source == "direct" else Config(**kwargs()).to_dict()
    if source in ("direct", "saved"):
        values[name] = 100
    elif source == "saved-request":
        saved = json.loads(values["allele_options_json"])
        saved["requested_options"][name] = 100
        values["allele_options_json"] = json.dumps(saved)
    else:
        values = kwargs(_allele_requested_options={name: 100})
    with pytest.raises(ValueError, match="private.*" + name):
        Config(**values)


def test_private_keys_reject_before_none_filter_and_saved_resolved_filter():
    import json

    with pytest.raises(ValueError, match="private"):
        Config(**kwargs(_primer_gc_hgc_max=None))
    values = Config(**kwargs()).to_dict()
    saved = json.loads(values["allele_options_json"])
    saved["_primer_gc_hgc_max"] = 100
    values["allele_options_json"] = json.dumps(saved)
    with pytest.raises(ValueError, match="private"):
        Config(**values)


def test_legitimate_internal_requests_and_worker_metadata_preserve_complete_presets():
    from primalscheme3.panel.coverage_discovery import discovery_profiles

    config = Config(**kwargs(_allele_requested_options=kwargs()))
    saved = config.to_dict() | {
        "discovery_workers_by_target_profile": {
            "target-1": {"normal": 1, "high-gc": 2}
        },
    }
    restored = Config(**saved)
    for current in (config, restored):
        profiles = discovery_profiles(current)
        assert (
            profiles["normal"].primer_size_min,
            profiles["normal"].primer_size_max,
            profiles["normal"].primer_gc_min,
            profiles["normal"].primer_gc_max,
        ) == (19, 36, 30, 55)
        assert (
            profiles["high-gc"].primer_size_min,
            profiles["high-gc"].primer_size_max,
            profiles["high-gc"].primer_gc_min,
            profiles["high-gc"].primer_gc_max,
        ) == (17, 30, 40, 65)
    assert (
        options().from_config(restored).to_dict()
        == options().from_config(config).to_dict()
    )
    # Execution counts are recomputed, not trusted from caller metadata.
    assert restored.discovery_workers_by_target_profile == {}
