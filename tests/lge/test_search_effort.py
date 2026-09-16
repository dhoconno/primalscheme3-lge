"""Versioned compute defaults are independent of scientific policy and overrides."""

import json

import pytest
from test_allele_options import invoke, kwargs
from typer.testing import CliRunner

from primalscheme3.cli import app
from primalscheme3.core.config import Config
from primalscheme3.panel.allele_options import AlleleOptions

STANDARD = {
    "optimizer_time_limit": 120.0,
    "optimizer_starts": 4,
    "optimizer_repair_rounds": 2,
    "work_construction_candidate_attempts": 2048,
    "work_families_per_refresh": 16,
}
QUALITY = {
    "optimizer_time_limit": 3600.0,
    "optimizer_starts": 8,
    "optimizer_repair_rounds": 3,
    "work_construction_candidate_attempts": 8192,
    "work_families_per_refresh": 32,
}


@pytest.mark.parametrize(
    "effort,expected",
    [(None, STANDARD), ("standard-v1", STANDARD), ("quality-v1", QUALITY)],
)
def test_programmatic_defaults_and_roundtrip(effort, expected):
    c = Config(**kwargs(**({"search_effort": effort} if effort else {})))
    o = AlleleOptions.from_config(c)
    assert o.search_effort == (effort or "standard-v1")
    assert {k: getattr(o, k) for k in expected} == expected
    assert Config(**c.to_dict()).to_dict() == c.to_dict()
    assert ("search_effort" in o.to_dict()["requested_options"]) == (effort is not None)
    assert o.preset == "allele-balanced-v1" and o.phase_scheduling == "serial"
    assert o.salvage == "off" and o.salvage_time_limit == 60
    assert (
        o.intended_product_policy == "exact-supported"
        and o.secondary_product_policy == "ordered-disjoint-intended-sites"
    )
    search = o.search_options()
    assert (search.time_limit, search.starts, search.repair_rounds) == tuple(
        expected[k]
        for k in ("optimizer_time_limit", "optimizer_starts", "optimizer_repair_rounds")
    )
    assert (
        search.work_limits.construction_candidate_attempts
        == expected["work_construction_candidate_attempts"]
    )
    assert (
        search.work_limits.families_per_refresh == expected["work_families_per_refresh"]
    )


@pytest.mark.parametrize("key,value", list(STANDARD.items()))
def test_explicit_old_default_wins_programmatic_and_cli(tmp_path, key, value):
    c = Config(**kwargs(search_effort="quality-v1", **{key: value}))
    assert getattr(c, key) == value
    flags = ("--search-effort", "quality-v1", "--" + key.replace("_", "-"), str(value))
    for extra in (flags, flags[2:] + flags[:2]):
        result, run = invoke(tmp_path, extra)
        assert result.exit_code == 0, result.output + repr(result.exception)
        o = AlleleOptions.from_config(run.call_args.kwargs["config"])
        assert getattr(o, key) == value
        assert o.to_dict()["requested_options"][key] == value
        assert {k: getattr(o, k) for k in QUALITY if k != key} == {
            k: v for k, v in QUALITY.items() if k != key
        }


def test_click_default_parameters_do_not_override_effort(tmp_path):
    result, run = invoke(tmp_path, ("--search-effort", "quality-v1"))
    assert result.exit_code == 0, result.output + repr(result.exception)
    o = AlleleOptions.from_config(run.call_args.kwargs["config"])
    assert {k: getattr(o, k) for k in QUALITY} == QUALITY
    assert not set(QUALITY) & o.to_dict()["requested_options"].keys()


def test_saved_values_are_not_recomputed_and_old_missing_is_standard():
    c = Config(
        **kwargs(
            search_effort="quality-v1",
            optimizer_time_limit=9,
            optimizer_starts=4,
            work_families_per_refresh=3,
        )
    )
    saved = c.to_dict()
    assert Config(**saved).to_dict() == saved
    # A saved snapshot is already resolved: even selecting another effort label
    # cannot silently re-resolve its serialized numeric controls.
    r = Config(**(saved | {"search_effort": "standard-v1"}))
    assert (
        r.optimizer_time_limit,
        r.optimizer_starts,
        r.work_families_per_refresh,
    ) == (9, 4, 3)
    old = Config(**kwargs()).to_dict()
    old.pop("search_effort")
    data = json.loads(old["allele_options_json"])
    data.pop("search_effort")
    old["allele_options_json"] = json.dumps(data)
    restored = Config(**old)
    assert restored.search_effort == "standard-v1"
    assert {k: getattr(restored, k) for k in STANDARD} == STANDARD


def test_capability_descriptor_and_independence():
    from primalscheme3.panel.allele_catalog_cache import discovery_settings
    from primalscheme3.panel.allele_validation import AlleleConstraintProfile
    from primalscheme3.panel.coverage_discovery import discovery_profiles
    from primalscheme3.panel.coverage_provenance import capabilities_document

    cap = capabilities_document()["alleleCoverage"]["searchEfforts"]
    assert cap == {
        "default": "standard-v1",
        "policies": {
            "standard-v1": {"id": "standard-v1", "defaults": STANDARD},
            "quality-v1": {"id": "quality-v1", "defaults": QUALITY},
        },
    }
    cap["policies"]["quality-v1"]["defaults"]["optimizer_starts"] = 100
    assert (
        capabilities_document()["alleleCoverage"]["searchEfforts"]["policies"][
            "quality-v1"
        ]["defaults"]
        == QUALITY
    )
    a = Config(**kwargs())
    b = Config(**kwargs(search_effort="quality-v1"))
    assert discovery_settings(a, discovery_profiles(a)) == discovery_settings(
        b, discovery_profiles(b)
    )
    assert AlleleConstraintProfile.from_config(
        a
    ) == AlleleConstraintProfile.from_config(b)


@pytest.mark.parametrize("effort", ["unknown", "quality-v2", ""])
def test_unknown_effort_rejected(tmp_path, effort):
    with pytest.raises(ValueError):
        Config(**kwargs(search_effort=effort))
    result, run = invoke(tmp_path, ("--search-effort", effort))
    assert result.exit_code == 2
    run.assert_not_called()


@pytest.mark.parametrize("algorithm", ["legacy", "coverage"])
def test_nonallele_explicit_effort_rejected(tmp_path, algorithm):
    with pytest.raises(ValueError):
        Config(selection_algorithm=algorithm, search_effort="standard-v1")
    path = tmp_path / "a.fa"
    path.write_text(">a\nAAAA\n")
    result = CliRunner().invoke(
        app,
        [
            "panel-create",
            "--msa",
            str(path),
            "--output",
            str(tmp_path / "out"),
            "--selection-algorithm",
            algorithm,
            "--search-effort",
            "standard-v1",
        ],
    )
    assert result.exit_code == 2 and "allele" in result.output


def test_direct_cli_function_explicit_historical_default_is_not_absent(tmp_path):
    from unittest.mock import patch

    from primalscheme3.cli import panel_create

    p = tmp_path / "a.fa"
    p.write_text(">a\nAAAA\n")
    with patch("primalscheme3.cli.panelcreate") as run:
        panel_create(
            msa=[p],
            output=tmp_path / "out",
            **kwargs(search_effort="quality-v1", optimizer_starts=4),
        )
    o = AlleleOptions.from_config(run.call_args.kwargs["config"])
    assert o.optimizer_starts == 4 and o.optimizer_time_limit == 3600
    assert o.to_dict()["requested_options"]["optimizer_starts"] == 4


def test_click_default_map_is_an_explicit_override(tmp_path):
    from unittest.mock import patch

    p = tmp_path / "a.fa"
    p.write_text(">a\nAAAA\n")
    argv = [
        "panel-create",
        "--msa",
        str(p),
        "--output",
        str(tmp_path / "out"),
        "--selection-algorithm",
        "allele-coverage",
        "--amplicon-size-min",
        "150",
        "--amplicon-size-max",
        "450",
        "--search-effort",
        "quality-v1",
    ]
    with patch("primalscheme3.cli.panelcreate") as run:
        r = CliRunner().invoke(
            app, argv, default_map={"panel-create": {"optimizer_starts": 4}}
        )
    assert r.exit_code == 0, r.output + repr(r.exception)
    o = AlleleOptions.from_config(run.call_args.kwargs["config"])
    assert o.optimizer_starts == 4 and o.optimizer_time_limit == 3600
    assert o.to_dict()["requested_options"]["optimizer_starts"] == 4


def test_quality_and_explicit_controls_have_identical_native_history():
    from test_allele_validation import dna, fixture

    from primalscheme3.panel.allele_search import search_allele_assignments
    from primalscheme3.panel.allele_validation import AlleleConstraintProfile
    from primalscheme3.panel.coverage_history import CoverageHistory

    cat, _, ledger = fixture(((0, 200), (300, 500)), seq=dna(seed=0))
    a = Config(**kwargs(search_effort="quality-v1"))
    b = Config(**kwargs(**QUALITY))
    ao, bo = map(AlleleOptions.from_config, (a, b))
    assert ao.search_options() == bo.search_options()
    results = []
    for c, o in ((a, ao), (b, bo)):
        results.append(
            search_allele_assignments(
                cat,
                AlleleConstraintProfile.from_config(c),
                o.search_options(),
                ledger,
                clock=lambda: 0,
                history=CoverageHistory(run_id="effort-equivalence"),
            )
        )
    a, b = results
    assert a.assignments == b.assignments and a.assignments
    assert (
        a.coverage == b.coverage
        and a.validation == b.validation
        and a.metadata == b.metadata
    )
    for stream in ("evidence", "assessments", "events", "snapshots"):
        assert getattr(a.history, stream) == getattr(b.history, stream)


def test_real_cli_effort_is_retained_in_native_receipt_and_options(tmp_path):
    from primalscheme3.panel.allele_inspection import audit_allele_bundle

    raw = tmp_path / "short.fa"
    raw.write_text(">short\n" + "A" * 30 + "\n")
    out = tmp_path / "panel"
    args = [
        "panel-create",
        "--msa",
        str(raw),
        "--output",
        str(out),
        "--selection-algorithm",
        "allele-coverage",
        "--amplicon-size",
        "200",
        "--amplicon-size-min",
        "150",
        "--amplicon-size-max",
        "250",
        "--search-effort",
        "quality-v1",
        "--optimizer-starts",
        "4",
        "--offline-plots",
    ]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output + repr(result.exception)
    config = json.loads((out / "config.json").read_text())
    o = json.loads(config["allele_options_json"])
    assert (
        o["search_effort"] == "quality-v1"
        and o["optimizer_starts"] == 4
        and o["optimizer_time_limit"] == 3600
    )
    assert o["requested_options"]["search_effort"] == "quality-v1"
    assert (
        o["requested_options"]["optimizer_starts"] == 4
        and "optimizer_time_limit" not in o["requested_options"]
    )
    receipt = json.loads((out / "panel-provenance.json").read_text())
    assert receipt["status"] == "success" and receipt["exitStatus"] == 0
    assert receipt["resolvedOptions"]["search_effort"] == "quality-v1"
    optimizer = json.loads((out / "panel-optimizer.json").read_text())
    assert optimizer["options"]["search_effort"] == "quality-v1"
    assert "search_effort" not in optimizer["profile"]
    assert audit_allele_bundle(out)["valid"]


@pytest.mark.parametrize("value", [[], {}, True, 1])
def test_effort_requires_a_known_string(value):
    with pytest.raises(ValueError, match="search_effort"):
        Config(**kwargs(search_effort=value))
