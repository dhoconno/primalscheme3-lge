#!/usr/bin/python3
import argparse
import inspect
import json
import pathlib
import sys
from importlib.metadata import version
from typing import Annotated

import click
import typer

# Module imports
from primalscheme3.core.config import (
    AmpliconSizeMetric,
    Config,
    MappingType,
    TerminalGapPolicy,
)
from primalscheme3.core.downsample import downsample_scheme
from primalscheme3.core.msa import parse_msa
from primalscheme3.core.primer_visual import bedfile_plot_html, primer_mismatch_heatmap
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.interaction.interaction import (
    visualise_interactions,
)
from primalscheme3.panel.coverage_provenance import capabilities_document
from primalscheme3.panel.legacy_salvage_cli import (
    LEGACY_SALVAGE_OPTION_NAMES,
    resolve_legacy_salvage_options,
)
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate
from primalscheme3.repair.repair import repair
from primalscheme3.replace.replace import ReplaceRunModes, replace

# Import main functions
from primalscheme3.scheme.scheme_main import schemecreate

## Commands are in the format of
# {pclass}-{mode}
# pclass = panel or scheme

# Example to create a scheme
# scheme-create

# To repair a scheme
# scheme-repair

# To create a panel
# panel-create


def check_path_is_file(value: str | pathlib.Path) -> pathlib.Path:
    if isinstance(value, str):
        value = pathlib.Path(value)
    if not value.is_file():
        raise argparse.ArgumentTypeError(f"No file found at: '{str(value.absolute())}'")
    return value


# Create the main app
app = typer.Typer(name="primalscheme3", no_args_is_help=True)


def check_output_dir(output: pathlib.Path, force: bool):
    if output.exists() and not force:
        raise typer.BadParameter(
            f"--output '{output}' directory already exists. Use --force to overwrite"
        )


def create_config_with_amplicon_bounds(params: dict) -> Config:
    """Only explicit creation flags opt in; old saved bounds retain their metric."""
    metric = (
        AmpliconSizeMetric.REFERENCE_SPAN
        if any(
            params.get(key) is not None
            for key in ("amplicon_size_min", "amplicon_size_max")
        )
        else AmpliconSizeMetric.LEGACY_PAIRING
    )
    if metric == AmpliconSizeMetric.REFERENCE_SPAN:
        if params.get("bedfile") is not None:
            raise typer.BadParameter(
                "reference-span amplicon bounds do not support imported primer pairs (--bedfile)"
            )
        if (
            params.get("region_bedfile") is not None
            or params.get("mode") == PanelRunModes.REGION_ONLY
        ):
            raise typer.BadParameter(
                "reference-span amplicon bounds require a whole-MSA equal or entropy panel without regions"
            )
    try:
        return Config(**params, amplicon_size_metric=metric)
    except ValueError as error:
        raise typer.BadParameter(str(error)) from error


def typer_callback_version(value: bool):
    if value:
        version_str = typer.style(
            version("primalscheme3"), fg=typer.colors.GREEN, bold=True
        )
        typer.echo("PrimalScheme3-LGE version: " + version_str)
        raise typer.Exit()


def typer_callback_capabilities(value: bool):
    if value:
        typer.echo(
            json.dumps(capabilities_document(), sort_keys=True, separators=(",", ":"))
        )
        raise typer.Exit()


@app.callback()
def primalscheme3(
    value: Annotated[bool, typer.Option] = typer.Option(
        False, "--version", callback=typer_callback_version
    ),
    capabilities_json: Annotated[bool, typer.Option] = typer.Option(
        False,
        "--capabilities-json",
        callback=typer_callback_capabilities,
        is_eager=True,
        help="Emit the machine-readable engine capability and runtime identity.",
    ),
):
    pass


@app.command(no_args_is_help=True)
def scheme_create(
    msa: Annotated[
        list[pathlib.Path],
        typer.Option(
            help="The MSA to design against. To use multiple MSAs, use multiple --msa flags. (--msa 1.fasta --msa 2.fasta)",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(
            help="The output directory",
            resolve_path=True,
        ),
    ],
    amplicon_size: Annotated[
        int,
        typer.Option(
            help="Nominal amplicon size; omitted bounds resolve to ±10 percent. Does not rank by target proximity.",
            min=100,
            max=2000,
        ),
    ] = Config.amplicon_size,
    amplicon_size_min: Annotated[
        int | None,
        typer.Option(
            help="Inclusive minimum full reference BED span, including primers. Enables reference-span pairing.",
            min=1,
        ),
    ] = None,
    amplicon_size_max: Annotated[
        int | None,
        typer.Option(
            help="Inclusive maximum full reference BED span, including primers. Enables reference-span pairing.",
            min=1,
        ),
    ] = None,
    bedfile: Annotated[
        pathlib.Path | None,
        typer.Option(
            help="An existing bedfile to add primers to",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ] = None,
    min_overlap: Annotated[
        int,
        typer.Option(help="min amount of overlap between primers", min=0),
    ] = Config.min_overlap,
    n_pools: Annotated[
        int, typer.Option(help="Number of pools to use", min=1)
    ] = Config.n_pools,
    dimer_score: Annotated[
        float,
        typer.Option(
            help="Threshold for dimer interaction",
        ),
    ] = Config.dimer_score,
    terminal_gap_policy: Annotated[
        TerminalGapPolicy,
        typer.Option(
            help="LGE custom policy: legacy uses upstream Rust; observed-only uses Python discovery, excluding terminal missing coverage without imputation. Internal gaps remain observations."
        ),
    ] = TerminalGapPolicy.LEGACY,
    min_base_freq: Annotated[
        float,
        typer.Option(help="Min freq to be included,[0<=x<=1]", min=0.0, max=1.0),
    ] = Config.min_base_freq,
    mapping: Annotated[
        MappingType,
        typer.Option(
            help="How should the primers in the bedfile be mapped",
        ),
    ] = Config.mapping.value,  # type: ignore
    circular: Annotated[
        bool, typer.Option(help="Should a circular amplicon be added")
    ] = Config.circular,
    backtrack: Annotated[
        bool, typer.Option(help="Should the algorithm backtrack")
    ] = Config.backtrack,
    ignore_n: Annotated[
        bool,
        typer.Option(help="Should N in the input genomes be ignored"),
    ] = Config.ignore_n,
    force: Annotated[
        bool, typer.Option(help="Override the output directory")
    ] = Config.force,
    input_bedfile: Annotated[
        pathlib.Path | None,
        typer.Option(
            help="Path to a primer.bedfile containing the pre-calculated primers"
        ),
    ] = Config.input_bedfile,
    high_gc: Annotated[bool, typer.Option(help="Use high GC primers")] = Config.high_gc,
    offline_plots: Annotated[
        bool,
        typer.Option(
            "--offline-plots/--online-plots",
            help="Offline plots includes 3Mb of dependencies, so they can be viewed offline",
        ),
    ] = True,
    use_matchdb: Annotated[
        bool,
        typer.Option(
            help="Create and use a mispriming database",
        ),
    ] = Config.use_matchdb,
    ncores: Annotated[
        int,
        typer.Option(help="Number of CPU cores to use during digestion", min=1),
    ] = Config.ncores,
    use_annealing: Annotated[
        bool,
        typer.Option(
            "--use-annealing/--use-tm",
            help="Using annealing proportion rather than Tm to calculate primers",
        ),
    ] = Config.use_annealing,
    # Downsample params
    downsample: Annotated[
        bool,
        typer.Option(
            help="EXPERIMENTAL: Reduce number of primers in a cloud by calculating inter-primercloud annealing",
            hidden=True,
        ),
    ] = Config.downsample,
    downsample_target: Annotated[
        float,
        typer.Option(
            help="EXPERIMENTAL: Ensure X proportion of primers have >= annealing (if using downsampling)",
            hidden=True,
        ),
    ] = Config.downsample_target,
):
    """
    Creates a tiling overlap scheme for each MSA file
    """
    # Update the config with CLI params
    config = create_config_with_amplicon_bounds(locals())

    # Check the output directory
    check_output_dir(output, force)

    # Set up the progress manager
    pm = ProgressManager()
    schemecreate(
        msa=msa,
        output_dir=output,
        config=config,
        pm=pm,
        force=force,
        input_bedfile=input_bedfile,
        offline_plots=offline_plots,
    )


@app.command(no_args_is_help=True)
def replace_primerpair(
    primername: Annotated[
        str, typer.Argument(help="The name of the primer to replace")
    ],
    primerbed: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primer to replace",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    msa: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The msa used to create the original primer scheme",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    amplicon_size_max: Annotated[
        int,
        typer.Option(
            help="The max size of an amplicon [100<=x<=2000]", min=100, max=2000
        ),
    ],
    amplicon_size_min: Annotated[
        int,
        typer.Option(
            help="The min size of an amplicon. [100<=x<=2000]", min=100, max=2000
        ),
    ],
    config_path: Annotated[
        pathlib.Path,
        typer.Option(
            help="The config.json used to create the original primer scheme",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(
            help="The output directory",
            resolve_path=True,
        ),
    ],
    mode: Annotated[
        ReplaceRunModes,
        typer.Option(
            help="Select what run mode",
        ),
    ] = ReplaceRunModes.ListAll.value,  # type: ignore
    force: Annotated[
        bool, typer.Option(help="Override the output directory")
    ] = Config.force,
    mask_old_sites: Annotated[
        bool,
        typer.Option(
            help="If True prevents replacement primers from spanning old primer regions."
        ),
    ] = True,
):
    """
    Replaces a primerpair in a bedfile
    """

    # Read in the config file
    with open(config_path) as file:
        _cfg: dict = json.load(file)
    config = Config(**_cfg)

    # Update the config with CLI params
    config.assign_kwargs(**locals())

    # Update the config with required
    config.in_memory_db = True  # Set up the db

    # Check the output directory
    check_output_dir(output, force)

    # Set up the progress manager
    pm = ProgressManager()

    replace(
        primername=primername,
        config=config,
        primerbed=primerbed,
        msapath=msa,
        pm=pm,
        output=output,
        force=force,
        mode=mode,
        mask_old_sites=mask_old_sites,
    )


@app.command(no_args_is_help=True)
def panel_create(
    msa: Annotated[
        list[pathlib.Path],
        typer.Option(
            help="Paths to the MSA files", exists=True, readable=True, resolve_path=True
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(
            help="The output directory",
            resolve_path=True,
        ),
    ],
    region_bedfile: Annotated[
        pathlib.Path | None,
        typer.Option(
            help="Path to the bedfile containing the wanted regions",
            readable=True,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ] = None,
    input_bedfile: Annotated[
        pathlib.Path | None,
        typer.Option(
            help="Path to a primer.bedfile containing the pre-calculated primers",
            readable=True,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ] = None,
    mode: Annotated[
        PanelRunModes | None,
        typer.Option(
            help="Select what run mode",
        ),
    ] = None,
    amplicon_size: Annotated[
        int,
        typer.Option(
            help="Nominal amplicon size; omitted bounds resolve to ±10 percent. Does not rank by target proximity.",
            min=100,
            max=2000,
        ),
    ] = Config.amplicon_size,
    amplicon_size_min: Annotated[
        int | None,
        typer.Option(
            help="Inclusive minimum full reference BED span, including primers. Enables reference-span pairing.",
            min=1,
        ),
    ] = None,
    amplicon_size_max: Annotated[
        int | None,
        typer.Option(
            help="Inclusive maximum full reference BED span, including primers. Enables reference-span pairing.",
            min=1,
        ),
    ] = None,
    n_pools: Annotated[
        int, typer.Option(help="Number of pools to use", min=1)
    ] = Config.n_pools,
    dimer_score: Annotated[
        float, typer.Option(help="Threshold for dimer interaction")
    ] = Config.dimer_score,
    terminal_gap_policy: Annotated[
        TerminalGapPolicy | None,
        typer.Option(
            help="LGE custom policy: legacy uses upstream Rust; observed-only uses Python discovery, excluding terminal missing coverage without imputation. Internal gaps remain observations."
        ),
    ] = None,
    min_base_freq: Annotated[
        float,
        typer.Option(help="Min freq to be included,[0<=x<=1]", min=0.0, max=1.0),
    ] = Config.min_base_freq,
    mapping: Annotated[
        MappingType,
        typer.Option(
            help="How should the primers in the bedfile be mapped",
        ),
    ] = Config.mapping.value,  # type: ignore
    max_amplicons: Annotated[
        int | None, typer.Option(help="Max number of amplicons to create", min=0)
    ] = None,
    max_amplicons_msa: Annotated[
        int | None, typer.Option(help="Max number of amplicons for each MSA", min=0)
    ] = None,
    max_amplicons_region_group: Annotated[
        int | None,
        typer.Option(help="Max number of amplicons for each region", min=1),
    ] = None,
    force: Annotated[bool, typer.Option(help="Override the output directory")] = False,
    high_gc: Annotated[bool, typer.Option(help="Use high GC primers")] = Config.high_gc,
    offline_plots: Annotated[
        bool,
        typer.Option(
            "--offline-plots/--online-plots",
            help="Offline plots includes 3Mb of dependencies, so they can be viewed offline",
        ),
    ] = True,
    use_matchdb: Annotated[
        bool,
        typer.Option(
            help="Create and use a mispriming database",
        ),
    ] = Config.use_matchdb,
    ncores: Annotated[
        int,
        typer.Option(help="Number of CPU cores to use during digestion", min=1),
    ] = Config.ncores,
    use_annealing: Annotated[
        bool,
        typer.Option(
            "--use-annealing/--use-tm",
            help="Using annealing proportion rather than Tm to calculate primers",
        ),
    ] = Config.use_annealing,
    # Downsample params
    downsample: Annotated[
        bool,
        typer.Option(
            help="EXPERIMENTAL: Reduce number of primers in a cloud by calculating inter-primercloud annealing",
            hidden=True,
        ),
    ] = Config.downsample,
    downsample_target: Annotated[
        float,
        typer.Option(
            help="EXPERIMENTAL: Ensure X proportion of primers have >= annealing (if using downsampling)",
            hidden=True,
        ),
    ] = Config.downsample_target,
    selection_algorithm: Annotated[
        str, typer.Option(help="Panel selector: legacy, coverage or allele-coverage")
    ] = "legacy",
    coverage_metric: Annotated[
        str | None,
        typer.Option(
            help="Metric: full-span, primer-trimmed, or observed-allele-primer-trimmed for allele mode"
        ),
    ] = None,
    coverage_target: Annotated[
        float | None,
        typer.Option(
            help="Per-target coverage objective; allele default .95, historical default .90"
        ),
    ] = None,
    optimizer_seed: Annotated[
        int, typer.Option(help="Deterministic coverage optimizer seed")
    ] = 0,
    optimizer_starts: Annotated[
        int | None,
        typer.Option(help="Positive number of bounded optimizer starts; standard default 4"),
    ] = None,
    optimizer_repair_rounds: Annotated[
        int | None,
        typer.Option(help="Nonnegative repair rounds per optimizer start; standard default 2"),
    ] = None,
    optimizer_time_limit: Annotated[
        float | None,
        typer.Option(help="Positive finite selector wall-time budget in seconds; standard default 120"),
    ] = None,
    mispriming_product_size: Annotated[
        int | None,
        typer.Option(
            help="Coverage supplied-MSA product bound; defaults to 2000 and must be positive"
        ),
    ] = None,
    preset: Annotated[
        str | None,
        typer.Option(
            help="Versioned allele preset: allele-balanced-v1",
            rich_help_panel="Allele science",
        ),
    ] = None,
    candidate_profiles: Annotated[
        str | None,
        typer.Option(
            help="Complete chemistry profiles: union, normal, high-gc",
            rich_help_panel="Allele science",
        ),
    ] = None,
    reuse_discovery: Annotated[
        pathlib.Path | None,
        typer.Option(
            exists=True, file_okay=False,
            help="Reuse a verified panel-cache artifact with identical discovery biology; copy all origin history into this output",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    search_effort: Annotated[
        str | None,
        typer.Option(
            help="Compute defaults: standard-v1 or quality-v1; explicit controls override. Science and scheduling unchanged.",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    phase_scheduling: Annotated[
        str | None,
        typer.Option(
            help="Phase scheduling: serial or reserved",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    variant_selection: Annotated[
        str | None,
        typer.Option(
            help="subsets or full-cloud controlled discovery ablation",
            rich_help_panel="Allele science",
        ),
    ] = None,
    allele_weighting: Annotated[
        str | None,
        typer.Option(
            help="Distinct observed allele classes: distinct-observed",
            rich_help_panel="Allele science",
        ),
    ] = None,
    discovery_length_mode: Annotated[
        str | None,
        typer.Option(
            help="first-compatible or all profile lengths per row/anchor",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    discovery_history: Annotated[
        str | None,
        typer.Option(
            help="Discovery history detail: compact (default) or full advanced diagnostics",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    specificity_terminal_k: Annotated[
        int | None,
        typer.Option(
            help="Common terminal seed length; allele default 17, one substitution, no indels",
            rich_help_panel="Allele science",
        ),
    ] = None,
    intended_product_policy: Annotated[
        str | None,
        typer.Option(
            help="exact-supported (default) or concrete-designated-sites; near-matches add no coverage",
            rich_help_panel="Allele science",
        ),
    ] = None,
    secondary_product_policy: Annotated[
        str | None,
        typer.Option(
            help="ordered-disjoint-intended-sites, ordered-disjoint-concrete-designated-sites/v1, or reject-secondary-products/v1",
            rich_help_panel="Allele science",
        ),
    ] = None,
    subset_beam_width: Annotated[
        int | None,
        typer.Option(
            help="Positive subset neighborhood beam width; default 16",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    subset_expansion_limit: Annotated[
        int | None,
        typer.Option(
            help="Positive deterministic subset materialization budget; default 256",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    exchange_width: Annotated[
        int | None,
        typer.Option(
            help="Maximum incumbent configurations per exchange, 0 through 2",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    salvage: Annotated[
        str | None,
        typer.Option(
            help="Experimental dimer salvage: off or bounded; default off",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    salvage_max_stages: Annotated[
        int | None,
        typer.Option(
            help="Maximum enabled salvage stages, 1 through 3",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    salvage_max_edges_per_pool: Annotated[
        int | None,
        typer.Option(
            help="Cumulative strict-violating physical edges per pool; default 8",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    salvage_max_oligos_per_pool: Annotated[
        int | None,
        typer.Option(
            help="Cumulative incident physical species per pool; default 4",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    salvage_time_limit: Annotated[
        float | None,
        typer.Option(
            help="Positive finite seconds per salvage tier; default 60",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    primary_tier: Annotated[
        str | None,
        typer.Option(
            help="strict or an enabled salvage-1, salvage-2, salvage-3 result",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    work_frontier_candidates: Annotated[
        int | None,
        typer.Option(
            help="Deterministic frontier candidates budget; default 64",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_construction_candidate_attempts: Annotated[
        int | None,
        typer.Option(
            help="Deterministic construction candidate attempts budget; default 2048",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_repair_candidate_probes_per_round: Annotated[
        int | None,
        typer.Option(
            help="Deterministic repair candidate probes per round budget; default 128",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_repair_neighborhoods_per_round: Annotated[
        int | None,
        typer.Option(
            help="Deterministic repair neighborhoods per round budget; default 256",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_repair_trials_per_round: Annotated[
        int | None,
        typer.Option(
            help="Deterministic repair trials per round budget; default 256",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_pool_lookahead_candidates: Annotated[
        int | None,
        typer.Option(
            help="Deterministic pool lookahead candidates budget; default 4",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_cleanup_moves_per_round: Annotated[
        int | None,
        typer.Option(
            help="Deterministic cleanup moves per round budget; default 64",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    work_families_per_refresh: Annotated[
        int | None,
        typer.Option(
            help="Deterministic families per refresh budget; default 16",
            rich_help_panel="Allele compute",
        ),
    ] = None,
    salvage_thresholds: Annotated[
        list[float] | None,
        typer.Option(
            "--salvage-threshold",
            help="Repeat a finite strictly decreasing cutoff below -26; default -28, -30, -32",
            rich_help_panel="Allele salvage",
        ),
    ] = None,
    legacy_salvage: Annotated[
        str | None,
        typer.Option(
            help="Legacy-only dimer salvage: off or bounded; default off",
            rich_help_panel="Legacy salvage",
        ),
    ] = None,
    legacy_salvage_thresholds: Annotated[
        list[float] | None,
        typer.Option(
            "--legacy-salvage-threshold",
            help="Repeat a finite decreasing legacy salvage cutoff; default -28, -30, -32",
            rich_help_panel="Legacy salvage",
        ),
    ] = None,
    legacy_salvage_floor: Annotated[
        float | None,
        typer.Option(help="Legacy salvage hard floor; default -32", rich_help_panel="Legacy salvage"),
    ] = None,
    legacy_salvage_max_edges_per_pool: Annotated[
        int | None,
        typer.Option(help="Cumulative relaxed conflict edges per pool; default 8", min=0, rich_help_panel="Legacy salvage"),
    ] = None,
    legacy_salvage_max_incident_species_per_pool: Annotated[
        int | None,
        typer.Option(help="Cumulative incident oligo species per pool; default 4", min=0, rich_help_panel="Legacy salvage"),
    ] = None,
    legacy_salvage_min_reference_gain: Annotated[
        int | None,
        typer.Option(help="Minimum new trimmed reference bases per addition; default 1", min=1, rich_help_panel="Legacy salvage"),
    ] = None,
    legacy_salvage_max_candidate_evaluations: Annotated[
        int | None,
        typer.Option(help="Cumulative candidate evaluations per pass; default 10000", min=1, rich_help_panel="Legacy salvage"),
    ] = None,
):
    """
    Creates a primer panel
    """
    # Parameter-source tracking keeps historical defaults out of allele requests.
    params = locals().copy()
    from primalscheme3.panel.allele_options import NEW_OPTION_NAMES

    context = click.get_current_context(silent=True)
    if context is not None:
        explicit = {
            name
            for name in params
            if context.get_parameter_source(name)
            not in (None, click.core.ParameterSource.DEFAULT)
        }
    else:
        defaults = inspect.signature(panel_create).parameters
        explicit = {
            name
            for name, value in params.items()
            if value is not None and value != defaults[name].default
        }
    if selection_algorithm not in {"legacy", "coverage", "allele-coverage"}:
        raise typer.BadParameter(
            "--selection-algorithm must be legacy, coverage or allele-coverage"
        )
    if selection_algorithm == "allele-coverage":
        requested = {
            name: params[name] for name in explicit if params[name] is not None
        }
        params = {name: value for name, value in params.items() if value is not None}
        if mispriming_product_size is not None:
            params["mismatch_product_size"] = mispriming_product_size
            requested["mismatch_product_size"] = requested.pop(
                "mispriming_product_size", mispriming_product_size
            )
        params.pop("mispriming_product_size", None)
        params.setdefault("mode", PanelRunModes.EQUAL)
        params["_allele_requested_options"] = requested
    else:
        for name, value in {
            "optimizer_starts": 4,
            "optimizer_repair_rounds": 2,
            "optimizer_time_limit": 120.0,
        }.items():
            if params[name] is None:
                params[name] = value
        supplied_new = sorted(
            name for name in NEW_OPTION_NAMES if params.get(name) is not None
        )
        if supplied_new:
            raise typer.BadParameter(
                "allele options require --selection-algorithm allele-coverage: "
                + ", ".join(supplied_new)
            )
        params["mode"] = mode if mode is not None else PanelRunModes.REGION_ONLY
        params["terminal_gap_policy"] = (
            terminal_gap_policy
            if terminal_gap_policy is not None
            else TerminalGapPolicy.LEGACY
        )
        params["coverage_metric"] = (
            coverage_metric if coverage_metric is not None else "full-span"
        )
        params["coverage_target"] = (
            coverage_target if coverage_target is not None else 0.90
        )
        if selection_algorithm == "legacy":
            defaults = {
                "coverage_metric": "full-span",
                "coverage_target": 0.90,
                "optimizer_seed": 0,
                "optimizer_starts": 4,
                "optimizer_repair_rounds": 2,
                "optimizer_time_limit": 120.0,
            }
            changed = [
                name for name, default in defaults.items() if params[name] != default
            ]
            if changed:
                raise typer.BadParameter(
                    "coverage optimizer options require --selection-algorithm coverage: "
                    + ", ".join(changed)
                )
            if mispriming_product_size not in {None, 0}:
                raise typer.BadParameter(
                    "nonzero --mispriming-product-size requires --selection-algorithm coverage"
                )
            if max_amplicons == 0 or max_amplicons_msa == 0:
                raise typer.BadParameter(
                    "legacy selection requires positive max-amplicons limits when supplied"
                )
            params["mismatch_product_size"] = 0
        else:
            params["mismatch_product_size"] = (
                2000 if mispriming_product_size is None else mispriming_product_size
            )
            if params["mode"] != PanelRunModes.EQUAL:
                raise typer.BadParameter(
                    "coverage selection supports only --mode equal"
                )
            if max_amplicons_region_group is not None:
                raise typer.BadParameter(
                    "coverage selection does not support --max-amplicons-region-group"
                )
        params.pop("mispriming_product_size", None)
    try:
        legacy_salvage_options = resolve_legacy_salvage_options(
            selection_algorithm=selection_algorithm,
            mode=params.get("mode", mode),
            mapping=mapping,
            region_bedfile=region_bedfile,
            input_bedfile=input_bedfile,
            config_input_bedfile=None,
            strict_cutoff=dimer_score,
            explicit=explicit & LEGACY_SALVAGE_OPTION_NAMES,
            legacy_salvage=legacy_salvage,
            legacy_salvage_thresholds=legacy_salvage_thresholds,
            legacy_salvage_floor=legacy_salvage_floor,
            legacy_salvage_max_edges_per_pool=legacy_salvage_max_edges_per_pool,
            legacy_salvage_max_incident_species_per_pool=legacy_salvage_max_incident_species_per_pool,
            legacy_salvage_min_reference_gain=legacy_salvage_min_reference_gain,
            legacy_salvage_max_candidate_evaluations=legacy_salvage_max_candidate_evaluations,
        )
    except ValueError as error:
        raise typer.BadParameter(str(error)) from error
    for name in LEGACY_SALVAGE_OPTION_NAMES:
        params.pop(name, None)
    config = create_config_with_amplicon_bounds(params)
    mode = PanelRunModes(params["mode"])

    # Check the output directory
    check_output_dir(output, force)

    # Set up the progress manager
    pm = ProgressManager()

    panelcreate(
        msa=msa,
        output_dir=output,
        region_bedfile=region_bedfile,
        input_bedfile=input_bedfile,
        mode=mode,
        config=config,
        pm=pm,
        max_amplicons=max_amplicons,
        max_amplicons_msa=max_amplicons_msa,
        max_amplicons_region_group=max_amplicons_region_group,
        force=force,
        offline_plots=offline_plots,
        executed_argv=list(sys.argv),
        legacy_salvage_options=legacy_salvage_options,
    )


@app.command(no_args_is_help=True)
def interactions(
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="Path to the bedfile",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    threshold: Annotated[
        float,
        typer.Option(
            help="Only show interactions more severe (Lower score) than this value",
        ),
    ] = -26.0,
):
    """
    Shows all the primer-primer interactions within a bedfile
    """
    visualise_interactions(bedfile, threshold)


@app.command(no_args_is_help=True)
def repair_mode(
    bedfile: Annotated[
        pathlib.Path,
        typer.Option(
            help="Path to the bedfile",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    msa: Annotated[
        pathlib.Path,
        typer.Option(
            help="An MSA, with the reference.fasta, aligned to any new genomes with mutations",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    config: Annotated[
        pathlib.Path,
        typer.Option(
            help="Path to the config.json",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path, typer.Option(help="The output directory", dir_okay=True)
    ],
    force: Annotated[bool, typer.Option(help="Override the output directory")] = False,
):
    """
    Repairs a primer scheme via adding more primers to account for new mutations
    """
    # Set up the progress manager
    pm = ProgressManager()

    repair(
        config_path=config,
        bedfile_path=bedfile,
        force=force,
        pm=pm,
        output_dir=output,
        msa_path=msa,
    )


@app.command(no_args_is_help=True)
def visualise_primer_mismatches(
    msa: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The MSA used to design the scheme",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(help="Output location of the plot", dir_okay=False, writable=True),
    ] = pathlib.Path("primer.html"),
    include_seqs: Annotated[
        bool,
        typer.Option(help="Reduces plot filesize, by excluding primer sequences"),
    ] = True,
    offline_plots: Annotated[
        bool,
        typer.Option(
            help="Includes 3Mb of dependencies into the plots, so they can be viewed offline"
        ),
    ] = True,
):
    """
    Visualise mismatches between primers and the input genomes
    """

    array, seqdict = parse_msa(msa)

    with open(output, "w") as outfile:
        outfile.write(
            primer_mismatch_heatmap(
                array=array,
                seqdict=seqdict,
                bedfile=bedfile,
                offline_plots=offline_plots,
                include_seqs=include_seqs,
            )
        )


@app.command(no_args_is_help=True)
def visualise_bedfile(
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    ref_id: Annotated[str, typer.Option(help="The reference genome ID")],
    ref_path: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(help="Output location of the plot", dir_okay=False, writable=True),
    ] = pathlib.Path("bedfile.html"),
):
    """
    Visualise the bedfile
    """
    import dnaio

    refs = []
    ref_genome: str | None = None

    # Find the wanted file
    with dnaio.open(ref_path, mode="r") as ref_file:
        for record in ref_file:
            refs.append(record.id)

            # See if wanted genome
            if record.id == ref_id:
                ref_genome = record.sequence

    if ref_genome is None:
        raise typer.BadParameter(
            f"Reference genome ID '{ref_id}' not found in '{ref_path}'. Options: {', '.join(refs)}"
        )

    with open(output, "w") as outfile:
        outfile.write(
            bedfile_plot_html(bedfile=bedfile, ref_name=ref_id, ref_seq=ref_genome)
        )


@app.command(no_args_is_help=True, hidden=True)
def downsample_existing_scheme(
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    downsample_target: Annotated[
        float,
        typer.Option(min=0, max=1),
    ] = Config.downsample_target,
    visualise: Annotated[
        bool,
        typer.Option(help="Will print graph visualisations of downsampling"),
    ] = False,
):
    """Will try and reduce number of primers in a cloud by calculating inter-primercloud annealing"""

    config = Config(**locals())
    config.downsample = True

    downsample_scheme(bedfile, config, visualise)


@app.command("panel-history")
def panel_history(
    bundle: Annotated[pathlib.Path, typer.Option(help="Source allele panel bundle")],
    output: Annotated[pathlib.Path, typer.Option(help="New result directory outside the bundle")],
    entity: str | None = None,
    target: Annotated[str | None, typer.Option(help="Target ID or reference name")] = None,
    region: Annotated[str | None, typer.Option(help="Zero-based half-open reference START:END")] = None,
    pool: Annotated[int | None, typer.Option(help="One-based pool number")] = None,
    stage: str | None = None,
    profile: str | None = None,
    lineage: bool = False,
    limit: int = 100,
    offset: int = 0,
):
    """Query generated, rejected and selected entities with bounded history pages."""
    from primalscheme3.panel.allele_inspection import run_inspection

    # Parsing belongs inside the receipt-producing workflow, including errors.
    result = run_inspection(
        "panel-history", bundle, output, argv=list(sys.argv), entity=entity,
        target=target, region=region, pool=pool, stage=stage, profile=profile,
        lineage=lineage, limit=limit, offset=offset,
    )
    typer.echo(str(output / "query.json"))
    if not result["valid"]:
        raise typer.Exit(1)


@app.command("panel-discovery-diagnose")
def panel_discovery_diagnose(
    bundle: Annotated[pathlib.Path, typer.Option(help="Source allele panel bundle")],
    output: Annotated[pathlib.Path, typer.Option(help="New replay result directory outside the bundle")],
    family_id: Annotated[str | None, typer.Option(help="Candidate family ID to replay")]=None,
    site_id: Annotated[str | None, typer.Option(help="Concrete site ID to replay")]=None,
):
    """Reconstruct detailed discovery at one stored site or family anchor."""
    from primalscheme3.panel.discovery_diagnostics import diagnose_discovery

    result = diagnose_discovery(
        bundle=bundle,
        output=output,
        family_id=family_id,
        site_id=site_id,
        argv=list(sys.argv),
    )
    typer.echo(str(output / "discovery-diagnostic.json"))
    if not result.get("valid", False):
        raise typer.Exit(1)


@app.command("panel-audit")
def panel_audit(
    bundle: Annotated[pathlib.Path, typer.Option(help="Source allele panel bundle")],
    output: Annotated[pathlib.Path, typer.Option(help="New audit directory outside the bundle")],
    tier: Annotated[str | None, typer.Option(help="Audit a named tier; default audits all published tiers")] = None,
):
    """Reparse stored original inputs and freshly audit scientific stage artifacts."""
    from primalscheme3.panel.allele_inspection import run_inspection

    result = run_inspection("panel-audit", bundle, output, argv=list(sys.argv), tier=tier)
    typer.echo(str(output / "validation.json"))
    if not result["valid"]:
        raise typer.Exit(1)


@app.command("panel-cache")
def panel_cache(
    bundle: Annotated[pathlib.Path, typer.Option(help="Completed allele panel whose discovery origins should be preserved")],
    output: Annotated[pathlib.Path, typer.Option(help="New portable discovery cache directory outside the bundle")],
):
    """Export verified discovery for repeated selector experiments without regeneration."""
    from primalscheme3.panel.allele_catalog_cache import export_panel_discovery_cache

    try:
        export_panel_discovery_cache(bundle, output, argv=list(sys.argv))
    except (Exception, KeyboardInterrupt) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    typer.echo(str(output / "manifest.json"))


if __name__ == "__main__":
    app()
