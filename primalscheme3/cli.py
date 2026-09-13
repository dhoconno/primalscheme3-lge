#!/usr/bin/python3
import argparse
import json
import pathlib
import sys
from importlib.metadata import version
from typing import Annotated

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
        PanelRunModes,
        typer.Option(
            help="Select what run mode",
        ),
    ] = PanelRunModes.REGION_ONLY.value,  # type: ignore
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
        str, typer.Option(help="Panel selector: legacy or coverage")
    ] = "legacy",
    coverage_metric: Annotated[
        str, typer.Option(help="Coverage interval: full-span or primer-trimmed")
    ] = "full-span",
    coverage_target: Annotated[
        float, typer.Option(help="Per-target coverage objective between zero and one")
    ] = 0.90,
    optimizer_seed: Annotated[
        int, typer.Option(help="Deterministic coverage optimizer seed")
    ] = 0,
    optimizer_starts: Annotated[
        int, typer.Option(help="Positive number of bounded optimizer starts")
    ] = 4,
    optimizer_repair_rounds: Annotated[
        int, typer.Option(help="Nonnegative repair rounds per optimizer start")
    ] = 2,
    optimizer_time_limit: Annotated[
        float, typer.Option(help="Positive finite selector wall-time budget in seconds")
    ] = 120.0,
    mispriming_product_size: Annotated[
        int | None,
        typer.Option(
            help="Coverage supplied-MSA product bound; defaults to 2000 and must be positive"
        ),
    ] = None,
):
    """
    Creates a primer panel
    """
    # Update the config with CLI params
    params = locals().copy()
    if selection_algorithm not in {"legacy", "coverage"}:
        raise typer.BadParameter("--selection-algorithm must be legacy or coverage")
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
                + ", ".join("--" + name.replace("_", "-") for name in changed)
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
        if mode != PanelRunModes.EQUAL:
            raise typer.BadParameter("coverage selection supports only --mode equal")
        if max_amplicons_region_group is not None:
            raise typer.BadParameter(
                "coverage selection does not support --max-amplicons-region-group"
            )
    params.pop("mispriming_product_size", None)
    config = create_config_with_amplicon_bounds(params)

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


if __name__ == "__main__":
    app()
