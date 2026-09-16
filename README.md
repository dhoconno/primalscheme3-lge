# PrimalScheme3-LGE — custom fork

**This is a separately maintained LGE custom fork, not an upstream PrimalScheme3 release.** It is based on [artic-network/primalscheme3](https://github.com/artic-network/primalscheme3) v3.3.0 and retains upstream authorship and GPL-3.0 licensing.

The development distribution version is **3.3.0+lge.4**; these development changes have not been released. Its executable remains `primalscheme3`, but `--version` identifies **PrimalScheme3-LGE**. Obtain published custom wheels and source from [GitHub releases](https://github.com/dhoconno/primalscheme3-lge/releases); installing the ordinary PyPI `primalscheme3` package does not select this fork.

## Development: allele-aware multi-MSA panels

`panel-create --selection-algorithm allele-coverage` combines normal and high-GC candidate discovery, selects compatible primer variants within an amplicon, and optimizes coverage after primer trimming with equal weight per distinct observed allele. The `allele-balanced-v1` preset targets 95% per MSA while retaining useful partial coverage. Strict dimer screening remains −26; optional bounded salvage exports separately audited alternatives and retains strict as primary by default.

The CLI exposes candidate profile, geometry, pool, specificity, subset-search, work-budget and salvage controls. Allele-aware discovery history defaults to compact summaries; `--discovery-history full` retains detailed discovery evidence for advanced comparison and debugging. `panel-history`, `panel-audit`, `panel-cache` and `panel-discovery-diagnose` query or reconstruct bounded details with exact reproducibility receipts. See [commands, defaults, interpretation and advanced controls](docs/allele-coverage.md). Legacy and prior coverage modes remain available.

Creation commands accept explicit `--amplicon-size-min` and `--amplicon-size-max` bounds on the full reference-coordinate BED span, including both primers. For example, target `200`, minimum `150`, and maximum `250` permit reported spans from 150 through 250 bases. See [the size contract and supported workflows](docs/lge-fork.md#explicit-amplicon-size-bounds).

The fork adds an explicit `--terminal-gap-policy` to `scheme-create` and `panel-create`:

- `legacy` preserves upstream Rust discovery and is the standalone CLI default.
- `observed-only` uses the patched Python discovery path, treating missing terminal observations as unavailable data while retaining observed internal sequence information.

The observed-only path supports **multiple worker processes** through `--ncores`. Work is bounded by requested workers, available CPUs and work chunks, and results are reconstructed in stable position order. Each worker has an alignment copy, so more workers use more memory. Small inputs can be dominated by process startup and downstream work.

Native configuration records the policy, discovery backend and effective worker counts. This is a backend change as well as a missing-data change; comparisons must not be described as changing only a denominator. Missing observations are not evidence that an allele matches a primer.

See [custom behavior, installation and validation](docs/lge-fork.md), [reproducible build script](scripts/build_lge_wheel.py), and [license](LICENSE). The upstream documentation below is retained for reference; its ordinary PyPI/source installation instructions refer to upstream unless the custom GitHub release is selected explicitly.

---

# Primalscheme3

[![CI](https://github.com/ChrisgKent/primalscheme3/actions/workflows/pytest.yml/badge.svg)](https://github.com/ChrisgKent/primalscheme3/actions/workflows/pytest.yml) [![Generic badge](https://img.shields.io/badge/DOI-10.1101/2024.12.20.629611-blue.svg)](https://doi.org/10.1101/2024.12.20.629611)

This is a command-line interface tool that generates a primer scheme from a Multiple Sequence Alignment (MSA) file, utilising degenerate primers to handle variation in the genomes.

## Installation

Via pip!
```
pip install primalscheme3
```

From source 
```         
git clone https://github.com/artic-network/primalscheme3
cd primalscheme3
uv sync
uv run primalscheme3
```

# `PrimalScheme3`

**Usage**:

```console
$ primalscheme3 [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--version`
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `interactions`: Shows all the primer-primer interactions...
* `panel-create`: Creates a primer panel
* `repair-mode`: Repairs a primer scheme via adding more...
* `scheme-create`: Creates a tiling overlap scheme for each...
* `scheme-replace`: Replaces a primerpair in a bedfile
* `visualise-bedfile`: Visualise the bedfile
* `visualise-primer-mismatches`: Visualise mismatches between primers and...

## `primalscheme3 interactions`

Shows all the primer-primer interactions within a bedfile

**Usage**:

```console
$ primalscheme3 interactions [OPTIONS] BEDFILE
```

**Arguments**:

* `BEDFILE`: Path to the bedfile  [required]

**Options**:

* `--threshold FLOAT`: Only show interactions more severe (Lower score) than this value  [default: -26.0]
* `--help`: Show this message and exit.

## `primalscheme3 panel-create`

Creates a primer panel

**Usage**:

```console
$ primalscheme3 panel-create [OPTIONS]
```

**Options**:

* `--msa PATH`: Paths to the MSA files  [required]
* `--output PATH`: The output directory  [required]
* `--region-bedfile FILE`: Path to the bedfile containing the wanted regions
* `--input-bedfile FILE`: Path to a primer.bedfile containing the pre-calculated primers
* `--mode [entropy|region-only|equal]`: Select what run mode  [default: region-only]
* `--amplicon-size INTEGER`: The size of an amplicon  [default: 400]
* `--n-pools INTEGER RANGE`: Number of pools to use  [default: 2; x>=1]
* `--dimer-score FLOAT`: Threshold for dimer interaction  [default: -26.0]
* `--min-base-freq FLOAT RANGE`: Min freq to be included,[0<=x<=1]  [default: 0.0; 0.0<=x<=1.0]
* `--mapping [first|consensus]`: How should the primers in the bedfile be mapped  [default: first]
* `--max-amplicons INTEGER RANGE`: Max number of amplicons to create  [x>=1]
* `--max-amplicons-msa INTEGER RANGE`: Max number of amplicons for each MSA  [x>=1]
* `--max-amplicons-region-group INTEGER RANGE`: Max number of amplicons for each region  [x>=1]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--high-gc / --no-high-gc`: Use high GC primers  [default: no-high-gc]
* `--offline-plots / --no-offline-plots`: Includes 3Mb of dependencies into the plots, so they can be viewed offline  [default: offline-plots]
* `--use-matchdb / --no-use-matchdb`: Create and use a mispriming database  [default: use-matchdb]
* `--help`: Show this message and exit.

## `primalscheme3 repair-mode`

Repairs a primer scheme via adding more primers to account for new mutations

**Usage**:

```console
$ primalscheme3 repair-mode [OPTIONS]
```

**Options**:

* `--bedfile PATH`: Path to the bedfile  [required]
* `--msa PATH`: An MSA, with the reference.fasta, aligned to any new genomes with mutations  [required]
* `--config PATH`: Path to the config.json  [required]
* `--output PATH`: The output directory  [required]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--help`: Show this message and exit.

## `primalscheme3 scheme-create`

Creates a tiling overlap scheme for each MSA file

**Usage**:

```console
$ primalscheme3 scheme-create [OPTIONS]
```

**Options**:

* `--msa PATH`: The MSA to design against. To use multiple MSAs, use multiple --msa flags. (--msa 1.fasta --msa 2.fasta)  [required]
* `--output PATH`: The output directory  [required]
* `--amplicon-size INTEGER`: The size of an amplicon. Min / max size are ± 10 percent [100<=x<=2000]  [default: 400]
* `--bedfile PATH`: An existing bedfile to add primers to
* `--min-overlap INTEGER RANGE`: min amount of overlap between primers  [default: 10; x>=0]
* `--n-pools INTEGER RANGE`: Number of pools to use  [default: 2; x>=1]
* `--dimer-score FLOAT`: Threshold for dimer interaction  [default: -26.0]
* `--min-base-freq FLOAT RANGE`: Min freq to be included,[0<=x<=1]  [default: 0.0; 0.0<=x<=1.0]
* `--mapping [first|consensus]`: How should the primers in the bedfile be mapped  [default: first]
* `--circular / --no-circular`: Should a circular amplicon be added  [default: no-circular]
* `--backtrack / --no-backtrack`: Should the algorithm backtrack  [default: no-backtrack]
* `--ignore-n / --no-ignore-n`: Should N in the input genomes be ignored  [default: no-ignore-n]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--input-bedfile PATH`: Path to a primer.bedfile containing the pre-calculated primers
* `--high-gc / --no-high-gc`: Use high GC primers  [default: no-high-gc]
* `--offline-plots / --no-offline-plots`: Includes 3Mb of dependencies into the plots, so they can be viewed offline  [default: offline-plots]
* `--use-matchdb / --no-use-matchdb`: Create and use a mispriming database  [default: use-matchdb]
* `--help`: Show this message and exit.

## `primalscheme3 scheme-replace`

Replaces a primerpair in a bedfile

**Usage**:

```console
$ primalscheme3 scheme-replace [OPTIONS] PRIMERNAME PRIMERBED MSA
```

**Arguments**:

* `PRIMERNAME`: The name of the primer to replace  [required]
* `PRIMERBED`: The bedfile containing the primer to replace  [required]
* `MSA`: The msa used to create the original primer scheme  [required]

**Options**:

* `--amplicon-size INTEGER`: The size of an amplicon. Use single value for ± 10 percent [100<=x<=2000]  [required]
* `--config PATH`: The config.json used to create the original primer scheme  [required]
* `--help`: Show this message and exit.

## `primalscheme3 visualise-bedfile`

Visualise the bedfile

**Usage**:

```console
$ primalscheme3 visualise-bedfile [OPTIONS] BEDFILE REF_PATH
```

**Arguments**:

* `BEDFILE`: The bedfile containing the primers  [required]
* `REF_PATH`: The bedfile containing the primers  [required]

**Options**:

* `--ref-id TEXT`: The reference genome ID  [required]
* `--output FILE`: Output location of the plot  [default: bedfile.html]
* `--help`: Show this message and exit.

## `primalscheme3 visualise-primer-mismatches`

Visualise mismatches between primers and the input genomes

**Usage**:

```console
$ primalscheme3 visualise-primer-mismatches [OPTIONS] MSA BEDFILE
```

**Arguments**:

* `MSA`: The MSA used to design the scheme  [required]
* `BEDFILE`: The bedfile containing the primers  [required]

**Options**:

* `--output FILE`: Output location of the plot  [default: primer.html]
* `--include-seqs / --no-include-seqs`: Reduces plot filesize, by excluding primer sequences  [default: include-seqs]
* `--offline-plots / --no-offline-plots`: Includes 3Mb of dependencies into the plots, so they can be viewed offline  [default: offline-plots]
* `--help`: Show this message and exit.
