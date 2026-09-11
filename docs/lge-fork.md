# PrimalScheme3-LGE custom fork

## Identity and scope

Public repository: https://github.com/dhoconno/primalscheme3-lge

Custom version: `3.3.0+lge.1`; release tag: `v3.3.0-lge.1`.

Upstream baseline: PrimalScheme3 v3.3.0, commit `dd13ec5cb1cf375f052640355c73101c0c4bf839`. Upstream retains its authorship and copyrights. The fork remains GPL-3.0 licensed. It is not an upstream release or endorsement.

The distribution/import/executable name remains `primalscheme3` for compatibility. The version probe emits `PrimalScheme3-LGE version: 3.3.0+lge.1`. The release wheel contains `lge-build.json` identifying the exact fork and upstream commits and build tooling.

## Selectable discovery behavior

Both `scheme-create` and `panel-create` accept `--terminal-gap-policy legacy|observed-only`. The standalone default is legacy. LGE explicitly selects observed-only by default, with a user option for legacy instead.

Legacy preserves upstream Rust discovery. Observed-only applies the terminal-gap changes through the existing Python discovery path, with compatibility fixes for the pinned native dependency. Internal gaps remain observations; terminal missing data are neither imputed nor counted as matches. Partially sequenced records retain their useful internal observations. This does not establish compatibility for unobserved alleles.

Native configuration identifies `rust-legacy` or `python-observed-only`. Experimental annealing/downsampling options unsupported by the Python path are rejected explicitly. Downstream pairing and scheme/panel construction remain shared. Different discovery backends may differ beyond terminal-gap frequency handling.

## Multicore execution

Observed-only supports `--ncores`. It uses spawned processes over bounded batches of candidate positions, sends primitive records across process boundaries, and reconstructs alternatives in stable order. One-worker execution uses the same canonical record conversion. Requested and effective counts are retained in configuration, including per-alignment counts for combined panels.

The worker count is capped by available processors and work. Each process has an alignment copy; memory requirements therefore increase with worker count. Exceptions, interruption and termination close the pool and reap its workers. LGE defaults to up to four workers and allows the user to change that count.

## Installation and rebuilding

Download the wheel and SHA256SUMS from the [custom GitHub release](https://github.com/dhoconno/primalscheme3-lge/releases/tag/v3.3.0-lge.1). Verify its checksum, then install that downloaded wheel into an isolated Python environment. Python 3.12.11 is the tested runtime and reproducible-build interpreter. The native dependency `primalschemers==0.1.13` is pinned because the Python adapter relies on its interface.

On Apple Silicon, LGE supplies Python and `primer3-py` through conda and installs the remaining pinned wheels into that environment. A wheel for every native dependency may not be available from PyPI for every platform; the custom pure-Python wheel does not bundle those native libraries.

Run `primalscheme3 --version` and `primalscheme3 scheme-create --help` to verify custom identity and options. Follow `scripts/build_lge_wheel.py --help` to rebuild from a clean committed checkout using the exact tooling in `requirements-build.txt`. The script builds from a git archive with a fixed source epoch, embeds source identity, and writes artifact evidence. Automatic upstream PyPI publishing and documentation deployment workflows are disabled in this fork; its CI exercises synthetic regression tests.

## Validation evidence and limits

The focused custom suite covers missing-data semantics, current native coordinate mapping, CLI/configuration dispatch, nonempty one/two/four-worker candidate equivalence, worker and consumer failures, and parent termination cleanup. The upstream Config regression suite also passes.

Nineteen local native runs on authentic human and macaque MHC inputs succeeded: twelve independent A/DPA1/DPB1/DQA1/DQB1/DRB1 cases, six combined DP/DQ/six-locus panels, and a repeated human A case with one worker. All inputs were complete supplied genomic records; biological completeness and laboratory assay performance were not established.

For the human A fixture, one worker took 5.57 seconds and four took 5.10 seconds in a single local Apple Silicon comparison. The 22 BED records had identical scientific fields and order; generated identifiers differed. These are individual measurements, not a general speedup claim. The six-locus panels completed in approximately 39 seconds for human and 31 seconds for macaque with four workers.

LGE preserves full scientific input/output provenance in its native result bundles. Local validation datasets, user paths, credentials and agent work records are not included in this public source repository.
