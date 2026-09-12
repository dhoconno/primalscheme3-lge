# PrimalScheme3-LGE custom fork

## Identity and scope

Public repository: https://github.com/dhoconno/primalscheme3-lge

Custom source version: `3.3.0+lge.2`; release tag when published: `v3.3.0-lge.2`.

Upstream baseline: PrimalScheme3 v3.3.0, commit `dd13ec5cb1cf375f052640355c73101c0c4bf839`. Upstream retains its authorship and copyrights. The fork remains GPL-3.0 licensed. It is not an upstream release or endorsement.

The distribution/import/executable name remains `primalscheme3` for compatibility. An installed wheel from this version emits `PrimalScheme3-LGE version: 3.3.0+lge.2`. The release wheel contains `lge-build.json` identifying the exact fork and upstream commits and build tooling.

## Selectable discovery behavior

Both `scheme-create` and `panel-create` accept `--terminal-gap-policy legacy|observed-only`. The standalone default is legacy. LGE explicitly selects observed-only by default, with a user option for legacy instead.

Legacy preserves upstream Rust discovery. Observed-only applies the terminal-gap changes through the existing Python discovery path, with compatibility fixes for the pinned native dependency. Internal gaps remain observations; terminal missing data are neither imputed nor counted as matches. Partially sequenced records retain their useful internal observations. This does not establish compatibility for unobserved alleles.

Native configuration identifies `rust-legacy` or `python-observed-only`. Experimental annealing/downsampling options unsupported by the Python path are rejected explicitly. Downstream pairing and scheme/panel construction remain shared. Different discovery backends may differ beyond terminal-gap frequency handling.

## Multicore execution

Observed-only supports `--ncores`. It uses spawned processes over bounded batches of candidate positions, sends primitive records across process boundaries, and reconstructs alternatives in stable order. One-worker execution uses the same canonical record conversion. Requested and effective counts are retained in configuration, including per-alignment counts for combined panels.

The worker count is capped by available processors and work. Each process has an alignment copy; memory requirements therefore increase with worker count. Exceptions, interruption and termination close the pool and reap its workers. LGE defaults to up to four workers and allows the user to change that count.

## Explicit amplicon size bounds

`scheme-create` and `panel-create` accept optional positive-integer `--amplicon-size-min` and `--amplicon-size-max` flags. Supplying either selects the persisted `amplicon_size_metric: reference-span` contract. An omitted individual bound resolves to the existing nominal target minus/plus 10 percent, truncated to an integer. The resolved values must satisfy `0 < minimum <= target <= maximum`; the nominal `--amplicon-size` target remains restricted to 100–2000. It sets defaults and does not add a preference for products closest to the target.

For example, `--amplicon-size 200 --amplicon-size-min 150 --amplicon-size-max 250` accepts reported full amplicon spans from 150 through 250 bases, inclusive. Length is exactly `amplicon.bed end - start`: the zero-based, half-open envelope from the earliest forward-primer start to the latest reverse-primer end on the selected reference. Both primers are included. With alternative oligo lengths this bounds the envelope reported for the cloud, not each shorter individual product or the product length on every aligned sequence. Primer-trimmed BED intervals keep their existing meaning.

Both discovery backends share this filter. Candidate search includes the reverse-primer extension before the exact span check, so products at the lower boundary remain eligible. Existing ordering, scoring and interaction checks are unchanged.

Reference-span design requires `--mapping first` (the default), measuring bases on the first input record after removing alignment gaps. Consensus mapping uses alignment columns and is rejected with this metric. Supported workflows are fresh linear schemes and whole-MSA panels in `--mode equal` or `--mode entropy`, without region BED files or imported primer pairs. Circular schemes, scheme `--bedfile`/`--input-bedfile`, panel `--input-bedfile`, and panel region modes/files are rejected before output creation when explicit bounds are supplied. Panel's default mode is region-only, so select equal or entropy for explicit bounds.

`replace-primerpair` also rejects a saved reference-span configuration before creating outputs; it cannot silently reinterpret these bounds through its legacy replacement path.

With neither new flag, configuration records `amplicon_size_metric: legacy-pairing`, and scheme pairing keeps its previous reverse-start minus earliest-forward-start calculation. Loading an older configuration with populated min/max fields does not opt into reference-span semantics; the new metric must be saved explicitly. Version lge.2 also fixes panel creation passing its maximum as its minimum: panel now uses the configured interval in either metric. This panel correction can change candidates even without new flags.

## Installation and rebuilding

Download a published wheel and SHA256SUMS from the [custom GitHub releases](https://github.com/dhoconno/primalscheme3-lge/releases). Verify its checksum, then install that downloaded wheel into an isolated Python environment. Python 3.12.11 is the tested runtime and reproducible-build interpreter. The native dependency `primalschemers==0.1.13` is pinned because the Python adapter relies on its interface.

On Apple Silicon, LGE supplies Python and `primer3-py` through conda and installs the remaining pinned wheels into that environment. A wheel for every native dependency may not be available from PyPI for every platform; the custom pure-Python wheel does not bundle those native libraries.

Run `primalscheme3 --version` and `primalscheme3 scheme-create --help` to verify custom identity and options. Follow `scripts/build_lge_wheel.py --help` to rebuild from a clean committed checkout using the exact tooling in `requirements-build.txt`. The script builds from a git archive with a fixed source epoch, embeds source identity, and writes artifact evidence. Automatic upstream PyPI publishing and documentation deployment workflows are disabled in this fork; its CI exercises synthetic regression tests.

## Validation evidence and limits

The focused custom suite covers missing-data semantics, current native coordinate mapping, CLI/configuration dispatch, nonempty one/two/four-worker candidate equivalence, worker and consumer failures, and parent termination cleanup. Synthetic coordinate tests cover inclusive full-span limits, alternative primer lengths, old configuration compatibility, CLI scope/range validation, and scheme/panel propagation. The upstream Config suite checks compatibility of existing defaults.

Before the lge.2 size changes, nineteen local native runs on authentic human and macaque MHC inputs succeeded: twelve independent A/DPA1/DPB1/DQA1/DQB1/DRB1 cases, six combined DP/DQ/six-locus panels, and a repeated human A case with one worker. All inputs were complete supplied genomic records; biological completeness and laboratory assay performance were not established. Those runs do not validate the new explicit size contract.

For the human A fixture, one worker took 5.57 seconds and four took 5.10 seconds in a single local Apple Silicon comparison. The 22 BED records had identical scientific fields and order; generated identifiers differed. These are individual measurements, not a general speedup claim. The six-locus panels completed in approximately 39 seconds for human and 31 seconds for macaque with four workers.

LGE preserves full scientific input/output provenance in its native result bundles. Local validation datasets, user paths, credentials and agent work records are not included in this public source repository.
