# Coverage panel selector: implementation scope

The user approved implementing the September 13 algorithm review. This first version adds an opt-in `coverage` selector for fresh linear whole-MSA equal panels with first-reference mapping and explicit reference-span sizing. `legacy` remains the default, unchanged. Region, imported-pair, circular, consensus and entropy optimization are rejected before scientific execution; their historical legacy modes remain available. Targeted discovery expansion and CP-SAT are later extensions, not prerequisites for this first selector.

Authoritative audit: `/Users/dho/Documents/lungfish-genome-explorer/docs/reports/2026-09-13-primalscheme-panel-algorithm-audit/report.md`.

## Runtime and integration

- Engine repository: `/Users/dho/Documents/lungfish-genome-explorer/.worktrees/primalscheme-panel-fork`, branch `codex/coverage-panel-optimizer`, based on `00eaa252446f01cabfeae71e10306a68cdb941d6`.
- LGE worktree: `/Users/dho/Documents/lungfish-genome-explorer/.worktrees/primalscheme-panel-optimizer`, branch `codex/primalscheme-panel-optimizer`.
- Develop in the fork's `.venv`, inherited dependencies from the installed conda runtime; never edit the installed managed package. New package identity `3.3.0+lge.3` is local and unpublished. LGE retains managed `3.3.0+lge.2` pins and accepts the exact new version via an explicit executable override/capability probe.
- All scientific outputs preserve exact argv, versions/build identity, explicit and resolved options, runtime identity, input/output paths, SHA-256 and byte sizes, wall time, status, useful stderr, and final durable/relative paths. Do not publish releases, update managed pins, or overwrite prior schemes.

## Public options and defaults

Native panel options: `--selection-algorithm legacy|coverage` (legacy); `--coverage-metric full-span|primer-trimmed` (full-span); `--coverage-target` (0.90); `--optimizer-seed` (0); `--optimizer-starts` (4); `--optimizer-repair-rounds` (2); `--optimizer-time-limit` (120 seconds); `--mispriming-product-size` (2000 for coverage, historical 0 for legacy). Existing size/GC/pool/core/count options remain authoritative. Coverage requires MatchDB enabled; an explicit disable flag is rejected rather than silently weakening the declared `panel-v1` constraint profile. Nondefault optimizer-only settings with legacy are rejected. LGE exposes the same option spellings except its existing pool/min-frequency spellings.

The newly exposed `--mispriming-product-size` CLI flag belongs to coverage mode. Coverage resolves an omitted value to 2000 and requires a positive integer. Legacy accepts its historical default, including an explicit zero, and rejects nonzero use of this new flag; existing programmatic legacy callers that set `Config.mismatch_product_size` retain their prior behavior. The reference-span minimum and maximum constrain the first-reference BED envelope. Row-specific products reconstructed from supplied MSA rows can differ with indels and remain diagnostic; they do not redefine the requested BED span.

The time limit bounds the selector's search, not discovery. A positive explicit wall cutoff can affect how many starts finish; scientific repeatability is guaranteed for fixed inputs/seed/work budgets when the wall cutoff is not reached, not across arbitrary time cutoffs or different hardware. Report completed starts, work counts and stop reason. A feasible empty solution must remain available if no valid candidate fits. Do not claim laboratory performance or global optimality.

## Representation

Immutable `Target`, `Candidate`, `Catalog`, `Assignment`, `SearchOptions`, `SearchResult` records separate science from mutable legacy `PrimerPair` objects. Target identity hashes the ordered alignment rows, excluding generated FASTA names. Preserve duplicate targets as deterministic occurrence IDs; never silently merge them. Candidate identity includes target identity, reference start/end and interior bounds, canonical actual 5'-3' forward/reverse oligo tuples. Export retains a mapping to native pair objects but the optimizer does not mutate them. All geometry is zero-based half-open. Reference lengths come from ungapped first rows, not alignment width.

Candidate records retain row-level intended binding support, jointly supported row IDs, unknown rows and per-row product spans. Compute this from the aligned source rows and reference-to-alignment anchors; cache per-anchor/cloud support because candidates share primers. Reverse oligos are actual reverse-complement primer sequences, not reference substrings. Terminal missing data is unknown, not imputed; at least one jointly supported observed row is required. For ambiguous IUPAC bases use explicit compatible-base matching and record uncertain support rather than treating unknown N as confirmed. The public coverage objective remains coordinate coverage; row support is diagnostic and a feasibility guard, not a claim of universal allele amplification.

Persist a deterministic compressed candidate catalogue with a schema, canonical target/oligo tables, candidate records, semantic digest, SHA-256, byte size, and source mapping. No top-N pruning in the first implementation. The caller records source discovery settings in catalogue metadata.

## Correctness and independent validation

New `panel-v1` compatibility rules are separate from legacy MatchDB and Multiplex state:

1. Candidate geometry and all oligos meet the requested bounds and existing thermochemical constraints. Use the existing thermochemical routines and native interaction kernel; do not lower thresholds to gain coverage.
2. Same-MSA full amplicon envelopes do not overlap within a pool. Pools must be within bounds, candidates unique, global/per-MSA caps respected.
3. Any oligo-cloud dimer conflict forbids same-pool coexistence. Cache pair conflicts symmetrically; immutable candidate clouds remain indivisible.
4. Specificity uses all supplied MSA rows with row identity, retaining cross-MSA hits and both orientations. Match terminal k-mers using the configured exact/single-mismatch rule. The unpublished `panel-v1` profile records specificity revision `intended-sites-v1`; its exact pure rules are in [coverage-specificity.md](coverage-specificity.md). Intrinsic checks exempt a candidate's fully observed, jointly supported, anchored intended products. Pair checks also exempt physical products proven intended for either of the two candidates and longer cross-combinations of the forward site of the left candidate and reverse site of the right candidate, only when their confirmed row-level products are disjoint, ordered and on the same target occurrence and row. Record these as allowed secondary products, never as absent. All other facing products with positive full ungapped-row span within the inclusive product bound are conflicts. Shared oligo enumeration origin cannot turn an already-proven intended product into an offsite product; exemptions never come from a third candidate or the global catalogue. Intrinsic validity remains mandatory, so adding a candidate cannot rescue prior incompatibility. Do not collapse hits from different alleles into a fabricated product. Describe the search space as supplied-MSA specificity, not whole-genome specificity.
5. Independent final validation rebuilds coverage and compatibility from the proposed assignments; it must not trust incremental search caches. It reports reference full-span and interior coverage separately, target deficits, joint/unknown support, allowed secondary products and violations. No failing solution is exported as successful.

`CompatibilityOracle` can memoize immutable chemistry/hit results for search. Its caches include all scientific settings; assignment-dependent rejection caches cannot survive a repool/removal. Reversible search state uses counted contributions or rebuilds. Do not reuse `Multiplex.remove_primerpair`.

## Search

Reuse the full frozen candidate set. Start with deterministic legacy-like and scarcity/deficit-aware constructions under the *new* compatibility rules. The search API can also consider a supplied fixed baseline assignment vector, including a historical legacy solution, only if independently valid under those rules. Native default execution does not rerun the unseeded, identity-hash-ordered legacy selector; label its deterministic baseline construction accurately. Do not promise no regression relative to an invalid or unsupplied legacy solution.

Use seeded starts with canonical semantic tie-breaks. Across target frontiers, prioritize uncovered target deficit and candidates with fewer alternatives; choose among feasible pools by balanced unique-oligo burden and remaining compatibility. Evaluate the objective lexicographically: minimize worst normalized target shortfall; minimize summed shortfall; maximize mean normalized unique coverage; minimize unique-oligo burden, amplicon count, pool imbalance and distance from requested target length. Canonical solution signature resolves final ties.

Keep the best feasible incumbent throughout. Each repair round examines undercovered targets, removes a bounded blocking neighborhood (one/two selected candidates), tries candidate swaps or repooling, and refills. Rebuild state when necessary. Accept only improvement under the declared objective. Inactive constraints cannot be ignored during repair. Use lazy symmetric conflicts, cached support/hits, interval unions/bitsets and cached pool oligo tuples rather than a dense C-squared graph. Return detailed start/move/evaluation counts and stop reason. Restrict work deterministically as well as by wall time; record any frontier/repair bounds as search limits, not candidate pruning or optimality guarantees.

## Native and LGE artifacts

Native `config.json` gains versioned `panel_optimizer` metadata only for coverage. Add `candidate-catalog.json.gz`, `panel-validation.json`, `panel-optimizer.json`, and `panel-provenance.json` under the native output. Include algorithm/objective/profile versions, all resolved options, candidate catalogue semantic/file digests, baseline/final objective vectors and per-locus metrics, seed, completed starts/repair/evaluation counts, budgets, stop reason, and independent validation success. Paths inside native records are relative to native output; executed argv is separately retained as historical fact.

Native `--capabilities-json` returns schema/version, exact package identity, supported algorithms/profile versions, and source/build/runtime identity. LGE validates this before running coverage and stores the probe as runtime evidence. Before bundle publication it verifies requested settings against config and the validation artifacts, hashes actual catalogue bytes, checks schema and successful status, and records the actual engine version. Missing or mismatched metadata blocks publication. Historical options decode with legacy defaults.

## Tests and delivery

Tests are synthetic interval/sequence or human/MHC only. Cover stable IDs under input permutation and object recreation, duplicated target handling, gapped first reference, partial/disjoint row support, asymmetric/cross-MSA hits, zero product limits, within-candidate products, reverse complement orientation, overlap across different pools, caps, deterministic seeded search, greedy counterexample repair, and timeout returning a valid incumbent. Independently enumerate tiny feasible graphs for an objective oracle. Compare incremental results to fresh validation after mutations.

Run existing 33 LGE-fork regressions; focused Swift CLI/options/publication checks; native HLA wide and narrow controls; and a deterministic interval-only 100-target scaling fixture. Label the last as synthetic algorithm scaling, not a representative 100-MSA wet-lab or biological benchmark. Preserve all benchmark provenance. No universal speedup or >90%-coverage promise is a release criterion; correctness, honest metadata and improvement of known feasible counterexamples are.

## Local use

The coverage selector is available only from the unpublished local `3.3.0+lge.3`
fork. Install it into this repository's existing virtual environment without resolving
or changing the inherited environment:

```sh
cd /Users/dho/Documents/lungfish-genome-explorer/.worktrees/primalscheme-panel-fork
.venv/bin/python -m pip install --no-deps --no-build-isolation -e .
.venv/bin/primalscheme3 --capabilities-json
```

LGE keeps the managed `3.3.0+lge.2` runtime for legacy designs. A coverage run must
use the local executable override and repeat `--msa` for every input:

```sh
cd /Users/dho/Documents/lungfish-genome-explorer/.worktrees/primalscheme-panel-optimizer
.build/debug/lungfish-cli primers design primalscheme3 \
  --msa A.lungfishmsa --msa B.lungfishmsa \
  --output HLA-coverage.lungfishprimeranalysis \
  --grouping combined \
  --primalscheme3-path ../primalscheme-panel-fork/.venv/bin/primalscheme3 \
  --amplicon-size 200 --amplicon-size-min 150 --amplicon-size-max 280 \
  --pool-count 2 --minimum-base-frequency 0.0 --core-count 4 \
  --terminal-gap-policy observed-only --dimer-score=-26.0 \
  --selection-algorithm coverage --coverage-metric full-span \
  --coverage-target 0.9 --optimizer-seed 0 --optimizer-starts 4 \
  --optimizer-repair-rounds 2 --optimizer-time-limit 120 \
  --mispriming-product-size 2000
```

The optimizer defaults shown above are full-span, target 0.9, seed 0, four starts,
two repair rounds, 120 seconds and a 2,000-base inclusive specificity-product
bound. The selector time applies only to search; discovery, final native validation,
LGE publication and bundle loading occur outside it. Coverage mode additionally
requires equal combined whole-MSA inputs, first-reference mapping, linear sequence,
MatchDB enabled, one or more pools and explicit resolved reference span bounds. Both
`legacy` and `observed-only` terminal-gap policies are supported. The native CLI
defaults to two pools and `legacy`; the LGE frontend defaults to two pools and
`observed-only`. The example and all three HLA benchmarks use two pools and
`observed-only`. Count caps default to unlimited, and omitted high-GC resolves false.
Legacy remains the default selector and retains its old argv and behavior.

## September 2026 benchmark evidence

The final software benchmark used nine human HLA MSA bundles with 20 rows each,
nominal size 200, two pools, default chemistry, supplied-MSA specificity, full-span
target 0.9, seed 0, four starts and two repair rounds per start. Native source was
clean commit `7b6148b6f80d70052fa2ed8aa3a44d4def3e8954`; LGE source was clean commit
`bc535d90f7c95eaae2c624ff87927a85de692bcd`. The wide 120-second run reached its
time limit after two starts. The wide 600-second run completed all four starts in
251.364 search seconds and returned exactly the same assignment vector and objective.
The narrow run completed all four starts in 91.653 search seconds.

| Locus | Wide 150–280 full / interior | Narrow 180–220 full / interior |
| --- | ---: | ---: |
| A | 42.7140% / 34.7905% | 38.6157% / 31.0565% |
| B | 62.6263% / 50.8724% | 55.3719% / 47.1993% |
| C | 46.7757% / 40.6903% | 56.4033% / 44.6866% |
| DPA1 | 72.9246% / 60.2810% | 58.7484% / 48.1481% |
| DPB1 | 77.9923% / 67.5676% | 58.9447% / 50.7079% |
| DQA1 | 94.9219% / 88.4115% | 90.1042% / 79.6875% |
| DQB1 | 91.4758% / 80.6616% | 56.1069% / 45.5471% |
| DRB1 | 68.9139% / 57.9276% | 64.9189% / 53.6829% |
| E | 44.3825% / 36.3045% | 68.8951% / 60.3528% |

The wide catalogue contained 50,969 candidates and selected 24 assignments. Its
120-second run evaluated 25,730,304 frontier entries and stopped at the wall limit;
the 600-second run evaluated 68,419,584 frontier entries, hit the declared per-unit
construction and repair limits, and completed with no objective improvement. The
narrow catalogue contained 16,376 candidates and selected 29 assignments; it
evaluated 26,345,472 frontier entries and completed after hitting its configured
construction and repair work limits. Narrow improved C and E relative to wide while
its worst locus, A, was lower. The lexicographic objective minimizes the worst
normalized shortfall before mean coverage and does not promise that every individual
target improves.

The corresponding LGE design wall times were 158.787 seconds for wide/120,
290.191 seconds for wide/600 and 124.484 seconds for narrow/120. Darwin
`/usr/bin/time -l` reported maximum resident-set high-water marks of 1,598,357,504,
2,295,857,152 and 1,165,656,064 bytes. These are cumulative whole-process maxima,
including discovery, validation and publication, rather than selector-only or
simultaneous-process memory. `PrimerAnalysisBundle.load` separately checked stored
bundle integrity; it did not rerun native biological validation. Independent Python
verification decompressed each final catalogue, checked every bundle/native
descriptor and stored input, reconstructed every assignment and BED row, and compared
literal full/interior bitset unions with validation. A separate root audit repeated
those checks for all three successful bundles.

The older audit reported wide C/E full-span coverage of 84.92%/77.34% with 44
selected amplicons and narrow C/E of 85.47%/73.63% with 45. Those historical runs
used legacy discovery/specificity (`D=0`) and a different selection system, while the
new results use `panel-v1`/`intended-sites-v1`, all supplied HLA rows and `D=2000`.
They are end-to-end results under materially different profiles, not a selector-only
speed or quality comparison and not evidence that the new selector is ready as a
default.

Exhaustive unary screening of the fixed wide catalogue under this exact profile found
382 of 1,354 C candidates intrinsically valid. Their literal union covers at most
913/1,101 bases (82.9246%) before pair, pool and count constraints, so the requested
90% C objective cannot be reached merely by more selector time or more pools. This is
a conditional ceiling for the current catalogue and `panel-v1`/
`intended-sites-v1`, not a biological impossibility bound. E's corresponding unary
union is 1,036/1,077 (96.1931%), but a unary bound is not an achievable scheme.
Candidate expansion, changed profiles or a later exact/CP-SAT formulation require
their own scientific and performance benchmarks.

The four configured starts are variants of one bounded heuristic. A comparison of
distinct selector approaches on the same frozen catalogue, profile and objective has
been requested as a separate experiment and remains pending; the current results do
not establish ensemble behavior.

## Allele-coverage benchmark fixture receipts

`scripts/benchmark_allele_coverage.py` freezes local MHC inputs before an
allele-aware benchmark. It reads the saved primer-analysis manifest, selects the
requested inputs by the manifest's human-readable `label`, verifies every selected
artifact's declared SHA-256 and byte size, and only then copies those bytes. The
snapshot contains `source-row-label-map.json` and `hash-manifest.json`. Repeated
labels remain separate records keyed by their original `sourceOccurrenceID`; UUID
ordering is never used to infer fixture identity. The runner refuses an existing
output directory and does not modify the source bundle.

The required command interface is:

```sh
.venv/bin/python scripts/benchmark_allele_coverage.py \
  --fixture-root /path/to/local-fixtures \
  --native-executable /path/to/primalscheme3 \
  --matrix /path/to/matrix.json \
  --output /new/immutable/evidence-directory
```

The matrix names `fixtureAnalysis`, the complete `inputLabels` selection, and
explicit comparison entries. The MHC matrix must distinguish these controls:

```json
{
  "schemaVersion": "allele-coverage-benchmark-matrix/v1",
  "fixtureAnalysis": "mhc-primal-scheme.lungfish/Analyses/MHC merged v1.lungfishprimeranalysis",
  "inputLabels": [
    "KIR2DL04.lungfishmsa", "KIR3DL10.lungfishmsa", "KIR3DS.lungfishmsa",
    "Mamu-A1.lungfishmsa", "Mamu-A2.lungfishmsa", "Mamu-A4.lungfishmsa",
    "Mamu-B.lungfishmsa", "Mamu-DPA.lungfishmsa", "Mamu-DQB.lungfishmsa",
    "Mamu-DRB.lungfishmsa", "Mamu-E.lungfishmsa"
  ],
  "comparisons": [
    {
      "id": "historical-independent-lge.2", "kind": "saved-artifact",
      "root": "mhc-primal-scheme.lungfish/Analyses/MHC v1.lungfishprimeranalysis",
      "manifest": "manifest.json"
    },
    {
      "id": "historical-combined-lge.2", "kind": "saved-artifact",
      "root": "mhc-primal-scheme.lungfish/Analyses/MHC merged v1.lungfishprimeranalysis",
      "manifest": "manifest.json"
    },
    {"id": "existing-lge.3", "kind": "pending-execution"},
    {
      "id": "upstream-original-3.3.0",
      "kind": "pending-execution"
    }
  ]
}
```

The upstream-original control is the clean upstream PrimalScheme repository
(`artic-network/primalscheme3`, version 3.3.0, commit prefix `60455e9`) in an
isolated environment. It is distinct from saved LGE legacy outputs. Historical
independent/combined lge.2 artifacts, the existing lge.3 result, and upstream
original use different scientific profiles and remain diagnostic comparisons until
they are re-evaluated under one declared metric and specificity contract.
`saved-artifact` entries require a manifest and freeze every artifact only after
checking its declared hash and size. Work with no saved result is
`pending-execution`; a name alone is never treated as a frozen control.

An executable row uses a measured identity probe, an expected identity, and a
complete resolved-option declaration. This is the native allele-aware baseline
interface (the executable path can instead be supplied by `--native-executable`):

```json
{
  "id": "native-allele-coverage",
  "kind": "native-execution",
  "workingDirectory": "/clean/tool/source/root",
  "identityProbeArgv": ["--capabilities-json"],
  "expectedToolIdentity": {
    "name": "primalscheme3",
    "version": "3.3.0+lge.4",
    "gitCommit": "FULL_OR_REVIEWED_COMMIT_PREFIX"
  },
  "argv": [
    "panel-create", "--mode", "equal", "{msa_args}",
    "--output", "{output}",
    "--amplicon-size", "200", "--amplicon-size-min", "150",
    "--amplicon-size-max", "250", "--n-pools", "2",
    "--min-base-freq", "0.0", "--mapping", "first", "--ncores", "4",
    "--terminal-gap-policy", "observed-only", "--dimer-score", "-26",
    "--selection-algorithm", "allele-coverage",
    "--candidate-profiles", "union",
    "--coverage-metric", "observed-allele-primer-trimmed",
    "--coverage-target", "0.95", "--allele-weighting", "distinct-observed",
    "--specificity-terminal-k", "17", "--mispriming-product-size", "2000",
    "--optimizer-seed", "0", "--optimizer-starts", "4",
    "--optimizer-repair-rounds", "2", "--optimizer-time-limit", "120",
    "--subset-beam-width", "16", "--subset-expansion-limit", "256",
    "--exchange-width", "2", "--salvage", "off", "--primary-tier", "strict"
  ],
  "optionsFullyResolved": true,
  "resolvedOptions": {
    "mode": "equal", "ampliconSize": 200, "ampliconSizeMin": 150,
    "ampliconSizeMax": 250, "poolCount": 2, "minimumBaseFrequency": 0.0,
    "mapping": "first", "coreCount": 4, "terminalGapPolicy": "observed-only",
    "dimerScore": -26, "selectionAlgorithm": "allele-coverage",
    "candidateProfiles": "union",
    "coverageMetric": "observed-allele-primer-trimmed", "coverageTarget": 0.95,
    "alleleWeighting": "distinct-observed", "specificityTerminalK": 17,
    "misprimingProductSize": 2000, "optimizerSeed": 0, "optimizerStarts": 4,
    "optimizerRepairRounds": 2, "optimizerTimeLimitSeconds": 120,
    "subsetBeamWidth": 16, "subsetExpansionLimit": 256, "exchangeWidth": 2,
    "salvage": "off", "primaryTier": "strict"
  },
  "optionContract": {
    "schemaVersion": "allele-coverage-resolved-options/v1",
    "toolVersion": "3.3.0+lge.4",
    "requiredResolvedOptionKeys": [
      "mode", "ampliconSize", "ampliconSizeMin", "ampliconSizeMax",
      "poolCount", "minimumBaseFrequency", "mapping", "coreCount",
      "terminalGapPolicy", "dimerScore", "selectionAlgorithm",
      "candidateProfiles", "coverageMetric", "coverageTarget",
      "alleleWeighting", "specificityTerminalK", "misprimingProductSize",
      "optimizerSeed", "optimizerStarts", "optimizerRepairRounds",
      "optimizerTimeLimitSeconds", "subsetBeamWidth", "subsetExpansionLimit",
      "exchangeWidth", "salvage", "primaryTier"
    ],
    "argvOptionMap": {
      "--mode": "mode", "--amplicon-size": "ampliconSize",
      "--amplicon-size-min": "ampliconSizeMin",
      "--amplicon-size-max": "ampliconSizeMax", "--n-pools": "poolCount",
      "--min-base-freq": "minimumBaseFrequency", "--mapping": "mapping",
      "--ncores": "coreCount", "--terminal-gap-policy": "terminalGapPolicy",
      "--dimer-score": "dimerScore", "--selection-algorithm": "selectionAlgorithm",
      "--candidate-profiles": "candidateProfiles", "--coverage-metric": "coverageMetric",
      "--coverage-target": "coverageTarget", "--allele-weighting": "alleleWeighting",
      "--specificity-terminal-k": "specificityTerminalK",
      "--mispriming-product-size": "misprimingProductSize",
      "--optimizer-seed": "optimizerSeed", "--optimizer-starts": "optimizerStarts",
      "--optimizer-repair-rounds": "optimizerRepairRounds",
      "--optimizer-time-limit": "optimizerTimeLimitSeconds",
      "--subset-beam-width": "subsetBeamWidth",
      "--subset-expansion-limit": "subsetExpansionLimit",
      "--exchange-width": "exchangeWidth", "--salvage": "salvage",
      "--primary-tier": "primaryTier"
    }
  },
  "actualResolvedOptions": {
    "path": "config.json",
    "jsonPath": ["panel_optimizer", "options"]
  }
}
```

`{msa_args}` expands to one `--msa PATH` pair for each verified snapshot MSA in
manifest order. `{output}` expands to a new path inside the run receipt directory.
The identity probe must return `tool`, `toolVersion`, `source`, and `runtime`; the
runner compares the measured name, version, and source commit with
`expectedToolIdentity` before scientific execution. Measured source identity must
include a nonempty Git commit and source digest. Measured runtime identity must
include the interpreter, complete kernel identity, and at least one versioned
runtime dependency. The version-bound `optionContract` requires every listed key
and compares each mapped argv value with `resolvedOptions` before launch. When
`actualResolvedOptions` is present, the runner hashes the native config and verifies
its resolved values after execution; any mismatch makes the run fail publication.
The same explicit `workingDirectory` is used for the identity probe and scientific
process so Python import resolution cannot select a neighboring checkout. Optional
`identityProbeArtifacts` are copied into the run and hashed in its receipt. For
value-taking flags, `argvOptionMap` maps the flag directly to a resolved key (or
uses `{"key": "name", "valueFromNextArg": true}`). Boolean switches use
`{"key": "name", "value": true}` or `false`, making their semantics explicit.

Each executed matrix row writes this minimum receipt shape, with additional schema,
workflow, shell-command and stdout fields allowed:

```python
receipt = {
    "argv": argv,
    "resolvedOptions": resolved,
    "inputArtifacts": inputs,
    "outputArtifacts": outputs,
    "sourceIdentity": source,
    "runtimeIdentity": runtime,
    "exitStatus": completed.returncode,
    "wallTimeSeconds": elapsed,
    "stderrPath": stderr_path,
}
```

Scientific input paths in an execution receipt resolve beneath the final snapshot
and are rehashed immediately before invocation. Original-bundle descriptors remain
separate source provenance. Paths in output records resolve beneath the final evidence directory. Receipts
retain exact argument arrays and reproducible shell rendering, workflow and version,
explicitly complete resolved options/defaults supplied by the matrix, measured tool
version/source/runtime identity, harness identity, input/output hashes and sizes, status, wall time,
and captured stderr. A failed subprocess still writes its receipt and stderr. A
launch failure writes the expanded per-run receipt. A fixture checksum failure
writes attempted manifest/input identity plus expected and observed hashes in runner
failure provenance and prevents scientific execution.

The deterministic abstract benchmark exercised the real search core on 100 synthetic
300-base interval targets, 300 candidates and 200 declared overlap edges. One pool,
one start and no repair selected one 260-base greedy interval per target
(26,000/30,000 bases); one repair round selected two disjoint 150-base intervals per
target (30,000/30,000) with graph and pool feasibility independently checked. The
runs took 0.420 and 13.053 seconds. This is abstract-oracle algorithm evidence, not
100 biological MSAs, primer chemistry or wet-lab scalability. Its `tracemalloc` peak
includes retained traced baseline objects because resetting the peak does not subtract
them, and process RSS is a cumulative maximum that includes validation before each
sample.

Full local commands, hashes, runtime/native-extension identity, input inventories,
resource measurements, stop reasons, work counts and final-byte checks are retained
under the ignored `outputs/task6-*` evidence roots. Published bundles preserve their
own complete provenance; the benchmark harness refuses to overwrite an existing
evidence directory.
