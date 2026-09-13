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
