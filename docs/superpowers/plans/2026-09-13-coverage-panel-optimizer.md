# Coverage Panel Optimizer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Implement an opt-in reproducible, validated joint candidate-selection/pool-repair algorithm and expose it through the LGE CLI.

**Architecture:** Frozen candidate records, a new versioned compatibility oracle and independent validator, a bounded multistart repair selector, then native/LGE integration. Legacy creation remains a separate default path. No dense all-candidate graph or new solver dependency.

**Tech Stack:** Python 3.12, NumPy, primalschemers/primer3, Typer, unittest/pytest; Swift 6 ArgumentParser and existing LGE workflow/provenance modules.

**Spec:** `/Users/dho/Documents/lungfish-genome-explorer/.worktrees/primalscheme-panel-fork/docs/designs/coverage-panel.md`.

## Global Constraints

- Work only in the two isolated repositories named in the spec; no installed package changes, published releases, managed pin updates or prior-output overwrites.
- Engine new identity `3.3.0+lge.3`; legacy managed identity `3.3.0+lge.2` remains supported by LGE.
- Native algorithm choices `legacy|coverage`, metric `full-span|primer-trimmed`; defaults and capability/provenance behavior are exactly as specified.
- Current scope: fresh whole-MSA equal panels, first-row mapping, explicit reference-span metric. Reject unsupported combinations before scientific output.
- All scientific artifacts must retain complete reproducibility provenance, correct sizes/hashes and durable/relative stored paths.
- Tests use synthetic or human/MHC inputs; no production-code test-only helpers.
- Implementers do not spawn agents. Record failing-test evidence, implementation, passing checks, commit and concerns in their report file; return only status/commit/test summary.

## Task 1: Immutable catalogue and coordinate/support model

**Files:** Create `primalscheme3/panel/coverage_types.py`, `primalscheme3/panel/coverage_catalog.py`, `tests/lge/test_coverage_catalog.py`.

**Interfaces:** Produce frozen `Target`, `Candidate`, `Catalog`, `Assignment` dataclasses; `build_catalog(msa_dict, config) -> Catalog`; `write_catalog(catalog, path) -> dict`. Catalog exposes ordered targets/candidates, ID maps, native-pair mapping, semantic digest. Candidate exposes id/target_id, full and interior reference intervals, canonical actual forward/reverse oligo tuples, jointly supported/unknown row IDs and observed row product spans. Target retains source MSA index and rows/reference mapping. Native-pair mapping is kept outside candidate equality/digests.

- [x] Write failing tests with hand-derived interval and row support expectations. Include two recreated native `FKmer/RKmer` objects yielding equal candidate IDs; input permutations yielding equal catalogue digest; a gapped reference whose denominator is ungapped length; partial rows with F-only/R-only support yielding no joint row; actual reverse oligo matching reverse-complement sequence; duplicate target inputs kept distinct.
- [x] Run `.venv/bin/python -m pytest tests/lge/test_coverage_catalog.py -q` and retain the failure.
- [x] Implement immutable canonical records and cached support reconstruction. Use `@dataclass(frozen=True)`; canonical JSON with sorted keys and compact separators feeds SHA-256. Example geometry test contract: `assert candidate.full_interval == (10, 50)` and `assert candidate.interior_interval == (18, 42)` for eight-base oligos at those anchors. Unknown-only support never enters `joint_rows`; a cloud may still report uncertainty for a row that another fully observed member confirms.
- [x] Write deterministic gzip with `mtime=0`, containing deduplicated targets/oligos and candidates; return semantic digest plus actual file SHA-256/size. Verify decoded rows against literal fixture data, and byte equality on two writes.
- [x] Run task tests and existing `tests/lge`; commit task changes.

## Task 2: Versioned compatibility oracle and independent validator

**Files:** Create `primalscheme3/panel/coverage_specificity.py`, `primalscheme3/panel/coverage_validation.py`, `tests/lge/test_coverage_validation.py`; consume Task 1 types without rewriting their contract.

**Interfaces:** `ConstraintProfile` records exact bounds, pool/count limits, thermochemical config and positive specificity product-size/kmer/fuzzy settings; `CompatibilityOracle(catalog, profile)` provides `candidate_valid(id)`, `conflict(id_a,id_b)` with reasons and counters. `validate_assignments(catalog, assignments, profile, metric, coverage_target) -> dict` rebuilds and returns schemaVersion, valid, violations, per-target full/interior metrics, support diagnostics. `Assignment(candidate_id: str, pool: int)` uses zero-based pool indices.

- [x] Tests first: intended same-row product accepted; disjoint same-target intended-site cross-combinations accepted in either order and reported as allowed secondary products; inward overlap, off-anchor, and cross-MSA products rejected regardless of insertion order; hits on different rows never form a product; unexpected within-candidate product rejected; zero/negative product bound rejected; true reverse orientation; exact boundary behavior; an overlapping interior in another pool remains covered after removing a candidate from a proposed assignment list. Explicit expected coverage: remaining interval `[30,70)` covers 40 bases, not 20. Shared-oligo enumeration must not invalidate a physical product intended for either evaluated candidate; no third candidate may rescue a conflict. See the literal counterexamples in `docs/designs/coverage-specificity.md`.
- [x] Run the focused file red. Implement complete row-indexed terminal-kmer hit sets and the exact `intended-sites-v1` exemptions in `docs/designs/coverage-specificity.md`. The maximum undesired-product size is an inclusive positive full ungapped-row span bound; plus 3-prime end must be at or before minus 3-prime end. Retain all source rows and cross-MSA hits. Ambiguous/unknown support semantics follow spec; no silent wildcard N support. Profile serialization/cache identity records the specificity revision and the allowed-secondary-product policy explicitly.
- [x] Implement pure candidate and pair constraints with symmetric cache keys, and an independent full validator with a fresh uncached evaluation path and fresh interval unions. Recheck bounds, caps, unique IDs, pools, all selected oligos and all active constraints. `assert not validate_assignments(...invalid...)['valid']` must catch a deliberately invalid pool assignment without relying on optimizer state.
- [x] Run task tests plus existing/catalogue regressions. Commit and document checker model limits, especially supplied-MSA rather than genomic background.

## Task 3: Bounded deterministic multistart selection and repair

**Files:** Create `primalscheme3/panel/coverage_search.py`, `tests/lge/test_coverage_search.py`; optional helper `coverage_objective.py` if needed for cohesion.

**Interfaces:** Frozen `SearchOptions(metric='full-span',coverage_target=.9,seed=0,starts=4,repair_rounds=2,time_limit=120)` and `SearchResult(assignments,metadata)`; `optimize_catalog(catalog,profile,options,baseline=()) -> SearchResult`. Objective is lexicographic as in spec, documented and serialized. Time source may be passed as an ordinary optional callable to the search function only if needed for deterministic timeout testing; do not create test-only production methods.

- [x] Write failing tests: the real catalogue/constraint-compatible equivalent of the audit's A/B greedy trap, three-pool repool opportunity, fixed-seed repeatability, input/candidate permutation invariance, count caps, sparse lazy conflict evaluation, no candidate duplicated, empty/unreachable target, and time budget returning validated incumbent.
- [x] Independently enumerate all assignments for a tiny literal candidate graph and calculate expected optimal coverage by set union in test code; the selected objective must match for specified tractable examples. No assertion generated by the production objective itself.
- [x] Implement canonical deterministic initial construction, seeded alternative starts, target deficit/scarcity ranking, balanced feasible pool choice and lazy memoized constraints. Use the complete catalogue with bounded frontier evaluation/search limits recorded in metadata; never imply a shortlist is the complete search.
- [x] Implement one/two-candidate bounded remove/repool/refill repair using fresh/counted state, preserving best validated objective. A valid baseline is retained; invalid legacy baseline is rejected with recorded reasons. Record starts, evaluations, repairs, objective history and stop reason (`completed` or `time-limit`). Wall time can interrupt only between bounded work units; explain timing/determinism limits.
- [x] Run all new tests and legacy fork regressions; commit.

## Task 4: Native CLI, artifacts, capabilities and local package

**Files:** Modify `primalscheme3/cli.py`, `primalscheme3/core/config.py`, `primalscheme3/panel/panel_main.py`, `pyproject.toml`; create focused `coverage_pipeline.py`/`coverage_provenance.py` as needed, `tests/lge/test_coverage_cli.py`, `tests/lge/test_coverage_publication.py`; update `docs/lge-fork.md`.

**Interfaces:** Native options/defaults exactly from spec. `--capabilities-json` emits schema1, actual toolVersion, algorithms, profile/optimizer schema versions, source/build and Python/runtime identity. Coverage `config.json.panel_optimizer` references versioned catalogue/validation/search/provenance artifacts and all resolved options. Task 5 consumes this verified JSON contract.

- [x] Test CLI legacy argv/default behavior and coverage unsupported-mode rejection before creation; test options finite/range validation and capability output. Legacy does not start new selector or alter historical sizing/constraints.
- [x] Route coverage after existing candidate generation into Tasks 1–3, then fresh publication to native BED/reference/plots with stable semantic amplicon naming. Do not mutate shared candidate records. Persist catalogue/validation/search/provenance files before return; independently invalid results fail closed.
- [x] Native provenance records executed argv, exact options/defaults, source file identity, runtime, input/output hashes/sizes, relative output paths, wall time and status/stderr. Candidate catalogue bytes and semantic digests are different named fields. Report version/algorithm/profile, metric/threshold, seed/budgets/completed work, baseline/final vector and validation validity.
- [x] Bump only local fork identity to `3.3.0+lge.3`, install editable with `.venv/bin/python -m pip install --no-deps --no-build-isolation -e .`, run CLI help/capabilities and synthetic native smoke. Do not publish or alter the managed conda installation.
- [x] Run fork regressions, focused CLI/publication tests and one human HLA command with complete provenance. Commit and report exact JSON contract for the LGE implementer.

## Task 5: LGE CLI and publication integration

**Files in the LGE worktree:** `Sources/LungfishWorkflow/PrimerDesign/PrimalScheme3DesignPipeline.swift`, `Sources/LungfishCLI/Commands/PrimerDesignCommand.swift`, `Tests/LungfishWorkflowTests/PrimalScheme3DesignPipelineTests.swift`, `Tests/LungfishWorkflowTests/PrimalScheme3PublicationTests.swift`, `Tests/LungfishCLITests/PrimerDesignCommandTests.swift`, `docs/features/primalscheme3-lge-fork.md`.

**Interfaces:** Map native selector options with backward-compatible Codable defaults. Legacy sends no new flags to managed lge.2. Coverage uses exact lge.3 capability negotiation via executable override; store probe/build evidence. Validate all Task4 metadata and actual catalogue/validation artifact content before publishing.

- [x] Add failing Swift tests for old payload decode, unchanged legacy argv, coverage-only combined validation and every invalid numeric/mode combination; CLI parse forwarding and resolved provenance fields.
- [x] Add publication tests that reject mismatched algorithm/profile/version/seed/budget/coverage/catalogue hash and failed or missing independent validation atomically. Acceptance records actual lge.3 and final stored paths; legacy lge.2 remains accepted.
- [x] Implement fields/arguments/capability handshake without changing registry/lock/managed pins or promising an unpublished managed install. Reject insufficient old runtime before scientific execution with a concrete executable-override instruction.
- [x] Run `scripts/setup-worktree.sh "$PWD"`, then focused `swift test --skip-update --filter 'PrimalScheme3DesignPipelineTests|PrimalScheme3PublicationTests|PrimerDesignCommandTests'`, build `swift build --skip-update --product lungfish-cli` and inspect help. Record environmental baseline failures separately and address task-related failures.
- [x] Commit LGE code/tests/docs and record built CLI path. GUI exposure is outside this CLI-first delivery; no partially functional UI controls.

## Task 6: End-to-end benchmark, documentation and final review

**Files:** Fork `scripts/benchmark_coverage.py`, synthetic benchmark tests, `docs/coverage-panel.md`; LGE report/documentation as needed. Scientific results stay in ignored `outputs/` with provenance.

- [x] Run deterministic 100-target interval-only scaling fixture with real search and declared synthetic conflict graph, verify final validity and greedy-trap improvement, record stage wall time/work/memory where available. Mark it synthetic, not a 100-MSA biological performance estimate.
- [x] Run human HLA combined coverage via built LGE CLI + local lge.3 override using 200 target,150–280 bounds. Compare original metrics under clearly stated different compatibility profiles; also run explicit180–220 coverage control. Never present profile differences as selector-only speed/quality improvements.
- [x] Independently verify bundle hashes/paths, selected BED spans and validation metrics; retain command/runtime/source/inputoutput manifests for every run, including failure logs if any.
- [x] Document install/override command, CLI examples, defaults, metrics/unknown allele support, legacy behavior, constraints/profile limitations, reproducibility and timeout scope, benchmark results and future candidate/CP-SAT extensions.
- [ ] Run final focused Python/Swift checks and Astra whole-branch review. Resolve all important findings through implementer/reviewer loops, then deliver both branches, binaries/artifacts, tests and limitations. Do not merge/push/publish without authorization.
