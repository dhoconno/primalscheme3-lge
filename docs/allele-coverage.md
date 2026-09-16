# Allele-aware PrimalScheme CLI controls

This development CLI is version 3.3.0+lge.4. Evaluation is ongoing; this document describes the implemented contract, not demonstrated improvement on a particular panel. Use an explicitly verified development executable (`primalscheme3 --version` and `--capabilities-json`).

## Starting point

```sh
primalscheme3 panel-create \
  --msa /absolute/path/target-a.fasta --msa /absolute/path/target-b.fasta \
  --output /absolute/path/new-panel \
  --selection-algorithm allele-coverage --preset allele-balanced-v1 \
  --amplicon-size 200 --amplicon-size-min 150 --amplicon-size-max 250 \
  --n-pools 2 --ncores 4
```

Repeat `--msa` for each target. Output must be a new directory. Supply the original aligned sequences, preserving unknown bases and alignment gaps. The preset combines normal and high-GC discovery, gives each distinct observed allele equal weight within its target, and aims for 95% coverage after primer trimming. Coverage is a predicted binding model, not measured amplification success.

## Scientific controls

| Option | Preset | What changing it does |
|---|---|---|
| `--candidate-profiles` | `union` | `normal` or `high-gc` restricts discovery; union retains profile provenance for shared candidates. |
| `--variant-selection` | `subsets` | Allows compatible variants to preserve a partially supported amplicon. `full-cloud` requires every eligible variant in a family and is useful as a comparison. |
| `--coverage-target` | `0.95` | Changes the objective's deficit penalty; it is an aspiration, not a requirement that discards lower-coverage panels. |
| `--amplicon-size`, `--amplicon-size-min`, `--amplicon-size-max` | Explicit bounds required | Defines reference-span geometry. Wider bounds admit more placement alternatives and increase candidate volume. |
| `--n-pools` | `2` | More physical reaction pools can separate conflicts but increase laboratory reactions. |
| `--min-base-freq` | `0` | Increasing it can remove rare allele-specific variants. Frequencies use a fixed anchored row cohort. |
| `--discovery-length-mode` | `first-compatible` | `all` examines all allowed lengths instead of stopping at the first individually compatible length per row, anchor and profile. This can be expensive. |
| `--specificity-terminal-k` | `17` | Changes the terminal seed used for supplied-row specificity screening. Comparisons must account for primer lengths and use a shared eligible catalog. |
| `--mispriming-product-size` | `2000` | Upper product length screened on supplied rows; positive inclusive bound. This mode does not permit disabling screening with zero. |
| `--intended-product-policy` | `exact-supported` | Opt-in `concrete-designated-sites` tolerates known potential products at the selected sites’ own concrete row footprints, without adding coverage. See below. |
| `--secondary-product-policy` | `ordered-disjoint-intended-sites` | Allows the declared class of products between ordered disjoint intended sites, without coverage credit. `reject-secondary-products/v1` rejects secondary products. |
| `--max-amplicons`, `--max-amplicons-msa` | Uncapped | Limits physical design size; caps may reduce achievable coverage. |

The strict dimer cutoff remains −26. The score is the native numerical score, not free energy or a probability. Chemistry is defined by complete versioned profiles; ad hoc temperature/salt overrides are not silently accepted. Both profiles use the same modeled temperature window. That compatibility does not establish laboratory performance.

## Search effort

`--optimizer-seed` (0), `--optimizer-starts` (4), `--optimizer-repair-rounds` (2), and `--optimizer-time-limit` (120 seconds) control search. `--subset-beam-width` (16) and `--subset-expansion-limit` (256) control how many variant subsets are considered. `--exchange-width` (2; maximum 2) controls replacement neighborhood size. Larger budgets may find improvements; they do not prove an optimum. Wall-clock budgets can interrupt different work on different machines.

Advanced deterministic work limits are also exposed:

| Option | Default |
|---|---:|
| `--work-frontier-candidates` | 64 |
| `--work-construction-candidate-attempts` | 2048 |
| `--work-repair-candidate-probes-per-round` | 128 |
| `--work-repair-neighborhoods-per-round` | 256 |
| `--work-repair-trials-per-round` | 256 |
| `--work-pool-lookahead-candidates` | 4 |
| `--work-cleanup-moves-per-round` | 64 |
| `--work-families-per-refresh` | 16 |

`--ncores` controls discovery workers. A verified reused catalog performs no discovery; reports record zero actual discovery workers even when more were requested.

## Optional salvage

Salvage is off by default. Enable `--salvage bounded`. The default exploratory ladder is −28, −30, −32, with at most 8 violating physical dimer edges and 4 incident oligo species per pool, measured cumulatively relative to strict −26. Each tier starts from an accepted incumbent and is freshly audited.

Controls: repeatable `--salvage-threshold`, `--salvage-max-stages` (up to 3), `--salvage-max-edges-per-pool`, `--salvage-max-oligos-per-pool`, and `--salvage-time-limit` (60 seconds per tier). Thresholds must be finite and strictly decreasing from −26. These are bounded experimental policies, not validated destructive-dimer cutoffs. Tighten exposure caps to restrict how many already selected oligos are put at additional modeled risk.

Strict remains the primary exported scheme unless `--primary-tier salvage-1` (or another enabled, successfully validated tier) is explicitly selected. Every completed tier and its measured tradeoffs remains in the output.

## Retained history, audits and reuse

The output retains generated sites, physical oligos, amplicon families, explored exact subsets, raw evidence, assessments, transitions, and stage snapshots. A failed candidate remains discoverable. “Not explored within budget” is distinct from “rejected by a measured constraint.” A compatible subset can remain selected after a conflicting member is removed, with the resulting allele-support loss recorded.

```sh
primalscheme3 panel-history --bundle /absolute/path/panel \
  --target REFERENCE_NAME --region 394:644 --limit 100 \
  --output /absolute/path/new-region-query
primalscheme3 panel-audit --bundle /absolute/path/panel \
  --output /absolute/path/new-audit
primalscheme3 panel-cache --bundle /absolute/path/panel \
  --output /absolute/path/new-discovery-cache
```

History regions are zero-based, half-open reference intervals; pool numbers are one-based. Use entity IDs for exact primer/family/subset histories, stage/profile filters, and offset/limit pagination. Queries and audits produce their own reproducibility receipts outside the source panel.

Add `--reuse-discovery /absolute/path/new-discovery-cache` to a compatible design command. Reuse verifies original target order/content, discovery settings, chemistry and scientific source/runtime fingerprints; incompatible caches fail closed. Copied origin history is explicitly distinguished from the new run's decisions. Reuse does not carry forward an old selection verdict.

All outputs preserve exact invocation, requested and resolved options, versions/runtime, stored paths, hashes/sizes, elapsed time and exit status. Original input paths are provenance origins; retained payloads remain usable after relocation.

### History database compatibility

New persistent histories use `primalscheme3.sqlite-history/v2`. Scientific record
IDs, compressed canonical payloads, record order, stage snapshots and exported
JSONL remain unchanged. The physical database stores entity and record relations
with integer foreign keys; canonical-ID SQL views preserve the reader contract.

Existing `primalscheme3.sqlite-history/v1` databases remain readable and
appendable in their original format. Opening, inspecting, auditing or reusing a
history never migrates it. A v1 cache origin may accompany a new v2 selection
history. Older tools that support only v1 must reject v2 rather than treating it
as v1; use a current reader for v2 outputs.

Both formats retain SQLite `synchronous=FULL`, `journal_mode=DELETE`, a 64 MiB
page cache and the default 1,000-record transaction batch. Stage completion,
explicit checkpoints and close commit durably. A killed writer can lose its
uncommitted batch; the committed prefix and completed snapshots remain valid.

For Python integrations, `SQLiteCoverageHistory(..., format_version=1)` explicitly
creates a new v1 history; the default is 2. The argument controls creation only:
an existing database always uses its recorded format. `rollback()` discards the
whole pending batch and reloads the durable prefix. Relation-write or commit
errors also roll back the pending batch. Callers must not continue using record
objects returned from a batch that was rolled back. The bounded entity-key cache
belongs to one writer and is cleared on rollback and close.


### Phase scheduling

`--phase-scheduling serial|reserved` defaults to `serial`. The opt-in `reserved`
policy (`initial-phase-reservations/v1`) reserves the initial cycle's remaining
time among full/normal seeds (20% combined), subset construction (40%), and
repair (40%). Disabled groups are omitted and weights normalized. Seed modes
and repair rounds divide their group's share equally. Unused time flows forward.
Within each initial repair round, preparation, cleanup, and exchange/refill get
20%, 20%, and 60%. Later starts use the remaining global budget in serial order.

Reservations are cooperative: an indivisible proposal batch or scientific
predicate may overrun a window. The hard global deadline and cancellation take
precedence at the next check, and the best validated incumbent remains retained.
This offers construction and exchange opportunities when operations fit; it does
not guarantee an improvement or completed repair. Serial remains the default
until same-catalog measurements support changing it.

Each stage records the resolved `scheduling_policy`, `phase_progress` (outcomes,
work/proposal deltas, before/after cursors and numeric objectives, accepted-repair
deltas, budgets, elapsed time and local overshoot), omitted
phases, the active phase at global stop, and global deadline overshoot. Outcomes
distinguish `work-cap`, `exhausted`, `no-eligible-work`, `phase-time-limit`, global
`time-limit`, and `cancelled`. Exhaustion refers to the named bounded phase or
family stream, never all possible variant subsets. `completed_seed_modes` means
returned bounded passes; `exhausted_seed_modes` separately identifies exhausted
seed streams. `fixed_work_completed` requires every requested bounded phase to
finish without any time or cancellation cutoff; reaching a deterministic work
cap is permitted. Phase lifecycle history events preserve status even without a
coverage improvement. These events introduce no synthetic primer entities.


### Optional designated-site screening

`--intended-product-policy concrete-designated-sites` changes only how potential
products at the intended locus are classified during supplied-row screening.
The default `exact-supported` requires a full-length exact F/R match on the same
row for an ordinary intended-product exemption.

The opt-in policy also permits a screened near-match when both oligos belong to
one complete selected configuration, hit their own designated projected full
footprints on the same target row in the correct orientation and order, and both
full footprints are observed A/C/G/T. Every endpoint must agree; a shifted repeat,
similar span, wrong orientation, other target, missing/N/IUPAC footprint, or a
terminal seed whose full primer extends beyond supplied sequence remains blocked.
Terminal screening still uses the configured k with at most one substitution;
it does not constrain mismatches outside that terminal segment.

These permitted near-matches add **zero coverage**. Coverage still requires
full-primer exact same-row support, with the original distinct-allele weighting.
The reference amplicon-size limits and positive inclusive screening product bound
are unchanged. Pair checks require a complete certificate from A or B; combining
one end from each configuration or borrowing a third cannot rescue a product.
The separate secondary-product rule continues to use exact-supported certificates.

Native profile/resolved options and stage constraints record
`intended_product_policy`; capability descriptors identify
`exact-supported/v1` and `concrete-designated-sites/v1`. Missing historical fields
mean `exact-supported`. Neither old artifacts nor the existing preset are silently
changed. This post-discovery control permits scientifically compatible catalog
reuse and changes the selection/validation profile identity.

Candidate/pair evidence and validation report `allowed_intended_products`
separately from `allowed_secondary_products`. Each newly permitted product records
its complete selected-site certificate, row footprints, observed templates, and
full-primer mismatch positions split into terminal/outside-terminal positions
(zero-based in synthesis orientation), with `coverage_credit: 0`. Validation's
`allowed_intended_product_count` counts contextual witnesses, which may repeat a
physical product across candidate/pair checks. Fresh audit reconstructs these
certificates from authoritative rows and checks the saved evidence and policy.
History inspection preserves those linked records.

This opt-in interpretation is not an amplification-probability or laboratory
specificity claim. Compare it explicitly against `exact-supported` using the same
catalog, search controls and source/runtime identity; a coverage increase would
not establish a global optimum or experimental performance.

### Live search observations

Allele panel creation writes `search-progress.jsonl` (schema
`primalscheme3.allele-search-progress/v1`). Each flushed record has a UTC timestamp,
stage, elapsed search time, active phase, incumbent ID, exact numeric objective
terms, per-target coverage, mean coverage, work counters and family cursors.
Search starts and incumbent improvements include configuration IDs and zero-based
pool IDs; the incumbent ID also identifies the canonical assignment tie breaker.
These IDs resolve against the retained configuration ledger. Subsequent records
reference that incumbent. Targets with no assessable classes have null coverage.

Phase transitions and improvements are recorded immediately. Cooperative ticks
emit a heartbeat after at least 30 seconds without another record. Elapsed and
latest-improvement times use the latest existing search-clock sample; expensive
indivisible operations can delay sampling and heartbeats. UTC timestamps record
emission time. Family counts describe explored proposal streams, not exhaustive
subset enumeration. In full-cloud mode, the unique-family count includes the
families visited by that mode's full-cloud seed stream.

All observations are **provisional oracle-checked**, not independently audited
scientific panels. Search completion can precede fresh validation or a fallback;
the published stage validation remains authoritative. The append-only log is
retained on failure and included in the final provenance output hashes. It is not
a checkpoint or resume format. Programmatic search callers can omit `observer`
to disable observations; the native pipeline enables them by default.
