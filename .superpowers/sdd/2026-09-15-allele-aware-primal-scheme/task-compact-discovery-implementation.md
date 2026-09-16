# Task 1: compact discovery history

## Implementation

Compact discovery is now the resolved native allele option default. The
low-level `build_variant_catalog` API keeps `history_detail="full"` as its
historical default; the allele pipeline passes the resolved compact/full mode
explicitly. Compact generation consumes the raw discovery records for the same
chemistry, mapping, site merge, and feasible-family enumeration. It retains
all concrete mapped sites, including chemically rejected or reference-rejected
sites, generated and accepting profile membership, observations, and family
members. It avoids row-origin payloads and per-site/family evidence,
assessments, and dispositions. Target/profile summaries report raw attempt,
concrete, mapped, accepted/rejected, nonsequence, and failure-reason counts.

Compact and full cache projections use the catalog's profile memberships. Full
projection continues to use indexed SQLite evidence bindings. Cache manifests
record detail and scope and reject explicit full-history reuse from compact
origins. Missing detail metadata remains the old full-history interpretation
only when the catalog is also an old artifact without a detail declaration.

Resolved options, optimizer history metadata, inspection scope, and capability
documents advertise compact/full policy. The CLI exposes
`--discovery-history compact|full` and registers
`panel-discovery-diagnose` after the diagnostic module became available.

## Focused validation

Commands run from the native worktree:

```text
./.venv/bin/python -m pytest -q tests/lge/test_compact_discovery.py
# 6 passed

./.venv/bin/python -m pytest -q tests/lge/test_allele_options.py
# 89 passed

./.venv/bin/python -m pytest -q tests/lge/test_allele_catalog_cache.py::test_compact_cache_projects_membership_without_reading_sqlite_history
# 1 passed

./.venv/bin/python -m pytest -q tests/lge/test_allele_inspection.py::test_history_inspection_discloses_compact_recorded_scope
# 1 passed

./.venv/bin/python -m pytest -q tests/lge/test_variant_discovery.py
# 28 passed
```

The complete cache module suite also passed in an isolated run (`12 passed`)
and the combined search/policy/cache run passed `65` tests before a separate
provenance identity drift failure caused by concurrent source edits in the
shared worktree. That failure is environmental and must be rerun after the
parallel diagnostic edits settle.

## Limitations and follow-up

Compact history deliberately cannot answer per-attempt chronology questions;
`panel-discovery-diagnose` reconstructs a bounded anchored slice with full
history and labels it reconstructed rather than historical. Catalog semantic
digests differ between compact and full modes because detail fields and
evidence bindings differ; scientific site/family IDs and profile membership
are the parity contract. Selector decision history and fresh final validation
remain unchanged.
