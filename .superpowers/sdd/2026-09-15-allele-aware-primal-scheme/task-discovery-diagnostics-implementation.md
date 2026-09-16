# Task 2: on-demand anchored discovery diagnostics

## Interface

`primalscheme3.panel.discovery_diagnostics.diagnose_discovery` accepts

```python
diagnose_discovery(
    *, bundle, output, family_id=None, site_id=None, argv=None
) -> dict
```

Exactly one entity ID is required. A family replay passes its original forward
and reverse anchor pair; a site replay passes only that site's strand and
anchor. The original target row corpus and the resolved source profile set are
used, with `history_detail="full"` passed to Task 1's
`build_variant_catalog` interface. No saved command string is executed.

## Output and validation

The output is a new directory containing `discovery-diagnostic.json`, a
gzip-serialized replay catalog, a SQLite full replay history, copied source
payload bytes, and `provenance.json`. The receipt uses the existing native
provenance finalizer and includes argv, resolved replay/config fields, source
and runtime identities, input/output hashes and sizes, wall time, status, and
stderr. Failed preflight, unknown entity, replay errors, and scientific
membership drift receive failure receipts; an existing output directory or an
output path inside the source bundle is rejected before mutation.

The report labels the result `reconstructed-anchored-discovery` and
`historicalTrace=false`. Scientific parity compares requested site/profile or
family/anchor/member/profile membership while deliberately ignoring evidence
IDs and event chronology. A parity mismatch is a failed diagnostic rather
than an explanation of an optimizer decision.

Source preflight requires a completed stable native panel provenance, verifies
stored input bytes and required optimizer/config/catalog/authoritative-target
artifact descriptors, reparses the stored inputs and compares authoritative
targets, and checks current source/runtime/native-kernel identities. It never
writes to the source bundle.

## Validation evidence

Focused command:

```text
./.venv/bin/python -m pytest -q tests/lge/test_discovery_diagnostics.py
```

Result: `6 passed`.

Static checks:

```text
ruff check primalscheme3/panel/discovery_diagnostics.py tests/lge/test_discovery_diagnostics.py
```

Result: clean.

The focused tests cover successful family and site anchor selection, exact-one
entity validation, unknown IDs, source-byte preservation, scientific
membership drift, builder failure, and failure provenance. End-to-end CLI
wiring is intentionally left to the parent integration task.
