# Task 1 report: immutable catalogue and coordinate/support model

## Status

Implemented the Task 1 catalogue layer on base `8b92381855008cce3c66ddc4a72a2e50c8c881ac`. This task changes only the new catalogue/type modules, their synthetic test, and this report. It does not change the CLI, package version, installed package, or publication state.

Commit: `Add immutable coverage candidate catalogue` (this report is part of that commit).

## Interfaces

`primalscheme3.panel.coverage_types` provides four frozen records:

- `Target(id, source_msa_index, occurrence, row_ids, rows, reference_sequence, reference_length, mapping, ref_to_alignment)`
- `Candidate(id, target_id, full_interval, interior_interval, forward_oligos, reverse_oligos, joint_rows, unknown_rows, row_product_spans)`
- `Catalog(targets, candidates, semantic_digest, resolved_config_json, source_mapping, _native_pair_items)` with cached read-only `target_by_id`, `candidate_by_id`, and `native_pair_by_candidate_id` maps
- `Assignment(candidate_id, pool)`

`build_catalog(msa_dict, config) -> Catalog` consumes the current native `MSA`, `PrimerPair`, `FKmer`, and `RKmer` APIs. `write_catalog(catalog, path) -> dict` writes schema `primalscheme3.coverage-catalog/v1` and returns `path`, `semantic_digest`, actual `file_sha256`, and actual `file_size`.

Target identities hash ordered alignment cells only. FASTA labels, target names, source MSA numbering, paths, resolved runtime configuration, and worker counts do not affect scientific IDs or the semantic digest. Duplicate alignments receive deterministic occurrence suffixes and remain separate targets. File metadata retains resolved configuration and the source-MSA-to-target mapping outside the hashed scientific payload.

`Target.mapping[i]` is the reference coordinate represented by alignment column `i`, or `None` for a reference gap. `Target.ref_to_alignment[r]` is the alignment column containing reference base `r`; the final entry is the exclusive column after the last reference base. Alignment rows are tuples of cells so native terminal `""` observations remain distinct from internal `"-"` gaps.

Candidate geometry is zero-based half-open. `interior_interval` contains the shared exclusive forward-end and inclusive reverse-start anchors. `full_interval` uses the longest actual oligo on each side. Native `RKmer.seqs()` strings are stored directly as actual 5-prime-to-3-prime reverse oligos; support reconstruction reverse-complements the row binding site, so reverse oligos are not reversed twice.

The private `_anchored_matches` contract for Task 2 is:

```text
((actual_oligo, zero_based_ungapped_row_start, zero_based_ungapped_row_end), ...), uncertain
```

It walks each oligo length from the shared anchor through row bases, skips alignment gap cells, preserves terminal missing cells as uncertainty, and retains the exact interval for every matching cloud length. Non-N IUPAC bases use explicit compatible-base sets but produce uncertainty, never confirmed support. `N` and terminal missing cells are unknown. A candidate row is jointly supported only when both clouds have exact anchored matches. `row_product_spans` retains every distinct matched forward/reverse interval combination in source row order.

Support is cached once per target and `(orientation, anchor, canonical oligo cloud)`. Catalogue ID maps are constructed once as `MappingProxyType` instances rather than rebuilt on access. Native pair objects are mapped by candidate ID outside candidate equality and every scientific digest.

The gzip writer uses canonical sorted-key compact JSON and `mtime=0` with no stored filename. Targets and candidate records are canonical; the oligo table is globally deduplicated by sequence. The writer hashes the exact compressed bytes after serialization and returns their byte count, satisfying the scientific-output provenance requirement for the catalogue artifact.

## TDD evidence

Initial RED, before production modules existed:

```text
.venv/bin/python -m pytest tests/lge/test_coverage_catalog.py -q
ERROR tests/lge/test_coverage_catalog.py
ModuleNotFoundError: No module named 'primalscheme3.panel.coverage_catalog'
exit 2
```

Reverse terminal-missing RED after the first minimal implementation:

```text
.venv/bin/python -m pytest tests/lge/test_coverage_catalog.py::test_joint_support_requires_both_exact_anchors_and_unknown_is_diagnostic -q
1 failed: reverse terminal-missing row absent from unknown_rows
exit 1
```

Review-driven RED for scientific ordering and hot-loop contracts:

```text
.venv/bin/python -m pytest tests/lge/test_coverage_catalog.py -q
3 failed, 7 passed
```

The three failures demonstrated that source renumbering changed target order/digest, reference-gap/indel rows lacked joint support, and ID-map properties rebuilt mapping proxies.

Focused GREEN after row-anchor walking, scientific/metadata separation, per-cloud caching, and cached maps:

```text
.venv/bin/python -m pytest tests/lge/test_coverage_catalog.py -q
10 passed in 0.39s
exit 0
```

Existing LGE regression GREEN:

```text
.venv/bin/python -m pytest tests/lge -q
43 passed, 2 warnings in 4.97s
exit 0
```

The two warnings are existing Typer deprecations for Click 9 APIs. Python bytecode compilation also completed successfully. Ruff was not available in this task's `.venv` (`No module named ruff`), so no lint-pass claim is made.

## Coverage and risks

The synthetic tests cover hand-derived `(10, 50)` full and `(18, 42)` interior intervals, actual reverse-primer orientation, native object recreation, source permutation and renumbering, duplicate target occurrences, gapped reference length/mapping, gapped and deletion-allele anchored matches, F-only/R-only disjoint support, N/IUPAC/terminal uncertainty, mixed cloud-length row spans, frozen records/read-only cached maps, canonical oligo serialization, deterministic gzip bytes, and actual output SHA-256/size.

This layer records candidates with no jointly supported row so the diagnostic state is explicit; later compatibility/validation code must enforce the approved feasibility rule that a selectable candidate has at least one jointly supported observed row. Source metadata is deliberately excluded from the semantic digest, so two catalogues can have the same scientific digest but different compressed file hashes when source numbering or runtime metadata differs. The support model is bounded to the current linear first-reference native MSA API and does not attempt generalized importer compatibility.
