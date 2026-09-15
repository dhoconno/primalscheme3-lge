"""Exact, reference-independent observed allele coverage. Never screening support."""
from __future__ import annotations
from statistics import mean
from primalscheme3.panel.coverage_types import (
    Target, ObservedAllele, OligoSite, BindingSupport, VariantCatalog,
    SelectionConfiguration, SupportedProduct, AlleleAssignment, ConfigurationLedger,
    ClassCoverage, TargetCoverage, CoverageSummary,
)
from primalscheme3.panel.coverage_catalog import _sequence_match, _reverse_observed


def canonical_observations(target: Target) -> tuple[ObservedAllele, ...]:
    groups = {}
    for row_id, row in zip(target.row_ids, target.rows, strict=True):
        cells = tuple(x.upper() for x in row)
        groups.setdefault(cells, []).append(row_id)
    result = []
    for cells, aliases in groups.items():
        reverse = tuple(i for i, cell in enumerate(cells) if cell not in ('', '-'))
        by_alignment = {column: pos for pos, column in enumerate(reverse)}
        result.append(ObservedAllele(
            target.id, cells, tuple(by_alignment.get(i) for i in range(len(cells))), reverse,
            tuple(pos for pos, column in enumerate(reverse) if cells[column] in 'ACGT'),
            tuple(sorted(aliases)), len(aliases)))
    return tuple(sorted(result, key=lambda x: x.id))


def binding_support(site: OligoSite, allele: ObservedAllele, *, row_id: str | None = None) -> BindingSupport:
    if site.target_id != allele.target_id:
        raise ValueError('site and allele belong to different targets')
    if row_id is not None and row_id not in allele.row_ids:
        raise ValueError('row alias does not belong to allele')
    reverse = site.strand == '-'
    cursor = site.alignment_anchor if reverse else site.alignment_anchor - 1
    step = 1 if reverse else -1
    selected = []
    unavailable = False
    while len(selected) < len(site.sequence):
        if cursor < 0 or cursor >= len(allele.cells) or allele.cells[cursor] == '':
            unavailable = True
            break
        if allele.cells[cursor] != '-':
            selected.append(cursor)
        cursor += step
    footprint = tuple(sorted(selected))
    observed = tuple(allele.cells[i] for i in footprint)
    if reverse:
        observed = _reverse_observed(observed)
    match = 'unknown' if unavailable else _sequence_match(site.sequence, observed)
    status = 'confirmed' if match == 'exact' else match
    positions = tuple(allele.alignment_to_row[i] for i in footprint)
    row_footprint = (positions[0], positions[-1] + 1) if positions and not unavailable else None
    return BindingSupport(site.id, allele.id, row_id or allele.row_ids[0], status,
                          'unavailable-footprint' if unavailable else 'full-primer-' + match,
                          footprint, row_footprint)


def configuration_products(catalog: VariantCatalog, config: SelectionConfiguration) -> tuple[SupportedProduct, ...]:
    family = catalog.family_by_id[config.family_id]
    if config.target_id != family.target_id or config.anchor_pair != family.anchor_pair:
        raise ValueError('configuration family/target/anchors disagree')
    products = []
    sides = []
    for ids, allowed, strand in ((config.forward_site_ids, family.forward_site_ids, '+'),
                                 (config.reverse_site_ids, family.reverse_site_ids, '-')):
        if not set(ids) <= set(allowed):
            raise ValueError('configuration sites outside family')
        sites = tuple(catalog.site_by_id[x] for x in ids)
        if any(x.target_id != config.target_id or x.strand != strand for x in sites):
            raise ValueError('configuration site target/strand mismatch')
        sides.append(sites)
    for allele in catalog.observations:
        if allele.target_id != config.target_id:
            continue
        fs = [(s, binding_support(s, allele)) for s in sides[0]]
        rs = [(s, binding_support(s, allele)) for s in sides[1]]
        for f, fm in fs:
            for r, rm in rs:
                if fm.status == rm.status == 'confirmed' and fm.row_footprint[1] <= rm.row_footprint[0]:
                    products.append(SupportedProduct(allele.id, f.id, r.id, fm.row_footprint, rm.row_footprint))
    return tuple(products)


def configuration_coverage(catalog: VariantCatalog, config: SelectionConfiguration) -> dict[str, frozenset[int]]:
    covered = {a.id: set() for a in catalog.observations if a.target_id == config.target_id}
    for product in configuration_products(catalog, config):
        covered[product.allele_id].update(range(*product.interior))
    return {key: frozenset(value) for key, value in covered.items()}


def coverage_utility(fractions: tuple[float, ...], *, goal: float = .95) -> float | None:
    if not 0 < goal <= 1:
        raise ValueError('coverage goal must be in (0, 1]')
    if any(not 0 <= x <= 1 for x in fractions):
        raise ValueError('coverage fractions must be in [0, 1]')
    return mean(x - max(0, goal-x)**2/goal for x in fractions) if fractions else None


def _intervals(positions):
    result = []
    for x in sorted(positions):
        if result and result[-1][1] == x:
            result[-1] = (result[-1][0], x + 1)
        else:
            result.append((x, x + 1))
    return tuple(result)


def allele_summary(catalog: VariantCatalog, assignments: tuple[AlleleAssignment, ...],
                   ledger: ConfigurationLedger, *, goal: float) -> CoverageSummary:
    coverage_utility((), goal=goal)
    if ledger.catalog_digest != catalog.semantic_digest:
        raise ValueError('ledger belongs to a different catalog')
    keys = [(a.configuration_id, a.pool) for a in assignments]
    if len(set(keys)) != len(keys):
        raise ValueError('duplicate exact assignment')
    covered = {a.id: set() for a in catalog.observations}
    for assignment in assignments:
        for allele_id, bases in configuration_coverage(catalog, ledger.configuration_by_id[assignment.configuration_id]).items():
            covered[allele_id].update(bases)
    classes = []
    for allele in catalog.observations:
        observed = set(allele.observed_positions)
        recovered = observed & covered[allele.id]
        classes.append(ClassCoverage(allele.target_id, allele.id, len(observed), len(recovered),
                                     len(allele.row_to_alignment)-len(observed), allele.cells.count(''),
                                     allele.cells.count('-'), len(recovered)/len(observed) if observed else None,
                                     _intervals(recovered)))
    targets = []
    for target in catalog.targets:
        rows = [x for x in classes if x.target_id == target.id]
        fractions = [x.fraction for x in rows if x.fraction is not None]
        targets.append(TargetCoverage(target.id, mean(fractions) if fractions else None,
                                      'assessable' if fractions else 'unassessable-target',
                                      len(fractions), len(rows)-len(fractions)))
    fractions = tuple(x.fraction for x in targets if x.fraction is not None)
    return CoverageSummary(tuple(classes), tuple(targets), coverage_utility(fractions, goal=goal),
                           mean(fractions) if fractions else None, min(fractions) if fractions else None,
                           sum(x > goal for x in fractions), sum(x >= goal for x in fractions), goal,
                           'ok' if fractions else 'no-assessable-targets')
