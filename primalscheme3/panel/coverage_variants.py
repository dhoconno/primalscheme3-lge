"""Exact retained-site configurations and bounded, lazy subset neighborhoods.

This proposes candidates, not scientifically validated placements. Each distinct
site subset remains distinct even when its current allele coverage is identical.
"""
from __future__ import annotations
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from collections import deque
import json

from .coverage_types import SelectionConfiguration, ConfigurationLedger
from .allele_coverage import configuration_products, configuration_coverage


class ConfigurationIneligible(ValueError):
    def __init__(self, reason):
        self.reason = reason
        super().__init__(reason)


def _bounds(catalog, min_size, max_size):
    resolved = json.loads(catalog.resolved_config_json)
    minimum = resolved.get('amplicon_size_min', 1) if min_size is None else min_size
    maximum = resolved.get('amplicon_size_max') if max_size is None else max_size
    if minimum < 1 or (maximum is not None and maximum < minimum):
        raise ValueError('invalid amplicon size bounds')
    return minimum, maximum


def make_configuration(catalog, family_id, forward_site_ids, reverse_site_ids, *,
                       min_size=None, max_size=None, profile_ids=None):
    """Rebuild geometry and same-allele support from the exact selected sites.

    None bounds resolve from catalog configuration (1/unbounded for bare test
    catalogs). Profile membership means acceptance by any complete named profile.
    """
    minimum, maximum = _bounds(catalog, min_size, max_size)
    family = catalog.family_by_id[family_id]
    fs, rs = tuple(sorted(set(forward_site_ids))), tuple(sorted(set(reverse_site_ids)))
    if not fs or not rs:
        raise ConfigurationIneligible('empty-side')
    sites = []
    for ids, allowed, strand, anchor in ((fs, family.forward_site_ids, '+', family.anchor_pair[0]),
                                        (rs, family.reverse_site_ids, '-', family.anchor_pair[1])):
        if not set(ids) <= set(allowed):
            raise ConfigurationIneligible('family-membership')
        for site_id in ids:
            site = catalog.site_by_id[site_id]
            if site.target_id != family.target_id or site.strand != strand or site.alignment_anchor != anchor:
                raise ConfigurationIneligible('site-anchor-membership')
            if not site.accepting_profile_ids or (profile_ids is not None and not set(profile_ids).intersection(site.accepting_profile_ids)):
                raise ConfigurationIneligible('no-accepting-profile')
            if site.mapping_failure or site.reference_footprint is None:
                raise ConfigurationIneligible('unavailable-reference-mapping')
            a, b = site.reference_footprint
            if not 0 <= a < b <= catalog.target_by_id[family.target_id].reference_length:
                raise ConfigurationIneligible('out-of-reference-bounds')
            target = catalog.target_by_id[family.target_id]
            column = anchor-1 if strand == '+' else anchor
            mapped = target.mapping[column] if 0 <= column < len(target.mapping) else None
            if mapped is None:
                raise ConfigurationIneligible('unavailable-reference-mapping')
            expected = (mapped+1-len(site.sequence), mapped+1) if strand == '+' else (mapped, mapped+len(site.sequence))
            if site.reference_footprint != expected:
                raise ConfigurationIneligible('inconsistent-reference-mapping')
            sites.append(site)
    fsites, rsites = sites[:len(fs)], sites[len(fs):]
    if family.anchor_pair[0] >= family.anchor_pair[1] or max(s.reference_footprint[1] for s in fsites) > min(s.reference_footprint[0] for s in rsites):
        raise ConfigurationIneligible('invalid-primer-geometry')
    envelope = (min(s.reference_footprint[0] for s in sites), max(s.reference_footprint[1] for s in sites))
    size = envelope[1] - envelope[0]
    if size < minimum or (maximum is not None and size > maximum):
        raise ConfigurationIneligible('amplicon-size')
    config = SelectionConfiguration(family.target_id, family.id, fs, rs, envelope, family.anchor_pair)
    products = configuration_products(catalog, config)
    masks = {a.id: set(a.observed_positions) for a in catalog.observations if a.target_id == family.target_id}
    if not any(any(x in masks[p.allele_id] for x in range(*p.interior)) for p in products):
        raise ConfigurationIneligible('no-useful-joint-support')
    return replace(config, supported_products=products)


@dataclass(frozen=True)
class SubsetLimits:
    beam_width: int = 16
    expansion_limit: int = 256

    def __post_init__(self):
        if self.beam_width < 1 or self.expansion_limit < 1:
            raise ValueError('subset work limits must be positive')


@dataclass(frozen=True)
class ProposalContext:
    pool: int | None = None
    context_digest: str = ''
    min_size: int | None = None
    max_size: int | None = None
    profile_ids: tuple[str, ...] | None = None


@dataclass(frozen=True)
class ProposalWitness:
    id: str
    site_ids: tuple[str, ...]


@dataclass(frozen=True)
class ProposalBatch:
    ledger: ConfigurationLedger
    new_ids: tuple[str, ...]
    expansion_count: int
    omitted_materialized_ids: tuple[str, ...]
    unexpanded_frontier_count: int
    truncated: bool
    stop_cause: str
    dispositions: object = field(default_factory=dict)
    # Frontier count counts pending lazy streams, not unenumerated subsets.
    frontier_count_kind: str = 'pending-neighborhood-streams'

    def __post_init__(self):
        object.__setattr__(self, 'dispositions', MappingProxyType(dict(self.dispositions)))


def propose_configurations(catalog, ledger, family_id, *, parent_ids=(), witnesses=(),
                           context=ProposalContext(), stage='strict', limits=SubsetLimits(), history=None):
    """Budget before materialization; batch ledger append once per call.

    Use a family-local ledger in large searches and freeze the global registry at
    checkpoints. No powerset is constructed. Beam omissions are explicit and may
    be reconsidered in a later call with those configurations as parents.
    """
    if ledger.catalog_digest != catalog.semantic_digest:
        raise ValueError('ledger belongs to a different catalog')
    family = catalog.family_by_id[family_id]
    kwargs = dict(min_size=context.min_size, max_size=context.max_size, profile_ids=context.profile_ids)
    selectable = set(catalog.selectable_site_ids)
    if context.profile_ids is not None:
        selectable = {s for s in selectable if set(catalog.site_by_id[s].accepting_profile_ids).intersection(context.profile_ids)}
    fs = tuple(sorted(set(family.forward_site_ids) & selectable))
    rs = tuple(sorted(set(family.reverse_site_ids) & selectable))
    parents = tuple(ledger.configuration_by_id[x] for x in sorted(set(parent_ids)))
    if any(p.family_id != family_id for p in parents):
        raise ValueError('parent belongs to a different family')
    witnesses = tuple(sorted(witnesses, key=lambda w: w.id))

    def deletions(parent):
        f, r = parent.forward_site_ids, parent.reverse_site_ids
        for witness in witnesses:
            for site in sorted(set(witness.site_ids) & set(f+r)):
                yield tuple(x for x in f if x != site), tuple(x for x in r if x != site), parent, (witness.id,), 'witness-deletion'

    full_configuration = None

    def seeds():
        for parent in parents:
            yield from deletions(parent)
        for parent in parents:
            yield from neighbors(parent)
        yield fs, rs, None, (), 'full-family'
        for f in fs:
            for r in rs:
                yield (f,), (r,), full_configuration, (), 'minimal-pair'

    def neighbors(parent):
        yield from deletions(parent)
        sides = (parent.forward_site_ids, parent.reverse_site_ids)
        for i, allowed in enumerate((fs, rs)):
            selected = sides[i]
            for old in selected:
                edited = list(sides); edited[i] = tuple(x for x in selected if x != old)
                yield *edited, parent, (), 'deletion'
            for new in allowed:
                if new in selected:
                    continue
                edited = list(sides); edited[i] = tuple(sorted(selected+(new,)))
                yield *edited, parent, (), 'addition'
                for old in selected:
                    edited = list(sides); edited[i] = tuple(sorted(tuple(x for x in selected if x != old)+(new,)))
                    yield *edited, parent, (), 'swap'

    streams = deque([(None, iter(seeds()))])
    seen = set()
    records, omitted, dispositions = [], [], {}
    existing = ledger.configuration_by_id
    event_ids = {entity: event.id for event in history.events for entity in event.entity_ids} if history else {}
    expanded = 0
    beam_used = 0
    while streams and expanded < limits.expansion_limit:
        try:
            f, r, parent, witness_ids, reason = next(streams[0][1])
        except StopIteration:
            streams.popleft(); beam_used = 0
            continue
        key = (tuple(sorted(f)), tuple(sorted(r)))
        if key in seen:
            continue
        seen.add(key); expanded += 1
        try:
            config = make_configuration(catalog, family_id, f, r, **kwargs)
        except ConfigurationIneligible as exc:
            if history is not None:
                history.emit(stage_id=stage, kind='configuration-ineligible', entity_ids=(family_id,),
                             changes={'selected_forward_site_ids': f, 'selected_reverse_site_ids': r,
                                      'reason': exc.reason, 'proposal_reason': reason,
                                      'parent_configuration_id': parent.id if parent else None,
                                      'witness_ids': witness_ids})
            continue
        if config.forward_site_ids == fs and config.reverse_site_ids == rs:
            full_configuration = config
        if config.id in existing:
            # An existing seed still supplies a neighborhood in this call.
            if beam_used < limits.beam_width:
                streams.append((config.id, iter(neighbors(config)))); beam_used += 1
            continue
        records.append(config)
        before = configuration_coverage(catalog, parent) if parent else {}
        after = configuration_coverage(catalog, config)
        delta = {a.id: len(after[a.id] & set(a.observed_positions)) - len(before.get(a.id, frozenset()) & set(a.observed_positions))
                 for a in catalog.observations if a.target_id == family.target_id}
        disposition = 'proposed'
        if beam_used >= limits.beam_width:
            omitted.append(config.id); disposition = 'not-explored/work-limit'
        else:
            streams.append((config.id, iter(neighbors(config)))); beam_used += 1
        dispositions[config.id] = disposition
        if history is not None:
            old = set(parent.forward_site_ids+parent.reverse_site_ids) if parent else set()
            new = set(f+r)
            parent_events = (event_ids[parent.id],) if parent and parent.id in event_ids else ()
            event = history.emit(stage_id=stage, kind='configuration-proposed', entity_ids=(config.id,),
                parent_event_ids=parent_events, changes={'parent_configuration_id': parent.id if parent else None,
                'added_site_ids': tuple(sorted(new-old)), 'removed_site_ids': tuple(sorted(old-new)),
                'witness_ids': witness_ids, 'proposal_reason': reason, 'coverage_delta': delta,
                'supported_allele_ids': tuple(sorted({p.allele_id for p in config.supported_products})),
                'gained_supported_allele_ids': tuple(sorted({p.allele_id for p in config.supported_products} - {p.allele_id for p in configuration_products(catalog, parent)})) if parent else tuple(sorted({p.allele_id for p in config.supported_products})),
                'lost_supported_allele_ids': tuple(sorted({p.allele_id for p in configuration_products(catalog, parent)} - {p.allele_id for p in config.supported_products})) if parent else (),
                'retained_profile_memberships': {s: catalog.site_by_id[s].accepting_profile_ids for s in sorted(new)},
                'disposition': disposition, 'pool': context.pool, 'context_digest': context.context_digest})
            event_ids[config.id] = event.id
    # Pending streams are unexpanded work; never claim an exact unseen powerset count.
    for configuration_id, _ in streams:
        if configuration_id is not None and configuration_id not in omitted:
            omitted.append(configuration_id)
            dispositions[configuration_id] = 'not-explored/work-limit'
            if history is not None:
                history.emit(stage_id=stage, kind='configuration-not-explored', entity_ids=(configuration_id,),
                             parent_event_ids=(event_ids[configuration_id],) if configuration_id in event_ids else (),
                             changes={'reason': 'expansion-limit', 'disposition': 'not-explored/work-limit'})
    pending = len(streams)
    truncated = bool(pending or omitted)
    cause = 'expansion-limit' if pending else ('beam-limit' if omitted else 'exhausted')
    if history is not None:
        history.emit(stage_id=stage, kind='subset-proposal-summary', entity_ids=(family_id,),
                     changes={'expansion_count': expanded, 'omitted_materialized_ids': tuple(omitted),
                              'pending_neighborhood_streams': pending, 'truncated': truncated, 'stop_cause': cause})
    return ProposalBatch(ledger.appended(*records), tuple(c.id for c in records), expanded,
                         tuple(omitted), pending, truncated, cause, dispositions)
