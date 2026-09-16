"""Variant-preserving profile union; independent of legacy PrimerCloud filtering."""
from collections import Counter, defaultdict
from copy import copy
from dataclasses import replace
from importlib.metadata import version
import numpy as np
from primalscheme3.core.parallel_discovery import discover
from primalscheme3.core.digestion import VARIANT_FREQUENCY_POLICY
from primalscheme3.core.variant_thermo import _FIELDS
from primalscheme3.panel.allele_coverage import canonical_observations
from primalscheme3.panel.coverage_types import (
    Target, OligoSite, CandidateFamily, VariantCatalog, canonical_json, semantic_id,
)


def discovery_profiles(config):
    """Apply each complete native preset, preserving user thermodynamic options."""
    result = {}
    for name, high in (('normal', False), ('high-gc', True)):
        profile = copy(config)
        profile.high_gc = high
        suffix = 'hgc' if high else 'default'
        for field in ('size_min', 'size_max', 'gc_min', 'gc_max'):
            part, bound = field.split('_')
            setattr(profile, 'primer_' + field, getattr(config, f'_primer_{part}_{suffix}_{bound}'))
        result[name] = profile
    return result


def variant_targets(msa_dict):
    """Freeze duplicate-content ordinals in caller input occurrence order."""
    counts, targets = Counter(), []
    for msa in msa_dict.values():
        rows = tuple(tuple(str(c).upper() for c in row) for row in msa.array)
        content = semantic_id('target-content', rows)
        occurrence = counts[content]
        counts[content] += 1
        identity = semantic_id('target-occurrence', (content, occurrence))
        mapping = tuple(None if x is None else int(x) for x in msa._mapping_array)
        reference = ''.join(c for c in rows[0] if c not in ('', '-'))
        reverse = tuple(i for i, x in enumerate(mapping) if x is not None)
        targets.append(Target(identity, int(msa.msa_index), occurrence,
                              tuple(semantic_id('row', (identity, i, row)) for i, row in enumerate(rows)),
                              rows, reference, len(reference), mapping,
                              reverse + ((reverse[-1]+1) if reverse else 0,)))
    return tuple(targets)


def _mapped_site(target, record):
    anchor, seq = record['index'], record['sequence']
    forward = record['direction'] == 'f'
    column = anchor-1 if forward else anchor
    mapped = target.mapping[column] if 0 <= column < len(target.mapping) else None
    interval, failure = None, 'unmapped-reference-anchor'
    if mapped is not None:
        interval = (mapped+1-len(seq), mapped+1) if forward else (mapped, mapped+len(seq))
        failure = None
        if interval[0] < 0 or interval[1] > target.reference_length:
            interval, failure = None, 'outside-reference'
    return OligoSite(target.id, seq, '+' if forward else '-', anchor, interval, failure)


HISTORY_ORIGIN_POLICY = 'grouped-row-origins/v1'
HISTORY_ASSESSMENT_POLICY = 'grouped-profile-conditions/v1'
_ORIGIN_FIELDS = frozenset(('row_index', 'alignment_footprint', 'support', 'row_frequency_weight'))


def discovery_condition_rows(query_result):
    """Yield lossless condition views with their persisted assessment context.

    Conditions are not standalone AssessmentRecords: source_assessment_id and
    origin_evidence_id resolve their immutable context and every row origin.
    Both live history queries and serialized query dictionaries are supported.
    Nonsequence groups have no chemistry conditions; their not-evaluated summary
    remains available in the original query's assessments.
    """
    def plain(record):
        return record.to_dict() if hasattr(record, 'to_dict') else record

    evidence = {item['id']: item for item in map(plain, query_result['evidence'])}
    for assessment in map(plain, query_result['assessments']):
        if assessment['check_name'] != 'profile-and-reference':
            continue
        for identity in assessment['evidence_ids']:
            item = evidence.get(identity)
            if item is None or item['measurement'] != 'row-enumeration':
                continue
            group = item['values']
            if group.get('origin_policy') != HISTORY_ORIGIN_POLICY:
                continue
            for name, condition in group['common']['checks'].items():
                yield {
                    'policy': HISTORY_ASSESSMENT_POLICY,
                    'source_assessment_id': assessment['id'],
                    'origin_evidence_id': identity,
                    'run_id': assessment['run_id'],
                    'stage_id': assessment['stage_id'],
                    'entity_ids': assessment['entity_ids'],
                    'pool': assessment['pool'],
                    'profile_id': assessment['profile_id'],
                    'context_digest': assessment['context_digest'],
                    'kernel_versions': assessment['kernel_versions'],
                    'check_name': name,
                    'reason': condition.get('reason', 'independent-profile-condition'),
                    'condition': condition,
                }


def _group_variant_records(records, target, row_content_digests):
    """Lossless factoring of identical outcomes, never a scientific filter.

    Every original record is exactly `common | origin` after removing the two
    catalog dependency aliases from origin. Unequal measurements, check outcomes,
    frequencies, lengths, reasons or enumeration boundaries cannot coalesce.
    """
    groups = {}
    for record in records:
        common = {k: v for k, v in record.items() if k not in _ORIGIN_FIELDS}
        key = canonical_json(common)
        group = groups.setdefault(key, {'origin_policy': HISTORY_ORIGIN_POLICY, 'common': common, 'origins': []})
        row_index = record['row_index']
        origin = {k: v for k, v in record.items() if k in _ORIGIN_FIELDS}
        origin.update(row_id=target.row_ids[row_index], row_content_digest=row_content_digests[row_index])
        group['origins'].append(origin)
    for key in sorted(groups, key=lambda k: (groups[k]['common']['direction'], groups[k]['common']['index'],
                                            groups[k]['common']['length'], groups[k]['common']['sequence'] or '', k)):
        group = groups[key]
        group['origins'] = tuple(sorted(group['origins'], key=lambda o: (o['row_index'], canonical_json(o))))
        yield group


def build_variant_catalog(targets, config, *, profiles=None, indexes=None, history=None,
                          length_mode=None, history_detail='full'):
    """Build diagnostic and selectable sites, then feasible hybrid anchor families.

    `targets` are authoritative Target records or an input-ordered MSA dictionary.
    `indexes` limits alignment anchors for controlled experiments. Discovery records
    a row-local chemistry-compatible length prefix, or explicit exhaustive lengths.
    Supply a persistent CoverageHistory for publication; absent history is in-memory.
    """
    length_mode = length_mode or getattr(config, 'discovery_length_mode', 'first-compatible')
    if length_mode not in ('first-compatible', 'all'):
        raise ValueError('discovery length mode must be first-compatible or all')
    if history_detail not in ('compact', 'full'):
        raise ValueError('history detail must be compact or full')
    if history is None:
        from primalscheme3.panel.coverage_history import CoverageHistory
        history = CoverageHistory(None, run_id='discovery')
    if isinstance(targets, dict):
        targets = variant_targets(targets)
    targets = tuple(targets)
    profiles = discovery_profiles(config) if profiles is None else profiles
    if not profiles:
        raise ValueError('at least one complete discovery profile is required')
    resolved = {'profiles': {name: {k: getattr(p, k) for k in _FIELDS}
                             for name, p in sorted(profiles.items())},
                'minimum_frequency': {name: p.min_base_freq for name, p in sorted(profiles.items())},
                'frequency_policy': VARIANT_FREQUENCY_POLICY,
                'history_origin_policy': HISTORY_ORIGIN_POLICY,
                'history_assessment_policy': HISTORY_ASSESSMENT_POLICY,
                'maximum_alignment_walk': {name: p.primer_max_walk for name, p in sorted(profiles.items())},
                'enumeration': 'row-anchor-profile-length/v2', 'discovery_length_mode': length_mode,
                'ambiguous_expansion_limit': 256,
                'amplicon_size_min': config.amplicon_size_min,
                'amplicon_size_max': config.amplicon_size_max,
                'history_detail': history_detail,
                'history_scope': ('compact-target-profile-summaries'
                                  if history_detail == 'compact' else 'full-attempt-origin-records'),
                'indexes': indexes}
    config.discovery_workers_by_target_profile = {}
    config.discovery_workers_by_msa = {}
    sites, dispositions = {}, {t.id: 'enumerated' for t in targets}
    stage = 'discovery'
    kernels = {name: version(name) for name in ('primalscheme3', 'primalschemers', 'primer3-py')}
    context_digest = semantic_id('discovery-context', resolved)
    for target in targets:
        # Bind every event to the authoritative aligned row without copying its
        # entire width into every retained length/anchor record. Hash once per
        # input row; the catalog preserves cells, including gap/missing symbols.
        row_content_digests = tuple(semantic_id('aligned-row-content/v1', row) for row in target.rows)
        observation_digest = semantic_id('target-row-corpus/v1', row_content_digests)
        for profile_id, profile in sorted(profiles.items()):
            history.emit(stage_id=stage, kind='enumeration-boundary', entity_ids=(target.id,),
                         changes={'profile_id': profile_id, 'policy': resolved, 'row_count': len(target.rows)})
            profile_indexes = indexes
            if indexes is not None:
                width = len(target.mapping)
                if any(i < 0 or i > width for side in indexes for i in side):
                    raise IndexError('alignment anchor outside target')
                profile_indexes = ([i for i in indexes[0] if i >= profile.primer_size_min],
                                   [i for i in indexes[1] if i + profile.primer_size_min <= width])
                history.emit(stage_id=stage, kind='profile-anchor-boundary', entity_ids=(target.id,),
                             changes={'profile_id': profile_id, 'requested': indexes,
                                      'enumerated': profile_indexes, 'reason': 'profile-minimum-length'})
            records, workers = discover(np.array(target.rows), profile, indexes=profile_indexes, variant_mode=True,
                                  discovery_length_mode=length_mode)
            config.discovery_workers_by_target_profile.setdefault(target.id, {})[profile_id] = workers
            config.discovery_workers_by_msa[str(target.source_msa_index)] = max(
                config.discovery_workers_by_target_profile[target.id].values())
            history.emit(stage_id=stage, kind='discovery-execution', entity_ids=(target.id,),
                         changes={'profile_id': profile_id, 'source_msa_index': target.source_msa_index,
                                  'requested_workers': profile.ncores, 'actual_workers': workers})
            if history_detail == 'compact':
                # Compact mode consumes raw discovery records directly.  It keeps
                # every concrete mapped site and profile membership, while
                # avoiding row-origin payloads, intrinsic evidence and per-site
                # assessments/events entirely.
                counts = Counter()
                failures = Counter()
                for record in records:
                    counts['raw_attempt_count'] += 1
                    sequence = record['sequence']
                    if not sequence:
                        counts['nonsequence_attempt_count'] += 1
                        failures[record['reason'] or 'unknown'] += 1
                        continue
                    counts['concrete_attempt_count'] += 1
                    site = _mapped_site(target, record)
                    if site is None:
                        counts['unmapped_attempt_count'] += 1
                        failures[record['reason'] or 'unmapped'] += 1
                        continue
                    previous = sites.get(site.id)
                    if site.mapping_failure is not None:
                        counts['mapped_failure_count'] += 1
                    chemistry_accepted = bool(record['accepted'])
                    accepted = chemistry_accepted and site.reference_footprint is not None
                    if accepted:
                        counts['accepted_concrete_attempt_count'] += 1
                    else:
                        counts['rejected_concrete_attempt_count'] += 1
                        reason = site.mapping_failure or record['reason']
                        failures[reason if reason and reason != 'pass' else 'rejected'] += 1
                    sites[site.id] = replace(
                        site,
                        # Match full mode: accepting membership records chemistry
                        # acceptance even when reference mapping later rejects the
                        # concrete site from family selection.
                        accepting_profile_ids=((profile_id,) if chemistry_accepted else ())
                        + (previous.accepting_profile_ids if previous else ()),
                        generated_profile_ids=(profile_id,)
                        + (previous.generated_profile_ids if previous else ()),
                        intrinsic_evidence_ids=(),
                    )
                history.emit(
                    stage_id=stage,
                    kind='compact-discovery-summary',
                    entity_ids=(target.id,),
                    changes={
                        'profile_id': profile_id,
                        'history_detail': history_detail,
                        **dict(counts),
                        'failed_attempts_by_reason': dict(sorted(failures.items())),
                    },
                )
                continue
            for group in _group_variant_records(records, target, row_content_digests):
                record = group['common']
                site = _mapped_site(target, record) if record['sequence'] else None
                payload_digest = semantic_id(HISTORY_ORIGIN_POLICY, group)
                entity_id = site.id if site else semantic_id('failed-footprint-group/v1',
                                                              (target.id, profile_id, payload_digest))
                evidence = history.record_evidence(entity_ids=(entity_id,), measurement='row-enumeration',
                    dependency_key={'target': target.id, 'observation_digest': observation_digest,
                                    'profile': profile_id, 'origin_payload_digest': payload_digest},
                    values=group, status='evaluated' if site else 'unknown')
                generated = history.emit(stage_id=stage, kind='generated', entity_ids=(entity_id,),
                                         changes={'profile_id': profile_id, 'origin_policy': HISTORY_ORIGIN_POLICY,
                                                  'origin_count': len(group['origins']),
                                                  'origin_evidence_id': evidence.id, 'reason': record['reason']})
                evidence_ids = [evidence.id]
                if site:
                    numeric = history.record_evidence(entity_ids=(site.species_id,), measurement='oligo-thermodynamics',
                        dependency_key={'sequence': site.sequence, 'kernels': kernels,
                                        'parameters': {k: getattr(profile, k) for k in
                                            ('mv_conc', 'dv_conc', 'dntp_conc', 'dna_conc',
                                             'use_annealing', 'primer_annealing_prop', 'primer_annealing_tempc')}},
                        values={name: check['value'] for name, check in record['checks'].items()
                                if name != 'minimum-frequency'}, status='evaluated')
                    evidence_ids.append(numeric.id)
                accepted = site is not None and record['accepted'] and site.reference_footprint is not None
                assessment = history.assess(stage_id=stage, entity_ids=(entity_id,), pool=None,
                    context_digest=context_digest, profile_id=profile_id,
                    kernel_versions=kernels, thresholds=resolved['profiles'][profile_id] | {
                        'minimum_frequency': profile.min_base_freq,
                        'assessment_policy': HISTORY_ASSESSMENT_POLICY},
                    check_name='profile-and-reference', outcome='pass' if accepted else ('fail' if site else 'not-evaluated'),
                    reason=(site.mapping_failure or record['reason']) if site else record['reason'],
                    evidence_ids=tuple(evidence_ids))
                if site:
                    previous = sites.get(site.id)
                    # Chemical acceptance is separate from reference eligibility.
                    chemistry = (profile_id,) if record['accepted'] else ()
                    sites[site.id] = replace(site,
                        accepting_profile_ids=chemistry + (previous.accepting_profile_ids if previous else ()),
                        generated_profile_ids=(profile_id,) + (previous.generated_profile_ids if previous else ()),
                        intrinsic_evidence_ids=tuple(evidence_ids) + (previous.intrinsic_evidence_ids if previous else ()))
                history.emit(stage_id=stage, kind='assessed', entity_ids=(entity_id,),
                             parent_event_ids=(generated.id,), assessment_ids=(assessment.id,),
                             changes={'selectable_in_profile': accepted,
                                      'assessment_policy': HISTORY_ASSESSMENT_POLICY,
                                      'failed_checks': tuple(name for name, check in record['checks'].items()
                                                             if check['outcome'] == 'fail')})
                dispositions[entity_id] = 'selectable' if accepted or dispositions.get(entity_id) == 'selectable' else ('rejected' if site else 'not-evaluated')
    families = []
    for target in targets:
        forward, reverse = defaultdict(list), defaultdict(list)
        for site in sites.values():
            if site.target_id == target.id and site.accepting_profile_ids and site.reference_footprint is not None:
                (forward if site.strand == '+' else reverse)[site.alignment_anchor].append(site)
        history.emit(stage_id=stage, kind='family-enumeration-boundary', entity_ids=(target.id,),
                     changes={'forward_anchor_count': len(forward), 'reverse_anchor_count': len(reverse),
                              'anchor_pair_count': len(forward)*len(reverse), 'subsets_materialized': 0})
        omitted = []
        for fa, fs in sorted(forward.items()):
            omitted_count = 0
            for ra, rs in sorted(reverse.items()):
                # Existence of a feasible member pair, not the longest envelope.
                feasible = fa <= ra and any(
                    f.reference_footprint[1] <= r.reference_footprint[0] and
                    config.amplicon_size_min <= r.reference_footprint[1]-f.reference_footprint[0] <= config.amplicon_size_max
                    for f in fs for r in rs)
                if not feasible:
                    omitted_count += 1
                    continue
                family = CandidateFamily(target.id, (fa, ra), tuple(s.id for s in fs), tuple(s.id for s in rs),
                    tuple((fp, rp) for fp in sorted({p for f in fs for p in f.accepting_profile_ids})
                          for rp in sorted({p for r in rs for p in r.accepting_profile_ids})))
                if history_detail == 'compact':
                    families.append(family)
                    continue
                generated = history.emit(stage_id=stage, kind='generated-family', entity_ids=(family.id,), changes={})
                evidence = history.record_evidence(entity_ids=(family.id,), measurement='family-geometry-existence',
                    dependency_key={'sites': sorted((s.id, s.reference_footprint) for s in fs+rs),
                                    'minimum': config.amplicon_size_min, 'maximum': config.amplicon_size_max},
                    values={'some_selection_feasible': feasible}, status='evaluated')
                history.emit(stage_id=stage, kind='family-geometry', entity_ids=(family.id,),
                             parent_event_ids=(generated.id,), changes={'feasible': feasible, 'evidence_id': evidence.id})
                dispositions[family.id] = 'selectable' if feasible else 'rejected-geometry'
                if feasible:
                    families.append(replace(family, geometry_evidence_ids=(evidence.id,)))
            omitted.append((fa, omitted_count))
        history.emit(stage_id=stage, kind='family-enumeration-complete', entity_ids=(target.id,),
                     changes={'unmaterialized_infeasible_anchor_pairs_by_forward_anchor': omitted,
                              'reason': 'no-member-pair-satisfies-reference-bounds',
                              'subsets_materialized': 0})
    catalog = VariantCatalog(targets, tuple(a for t in targets for a in canonical_observations(t)),
                             tuple(sites.values()), tuple(families), canonical_json(resolved),
                             tuple((t.source_msa_index, t.id) for t in targets))
    history.complete_stage(stage_id=stage, dispositions=dispositions, catalog_digest=catalog.semantic_digest,
                           ledger_digest=semantic_id('empty-ledger', ()))
    return catalog
