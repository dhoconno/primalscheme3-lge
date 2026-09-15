"""Variant evidence precedes all cloud-level filters."""
from dataclasses import replace
import numpy as np
import pytest
from primalscheme3.core.config import Config
from primalscheme3.core.parallel_discovery import discover
from primalscheme3.panel.coverage_types import Target

SHARED = 'CAACGGCGGACTTTATTGTATCTCC'
HIGH = 'AATCGGCCGACTGTACGCGA'
NORMAL = 'CATGCTCATTAGGTATATCTTTCAATAAGTTGCA'


def target(rows):
    width = max(map(len, rows))
    rows = tuple(tuple(row.rjust(width, '-')) for row in rows)
    ref = ''.join(x for x in rows[0] if x not in ('', '-'))
    mapping, position = [], 0
    for x in rows[0]:
        mapping.append(position if x not in ('', '-') else None)
        position += x not in ('', '-')
    reverse = tuple(i for i, x in enumerate(mapping) if x is not None)
    return Target('t', 0, 0, tuple(str(i) for i in range(len(rows))), rows,
                  ref, len(ref), tuple(mapping), reverse + (width,))


def test_row_variants_survive_failing_sibling_and_preserve_failed_evidence():
    t = target([SHARED, 'A' * len(SHARED), HIGH])
    records, _ = discover(np.array(t.rows), Config(), indexes=([len(t.rows[0])], []), variant_mode=True, discovery_length_mode='all')
    assert any(r['sequence'] == SHARED and r['accepted'] for r in records)
    assert any(r['sequence'] == 'A' * len(SHARED) and not r['accepted'] for r in records)
    assert all('checks' in r for r in records if r['sequence'])


def test_union_profiles_merge_shared_sites_without_rectangular_chemistry():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    t = target([NORMAL, SHARED, HIGH])
    catalog = build_variant_catalog((t,), Config(), length_mode='all', indexes=([len(t.rows[0])], []))
    sites = {s.sequence: s for s in catalog.sites}
    assert sites[NORMAL].accepting_profile_ids == ('normal',)
    assert sites[HIGH].accepting_profile_ids == ('high-gc',)
    assert sites[SHARED].accepting_profile_ids == ('high-gc', 'normal')


def test_gap_walk_uses_original_anchor_and_missing_is_diagnostic():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    cells = tuple(SHARED[:8] + '-' + SHARED[8:])
    t = target([''.join(cells), ''.join(cells)])
    t = replace(t, rows=(cells, ('',) * len(cells)))
    catalog = build_variant_catalog((t,), Config(), length_mode='all', indexes=([len(cells)], []))
    site = next(s for s in catalog.sites if s.sequence == SHARED)
    assert site.alignment_anchor == len(cells)
    assert site.reference_footprint == (0, len(SHARED))


def test_hybrid_family_retains_short_member_when_longer_mapping_fails():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    # The longer normal-only site falls outside the shorter first reference.
    t = target([SHARED + 'ATATATATAT' + HIGH, NORMAL + 'ATATATATAT' + HIGH])
    end = len(t.rows[0]) - len(HIGH) - 10
    config = Config(amplicon_size=54, amplicon_size_min=50, amplicon_size_max=65)
    catalog = build_variant_catalog((t,), config, length_mode='all', indexes=([end], [end + 10]))
    short = next(s for s in catalog.sites if s.sequence == SHARED and s.strand == '+')
    long = next(s for s in catalog.sites if s.sequence == NORMAL and s.strand == '+')
    assert short.reference_footprint is not None
    assert long.reference_footprint is None
    assert any(short.id in f.forward_site_ids for f in catalog.families)
    assert all(long.id not in f.forward_site_ids for f in catalog.families)


def test_discovery_worker_count_does_not_change_records():
    t = target([SHARED * 4])
    one = discover(np.array(t.rows), Config(ncores=1), variant_mode=True, discovery_length_mode='all')[0]
    two = discover(np.array(t.rows), Config(ncores=2), variant_mode=True, discovery_length_mode='all')[0]
    assert one == two


def test_self_dimer_measurement_agrees_with_native_boolean_at_boundary():
    from primalschemers import do_pool_interact
    from primalscheme3.core.variant_thermo import variant_measurements
    for seq in (SHARED, HIGH, NORMAL, 'ACGTACGTACGTACGTACGT'):
        score = variant_measurements(seq, Config())['self-dimer']['value']
        for cutoff in (score - 1e-6, score, score + 1e-6, -26):
            assert (score <= cutoff) == do_pool_interact([seq.encode()], [seq.encode()], cutoff)


def test_failed_sites_and_aliases_have_history_and_no_sibling_frequency_veto(tmp_path):
    from primalscheme3.panel.coverage_history import CoverageHistory
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    history = CoverageHistory(tmp_path, run_id='test')
    t = target([SHARED, SHARED, 'A'*len(SHARED)])
    catalog = build_variant_catalog((t,), Config(min_base_freq=.5), length_mode='all', indexes=([len(SHARED)], []), history=history)
    shared = next(s for s in catalog.sites if s.sequence == SHARED)
    failed = next(s for s in catalog.sites if s.sequence == 'A'*len(SHARED))
    assert shared.id in catalog.selectable_site_ids
    assert failed.id not in catalog.selectable_site_ids
    assert history.snapshots[-1].completeness == 'complete'
    for site in (shared, failed):
        events = [e for e in history.events if site.id in e.entity_ids]
        assert events[0].kind == 'generated'
        assert any(e.kind == 'assessed' for e in events)
    origins = [e.values['origins'] for e in history.evidence
               if e.measurement == 'row-enumeration' and shared.id in e.entity_ids]
    assert any({o['row_id'] for o in group} == {'0', '1'} for group in origins)


def test_high_gc_short_anchor_does_not_fail_normal_index_validation():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    t = target(['GCGTACGCGTACGCGTA'])
    catalog = build_variant_catalog((t,), Config(), length_mode='all', indexes=([17], [0]))
    assert catalog.sites
    assert all(s.generated_profile_ids == ('high-gc',) for s in catalog.sites)


def test_cloud_internal_dimer_does_not_erase_either_chemical_variant():
    from primalschemers import do_pool_interact
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    a, b = 'TCGCTGTTAGCATTTCAGTCTCCT', 'AGGTGTAGGAAAGTGGAGGAGACC'
    assert do_pool_interact([a.encode()], [b.encode()], -26)
    catalog = build_variant_catalog((target([a, b]),), Config(), length_mode='all', indexes=([24], []))
    selected = {s.sequence for s in catalog.sites if s.id in catalog.selectable_site_ids}
    assert {a, b} <= selected


def test_ambiguous_expansion_has_unknown_support_and_omission_boundary():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.allele_coverage import binding_support
    from primalscheme3.panel.coverage_history import CoverageHistory
    t = target([SHARED[:-1] + 'N', 'N' * len(SHARED)])
    history = CoverageHistory(None, run_id='ambiguity')
    catalog = build_variant_catalog((t,), Config(), length_mode='all', indexes=([len(SHARED)], []), history=history)
    site = next(s for s in catalog.sites if s.sequence == SHARED)
    assert all(binding_support(site, a).status == 'unknown' for a in catalog.observations)
    assert any(e.values.get('common', {}).get('reason') == 'ambiguous-expansion-limit' for e in history.evidence)


def test_union_catalog_identity_does_not_depend_on_worker_count_or_profile_order():
    from copy import copy
    from primalscheme3.panel.coverage_discovery import build_variant_catalog, discovery_profiles
    t = target([SHARED * 4])
    profiles = discovery_profiles(Config())
    # One length per profile keeps the history fixture small; >128 tasks still
    # uses two real spawned workers, with actual thermodynamic measurements.
    for p in profiles.values():
        p.primer_size_max = p.primer_size_min
    one = build_variant_catalog((t,), Config(), profiles=profiles)
    reverse = {name: copy(p) for name, p in reversed(list(profiles.items()))}
    for p in reverse.values():
        p.ncores = 2
    two = build_variant_catalog((t,), Config(), profiles=reverse)
    assert one.semantic_digest == two.semantic_digest
    assert one.sites == two.sites


def test_cross_profile_rectangle_is_not_a_complete_accepting_profile():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    # Its 32 nt length requires normal, its 62.5% GC requires high-gc.
    sequence = 'GCGCAT' * 5 + 'AA'
    assert len(sequence) == 32
    catalog = build_variant_catalog((target([sequence]),), Config(), length_mode='all', indexes=([32], []))
    site = next(s for s in catalog.sites if s.sequence == sequence)
    assert site.accepting_profile_ids == ()
    assert site.id not in catalog.selectable_site_ids


def test_first_compatible_is_row_local_and_records_unexamined_longer_lengths():
    t = target([SHARED, 'A' * len(SHARED)])
    records, _ = discover(np.array(t.rows), Config(), indexes=([len(SHARED)], []), variant_mode=True)
    first_row = [r for r in records if r['row_index'] == 0]
    passing = [r for r in first_row if r['sequence'] and r['accepted']]
    assert len(passing) == 1
    boundary = next(r for r in first_row if r['reason'] == 'first-compatible-length-stop')
    assert boundary['unexamined_length_interval'] == (passing[0]['length'] + 1, 37)
    assert all(r['length'] <= passing[0]['length'] for r in first_row if r['sequence'])
    # A failing sibling still evaluates its entire available footprint.
    assert any(r['row_index'] == 1 and r['length'] == len(SHARED) and r['sequence'] for r in records)


def test_explicit_exhaustive_length_mode_retains_later_compatible_lengths():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    t = target([SHARED])
    first = build_variant_catalog((t,), Config(), indexes=([len(SHARED)], []))
    all_lengths = build_variant_catalog((t,), Config(), indexes=([len(SHARED)], []), length_mode='all')
    assert len(first.sites) < len(all_lengths.sites)
    assert '"discovery_length_mode":"first-compatible"' in first.resolved_config_json
    assert '"discovery_length_mode":"all"' in all_lengths.resolved_config_json


def test_family_geometry_recomputes_envelope_without_longest_valid_member():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    t = target([NORMAL + 'ATATATATAT' + HIGH, SHARED + 'ATATATATAT' + HIGH])
    end = len(NORMAL)
    catalog = build_variant_catalog((t,), Config(amplicon_size=55, amplicon_size_min=50, amplicon_size_max=60),
                                    length_mode='all', indexes=([end], [end + 10]))
    short = next(s for s in catalog.sites if s.sequence == SHARED and s.strand == '+')
    long = next(s for s in catalog.sites if s.sequence == NORMAL and s.strand == '+')
    assert short.reference_footprint is not None and long.reference_footprint is not None
    family = next(f for f in catalog.families if short.id in f.forward_site_ids and long.id in f.forward_site_ids)
    assert ('normal', 'high-gc') in family.discovery_profile_combinations
    reverse = next(catalog.site_by_id[x] for x in family.reverse_site_ids
                   if len(catalog.site_by_id[x].sequence) == len(HIGH))
    assert reverse.reference_footprint[1] - long.reference_footprint[0] > 60
    assert 50 <= reverse.reference_footprint[1] - short.reference_footprint[0] <= 60


def test_numeric_chemistry_evidence_is_shared_across_sites_and_profiles():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.coverage_history import CoverageHistory
    history = CoverageHistory(None, run_id='shared')
    catalog = build_variant_catalog((target([SHARED * 2]),), Config(), length_mode='all',
                                    indexes=([len(SHARED), len(SHARED) * 2], []), history=history)
    shared_sites = [s for s in catalog.sites if s.sequence == SHARED]
    evidence = [e for e in history.evidence if e.measurement == 'oligo-thermodynamics'
                and e.dependency_key['sequence'] == SHARED]
    assert len(shared_sites) == 2
    assert len(evidence) == 1
    assert all(evidence[0].id in s.intrinsic_evidence_ids for s in shared_sites)
    assert evidence[0].entity_ids == (shared_sites[0].species_id,)


def frequency_config(minimum=.75):
    config = Config(min_base_freq=minimum)
    config.primer_size_min = config.primer_size_max = len(SHARED)
    return config


def test_frequency_excludes_terminal_missing_rows_at_fixed_footprint():
    array = np.array([list(SHARED), [''] * len(SHARED)])
    records, _ = discover(array, frequency_config(), indexes=([len(SHARED)], []), variant_mode=True)
    shared = next(r for r in records if r['sequence'] == SHARED)
    assert shared['accepted']
    assert shared['frequency'] == 1
    assert shared['frequency_evidence']['numerator'] == 1
    assert shared['frequency_evidence']['denominator'] == 1
    assert shared['frequency_evidence']['excluded_terminal_row_indexes'] == (1,)


def test_frequency_apportions_ambiguous_row_at_fixed_footprint():
    array = np.array([list(SHARED[:-1] + 'N')])
    records, _ = discover(array, frequency_config(), indexes=([len(SHARED)], []), variant_mode=True)
    alternatives = [r for r in records if r['sequence']]
    assert len(alternatives) == 4
    assert all(r['frequency'] == .25 and not r['accepted'] for r in alternatives)
    assert all(r['support'] == 'unknown' for r in alternatives)
    assert sum(r['row_frequency_weight'] for r in alternatives) == 1


def test_frequency_shared_alternatives_sum_fractional_row_weights():
    array = np.array([list(SHARED), list(SHARED[:-1] + 'N'), [''] * len(SHARED)])
    records, _ = discover(array, frequency_config(), indexes=([len(SHARED)], []), variant_mode=True)
    shared_records = [r for r in records if r['sequence'] == SHARED]
    assert all(r['frequency'] == .625 and not r['accepted'] for r in shared_records)
    for r in shared_records:
        evidence = r['frequency_evidence']
        assert evidence['numerator'] == 1.25 and evidence['denominator'] == 2
        assert evidence['row_weights'] == ((0, 1.0), (1, .25))
        assert evidence['policy'] == 'observed-only-anchored-length-row-mass/v2'


def test_failed_prefixes_and_exhaustive_lengths_do_not_multiply_row_mass():
    array = np.array([list(SHARED), list(SHARED[:-1] + 'N')])
    by_mode = {}
    for mode in ('first-compatible', 'all'):
        records, _ = discover(array, Config(min_base_freq=.75), indexes=([len(SHARED)], []),
                              variant_mode=True, discovery_length_mode=mode)
        by_mode[mode] = {r['sequence']: r['frequency_evidence'] for r in records if r['sequence']}
        for row in (0, 1):
            for length in {r['length'] for r in records if r['sequence']}:
                mass = sum(r['row_frequency_weight'] for r in records
                           if r['sequence'] and r['row_index'] == row and r['length'] == length)
                assert 0 <= mass <= 1
    assert len(by_mode['all']) > len(by_mode['first-compatible'])
    for sequence, evidence in by_mode['first-compatible'].items():
        assert evidence == by_mode['all'][sequence]
    # Exhaustive later compatible sequence has full reconstructed support even
    # though first-compatible discovery stopped that row at a shorter length.
    assert by_mode['all'][SHARED]['numerator'] == 1.25
    assert by_mode['all'][SHARED]['denominator'] == 2


def test_capped_ambiguous_row_retains_denominator_without_invented_support():
    array = np.array([list(SHARED), list('N' * len(SHARED))])
    records, _ = discover(array, frequency_config(), indexes=([len(SHARED)], []), variant_mode=True)
    shared = next(r for r in records if r['sequence'] == SHARED)
    assert shared['frequency'] == .5 and not shared['accepted']
    assert shared['frequency_evidence']['denominator'] == 2
    assert shared['frequency_evidence']['unallocated_row_indexes'] == (1,)
    assert any(r['reason'] == 'ambiguous-expansion-limit' for r in records)


def test_frequency_denominator_is_footprint_specific_and_reverse_aware():
    from primalscheme3.core.seq_functions import reverse_complement
    row = list(reverse_complement(SHARED))
    partial = row.copy()
    partial[-1] = ''
    config = Config(min_base_freq=.75)
    config.primer_size_min = len(SHARED) - 1
    config.primer_size_max = len(SHARED)
    records, _ = discover(np.array([row, partial]), config, indexes=([], [0]),
                          variant_mode=True, discovery_length_mode='all')
    full = next(r for r in records if r['sequence'] == SHARED)
    assert full['frequency'] == 1 and full['accepted']
    assert full['frequency_evidence']['denominator'] == 1
    assert full['frequency_evidence']['excluded_terminal_row_indexes'] == (1,)
    short = next(r for r in records if r['sequence'] and r['length'] == len(SHARED) - 1)
    assert short['frequency'] == 1
    assert short['frequency_evidence']['denominator'] == 2


def test_row_dependencies_use_stable_catalog_resolvable_content_digests():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.coverage_history import CoverageHistory
    from primalscheme3.panel.coverage_types import canonical_json, semantic_id

    def discover_prefix(row, run_id):
        history = CoverageHistory(None, run_id=run_id)
        catalog = build_variant_catalog((target([row]),), Config(),
                                        indexes=([len(SHARED)], []), history=history)
        records = tuple(e for e in history.evidence if e.measurement == 'row-enumeration')
        assert records
        for evidence in records:
            dependency = evidence.dependency_key
            assert 'row' not in dependency
            stored = catalog.target_by_id[dependency['target']]
            for origin in evidence.values['origins']:
                row_index = stored.row_ids.index(origin['row_id'])
                assert origin['row_content_digest'] == semantic_id('aligned-row-content/v1', stored.rows[row_index])
        return records

    original = discover_prefix(SHARED + 'A' * 100, 'first')
    repeated = discover_prefix(SHARED + 'A' * 100, 'repeat')
    distal_change = discover_prefix(SHARED + 'A' * 99 + 'C', 'changed')
    long_row = discover_prefix(SHARED + 'A' * 10000, 'long')
    assert original == repeated
    # The changed distal base does not alter measured prefix values, but must
    # invalidate the source dependency even with the same source/row aliases.
    assert [e.values['common'] for e in original] == [e.values['common'] for e in distal_change]
    assert {e.id for e in original}.isdisjoint(e.id for e in distal_change)
    assert {o['row_content_digest'] for e in original for o in e.values['origins']}.isdisjoint(
        o['row_content_digest'] for e in distal_change for o in e.values['origins'])
    assert [len(canonical_json(e.dependency_key)) for e in original] == [
        len(canonical_json(e.dependency_key)) for e in long_row]


def test_grouped_occurrence_evidence_expands_to_every_primitive_outcome():
    import json
    from primalscheme3.panel.coverage_discovery import build_variant_catalog, discovery_profiles
    from primalscheme3.panel.coverage_history import CoverageHistory
    from primalscheme3.panel.coverage_types import canonical_json
    gapped = SHARED[:8] + '-' + SHARED[8:]
    t = target([gapped, gapped, gapped[:-1] + 'N', 'A' * len(gapped), 'N' * len(gapped), gapped])
    t = replace(t, rows=t.rows[:-1] + (('',) * len(gapped),))
    indexes = ([len(gapped)], [0])
    for mode in ('first-compatible', 'all'):
        history = CoverageHistory(None, run_id='grouped')
        config = Config(min_base_freq=.25)
        catalog = build_variant_catalog((t,), config, indexes=indexes, length_mode=mode, history=history)
        assert json.loads(catalog.resolved_config_json)['history_origin_policy'] == 'grouped-row-origins/v1'
        for name, profile in discovery_profiles(config).items():
            primitive, _ = discover(np.array(t.rows), profile, indexes=indexes, variant_mode=True, discovery_length_mode=mode)
            groups = [e for e in history.evidence if e.measurement == 'row-enumeration'
                      and e.dependency_key['profile'] == name]
            expanded = []
            for group in groups:
                assert group.values['origin_policy'] == 'grouped-row-origins/v1'
                for origin in group.values['origins']:
                    assert origin['row_id'] == t.row_ids[origin['row_index']]
                    raw_origin = {k: v for k, v in origin.items() if k not in ('row_id', 'row_content_digest')}
                    expanded.append(dict(group.values['common']) | raw_origin)
            assert sorted(map(canonical_json, expanded)) == sorted(map(canonical_json, primitive))
            assert len(groups) < len(primitive)
        assert any(len(e.values['origins']) >= 2 for e in history.evidence if e.measurement == 'row-enumeration')


def test_grouping_retains_different_common_frequency_decisions():
    from copy import deepcopy
    from primalscheme3.panel.coverage_discovery import _group_variant_records
    from primalscheme3.panel.coverage_types import semantic_id
    t = target([SHARED, SHARED])
    records, _ = discover(np.array(t.rows), frequency_config(), indexes=([len(SHARED)], []), variant_mode=True)
    concrete = [deepcopy(r) for r in records if r['sequence'] == SHARED]
    assert len(concrete) == 2
    concrete[1]['checks']['minimum-frequency']['outcome'] = 'fail'
    concrete[1]['accepted'] = False
    digests = tuple(semantic_id('aligned-row-content/v1', row) for row in t.rows)
    groups = list(_group_variant_records(concrete, t, digests))
    assert len(groups) == 2
    assert {group['common']['accepted'] for group in groups} == {True, False}


def test_compacted_catalog_retains_precompaction_scientific_projection():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.coverage_types import semantic_id
    t = target([NORMAL + 'ATATATATAT' + HIGH, SHARED + 'ATATATATAT' + HIGH, SHARED + 'ATATATATAT' + HIGH])
    catalog = build_variant_catalog((t,), Config(amplicon_size=55, amplicon_size_min=50, amplicon_size_max=60),
                                    length_mode='all', indexes=([len(NORMAL)], [len(NORMAL) + 10]))
    projection = {'targets': catalog.targets, 'observations': catalog.observations,
                  'sites': [dict(s.to_dict(), intrinsic_evidence_ids=[]) for s in catalog.sites],
                  'families': [dict(f.to_dict(), geometry_evidence_ids=[]) for f in catalog.families]}
    # Captured from the preceding ungrouped implementation on this exact fixture.
    assert semantic_id('discovery-science-projection', projection) == 'discovery-science-projection-c1ecfb4b672475dbc9fcf94d257f83cc237303a525a3ba9c3f381ce9381c2a2f'


def test_duplicate_rows_expand_origins_without_multiplying_history_decisions():
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.coverage_history import CoverageHistory
    counts = []
    for copies in (1, 12):
        history = CoverageHistory(None, run_id='duplicates')
        build_variant_catalog((target([SHARED] * copies),), Config(), indexes=([len(SHARED)], []), history=history)
        occurrences = [e for e in history.evidence if e.measurement == 'row-enumeration']
        assert occurrences and all(len(e.values['origins']) == copies for e in occurrences)
        assert all(len({o['row_id'] for o in e.values['origins']}) == copies for e in occurrences)
        counts.append((len(occurrences), len(history.assessments), len(history.events)))
    assert counts[0] == counts[1]


@pytest.mark.parametrize("length_mode", ("first-compatible", "all"))
def test_compact_condition_queries_preserve_every_primitive_check_after_reload(tmp_path, length_mode):
    from primalscheme3.panel.coverage_discovery import (
        build_variant_catalog, discovery_profiles, discovery_condition_rows,
    )
    from primalscheme3.panel.coverage_history import SQLiteCoverageHistory
    from primalscheme3.panel.coverage_types import canonical_json
    t = target([SHARED, SHARED, SHARED[:-1] + 'N', 'A' * len(SHARED), 'N' * len(SHARED)])
    config = Config(min_base_freq=.25)
    path = tmp_path / 'history'
    history = SQLiteCoverageHistory(path, run_id='conditions')
    build_variant_catalog((t,), config, indexes=([len(SHARED)], []), history=history, length_mode=length_mode)
    groups = [e for e in history.evidence if e.measurement == 'row-enumeration']
    assert len(history.assessments) == len(groups)
    query = history.query(stage_id='discovery')
    before = list(discovery_condition_rows(query))
    import json
    serialized = json.loads(canonical_json({k: [r.to_dict() for r in v] for k, v in query.items()}))
    assert canonical_json(before) == canonical_json(list(discovery_condition_rows(serialized)))
    history.close()
    history = SQLiteCoverageHistory(path, run_id='conditions')
    assert before == list(discovery_condition_rows(history.query(stage_id='discovery')))
    for name, profile in discovery_profiles(config).items():
        raw, _ = discover(np.array(t.rows), profile, indexes=([len(SHARED)], []), variant_mode=True, discovery_length_mode=length_mode)
        query = history.query(stage_id='discovery', profile_id=name)
        evidence = {e.id: e for e in query['evidence']}
        assessments = {a.id: a for a in query['assessments']}
        expanded = []
        for condition in discovery_condition_rows(query):
            assessment = assessments[condition['source_assessment_id']]
            assert condition['stage_id'] == 'discovery'
            assert condition['profile_id'] == name
            assert condition['kernel_versions'] == dict(assessment.kernel_versions)
            assert condition['context_digest'] == assessment.context_digest
            assert assessment.thresholds['minimum_frequency'] == .25
            group = evidence[condition['origin_evidence_id']].values
            for origin in group['origins']:
                expanded.append((origin['row_index'], group['common']['sequence'],
                                 condition['check_name'], condition['condition']))
        expected = [(r['row_index'], r['sequence'], k, v) for r in raw for k, v in r['checks'].items()]
        assert sorted(map(canonical_json, expanded)) == sorted(map(canonical_json, expected))
    assert any(c['condition']['outcome'] == 'not-evaluated' for c in before)
    assert any(c['condition']['outcome'] == 'fail' for c in before)
    assert any(a.outcome == 'not-evaluated' for a in history.assessments)
    history.close()


def test_discovery_records_actual_workers_without_changing_scientific_identity():
    import json
    from primalscheme3.panel.coverage_discovery import build_variant_catalog
    from primalscheme3.panel.coverage_history import CoverageHistory
    t = target([SHARED])
    config = Config(selection_algorithm='allele-coverage', amplicon_size_min=360, ncores=4)
    history = CoverageHistory(None, run_id='workers')
    catalog = build_variant_catalog((t,), config, indexes=([len(SHARED)], []), history=history)
    assert config.discovery_workers_by_target_profile == {'t': {'high-gc': 1, 'normal': 1}}
    assert config.discovery_workers_by_msa == {'0': 1}
    assert config.discovery_core_count == 1
    execution = [e for e in history.events if e.kind == 'discovery-execution']
    assert len(execution) == 2
    assert all(e.changes['requested_workers'] == 4 and e.changes['actual_workers'] == 1 for e in execution)
    resolved = json.loads(catalog.resolved_config_json)
    assert resolved['history_assessment_policy'] == 'grouped-profile-conditions/v1'
    other = Config(selection_algorithm='allele-coverage', amplicon_size_min=360, ncores=1)
    same = build_variant_catalog((t,), other, indexes=([len(SHARED)], []))
    assert catalog.semantic_digest == same.semantic_digest
