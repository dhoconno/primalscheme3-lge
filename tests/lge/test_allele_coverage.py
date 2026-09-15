import json
import pytest
from primalscheme3.panel import coverage_types as types
from allele_fixtures import target, literal_coverage


def api():
    from primalscheme3.panel import allele_coverage
    return allele_coverage


def test_reference_independent_fraction():
    assert literal_coverage(set(range(10)), [(2, 8)]) == 0.6
    assert literal_coverage(set(), []) is None


def test_duplicate_observations_and_masks():
    a = api()
    t = target('AACCGG', 'AACCGG', ('', 'A', 'C', 'C', 'G', 'G'))
    obs = a.canonical_observations(t)
    assert sorted(x.multiplicity for x in obs) == [1, 2]
    original = a.canonical_observations(target('AACCGG'))[0]
    duplicate = next(x for x in obs if x.multiplicity == 2)
    assert original.id == duplicate.id
    assert duplicate == types.ObservedAllele.from_dict(json.loads(json.dumps(duplicate.to_dict())))


def test_row_maps_and_support_walk():
    a = api()
    t = target('AA--CCNGTT', ('', 'A', 'T', '-', 'C', 'C', 'N', 'G', 'T', ''))
    obs = next(x for x in a.canonical_observations(t) if x.cells[0] == '')
    assert obs.row_to_alignment == (1, 2, 4, 5, 6, 7, 8)
    assert obs.observed_positions == (0, 1, 2, 3, 5, 6)
    f = types.OligoSite(t.id, 'ATC', '+', 5, (0, 3))
    s = a.binding_support(f, obs)
    assert s.status == 'confirmed'
    assert s.alignment_footprint == (1, 2, 4)
    assert s.row_footprint == (0, 3)
    assert a.binding_support(types.OligoSite(t.id, 'AAC', '+', 5, (0, 3)), obs).status == 'mismatch'
    assert a.binding_support(types.OligoSite(t.id, 'AG', '+', 8, (4, 6)), obs).status == 'unknown'
    assert a.binding_support(types.OligoSite(t.id, 'AA', '+', 2, (0, 2)), obs).status == 'unknown'


def catalog_and_configs():
    a = api()
    t = target('AACCGGTT', 'TTCCGGAA')
    f = types.OligoSite(t.id, 'AA', '+', 2, (0, 2))
    wrong_r = types.OligoSite(t.id, 'TT', '-', 6, (6, 8))
    right_r = types.OligoSite(t.id, 'AA', '-', 6, (6, 8))
    family = types.CandidateFamily(t.id, (2, 6), (f.id,), (wrong_r.id, right_r.id))
    catalog = types.VariantCatalog((t,), a.canonical_observations(t), (f, wrong_r, right_r), (family,))
    bad = types.SelectionConfiguration(t.id, family.id, (f.id,), (wrong_r.id,), (0, 8), (2, 6))
    good = types.SelectionConfiguration(t.id, family.id, (f.id,), (right_r.id,), (0, 8), (2, 6))
    return catalog, bad, good


def test_same_class_joint_support_and_split_pools():
    a = api()
    catalog, bad, good = catalog_and_configs()
    assert all(not x for x in a.configuration_coverage(catalog, bad).values())
    coverage = a.configuration_coverage(catalog, good)
    assert sorted(map(len, coverage.values())) == [0, 4]
    ledger = types.ConfigurationLedger(catalog.semantic_digest, (bad, good))
    assignments = (types.AlleleAssignment(bad.id, 0, 'strict'), types.AlleleAssignment(good.id, 1, 'strict'))
    summary = a.allele_summary(catalog, assignments, ledger, goal=.95)
    assert summary.mean_coverage == .25
    assert summary.targets[0].fraction == .25
    assert catalog == types.VariantCatalog.from_dict(catalog.to_dict())
    assert ledger == types.ConfigurationLedger.from_dict(ledger.to_dict())
    assert bad.id != good.id


def test_union_insertions_and_unknown_denominator():
    a = api()
    t = target('AA--CCNGTT', 'AATTCCNGTT')
    f = types.OligoSite(t.id, 'AA', '+', 2, (0, 2))
    r = types.OligoSite(t.id, 'AA', '-', 8, (6, 8))
    fam = types.CandidateFamily(t.id, (2, 8), (f.id,), (r.id,))
    c = types.SelectionConfiguration(t.id, fam.id, (f.id,), (r.id,), (0, 8), (2, 8))
    cat = types.VariantCatalog((t,), a.canonical_observations(t), (f, r), (fam,))
    ledger = types.ConfigurationLedger(cat.semantic_digest, (c,))
    summary = a.allele_summary(cat, (types.AlleleAssignment(c.id, 0, 'strict'), types.AlleleAssignment(c.id, 1, 'strict')), ledger, goal=.95)
    assert sorted(x.fraction for x in summary.classes) == [3/7, 5/9]
    for obs in cat.observations:
        positions = a.configuration_coverage(cat, c)[obs.id]
        expected = literal_coverage(obs.observed_positions, [(2, len(obs.row_to_alignment)-2)])
        assert len(positions & set(obs.observed_positions))/len(obs.observed_positions) == expected


def test_soft_utility_and_empty_outcomes():
    a = api()
    assert a.coverage_utility((.5, .8)) > a.coverage_utility((0., .8))
    assert a.coverage_utility((.96, .8)) > a.coverage_utility((.95, .8))
    assert a.coverage_utility(()) is None
    with pytest.raises(ValueError): a.coverage_utility((.2,), goal=0)
    t = target('NNN')
    cat = types.VariantCatalog((t,), a.canonical_observations(t), (), ())
    s = a.allele_summary(cat, (), types.ConfigurationLedger(cat.semantic_digest, ()), goal=.95)
    assert s.status == 'no-assessable-targets'
    assert s.mean_coverage is None and s.classes[0].fraction is None
    assert s.targets[0].status == 'unassessable-target'


def test_catalog_cached_interfaces_and_profile_aliases():
    catalog, _, _ = catalog_and_configs()
    assert catalog.site_by_id is catalog.site_by_id
    assert catalog.family_by_id is catalog.family_by_id
    site = types.OligoSite('t', 'acgt', '+', 4, (0, 4), accepting_profile_ids=('normal',))
    assert site.accepting_profiles == ('normal',)
    cat = types.VariantCatalog((), (), (site,), ())
    assert cat.selectable_site_ids == (site.id,)
    assert site.id == types.OligoSite('t', 'ACGT', '+', 4, (0, 4), accepting_profile_ids=('high-gc',)).id
    relocated = types.VariantCatalog((), (), (site,), (), source_mapping=((0, '/other/path'),))
    assert cat.semantic_digest == relocated.semantic_digest


def test_reference_metric_version_and_duplicate_alias_invariance():
    a = api()
    catalog, _, good = catalog_and_configs()
    s = a.allele_summary(catalog, (), types.ConfigurationLedger(catalog.semantic_digest, (good,)), goal=.95)
    assert s.metric_id == 'observed-allele-primer-trimmed/v1'
    one = a.canonical_observations(target('AACCGG'))[0]
    two = a.canonical_observations(target('AACCGG', 'AACCGG', aliases=('elsewhere', 'alias')))[0]
    assert one.id == two.id


def test_partial_overlapping_interiors_form_union_and_threshold_counts():
    a = api()
    t = target('AAAAAAAAAA')
    fs = (types.OligoSite(t.id,'A','+',1,(0,1)), types.OligoSite(t.id,'A','+',3,(2,3)))
    rs = (types.OligoSite(t.id,'T','-',6,(6,7)), types.OligoSite(t.id,'T','-',9,(9,10)))
    families = tuple(types.CandidateFamily(t.id,(f.alignment_anchor,r.alignment_anchor),(f.id,),(r.id,)) for f,r in zip(fs,rs))
    cs = tuple(types.SelectionConfiguration(t.id,fam.id,(f.id,),(r.id,),(f.reference_footprint[0],r.reference_footprint[1]),fam.anchor_pair) for fam,f,r in zip(families,fs,rs))
    cat = types.VariantCatalog((t,),a.canonical_observations(t),fs+rs,families)
    ledger = types.ConfigurationLedger(cat.semantic_digest,cs)
    assigned = tuple(types.AlleleAssignment(c.id,i,'strict') for i,c in enumerate(cs))
    s = a.allele_summary(cat,assigned,ledger,goal=.8)
    assert s.mean_coverage == .8
    assert s.classes[0].covered_intervals == ((1,9),)
    assert s.at_least_goal == 1 and s.strictly_above_goal == 0
    assert a.allele_summary(cat,assigned[:1],ledger,goal=.8).mean_coverage == .5
    with pytest.raises(ValueError,match='duplicate'):
        a.allele_summary(cat, assigned+assigned[:1],ledger,goal=.8)


@pytest.mark.parametrize('field,value', [
    ('observed_positions', [2, 3, 4, 5]),
    ('alignment_to_row', [7, 6, 5, 4, 3, 2, 1, 0]),
    ('row_to_alignment', [7, 6, 5, 4, 3, 2, 1, 0]),
    ('multiplicity', 0), ('row_ids', []),
])
def test_observation_rejects_incoherent_derived_payload(field, value):
    catalog, _, _ = catalog_and_configs()
    payload = catalog.to_dict()
    payload['observations'][0][field] = value
    with pytest.raises(ValueError):
        types.VariantCatalog.from_dict(payload)


def test_nested_records_reject_stale_semantic_ids():
    catalog, _, config = catalog_and_configs()
    payload = catalog.to_dict()
    payload['sites'][0]['id'] = 'forged'
    with pytest.raises(ValueError, match='identity'):
        types.VariantCatalog.from_dict(payload)
    ledger = types.ConfigurationLedger(catalog.semantic_digest, (config,)).to_dict()
    ledger['configurations'][0]['id'] = 'forged'
    with pytest.raises(ValueError, match='identity'):
        types.ConfigurationLedger.from_dict(ledger)
