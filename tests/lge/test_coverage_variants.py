import pytest
from allele_fixtures import target
from primalscheme3.panel.coverage_types import OligoSite, CandidateFamily, VariantCatalog, ConfigurationLedger
from primalscheme3.panel.allele_coverage import canonical_observations, configuration_coverage
from primalscheme3.panel.coverage_history import CoverageHistory


def api():
    from primalscheme3.panel import coverage_variants
    return coverage_variants


def fixture():
    t = target('AACCCCTTAA', 'AACCCCGGAA', 'AACCCCCCAA')
    f = OligoSite(t.id, 'AA', '+', 2, (0, 2), accepting_profile_ids=('normal',))
    rs = tuple(OligoSite(t.id, s, '-', 6, (6, 8), accepting_profile_ids=(p,))
               for s, p in [('AA', 'normal'), ('CC', 'high-gc'), ('GG', 'normal')])
    fam = CandidateFamily(t.id, (2, 6), (f.id,), tuple(r.id for r in rs))
    return VariantCatalog((t,), canonical_observations(t), (f,)+rs, (fam,)), fam, f, rs


def test_subset_rebuilds_support_and_retains_partial_classes():
    v = api(); cat, fam, f, rs = fixture()
    full = v.make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs))
    sub = v.make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs[:2]))
    assert full.full_interval == sub.full_interval == (0, 8)
    assert sorted(map(len, configuration_coverage(cat, sub).values())) == [0, 4, 4]
    assert len(full.supported_products) == 3 and len(sub.supported_products) == 2
    assert all(p.interior == (2, 6) for p in sub.supported_products)
    assert sub.reverse_site_ids == tuple(sorted(r.id for r in rs[:2]))
    with pytest.raises(v.ConfigurationIneligible, match='profile'):
        v.make_configuration(cat, fam.id, (f.id,), (rs[1].id,), profile_ids=('normal',))


def test_rejects_no_joint_support_and_recomputes_selected_envelope():
    v = api(); cat, fam, f, rs = fixture()
    from dataclasses import replace
    wrong = replace(f, sequence='TT')
    long = replace(rs[0], sequence='TTAA', reference_footprint=(6, 10))
    fam = replace(fam, forward_site_ids=(f.id, wrong.id), reverse_site_ids=(rs[0].id, long.id))
    cat = replace(cat, sites=cat.sites+(wrong, long), families=(fam,))
    with pytest.raises(v.ConfigurationIneligible, match='joint-support'):
        v.make_configuration(cat, fam.id, (wrong.id,), (rs[0].id,))
    with pytest.raises(v.ConfigurationIneligible, match='size'):
        v.make_configuration(cat, fam.id, (f.id,), (rs[0].id, long.id), max_size=8)
    assert v.make_configuration(cat, fam.id, (f.id,), (rs[0].id,), max_size=8).full_interval == (0, 8)
    with pytest.raises(v.ConfigurationIneligible, match='membership'):
        v.make_configuration(cat, fam.id, (f.id,), (rs[1].id,))


def test_witness_deletion_records_lineage_and_loss_and_keeps_alternatives():
    v = api(); cat, fam, f, rs = fixture()
    full = v.make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs))
    ledger = ConfigurationLedger(cat.semantic_digest, (full,)); history = CoverageHistory(run_id='test')
    batch = v.propose_configurations(cat, ledger, fam.id, parent_ids=(full.id,),
        witnesses=(v.ProposalWitness('conflicting-third', (rs[2].id,)),), history=history)
    desired = v.make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs[:2]))
    assert desired.id in batch.new_ids
    event = next(e for e in history.events if desired.id in e.entity_ids and e.kind == 'configuration-proposed')
    assert event.changes['parent_configuration_id'] == full.id
    assert event.changes['removed_site_ids'] == (rs[2].id,)
    assert 'conflicting-third' in event.changes['witness_ids']
    assert sorted(event.changes['coverage_delta'].values()) == [-4, 0, 0]
    assert ledger.configurations == (full,)
    assert len([c for c in batch.ledger.configurations if len(c.reverse_site_ids) == 1]) == 3


def test_bounded_generation_reports_materialized_omissions():
    v = api(); cat, fam, f, rs = fixture()
    batch = v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, ()), fam.id,
        limits=v.SubsetLimits(beam_width=1, expansion_limit=3), history=CoverageHistory(run_id='bounded'))
    assert batch.expansion_count == 3
    assert len(batch.new_ids) <= 3
    assert batch.truncated and batch.unexpanded_frontier_count > 0
    assert batch.omitted_materialized_ids
    assert all(batch.dispositions[x] == 'not-explored/work-limit' for x in batch.omitted_materialized_ids)
    again = v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, ()), fam.id,
        limits=v.SubsetLimits(beam_width=1, expansion_limit=3))
    assert again.new_ids == batch.new_ids


def test_defaults_resolve_catalog_bounds_and_unexpanded_new_ids_are_explicit():
    from dataclasses import replace
    v = api(); cat, fam, f, rs = fixture()
    cat = replace(cat, resolved_config_json='{"amplicon_size_min": 9, "amplicon_size_max": 10}')
    with pytest.raises(v.ConfigurationIneligible, match='size'):
        v.make_configuration(cat, fam.id, (f.id,), (rs[0].id,))
    h = CoverageHistory(run_id='one')
    batch = v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, ()), fam.id,
        context=v.ProposalContext(min_size=1), limits=v.SubsetLimits(expansion_limit=1), history=h)
    assert batch.omitted_materialized_ids == batch.new_ids
    assert any(e.kind == 'configuration-not-explored' for e in h.events)


def test_reference_touching_primers_allow_observed_insertion_interior():
    v = api()
    t = target('AA--TT', 'AACCTT')
    f = OligoSite(t.id, 'AA', '+', 2, (0, 2), accepting_profile_ids=('normal',))
    r = OligoSite(t.id, 'AA', '-', 4, (2, 4), accepting_profile_ids=('normal',))
    fam = CandidateFamily(t.id, (2, 4), (f.id,), (r.id,))
    cat = VariantCatalog((t,), canonical_observations(t), (f, r), (fam,))
    c = v.make_configuration(cat, fam.id, (f.id,), (r.id,))
    assert sorted(map(len, configuration_coverage(cat, c).values())) == [0, 2]


def test_large_family_is_lazy_and_keeps_equal_coverage_conflict_signatures():
    v = api()
    t = target('A'*80)
    fs = tuple(OligoSite(t.id, 'A'*n, '+', 30, (30-n, 30), accepting_profile_ids=('normal',)) for n in range(1, 26))
    rs = tuple(OligoSite(t.id, 'T'*n, '-', 40, (40, 40+n), accepting_profile_ids=('normal',)) for n in range(1, 26))
    fam = CandidateFamily(t.id, (30, 40), tuple(s.id for s in fs), tuple(s.id for s in rs))
    cat = VariantCatalog((t,), canonical_observations(t), fs+rs, (fam,))
    batch = v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, ()), fam.id)
    assert batch.expansion_count == 256 and len(batch.new_ids) == 256
    assert batch.truncated and batch.stop_cause == 'expansion-limit'
    assert len({tuple(configuration_coverage(cat, c).values()) for c in batch.ledger.configurations}) == 1
    assert len({c.forward_site_ids+c.reverse_site_ids for c in batch.ledger.configurations}) == 256


def test_minimal_pairs_link_to_full_and_existing_nonfull_parent_expands():
    v = api(); cat, fam, f, rs = fixture()
    full = v.make_configuration(cat, fam.id, (f.id,), tuple(r.id for r in rs))
    h = CoverageHistory(run_id='lineage')
    batch = v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, ()), fam.id, history=h)
    minimum = next(e for e in h.events if e.kind == 'configuration-proposed' and e.changes['proposal_reason'] == 'minimal-pair')
    assert minimum.changes['parent_configuration_id'] == full.id
    assert minimum.parent_event_ids
    parent = v.make_configuration(cat, fam.id, (f.id,), (rs[0].id, rs[1].id))
    h = CoverageHistory(run_id='neighbors')
    v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, (parent,)), fam.id,
        parent_ids=(parent.id,), history=h)
    assert any(e.changes.get('parent_configuration_id') == parent.id for e in h.events)


@pytest.mark.parametrize('change,reason', [
    ({'accepting_profile_ids': ()}, 'profile'),
    ({'mapping_failure': 'missing'}, 'mapping'),
    ({'reference_footprint': (6, 11)}, 'bounds'),
    ({'alignment_anchor': 7}, 'anchor'),
])
def test_retained_site_validation(change, reason):
    from dataclasses import replace
    v = api(); cat, fam, f, rs = fixture()
    r = replace(rs[0], **change)
    fam = replace(fam, reverse_site_ids=(r.id,))
    cat = replace(cat, sites=(f, r), families=(fam,))
    with pytest.raises(v.ConfigurationIneligible, match=reason):
        v.make_configuration(cat, fam.id, (f.id,), (r.id,))


def test_rejects_stale_reference_footprint_even_when_in_bounds():
    from dataclasses import replace
    v = api(); cat, fam, f, rs = fixture()
    stale = replace(rs[0], reference_footprint=(7, 9))
    fam = replace(fam, reverse_site_ids=(stale.id,))
    cat = replace(cat, sites=(f, stale), families=(fam,))
    with pytest.raises(v.ConfigurationIneligible, match='mapping'):
        v.make_configuration(cat, fam.id, (f.id,), (stale.id,))


def test_parent_neighborhood_materializes_additions_and_swaps():
    v = api(); cat, fam, f, rs = fixture()
    parent = v.make_configuration(cat, fam.id, (f.id,), (rs[0].id,))
    h = CoverageHistory(run_id='edits')
    v.propose_configurations(cat, ConfigurationLedger(cat.semantic_digest, (parent,)), fam.id,
                             parent_ids=(parent.id,), history=h)
    events = [e for e in h.events if e.kind == 'configuration-proposed']
    assert {'addition', 'swap'} <= {e.changes['proposal_reason'] for e in events}
    swap = next(e for e in events if e.changes['proposal_reason'] == 'swap')
    assert swap.changes['removed_site_ids'] == (rs[0].id,)
    assert len(swap.changes['added_site_ids']) == 1
    assert len(swap.changes['lost_supported_allele_ids']) == 1
    assert len(swap.changes['gained_supported_allele_ids']) == 1
