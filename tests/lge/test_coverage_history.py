import json
import pytest


def api():
    from primalscheme3.panel import coverage_history
    return coverage_history


def test_contextual_assessments_are_not_global_rejection(tmp_path):
    h = api().CoverageHistory(tmp_path, run_id='run')
    evidence = h.record_evidence(entity_ids=('site',), measurement='dimer',
        dependency_key={'kernel':'v1','sequences':['AA','CC']}, values={'score':-27.}, status='measured')
    failed = h.assess(stage_id='strict', entity_ids=('site',), pool=0, context_digest='pool0',
        profile_id='normal', kernel_versions={'dimer':'v1'}, thresholds={'dimer':-26},
        check_name='dimer', outcome='fail', reason='threshold', evidence_ids=(evidence.id,))
    passed = h.assess(stage_id='strict', entity_ids=('site',), pool=1, context_digest='pool1',
        profile_id='normal', kernel_versions={'dimer':'v1'}, thresholds={'dimer':-26},
        check_name='dimer', outcome='pass', reason='no-conflict')
    unperformed = h.assess(stage_id='strict', entity_ids=('site',), pool=0, context_digest='pool0',
        profile_id='high-gc', kernel_versions={}, thresholds={}, check_name='specificity',
        outcome='not-evaluated', reason='short-circuit')
    assert h.query(entity_id='site', pool=0)['assessments'] == (failed, unperformed)
    assert h.query(entity_id='site', pool=1)['assessments'] == (passed,)
    assert h.query(profile_id='high-gc')['assessments'][0].outcome == 'not-evaluated'
    with pytest.raises(TypeError): evidence.values['score'] = 0
    assert evidence == api().IntrinsicEvidence.from_dict(json.loads(json.dumps(evidence.to_dict())))
    with pytest.raises(ValueError):
        h.assess(stage_id='strict',entity_ids=('site',),pool=0,context_digest='x',profile_id='normal',
                 kernel_versions={},thresholds={},check_name='bad',outcome='skipped',reason='x')


def test_lineage_snapshots_and_reload(tmp_path):
    a = api()
    h = a.CoverageHistory(tmp_path, run_id='run')
    parent = h.emit(stage_id='discovery',kind='generated',entity_ids=('parent',))
    child = h.emit(stage_id='strict',kind='proposed',entity_ids=('child',),parent_event_ids=(parent.id,),
                   changes={'removed_sites':['bad'], 'coverage_delta':{'allele':-.1}})
    snapshot = h.complete_stage(stage_id='strict',dispositions={'parent':'feasible-but-not-selected','child':'selected-strict'},
                                catalog_digest='catalog',ledger_digest='ledger')
    assert snapshot.completeness == 'complete'
    assert snapshot.dispositions['parent'] != 'rejected'
    assert h.query(entity_id='child', include_lineage=True)['events'] == (parent, child)
    again = a.CoverageHistory(tmp_path, run_id='run')
    assert again.events == h.events and again.snapshots == h.snapshots
    assert again.last_complete_stage == snapshot
    third = again.emit(stage_id='salvage-1', kind='reconsidered', entity_ids=('parent',), parent_event_ids=(parent.id,))
    assert third.sequence_number == 2
    failed = again.complete_stage(stage_id='salvage-1', dispositions={}, catalog_digest='catalog',
                                  ledger_digest='ledger', completeness='failed', prior_snapshot_id=snapshot.id)
    assert again.last_complete_stage == snapshot
    assert failed.prior_snapshot_id == snapshot.id
    exported = again.export(tmp_path/'export')
    assert exported == again.export(tmp_path/'export2')


def test_complete_snapshot_cannot_omit_generated_entities(tmp_path):
    h = api().CoverageHistory(tmp_path, run_id='run')
    h.emit(stage_id='discovery',kind='generated',entity_ids=('unexamined',))
    with pytest.raises(ValueError, match='disposition'):
        h.complete_stage(stage_id='strict',dispositions={},catalog_digest='c',ledger_digest='l')
    s = h.complete_stage(stage_id='strict',dispositions={'unexamined':'search-limit-not-explored'},catalog_digest='c',ledger_digest='l')
    assert h.complete_stage(stage_id='salvage',dispositions={},catalog_digest='c',ledger_digest='l',prior_snapshot_id=s.id).completeness == 'complete'


def test_evidence_reuse_uses_complete_dependency_and_corrupt_tail_fails(tmp_path):
    a = api()
    h = a.CoverageHistory(tmp_path, run_id='run')
    kw = dict(entity_ids=('site',), measurement='support', values={'exact':True}, status='measured')
    one = h.record_evidence(dependency_key={'observations':'one'}, **kw)
    two = h.record_evidence(dependency_key={'observations':'two'}, **kw)
    assert one.id != two.id
    assert h.record_evidence(dependency_key={'observations':'one'}, **kw).id == one.id
    assert len(h.evidence) == 2
    with (tmp_path/'events.jsonl').open('a') as handle: handle.write('{broken')
    with pytest.raises(ValueError, match='incomplete|corrupt'):
        a.CoverageHistory(tmp_path, run_id='run')
