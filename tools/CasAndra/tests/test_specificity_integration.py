"""Meaningful integration guards using frozen old evidence, without fresh truth."""
import copy
import json
from pathlib import Path

import pytest

from casandra.api import Analyzer, DEFAULT_MODEL_BUNDLE
from casandra.specificity import SpecificityOverlay
from casandra._specificity_v3 import validate_auxiliary_batch
from casandra.model import normalize_bundle_config, normalize_protein_card, ModelIntegrityError

FIXTURES = Path(__file__).resolve().parent / 'fixtures'

def rows(p):return [json.loads(s) for s in p.read_text().splitlines() if s]

@pytest.fixture(scope='module')
def model():
    return Analyzer.from_bundle(DEFAULT_MODEL_BUNDLE).protein_model


def load_with(model,configuration):
    return SpecificityOverlay.from_protein_directory(model.directory,configuration,
        model.metadata,model.card,model.hmms)


def test_explicit_registry_and_search_space_identity(model):
    assert model.specificity.manifest['search_spaces']=={'core':361,'inherited':370,'new':387}
    assert len(model.hmms)==361 and len(model.specificity.inherited_hmms)==9 and len(model.specificity.new_hmms)==26
    assert set(model.specificity.card['positive_profile_registry'])=={p for p,m in model.metadata.items() if m['intended_use']=='primary_search'}


def test_tampered_manifest_fails_before_auxiliary_deserialization(model,monkeypatch):
    configuration=copy.deepcopy(model.card['specificity_overlay']);configuration['manifest_sha256']='0'*64
    def forbidden(*args,**kwargs):pytest.fail('Must check digest before loading any auxiliary HMM')
    monkeypatch.setattr('casandra.specificity.pyhmmer.plan7.HMMFile',forbidden)
    with pytest.raises(ValueError,match='checksum'):load_with(model,configuration)


def test_unknown_overlay_schema_fails_closed(model):
    configuration=copy.deepcopy(model.card['specificity_overlay']);configuration['schema_version']=999
    with pytest.raises(ValueError,match='schema'):load_with(model,configuration)
    configuration['schema_version']=True
    with pytest.raises(ValueError,match='schema'):load_with(model,configuration)


def test_required_bundle_format_normalizes_only_recognized_contracts(model):
    path=DEFAULT_MODEL_BUNDLE
    config=json.loads((path/'config.json').read_text());card=json.loads((path/'protein/model_card.json').read_text())
    assert 'protein_model_directory' not in config and 'report_evalue' not in card and 'detection_threshold' not in card
    assert normalize_bundle_config(config)['protein_model_directory']=='protein'
    assert normalize_protein_card(card)['report_evalue']==1000
    bad=copy.deepcopy(config);bad['specificity_runtime_contract']['schema_version']=True
    with pytest.raises(ModelIntegrityError,match='contract'):normalize_bundle_config(bad)
    bad=copy.deepcopy(config);bad['protein_model_directory']='protein'
    with pytest.raises(ModelIntegrityError,match='Ambiguous'):normalize_bundle_config(bad)
    bad=copy.deepcopy(card);bad.pop('specificity_overlay')
    with pytest.raises(ModelIntegrityError,match='incomplete'):normalize_protein_card(bad)
    bad=copy.deepcopy(card);bad['required_specificity_core']['detection_threshold']=True
    with pytest.raises(ModelIntegrityError,match='settings'):normalize_protein_card(bad)


def test_manifest_cannot_escape_verified_protein_directory(model):
    configuration=copy.deepcopy(model.card['specificity_overlay']);configuration['manifest']='../manifest.json'
    with pytest.raises(ValueError,match='escapes'):load_with(model,configuration)


def test_boolean_or_mismatched_auxiliary_z_is_not_accepted(model):
    for bad in [True,361,370,None,'387']:
        with pytest.raises(ValueError,match='space'):validate_auxiliary_batch({'Z':bad,'hits':[]},model.specificity.card)


def test_unknown_auxiliary_profile_is_not_accepted(model):
    with pytest.raises(ValueError,match='Unregistered'):
        validate_auxiliary_batch({'Z':387,'hits':[{'profile':'invented','score':100,'evalue':0}]},model.specificity.card)


def test_canonical_results_and_raw_evidence_are_not_mutated(model,monkeypatch):
    path=FIXTURES/'specificity';queries=rows(path/'input.jsonl')
    original=rows(path/'evidence/production_predictions.jsonl')
    raw=rows(path/'evidence/ordered_hits.jsonl')
    scans=[{'specificity_core_evidence':r} for r in raw]
    old=rows(path/'evidence/previous_auxiliary_hits.jsonl');new=rows(path/'evidence/new_auxiliary_hits.jsonl')
    initial=copy.deepcopy((queries,original,scans,old,new))
    monkeypatch.setattr(SpecificityOverlay,'_scan',lambda self,qs,hmms,Z,threads,**kw:old if Z==370 else new)
    outputs=model.specificity.apply_records(queries,original,scans,threads=2)
    assert (queries,original,scans,old,new)==initial
    # These frozen inputs originated with uncapped V3. Cap15 can retain calls
    # that V3 withheld; verify the shared cumulative contract, not V3 labels.
    assert len(outputs) == len(original)
    assert all(not actual['is_cas'] or baseline['is_cas']
               for actual, baseline in zip(outputs, original))
    assert any(r['abstained_from_baseline_Cas_call'] for r in outputs)
    for r in outputs:
        assert r['result']==(r['cas_family'] if r['is_cas'] else 'no cas')
        assert r['specificity_evidence']['confidence_probability'] is None
        if not r['is_cas']:
            assert r['cas_family'] is None and all(x is None for x in r['classification'].values())
        if not r['hard_negative_hit_present']:
            assert r['raw_best_negative_score'] is None and r['raw_best_negative_profile'] is None
            assert r['hard_negative_profile_score']==-200


def test_margin_bundle_changes_only_enabled_advantages():
    before=json.loads((FIXTURES/'policies/v3.json').read_text())
    after=json.loads((FIXTURES/'policies/margin10.json').read_text())
    after=copy.deepcopy(after);after['id']=before['id'];after.pop('v5_scope')
    changed=0
    for pid,rule in before['profiles'].items():
        for group,comp in rule['competition'].items():
            actual=after['profiles'][pid]['competition'][group]
            expected=min(comp['minimum_score_advantage'],10) if comp['enabled'] else comp.get('minimum_score_advantage')
            assert actual.get('minimum_score_advantage')==expected
            if comp['enabled']:
                changed+=actual['minimum_score_advantage']!=comp['minimum_score_advantage']
                actual['minimum_score_advantage']=comp['minimum_score_advantage']
    assert after==before and changed==274


def test_cap15_floor_maps_match_the_frozen_selection():
    directory = DEFAULT_MODEL_BUNDLE / 'protein/specificity'
    provenance = json.loads((directory / 'floor_cap_provenance.json').read_text())
    policies = {name: json.loads((directory / name).read_text())
                for name in ('inherited_floor.json', 'repair_policy.json')}
    assert provenance['cap'] == 15
    assert provenance['HMM_or_runtime_changes'] is False
    assert provenance['both_floor_maps']
    for item in provenance['both_floor_maps']:
        expected = min(item['before'], 15) if item['before'] is not None else None
        assert item['after'] == expected
        profile = policies[item['file']]['profiles'][item['profile_id']]
        assert profile.get(item['field']) == expected
