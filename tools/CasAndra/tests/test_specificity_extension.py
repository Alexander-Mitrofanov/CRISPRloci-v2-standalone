"""Boundary and cumulative semantics for a bank interface with no fitted rules."""
import copy
import hashlib
import json
from pathlib import Path

import pytest
import pyhmmer

from casandra.api import DEFAULT_MODEL_BUNDLE

from casandra.specificity_extension import AdditionalSpecificityBank, apply_extension, validate_extension_card


def fixture():
    metadata={'P1':{'intended_use':'primary_search'},'P2':{'intended_use':'primary_search'}}
    parent={'mode':'margin10','repair_policy_sha256':'a'*64,
        'inherited_floor_sha256':'b'*64,'inherited_competition_sha256':'c'*64}
    card={'schema_version':1,'id':'synthetic-interface-fixture','auxiliary_Z':362,
        'scan_parameters':{'E':1000.0,'domE':1000.0,'incE':1000.0,'incdomE':1000.0},
        'maximum_competitor_evalue':1e-5,'allowed_advantages':[10,20,40,60],
        'positive_profile_registry':['P1','P2'],'parent':parent,
        'negative_profiles':{'N1':{'group_id':'synthetic','allowed_positive_profiles':['P1']}},
        'rules':{'P1':{'N1':{'enabled':True,'minimum_score_advantage':10,'current_positive_guard_n':1}},'P2':{}}}
    core={'sequence_id':'q','protein_sha256':'d'*64,'is_cas':True,'best_positive_profile':'P1',
        'positive_profile_score':20.0,'hard_negative_profile_score':-200.0,'score_margin':220.0,
        'cas_family':'Cas3','classification':{'class':'1','type':'I','subtype':'I-E'},'result':'Cas3'}
    batch={'sequence_id':'q','protein_sha256':'d'*64,'Z':362.0,'hits':[{'profile':'N1','score':30.0,'evalue':1e-5}]}
    return card,metadata,parent,core,batch


def test_inclusive_thresholds_and_unchanged_inputs():
    card,metadata,parent,core,batch=fixture();validate_extension_card(card,metadata,parent)
    initial=copy.deepcopy((card,core,batch))
    actual=apply_extension(core,core,batch,card)
    assert not actual['is_cas'] and actual['result']=='no cas'
    assert actual['cas_family'] is None and all(v is None for v in actual['classification'].values())
    assert actual['abstained_from_baseline_Cas_call']
    assert actual['positive_profile_score']==20 and actual['hard_negative_profile_score']==-200 and actual['score_margin']==220
    assert (card,core,batch)==initial
    assert actual['additional_specificity_evidence'][0]['parent_call']['classification']==core['classification']
    assert actual['additional_specificity_evidence'][0]['checks'][0]['advantage_minus_threshold']==0


@pytest.mark.parametrize('field,value',[('score',float('nan')),('score',float('inf')),('score',True),('evalue',-1),('evalue',float('nan')),('evalue',True)])
def test_malformed_hit_numbers_rejected(field,value):
    card,_,_,core,batch=fixture();batch['hits'][0][field]=value
    with pytest.raises(ValueError,match='finite'):apply_extension(core,core,batch,card)


@pytest.mark.parametrize('value',[None,True,'362',361,float('inf')])
def test_bad_actual_search_space_rejected(value):
    card,_,_,core,batch=fixture();batch['Z']=value
    with pytest.raises(ValueError,match='search space'):apply_extension(core,core,batch,card)


def test_unknown_duplicate_or_mismatched_query_rejected():
    card,_,_,core,batch=fixture();batch['hits'][0]['profile']='UNKNOWN'
    with pytest.raises(ValueError,match='Unregistered'):apply_extension(core,core,batch,card)
    card,_,_,core,batch=fixture();batch['hits']*=2
    with pytest.raises(ValueError,match='Duplicate'):apply_extension(core,core,batch,card)
    card,_,_,core,batch=fixture();batch['sequence_id']='wrong'
    with pytest.raises(ValueError,match='identity'):apply_extension(core,core,batch,card)


def test_disabled_pair_no_mask_and_upstream_rejection_cannot_veto_or_rescue():
    card,metadata,parent,core,batch=fixture()
    rule=card['rules']['P1']['N1'];rule.update(enabled=False,minimum_score_advantage=None,current_positive_guard_n=0)
    validate_extension_card(card,metadata,parent)
    assert apply_extension(core,core,batch,card)['is_cas']
    card,_,_,core,batch=fixture();core['best_positive_profile']='P2'
    assert apply_extension(core,core,batch,card)['is_cas']
    card,_,_,core,batch=fixture();inherited=copy.deepcopy(core)
    inherited.update(is_cas=False,cas_family=None,result='no cas',classification={'class':None,'type':None,'subtype':None})
    actual=apply_extension(core,inherited,batch,card)
    assert not actual['is_cas'] and actual['abstained_from_baseline_Cas_call']
    assert not actual['additional_specificity_evidence'][0]['applied_extra_veto']
    core=copy.deepcopy(inherited)
    assert not apply_extension(core,inherited,batch,card)['abstained_from_baseline_Cas_call']


def test_guardless_or_unmasked_enabled_rules_rejected():
    card,metadata,parent,_,_=fixture();card['rules']['P1']['N1']['current_positive_guard_n']=0
    with pytest.raises(ValueError,match='witness'):validate_extension_card(card,metadata,parent)
    card,metadata,parent,_,_=fixture();card['negative_profiles']['N1']['allowed_positive_profiles']=['P2']
    with pytest.raises(ValueError,match='mask'):validate_extension_card(card,metadata,parent)
    card,metadata,parent,_,_=fixture();card['rules']['P1']['N1']['minimum_score_advantage']=15
    with pytest.raises(ValueError,match='ladder'):validate_extension_card(card,metadata,parent)


def test_parent_identity_and_registry_required():
    card,metadata,parent,_,_=fixture();other={**parent,'repair_policy_sha256':'0'*64}
    with pytest.raises(ValueError,match='inherited'):validate_extension_card(card,metadata,other)
    card,metadata,parent,_,_=fixture();card['negative_profiles']['N1']['allowed_positive_profiles']=['UNKNOWN']
    with pytest.raises(ValueError,match='masks'):validate_extension_card(card,metadata,parent)


def test_two_stages_preserve_the_first_rejection_and_evidence():
    card,_,_,core,batch=fixture();first=apply_extension(core,core,batch,card)
    card=copy.deepcopy(card);card['id']='second-fixture';batch={**batch,'hits':[]}
    second=apply_extension(core,first,batch,card)
    assert not second['is_cas'] and len(second['additional_specificity_evidence'])==2
    assert second['additional_specificity_evidence'][0]==first['additional_specificity_evidence'][0]
    assert second['additional_specificity_evidence'][1]['reason']=='upstream_rejection_preserved'


def test_serialized_extension_loader_checks_hash_and_registry(tmp_path,monkeypatch):
    card,metadata,parent,_,_=fixture()
    source=DEFAULT_MODEL_BUNDLE/'protein/specificity/inherited_profiles.hmm'
    with pyhmmer.plan7.HMMFile(source) as stream:hmm=stream.read().copy()
    hmm.name='N1'
    with (tmp_path/'named-bank.hmm').open('wb') as stream:hmm.write(stream)
    (tmp_path/'registry.json').write_text(json.dumps({'N1':{'fixture_only':True}}))
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    card.update(bank_hmm_sha256=digest(tmp_path/'named-bank.hmm'),bank_metadata_sha256=digest(tmp_path/'registry.json'))
    (tmp_path/'candidate.json').write_text(json.dumps(card))
    manifest={'schema_version':1,'files':{'policy':'candidate.json','metadata':'registry.json','hmm':'named-bank.hmm'},
        'artifacts':{p.name:{'size':p.stat().st_size,'sha256':digest(p)} for p in tmp_path.iterdir()}}
    path=tmp_path/'bank-manifest.json';path.write_text(json.dumps(manifest))
    pointer={'manifest':path.name,'manifest_sha256':digest(path)}
    loaded=AdditionalSpecificityBank.load(tmp_path,pointer,metadata,parent)
    assert len(loaded.hmms)==1 and loaded.card['auxiliary_Z']==362 and loaded.manifest_sha256==digest(path)
    (tmp_path/'named-bank.hmm').write_bytes(b'altered fixture')
    def forbidden(*args,**kwargs):pytest.fail('Tampered extension must fail before HMM deserialization')
    monkeypatch.setattr('casandra.specificity_extension.pyhmmer.plan7.HMMFile',forbidden)
    with pytest.raises(ValueError,match='checksum'):AdditionalSpecificityBank.load(tmp_path,pointer,metadata,parent)
