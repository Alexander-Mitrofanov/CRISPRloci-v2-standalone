#!/usr/bin/env python3
"""Auxiliary-reference veto, isolated from original CasAndra implementation."""
import copy
import math
from casandra._specificity_prior import apply, _finite_number

def apply_competition(row,hits,card):
    # The auxiliary bank was scanned at this fixed search space. Do not label
    # evidence from a mismatched or malformed card as if that scan were used.
    if type(card.get('Z')) is not int or card['Z'] != 370:
        raise ValueError('Reference competition requires the pinned integer Z=370')
    result=copy.deepcopy(row)
    references=card['allowed_references']
    if not isinstance(references,list) or any(not isinstance(x,str) or not x.strip() for x in references) or len(set(references))!=len(references):
        raise ValueError('allowed_references must be a list of unique nonempty identifiers')
    allowed=set(references)
    margin=_finite_number(card['minimum_score_advantage'],'score advantage',minimum=0)
    maximum=_finite_number(card['maximum_evalue'],'E-value maximum',strictly_positive=True)
    score=_finite_number(row['positive_profile_score'],'positive score')
    eligible=[]
    for h in hits:
        _finite_number(h['score'],'reference score');_finite_number(h['evalue'],'reference E-value',minimum=0)
        if h['reference_id'] in allowed and h['evalue']<=maximum:
            eligible.append(h)
    best=max(eligible,key=lambda h:h['score']) if eligible else None
    veto=bool(row['is_cas'] and best is not None and best['score']>=score+margin)
    result['reference_competition_evidence']={'policy_id':card['id'],'baseline_is_cas':row['is_cas'],
        'baseline_classification':copy.deepcopy(row['classification']),'baseline_cas_family':row.get('cas_family'),
        'veto':veto,'best_eligible_reference':copy.deepcopy(best),'score_advantage':best['score']-score if best else None,
        'Z':370,'minimum_score_advantage':margin,'maximum_evalue':maximum,
        'confidence_probability':None,'scope':card['scope']}
    if veto:
        result['is_cas']=False;result['classification']={'class':None,'type':None,'subtype':None}
        result['cas_family']=None;result['result']='no cas'
    else:result['result']=row.get('cas_family') if row['is_cas'] else 'no cas'
    return result
