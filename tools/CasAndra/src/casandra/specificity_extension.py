"""Versioned interface for a separately scanned, curated competitor bank.

No enabled bank or scientific thresholds are supplied by this module. A bank
must carry its own immutable registry, policy, parent hashes and search space.
The original 361/370/387 scans and their decisions remain unchanged.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import pyhmmer


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def _number(value, name, *, nonnegative=False):
    if type(value) not in (int, float) or not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError('Invalid finite ' + name)
    return float(value)


def _inside(base, relative):
    if not isinstance(relative, str):
        raise ValueError('Extension artifact path must be text')
    result = (base/relative).resolve()
    if not result.is_relative_to(base) or not result.is_file():
        raise ValueError('Extension artifact missing or outside verified directory')
    return result


def validate_extension_card(card, core_metadata, parent):
    """Validate executable masks and prototype rules without inferring truth."""
    if not isinstance(card, dict) or type(card.get('schema_version')) is not int or card['schema_version'] != 1:
        raise ValueError('Unknown additional specificity card schema')
    if not isinstance(card.get('id'), str) or not card['id']:
        raise ValueError('Additional specificity policy ID is missing')
    registry = card.get('negative_profiles')
    if not isinstance(registry, dict) or not registry or not all(isinstance(k, str) and k for k in registry):
        raise ValueError('Additional bank needs a named negative profile registry')
    expected_z = 361 + len(registry)
    if type(card.get('auxiliary_Z')) is not int or card['auxiliary_Z'] != expected_z:
        raise ValueError('Additional bank search space must equal 361 plus its HMM count')
    if card.get('scan_parameters') != {'E':1000.0,'domE':1000.0,'incE':1000.0,'incdomE':1000.0}:
        raise ValueError('Additional bank reporting parameters are not the declared contract')
    if _number(card.get('maximum_competitor_evalue'), 'maximum E-value', nonnegative=True) != 1e-5:
        raise ValueError('Additional bank significance contract mismatch')
    if card.get('allowed_advantages') != [10,20,40,60]:
        raise ValueError('Additional bank advantage ladder mismatch')
    primary = {p for p,m in core_metadata.items() if m['intended_use'] == 'primary_search'}
    positive = card.get('positive_profile_registry')
    if not isinstance(positive, list) or len(positive) != len(set(positive)) or set(positive) != primary:
        raise ValueError('Additional bank positive profile registry mismatch')
    if card.get('parent') != parent:
        raise ValueError('Additional bank was specified for a different inherited policy')
    for negative, entry in registry.items():
        if not isinstance(entry, dict) or not isinstance(entry.get('group_id'), str) or not entry['group_id']:
            raise ValueError('Additional negative profile group is missing')
        mask = entry.get('allowed_positive_profiles')
        if (not isinstance(mask, list) or not mask or not all(isinstance(x, str) for x in mask)
                or len(mask) != len(set(mask)) or not set(mask).issubset(primary)):
            raise ValueError('Additional bank requires explicit valid positive-profile masks')
    rules = card.get('rules')
    if not isinstance(rules, dict) or set(rules) != primary:
        raise ValueError('Additional bank rules require every primary profile, including fallbacks')
    for positive, pairs in rules.items():
        if not isinstance(pairs, dict) or set(pairs)-set(registry):
            raise ValueError('Unknown additional negative profile in a rule')
        for negative, rule in pairs.items():
            if not isinstance(rule, dict) or type(rule.get('enabled')) is not bool:
                raise ValueError('Additional competitor rule requires a Boolean enabled flag')
            count = rule.get('current_positive_guard_n')
            if type(count) is not int or count < 0:
                raise ValueError('Additional competitor rule requires a nonnegative current witness count')
            if rule['enabled']:
                if count == 0:
                    raise ValueError('An enabled additional competitor requires a current positive witness')
                if positive not in registry[negative]['allowed_positive_profiles']:
                    raise ValueError('Enabled additional competitor lies outside its biological mask')
                advantage = _number(rule.get('minimum_score_advantage'), 'score advantage', nonnegative=True)
                if advantage not in card['allowed_advantages']:
                    raise ValueError('Additional competitor advantage is outside the declared ladder')
            elif rule.get('minimum_score_advantage') is not None:
                raise ValueError('Disabled additional competitor must have no executable threshold')


def apply_extension(core, inherited, batch, card):
    """Apply one monotone extra veto, retaining all scores and original heads."""
    z = batch.get('Z')
    if type(z) not in (int, float) or not math.isfinite(z) or z != card['auxiliary_Z']:
        raise ValueError('Actual additional competitor search space mismatch')
    if batch.get('sequence_id') != core.get('sequence_id') or inherited.get('sequence_id') != core.get('sequence_id'):
        raise ValueError('Additional competitor query identity mismatch')
    if 'protein_sha256' in core and batch.get('protein_sha256') != core['protein_sha256']:
        raise ValueError('Additional competitor protein identity mismatch')
    hits = batch.get('hits')
    if not isinstance(hits, list):
        raise ValueError('Additional competitor hits must be a list')
    seen=set()
    for hit in hits:
        if not isinstance(hit, dict) or hit.get('profile') not in card['negative_profiles']:
            raise ValueError('Unregistered additional negative profile')
        if hit['profile'] in seen:
            raise ValueError('Duplicate additional profile hit for one query')
        seen.add(hit['profile'])
        _number(hit.get('score'), 'additional score')
        _number(hit.get('evalue'), 'additional E-value', nonnegative=True)
    if type(core.get('is_cas')) is not bool or type(inherited.get('is_cas')) is not bool:
        raise ValueError('Additional policy requires Boolean upstream calls')
    if inherited['is_cas'] and not core['is_cas']:
        raise ValueError('Additional stage cannot rescue an original non-Cas call')
    positive = core.get('best_positive_profile')
    if positive is not None and positive not in card['rules']:
        raise ValueError('Unknown original positive profile in additional stage')
    score = _number(core.get('positive_profile_score'), 'original positive score')
    checks=[]
    for hit in hits:
        negative = hit['profile']
        pair = card['rules'].get(positive, {}).get(negative)
        masked = positive in card['negative_profiles'][negative]['allowed_positive_profiles']
        enabled = bool(pair and pair['enabled'] and masked)
        advantage = hit['score'] - score
        threshold = pair['minimum_score_advantage'] if enabled else None
        qualifies = enabled and hit['evalue'] <= card['maximum_competitor_evalue'] and advantage >= threshold
        checks.append({'negative_profile':negative,'biological_mask_matches':masked,'pair_enabled':enabled,
            'score_advantage':advantage,'minimum_score_advantage':threshold,
            'advantage_minus_threshold':advantage-threshold if enabled else None,
            'maximum_competitor_evalue':card['maximum_competitor_evalue'],
            'evalue':hit['evalue'],'qualifies_for_veto':qualifies})
    vetoes=[r for r in checks if r['qualifies_for_veto']]
    result=copy.deepcopy(inherited)
    before=copy.deepcopy({k:inherited.get(k) for k in ['is_cas','cas_family','classification','result','decision_status']})
    veto=bool(inherited['is_cas'] and vetoes)
    if veto:
        result['is_cas']=False
        result['cas_family']=None
        result['classification']={'class':None,'type':None,'subtype':None}
        result['decision_status']='abstained_by_additional_specificity'
    result['result']=result['cas_family'] if result['is_cas'] else 'no cas'
    result['abstained_from_baseline_Cas_call']=bool(core['is_cas'] and not result['is_cas'])
    result.setdefault('additional_specificity_evidence', []).append({
        'schema_version':1,'policy_id':card['id'],'parent_call':before,'applied_extra_veto':veto,
        'eligible_veto_profiles':[r['negative_profile'] for r in vetoes],
        'checks':checks,'search':copy.deepcopy(batch),'confidence_probability':None,
        'reason':'additional_competitor_veto' if veto else 'upstream_rejection_preserved' if not inherited['is_cas'] else 'no_enabled_competitor_veto'})
    return result


@dataclass(frozen=True)
class AdditionalSpecificityBank:
    directory: Path
    manifest: dict
    card: dict
    hmms: tuple
    manifest_sha256: str

    @classmethod
    def load(cls, base, pointer, core_metadata, parent):
        base=Path(base).resolve()
        if not isinstance(pointer, dict):
            raise ValueError('Additional bank pointer must be an object')
        path=_inside(base,pointer.get('manifest'))
        if _sha(path) != pointer.get('manifest_sha256'):
            raise ValueError('Additional bank manifest checksum mismatch')
        manifest=json.loads(path.read_text())
        if type(manifest.get('schema_version')) is not int or manifest['schema_version'] != 1:
            raise ValueError('Unknown additional bank manifest schema')
        artifacts=manifest.get('artifacts')
        roles=manifest.get('files',{'policy':'policy.json','hmm':'profiles.hmm','metadata':'metadata.json'})
        if (not isinstance(roles,dict) or set(roles) != {'policy','hmm','metadata'}
                or not all(isinstance(v,str) for v in roles.values())
                or len(set(roles.values())) != 3 or not isinstance(artifacts,dict)
                or not set(roles.values()).issubset(artifacts)):
            raise ValueError('Additional bank manifest is incomplete')
        directory=path.parent
        for relative,pin in artifacts.items():
            artifact=_inside(directory,relative)
            if (not isinstance(pin,dict) or type(pin.get('size')) is not int
                    or artifact.stat().st_size != pin['size'] or _sha(artifact) != pin.get('sha256')):
                raise ValueError('Additional bank artifact checksum mismatch')
        card=json.loads((directory/roles['policy']).read_text())
        validate_extension_card(card,core_metadata,parent)
        if _sha(directory/roles['hmm']) != card.get('bank_hmm_sha256') or _sha(directory/roles['metadata']) != card.get('bank_metadata_sha256'):
            raise ValueError('Additional bank policy/artifact identity mismatch')
        metadata=json.loads((directory/roles['metadata']).read_text())
        if not isinstance(metadata,dict) or set(metadata) != set(card['negative_profiles']):
            raise ValueError('Additional bank metadata registry mismatch')
        with pyhmmer.plan7.HMMFile(directory/roles['hmm']) as handle:
            hmms=tuple(handle)
        ids=[h.name.decode() if isinstance(h.name,bytes) else str(h.name) for h in hmms]
        if len(ids) != len(set(ids)) or set(ids) != set(card['negative_profiles']):
            raise ValueError('Additional HMM registry/count mismatch')
        return cls(directory,manifest,card,hmms,_sha(path))
