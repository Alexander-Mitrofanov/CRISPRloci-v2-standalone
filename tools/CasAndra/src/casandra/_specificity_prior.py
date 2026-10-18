#!/usr/bin/env python3
"""Experimental read-only acceptance policies; these are not production models."""
import copy
import math


def _finite_number(value, name, *, minimum=None, strictly_positive=False):
    """Reject malformed evidence instead of silently accepting/ranking it."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name + ' must be a finite number')
    if minimum is not None and value < minimum:
        raise ValueError(name + ' is below its allowed minimum')
    if strictly_positive and value <= 0:
        raise ValueError(name + ' must be positive')
    return value


def _validate_numeric_inputs(row, policy):
    for key in ('positive_profile_score', 'hard_negative_profile_score', 'score_margin'):
        if key in row:
            _finite_number(row[key], key)
    if row.get('positive_evalue') is not None:
        _finite_number(row['positive_evalue'], 'positive_evalue', minimum=0)
    kind = policy['kind']
    if kind == 'score':
        _finite_number(policy['floor'], 'score floor')
    elif kind == 'evalue':
        _finite_number(policy['maximum'], 'E-value maximum', strictly_positive=True)
    elif kind == 'profile_floor':
        rule = policy['profiles'].get(row.get('best_positive_profile'))
        if rule and rule.get('floor') is not None:
            _finite_number(rule['floor'], 'profile score floor')
    elif kind != 'baseline':
        raise ValueError('Unknown policy kind: ' + str(kind))


def decision(row, policy):
    _validate_numeric_inputs(row, policy)
    if not row['is_cas']:
        return False, 'baseline_negative'
    kind = policy['kind']
    if kind == 'baseline':
        return True, 'baseline_positive'
    if kind == 'score':
        passed = row['positive_profile_score'] >= policy['floor']
        return passed, 'score_pass' if passed else 'absolute_score_below_floor'
    if kind == 'evalue':
        value = row.get('positive_evalue')
        passed = value is not None and math.isfinite(value) and value <= policy['maximum']
        return passed, 'evalue_pass' if passed else 'positive_evalue_fails'
    if kind == 'profile_floor':
        rule = policy['profiles'].get(row['best_positive_profile'])
        if not rule or rule.get('floor') is None:
            return True, 'uncalibrated_profile_baseline_fallback'
        passed = row['positive_profile_score'] >= rule['floor']
        return passed, 'profile_floor_pass' if passed else 'absolute_score_below_profile_floor'
    raise ValueError('Unknown policy kind: ' + str(kind))


def apply(row, policy):
    result = copy.deepcopy(row)
    accepted, reason = decision(row, policy)
    result['is_cas'] = accepted
    result['policy_evidence'] = {'policy_id': policy['id'], 'baseline_is_cas': row['is_cas'],
        'baseline_classification': copy.deepcopy(row['classification']), 'baseline_cas_family': row.get('cas_family'),
        'decision_reason': reason, 'scope': 'experimental_development_candidate',
        'confidence_probability': None}
    if accepted:
        result['result'] = row.get('cas_family')
    else:
        result['classification'] = {'class': None, 'type': None, 'subtype': None}
        result['cas_family'] = None
        result['result'] = 'no cas'
    return result


def ranking_score(row, policy):
    """Signed acceptance margin for ranking, never a probability/calibration claim."""
    _validate_numeric_inputs(row, policy)
    if row.get('best_positive_profile') is None:
        return -1e6
    base = float(row['score_margin'])
    kind = policy['kind']
    if kind == 'baseline': return base
    if kind == 'score': return min(base, row['positive_profile_score'] - policy['floor'])
    if kind == 'profile_floor':
        rule = policy['profiles'].get(row['best_positive_profile'])
        return min(base, row['positive_profile_score'] - rule['floor']) if rule and rule.get('floor') is not None else base
    if kind == 'evalue':
        e = row.get('positive_evalue')
        return min(base, -math.log10(max(e,1e-300)) + math.log10(policy['maximum'])) if e is not None else -1e6
    raise ValueError(kind)
