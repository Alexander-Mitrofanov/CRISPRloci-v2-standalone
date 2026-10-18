"""Versioned acceptance overlay; original scores and hierarchy are never refitted here."""
import copy
import fnmatch
import math


def number(value, name):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError(name+' must be a finite number')
    return float(value)


def coverage(row):
    h, q = row.get('hmm_coverage'), row.get('query_coverage')
    if h is None or q is None:
        return None
    h, q = number(h, 'HMM coverage'), number(q, 'query coverage')
    if not (0 <= h <= 1 and 0 <= q <= 1):
        raise ValueError('Coverage must be in [0,1]')
    return max(h, q)


def allowed_family(family, patterns):
    return family is not None and any(fnmatch.fnmatchcase(family, pattern) for pattern in patterns)


def validate_card(card):
    """Call once when loading a pinned card, before processing a batch."""
    registry = card['positive_profile_registry']
    if not isinstance(registry, list) or len(set(registry)) != len(registry) or not all(isinstance(x, str) and x for x in registry):
        raise ValueError('Invalid positive-profile registry')
    if type(card['auxiliary_Z']) is not int or card['auxiliary_Z'] != 361+len(card['negative_profiles']):
        raise ValueError('Auxiliary Z must equal361+K')
    if number(card['maximum_competitor_evalue'], 'competitor E cutoff') <= 0:
        raise ValueError('E cutoff must be positive')
    known_groups = set()
    for neg in card['negative_profiles'].values():
        known_groups.add(neg['family_group'])
        patterns = neg['allowed_positive_family_patterns']
        ids = neg['allowed_positive_profile_ids']
        if not isinstance(patterns, list) or not patterns or not all(isinstance(x,str) and x for x in patterns):
            raise ValueError('Family masks must be a nonempty list of strings')
        if not isinstance(ids,list) or not ids or not set(ids).issubset(registry):
            raise ValueError('Invalid allowed winning-profile IDs')
        if neg['family_group'].startswith('RNA_processing_') and (patterns != ['Cas14'] or set(ids) != {'C25_Cas14_1'}):
            raise ValueError('RNA-processing competition must bind to the new type-VII Cas14 profile')
    for pid, rule in card['profiles'].items():
        if pid not in registry:
            raise ValueError('Unknown positive profile in rules')
        if rule.get('score_floor') is not None:
            number(rule['score_floor'], 'score floor')
        if rule.get('coverage_floor') is not None and not 0 <= number(rule['coverage_floor'], 'coverage floor') <= 1:
            raise ValueError('Coverage floor must be in [0,1]')
        for group, comp in rule.get('competition', {}).items():
            if group not in known_groups or type(comp.get('enabled')) is not bool:
                raise ValueError('Invalid competition rule')
            if comp['enabled'] and number(comp['minimum_score_advantage'], 'score advantage') < 0:
                raise ValueError('Competition advantage must be nonnegative')


def validate_auxiliary_batch(batch, card):
    if type(batch.get('Z')) not in (int,float) or batch['Z'] != card['auxiliary_Z']:
        raise ValueError('Auxiliary evidence search space mismatch')
    for hit in batch['hits']:
        if hit['profile'] not in card['negative_profiles']:
            raise ValueError('Unregistered auxiliary profile')
        number(hit['score'], 'negative score')
        if number(hit['evalue'], 'negative E-value') < 0:
            raise ValueError('Negative E-value cannot be negative')
    return batch['hits']


def evaluate(row, auxiliary_batch, card, arm):
    if arm not in {'support_gate', 'competitor_gate', 'intersection'}:
        raise ValueError('Unknown v3 arm')
    if type(row['is_cas']) is not bool:
        raise ValueError('is_cas must be Boolean')
    if not row['is_cas']:
        return {'accepted': False, 'status': 'baseline_negative', 'reasons': ['baseline_negative'],
                'coverage_evaluable': coverage(row) is not None, 'competing_hits': []}
    score = number(row['positive_profile_score'], 'positive score')
    auxiliary_hits = validate_auxiliary_batch(auxiliary_batch, card)
    pid = row['best_positive_profile']
    cov = coverage(row)
    rule = card['profiles'].get(pid, {})
    reasons = []
    notes = []
    matches = []
    if arm in {'support_gate', 'intersection'}:
        floor = rule.get('score_floor')
        cov_floor = rule.get('coverage_floor')
        if floor is not None and score < number(floor, 'score floor'):
            reasons.append('positive_score_below_profile_floor')
        if cov_floor is not None:
            if cov is None:
                notes.append('coverage_gate_unevaluable_missing_domain_evidence')
            elif cov < number(cov_floor, 'coverage floor'):
                reasons.append('alignment_support_below_profile_floor')
        if floor is None and cov_floor is None:
            notes.append('support_gate_explicit_baseline_fallback')
    if arm in {'competitor_gate', 'intersection'}:
        for hit in auxiliary_hits:
            profile = hit['profile']
            neg = card['negative_profiles'][profile]
            group = neg['family_group']
            if pid not in neg['allowed_positive_profile_ids']:
                continue
            if not allowed_family(row.get('cas_family'), neg['allowed_positive_family_patterns']):
                continue
            competition = rule.get('competition', {}).get(group)
            if not competition or not competition.get('enabled'):
                continue
            evalue = number(hit['evalue'], 'negative E-value')
            neg_score = number(hit['score'], 'negative score')
            if evalue < 0:
                raise ValueError('Negative E-value cannot be negative')
            advantage = number(competition['minimum_score_advantage'], 'competition advantage')
            cutoff = number(card['maximum_competitor_evalue'], 'negative E-value cutoff')
            if evalue <= cutoff and neg_score >= score+advantage:
                matches.append({'profile': profile, 'family_group': group, 'score': neg_score,
                                'evalue': evalue, 'score_advantage': neg_score-score,
                                'required_advantage': advantage,
                                'scope': competition['scope']})
        if matches:
            reasons.append('stronger_supported_alternative_function_profile')
    accepted = not reasons
    status = 'accepted' if accepted else ('uncertain_Cas_association' if matches else 'insufficient_profile_support')
    return {'accepted': accepted, 'status': status, 'reasons': reasons or ['gates_passed'],
            'notes': notes, 'coverage_evaluable': cov is not None,
            'max_hmm_or_query_coverage': cov, 'competing_hits': matches}


def apply(row, auxiliary_batch, card, arm):
    result = copy.deepcopy(row)
    decision = evaluate(row, auxiliary_batch, card, arm)
    result['is_cas'] = decision['accepted']
    result['decision_status'] = decision['status']
    result['abstained_from_baseline_Cas_call'] = row['is_cas'] and not decision['accepted']
    result['repair_evidence'] = {'candidate_id': card['id'], 'arm': arm,
        'baseline_is_cas': row['is_cas'], 'baseline_cas_family': row.get('cas_family'),
        'baseline_classification': copy.deepcopy(row['classification']),
        'baseline_result': row.get('result'), 'positive_profile': row['best_positive_profile'],
        'hmm_coverage': row.get('hmm_coverage'), 'query_coverage': row.get('query_coverage'),
        'decision': decision, 'auxiliary_Z': card['auxiliary_Z'], 'confidence_probability': None,
        'scope': 'local_research_candidate; abstentions count as missed positive recovery'}
    if not decision['accepted']:
        result['classification'] = {'class': None, 'type': None, 'subtype': None}
        result['cas_family'] = None
        result['result'] = 'no cas' if not row['is_cas'] else 'uncertain Cas assignment'
    return result


def apply_cumulative(core, previous_combined, auxiliary_batch, card, arm):
    """Both new-stage evidence and total recovery remain relative to the original core."""
    if core['sequence_id'] != previous_combined['sequence_id']:
        raise ValueError('Core and inherited prediction IDs differ')
    if previous_combined['is_cas'] and not core['is_cas']:
        raise ValueError('Inherited policy must not rescue a core noncall')
    result = apply(core, auxiliary_batch, card, arm)
    new_stage = result['is_cas']
    result['is_cas'] = previous_combined['is_cas'] and new_stage
    result['repair_evidence']['new_stage_accepted'] = new_stage
    result['repair_evidence']['inherited_combined_accepted'] = previous_combined['is_cas']
    result['repair_evidence']['inherited_policy_evidence'] = copy.deepcopy(previous_combined.get('policy_evidence'))
    result['repair_evidence']['inherited_reference_evidence'] = copy.deepcopy(previous_combined.get('reference_competition_evidence'))
    result['abstained_from_baseline_Cas_call'] = core['is_cas'] and not result['is_cas']
    if core['is_cas'] and not previous_combined['is_cas']:
        result['decision_status'] = 'inherited_combined_rejection'
    if not result['is_cas']:
        result['classification'] = {'class': None, 'type': None, 'subtype': None}
        result['cas_family'] = None
        result['result'] = 'uncertain Cas assignment' if core['is_cas'] else 'no cas'
    return result
