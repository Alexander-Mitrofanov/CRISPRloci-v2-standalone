import math

import numpy as np
import pandas as pd


BIAS_PROVENANCE_COLUMNS = (
    'alpha_bias_correction_status',
    'rho_bias_correction_status',
    'rho_bias_failure_policy',
    'rho_bias_correction_warning_code',
    'rho_bias_correction_warning',
)


def compare_likelihoods_for_orientation(protocol, protocol_reversed_raw, protocol_boring=None,
                                        protocol_boring_reversed=None,
                                        decision_boundary=0, dict_trees=None):
    """

    :param dict_trees:
    :param protocol:
    :param protocol_reversed_raw:
    :param protocol_boring:
    :param protocol_boring_reversed:
    :param decision_boundary:
    :return:
    """
    if not math.isfinite(decision_boundary) or decision_boundary < 0:
        raise ValueError("decision_boundary must be finite and non-negative")

    if dict_trees is not None:
        if not isinstance(dict_trees, dict):
            raise ValueError('dict_trees must contain forward and reversed tree maps')
        dict_tree_forward = dict_trees.get('forward')
        dict_tree_reversed = dict_trees.get('reversed')
        if not isinstance(dict_tree_forward, dict) or not isinstance(dict_tree_reversed, dict):
            raise ValueError('dict_trees must contain forward and reversed tree maps')
    else:
        dict_tree_forward = None
        dict_tree_reversed = None

    if protocol_boring is not None:
        if not protocol_boring.empty:
            protocol = pd.concat([protocol, protocol_boring])
        if protocol_boring_reversed is not None and not protocol_boring_reversed.empty:
            protocol_reversed_raw = pd.concat([protocol_reversed_raw, protocol_boring_reversed])

    required_columns = {'lh_idm', 'ln_lh_idm', 'lh_bdm', 'ln_lh_bdm'}
    for label, frame in [('forward', protocol), ('reverse', protocol_reversed_raw)]:
        missing_columns = required_columns.difference(frame.columns)
        if missing_columns:
            raise ValueError(f'{label} protocol is missing required columns: {sorted(missing_columns)}')
        if not frame.index.is_unique:
            raise ValueError(f'{label} protocol contains duplicate group names')
        if not np.isfinite(frame[['ln_lh_idm', 'ln_lh_bdm']].to_numpy(dtype=float)).all():
            raise ValueError(f'{label} protocol contains non-finite log-likelihoods')

    missing_reverse = protocol.index.difference(protocol_reversed_raw.index)
    extra_reverse = protocol_reversed_raw.index.difference(protocol.index)
    if len(missing_reverse) or len(extra_reverse):
        raise ValueError("Forward and reverse protocols must contain the same groups")
    protocol_reversed_raw = protocol_reversed_raw.reindex(protocol.index)
    if dict_trees is not None:
        missing_forward_trees = protocol.index.difference(dict_tree_forward)
        missing_reverse_trees = protocol.index.difference(dict_tree_reversed)
        if len(missing_forward_trees) or len(missing_reverse_trees):
            raise ValueError('Orientation tree maps must contain every reconstructed group')

    reversed_columns = {
        'lh_idm': 'reversed_lh_idm',
        'ln_lh_idm': 'reversed_ln_lh_idm',
        'lh_bdm': 'reversed_lh_bdm',
        'ln_lh_bdm': 'reversed_ln_lh_bdm',
    }
    reversed_columns.update({
        field: f'reversed_{field}'
        for field in BIAS_PROVENANCE_COLUMNS
        if field in protocol_reversed_raw.columns
    })
    protocol_reversed = protocol_reversed_raw.rename(columns=reversed_columns)

    new_protocol = pd.merge(protocol, protocol_reversed[list(reversed_columns.values())],
                            left_index=True, right_index=True, how='left')

    new_protocol['ln_lh_idm - reversed_ln_lh_idm'] = new_protocol['ln_lh_idm'] - new_protocol['reversed_ln_lh_idm']
    new_protocol['ln_lh_bdm - reversed_ln_lh_bdm'] = new_protocol['ln_lh_bdm'] - new_protocol['reversed_ln_lh_bdm']
    with np.errstate(over='ignore', under='ignore'):
        new_protocol['lh_idm / reversed_lh_idm'] = np.exp(
            new_protocol['ln_lh_idm - reversed_ln_lh_idm']
        )
        new_protocol['lh_bdm / reversed_lh_bdm'] = np.exp(
            new_protocol['ln_lh_bdm - reversed_ln_lh_bdm']
        )

    new_protocol['recommend reversing array'] = [decision_fct(v, decision_boundary=decision_boundary)
                                                 for v in new_protocol['ln_lh_bdm - reversed_ln_lh_bdm']]

    oriented_protocol = protocol.copy()
    protocol_boring_index = set() if protocol_boring is None or protocol_boring.empty else set(protocol_boring.index)
    protocol_boring_reversed_index = set() if protocol_boring_reversed is None or protocol_boring_reversed.empty \
        else set(protocol_boring_reversed.index)

    ls_boring = []
    ls_recommendation = []
    ls_pred_orient = []
    ls_decision_bdy = []
    ls_reversed_lh_0, ls_reversed_lh_1, ls_reversed_ln_lh_0, ls_reversed_ln_lh_1 = [], [], [], []
    dict_trees_final = {}
    for i, n in enumerate(protocol.index):
        ls_decision_bdy.append(decision_boundary)
        recommendation = new_protocol.loc[n, 'recommend reversing array']
        if recommendation == True:  # noqa: E712 - tri-state column also contains the string "nd"
            oriented_protocol.loc[n, :] = protocol_reversed_raw.loc[n, :]
            ls_reversed_lh_0.append(new_protocol.loc[n, 'lh_idm'])
            ls_reversed_lh_1.append(new_protocol.loc[n, 'lh_bdm'])
            ls_reversed_ln_lh_0.append(new_protocol.loc[n, 'ln_lh_idm'])
            ls_reversed_ln_lh_1.append(new_protocol.loc[n, 'ln_lh_bdm'])
            ls_pred_orient.append('Reverse')
            ls_recommendation.append(True)
            if n in protocol_boring_reversed_index:
                ls_boring.append(n)

            if dict_trees:
                dict_trees_final[n] = dict_tree_reversed[n]
        else:
            if n in protocol_boring_index:
                ls_boring.append(n)
            ls_reversed_lh_0.append(protocol_reversed_raw.loc[n, 'lh_idm'])
            ls_reversed_lh_1.append(protocol_reversed_raw.loc[n, 'lh_bdm'])
            ls_reversed_ln_lh_0.append(protocol_reversed_raw.loc[n, 'ln_lh_idm'])
            ls_reversed_ln_lh_1.append(protocol_reversed_raw.loc[n, 'ln_lh_bdm'])
            ls_pred_orient.append('Forward' if recommendation == False else 'ND')  # noqa: E712
            ls_recommendation.append(recommendation)

            if dict_trees:
                dict_trees_final[n] = dict_tree_forward[n]
    oriented_protocol['recommend reversing array'] = ls_recommendation
    oriented_protocol['predicted orientation'] = ls_pred_orient
    oriented_protocol['reversed_lh_idm'] = ls_reversed_lh_0
    oriented_protocol['reversed_lh_bdm'] = ls_reversed_lh_1
    oriented_protocol['reversed_ln_lh_idm'] = ls_reversed_ln_lh_0
    oriented_protocol['reversed_ln_lh_bdm'] = ls_reversed_ln_lh_1
    oriented_protocol['lh_idm / reversed_lh_idm'] = new_protocol['lh_idm / reversed_lh_idm']
    oriented_protocol['ln_lh_idm - reversed_ln_lh_idm'] = new_protocol['ln_lh_idm - reversed_ln_lh_idm']
    oriented_protocol['lh_bdm / reversed_lh_bdm'] = new_protocol['lh_bdm / reversed_lh_bdm']
    oriented_protocol['confidence threshold'] = ls_decision_bdy
    oriented_protocol['ln_lh_bdm - reversed_ln_lh_bdm'] = new_protocol['ln_lh_bdm - reversed_ln_lh_bdm']
    oriented_protocol_w_trivial = oriented_protocol.copy()
    oriented_protocol = oriented_protocol.drop(ls_boring, axis=0, errors='ignore')

    # print(oriented_protocol)
    return new_protocol, oriented_protocol, oriented_protocol_w_trivial, dict_trees_final


def decision_fct(ln_lh_diff, decision_boundary=0):
    if not math.isfinite(ln_lh_diff):
        raise ValueError('log-likelihood difference must be finite')
    if not math.isfinite(decision_boundary) or decision_boundary < 0:
        raise ValueError('decision_boundary must be finite and non-negative')
    decision = 'nd'
    if ln_lh_diff < -decision_boundary:
        decision = True
    elif ln_lh_diff > decision_boundary:
        decision = False
    return decision
