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
    """Compare aligned group fits and select a fit/tree with a tri-state decision.

    Log differences and ratios always describe input order versus its reverse;
    the selected protocol's primary metrics describe the recommended fit. ND
    retains the input-order fit. Trivial groups are removed only from the second
    return value, according to their selected orientation's classification.
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

    recommendation = new_protocol['recommend reversing array']
    reverse_selected = recommendation.eq(True)
    # ND retains the forward fit and tree. Select whole rows in one operation
    # so rates, event diagnostics and effective bias provenance stay together.
    oriented_protocol = protocol.where(~reverse_selected, protocol_reversed_raw, axis=0)
    oriented_protocol['recommend reversing array'] = recommendation
    oriented_protocol['predicted orientation'] = [
        'Reverse' if value == True else 'Forward' if value == False else 'ND'  # noqa: E712 - tri-state values
        for value in recommendation
    ]
    for field in ('lh_idm', 'lh_bdm', 'ln_lh_idm', 'ln_lh_bdm'):
        oriented_protocol[f'reversed_{field}'] = protocol_reversed_raw[field].where(
            ~reverse_selected, protocol[field]
        )
    for field in (
        'lh_idm / reversed_lh_idm',
        'ln_lh_idm - reversed_ln_lh_idm',
        'lh_bdm / reversed_lh_bdm',
    ):
        oriented_protocol[field] = new_protocol[field]
    oriented_protocol['confidence threshold'] = decision_boundary
    oriented_protocol['ln_lh_bdm - reversed_ln_lh_bdm'] = new_protocol['ln_lh_bdm - reversed_ln_lh_bdm']

    boring_forward = () if protocol_boring is None else protocol_boring.index
    boring_reverse = () if protocol_boring_reversed is None else protocol_boring_reversed.index
    selected_trivial = (
        (~reverse_selected & protocol.index.isin(boring_forward))
        | (reverse_selected & protocol.index.isin(boring_reverse))
    )
    dict_trees_final = {} if dict_trees is None else {
        name: (dict_tree_reversed if reverse_selected.loc[name] else dict_tree_forward)[name]
        for name in protocol.index
    }
    return new_protocol, oriented_protocol.loc[~selected_trivial].copy(), oriented_protocol.copy(), dict_trees_final


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
