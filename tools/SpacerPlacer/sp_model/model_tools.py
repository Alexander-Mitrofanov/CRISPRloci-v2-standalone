import numpy as np
import copy
import scipy.optimize
import pandas as pd
import collections
import math
import os

from sp_model.helpers import misc

EPS = 1e-40
BIG_EPS = 1e20
LOG_EPS = math.log(EPS)


def minimize_fct(lh_fct, start_values, bounds=None, method='L-BFGS-B', **kwargs):
    result = scipy.optimize.minimize(lh_fct, start_values, bounds=bounds, method=method, **kwargs)
    return validate_optimization_result(result, method)


def minimize_scalar_fct(lh_fct, bounds=None, method='bounded', **kwargs):
    result = scipy.optimize.minimize_scalar(lh_fct, bounds=bounds, method=method, **kwargs)
    return validate_optimization_result(result, method)


def validate_optimization_result(result, method):
    """Reject failed or non-finite optimizations instead of emitting plausible-looking output."""
    parameters = np.asarray(result.x)
    if not result.success or not np.isfinite(result.fun) or not np.all(np.isfinite(parameters)):
        raise RuntimeError(
            f'Optimization with {method} failed: {result.message}; '
            f'objective={result.fun!r}, parameters={result.x!r}'
        )
    return result


def compute_lh_ratio_of_multiple_trees(ls_data, method='L-BFGS-B', filter_by=False, lambdifyed_lh_fct=None):
    def lh_fct_0(x):
        return -sum([compute_tree_lh(x, 1.0, d[0], d[1], d[2], log=True) for d in ls_data])

    if lambdifyed_lh_fct is None:
        def lh_fct_1(x):
            return -sum([compute_tree_lh(x[0], x[1], d[0], d[1], d[2], log=True)
                                           for d in ls_data])
    else:
        def lh_fct_1(x):
            return -sum([compute_tree_lh_for_given_lh_fct(x[0], x[1], d[0], d[1], d[2], log=True,
                                                                            lambdifyed_lh_fct=lambdifyed_lh_fct)
                                           for d in ls_data])

    start_values_1 = np.array([0.5, 1.5])
    bnds_0 = (0.0, 10000)
    bnds_1 = ((0.0, None), (1.0, None))

    max_0 = minimize_scalar_fct(lh_fct_0, bounds=bnds_0, method='bounded')
    max_1 = minimize_fct(lh_fct_1, start_values_1, bounds=bnds_1, method=method)

    lh_0 = np.exp(-max_0.fun)
    lh_1 = np.exp(-max_1.fun)

    return lh_0, lh_1, max_0, max_1


def compute_tree_lh(loss_rate, alpha, ls_bl, nb_survivors, nb_max_length_losses, log=True,
                    rec_tree=None, gains=None):
    total_lh = 0.0
    scalar_loss_rate = float(loss_rate)
    if rec_tree is None:
        for bl, n, mll in zip(ls_bl, nb_survivors, nb_max_length_losses):
            loss_prob = max(-math.expm1(-scalar_loss_rate * bl), 0.0)
            total_lh += sum_log_block_likelihoods(loss_prob, alpha, mll) - scalar_loss_rate * bl * n
    else:
        ls_tree_ln_u_lh = rec_tree.compute_unobserved_lh(loss_rate)
        for i, (bl, n, mll, g) in enumerate(zip(ls_bl, nb_survivors, nb_max_length_losses, gains[1:])):
            loss_prob = max(-math.expm1(-scalar_loss_rate * bl), 0.0)
            branch_lh = sum_log_block_likelihoods(loss_prob, alpha, mll) - scalar_loss_rate * bl * n
            ln_u_lh = ls_tree_ln_u_lh[i + 1]
            if scalar_loss_rate != 0 and ln_u_lh is not None and g:
                branch_lh -= log1mexp(ln_u_lh) * g
            total_lh += branch_lh
        if scalar_loss_rate != 0 and gains[0]:
            total_lh -= log1mexp(ls_tree_ln_u_lh[0]) * gains[0]
    return total_lh if log else math.exp(total_lh)


def log1mexp(x):
    """
    Computes log(1 - exp(x)) in a numerically stable way. Note, this is only viable for abs(x) < 1 and x < 0.
    log(2) \approx 0.693 seems to be the optimal cutoff (for highest precision).
    :param x:
    :return:
    """
    # print('log1mexp1', np.log1p(-np.exp(x)))
    # print('log1mexp2', np.log(-np.expm1(x)))
    if x > 0:
        raise ValueError('log1mexp is defined only for values less than or equal to zero')
    if x == 0:
        return -np.inf
    if x < -0.693:
        return math.log1p(-math.exp(x))
    return math.log(-math.expm1(x))


def compute_tree_lh_for_given_lh_fct(loss_rate, alpha, ls_bl, nb_survivors, nb_max_length_losses, log=True,
                                     lambdifyed_lh_fct=None):
    """
    Computes blm-lh of a tree.
    :param lambdifyed_lh_fct:
    :param loss_rate:
    :param alpha:
    :param ls_bl:
    :param nb_survivors:
    :param nb_max_length_losses:
    :param log:
    :return:
    """
    tree_lh = []
    for bl, n, mll in zip(ls_bl, nb_survivors, nb_max_length_losses):
        max_mll = max(mll) if mll else 0
        ls_f = np.log([max(f(bl, loss_rate, alpha), EPS) for f in lambdifyed_lh_fct[:max_mll + 1]])
        ls_f[np.isnan(ls_f) | np.isinf(ls_f)] = -BIG_EPS
        ln_keep = - loss_rate * bl * n
        # Handle blocks that are longer than the available likelihood functions. Divide them in max length blocks
        # + Rest.
        if max_mll + 1 > len(ls_f):
            if len(ls_f) < 2:
                raise ValueError('At least likelihood functions for block lengths 0 and 1 are required')
            ls_ln_mult_loss_lh = []
            for k_i in mll:
                block_size = len(ls_f) - 1
                quotient, remainder = divmod(k_i, block_size)
                val = quotient * ls_f[-1]
                if remainder:
                    val += ls_f[remainder]
                ls_ln_mult_loss_lh.append(val)
        else:
            ls_ln_mult_loss_lh = [ls_f[k_i] for k_i in mll]
        ln_mult_loss_lh = sum(ls_ln_mult_loss_lh)
        tree_lh.append(ln_mult_loss_lh + ln_keep)
    return sum(tree_lh) if log else np.exp(sum(tree_lh))


def eff_mult_loss_lh(loss_prob, alpha, k_i, log=True):
    """Evaluate the geometric block-loss recurrence in linear time.

    The recurrence has the exact closed form
    ``(loss_prob / alpha) * (1 - 1 / alpha + loss_prob / alpha) ** (k - 1)``.
    Computing it directly removes the quadratic convolution formerly performed
    for every likelihood evaluation.
    """
    if k_i < 0:
        raise ValueError('k_i must be non-negative')
    if k_i == 0:
        return np.empty(0, dtype=float) if log else []

    first = loss_prob / alpha
    continuation = 1 - 1 / alpha + first
    if not math.isfinite(first) or not math.isfinite(continuation):
        return np.full(k_i, LOG_EPS, dtype=float) if log else [EPS] * k_i

    if log:
        values = np.full(k_i, LOG_EPS, dtype=float)
        if first > 0:
            log_first = math.log(first)
            values[0] = max(LOG_EPS, log_first)
            if k_i > 1 and continuation > 0:
                exponents = np.arange(1, k_i, dtype=float)
                values[1:] = np.maximum(LOG_EPS, log_first + exponents * math.log(continuation))
        return values
    return [max(EPS, first * continuation ** index) for index in range(k_i)]


def sum_log_block_likelihoods(loss_prob, alpha, block_lengths):
    """Sum geometric block-loss log likelihoods without allocating an array."""
    if not block_lengths:
        return 0.0
    first = loss_prob / alpha
    if not math.isfinite(first) or first <= 0:
        return len(block_lengths) * LOG_EPS
    log_first = math.log(first)
    continuation = 1 - 1 / alpha + first
    if not math.isfinite(continuation):
        return len(block_lengths) * LOG_EPS
    if continuation <= 0:
        return sum(log_first if length == 1 else LOG_EPS for length in block_lengths)
    log_continuation = math.log(continuation)
    longest_block = max(block_lengths)
    if log_first + (longest_block - 1) * log_continuation >= LOG_EPS:
        block_count = len(block_lengths)
        return block_count * log_first + (sum(block_lengths) - block_count) * log_continuation
    return sum(
        max(LOG_EPS, log_first + (length - 1) * log_continuation)
        for length in block_lengths
    )


def mult_loss_lh(loss_prob, alpha, geom_k, log=True):
    """Compute block-loss probabilities iteratively without recursion limits."""
    probabilities = []
    for index, geometric_probability in enumerate(geom_k):
        convolution = sum(
            geom_k[offset] * probabilities[index - offset - 1]
            for offset in range(index)
        )
        probabilities.append(loss_prob * (geometric_probability + convolution))
    return probabilities


def geom(alpha, k):
    if alpha == 1:
        if k == 1:
            return 1
        else:
            return 0
    p = 1 / alpha
    return (1 - p) ** (k - 1) * p


def create_df_alignment(ls_arrays, ls_names, topological_order=None):
    """
    I hate this function, should remove it, because it isn't beneficial except for printing the alignment.
    But I'm also using dfs in other functions so...
    :param ls_arrays:
    :param ls_names:
    :param topological_order:
    :return:
    """
    if len(ls_arrays) != len(ls_names):
        raise ValueError('Spacer arrays and names must have the same length')
    if not ls_arrays:
        raise ValueError('The spacer alignment is empty!')

    if topological_order:
        spacer_sets = [set(array) for array in ls_arrays]
        rows = [
            [spacer if spacer in present else '-' for spacer in topological_order]
            for present in spacer_sets
        ]
        return pd.DataFrame(rows, index=ls_names, columns=topological_order)

    if len(ls_arrays[0]) == 0:
        raise ValueError('The spacer alignment is empty!')
    return pd.DataFrame(
        ls_arrays,
        index=ls_names,
        columns=list(range(1, len(ls_arrays[0]) + 1)),
    )


def remove_unobserved_spacer(df_aligned, tree, gain_dict, loss_dict, gain_loss_event_dict=None, logger=None):
    """
    Could also only search for only '-' columns in alignment, should be cheaper.
    :param tree:
    :param gain_loss_event_dict:
    :param df_aligned:
    :param top_order:
    :param gain_dict:
    :param loss_dict:
    :return:
    """
    top_order = []
    top_order_col = []
    for col in df_aligned.columns:
        unique_val = next((item for item in df_aligned[col] if item != '-'), None)
        if unique_val is not None:
            top_order.append(unique_val)
            top_order_col.append(col)

    # Determines the nodes that are considered to observe spacers, i.e. the leafs in most cases.
    considered_nodes = [node.name for node in tree.get_terminals()]
    observed_spacers = set()
    for name in considered_nodes:
        observed_spacers = observed_spacers.union(set(df_aligned.loc[name]))
    observed_spacers = observed_spacers - {'-'}

    new_top_order = []
    col_to_drop = []
    for s, col in zip(top_order, top_order_col):
        if s not in observed_spacers:
            col_to_drop.append(col)
        else:
            new_top_order.append(s)
    new_df = df_aligned.drop(col_to_drop, axis=1)
    if logger:
        logger.info(f'Original number of spacers: {df_aligned.shape[1]} ; '
                    f'Number of spacers w/o unobserved: {new_df.shape[1]}')
    else:
        print('Original number of spacers: ', df_aligned.shape[1],
              'Number of spacers w/o unobserved: ', new_df.shape[1])
    new_gain_dict = {}
    new_loss_dict = {}
    for key in gain_dict.keys():
        new_gain_dict[key] = [gain for gain in gain_dict[key] if gain in observed_spacers]
        new_loss_dict[key] = [loss for loss in loss_dict[key] if loss in observed_spacers]

    if gain_loss_event_dict:
        for ls_events in gain_loss_event_dict.values():
            retained_events = []
            for event in ls_events:
                event[2] = [s for s in event[2] if s in observed_spacers]
                event[4] = [s for s in event[4] if s in observed_spacers]
                if event[2]:
                    retained_events.append(event)
            ls_events[:] = retained_events

    return new_df, new_top_order, new_gain_dict, new_loss_dict


def find_wrong_order_gains_dict(rec_gains, tree, save_to_tree=True, given_upgraph_order=None):
    """
    Compares previously gained spacers to known younger spacers (upgraph order) of each gain at current node. If there
    is an intersection --> contradiction.
    :param rec_gains:
    :param tree:
    :param save_to_tree:
    :param given_upgraph_order:
    :return:
    """
    contradictions = {}
    gained_spacers = {}
    for c in tree.find_clades(order='preorder'):
        if c.up is None:
            gained_spacers[c.name] = set(rec_gains[c.name])
            continue

        for gain in rec_gains[c.name]:
            if len(gained_spacers[c.up.name].intersection(given_upgraph_order[gain])) > 0:
                contradictions[c.name] = [gain] + contradictions.get(c.name, [])
        gained_spacers[c.name] = gained_spacers[c.up.name].union(set(rec_gains[c.name]))
        if save_to_tree:
            c.contra_gain = contradictions.get(c.name, [])
    return contradictions


def determine_order_dict(leaf_arrays, save_path=None, plot_order=True, exclude_self=True):
    """Return older/younger spacer sets, a stable order, and its adjacency matrix.

    Spacer identities follow first observation. Topological ties retain the
    historical LIFO order because reconstruction numbering depends on it.
    Reachability is iterative so long valid spacer chains need no recursion.
    """
    unique_spacers = list(dict.fromkeys(spacer for array in leaf_arrays for spacer in array))
    positions = {spacer: pos for pos, spacer in enumerate(unique_spacers)}
    adj_matrix = np.zeros((len(positions), len(positions)))
    for array in leaf_arrays:
        for first, second in zip(array, array[1:]):
            adj_matrix[positions[first], positions[second]] = 1
    label_dict = dict(enumerate(unique_spacers))
    roots = np.flatnonzero(~adj_matrix.any(axis=0)).tolist()
    topological_order = misc.topological_sort(roots.copy(), adj_matrix.copy(), label_dict)

    subgraph_order_dict = {pos: set() for pos in label_dict}
    for root in roots:
        get_subgraph_sets(root, adj_matrix, subgraph_order_dict)
    if exclude_self:
        for node, descendants in subgraph_order_dict.items():
            descendants.discard(node)
    named_subgraph_order_dict = {
        label_dict[node]: {label_dict[child] for child in descendants}
        for node, descendants in subgraph_order_dict.items()
    }
    upgraph_order_dict = {spacer: set() for spacer in topological_order}
    for parent, descendants in named_subgraph_order_dict.items():
        for child in descendants:
            upgraph_order_dict[child].add(parent)

    if plot_order:
        from sp_model.spacer_visualization import visualization as vis

        folder, filename = os.path.split(save_path)
        vis.plot_order_w_graphviz(
            adj_matrix, label_dict=label_dict, do_show=False, save_folder=folder,
            graph_name=filename.split('.')[0], file_name=filename + '.dot', color_dict=None,
        )
    return (named_subgraph_order_dict, upgraph_order_dict), topological_order, (adj_matrix, label_dict)


def get_subgraph_sets(n, adj_m, save_dict):
    """Memoize reachable nodes including self without recursive Python calls.

    Preserve populated caller caches and the input matrix. A back edge raises
    the same cycle exception used by the topological order validator.
    """
    pending = [(n, False)]
    active = set()
    while pending:
        node, expanded = pending.pop()
        if save_dict[node]:
            continue
        children = np.flatnonzero(adj_m[node] == 1).tolist()
        if expanded:
            descendants = {node}
            for child in children:
                descendants.update(save_dict[child])
            save_dict[node] = descendants
            active.remove(node)
        else:
            if node in active:
                raise misc.GraphContainsCycleException('Spacer arrays contain (at least one) cycle!')
            active.add(node)
            pending.append((node, True))
            pending.extend((child, False) for child in reversed(children) if not save_dict[child])
    return save_dict[n]


def get_fix_mismatches_df_alignment(df, reverse_walkthrough=False):
    """
    Fixes errors produced by mafft. There are columns with multiple spacers, that need to be separated.
    More problematic are artificial duplications produced by mafft, i.e. columns with the same name that actually can
    be combined, if the PSIO is considered.
    :param df:
    :return:
    """
    if df.empty:
        raise ValueError('The spacer alignment is empty!')
    mismatches = []
    count_mismatches = 0

    ls_resolved_mismatches_col = []
    series_name = 0
    old_uniques = set()
    set_resolved_mismatch_names = set()
    # Separate columns with multiple spacers. They are moved accordingly to what is to the left or right of them, to not
    # introduce additional artificial duplicates.

    # if reverse_walkthrough:
    #     df = df[df.columns[::-1]]

    for idx in range(df.shape[1]):
        series = df.iloc[:, idx]
        ordered_uniques = list(dict.fromkeys(value for value in series if value != '-'))
        uniques = set(ordered_uniques)

        next_uniques = set(df.iloc[:, idx + 1]) - {'-'} if idx + 1 < df.shape[1] else set()
        if len(uniques) > 1:
            mismatches.append(tuple(ordered_uniques))
            count_mismatches += 1
            inter_w_old = uniques.intersection(old_uniques)
            inter_w_next = uniques.intersection(next_uniques)
            inter_all = inter_w_old.intersection(inter_w_next)
            remaining_uniques = uniques - inter_w_old - inter_w_next
            inter_w_old -= inter_all
            inter_w_next -= inter_all
            ls_first = []
            ls_last = []
            ls_middle = []

            for s in (value for value in ordered_uniques if value in remaining_uniques):
                resolved_mismatches_col = series.where(series.eq(s), '-').rename(series_name)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_last.append(resolved_mismatches_col)
            for s in (value for value in ordered_uniques if value in inter_all):
                resolved_mismatches_col = series.where(series.eq(s), '-').rename(series_name)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_middle.append(resolved_mismatches_col)
            for s in (value for value in ordered_uniques if value in inter_w_next):
                resolved_mismatches_col = series.where(series.eq(s), '-').rename(series_name)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_last.append(resolved_mismatches_col)
            for s in (value for value in ordered_uniques if value in inter_w_old):
                resolved_mismatches_col = series.where(series.eq(s), '-').rename(series_name)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_first.append(resolved_mismatches_col)
            ls_resolved_mismatches_col += ls_first + ls_middle + ls_last

        else:
            series.name = series_name
            series_name += 1
            ls_resolved_mismatches_col.append(series)

        old_uniques = set(ls_resolved_mismatches_col[-1]) - {'-'}
    resolved_mismatches_df = pd.concat(ls_resolved_mismatches_col, axis=1)

    # if reverse_walkthrough:
    #     resolved_mismatches_df = resolved_mismatches_df[resolved_mismatches_df.columns[::-1]]

    candidates = [(col, next(iter(set(resolved_mismatches_df[col]) - {'-'})))
                  for col in resolved_mismatches_df.columns if len(set(resolved_mismatches_df[col]) - {'-'}) != 0
                  ]

    duplicates = [item for item, count in collections.Counter([val[1]
                                                               for val in candidates]).items()
                  if count > 1]

    # Check if duplicates can be joined by looking at the order between them.
    dict_of_ls_cols_to_join = dict()
    min_idx = resolved_mismatches_df.shape[1]
    max_idx = 0
    # min_col, max_col = None, None
    for dup in duplicates:
        columns_not_to_join = []
        all_cols = [(i, col) for i, (col, val) in enumerate(candidates) if val == dup]
        new_min_idx, new_min_col = min(all_cols, key=lambda x: x[0])
        new_max_idx, new_max_col = max(all_cols, key=lambda x: x[0])
        if new_min_idx < min_idx:
            min_idx = new_min_idx
            # min_col = list(resolved_mismatches_df.columns)[min_idx]
        if new_max_idx > max_idx:
            max_idx = new_max_idx
            # max_col = list(resolved_mismatches_df.columns)[max_idx]
        in_between_cols = copy.deepcopy(resolved_mismatches_df.iloc[:, min_idx:max_idx + 1])
        ls_renamed_dup = dict()
        columns = set(in_between_cols.columns)
        for d in duplicates:
            # print('d', d)
            dup_counter = 0
            # columns that can be joined w/o issue
            joinables = dict_of_ls_cols_to_join.get(d, [])
            # columns that can be joined w/o issue and are actually in between columns,
            same_name = []
            # new names for same_name spacers
            ls_dup_name = []
            i = 0
            for ls in joinables:
                same_name.append(set([a for a in ls if a in columns]))
                ls_dup_name += [i]
                i += 1
            # Give joinable spacers the same name
            if same_name:
                for j, s in enumerate(same_name):
                    for col in s:
                        for idx in in_between_cols.index:
                            val = in_between_cols.loc[idx, col]
                            in_between_cols.loc[idx, col] = d + '_' + str(ls_dup_name[j]) if val == d else val
            else:
                # Give not joinables individual names
                for col in in_between_cols.columns:
                    # if col not in skip_renaming and d in set(in_between_cols[col]):
                    if d in set(in_between_cols[col]):
                        # print(col)
                        for idx in in_between_cols.index:
                            val = in_between_cols.loc[idx, col]
                            in_between_cols.loc[idx, col] = d + '_' + str(dup_counter) if val == d else val
                        # Remember duplications that need be transformed back.
                        if d == dup:
                            ls_renamed_dup[d + '_' + str(dup_counter)] = col
                        dup_counter += 1
        ls_aligned_arrays = [[value for value in ls if value != '-'] for ls in in_between_cols.values]
        # try:
        (named_subgraph_order_dict, upgraph_order_dict), _, _ = determine_order_dict(ls_aligned_arrays,
                                                                                     save_path=os.path.join(
                                                                                         '../0_result_folder',
                                                                                         'testing'),
                                                                                     plot_order=False,
                                                                                     )

        # Find if up and down order of one individual overlap with any other duplicate candidate,
        # if so they are divided and can't be joined!
        # Should implement this with permutation invariance, i.e. only one of (1, 2) (2, 1).
        for rd, col in ls_renamed_dup.items():
            set_wo_rd = set(ls_renamed_dup.keys()) - {rd}
            subgraph_order_rd = named_subgraph_order_dict[rd]
            upgraph_order_rd = upgraph_order_dict[rd]
            inter_subgraph = set_wo_rd.intersection(subgraph_order_rd)
            inter_upgraph = set_wo_rd.intersection(upgraph_order_rd)
            set_d_not_to_join = inter_upgraph.union(inter_subgraph)
            col_not_to_join = [ls_renamed_dup[s] for s in set_d_not_to_join]

            # col is not allowed to be joined with c.
            columns_not_to_join.extend((col, other_column) for other_column in col_not_to_join)

        ls_usable_cols = list(ls_renamed_dup.values())
        # Put join candidates together for maximum length connected sets, i.e.
        # don't allow not joinable pairs from above
        while ls_usable_cols:
            new_ls = []
            to_keep = []
            for i in range(len(ls_usable_cols)):
                candidate_col = ls_usable_cols[i]
                addable = all([True if (candidate_col, c) not in columns_not_to_join else False for c in new_ls])
                if addable or not new_ls:
                    new_ls.append(candidate_col)
                else:
                    to_keep.append(candidate_col)
            dict_of_ls_cols_to_join[dup] = dict_of_ls_cols_to_join.get(dup, []) + [new_ls]
            ls_usable_cols = to_keep

    top_order = []
    for col in resolved_mismatches_df.columns:
        unique_col = pd.unique(pd.Series([a for a in resolved_mismatches_df[col] if a != '-']))
        if len(unique_col) > 0:
            if len(unique_col) > 1:
                raise Exception('Mismatch in column' + str(col))
            else:
                val = list(unique_col)
        else:
            val = []
        top_order += val

    mask = [True]
    series = pd.Series(top_order)
    dups = []
    count_dups = 0
    while any(mask):
        mask = series.duplicated(keep='last')
        series = series[mask]
        unique_series = pd.unique(series)
        dups += (list(unique_series))
        count_dups += len(unique_series)
    return resolved_mismatches_df, top_order, dict_of_ls_cols_to_join, (count_dups, count_mismatches, mismatches, dups)


def differences_child_parent(child_profile, parent_profile):
    diff_p_ch = child_profile - parent_profile
    # should be careful this relies on maximal length of the profile (for the labeling in the dictionary)
    len_differences = len(diff_p_ch)
    gains = [(len_differences - idx) for idx, val in enumerate(diff_p_ch) if val == 1]
    losses = [(len_differences - idx) for idx, val in enumerate(diff_p_ch) if val == -1]
    gains = list(reversed(gains))
    losses = list(reversed(losses))
    return gains, losses


def differences_child_parent_pos(child_profile, parent_profile):
    diff_p_ch = child_profile - parent_profile
    gains = [idx for idx, val in enumerate(diff_p_ch) if val == 1]
    losses = [idx for idx, val in enumerate(diff_p_ch) if val == -1]
    return gains, losses


def seq2prof(seq, profile_map=None):
    """
    Converts presence absence to likelihood profile (default guaranteed 0 or 1). This goes solely by position.
    1-d -> 2-d array
    :param seq:
    :param profile_map:
    :return:
    """
    if profile_map is None:
        profile_map = {0: [1.0, 0], 1: [0, 1.0]}
    return np.array([profile_map[k] for k in seq])


def prof2seq(profile):
    """
    Converts likelihood profile to presence absence profile. This goes solely by position.
    2-d -> 1-d array
    :param profile:
    :return:
    """
    if len(profile.shape) == 1:
        length = len(profile)
        seq = []
        for i, p in enumerate(profile):
            if p == 1:
                seq.append(length - i)
        seq = np.array(seq)
    else:
        seq = profile.argmax(axis=1).astype(int)
    return seq


def spacer_array_to_binary(spacer_array, alignment=None):
    # if 0 in spacer_array:
    #     raise ValueError('Spacer array contains 0!')
    return np.isin(alignment, test_elements=spacer_array).astype(int)
