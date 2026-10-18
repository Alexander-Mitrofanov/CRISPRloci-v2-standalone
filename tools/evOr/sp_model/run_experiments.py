import logging
from concurrent.futures import ThreadPoolExecutor
import os
import pandas as pd
import copy
import pickle
import numpy as np
import json
from functools import lru_cache
import warnings

from sp_model import mafft_integration
from sp_model.helpers import misc, stats
from sp_model.orientation import orientation_tools
from sp_model.experiments_tools import load_align_pickled_data, load_arrays_from_pkl, align_crispr_groups_for_sim_as_rec, pooling_for_parameter_estimation, \
    load_align_single_fasta, integer_group_labels, load_aligned_groups_json, write_aligned_groups_json, \
    validate_aligned_groups_inputs
from sp_model.run_reconstruction import run_reconstruction

WORK_PATH = os.path.join('data', 'simulation_alignment')
LH_FCT_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'sp_model',
                           'bdm_likelihood_computation', 'sympy_bdm_lh_fct',
                           '230329_death_lh_up_to_68_lambdifyed.pickle')
LH_RAW_FCT_PATH = os.path.join(os.path.dirname(LH_FCT_PATH), '230329_death_lh_up_to_68.pickle')
# LH_FCT_PATH = os.path.join('../additional_data', '0_lh',
#                            '230329_death_lh_up_to_68_lambdifyed.pickle'
#                            )


@lru_cache(maxsize=1)
def load_ode_likelihood_functions():
    """Load compatible cached callables or regenerate them from portable SymPy data."""
    from sp_model.bdm_likelihood_computation import symbolic_lh_computation

    try:
        likelihood_functions = symbolic_lh_computation.load_pickled_lambdified_lh_fct(LH_FCT_PATH)
        validation_value = likelihood_functions[1](0.1, 0.2, 1.5)
        if not np.isfinite(float(np.asarray(validation_value).squeeze())):
            raise ValueError('cached likelihood validation returned a non-finite value')
        return likelihood_functions
    except (IndexError, TypeError, ValueError):
        # Dill serializes Python code objects, which are not portable across all
        # interpreter versions. The raw SymPy matrix is portable but slower to
        # lambdify; cache the regenerated callables for the rest of the process.
        warnings.warn(
            'Bundled ODE callable cache is incompatible with this Python version; '
            'regenerating it from the portable SymPy matrix for this process.',
            RuntimeWarning,
            stacklevel=2,
        )
        return symbolic_lh_computation.load_lambdify_lh_fct(
            LH_RAW_FCT_PATH, save_lambdified_lh_fct=False
        )


def _resolve_likelihood_function(likelihood):
    """Resolve CLI ODE aliases while retaining simple/custom likelihood inputs."""
    if isinstance(likelihood, str) and likelihood.lower() in {'ode_based', 'precise'}:
        return load_ode_likelihood_functions()
    return likelihood


def _write_orientation_outputs(save_path, comparison, oriented, forward_trees, reverse_trees, selected_trees):
    """Write the common CLI contract; legacy workflows may also retain their filenames."""
    detailed = os.path.join(save_path, 'detailed_results')
    additional = os.path.join(save_path, 'additional_data')
    os.makedirs(detailed, exist_ok=True)
    os.makedirs(additional, exist_ok=True)
    comparison.to_csv(os.path.join(detailed, '0_detailed_results_w_orientation.csv'))
    comparison.to_pickle(os.path.join(detailed, '0_detailed_results_w_orientation.pkl'))
    oriented.to_csv(os.path.join(detailed, '0_detailed_results_wo_trivial_groups.csv'))
    oriented.to_pickle(os.path.join(detailed, '0_detailed_results_wo_trivial_groups.pkl'))
    tree_maps = {'forward': forward_trees, 'reverse': reverse_trees, 'selected': selected_trees}
    with open(os.path.join(additional, 'orientation_trees.json'), 'w', encoding='utf-8') as tree_file:
        json.dump(tree_maps, tree_file, indent=2)


def _require_reconstructed_protocols(protocol, trivial_protocol, forward_save_path):
    if protocol.empty and trivial_protocol.empty:
        skipped_path = os.path.join(forward_save_path, '0_protocol_skipped.csv')
        raise ValueError(
            'No groups could be reconstructed. Check input groups and tree names; '
            f'see the skipped-group report at {skipped_path}'
        )


def run_multiple_groups(ls_data_path, save_path, rec_parameter_dict=None, lh_fct=None, logger=None, plot_tree=True,
                        do_show=False, combine_non_unique_arrays=False, determine_orientation=True,
                        orientation_decision_boundary=5,
                        tree_path=None, plot_order=True, significance_level=0.05, extend_branches=False,
                        tree_distance_fct='likelihood', tree_construction_method='nj',
                        tree_lh_fct=None, tree_insertion_rate=None, tree_deletion_rate=None, tree_alpha=None,
                        alpha_bias_correction=True, rho_bias_correction=True,
                        rho_bias_failure_policy='fail',
                        seed=None,
                        save_reconstructed_events=False,
                        dpi=90,
                        figsize_rec=(None, None, 'px'),
                        show_spacer_name_row=True,
                        alignment_workers=1,
                        aligned_groups_path=None,
                        ):
    if isinstance(alignment_workers, bool) or not isinstance(alignment_workers, (int, np.integer)) or alignment_workers < 1:
        raise ValueError('alignment_workers must be a positive integer')
    paths_by_group = {}
    for data_path in ls_data_path:
        group_name = os.path.splitext(os.path.basename(data_path))[0]
        if group_name in paths_by_group:
            raise ValueError(f'Duplicate input group name: {group_name}')
        paths_by_group[group_name] = data_path
    if not paths_by_group:
        raise ValueError('At least one input group is required')
    if seed is not None:
        np.random.seed(seed)
    if logger is None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        logger = misc.create_logger('multiple groups reconstruction', logging.INFO, outfile=os.path.join(save_path,
                                                                                                      '0_logger.log'))
    if rec_parameter_dict is None:
        rec_parameter_dict = {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 1e-1}
    dict_crispr_groups = {}
    save_data_folder = os.path.join(save_path, 'additional_data')
    logger.info('Loading and aligning the data...')
    def align_group(index_and_path):
        i, data_path = index_and_path
        group_name = os.path.splitext(os.path.basename(data_path))[0]
        ind_save_path = os.path.join(save_data_folder, 'work_folder', group_name)
        os.makedirs(ind_save_path, exist_ok=True)
        logger.info(f'{i + 1} / {len(ls_data_path)} aligning group {group_name}... \n')
        crispr_group = load_align_single_fasta(data_path, ind_save_path,
                                               group_name=group_name,
                                               logger=logger,
                                               mafft_options=None,
                                               save_path=None,
                                               seed=seed)
        return group_name, crispr_group

    if aligned_groups_path is not None:
        logger.info('Reusing validated alignments from %s', aligned_groups_path)
        cached_groups = load_aligned_groups_json(aligned_groups_path)
        validate_aligned_groups_inputs(cached_groups, paths_by_group)
        dict_crispr_groups.update(cached_groups)
    else:
        indexed_paths = list(enumerate(ls_data_path))
        if alignment_workers == 1 or len(indexed_paths) < 2:
            aligned_groups = [align_group(item) for item in indexed_paths]
        else:
            with ThreadPoolExecutor(max_workers=min(alignment_workers, len(indexed_paths))) as pool:
                aligned_groups = list(pool.map(align_group, indexed_paths))
        dict_crispr_groups.update(aligned_groups)
    print_save_path = os.path.join(save_data_folder, 'work_folder')
    logger.debug(f'Saving CRISPR groups to: {print_save_path}')
    if not os.path.exists(save_data_folder):
        os.makedirs(save_data_folder)
    write_aligned_groups_json(dict_crispr_groups, os.path.join(save_data_folder, 'aligned_groups.json'))
    with open(os.path.join(save_data_folder, 'dict_crispr_groups.pkl'), 'wb') as group_file:
        pickle.dump(dict_crispr_groups, group_file)

    lh_fct = _resolve_likelihood_function(lh_fct)
    tree_lh_fct = _resolve_likelihood_function(tree_lh_fct)
    dict_crispr_groups_for_reverse = copy.deepcopy(dict_crispr_groups) if determine_orientation else {}
    logger.info('Beginning reconstruction of (forward oriented) group(s)...\n')
    df_rec_protocol, \
        df_rec_protocol_trivial, _, \
        dict_provided_numbering, dict_trees_forward, \
        dict_provided_aligned_arrays, \
        dict_provided_duplicated_spacers = run_reconstruction(rec_parameter_dict,
                                                              dict_crispr_groups,
                                                              save_path=os.path.join(
                                                                  save_path,
                                                                  '0_forward'),
                                                              plot_tree=plot_tree,
                                                              lh_fct=lh_fct,
                                                              logfile_path=None,
                                                              logger=logger, do_show=do_show,
                                                              combine_non_unique_arrays=combine_non_unique_arrays,
                                                              tree_path=tree_path,
                                                              use_provided_trees=True if tree_path is not None
                                                              else False,
                                                              tree_save_path=None,
                                                              hide_unobserved_spacers=False,
                                                              plot_order=plot_order,
                                                              finite_bl_at_root=True,
                                                              group_by=None,
                                                              significance_level=significance_level,
                                                              extend_branches=extend_branches,
                                                              tree_lh_fct=tree_lh_fct,
                                                              tree_construction_method=tree_construction_method,
                                                              tree_distance_function=tree_distance_fct,
                                                              tree_gain_rate=tree_insertion_rate,
                                                              tree_loss_rate=tree_deletion_rate,
                                                              tree_alpha=tree_alpha,
                                                              alpha_bias_correction=alpha_bias_correction,
                                                              rho_bias_correction=rho_bias_correction,
                                                              rho_bias_failure_policy=rho_bias_failure_policy,
                                                              core_genome_trees=False,
                                                              metadata=False,
                                                              save_reconstructed_events=save_reconstructed_events,
                                                              dpi=dpi,
                                                              figsize_rec=figsize_rec,
                                                              show_spacer_name_row=show_spacer_name_row,
                                                              )
    _require_reconstructed_protocols(df_rec_protocol, df_rec_protocol_trivial,
                                     os.path.join(save_path, '0_forward'))
    dict_trees = dict_trees_forward
    if not df_rec_protocol.empty:
        df_rec_protocol = df_rec_protocol.set_index('name')
    if not df_rec_protocol_trivial.empty:
        df_rec_protocol_trivial = df_rec_protocol_trivial.set_index('name')

    df_results_wo_details = pd.concat([df_rec_protocol, df_rec_protocol_trivial], axis=0)

    for gn, provided_numbering in dict_provided_numbering.items():
        dict_spacer_names_to_numbers = provided_numbering[0]
        p = os.path.join(save_data_folder, 'spacer_fasta', gn + '_spacer_name_to_sp_number.json')
        if not os.path.exists(os.path.dirname(p)):
            os.makedirs(os.path.dirname(p))
        with open(p, 'w') as f:
            json.dump(dict_spacer_names_to_numbers, f)
    if determine_orientation:
        logger.info('Beginning reconstruction of reversed group(s)...\n')
        if not os.path.exists(os.path.join(save_path, '0_reversed')):
            os.makedirs(os.path.join(save_path, '0_reversed'))

        for gn, crispr_group in dict_crispr_groups_for_reverse.items():
            for ca in crispr_group.crispr_dict.values():
                ca.spacer_array = ca.spacer_array[::-1]
        dict_provided_aligned_arrays = {k: v[v.columns[::-1]] for k, v in dict_provided_aligned_arrays.items()}

        df_rec_protocol_reversed, \
            df_rec_protocol_reversed_trivial, \
            _, _, dict_trees_reversed, \
            _, _ = run_reconstruction(rec_parameter_dict,
                                      dict_crispr_groups_for_reverse,
                                      save_path=os.path.join(
                                          save_path,
                                          '0_reversed'),
                                      plot_tree=plot_tree,
                                      lh_fct=lh_fct,
                                      logfile_path=None,
                                      logger=logger, do_show=do_show,
                                      combine_non_unique_arrays=combine_non_unique_arrays,
                                      tree_path=tree_path,
                                      use_provided_trees=True if tree_path is not None
                                      else False,
                                      tree_save_path=None,
                                      hide_unobserved_spacers=False,
                                      plot_order=plot_order,
                                      finite_bl_at_root=True,
                                      group_by=None,
                                      significance_level=significance_level,
                                      extend_branches=extend_branches,
                                      tree_lh_fct=tree_lh_fct,
                                      tree_construction_method=tree_construction_method,
                                      tree_distance_function=tree_distance_fct,
                                      tree_gain_rate=tree_insertion_rate,
                                      tree_loss_rate=tree_deletion_rate,
                                      tree_alpha=tree_alpha,
                                      alpha_bias_correction=alpha_bias_correction,
                                      rho_bias_correction=rho_bias_correction,
                                      rho_bias_failure_policy=rho_bias_failure_policy,
                                      core_genome_trees=False,
                                      provided_numbering=dict_provided_numbering,
                                      provided_aligned_arrays=dict_provided_aligned_arrays,
                                      provided_dict_duplicated_spacers=dict_provided_duplicated_spacers,
                                      metadata=False,
                                      save_reconstructed_events=save_reconstructed_events,
                                      dpi=dpi,
                                      figsize_rec=figsize_rec,
                                      show_spacer_name_row=show_spacer_name_row,)
        dict_trees = {'forward': dict_trees_forward,
                      'reversed': dict_trees_reversed}
        if not df_rec_protocol_reversed.empty:
            df_rec_protocol_reversed = df_rec_protocol_reversed.set_index('name')
        if not df_rec_protocol_reversed_trivial.empty:
            df_rec_protocol_reversed_trivial = df_rec_protocol_reversed_trivial.set_index('name')

        (df_rec_protocol, df_oriented_rec_protocol,
         df_oriented_rec_protocol_w_trivial, dict_trees) = orientation_tools.compare_likelihoods_for_orientation(
            df_rec_protocol,
            df_rec_protocol_reversed,
            df_rec_protocol_trivial,
            df_rec_protocol_reversed_trivial,
            decision_boundary=orientation_decision_boundary, dict_trees=dict_trees)

        _write_orientation_outputs(
            save_path, df_rec_protocol, df_oriented_rec_protocol,
            dict_trees_forward, dict_trees_reversed, dict_trees,
        )
        df_results_wo_details = df_oriented_rec_protocol_w_trivial

    df_results_wo_details = df_results_wo_details.drop(
        columns=['relative deletion positions', 'nb of existent spacers',
                 'all max block deletion lengths',
                 'all max block deletion lengths (normalized)',
                 'reversed_lh_idm', 'reversed_lh_bdm',
                 'lh_idm / reversed_lh_idm',
                 'lh_bdm / reversed_lh_bdm',
                 'lh_idm', 'lh_bdm',
                 'excluding young spacers deletion_rate_idm',
                 'excluding young spacers deletion_rate_bdm',
                 'excluding young spacers alpha_bdm',
                 'containing old spacers deletion_rate_idm',
                 'containing old spacers deletion_rate_bdm',
                 'containing old spacers alpha_bdm', 'fes +m presence',
                 'les -m presence',
                 'ls -m presence', 'fs +m presence', 'global fes +m presence',
                 'global les -m presence', 'global ls -m presence',
                 # 'combined array names',
                 ],
        errors='ignore')

    df_results_wo_details.to_csv(os.path.join(save_path, '_'.join(['0_results']) + '.csv'))

    # write_summary()
    # summary = {**summary_forward, **summary_reversed}

    tree_save_path = os.path.join(save_path, 'additional_data')
    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'dict_nwk_trees.json'), 'w') as f:
            json.dump(dict_trees, f)

    return df_results_wo_details


def run_pickled_data(rec_parameter_dict, data_path, save_path=None, plot_tree=True,
                     logfile_path=None,
                     do_show=False,
                     combine_non_unique_arrays=True,
                     tree_path=None,
                     tree_save_path=None,
                     hide_unobserved_spacers=False,
                     selection_fct=None,
                     plot_order=True,
                     lh_fct=None,
                     significance_level=0.05,
                     finite_bl_at_root=True,
                     group_by=None,
                     alignment_options_dict=None,
                     extend_branches=False,
                     tree_distance_function='breakpoint',
                     tree_gain_rate=None,
                     tree_loss_rate=None,
                     tree_alpha=None,
                     determine_orientation=False,
                     orient_boundary=0,
                     alpha_bias_correction=False,
                     rho_bias_correction=False,
                     rho_bias_failure_policy='fail',
                     core_genome_trees=False,
                     save_reconstructed_events=False,
                     dpi=90,
                     figsize_rec=(None, None, 'px'),
                     tree_lh_fct=None,
                     tree_construction_method='nj',
                     show_spacer_name_row=True,
                     ):
    if isinstance(logfile_path, (str, os.PathLike)):
        logger = misc.create_logger('pickled data based reconstruction', logging.INFO, outfile=logfile_path)
    elif logfile_path is None:
        logger = misc.create_logger('pickled data based reconstruction', logging.INFO)
    elif isinstance(logfile_path, logging.Logger):
        logger = logfile_path
    else:
        raise TypeError('logfile_path must be a path, logging.Logger, or None')

    if data_path.split('.')[-1] in ['pkl', 'pickle']:
        logger.info(f'Loading data from pickle {data_path}')
        dict_crispr_groups = load_arrays_from_pkl(data_path)
        logger.info(f'Number of groups before reconstruction: {len(dict_crispr_groups)}')
        nb_arrays = sum([len(g.crispr_dict) for g in dict_crispr_groups.values()])
        logger.info(f'Number of arrays before reconstruction: {nb_arrays}')
    else:
        logger.info(f'Loading and aligning from files {data_path}')
        dict_crispr_groups = load_align_pickled_data(data_path, mafft_options=alignment_options_dict['mafft_options'],
                                                     exclude_files=alignment_options_dict['exclude_files'],
                                                     save_path=alignment_options_dict['save_path'])

    lh_fct = _resolve_likelihood_function(lh_fct)
    tree_lh_fct = _resolve_likelihood_function(tree_lh_fct)
    dict_crispr_groups_for_reverse = copy.deepcopy(dict_crispr_groups) if determine_orientation else {}
    (df_rec_protocol, df_rec_protocol_boring, (dict_data_for_group_by, lh_fct), dict_provided_numbering,
     dict_trees_forward,
     dict_provided_aligned_arrays,
     dict_provided_duplicated_spacers) = run_reconstruction(rec_parameter_dict,
                                                            dict_crispr_groups,
                                                            save_path=os.path.join(
                                                                save_path,
                                                                '0_forward'),
                                                            plot_tree=plot_tree,
                                                            lh_fct=lh_fct,
                                                            logfile_path=logfile_path,
                                                            logger=logger, do_show=do_show,
                                                            combine_non_unique_arrays=combine_non_unique_arrays,
                                                            tree_path=tree_path,
                                                            use_provided_trees=True if tree_path is not None else False,
                                                            tree_save_path=os.path.join(
                                                                tree_save_path,
                                                                '0_forward'),
                                                            hide_unobserved_spacers=hide_unobserved_spacers,
                                                            selection_fct=selection_fct,
                                                            plot_order=plot_order,
                                                            finite_bl_at_root=finite_bl_at_root,
                                                            group_by=group_by,
                                                            significance_level=significance_level,
                                                            extend_branches=extend_branches,
                                                            tree_distance_function=tree_distance_function,
                                                            tree_lh_fct=tree_lh_fct,
                                                            tree_construction_method=tree_construction_method,
                                                            tree_gain_rate=tree_gain_rate,
                                                            tree_loss_rate=tree_loss_rate,
                                                            tree_alpha=tree_alpha,
                                                            alpha_bias_correction=alpha_bias_correction,
                                                            rho_bias_correction=rho_bias_correction,
                                                            rho_bias_failure_policy=rho_bias_failure_policy,
                                                            core_genome_trees=core_genome_trees,
                                                            save_reconstructed_events=save_reconstructed_events,
                                                            dpi=dpi,
                                                            figsize_rec=figsize_rec,
                                                            show_spacer_name_row=show_spacer_name_row,)
    dict_trees = dict_trees_forward
    _require_reconstructed_protocols(df_rec_protocol, df_rec_protocol_boring,
                                     os.path.join(save_path, '0_forward'))
    if not df_rec_protocol.empty:
        df_rec_protocol = df_rec_protocol.set_index('name')
    if not df_rec_protocol_boring.empty:
        df_rec_protocol_boring = df_rec_protocol_boring.set_index('name')
    if determine_orientation:
        if not os.path.exists(os.path.join(save_path, '0_reversed')):
            os.makedirs(os.path.join(save_path, '0_reversed'))

        for crispr_group in dict_crispr_groups_for_reverse.values():
            for ca in crispr_group.crispr_dict.values():
                ca.spacer_array = ca.spacer_array[::-1]
        dict_provided_aligned_arrays = {k: v[v.columns[::-1]] for k, v in dict_provided_aligned_arrays.items()}

        (df_rec_protocol_reversed,
         df_rec_protocol_reversed_boring,
         (dict_data_for_group_by_reversed, lh_fct), _, dict_trees_reversed,
         _, _) = run_reconstruction(rec_parameter_dict,
                                    dict_crispr_groups_for_reverse,
                                    save_path=os.path.join(
                                        save_path,
                                        '0_reversed'),
                                    plot_tree=plot_tree,
                                    lh_fct=lh_fct,
                                    logfile_path=logfile_path,
                                    logger=logger,
                                    do_show=do_show,
                                    combine_non_unique_arrays=combine_non_unique_arrays,
                                    tree_path=tree_path,
                                    use_provided_trees=True
                                    if tree_path is not None else False,
                                    tree_save_path=os.path.join(
                                        tree_save_path,
                                        '0_reversed'),
                                    hide_unobserved_spacers=hide_unobserved_spacers,
                                    selection_fct=selection_fct,
                                    plot_order=plot_order,
                                    finite_bl_at_root=finite_bl_at_root,
                                    group_by=group_by,
                                    significance_level=significance_level,
                                    extend_branches=extend_branches,
                                    tree_distance_function=tree_distance_function,
                                    tree_lh_fct=tree_lh_fct,
                                    tree_construction_method=tree_construction_method,
                                    tree_gain_rate=tree_gain_rate,
                                    tree_loss_rate=tree_loss_rate,
                                    tree_alpha=tree_alpha,
                                    provided_numbering=dict_provided_numbering,
                                    alpha_bias_correction=alpha_bias_correction,
                                    rho_bias_correction=rho_bias_correction,
                                    rho_bias_failure_policy=rho_bias_failure_policy,
                                    core_genome_trees=core_genome_trees,
                                    provided_aligned_arrays=dict_provided_aligned_arrays,
                                    provided_dict_duplicated_spacers=dict_provided_duplicated_spacers,
                                    save_reconstructed_events=save_reconstructed_events,
                                    dpi=dpi,
                                    figsize_rec=figsize_rec,
                                    show_spacer_name_row=show_spacer_name_row,)
        dict_trees = {'forward': dict_trees_forward,
                      'reversed': dict_trees_reversed}
        if not df_rec_protocol_reversed.empty:
            df_rec_protocol_reversed = df_rec_protocol_reversed.set_index('name')
        if not df_rec_protocol_reversed_boring.empty:
            df_rec_protocol_reversed_boring = df_rec_protocol_reversed_boring.set_index('name')

        (df_rec_protocol, df_oriented_rec_protocol,
         df_oriented_rec_protocol_w_trivial, dict_trees) = orientation_tools.compare_likelihoods_for_orientation(
            df_rec_protocol,
            df_rec_protocol_reversed,
            df_rec_protocol_boring,
            df_rec_protocol_reversed_boring,
            decision_boundary=orient_boundary, dict_trees=dict_trees)

        if group_by is not None:
            logger.info(f'Grouping trees for parameter estimation by: {group_by}')
            if isinstance(group_by, (int, np.integer)):
                val_group_by = integer_group_labels(df_oriented_rec_protocol.shape[0], group_by)
            elif isinstance(group_by, pd.Series):
                # solve this issue here pooling and group_by do not agree in dimensions!!!!
                val_group_by = [group_by[idx] for idx in df_oriented_rec_protocol.index]
            else:
                val_group_by = df_oriented_rec_protocol[group_by]
            ls_data_for_lh_ratio = [
                dict_data_for_group_by_reversed[idx]
                if df_oriented_rec_protocol.loc[idx, 'predicted orientation'] == 'Reverse'
                else dict_data_for_group_by[idx]
                for idx in df_oriented_rec_protocol.index
            ]
            group_names, groups_lh_0, groups_lh_1, \
                groups_result_0, groups_result_1 = pooling_for_parameter_estimation(ls_data_for_lh_ratio,
                                                                                    group_by=val_group_by,
                                                                                    give_lh_fct=lh_fct)

            ls_loss_rates_0 = [g.x for g in groups_result_0]
            ls_loss_rates_1 = [g.x[0] for g in groups_result_1]
            ls_alphas_1 = [g.x[1] for g in groups_result_1]
            ls_ln_lh_ratio = [2 * (g_0.fun - g_1.fun) for g_0, g_1 in zip(groups_result_0, groups_result_1)]
            ls_test_results = []
            ls_confidences = []
            ls_sig_values = []
            for llh in ls_ln_lh_ratio:
                test_result, confidence = stats.test_significance_ratio_chi2(llh, significance_level)
                ls_test_results.append(test_result)
                ls_confidences.append(confidence)
                ls_sig_values.append(significance_level)

            df_group_protocol = pd.DataFrame(list(zip(group_names, groups_lh_0, groups_lh_1,
                                                      ls_loss_rates_0, ls_loss_rates_1,
                                                      ls_alphas_1, ls_ln_lh_ratio, ls_test_results,
                                                      ls_confidences, ls_sig_values,
                                                      )),
                                             columns=['group_by', 'lh_idm', 'lh_bdm', 'deletion_rate_idm',
                                                      'deletion_rate_bdm',
                                                      'alpha_bdm', 'test statistic (-2*ln_lh_ratio)',
                                                      'test result', 'chi2_quantile',
                                                      'significance level'])
            df_group_protocol.to_csv(os.path.join(save_path, '0_final_group_by_protocol.csv'))

        _write_orientation_outputs(
            save_path, df_rec_protocol, df_oriented_rec_protocol,
            dict_trees_forward, dict_trees_reversed, dict_trees,
        )
        df_rec_protocol.to_csv(os.path.join(save_path, '_'.join(['0_protocol_w_orientation']) + '.csv'))
        df_rec_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_protocol_w_orientation']) + '.pkl'))
        df_oriented_rec_protocol.to_csv(os.path.join(save_path, '0_final_oriented_protocol.csv'))
        df_oriented_rec_protocol.to_pickle(os.path.join(save_path, '0_final_oriented_protocol.pkl'))
    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'final_dict_nwk_trees.json'), 'w') as f:
            json.dump(dict_trees, f)
    return


def run_simulation_and_reconstruction(sim_parameter_dict, rec_parameter_dict, save_path=None,
                                      rec_save_path=None, plot_tree=True,
                                      sim_save_path=None, sim_plot_tree=True,
                                      sim_logfile_path=None,
                                      logfile_path=None,
                                      do_show=False,
                                      sim_as_rec=False,
                                      hide_unobserved_spacers=False,
                                      selection_fct=None,
                                      plot_order=True,
                                      significance_level=0.05,
                                      alignment_options=None,
                                      load_sim_from_pkl=False,
                                      finite_bl_at_root=True,
                                      sim_simple_alignment=False,
                                      group_by=None,
                                      extend_branches=False,
                                      use_sim_tree=True,
                                      tree_distance_function='breakpoint',
                                      tree_gain_rate=None,
                                      tree_loss_rate=None,
                                      tree_alpha=None,
                                      tree_save_path=None,
                                      determine_orientation=False,
                                      randomize_orientation=False,
                                      orient_decision_boundary=0,
                                      alpha_bias_correction=False,
                                      rho_bias_correction=False,
                                      dpi=90,
                                      figsize_rec=(None, None, 'px'),
                                      ):
    """
    This function runs a simulation and a reconstruction and compares the results.
    :param figsize_rec:
    :param rho_bias_correction:
    :param alpha_bias_correction:
    :param sim_parameter_dict:
    :param rec_parameter_dict:
    :param save_path:
    :param rec_save_path:
    :param plot_tree:
    :param sim_save_path:
    :param sim_plot_tree:
    :param sim_logfile_path:
    :param logfile_path:
    :param do_show:
    :param sim_as_rec:
    :param hide_unobserved_spacers:
    :param selection_fct:
    :param plot_order:
    :param significance_level:
    :param alignment_options:
    :param load_sim_from_pkl:
    :param finite_bl_at_root:
    :param sim_simple_alignment:
    :param group_by:
    :param extend_branches:
    :param use_sim_tree:
    :param tree_distance_function:
    :param tree_gain_rate:
    :param tree_loss_rate:
    :param tree_alpha:
    :param tree_save_path:
    :param determine_orientation:
    :param randomize_orientation:
    :param orient_decision_boundary:
    :param load_sim_param_from_pkl:
    :return:
    """
    logger = misc.create_logger('Simulation -> Reconstruction', logging.INFO, outfile=logfile_path)

    if load_sim_from_pkl:
        logger.info(f'Loading simulations from {sim_save_path}...')
        with open(os.path.join(sim_save_path, 'dict_crispr_groups.pkl'), 'rb') as f:
            dict_crispr_groups = pickle.load(f)
        df_sim_protocol = pd.read_pickle(os.path.join(sim_save_path, '_'.join(['0_sim_protocol']) + '.pkl'))
    else:
        logger.info('Starting simulations...')
        from sp_model.run_simulation import run_simulation

        df_sim_protocol, dict_crispr_groups, ls_sim_m = run_simulation(sim_parameter_dict, save_path=sim_save_path,
                                                                       plot_tree=sim_plot_tree,
                                                                       logfile_path=sim_logfile_path,
                                                                       sim_as_rec=sim_as_rec,
                                                                       randomize_orientation=randomize_orientation,
                                                                       )
        # Alignment
        if sim_as_rec or sim_simple_alignment:
            logger.info('Simulation is used as Reconstruction or sim_simple_alignment was chosen -> '
                        'alignment is done based on simulated top_order.')
            dict_crispr_groups = align_crispr_groups_for_sim_as_rec(dict_crispr_groups)
        else:
            logger.info('Starting alignment...')
            dict_crispr_groups = mafft_integration.align_crispr_groups(WORK_PATH, dict_crispr_groups,
                                                                       mafft_options=alignment_options,
                                                                       logger=logger, )

        if sim_save_path:
            with open(os.path.join(sim_save_path, 'dict_crispr_groups.pkl'), 'wb') as f:
                pickle.dump(dict_crispr_groups, f)

    # If we want to use simulated parameters for grouping.
    if group_by is not None:
        if isinstance(group_by, str):
            if group_by.split('|')[0] == 'sim':
                group_by = df_sim_protocol.set_index('name')[''.join(group_by.split('|')[1:])]
    logger.info('Starting reconstruction...')

    if not os.path.exists(os.path.join(rec_save_path, '0_forward')):
        os.makedirs(os.path.join(rec_save_path, '0_forward'))

    dict_crispr_groups_for_reverse = copy.deepcopy(dict_crispr_groups) if determine_orientation else {}

    (df_rec_protocol, df_rec_protocol_boring, (
        dict_data_for_group_by, lh_fct),
     dict_provided_numbering, dict_trees_forward,
     dict_provided_aligned_arrays,
     dict_provided_duplicated_spacers) = run_reconstruction(
        rec_parameter_dict,
        dict_crispr_groups,
        save_path=os.path.join(
            rec_save_path,
            '0_forward'),
        plot_tree=plot_tree,
        logfile_path=logfile_path,
        do_show=do_show,
        combine_non_unique_arrays=False,
        tree_path=os.path.join(
            sim_save_path,
            'dict_nwk_trees.pkl'),
        use_provided_trees=use_sim_tree,
        tree_save_path=os.path.join(
            tree_save_path,
            '0_forward') if tree_save_path is not None else None,
        sim_as_rec=sim_as_rec,
        hide_unobserved_spacers=hide_unobserved_spacers,
        selection_fct=selection_fct,
        plot_order=plot_order,
        significance_level=significance_level,
        finite_bl_at_root=finite_bl_at_root,
        logger=logger,
        group_by=group_by,
        extend_branches=extend_branches,
        tree_distance_function=tree_distance_function,
        tree_gain_rate=tree_gain_rate,
        tree_loss_rate=tree_loss_rate,
        tree_alpha=tree_alpha,
        alpha_bias_correction=alpha_bias_correction,
        rho_bias_correction=rho_bias_correction,
        dpi=dpi,
        figsize_rec=figsize_rec,
    )

    dict_trees = dict_trees_forward
    _require_reconstructed_protocols(df_rec_protocol, df_rec_protocol_boring,
                                     os.path.join(rec_save_path, '0_forward'))
    if not df_rec_protocol.empty:
        df_rec_protocol = df_rec_protocol.set_index('name')
    if not df_rec_protocol_boring.empty:
        df_rec_protocol_boring = df_rec_protocol_boring.set_index('name')

    if determine_orientation:
        if not os.path.exists(os.path.join(rec_save_path, '0_reversed')):
            os.makedirs(os.path.join(rec_save_path, '0_reversed'))

        for crispr_group in dict_crispr_groups_for_reverse.values():
            for ca in crispr_group.crispr_dict.values():
                ca.spacer_array = ca.spacer_array[::-1]
        dict_provided_aligned_arrays = {k: v[v.columns[::-1]] for k, v in dict_provided_aligned_arrays.items()}

        (df_rec_protocol_reversed,
         df_rec_protocol_reversed_boring,
         (dict_data_for_group_by_reversed, lh_fct), _, dict_tree_reversed,
         _, _) = run_reconstruction(rec_parameter_dict,
                                    dict_crispr_groups_for_reverse,
                                    save_path=os.path.join(
                                        rec_save_path,
                                        '0_reversed'),
                                    plot_tree=plot_tree,
                                    logfile_path=logfile_path,
                                    do_show=do_show,
                                    combine_non_unique_arrays=False,
                                    tree_path=os.path.join(
                                        sim_save_path,
                                        'dict_nwk_trees.pkl'),
                                    use_provided_trees=use_sim_tree,
                                    tree_save_path=os.path.join(
                                        tree_save_path,
                                        '0_reversed') if tree_save_path is not None else None,
                                    sim_as_rec=sim_as_rec,
                                    hide_unobserved_spacers=hide_unobserved_spacers,
                                    selection_fct=selection_fct,
                                    plot_order=plot_order,
                                    significance_level=significance_level,
                                    finite_bl_at_root=finite_bl_at_root,
                                    logger=logger,
                                    group_by=group_by,
                                    extend_branches=extend_branches,
                                    tree_distance_function=tree_distance_function,
                                    tree_gain_rate=tree_gain_rate,
                                    tree_loss_rate=tree_loss_rate,
                                    tree_alpha=tree_alpha,
                                    alpha_bias_correction=alpha_bias_correction,
                                    rho_bias_correction=rho_bias_correction,
                                    provided_aligned_arrays=dict_provided_aligned_arrays,
                                    provided_dict_duplicated_spacers=dict_provided_duplicated_spacers,
                                    dpi=dpi,
                                    figsize_rec=figsize_rec,
                                    )
        dict_trees = {'forward': dict_trees_forward, 'reversed': dict_tree_reversed}

        if not df_rec_protocol_reversed.empty:
            df_rec_protocol_reversed = df_rec_protocol_reversed.set_index('name')
        if not df_rec_protocol_reversed_boring.empty:
            df_rec_protocol_reversed_boring = df_rec_protocol_reversed_boring.set_index('name')

        (df_rec_protocol, df_oriented_rec_protocol,
         _, dict_trees) = orientation_tools.compare_likelihoods_for_orientation(
            df_rec_protocol,
            df_rec_protocol_reversed,
            df_rec_protocol_boring,
            df_rec_protocol_reversed_boring,
            decision_boundary=orient_decision_boundary, dict_trees=dict_trees)

        if group_by is not None:
            logger.info(f'Grouping trees for parameter estimation by: {group_by}')
            if isinstance(group_by, (int, np.integer)):
                val_group_by = integer_group_labels(df_oriented_rec_protocol.shape[0], group_by)
            elif isinstance(group_by, pd.Series):
                val_group_by = [group_by[idx] for idx in df_oriented_rec_protocol.index]
            else:
                val_group_by = df_oriented_rec_protocol[group_by]
            ls_data_for_lh_ratio = []
            filtered_val_group_by = []
            for i, idx in enumerate(df_oriented_rec_protocol.index):
                if df_oriented_rec_protocol.loc[idx, 'predicted orientation'] == 'Reverse':
                    if idx in dict_data_for_group_by_reversed:
                        ls_data_for_lh_ratio.append(dict_data_for_group_by_reversed[idx])
                        filtered_val_group_by.append(val_group_by[i])
                else:
                    if idx in dict_data_for_group_by:
                        ls_data_for_lh_ratio.append(dict_data_for_group_by[idx])
                        filtered_val_group_by.append(val_group_by[i])

            group_names, groups_lh_0, groups_lh_1, \
                groups_result_0, groups_result_1 = pooling_for_parameter_estimation(ls_data_for_lh_ratio,
                                                                                    group_by=filtered_val_group_by,
                                                                                    give_lh_fct=lh_fct)

            ls_loss_rates_0 = [g.x for g in groups_result_0]
            ls_loss_rates_1 = [g.x[0] for g in groups_result_1]
            ls_alphas_1 = [g.x[1] for g in groups_result_1]
            ls_ln_lh_ratio = [2 * (g_0.fun - g_1.fun) for g_0, g_1 in zip(groups_result_0, groups_result_1)]
            ls_test_results = []
            ls_confidences = []
            ls_sig_values = []
            for llh in ls_ln_lh_ratio:
                test_result, confidence = stats.test_significance_ratio_chi2(llh, significance_level)
                ls_test_results.append(test_result)
                ls_confidences.append(confidence)
                ls_sig_values.append(significance_level)

            df_group_protocol = pd.DataFrame(list(zip(group_names, groups_lh_0, groups_lh_1,
                                                      ls_loss_rates_0, ls_loss_rates_1,
                                                      ls_alphas_1, ls_ln_lh_ratio, ls_test_results,
                                                      ls_confidences, ls_sig_values,
                                                      )),
                                             columns=['group_by', 'lh_0', 'lh_1', 'loss_rate_0', 'loss_rate_1',
                                                      'alpha_1', '-2*ln_lh_ratio', 'test result', 'test_confidence',
                                                      'significance value'])
            df_group_protocol.to_csv(os.path.join(save_path, '0_final_group_by_protocol.csv'))

        df_rec_protocol.to_csv(os.path.join(rec_save_path, '_'.join(['0_protocol_orientation']) + '.csv'))
        df_rec_protocol.to_pickle(os.path.join(rec_save_path, '_'.join(['0_protocol_orientation']) + '.pkl'))
        df_oriented_rec_protocol.to_csv(os.path.join(rec_save_path, '0_final_oriented_protocol.csv'))
        df_oriented_rec_protocol.to_pickle(os.path.join(rec_save_path, '0_final_oriented_protocol.pkl'))

    # trying concat for now, might want to merge/join, especially, if keys overlap
    df_sim_protocol = df_sim_protocol.set_index('name')
    if not determine_orientation:
        df_rec_protocol = pd.concat([df_rec_protocol, df_rec_protocol_boring], axis=0)

    df_rec_protocol = df_oriented_rec_protocol if determine_orientation else df_rec_protocol
    df_sim_rec_protocol = pd.concat([df_sim_protocol, df_rec_protocol], axis=1, join='inner')

    protocol_path = os.path.join(save_path, '_'.join(['0_sim_rec_protocol']))
    logger.info(f'Finished simulation and reconstruction, saving to {protocol_path}')
    df_sim_rec_protocol.to_csv(os.path.join(save_path, '_'.join(['0_sim_rec_protocol']) + '.csv'))
    df_sim_rec_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_sim_rec_protocol']) + '.pkl'))

    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'final_dict_nwk_trees.pkl'), 'wb') as f:
            pickle.dump(dict_trees, f)
