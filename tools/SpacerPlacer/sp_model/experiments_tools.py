import copy
import json

import pandas as pd
import numpy as np
import os
import pickle

from sp_model.helpers import misc, import_data
from sp_model import model_tools
from sp_model.coalescent_tree_sim import mult_rnd_coalescent_tree, rnd_coalescent_tree
from sp_model import distance_tree_constructor
from sp_model import mafft_integration
from sp_model.data_classes.crisprdata import CRISPRGroup, CRISPRArray, Species
from sp_model.data_classes.advanced_tree import AdvancedTree

GROUPS_TO_EXCLUDE = ['g_39', 'g_284', 'g_690', 'g_140', 'g_719', 'g_931', 'g_598',
                     'g_942']


def integer_group_labels(row_count, group_size):
    """Return stable chunk labels, folding a short remainder into the last full chunk."""
    if isinstance(row_count, bool) or not isinstance(row_count, (int, np.integer)) or row_count < 0:
        raise ValueError('row_count must be a non-negative integer')
    if isinstance(group_size, bool) or not isinstance(group_size, (int, np.integer)) or group_size <= 0:
        raise ValueError('group_size must be a positive integer')
    if row_count == 0:
        return []

    full_group_count = max(row_count // group_size, 1)
    return [min(index // group_size, full_group_count - 1) for index in range(row_count)]


def load_align_pickled_data(path, mafft_options=None, exclude_files=None, save_path=None):
    """
    :param path:
    :param mafft_options:
    :param exclude_files:
    :param save_path:
    :return:
    """
    ls_files = []
    for entry in sorted(os.listdir(path)):
        if os.path.isdir(os.path.join(path, entry)):
            continue
        ls_files.append(entry)
    dict_crispr_groups = {}
    if exclude_files:
        ls_files = [file for file in ls_files if file not in exclude_files]
    for file in ls_files:
        ls_arrays, ls_names, head = mafft_integration.pickled_group_to_mafft(path, [file],
                                                                             options=mafft_options)
        ls_crispr_arrays = []
        for name, array in zip(ls_names, ls_arrays):
            crispr_array = CRISPRArray(name, '', None, 4, head[1],
                                       'chromosome', 'Bacteria', head[0], array, head[2],
                                       cas_genes=[], all_cas_types=[], all_cas_genes=[],
                                       species=Species('', 'Bacteria', 'sim', 'sim', 'sim', 'sim',
                                                                      'sim', 'sim'),
                                       species_fc='sim')
            ls_crispr_arrays.append(crispr_array)
        crispr_group = CRISPRGroup(head[1], ls_crispr_arrays, name=file)
        crispr_group.convert_arrays_from_mafft_fmt()
        dict_crispr_groups[file] = crispr_group

    if save_path:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        with open(os.path.join(save_path, '0_aligned_crisprgroups.pkl'), 'wb') as f:
            pickle.dump(dict_crispr_groups, f)
    return dict_crispr_groups


def give_count_duplicates_in_list(ls):
    """
    Returns a list of all duplicates in the list
    :param ls:
    :return:
    """
    dict_duplicates = {}
    ls_unique = []
    for item in ls:
        if item in ls_unique:
            dict_duplicates[item] = dict_duplicates.get(item, 1) + 1
        else:
            ls_unique.append(item)
    return dict_duplicates


def name_inner_nodes_if_unnamed(tree):
    names = [node.name for node in tree.find_clades() if node.name is not None]
    if len(names) != len(set(names)):
        raise ValueError('Tree node names must be unique')
    used_names = set(names)
    n_i = 0
    for node in tree.get_nonterminals():
        if node.name is None:
            while 'in_' + str(n_i) in used_names:
                n_i += 1
            node.name = 'in_' + str(n_i)
            used_names.add(node.name)
            n_i += 1
    return tree


def tree_handling(tree, crispr_group, name_inner_nodes=False, bl_eps=0):
    dict_duplicates = give_count_duplicates_in_list(crispr_group.get_ls_acc_nums())
    duplicate_to_ls_arrays = {key: crispr_group.get_ls_arrays_by_acc_num(key) for key in dict_duplicates.keys()}
    # acc_num: name of arrays that need to be pruned
    dict_to_prune = {}
    # acc_num to names of arrays that need to use the same leaf... (or combine them?)
    dict_acc_num_ls_duplicates = {}
    for key, ls_ca in duplicate_to_ls_arrays.items():
        ls_ca_w_o_gaps = [[a for a in x.spacer_array if a != '-'] for x in ls_ca]
        if all([ls_ca_w_o_gaps[0] == x for x in ls_ca_w_o_gaps]):
            dict_acc_num_ls_duplicates[key] = [x.name for x in ls_ca]
        else:
            dict_to_prune[key] = [x.name for x in ls_ca]

    ls_crispr_arrays_to_delete = []
    for name, ca in crispr_group.crispr_dict.items():
        for leaf in tree.get_terminals():
            if leaf.name == ca.acc_num:
                if leaf.name in dict_to_prune.keys():
                    # print('pruning', leaf.name)
                    if len(list(tree.get_terminals())) > 1:
                        tree.prune(leaf)
                    ls_crispr_arrays_to_delete.append(ca.acc_num)
                # This splits the tree leafs for equal CRISPR arrays in the same strain...
                elif leaf.name in dict_acc_num_ls_duplicates.keys():
                    leaf.split(n=len(dict_acc_num_ls_duplicates[leaf.name]), branch_length=bl_eps)
                    for child in leaf.clades:
                        child.name = dict_acc_num_ls_duplicates[leaf.name].pop()
                else:
                    leaf.name = name
    for acc_num in ls_crispr_arrays_to_delete:
        crispr_group.drop_array_by_acc_num(acc_num)
    # Naming inner nodes...
    if name_inner_nodes:
        tree = name_inner_nodes_if_unnamed(tree)
    tree = AdvancedTree(tree, rooted=True, model_name=tree.model_name)
    return tree, crispr_group


def load_align_single_fasta(data_path, work_path, mafft_options=None, group_name='g', logger=None, save_path=None,
                            seed=None):
    dict_arrays = read_spacer_arrays(data_path)
    ls_array_names = list(dict_arrays)
    ls_spacer_arrays = list(dict_arrays.values())
    if len(ls_array_names) < 2:
        raise ValueError(f"Group '{group_name}' must contain at least two arrays")
    ls_crispr_arrays = []
    for name, array in zip(ls_array_names, ls_spacer_arrays):
        crispr_array = CRISPRArray(name, None, None, None,
                                   None, None, None, None, array,
                                   None,
                                   cas_genes=None, all_cas_types=None, all_cas_genes=None,
                                   species=None,
                                   species_fc=None)
        ls_crispr_arrays.append(crispr_array)
    crispr_group = CRISPRGroup(None, ls_crispr_arrays, name=group_name)
    dict_crispr_group = {crispr_group.name: crispr_group}
    dict_crispr_group = mafft_integration.align_crispr_groups(work_path, dict_crispr_group, mafft_options=mafft_options,
                                                              logger=logger, seed=seed)
    if save_path:
        with open(save_path, 'wb') as f:
            pickle.dump(dict_crispr_group, f)
    return list(dict_crispr_group.values())[0]


def read_spacer_arrays(data_path):
    """Read observed spacer IDs without alignment or executable data loading."""
    dict_fasta = read_fasta(data_path)
    arrays = {}
    for key, value in dict_fasta.items():
        array_name = key.split(',', maxsplit=1)[0].strip()
        if not array_name:
            raise ValueError(f"FASTA record '{key}' has an empty array name")
        if array_name in arrays:
            raise ValueError(f"Duplicate array name after metadata removal: {array_name}")
        spacer_array = [spacer.strip() for spacer in value.split(',')]
        if any(not spacer for spacer in spacer_array):
            raise ValueError(f"Array '{array_name}' contains an empty spacer ID")
        if any(any(character.isspace() for character in spacer) for spacer in spacer_array):
            raise ValueError(f"Array '{array_name}' contains whitespace in a spacer ID")
        if any(spacer in {'-', '--'} for spacer in spacer_array):
            raise ValueError(f"Array '{array_name}' contains a reserved alignment gap as a spacer ID")
        arrays[array_name] = spacer_array
    return arrays


def write_aligned_groups_json(dict_crispr_groups, output_path):
    """Persist aligned groups in a portable, non-executable interchange format."""
    document = {
        'schema': {'name': 'spacerplacer-aligned-groups', 'version': 1},
        'groups': [
            {
                'name': group.name,
                'repeat': group.repeat,
                'arrays': [
                    {'name': name, 'spacers': crispr_array.spacer_array}
                    for name, crispr_array in group.crispr_dict.items()
                ],
            }
            for group in dict_crispr_groups.values()
        ],
    }
    with open(output_path, 'w', encoding='utf-8') as output_file:
        json.dump(document, output_file, indent=2, allow_nan=False)
        output_file.write('\n')


def load_aligned_groups_json(input_path):
    """Load and validate a portable alignment cache produced by SpacerPlacer."""
    try:
        with open(input_path, 'r', encoding='utf-8') as input_file:
            document = json.load(input_file)
    except json.JSONDecodeError as error:
        raise ValueError(f'Invalid aligned-groups JSON: {input_path}') from error

    schema = {
        'name': 'spacerplacer-aligned-groups',
        'version': 1,
    }
    if not isinstance(document, dict) or document.get('schema') != schema:
        raise ValueError('Unsupported aligned-groups JSON schema')
    group_entries = document.get('groups')
    if not isinstance(group_entries, list) or not group_entries:
        raise ValueError('Aligned-groups JSON must contain at least one group')

    groups = {}
    for group_entry in group_entries:
        if not isinstance(group_entry, dict):
            raise ValueError('Every aligned group must be a JSON object')
        group_name = group_entry.get('name')
        if not isinstance(group_name, str) or not group_name.strip():
            raise ValueError('Every aligned group must have a non-empty string name')
        if group_name in groups:
            raise ValueError(f'Duplicate aligned group name: {group_name}')
        repeat = group_entry.get('repeat')
        if repeat is not None and not isinstance(repeat, str):
            raise ValueError(f'Invalid repeat for aligned group {group_name}')
        array_entries = group_entry.get('arrays')
        if not isinstance(array_entries, list) or len(array_entries) < 2:
            raise ValueError(f'Aligned group {group_name} must contain at least two arrays')

        arrays = []
        array_names = set()
        alignment_length = None
        for array_entry in array_entries:
            if not isinstance(array_entry, dict):
                raise ValueError(f'Invalid array entry in aligned group {group_name}')
            array_name = array_entry.get('name')
            spacers = array_entry.get('spacers')
            if not isinstance(array_name, str) or not array_name.strip() or array_name in array_names:
                raise ValueError(f'Invalid or duplicate array name in aligned group {group_name}')
            if (
                not isinstance(spacers, list)
                or any(
                    not isinstance(spacer, str) or not spacer
                    for spacer in spacers
                )
            ):
                raise ValueError(f'Invalid spacer alignment for {group_name}/{array_name}')
            if alignment_length is None:
                alignment_length = len(spacers)
            elif len(spacers) != alignment_length:
                raise ValueError(f'Ragged spacer alignment in group {group_name}')
            array_names.add(array_name)
            arrays.append(CRISPRArray(
                array_name, None, None, None, None, None, None, None, spacers, None,
                cas_genes=None, all_cas_types=None, all_cas_genes=None, species=None, species_fc=None,
            ))
        groups[group_name] = CRISPRGroup(repeat, arrays, name=group_name)
    return groups


def read_fasta(data_path):
    """Read spacer-array FASTA, supporting wrapped records and strict validation."""
    records = {}
    current_header = None
    sequence_lines = []

    def store_record():
        if current_header is None:
            return
        sequence = "".join(sequence_lines).strip()
        if not sequence:
            raise ValueError(f"FASTA record '{current_header}' has no spacer IDs")
        if current_header in records:
            raise ValueError(f"Duplicate FASTA header: {current_header}")
        records[current_header] = sequence

    with open(data_path, "r", encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store_record()
                current_header = line[1:].strip()
                sequence_lines = []
                if not current_header:
                    raise ValueError(f"Empty FASTA header on line {line_number}")
            else:
                if current_header is None:
                    raise ValueError(
                        f"Sequence data appears before the first FASTA header on line {line_number}"
                    )
                sequence_lines.append(line)
    store_record()
    if not records:
        raise ValueError(f"No FASTA records found in {data_path}")
    return records


def align_crispr_groups_for_sim_as_rec(dict_crispr_groups):
    for crispr_group in dict_crispr_groups.values():
        dict_arrays_as_list = crispr_group.get_arrays_as_lists()
        top_order = crispr_group.top_order
        ls_names = list(dict_arrays_as_list.keys())
        ls_arrays = list(dict_arrays_as_list.values())

        ls_aligned_arrays = []
        for array in ls_arrays:
            aligned_array = [str(s) if str(s) in array else '-' for s in top_order]
            ls_aligned_arrays.append(aligned_array)
        crispr_group.update_spacer_arrays_by_ls_arrays_as_list(ls_names, ls_aligned_arrays)
    return dict_crispr_groups


def filter_fct(dict_protocol):
    """
    :param dict_protocol:
    :return:
    """
    # return True
    return True if dict_protocol['all max block deletion lengths'] else False


def construct_tree(ls_array_names, ls_arrays, group_name, logger=None, distance_fct='breakpoint', tree_save_path=None,
                   gain_rate=None, loss_rate=None, alpha=None, provided_lh_fct=None, tree_construction_method='upgma'):
    """
    Constructs an upgma tree based on different distance functions.
    :param tree_construction_method:
    :param ls_array_names:
    :param ls_arrays:
    :param group_name:
    :param logger:
    :param distance_fct:
    :param tree_save_path:
    :param gain_rate:
    :param loss_rate:
    :param alpha:
    :param provided_lh_fct:
    :return:
    """
    if distance_fct == 'likelihood':
        distance = distance_tree_constructor.LikelihoodDistance(gain_rate, loss_rate, alpha,
                                                                provided_lh_fct=provided_lh_fct, skip_letters={'-'})
    else:
        raise NotImplementedError('Distance function not implemented!')
    tree_constructor = distance_tree_constructor.FixedDistanceTreeConstructor(distance_calculator=distance,
                                                                              method=tree_construction_method)
    # tree_constructor = DistanceTreeConstructor(distance_calculator=distance, method=tree_construction_method)

    if logger is not None:
        logger.info(f'Constructing tree of {group_name}...')
    tree = tree_constructor.build_tree([ls_array_names, ls_arrays])


    # for clade in tree.find_clades():
    #     print(clade.name)
    #     print(clade.branch_length)
    if tree_save_path:
        if not os.path.exists(tree_save_path):
            os.makedirs(tree_save_path)
        save_name = os.path.join(tree_save_path, group_name + '_tree.nwk')
        logger.info(f'Tree saved to: {save_name}')
        misc.write_tree_to_file(tree, save_name)
    return tree


def multiple_selection_fcts(array, selection_criteria, selection_parameters=None):
    skip = False
    reason = None
    if not selection_criteria:
        skip = False
    else:
        for selection_criterium, selection_parameter in zip(selection_criteria, selection_parameters):
            skip, reason = selection_fct(array, selection_criterium, selection_parameter=selection_parameter)
            if skip:
                return skip, reason
    return skip, reason


def selection_fct(group, selection_criterium, selection_parameter=None):
    if selection_criterium is None:
        skip = False
    elif selection_criterium == 'maximum_unique_array_number':
        ls_arrays = []
        for array in group.crispr_dict.values():
            if array.spacer_array in ls_arrays:
                continue
            else:
                ls_arrays.append(array.spacer_array)
        skip = True if len(ls_arrays) > selection_parameter else False
    elif selection_criterium == 'maximum_array_number':
        skip = True if len(group.crispr_dict) > selection_parameter else False
    elif selection_criterium == 'minimum_array_number':
        skip = True if len(group.crispr_dict) < selection_parameter else False
    elif selection_criterium == 'empty_arrays':
        skip = True
        for array in group.crispr_dict.values():
            if set(array.spacer_array) not in [{'-'}, set()]:
                skip = False
                break
    elif selection_criterium == 'crispr_type':
        if group.cas_type.lower() in selection_parameter:
            skip = False
        else:
            skip = True
    elif selection_criterium == 'selection of groups':
        if group.name in selection_parameter:
            skip = False
        else:
            skip = True
    elif selection_criterium == 'selection not in groups':
        if group.name in selection_parameter:
            skip = True
        else:
            skip = False
    else:
        raise NotImplementedError('Selection function not implemented!')
    reason = selection_criterium
    return skip, reason


def load_arrays_from_pkl(data_path):
    dict_crispr_arrays = import_data.load_crispr_arrays_from_pickle(os.path.join(data_path))
    return dict_crispr_arrays


def expand_sim_parameters_based_on_import(parameter_dict):
    """
    needs import_path for parameters, trees_path for trees, and repeat_runs for number of repeated runs.
    :param parameter_dict:
    :return:
    """
    import_path = parameter_dict['import_path']
    with open(os.path.join(import_path), 'rb') as f:
        df_imported_parameters = pd.read_pickle(f)
        ls_names = df_imported_parameters.index.values
        ls_names = [name for name in ls_names if name not in GROUPS_TO_EXCLUDE]
        df_imported_parameters = df_imported_parameters.loc[ls_names]
        # df_imported_parameters = df_imported_parameters.set_index('name')

    ls_avg_array_length = list(df_imported_parameters['avg array length'].values)
    if parameter_dict['deletion_model'] == 'independent':
        ls_gtr_loss_rate = list(df_imported_parameters['loss_rate_0'].values)
        ls_deletion_model = ['independent'] * len(ls_gtr_loss_rate)
        ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_length,
                                                                  ls_gtr_loss_rate,
                                                                  'independent'
                                                                  )
        ls_deletion_model_parameter = [None] * len(ls_gtr_loss_rate)
    elif parameter_dict['deletion_model'] == 'block':
        ls_gtr_loss_rate = list(df_imported_parameters['loss_rate_1'].values)
        ls_deletion_model = ['block'] * len(ls_gtr_loss_rate)
        ls_deletion_model_parameter = list(df_imported_parameters['alpha_1'].values)
        ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_length,
                                                                  ls_gtr_loss_rate,
                                                                  'block',
                                                                  ls_deletion_model_parameter,
                                                                  )
    extended_ls_tree_names = []
    extended_ls_trees = []
    extended_ls_gtr_gain_rate = []
    extended_ls_gtr_loss_rate = []
    extended_ls_deletion_model = []
    extended_ls_deletion_model_parameter = []
    extended_ls_repeat_runs = []
    repeat_runs = parameter_dict['repeat_runs']
    extended_ls_sim_seed = list(range(parameter_dict['sim_seed'],
                                      parameter_dict['sim_seed'] + len(ls_names) * repeat_runs))

    nb_leafs_col = df_imported_parameters['nb leafs (after combining)']
    tree_sim_per_repeat = bool(parameter_dict.get('tree_sim_per_repeat_run', False))
    if parameter_dict.get('trees_path', False):
        trees_path = parameter_dict['trees_path']
        with open(trees_path, 'rb') as f:
            dict_trees = pickle.load(f)
    elif not tree_sim_per_repeat:
        ls_tree_sim_seed = [parameter_dict['tree_seed'] + i for i in range(len(ls_names))] \
            if 'tree_seed' in parameter_dict else [None] * len(ls_names)
        dict_trees = {ls_names[i]: rnd_coalescent_tree(nb_leafs_col.iloc[i], rnd_bl=True, seed=ls_tree_sim_seed[i],
                                                       tree_name=ls_names[i]) for i in range(len(ls_names))}
    else:
        dict_trees = None
        ls_tree_sim_seed = [parameter_dict['tree_seed'] + i for i in range(repeat_runs * len(ls_names))] \
            if 'tree_seed' in parameter_dict else [None] * (repeat_runs * len(ls_names))

    if dict_trees is not None:
        missing_trees = set(ls_names).difference(dict_trees)
        if missing_trees:
            raise ValueError(f'Missing trees for imported groups: {sorted(missing_trees)}')

    for r in range(repeat_runs):
        new_ls_tree_names = [name + '-' + str(r) for name in ls_names]
        extended_ls_tree_names += new_ls_tree_names
        if dict_trees is not None:
            new_ls_trees = [
                import_data.load_single_tree_from_string(tree) if isinstance(tree, str) else copy.deepcopy(tree)
                for tree in (dict_trees[name] for name in ls_names)
            ]
        else:
            new_ls_trees = [rnd_coalescent_tree(nb_leafs_col.loc[name],
                                                rnd_bl=True,
                                                seed=ls_tree_sim_seed[r * len(ls_names) + i],
                                                tree_name=new_ls_tree_names[i])
                            for i, name in enumerate(ls_names)]
        for i, t in enumerate(new_ls_trees):
            t.name = new_ls_tree_names[i]
        extended_ls_repeat_runs += [r] * len(new_ls_tree_names)
        extended_ls_trees += new_ls_trees
        extended_ls_gtr_gain_rate += ls_gtr_gain_rate
        extended_ls_gtr_loss_rate += ls_gtr_loss_rate
        extended_ls_deletion_model += ls_deletion_model
        extended_ls_deletion_model_parameter += ls_deletion_model_parameter

    return extended_ls_tree_names, extended_ls_trees, extended_ls_gtr_gain_rate, extended_ls_gtr_loss_rate, \
        extended_ls_deletion_model, extended_ls_deletion_model_parameter, extended_ls_sim_seed, extended_ls_repeat_runs


def expand_sim_parameters(parameter_dict):
    """
    Should include: nb_trees, nb_leafs, deletion_model (string or list), gtr_loss_rate, deletion_model
    tree_seed for consistency
    deletion_model_parameter if used
    parameters_path if used
    lr_dependent_on_length_alpha desired_length if used then no gtr_gain_rate
    trees_path if used then no rnd_bl
    :param parameter_dict:
    :return: dictionary with trees, gtr_gain/loss_rate, deletion_model, deletion_model_parameter
    """
    nb_trees = parameter_dict['nb_trees']
    nb_leafs = parameter_dict['nb_leafs']
    ls_tree_names = None

    sim_seed = parameter_dict.get('sim_seed', None)

    if parameter_dict.get('parameters_path', False):
        parameters_path = parameter_dict['parameters_path']
        df_param = pd.read_pickle(parameters_path)
        ls_tree_names = list(df_param['name'])
        ls_avg_array_len = list(df_param['avg array length'])
        if parameter_dict['deletion_model'] == 'independent':
            ls_gtr_loss_rate = list(df_param['loss_rate_0'])
            ls_deletion_model = ['independent'] * len(ls_gtr_loss_rate)
            ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_len,
                                                                      ls_gtr_loss_rate,
                                                                      'independent'
                                                                      )
            ls_deletion_model_parameter = [None] * len(ls_gtr_gain_rate)
        elif parameter_dict['deletion_model'] == 'block':
            ls_gtr_loss_rate = list(df_param['loss_rate_1'])
            ls_deletion_model_parameter = list(df_param['alpha_1'])
            ls_deletion_model = ['block'] * len(ls_gtr_loss_rate)
            ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_len,
                                                                      ls_gtr_loss_rate,
                                                                      'block',
                                                                      ls_deletion_model_parameter,
                                                                      )
        else:
            raise NotImplementedError('Model not implemented for importing parameters!')
        ls_sim_seed = [None] * len(ls_gtr_loss_rate) if sim_seed is None else list(
            range(sim_seed, sim_seed + len(ls_gtr_loss_rate)))
    else:
        deletion_model = parameter_dict['deletion_model']
        deletion_model_parameter = parameter_dict.get('deletion_model_parameter', None)
        if parameter_dict['lr_based_on_gr_avg_array_length']:
            gtr_gain_rate = parameter_dict['gtr_gain_rate']
            gtr_loss_rate = misc.keep_len_constant(gtr_gain_rate, parameter_dict['avg_array_length'], deletion_model,
                                                   ls_alpha=deletion_model_parameter)
        else:
            gtr_gain_rate = parameter_dict['gtr_gain_rate']
            gtr_loss_rate = parameter_dict['gtr_loss_rate']
        nb_trees = parameter_dict['nb_trees']
        nb_leafs = parameter_dict['nb_leafs']
        ls_gtr_gain_rate = []
        ls_gtr_loss_rate = []
        ls_deletion_model = []
        ls_deletion_model_parameter = []
        ls_sim_seed = []
        for i, (tree_count, leaf_count) in enumerate(zip(nb_trees, nb_leafs)):
            if isinstance(gtr_gain_rate, list):
                ls_gtr_gain_rate += list(np.repeat(gtr_gain_rate[i], tree_count))
                ls_gtr_loss_rate += list(np.repeat(gtr_loss_rate[i], tree_count))
            else:
                ls_gtr_gain_rate += list(np.repeat(gtr_gain_rate, tree_count))
                ls_gtr_loss_rate += list(np.repeat(gtr_loss_rate, tree_count))
            ls_deletion_model += list(np.repeat(deletion_model[i], tree_count))
            ls_sim_seed += [None] * tree_count if sim_seed is None else list(range(sim_seed, sim_seed + tree_count))
            if deletion_model_parameter:
                ls_deletion_model_parameter += list(np.repeat(deletion_model_parameter[i], tree_count))
            else:
                ls_deletion_model_parameter += [None] * tree_count

    if parameter_dict.get('trees_path', False):
        if ls_tree_names is None:
            ls_tree_names = parameter_dict['tree_file_names']
        trees_format = parameter_dict.get('trees_format', 'nwk').lower()
        if trees_format in {'nwk', 'newick'}:
            ls_trees = import_data.load_real_data_tree(
                parameter_dict['trees_path'],
                [n + '_tree.nwk' for n in ls_tree_names],
            )
        elif trees_format in {'pkl', 'pickle'}:
            with open(parameter_dict['trees_path'], 'rb') as f:
                ls_trees = pickle.load(f)
        else:
            raise ValueError(f'Unsupported trees_format: {trees_format}')
        if len(ls_trees) != len(ls_tree_names):
            raise ValueError('The number of loaded trees must match tree_file_names')
        for n, t in zip(ls_tree_names, ls_trees):
            t.name = n
    else:
        ls_trees = []
        ls_tree_names = []
        tree_seed = parameter_dict.get('tree_seed', None)
        for tree_count, leaf_count in zip(nb_trees, nb_leafs):
            ls_trees += mult_rnd_coalescent_tree(tree_count, leaf_count, rnd_bl=parameter_dict['rnd_bl'], seed=tree_seed,
                                                 tree_name='t')
        ls_tree_names += ['_'.join([parameter_dict['model_name'], str(i)]) for i, t in enumerate(ls_trees)]

    return ls_tree_names, ls_trees, ls_gtr_gain_rate, ls_gtr_loss_rate, ls_deletion_model, \
        ls_deletion_model_parameter, ls_sim_seed


def pooling_for_parameter_estimation(ls_data, group_by=None, give_lh_fct=None):
    if group_by is None:
        raise ValueError('group_by is required for pooled parameter estimation')
    if len(group_by) != len(ls_data):
        raise ValueError('group_by and likelihood data must have the same length')
    grouped_indices = {}
    for index, label in enumerate(group_by):
        missing = pd.isna(label)
        if not isinstance(missing, (bool, np.bool_)):
            raise ValueError('group labels must be scalar values')
        if missing:
            raise ValueError('group labels must not contain missing values')
        try:
            grouped_indices.setdefault(label, []).append(index)
        except TypeError as error:
            raise ValueError('group labels must be hashable scalar values') from error

    values = np.empty(len(grouped_indices), dtype=object)
    values[:] = list(grouped_indices)
    ls_lh_0, ls_lh_1, ls_opt_result_0, ls_opt_result_1 = [], [], [], []
    for group_idx in grouped_indices.values():
        group_data = [ls_data[idx] for idx in group_idx]
        lh_0, lh_1, opt_result_0, opt_result_1 = model_tools.compute_lh_ratio_of_multiple_trees(group_data,
                                                                                                method='Nelder-Mead',
                                                                                                filter_by=False,
                                                                                                lambdifyed_lh_fct=give_lh_fct)
        ls_lh_0.append(lh_0)
        ls_lh_1.append(lh_1)
        ls_opt_result_0.append(opt_result_0)
        ls_opt_result_1.append(opt_result_1)
    return values, ls_lh_0, ls_lh_1, ls_opt_result_0, ls_opt_result_1
