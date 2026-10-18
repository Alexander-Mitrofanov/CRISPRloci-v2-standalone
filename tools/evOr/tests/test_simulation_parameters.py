import numpy as np
import pandas as pd
import pytest

from sp_model.coalescent_tree_sim import rnd_coalescent_tree
from sp_model.experiments_tools import expand_sim_parameters, expand_sim_parameters_based_on_import


def _tree_signature(tree):
    return tuple(
        (len(clade.clades), None if clade.branch_length is None else round(clade.branch_length, 12))
        for clade in tree.find_clades(order='preorder')
    )


def _imported_parameter_file(tmp_path):
    path = tmp_path / 'parameters.pkl'
    pd.DataFrame(
        {
            'avg array length': [8.0, 10.0],
            'loss_rate_0': [0.1, 0.2],
            'nb leafs (after combining)': [5, 5],
        },
        index=['a', 'b'],
    ).to_pickle(path)
    return path


def test_seeded_coalescent_generation_does_not_mutate_global_numpy_rng():
    np.random.seed(1234)
    expected_next = np.random.random()
    np.random.seed(1234)
    first = rnd_coalescent_tree(6, rnd_bl=True, seed=7)
    actual_next = np.random.random()
    second = rnd_coalescent_tree(6, rnd_bl=True, seed=7)
    assert actual_next == expected_next
    assert _tree_signature(first) == _tree_signature(second)


@pytest.mark.parametrize('per_repeat', [False, True])
def test_imported_simulation_tree_repeat_policy(tmp_path, per_repeat):
    parameters = {
        'import_path': _imported_parameter_file(tmp_path),
        'deletion_model': 'independent',
        'repeat_runs': 2,
        'sim_seed': 10,
        'tree_seed': 20,
        'tree_sim_per_repeat_run': per_repeat,
    }
    names, trees, *_, repeat_runs = expand_sim_parameters_based_on_import(parameters)
    assert names == ['a-0', 'b-0', 'a-1', 'b-1']
    assert repeat_runs == [0, 0, 1, 1]
    assert len({id(tree) for tree in trees}) == 4
    same_group_signatures_equal = _tree_signature(trees[0]) == _tree_signature(trees[2])
    assert same_group_signatures_equal is (not per_repeat)


def test_named_newick_tree_loading_uses_mapping_lookup_and_explicit_format(tmp_path):
    (tmp_path / 'custom_tree.nwk').write_text('(a:1,b:1)root;\n', encoding='utf-8')
    parameters = {
        'nb_trees': [1],
        'nb_leafs': [2],
        'deletion_model': ['independent'],
        'deletion_model_parameter': None,
        'lr_based_on_gr_avg_array_length': False,
        'gtr_gain_rate': [1.0],
        'gtr_loss_rate': [0.1],
        'sim_seed': 5,
        'trees_path': tmp_path,
        'tree_file_names': ['custom'],
        'trees_format': 'newick',
    }
    names, trees, *_ = expand_sim_parameters(parameters)
    assert names == ['custom']
    assert [leaf.name for leaf in trees[0].get_terminals()] == ['a', 'b']


def test_named_tree_loading_rejects_unknown_format(tmp_path):
    parameters = {
        'nb_trees': [1],
        'nb_leafs': [2],
        'deletion_model': ['independent'],
        'deletion_model_parameter': None,
        'lr_based_on_gr_avg_array_length': False,
        'gtr_gain_rate': [1.0],
        'gtr_loss_rate': [0.1],
        'sim_seed': 5,
        'trees_path': tmp_path,
        'tree_file_names': ['custom'],
        'trees_format': 'xml',
    }
    with pytest.raises(ValueError, match='Unsupported trees_format'):
        expand_sim_parameters(parameters)
