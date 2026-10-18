"""Independent controls for spacer order and reusable result summaries."""

import copy
import random

import numpy as np
import pandas as pd
import pytest
from Bio.Phylo.Newick import Clade

from sp_model import model_tools
from sp_model.helpers.misc import GraphContainsCycleException, topological_sort
from sp_model.summary import compose_summary_dict, write_summary
from sp_model.reconstruction_tree import ReconstructionTree


def _reference_topological_sort(roots, adjacency):
    """Historical LIFO traversal, including its matrix mutation contract."""
    order = []
    while roots:
        node = roots.pop()
        order.append(node)
        for child in range(adjacency.shape[1]):
            if adjacency[node, child] == 1:
                adjacency[node, child] = 0
                if sum(adjacency[:, child]) == 0:
                    roots.append(child)
    if np.any(adjacency):
        raise GraphContainsCycleException
    return order


def _reachable(start, adjacency):
    found = set()
    pending = [start]
    while pending:
        node = pending.pop()
        if node not in found:
            found.add(node)
            pending.extend(np.flatnonzero(adjacency[node]))
    return found


@pytest.mark.parametrize('exclude_self', [False, True])
def test_order_graph_matches_historical_ties_and_independent_reachability(exclude_self):
    generator = random.Random(2357)
    for _ in range(100):
        names = [f'spacer.{index}' for index in range(generator.randrange(1, 20))]
        arrays = [[name] for name in names]
        arrays.extend(
            [first, second]
            for index, first in enumerate(names)
            for second in names[index + 1:]
            if generator.random() < 0.2
        )
        (descendants, ancestors), order, (adjacency, labels) = model_tools.determine_order_dict(
            arrays, plot_order=False, exclude_self=exclude_self,
        )
        roots = [node for node in range(len(names)) if not np.any(adjacency[:, node])]
        expected_order = _reference_topological_sort(roots, adjacency.copy())
        assert order == [labels[node] for node in expected_order]
        assert list(descendants) == names
        assert list(ancestors) == order
        for node, name in labels.items():
            expected = _reachable(node, adjacency)
            if exclude_self:
                expected.remove(node)
            assert descendants[name] == {labels[child] for child in expected}
            assert ancestors[name] == {parent for parent in names if name in descendants[parent]}


@pytest.mark.parametrize('edges', [[], [(0, 1), (0, 2), (1, 3), (2, 3)], [(0, 1), (1, 0)]])
def test_topological_sort_preserves_mutations_and_cycle_errors(edges):
    adjacency = np.zeros((4, 4))
    for first, second in edges:
        adjacency[first, second] = 1
    actual = adjacency.copy()
    expected = adjacency.copy()
    actual_roots = [node for node in range(4) if not np.any(adjacency[:, node])]
    expected_roots = actual_roots.copy()
    labels = {node: str(node) for node in range(4)}
    if edges == [(0, 1), (1, 0)]:
        with pytest.raises(GraphContainsCycleException):
            _reference_topological_sort(expected_roots, expected)
        with pytest.raises(GraphContainsCycleException):
            topological_sort(actual_roots, actual, labels)
    else:
        assert topological_sort(actual_roots, actual, labels) == [
            labels[node] for node in _reference_topological_sort(expected_roots, expected)
        ]
    assert np.array_equal(actual, expected)
    assert actual_roots == expected_roots


def test_order_graph_supports_chains_beyond_python_recursion_limit():
    names = [str(index) for index in range(1200)]
    (descendants, ancestors), order, _ = model_tools.determine_order_dict([names], plot_order=False)
    assert order == names
    assert descendants[names[0]] == set(names[1:])
    assert ancestors[names[-1]] == set(names[:-1])


def test_legacy_subgraph_helper_is_iterative_and_preserves_memoized_sets():
    size = 1200
    adjacency = np.zeros((size, size))
    adjacency[np.arange(size - 1), np.arange(1, size)] = 1
    cached = {node: set() for node in range(size)}
    cached[size - 1] = {size - 1}
    final_cached_set = cached[size - 1]
    assert model_tools.get_subgraph_sets(0, adjacency, cached) == set(range(size))
    assert cached[size - 1] is final_cached_set
    assert np.count_nonzero(adjacency) == size - 1


def test_legacy_subgraph_helper_rejects_cycles():
    adjacency = np.array([[0, 1], [1, 0]])
    with pytest.raises(GraphContainsCycleException):
        model_tools.get_subgraph_sets(0, adjacency, {0: set(), 1: set()})


def test_order_graph_rejects_self_and_multi_node_cycles():
    for arrays in ([['a', 'a']], [['a', 'b'], ['b', 'a']]):
        with pytest.raises(GraphContainsCycleException):
            model_tools.determine_order_dict(arrays, plot_order=False)


@pytest.mark.parametrize(('statistic_key', 'significance_key'), [
    ('test_statistic (-2*ln_lh_ratio)', 'significance level'),
    ('test statistic (-2*(ln_lh_ratio)', 'significance_level'),
])
def test_summary_preserves_native_and_legacy_likelihood_ratio_fields(tmp_path, statistic_key, significance_key):
    parameters = {'input_path': 'in.fa', 'output_path': 'out', 'determine_orientation': False,
                  'combine_non_unique_arrays': True}
    fields = {'ln_lh_idm': -5.0, 'ln_lh_bdm': -3.25, statistic_key: 3.5, significance_key: 0.05}
    results = pd.DataFrame([{**fields, 'nb of reconstructed deletions': 2}], index=['group'])
    original_results = results.copy(deep=True)
    summary = compose_summary_dict(results, parameters)
    assert summary['ind_groups']['group']['idm_vs_bdm_test'] == fields
    original_summary = copy.deepcopy(summary)
    for name in ('first.txt', 'second.txt'):
        output = tmp_path / name
        write_summary(summary, output)
        assert f'{statistic_key}: 3.5' in output.read_text()
        assert f'{significance_key}: 0.05' in output.read_text()
    assert summary == original_summary
    pd.testing.assert_frame_equal(results, original_results)
    assert compose_summary_dict(results, parameters) == summary


@pytest.mark.parametrize('combine_non_unique_arrays', [False, True])
def test_summary_can_be_composed_and_rendered_repeatedly_without_mutating_inputs(
    tmp_path, combine_non_unique_arrays,
):
    parameters = {
        'input_path': 'input.fa',
        'output_path': 'results',
        'determine_orientation': True,
        'combine_non_unique_arrays': combine_non_unique_arrays,
    }
    results = pd.DataFrame({
        'nb of reconstructed deletions': [0],
        'tree length': [1.5],
        'nb of leafs (before combining non-uniques)': [3],
        'nb of leafs (after combining non-uniques)': [2],
        'predicted orientation': ['ND'],
    }, index=['group'])
    original_parameters = parameters.copy()
    summary = compose_summary_dict(results, parameters)
    assert parameters == original_parameters
    assert compose_summary_dict(results, parameters) == summary
    original_summary = copy.deepcopy(summary)
    first = tmp_path / 'first.txt'
    second = tmp_path / 'second.txt'
    write_summary(summary, first)
    write_summary(summary, second)
    assert summary == original_summary
    assert first.read_bytes() == second.read_bytes()
    text = first.read_text()
    assert 'WARNING: No deletions were reconstructed' in text
    assert 'predicted orientation: ND' in text
    if combine_non_unique_arrays:
        assert 'nb of leafs (before combining non-uniques): 3' in text
        assert 'nb of leafs (after combining non-uniques): 2' in text
    else:
        assert 'nb of leafs: 3' in text
        assert 'combining non-uniques' not in text


@pytest.mark.parametrize(('event_attribute', 'event_index'), [
    ('rec_duplications_dict', 0),
    ('rec_rearrangements_dict', 1),
    ('rec_reacquisition_dict', 2),
    ('rec_indep_gain_dict', 3),
    ('rec_other_dup_events_dict', 4),
])
def test_independent_acquisition_count_excludes_repeated_duplicate_candidates(event_attribute, event_index):
    root = Clade(name='root', branch_length=0, clades=[
        Clade(name='first', branch_length=1), Clade(name='second', branch_length=1),
    ])
    tree = ReconstructionTree(
        {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 0.1}, None, root, True,
    )
    tree.duplicated_spacers = {'duplicate': [2, 3]}
    tree.rec_gain_dict = {'root': [1, 2, 3], 'first': [1, 2], 'second': [1, 3]}
    setattr(tree, event_attribute, {'first': [2]})
    original_gains = copy.deepcopy(tree.rec_gain_dict)
    expected = [0] * 5 + [2]
    expected[event_index] = 1
    counts = tree.get_rec_duplication_type_counts()
    assert counts == tuple(expected)
    assert all(type(count) is int for count in counts)
    assert tree.rec_gain_dict == original_gains
    assert tree.duplicated_spacers == {'duplicate': [2, 3]}
