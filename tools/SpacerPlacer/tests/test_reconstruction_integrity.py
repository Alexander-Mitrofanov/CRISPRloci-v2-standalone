import json
import logging
import pickle

import pandas as pd
import pytest
from Bio.Phylo.Newick import Clade, Tree

import spacerplacer
from sp_model import mafft_integration, run_experiments
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup
from sp_model.experiments_tools import name_inner_nodes_if_unnamed, write_aligned_groups_json
from sp_model.helpers.import_data import load_single_tree_from_string
from sp_model.run_reconstruction import run_reconstruction


def _group(name='group', arrays=None):
    if arrays is None:
        arrays = {'a': ['1', '2', '3'], 'b': ['1', '3'], 'c': ['1', '2']}
    return CRISPRGroup(None, [
        CRISPRArray(array_name, None, None, None, None, None, None, None, spacers, None)
        for array_name, spacers in arrays.items()
    ], name=name)


def _capture_reconstruction_calls(monkeypatch):
    calls = []

    def reconstruct(parameters, groups, **kwargs):
        calls.append((groups, kwargs))
        return pd.DataFrame({'name': ['group']}), pd.DataFrame(), ({}, None), {}, {'group': '(a:1,b:1);'}, {}, {}

    def compare(forward, reverse, *args, **kwargs):
        return forward, forward.copy(), forward.copy(), {'group': '(a:1,b:1);'}

    monkeypatch.setattr(run_experiments, 'run_reconstruction', reconstruct)
    monkeypatch.setattr(run_experiments.orientation_tools, 'compare_likelihoods_for_orientation', compare)
    return calls


@pytest.mark.parametrize('determine_orientation', [False, True])
@pytest.mark.parametrize('save_trees', [False, True])
def test_trusted_pickle_preserves_optional_tree_save_paths(tmp_path, monkeypatch, determine_orientation, save_trees):
    source = tmp_path / 'groups.pkl'
    group = _group(arrays={'a': ['1', '2'], 'b': ['1']})
    with source.open('wb') as handle:
        pickle.dump({'group': group}, handle)
    output = tmp_path / 'result'
    output.mkdir()
    tree_root = tmp_path / 'trees' if save_trees else None
    if tree_root is not None:
        tree_root.mkdir()
    calls = _capture_reconstruction_calls(monkeypatch)

    run_experiments.run_pickled_data(
        {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 0.1}, str(source), save_path=output,
        tree_save_path=tree_root, determine_orientation=determine_orientation, plot_tree=False, plot_order=False,
    )

    orientations = ['0_forward', '0_reversed'] if determine_orientation else ['0_forward']
    assert [options['tree_save_path'] for _, options in calls] == [
        str(tree_root / orientation) if save_trees else None for orientation in orientations
    ]
    assert calls[0][0]['group'].get_arrays_as_lists() == group.get_arrays_as_lists()
    if determine_orientation:
        assert calls[1][0]['group'].get_arrays_as_lists() == {'a': ['2', '1'], 'b': ['1']}
    if save_trees:
        assert json.loads((tree_root / 'final_dict_nwk_trees.json').read_text()) == {'group': '(a:1,b:1);'}


@pytest.mark.parametrize('show_spacer_name_row', [False, True])
def test_both_orientations_preserve_spacer_name_plot_option(tmp_path, monkeypatch, show_spacer_name_row):
    monkeypatch.setattr(run_experiments, 'load_align_single_fasta', lambda *args, **kwargs: _group())
    calls = _capture_reconstruction_calls(monkeypatch)
    run_experiments.run_multiple_groups(
        [tmp_path / 'group.fa'], tmp_path / 'result', plot_tree=True, plot_order=False,
        show_spacer_name_row=show_spacer_name_row,
    )
    assert [options.get('show_spacer_name_row', True) for _, options in calls] == [show_spacer_name_row] * 2


@pytest.mark.parametrize('option', [
    '--insertion_rate', '--deletion_rate', '--tree_insertion_rate', '--tree_deletion_rate',
    '--tree_alpha', '--extend_branches', '--orientation_decision_boundary', '--significance_level',
])
@pytest.mark.parametrize('value', ['nan', 'inf', '-inf'])
def test_cli_rejects_non_finite_scientific_options(option, value):
    with pytest.raises(SystemExit):
        spacerplacer.parse_args(['input.fa', 'output', f'{option}={value}'])


def test_dotted_groups_reach_reconstruction_separately(tmp_path, monkeypatch):
    paths = [tmp_path / 'isolate.one.fa', tmp_path / 'isolate.two.fa']
    aligned = []

    def align(data_path, work_path, *, group_name, **kwargs):
        aligned.append((group_name, work_path))
        return _group(group_name)

    def reconstruct(parameters, groups, **kwargs):
        assert set(groups) == {'isolate.one', 'isolate.two'}
        raise RuntimeError('reconstruction reached both groups')

    monkeypatch.setattr(run_experiments, 'load_align_single_fasta', align)
    monkeypatch.setattr(run_experiments, 'run_reconstruction', reconstruct)
    with pytest.raises(RuntimeError, match='reached both groups'):
        run_experiments.run_multiple_groups(paths, tmp_path / 'result', alignment_workers=2)
    assert len({path for _, path in aligned}) == 2


def test_duplicate_group_stems_fail_before_any_alignment(tmp_path, monkeypatch):
    def align(*args, **kwargs):
        raise AssertionError('ambiguous inputs must fail before alignment')

    monkeypatch.setattr(run_experiments, 'load_align_single_fasta', align)
    with pytest.raises(ValueError, match='Duplicate input group name'):
        run_experiments.run_multiple_groups(
            [tmp_path / 'group.fa', tmp_path / 'group.fasta'], tmp_path / 'result', alignment_workers=2,
        )


@pytest.mark.parametrize('mutation', ['group', 'array_name', 'spacer', 'order', 'extra_array'])
def test_alignment_cache_must_match_observed_input_arrays(tmp_path, mutation):
    source = tmp_path / 'group.fa'
    source.write_text('>a, metadata\n1, 2\n>b\n1\n', encoding='utf-8')
    name = 'other' if mutation == 'group' else 'group'
    arrays = {'a': ['1', '2'], 'b': ['1', '-']}
    if mutation == 'array_name':
        arrays['other'] = arrays.pop('a')
    elif mutation == 'spacer':
        arrays['a'][1] = '3'
    elif mutation == 'order':
        arrays['a'].reverse()
    elif mutation == 'extra_array':
        arrays['c'] = ['1', '-']
    cache = tmp_path / 'cache.json'
    write_aligned_groups_json({name: _group(name, arrays)}, cache)
    with pytest.raises(ValueError, match='cache'):
        run_experiments.run_multiple_groups([source], tmp_path / 'result', aligned_groups_path=cache)


def test_valid_alignment_cache_and_extended_tree_export_match_the_fit(tmp_path):
    source = tmp_path / 'group.FA'
    source.write_text('>a\n1,2,3\n>b\n1,3\n>c\n1,2\n', encoding='utf-8')
    group = _group(arrays={'a': ['1', '2', '3'], 'b': ['1', '-', '3'], 'c': ['1', '2', '-']})
    cache = tmp_path / 'cache.json'
    write_aligned_groups_json({'group': group}, cache)
    trees = tmp_path / 'trees.JSON'
    trees.write_text(json.dumps({'group': '(a:1,b:1,c:1);'}), encoding='utf-8')
    output = tmp_path / 'result'
    arguments = spacerplacer.parse_args([str(source), str(output)])
    paths, _ = spacerplacer.check_and_parse_input_data(
        arguments.input_type, str(source), str(output), None, logging.getLogger('integrity'),
    )
    assert paths == [str(source)]
    results = run_experiments.run_multiple_groups(
        paths, output, aligned_groups_path=cache, tree_path=trees,
        determine_orientation=False, plot_tree=False, plot_order=False,
        rho_bias_correction=False, alpha_bias_correction=False,
        extend_branches=0.5, save_reconstructed_events=True,
    )
    exported = json.loads((output / 'additional_data' / 'dict_nwk_trees.json').read_text())
    tree = load_single_tree_from_string(exported['group'])
    assert tree.total_branch_length() == pytest.approx(4.5)
    assert results.loc['group', 'tree length'] == pytest.approx(tree.total_branch_length())
    event_tree = load_single_tree_from_string(
        (output / '0_forward' / 'reconstructed_events' / 'group_rec' / 'group_rec_tree.nwk').read_text(),
    )
    assert event_tree.total_branch_length() == pytest.approx(tree.total_branch_length())
    assert {node.name for node in tree.find_clades()} == {node.name for node in event_tree.find_clades()}


@pytest.mark.parametrize('document', [{}, [], {'group': None}, {'group': ''}])
def test_explicit_tree_mapping_cannot_silently_fall_back_to_inference(tmp_path, document):
    tree_path = tmp_path / 'trees.json'
    tree_path.write_text(json.dumps(document), encoding='utf-8')
    with pytest.raises(ValueError, match='non-empty mapping'):
        run_reconstruction(
            {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 0.1}, {'group': _group()},
            save_path=tmp_path / 'result', tree_path=tree_path, plot_tree=False, plot_order=False,
        )


@pytest.mark.parametrize('newick', [
    '(a:1,b:1);\n(a:2,b:2);', '(a:-1,b:1);', '(a,b:1);', '(a:1,a:1);', '(:1,b:1);',
    '(a:nan,b:1);', '(a:inf,b:1);',
])
def test_single_tree_rejects_ambiguous_or_invalid_phylogenies(newick):
    with pytest.raises(ValueError):
        load_single_tree_from_string(newick)


def test_supplied_tree_leaf_identity_is_checked_before_reconstruction(tmp_path):
    tree_path = tmp_path / 'trees.json'
    tree_path.write_text(json.dumps({'group': '(a:1,b:1,unrelated:1);'}), encoding='utf-8')
    with pytest.raises(ValueError, match='leaves disagree'):
        run_reconstruction(
            {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 0.1}, {'group': _group()},
            save_path=tmp_path / 'result', tree_path=tree_path, plot_tree=False, plot_order=False,
        )


def test_internal_node_names_cannot_collide_with_existing_leaf_or_internal_names():
    tree = Tree(root=Clade(clades=[
        Clade(name='in_0'), Clade(name='in_1', clades=[Clade(name='leaf')]),
    ]))
    name_inner_nodes_if_unnamed(tree)
    assert tree.root.name == 'in_2'
    names = [node.name for node in tree.find_clades()]
    assert len(names) == len(set(names))
    tree.root.name = 'leaf'
    with pytest.raises(ValueError, match='unique'):
        name_inner_nodes_if_unnamed(tree)


@pytest.mark.parametrize(('content', 'message'), [
    ('>a\n01 02\n', 'array names disagree'),
    ('>a\n01 02\n>a\n01 --\n', 'array names disagree'),
    ('>a\n02 01\n>b\n01 --\n', 'changed the observed'),
    ('>a\n01 03\n>b\n01 --\n', 'unknown spacer codes'),
])
def test_mafft_must_preserve_array_identity_and_each_observed_spacer(tmp_path, monkeypatch, content, message):
    def fake_mafft(file_path, output_path, *args, **kwargs):
        with open(output_path + '_output.txt', 'w', encoding='utf-8') as output_file:
            output_file.write(content)

    monkeypatch.setattr(mafft_integration, 'run_mafft', fake_mafft)
    group = _group(arrays={'a': ['1', '2'], 'b': ['1']})
    with pytest.raises(ValueError, match=message):
        mafft_integration.align_crispr_groups(tmp_path, {'group': group})


def test_mafft_can_reorder_records_while_preserving_their_spacer_identity(tmp_path, monkeypatch):
    def fake_mafft(file_path, output_path, *args, **kwargs):
        with open(output_path + '_output.txt', 'w', encoding='utf-8') as output_file:
            output_file.write('>b\n01 --\n>a\n01 02\n')

    monkeypatch.setattr(mafft_integration, 'run_mafft', fake_mafft)
    group = _group(arrays={'a': ['1', '2'], 'b': ['1']})
    mafft_integration.align_crispr_groups(tmp_path, {'group': group})
    assert group.get_arrays_as_lists() == {'a': ['1', '2'], 'b': ['1', '-']}
