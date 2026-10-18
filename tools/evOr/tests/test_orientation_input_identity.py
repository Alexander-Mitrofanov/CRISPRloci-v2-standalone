"""Scientific input identity and fitted-tree provenance regressions."""

import json
import logging

import pytest

import crispr_evor
import spacerplacer
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup
from sp_model.experiments_tools import (
    load_aligned_groups_json,
    name_inner_nodes_if_unnamed,
    read_spacer_arrays,
    validate_aligned_groups_inputs,
    write_aligned_groups_json,
)
from sp_model import run_experiments
from sp_model.helpers.import_data import load_single_tree_from_string
from sp_model.run_experiments import run_multiple_groups
from sp_model.run_reconstruction import run_reconstruction


def _group(name):
    arrays = [
        CRISPRArray(array_name, None, None, None, None, None, None, None, spacers, None)
        for array_name, spacers in (
            ('a', ['s1', 's2', 's3']),
            ('b', ['-', 's2', 's3']),
            ('c', ['s1', '-', 's3']),
        )
    ]
    return CRISPRGroup(None, arrays, name=name)


def _input(tmp_path, name='group'):
    path = tmp_path / f'{name}.fa'
    path.write_text('>a, metadata\ns1,s2,s3\n>b\ns2,s3\n>c\ns1,s3\n', encoding='utf-8')
    return path


@pytest.mark.parametrize('option', ['--insertion_rate', '--extend_branches', '--orientation_decision_boundary'])
@pytest.mark.parametrize('value', ['nan', 'inf', '-inf'])
def test_cli_rejects_nonfinite_scientific_values(option, value):
    with pytest.raises(SystemExit):
        spacerplacer.parse_args(['input.fa', 'output', f'{option}={value}'])


def test_alignment_cache_requires_actual_input_group_array_and_spacer_identity(tmp_path):
    path = _input(tmp_path)
    group = _group('group')
    groups = {'group': group}
    validate_aligned_groups_inputs(groups, {'group': path})
    with pytest.raises(ValueError, match='exactly the input groups'):
        validate_aligned_groups_inputs(groups, {'other': path})
    group.crispr_dict['b'].spacer_array.reverse()
    with pytest.raises(ValueError, match='spacer order or identities'):
        validate_aligned_groups_inputs(groups, {'group': path})
    group.crispr_dict['b'].spacer_array.reverse()
    group.crispr_dict['renamed'] = group.crispr_dict.pop('b')
    with pytest.raises(ValueError, match='array identities'):
        validate_aligned_groups_inputs(groups, {'group': path})


def test_workflow_rejects_unrelated_cache_and_duplicate_group_stems(tmp_path):
    path = _input(tmp_path)
    cache = tmp_path / 'cache.json'
    write_aligned_groups_json({'other': _group('other')}, cache)
    with pytest.raises(ValueError, match='exactly the input groups'):
        run_multiple_groups([path], tmp_path / 'wrong-cache', aligned_groups_path=cache)
    with pytest.raises(ValueError, match='Duplicate input group name'):
        run_multiple_groups([path, tmp_path / 'group.fna'], tmp_path / 'duplicate')


@pytest.mark.parametrize('spacers', [[], ['-', '-'], ['bad id', 's3']])
def test_cache_rejects_empty_all_gap_and_whitespace_arrays(tmp_path, spacers):
    path = tmp_path / 'cache.json'
    write_aligned_groups_json({'group': _group('group')}, path)
    document = json.loads(path.read_text(encoding='utf-8'))
    document['groups'][0]['arrays'][0]['spacers'] = spacers
    path.write_text(json.dumps(document), encoding='utf-8')
    with pytest.raises(ValueError, match='Invalid spacer alignment'):
        load_aligned_groups_json(path)


def test_input_does_not_treat_a_gap_as_an_observed_spacer(tmp_path):
    path = _input(tmp_path)
    path.write_text('>a\ns1,-,s3\n>b\ns2,s3\n', encoding='utf-8')
    with pytest.raises(ValueError, match='reserved alignment gap'):
        read_spacer_arrays(path)


@pytest.mark.parametrize('trees', [{}, [], {'group': '(a:-0.1,b:0.1,c:0.1);'}, {'group': '(a,b:0.1,c:0.1);'}])
def test_provided_trees_fail_before_fitting_when_empty_or_invalid(tmp_path, trees):
    path = tmp_path / 'trees.json'
    path.write_text(json.dumps(trees), encoding='utf-8')
    with pytest.raises(ValueError, match='Provided trees|branch lengths'):
        run_reconstruction(
            {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 0.1},
            {'group': _group('group')},
            save_path=tmp_path / 'result',
            tree_path=path,
            plot_tree=False,
            plot_order=False,
        )


@pytest.mark.parametrize('use_cache', [False, True])
def test_both_orientation_fits_keep_dotted_groups_and_export_actual_extended_trees(tmp_path, monkeypatch, use_cache):
    names = ['family.a', 'family.b']
    paths = [_input(tmp_path, name) for name in names]
    cache = tmp_path / 'cache.json'
    write_aligned_groups_json({name: _group(name) for name in names}, cache)
    if not use_cache:
        # Supply a known alignment while exercising the input naming path;
        # reconstruction, reverse comparison, and all exports remain real.
        monkeypatch.setattr(run_experiments, 'load_align_single_fasta',
                            lambda _data_path, _work_path, *, group_name, **_kwargs: _group(group_name))
    trees = tmp_path / 'trees.json'
    trees.write_text(json.dumps({name: '(a:0.1,b:0.1,c:0.1);' for name in names}), encoding='utf-8')
    output = tmp_path / 'result'
    results = run_multiple_groups(
        paths,
        output,
        tree_path=trees,
        aligned_groups_path=cache if use_cache else None,
        extend_branches=0.2,
        determine_orientation=True,
        plot_tree=False,
        plot_order=False,
        alpha_bias_correction=False,
        rho_bias_correction=False,
        save_reconstructed_events=True,
        logger=logging.getLogger('orientation-input-identity-test'),
    )
    assert set(results.index) == set(names)
    metadata_path = crispr_evor.write_orientation_json(
        output, 5, tree_path=trees, alpha_bias_correction=False, rho_bias_correction=False,
    )
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    tree_maps = json.loads((output / 'additional_data' / 'orientation_trees.json').read_text(encoding='utf-8'))
    for entry in metadata['results']:
        name = entry['group']
        for policy, folder in [('forward', '0_forward'), ('reverse', '0_reversed')]:
            tree_string = tree_maps[policy][name]
            tree = load_single_tree_from_string(tree_string)
            assert all(leaf.branch_length == pytest.approx(0.3) for leaf in tree.get_terminals())
            assert all(node.name for node in tree.get_nonterminals())
            event_tree = output / folder / 'reconstructed_events' / f'{name}_rec' / f'{name}_rec_tree.nwk'
            assert event_tree.read_text(encoding='utf-8') == tree_string
        selected = 'reverse' if entry['prediction'] == 'Reverse' else 'forward'
        assert tree_maps['selected'][name] == tree_maps[selected][name]


def test_generated_internal_names_do_not_collide_with_observed_leaf_names():
    tree = load_single_tree_from_string('(in_0:0.1,b:0.1);')
    name_inner_nodes_if_unnamed(tree)
    assert tree.root.name == 'in_1'
    assert len({node.name for node in tree.find_clades()}) == 3
    duplicate = load_single_tree_from_string('(a:0.1,b:0.1)a;')
    with pytest.raises(ValueError, match='node names must be unique'):
        name_inner_nodes_if_unnamed(duplicate)
