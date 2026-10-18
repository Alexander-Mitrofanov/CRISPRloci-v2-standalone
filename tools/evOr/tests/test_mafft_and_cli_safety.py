
import numpy as np
import pytest

import spacerplacer
from sp_model import run_experiments
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup
from sp_model.experiments_tools import (
    load_aligned_groups_json,
    write_aligned_groups_json,
)
from sp_model.mafft_integration import (
    DEFAULT_MAFFT_OPTIONS,
    convert_mafft_output,
    convert_spacer_arrays_to_random_pseudo_dna,
)


def test_cli_numeric_validation_and_version_parser():
    args = spacerplacer.parse_args([
        'input.fa', 'output', '--significance_level', '0.5', '--alignment_workers', '3'
    ])
    assert args.significance_level == 0.5
    assert args.alignment_workers == 3
    with pytest.raises(SystemExit):
        spacerplacer.parse_args(['input.fa', 'output', '--insertion_rate', '0'])
    with pytest.raises(SystemExit):
        spacerplacer.parse_args(['input.fa', 'output', '--cluster_spacers_max_distance', '-1'])
    with pytest.raises(SystemExit):
        spacerplacer.parse_args(['input.fa', 'output', '--alignment_workers', '0'])
    with pytest.raises(ValueError, match='positive integer'):
        run_experiments.run_multiple_groups([], 'output', alignment_workers=0)


def test_path_validation_blocks_pickle_nested_and_stale_output(tmp_path):
    input_file = tmp_path / 'input.pkl'
    input_file.write_bytes(b'not actually loaded')
    output = tmp_path / 'output'
    args = spacerplacer.parse_args([str(input_file), str(output), '--input_type', 'pickled'])
    with pytest.raises(ValueError, match='disabled'):
        spacerplacer.validate_paths_and_security(args)

    input_directory = tmp_path / 'input'
    input_directory.mkdir()
    nested = input_directory / 'results'
    args = spacerplacer.parse_args([str(input_directory), str(nested)])
    with pytest.raises(ValueError, match='inside'):
        spacerplacer.validate_paths_and_security(args)

    clean_input = tmp_path / 'clean.fa'
    clean_input.write_text('>a\n1\n', encoding='utf-8')
    output.mkdir()
    (output / 'stale.txt').write_text('stale', encoding='utf-8')
    args = spacerplacer.parse_args([str(clean_input), str(output)])
    with pytest.raises(FileExistsError, match='not empty'):
        spacerplacer.validate_paths_and_security(args)


def test_mafft_output_parser_rejects_malformed_and_ragged_files(tmp_path):
    before_header = tmp_path / 'before.txt'
    before_header.write_text('01 02\n', encoding='utf-8')
    with pytest.raises(ValueError, match='before a header'):
        convert_mafft_output(before_header)

    ragged = tmp_path / 'ragged.txt'
    ragged.write_text('>a\n01 02\n>b\n01\n', encoding='utf-8')
    with pytest.raises(ValueError, match='not equally aligned'):
        convert_mafft_output(ragged)


def test_large_alphabet_pseudo_encoding_is_seeded_and_does_not_mutate_defaults():
    arrays = [[index] for index in range(249)]
    defaults_before = list(DEFAULT_MAFFT_OPTIONS)
    first = convert_spacer_arrays_to_random_pseudo_dna(arrays, rng=np.random.default_rng(42))
    second = convert_spacer_arrays_to_random_pseudo_dna(arrays, rng=np.random.default_rng(42))
    assert first == second
    assert DEFAULT_MAFFT_OPTIONS == defaults_before
    assert len(first[1]) == 249


def test_ode_loader_cache_can_be_cleared_without_loading_large_artifact():
    run_experiments.load_ode_likelihood_functions.cache_clear()
    assert run_experiments.load_ode_likelihood_functions.cache_info().currsize == 0


def test_aligned_group_json_cache_round_trips_without_pickle(tmp_path):
    arrays = [
        CRISPRArray(
            name, None, None, None, None, None, None, None, spacers, None,
            cas_genes=None, all_cas_types=None, all_cas_genes=None, species=None, species_fc=None,
        )
        for name, spacers in [('first', ['1', '-', '2']), ('second', ['1', '3', '-'])]
    ]
    groups = {'group': CRISPRGroup(None, arrays, name='group')}
    cache_path = tmp_path / 'aligned_groups.json'

    write_aligned_groups_json(groups, cache_path)
    restored = load_aligned_groups_json(cache_path)

    assert list(restored) == ['group']
    assert restored['group'].get_arrays_as_lists() == groups['group'].get_arrays_as_lists()

    invalid_path = tmp_path / 'invalid.json'
    invalid_path.write_text('{}', encoding='utf-8')
    with pytest.raises(ValueError, match='schema'):
        load_aligned_groups_json(invalid_path)
