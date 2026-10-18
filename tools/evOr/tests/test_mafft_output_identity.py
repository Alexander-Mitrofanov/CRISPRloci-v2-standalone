"""Native alignment output must retain each submitted array and spacer in order."""

from pathlib import Path

import pytest

from sp_model import mafft_integration
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup


def _group(arrays):
    return CRISPRGroup(None, [
        CRISPRArray(name, None, None, None, None, None, None, None, spacers, None)
        for name, spacers in arrays.items()
    ], name='group')


@pytest.mark.parametrize('output, message', [
    ('>a\n01 03\n>b\n02 03\n', 'spacer order or identities'),
    ('>a\n02 01 03\n>b\n-- 02 03\n', 'spacer order or identities'),
    ('>a\n01 04 03\n>b\n-- 02 03\n', 'Unknown MAFFT spacer token'),
    ('>a\n01 02 03\n>a\n01 02 03\n', 'array identities'),
    ('>a\n01 02 03\n', 'array identities'),
    ('>a\n01 02 03\n>other\n-- 02 03\n', 'array identities'),
])
def test_direct_alignment_rejects_corrupt_output_before_export_or_mutation(tmp_path, monkeypatch, output, message):
    original = {'a': ['s1', 's2', 's3'], 'b': ['s2', 's3']}
    group = _group(original)

    def fake_mafft(_input_path, output_path, *_args, **_kwargs):
        Path(output_path + '_output.txt').write_text(output, encoding='utf-8')

    monkeypatch.setattr(mafft_integration, 'run_mafft', fake_mafft)
    with pytest.raises(ValueError, match=message):
        mafft_integration.align_crispr_groups(tmp_path, {'group': group})

    assert group.get_arrays_as_lists() == original
    assert not (tmp_path / 'group_finished.txt').exists()


def test_direct_alignment_accepts_reordered_rows_and_gaps_in_original_group_order(tmp_path, monkeypatch):
    group = _group({'a': ['s1', 's2', 's3'], 'b': ['s2', 's3']})

    def fake_mafft(_input_path, output_path, *_args, **_kwargs):
        Path(output_path + '_output.txt').write_text('>b\n-- -- 02 03\n>a\n01 -- 02 03\n', encoding='utf-8')

    monkeypatch.setattr(mafft_integration, 'run_mafft', fake_mafft)
    mafft_integration.align_crispr_groups(tmp_path, {'group': group})

    assert list(group.crispr_dict) == ['a', 'b']
    assert group.get_arrays_as_lists() == {'a': ['s1', '-', 's2', 's3'], 'b': ['-', '-', 's2', 's3']}
    assert (tmp_path / 'group_finished.txt').is_file()


@pytest.mark.parametrize('corruption', [None, 'lost', 'reordered', 'unknown', 'duplicate-header'])
def test_pseudo_alignment_checks_decoded_identity_before_export_or_mutation(tmp_path, monkeypatch, corruption):
    # More than 248 identities selects the real pseudo-DNA encoder and decoder.
    spacers = [f's{index}' for index in range(249)]
    original = {'a': spacers, 'b': spacers[1:]}
    group = _group(original)

    def fake_mafft(input_path, output_path, *_args, **_kwargs):
        lines = Path(input_path + '.txt').read_text(encoding='utf-8').splitlines()
        # Native text-to-hex output drops the prefix written in the input.
        first = [token.removeprefix('0x') for token in lines[1].split()]
        second = [token.removeprefix('0x') for token in lines[3].split()]
        spacer_width = 27  # Twenty pseudo bases plus the seven-base repeat.
        second = ['--'] * spacer_width + second
        if corruption == 'lost':
            first[7:7 + spacer_width] = ['--'] * spacer_width
        elif corruption == 'reordered':
            first[7:34], first[34:61] = first[34:61], first[7:34]
        elif corruption == 'unknown':
            first[7] = '00'
        name = 'a' if corruption == 'duplicate-header' else 'b'
        Path(output_path + '_output.txt').write_text(
            f'>{name}\n{" ".join(second)}\n>a\n{" ".join(first)}\n', encoding='utf-8',
        )

    monkeypatch.setattr(mafft_integration, 'run_mafft', fake_mafft)
    if corruption is None:
        mafft_integration.align_crispr_groups(tmp_path, {'group': group}, seed=17)
        assert list(group.crispr_dict) == ['a', 'b']
        assert group.get_arrays_as_lists() == {'a': spacers, 'b': ['-', *spacers[1:]]}
        assert (tmp_path / 'group_finished.txt').is_file()
    else:
        with pytest.raises(ValueError, match='array identities|spacer order or identities|equally aligned'):
            mafft_integration.align_crispr_groups(tmp_path, {'group': group}, seed=17)
        assert group.get_arrays_as_lists() == original
        assert not (tmp_path / 'group_finished.txt').exists()
