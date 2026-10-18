"""Both hypotheses receive the same legacy spacer-label plotting choice."""

import pandas as pd
import pytest

import spacerplacer
from sp_model import run_experiments
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup


@pytest.mark.parametrize('trusted_pickle', [False, True])
@pytest.mark.parametrize('show_names', [None, False])
def test_spacer_label_choice_reaches_both_orientation_fits(tmp_path, monkeypatch, trusted_pickle, show_names):
    arrays = [CRISPRArray(name, None, None, None, None, None, None, None, ['s1', 's2'], None)
              for name in ('a', 'b')]
    group = CRISPRGroup(None, arrays, name='group')
    calls = []

    class FitsIntercepted(Exception):
        pass

    def fake_reconstruction(*_args, **kwargs):
        calls.append(kwargs)
        if len(calls) == 2:
            raise FitsIntercepted
        return pd.DataFrame([{'name': 'group'}]), pd.DataFrame(), ({}, None), {}, {}, {}, {}

    monkeypatch.setattr(run_experiments, 'run_reconstruction', fake_reconstruction)
    options = {} if show_names is None else {'show_spacer_name_row': show_names}
    output = tmp_path / 'output'
    if trusted_pickle:
        monkeypatch.setattr(run_experiments, 'load_arrays_from_pkl', lambda _path: {'group': group})
        with pytest.raises(FitsIntercepted):
            run_experiments.run_pickled_data(
                {}, 'trusted.pkl', save_path=output, tree_save_path=output,
                determine_orientation=True, **options,
            )
    else:
        monkeypatch.setattr(run_experiments, 'load_align_single_fasta', lambda *_args, **_kwargs: group)
        with pytest.raises(FitsIntercepted):
            run_experiments.run_multiple_groups(['group.fa'], output, **options)

    assert [call['show_spacer_name_row'] for call in calls] == [show_names is not False] * 2


def test_legacy_spacer_label_flag_hides_names_and_help_describes_it(capsys):
    assert spacerplacer.parse_args(['input.fa', 'output']).show_spacer_name_row is True
    assert spacerplacer.parse_args(['input.fa', 'output', '--show_spacer_name_row']).show_spacer_name_row is False
    with pytest.raises(SystemExit):
        spacerplacer.parse_args(['--help'])
    assert 'original names are hidden' in ' '.join(capsys.readouterr().out.split())
