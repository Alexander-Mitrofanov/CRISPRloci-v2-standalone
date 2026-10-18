"""Turner2004 evidence must survive later global parameter-table replacement."""

import pytest

from crisprrepeat.folding import Folder, Settings

SETTINGS = [Settings(), Settings(temperature=20, window=50, span=20, no_lonely_pairs=True),
            Settings(temperature=65, window=30, span=10)]
REPEATS = ["GCGCGAAAACGCGC", "GUUCAACACCCUCUUUUCCCCGUCAGGGGACUGAAAC"]
CONTEXTS = [sequence + "AAAA" + sequence for sequence in REPEATS] + [REPEATS[1] * 5]


@pytest.mark.parametrize("settings", SETTINGS)
@pytest.mark.parametrize("method,sequences", [("isolated", REPEATS), ("context", CONTEXTS)])
def test_folds_preserve_owned_turner_parameters_after_global_replacement(settings, method, sequences):
    baseline = Folder(settings)
    expected = [getattr(baseline, method)(sequence) for sequence in sequences]
    folder = Folder(settings)
    try:
        for sequence, result in zip(sequences, expected):
            folder.RNA.params_load_DNA_Mathews2004()
            ambient = folder.RNA.param(folder.model()).param_file
            assert getattr(folder, method)(sequence) == result
            assert folder.RNA.param(folder.model()).param_file == ambient
        assert folder.provenance == baseline.provenance
    finally:
        folder.RNA.params_load_RNA_Turner2004()


def test_different_settings_remain_independent_with_interleaved_native_folds():
    baselines = [Folder(settings) for settings in SETTINGS]
    expected = [(baseline.isolated(REPEATS[1]), baseline.context(CONTEXTS[-1])) for baseline in baselines]
    folders = [Folder(settings) for settings in SETTINGS]
    try:
        for folder, (isolated, context) in zip(reversed(folders), reversed(expected)):
            folder.RNA.params_load_RNA_Turner1999()
            assert folder.context(CONTEXTS[-1]) == context
            folder.RNA.params_load_DNA_Mathews2004()
            assert folder.isolated(REPEATS[1]) == isolated
    finally:
        folders[0].RNA.params_load_RNA_Turner2004()


@pytest.mark.parametrize("method,sequence", [("isolated", REPEATS[1]), ("context", CONTEXTS[-1])])
def test_global_replacement_during_mfe_or_local_callbacks_preserves_results(monkeypatch, method, sequence):
    baseline = Folder(SETTINGS[1])
    expected = getattr(baseline, method)(sequence)
    folder = Folder(SETTINGS[1])
    native_compound = folder.RNA.fold_compound

    class NativeResult:
        def __init__(self, *args):
            self.compound = native_compound(*args)

        def __getattr__(self, name):
            return getattr(self.compound, name)

        def mfe(self):
            result = self.compound.mfe()
            folder.RNA.params_load_DNA_Mathews2004()
            return result

        def probs_window(self, length, options, callback, data):
            changed = False

            def replace_global_parameters(*args):
                nonlocal changed
                if not changed:
                    folder.RNA.params_load_DNA_Mathews2004()
                    changed = True
                callback(*args)

            return self.compound.probs_window(length, options, replace_global_parameters, data)

    monkeypatch.setattr(folder.RNA, "fold_compound", NativeResult)
    try:
        assert getattr(folder, method)(sequence) == expected
    finally:
        folder.RNA.params_load_RNA_Turner2004()
