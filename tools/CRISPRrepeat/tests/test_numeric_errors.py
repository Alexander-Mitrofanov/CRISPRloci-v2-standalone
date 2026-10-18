"""Unavailable or invalid native values must never be interpreted as zero."""

import pytest

from crisprrepeat.folding import Folder, Settings


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -0.01, 1.2])
def test_global_matrix_rejects_invalid_values_before_filtering(monkeypatch, invalid):
    class NativeResult:
        def params_subst(self, parameters):
            self.params = parameters

        def exp_params_subst(self, parameters):
            self.exp_params = parameters

        def mfe(self):
            return "(...)", 0.0

        def exp_params_rescale(self, mfe):
            pass

        def pf(self):
            return ".....", 0.0

        def bpp(self):
            matrix = [[0.0] * 6 for _ in range(6)]
            matrix[1][5] = invalid
            return matrix

    folder = Folder()
    monkeypatch.setattr(folder.RNA, "fold_compound", lambda *args: NativeResult())
    with pytest.raises(ValueError, match="Invalid probability"):
        folder.isolated("GAAAC")
    assert not folder.cache


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -0.01, 1.2, None])
def test_local_callback_rejects_invalid_and_unavailable_values(monkeypatch, invalid):
    folder = Folder()

    class NativeResult:
        def params_subst(self, parameters):
            self.params = parameters

        def exp_params_subst(self, parameters):
            self.exp_params = parameters

        def probs_window(self, length, options, callback, data):
            for i in range(1, 6):
                row = [0.0] * 6
                if i == 1:
                    row[5] = invalid
                callback(row, 5, i, 5, folder.RNA.PROBS_WINDOW_BPP, data)
                callback([None, 1.0], 1, i, 1, folder.RNA.PROBS_WINDOW_UP, data)
            return 1

    monkeypatch.setattr(folder.RNA, "fold_compound", lambda *args: NativeResult())
    with pytest.raises(ValueError, match="Incomplete local probabilities"):
        folder.context("GAAAC")


@pytest.mark.parametrize("kwargs", [
    {"temperature": True}, {"window": 150.5}, {"span": 80.0}, {"no_lonely_pairs": 1},
])
def test_settings_reject_wrong_types(kwargs):
    with pytest.raises(ValueError):
        Settings(**kwargs)
