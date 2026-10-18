"""Fail closed at native entry and exit points before publishing evidence."""

import pytest

from crisprrepeat.folding import Folder


def test_reported_loop_size_is_independent_of_vienna_global_defaults(monkeypatch):
    folder = Folder()
    monkeypatch.setattr(folder.RNA.cvar, "min_loop_size", 0)
    assert folder.model().min_loop_size == folder.provenance["minimum_loop_size"] == 3
    assert folder.isolated("GAAAC")["structure"] == "....."


@pytest.mark.parametrize("method,sequence", [
    ("isolated", ""), ("isolated", "A" * 201), ("context", ""), ("context", "A" * 5001),
])
def test_direct_folding_api_bounds_sequences_before_native_allocation(monkeypatch, method, sequence):
    folder = Folder()
    monkeypatch.setattr(folder.RNA, "fold_compound", lambda *a: pytest.fail("Invalid substrate reached ViennaRNA"))
    with pytest.raises(ValueError, match="nucleotides"):
        getattr(folder, method)(sequence)


@pytest.mark.parametrize("structure,mfe", [
    (".....", float("nan")), (".....", float("inf")), ("....", 0.0), ("(.)..", 0.0), ("(...)", 0.0),
])
def test_bad_mfe_is_rejected_before_partition_function_rescaling(monkeypatch, structure, mfe):
    class NativeResult:
        def params_subst(self, parameters):
            self.params = parameters

        def exp_params_subst(self, parameters):
            self.exp_params = parameters

        def mfe(self):
            return structure, mfe

        def exp_params_rescale(self, energy):
            pytest.fail("Invalid MFE result reached native rescaling")

    folder = Folder()
    monkeypatch.setattr(folder.RNA, "fold_compound", lambda *a: NativeResult())
    with pytest.raises(ValueError):
        folder.isolated("AAAAA")
    assert not folder.cache


@pytest.mark.parametrize("fault", ["shifted_rows", "duplicate_pair_row", "duplicate_accessibility", "span"])
def test_invalid_native_local_row_identity_and_span_cannot_complete(monkeypatch, fault):
    folder = Folder()
    sequence = "GCGCGAAAACGCGC"
    n = len(sequence)

    class NativeResult:
        def params_subst(self, parameters):
            self.params = parameters

        def exp_params_subst(self, parameters):
            self.exp_params = parameters

        def probs_window(self, length, options, callback, data):
            for i in range(1, n + 1):
                callback([None, 1.0], 1, i, 1, folder.RNA.PROBS_WINDOW_UP, data)
                row = [0.0] * (n + 1)
                if fault == "span" and i == 1:
                    row[n] = 0.5
                callback(row, n, i - 1 if fault == "shifted_rows" else i, n,
                         folder.RNA.PROBS_WINDOW_BPP, data)
                if fault == "duplicate_pair_row" and i == 1:
                    callback(row, n, i, n, folder.RNA.PROBS_WINDOW_BPP, data)
                if fault == "duplicate_accessibility" and i == 1:
                    callback([None, 1.0], 1, i, 1, folder.RNA.PROBS_WINDOW_UP, data)
            return 1

    monkeypatch.setattr(folder.RNA, "fold_compound", lambda *a: NativeResult())
    if fault == "span":
        from crisprrepeat.folding import Settings
        folder.settings = Settings(window=10, span=5)
    with pytest.raises(ValueError, match="Incomplete local probabilities"):
        folder.context(sequence)


@pytest.mark.parametrize("global_name,value", [
    ("max_bp_span", 5), ("window_size", 10), ("energy_set", 1), ("tetra_loop", 0), ("logML", 1),
    ("betaScale", 2.0), ("pf_smooth", 0), ("sfact", 2.0), ("backtrack", 0),
])
def test_native_folds_are_independent_of_remaining_vienna_global_defaults(monkeypatch, global_name, value):
    sequence = "GCGCGAAAACGCGC"
    baseline = Folder()
    isolated = baseline.isolated(sequence)
    context = baseline.context(sequence + "AAAA" + sequence)
    model = baseline.model()
    monkeypatch.setattr(baseline.RNA.cvar, global_name, value)
    folder = Folder()
    actual = folder.model()
    for field in ("max_bp_span", "window_size", "energy_set", "special_hp", "logML",
                  "betaScale", "pf_smooth", "sfact", "backtrack"):
        assert getattr(actual, field) == getattr(model, field)
    assert folder.provenance == baseline.provenance
    assert folder.isolated(sequence) == isolated
    assert folder.context(sequence + "AAAA" + sequence) == context
