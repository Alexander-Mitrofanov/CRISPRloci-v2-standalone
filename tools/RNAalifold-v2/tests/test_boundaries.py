"""Native reproducibility and source/alignment boundaries, independent of the server."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

from rna_family_tool import core
from rna_family_tool.cli import run

FIXTURE = Path(__file__).parent / "data/family.fasta"
SEQUENCES = [
    "GGGAAACCCGGGAAACCCGGGAAACCCGGGAAACCC",
    "GGGAAACCCGGGAAGCCCGGGAAACCCGGGAAACCC",
    "CCGAAACGGGGGAAACCCGGGAAACCCGGGAAACCC",
]


@pytest.fixture
def native_oracle():
    """Use fresh, unmodified ViennaRNA defaults as the pre-refinement oracle."""
    pytest.importorskip("RNA")
    script = '''
import json, RNA, sys
assert RNA.__version__ == "2.7.0"
RNA.params_load_RNA_Turner2004()
fc = RNA.fold_compound(json.loads(sys.stdin.read()), RNA.md(), RNA.OPTION_MFE | RNA.OPTION_PF)
structure, mfe = fc.mfe()
fc.exp_params_rescale(mfe)
_, ensemble = fc.pf()
print(json.dumps([structure, mfe, ensemble, fc.bpp()]))
'''
    result = subprocess.run([sys.executable, "-c", script], input=json.dumps(SEQUENCES),
                            text=True, capture_output=True, check=True)
    return json.loads(result.stdout)


def assert_native_result(result, oracle):
    structure, mfe, ensemble, matrix = oracle
    assert result["structure"] == structure == "(((...)))(((...((((((...))))))...)))"
    assert result["mfe_score_kcal_mol"] == mfe == -18.6200008392334
    assert result["ensemble_free_energy_kcal_mol"] == ensemble == -18.721271514892578
    mfe_pairs = {(1, 9), (2, 8), (3, 7), (10, 36), (11, 35), (12, 34),
                 (16, 30), (17, 29), (18, 28), (19, 27), (20, 26), (21, 25)}
    expected = [{"left": left, "right": right, "probability": matrix[left][right]}
                for left in range(1, 37) for right in range(left + 1, 37)
                if matrix[left][right] >= 0.01 or (left, right) in mfe_pairs]
    assert result["probabilities"] == expected


@pytest.mark.parametrize("name,value", [("dangles", 0), ("min_loop_size", 8), ("energy_set", 2),
                                       ("betaScale", 0.2), ("backtrack_type", "C"), ("max_bp_span", 10)])
def test_fold_ignores_ambient_model_defaults(monkeypatch, native_oracle, name, value):
    import RNA
    monkeypatch.setattr(RNA.cvar, name, value)
    assert_native_result(core.fold(SEQUENCES), native_oracle)
    assert getattr(RNA.cvar, name) == value


def test_fold_owns_thermodynamic_and_boltzmann_tables(monkeypatch, native_oracle):
    import RNA
    original = RNA.fold_compound

    def foreign_tables(*args):
        if not args[2] & RNA.OPTION_EVAL_ONLY:
            RNA.params_load_DNA_Mathews2004()
        return original(*args)

    monkeypatch.setattr(RNA, "fold_compound", foreign_tables)
    try:
        assert_native_result(core.fold(SEQUENCES), native_oracle)
    finally:
        RNA.params_load_RNA_Turner2004()


@pytest.mark.parametrize("options", [{"molecule": "protein"}, {"threads": 0}, {"threads": 9},
                                      {"threads": True}, {"timeout_seconds": 0},
                                      {"timeout_seconds": 14401}, {"timeout_seconds": 1.5}])
def test_library_run_checks_bounds_before_native_work(tmp_path, monkeypatch, options):
    def unexpected_native(*args, **kwargs):
        pytest.fail("Invalid options reached native execution")

    monkeypatch.setattr(subprocess, "run", unexpected_native)
    output = tmp_path / "invalid"
    with pytest.raises(ValueError):
        run(FIXTURE, output, **{"molecule": "DNA", **options})
    assert not output.exists()


@pytest.mark.parametrize("payload", [b">a\n" + b"A" * 30 + b"\n>b\n" + b"C" * 30 + b"\n",
                                     b">a\n" + b"T" * 30 + b"\n>b\n" + b"C" * 30 + b"\n"])
def test_parser_rejects_unknown_molecule(payload):
    with pytest.raises(ValueError, match="molecule"):
        core.fasta(payload, molecule="protein")


def test_alignment_preserves_gaps_source_order_and_membership():
    groups = [{"id": "family000001", "member_ids": ["a", "c"], "sequence": "AU"},
              {"id": "family000002", "member_ids": ["b"], "sequence": "CAU"}]
    result = core.alignment_rows(groups, {"family000001": "A-T", "family000002": "CAT"})
    assert result == [{**groups[0], "sequence": "A-U"}, {**groups[1], "sequence": "CAU"}]
    assert groups[0]["sequence"] == "AU"
    assert core.columns([row["sequence"] for row in result])[0] == "NAU"


@pytest.mark.parametrize("aligned", [
    {"family000002": "CAT", "family000001": "A-T"},
    {"family000001": "A-T"},
    {"family000001": "A-T", "family000002": "CA"},
    {"family000001": "A-C", "family000002": "CAT"},
    {"family000001": "A--T", "family000002": "C-AT"},
    {"family000001": "A" + "-" * 999 + "T", "family000002": "CA" + "-" * 998 + "T"},
])
def test_alignment_rejects_identity_shape_source_and_gap_faults(aligned):
    groups = [{"id": "family000001", "member_ids": ["a"], "sequence": "AU"},
              {"id": "family000002", "member_ids": ["b"], "sequence": "CAU"}]
    with pytest.raises(ValueError):
        core.alignment_rows(groups, aligned)


@pytest.mark.parametrize("alignment", [[], ["A" * 30], ["A" * 30, "C" * 31],
                                        ["", ""], ["A" * 1001, "C" * 1001], ["A" * 30] * 51])
def test_fold_rejects_unsafe_shapes_before_native(monkeypatch, alignment):
    RNA = pytest.importorskip("RNA")
    monkeypatch.setattr(RNA, "fold_compound", lambda *args: pytest.fail("Invalid shape reached native folding"))
    with pytest.raises(ValueError, match="shape"):
        core.fold(alignment)


@pytest.mark.parametrize("fault", ["short_structure", "unsupported_structure", "unbalanced_structure",
                                  "mfe_nan", "ensemble_nan", "matrix_shape", "probability_nan",
                                  "probability_negative", "probability_above_one"])
def test_invalid_native_evidence_is_rejected(monkeypatch, fault):
    RNA = pytest.importorskip("RNA")

    class Compound:
        params = exp_params = None

        def params_subst(self, parameters):
            pass

        def exp_params_subst(self, parameters):
            pass

        def mfe(self):
            structure = "." * 36
            if fault == "short_structure":
                structure = "." * 35
            if fault == "unsupported_structure":
                structure = "[" + "." * 34 + "]"
            if fault == "unbalanced_structure":
                structure = "(" + "." * 35
            return structure, float("nan") if fault == "mfe_nan" else 0.0

        def exp_params_rescale(self, mfe):
            assert fault not in {"short_structure", "unsupported_structure", "unbalanced_structure", "mfe_nan"}

        def pf(self):
            return "", float("nan") if fault == "ensemble_nan" else 0.0

        def bpp(self):
            if fault == "matrix_shape":
                return [[0.0]]
            matrix = [[0.0] * 37 for _ in range(37)]
            matrix[1][9] = {"probability_nan": float("nan"), "probability_negative": -0.01,
                            "probability_above_one": 1.01}.get(fault, 0.0)
            return matrix

    monkeypatch.setattr(RNA, "fold_compound", lambda *args: Compound())
    with pytest.raises(ValueError):
        core.fold(SEQUENCES)


def test_native_timeout_does_not_install_report(tmp_path, monkeypatch):
    def timeout(command, **kwargs):
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, "v7.505", "")
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", timeout)
    output = tmp_path / "timeout"
    with pytest.raises(subprocess.TimeoutExpired):
        run(FIXTURE, output, molecule="DNA", timeout_seconds=1)
    assert not (output / "report.json").exists()


def test_compensatory_count_matches_observed_sequence_pairs():
    from itertools import combinations
    rows = [left + right for left in "ACGU-" for right in "ACGU-"] * 2
    observed = sum(a in {"AU", "UA", "GC", "CG", "GU", "UG"} and
                   b in {"AU", "UA", "GC", "CG", "GU", "UG"} and
                   a[0] != b[0] and a[1] != b[1] for a, b in combinations(rows, 2))
    result = core.pair_evidence(rows, 1, 2, 0.5)
    assert result["compensatory_sequence_pairs"] == observed
    assert result["canonical_sequences"] + result["noncanonical_sequences"] + result["gapped_sequences"] == 50
