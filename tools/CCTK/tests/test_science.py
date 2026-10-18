import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from PIL import Image
import pytest

from cctk import (
    sequence_operations as seq,
    array_parsimony as pars,
    CRISPRdiff as diff,
    colour_schemes as colors,
    tree_operations as trees,
)
from cctk.cli import main

FIX = Path(__file__).parent / "fixtures"
ORACLE = json.loads((FIX / "oracle.json").read_text())


@pytest.mark.parametrize("case", ORACLE["pairs"])
def test_original_pairwise_alignment_and_modules(case):
    assert list(seq.needle(case["a"], case["b"])) == case["alignment"]
    arrays = pars.find_modules(pars.Array("a", case["a"]), pars.Array("b", case["b"]))
    result = [
        {
            "aligned": a.aligned,
            "modules": [
                {"type": m.type, "spacers": m.spacers, "indices": m.indices} for m in a.modules
            ],
        }
        for a in arrays
    ]
    assert result == case["modules"]


def test_original_order_colour_and_internal_labels():
    from cctk.io import read_arrays

    arrays = read_arrays(FIX / "arrays.txt")
    assert diff.decide_array_order_global_best(arrays) == ORACLE["order"]
    assert diff.get_list_score(arrays, list(arrays)) == ORACLE["score"]
    assert [list(c) for c in colors.choose_col_scheme(12, 2)] == ORACLE["colors"]
    assert trees.create_internal_node_ids(40) == ORACLE["nodes"]


def test_original_network_and_deliberate_singletons(tmp_path):
    assert main(["network", "-i", str(FIX / "arrays.txt"), "-o", str(tmp_path)]) == 0
    assert (tmp_path / "Array_network.txt").read_text() == (FIX / "Array_network.txt").read_text()
    assert (tmp_path / "Array_clusters.txt").read_text().splitlines() == ["A B C", "isolated"]
    assert (
        main(["network", "-i", str(FIX / "arrays.txt"), "-o", str(tmp_path), "--min-shared", "2"])
        == 0
    )
    assert (tmp_path / "Array_clusters.txt").read_text().splitlines() == ["A B", "C", "isolated"]


def test_real_crisprdiff_pixel_parity(tmp_path):
    out = tmp_path / "native.png"
    assert (
        main(
            [
                "crisprdiff",
                "-a",
                str(FIX / "arrays.txt"),
                "-o",
                str(out),
                "--preordered",
                "--seed",
                "2",
            ]
        )
        == 0
    )
    assert np.array_equal(np.asarray(Image.open(out)), np.asarray(Image.open(FIX / "diff.png")))


def test_original_tree_reconstruction(tmp_path, capsys):
    out = tmp_path / "tree.txt"
    assert (
        main(
            [
                "crisprtree",
                "-a",
                str(FIX / "arrays.txt"),
                "--output-arrays",
                str(out),
                "--seed",
                "3",
                "-r",
                "2",
                "-x",
                "A",
                "B",
                "C",
            ]
        )
        == 0
    )
    assert out.read_text() == (FIX / "tree_arrays.txt").read_text()
    assert "((B:0,C:11):1,A:52):0;" in capsys.readouterr().out


def test_original_simulation(tmp_path):
    assert main(["evolve", "-n", "12", "-s", "3", "-o", str(tmp_path)]) == 0
    assert (
        next(tmp_path.glob("evolved_arrays*")).read_text()
        == (FIX / "evolution_arrays.txt").read_text()
    )
    assert (
        next(tmp_path.glob("evolved_tree_file*")).read_text()
        == (FIX / "evolution_tree.nwk").read_text()
    )


def test_original_constrain(tmp_path, capsys):
    assert (
        main(
            [
                "constrain",
                "-a",
                str(FIX / "arrays.txt"),
                "-t",
                str(FIX / "constrain.nwk"),
                "-g",
                str(FIX / "mapping.txt"),
                "-o",
                str(tmp_path / "tree.png"),
            ]
        )
        == 0
    )
    assert "((genA.A:5,genB.B:5)Anc_a:5,genC.C:10)Anc_b:0;" in capsys.readouterr().out


def test_no_mutable_array_alias():
    first, second = pars.Array("a"), pars.Array("b")
    first.spacers.append("a")
    assert second.spacers == []
    source = ["a"]
    first = pars.Array("a", source)
    first.spacers.append("b")
    assert source == ["a"]


def test_zero_branches_preserved():
    import dendropy

    tree = dendropy.Tree.get(data="(a:0,b:0):0;", schema="newick")
    trees.scale_branches(tree)
    assert all(n.edge_length == 0 for n in tree)


def test_order_work_explicitly_bounded():
    with pytest.raises(ValueError, match="eight"):
        diff.decide_array_order_global_best({str(i): ["1"] for i in range(9)})


@pytest.mark.parametrize("score", [float("nan"), float("inf"), -float("inf"), 1e308, 10**400])
@pytest.mark.parametrize("field", ["match", "mismatch", "gap"])
def test_alignment_rejects_scores_that_cannot_trace_back(score, field):
    from cctk.errors import CCTKError

    with pytest.raises(CCTKError, match="scores"):
        seq.needle("AA", "AA", **{field: score})


@pytest.mark.parametrize("source,expected", [
    ("ARY", "RYT"), ("kmbvdh", "dhbvkm"), ("SWN-", "-NWS"),
    ("ACGTRYSWKMBDHVN", "NBDHVKMWSRYACGT"),
])
def test_reverse_complement_covers_the_admitted_fasta_alphabet(source, expected):
    assert seq.rev_comp(source) == expected
    assert seq.rev_comp(expected) == source


def test_deletion_only_rejected_without_loop(tmp_path):
    assert (
        main(["evolve", "-n", "2", "-a", "0", "-t", "0", "-d", "1", "-i", "1", "-o", str(tmp_path)])
        == 2
    )


def test_parsimony_polytomy_rejected(tmp_path):
    tree = tmp_path / "poly.nwk"
    tree.write_text("(genA:1,genB:1,genC:1);")
    assert (
        main(
            [
                "constrain",
                "-a",
                str(FIX / "arrays.txt"),
                "-t",
                str(tree),
                "-g",
                str(FIX / "mapping.txt"),
                "-o",
                str(tmp_path / "tree.png"),
            ]
        )
        == 2
    )


@pytest.mark.parametrize(
    "command",
    [
        "compare",
        "network",
        "crisprdiff",
        "crisprtree",
        "constrain",
        "evolve",
        "blast",
        "spacerblast",
        "minced",
        "quickrun",
    ],
)
def test_installed_commands_help(command):
    result = subprocess.run(
        [sys.executable, "-m", "cctk", command, "--help"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


def test_exhaustive_order_matches_independent_set_objective():
    from itertools import permutations

    arrays = {"b": ["1", "1", "2"], "a": ["2", "3"], "c": ["3"], "d": ["4"]}

    def score(order):
        return sum(len(set(arrays[a]) & set(arrays[b])) for a, b in zip(order, order[1:]))

    candidates = list(permutations(arrays))
    optimum = max(score(order) for order in candidates)
    first_optimum = next(order for order in candidates if score(order) == optimum)
    assert diff.decide_array_order_global_best(arrays) == list(first_optimum)
    assert all(diff.get_list_score(arrays, order) == score(order) for order in candidates)
    assert arrays["b"] == ["1", "1", "2"]


def test_local_order_keeps_original_seeded_positive_sharing_walk():
    arrays = {str(i): [str(n) for n in range(i, i + 20)] for i in range(8)}
    diff.seed(42)
    assert diff.decide_array_order_local_best(arrays, 15, 2) == (
        ["7", "5", "3", "4", "6", "1", "2", "0"], 125
    )


def test_local_order_handles_groups_routed_by_the_cli(tmp_path):
    arrays = {str(i): [str(i), str(i + 1)] for i in range(9)}
    diff.seed(2)
    order, score = diff.decide_array_order_local_best(arrays, 10, 2)
    assert sorted(order) == sorted(arrays)
    assert score == sum(len(set(arrays[a]) & set(arrays[b])) for a, b in zip(order, order[1:]))
    source = tmp_path / "arrays.txt"
    source.write_text("".join(f"{a} {' '.join(spacers)}\n" for a, spacers in arrays.items()))
    output = tmp_path / "optimized.png"
    assert main([
        "crisprdiff", "-a", str(source), "-o", str(output), "--optimize-order", "--iterations", "1"
    ]) == 0
    assert output.is_file()


@pytest.mark.parametrize("arrays", [{"only": ["x", "x"]}, {"b": ["x"], "a": ["y"]}])
def test_local_order_without_shared_neighbours_has_supplied_fallback(arrays):
    assert diff.decide_array_order_local_best(arrays, 2, 1) == (list(arrays), 0)
    order = list(arrays)[:1]
    assert diff.shuffle_random_arrays(order) is order


@pytest.mark.parametrize("reps,nits", [(0, 1), (1, 0), (1001, 1), (1, True)])
def test_local_order_rejects_unbounded_or_empty_search(reps, nits):
    with pytest.raises(ValueError, match="repetitions"):
        diff.decide_array_order_local_best({"a": ["x"], "b": ["x"]}, reps, nits)


def test_network_matches_independent_set_metrics_in_pair_order():
    from itertools import combinations
    from types import SimpleNamespace

    records = [
        SimpleNamespace(id=ident, spacers=values, repeat_id=kind)
        for ident, values, kind in [
            ("b", ["1", "1", "2"], "I"), ("a", ["2", "3"], "II"),
            ("c", ["1", "3", "3"], "III"), ("empty", [], ""), ("alone", ["4"], "V")
        ]
    ]
    expected = []
    for a, b in combinations(records, 2):
        first, second = set(a.spacers), set(b.spacers)
        shared = first & second
        if shared:
            expected.append((a.id, b.id, len(shared), len(shared) / len(first | second),
                             a.repeat_id, b.repeat_id, shared))
    assert [(e.a, e.b, e.nshared, e.jaccard, e.a_type, e.b_type, e.shared_spacers)
            for e in seq.build_network(records)] == expected
    empty = seq.NetworkEdge(records[3], records[3])
    assert empty.nshared == 0 and empty.jaccard == 0.0
    assert records[0].spacers == ["1", "1", "2"]


@pytest.mark.parametrize("spacers", [[], ["s"]])
def test_simulated_deletion_rejects_impossible_draw_states(spacers):
    from cctk import evolve
    from cctk.errors import CCTKError

    with pytest.raises(CCTKError, match="at least two"):
        evolve.do_deletion(pars.Array("array", spacers))


def test_simulated_trailer_loss_rejects_an_empty_array():
    from cctk import evolve
    from cctk.errors import CCTKError

    with pytest.raises(CCTKError, match="at least one"):
        evolve.do_trailer_loss(pars.Array("array", []))


def test_simulated_deletion_preserves_remaining_occurrences_and_draws(monkeypatch):
    from cctk import evolve

    draws = iter([3, 1])
    calls = []

    def pick(maximum, mean, stddev):
        calls.append((maximum, mean, stddev))
        return next(draws)

    monkeypatch.setattr(evolve, "pick_normal_index", pick)
    original = ["x", "y", "x", "z", "x"]
    array = pars.Array("array", original)
    assert evolve.do_deletion(array).spacers == ["x", "z", "x"]
    assert original == ["x", "y", "x", "z", "x"]
    assert calls == [(5, 2.5, 1.25), (5, 2.5, 1.25)]


def test_simulation_trailer_loss_exhaustion_is_a_domain_error(tmp_path, capsys):
    assert main([
        "evolve", "-n", "2", "-i", "1", "-a", "0", "-d", "0", "-t", "1", "-l", "100",
        "-o", str(tmp_path)
    ]) == 2
    assert "trailer loss requires at least one" in capsys.readouterr().err
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("annotations", ["A\n", "A I extra\n", "A I\nA II\n"])
def test_malformed_network_type_annotations_are_domain_errors(tmp_path, annotations):
    from cctk import file_handling
    from cctk.errors import CCTKError

    types = tmp_path / "types.txt"
    types.write_text(annotations)
    with pytest.raises(CCTKError, match="two columns and unique"):
        file_handling.read_array_types_file(str(types))
    assert main([
        "network", "-i", str(FIX / "arrays.txt"), "-t", str(types), "-o", str(tmp_path / "out")
    ]) == 2
    assert not (tmp_path / "out").exists()


def test_type_annotation_blanks_and_order_are_retained(tmp_path):
    from cctk import file_handling

    types = tmp_path / "types.txt"
    types.write_text("\nB II\n\nA I\n")
    assert list(file_handling.read_array_types_file(str(types)).items()) == [("B", "II"), ("A", "I")]
