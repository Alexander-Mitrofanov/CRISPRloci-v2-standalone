"""Trusted legacy pickle inputs must honor the same orientation CLI contract."""

import json
import pickle

import pytest

import crispr_evor
import spacerplacer
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup
from sp_model import run_experiments


@pytest.mark.parametrize('hide_names, combine_arrays', [(False, False), (True, True)])
def test_pickled_cli_forwards_tree_configuration(tmp_path, monkeypatch, hide_names, combine_arrays):
    input_path = tmp_path / "trusted.pkl"
    input_path.write_bytes(b"unused: reconstruction is intercepted")
    calls = []
    monkeypatch.setattr(spacerplacer, "run_pickled_data", lambda *args, **kwargs: calls.append(kwargs))
    flags = ['--show_spacer_name_row'] if hide_names else []
    flags += ['--combine_non_unique_arrays'] if combine_arrays else []
    spacerplacer.main([
        str(input_path), str(tmp_path / "output"), "--input_type", "pickled",
        "--allow_unsafe_pickle", "--tree_lh_fct", "ode_based", "--tree_construction_method", "upgma",
        *flags,
    ])
    assert calls[0]["tree_distance_function"] == "likelihood"
    assert calls[0]["tree_lh_fct"] == "ode_based"
    assert calls[0]["tree_construction_method"] == "upgma"
    assert calls[0]["combine_non_unique_arrays"] is combine_arrays
    assert calls[0]["show_spacer_name_row"] is not hide_names


@pytest.mark.parametrize("provided_tree", [False, True])
def test_trusted_pickled_cli_exports_valid_orientation_and_actual_trees(tmp_path, provided_tree):
    arrays = [
        CRISPRArray(name, None, None, None, None, None, None, None, spacers, None)
        for name, spacers in [("a", ["1", "2", "3", "4"]), ("b", ["-", "2", "3", "4"])]
    ]
    group = CRISPRGroup(None, arrays, name="group.with.dots")
    input_path = tmp_path / "trusted.pkl"
    with input_path.open("wb") as file:
        pickle.dump({group.name: group}, file)
    output = tmp_path / "output"
    arguments = [
        str(input_path), str(output), "--input_type", "pickled", "--allow_unsafe_pickle",
        "--no_plot_reconstruction", "--no_plot_order_graph", "--no_alpha_bias_correction",
        "--no_rho_bias_correction",
    ]
    if provided_tree:
        tree_path = tmp_path / "trees.json"
        tree_path.write_text(json.dumps({group.name: "(a:0.2,b:0.3):0;"}), encoding="utf-8")
        arguments.extend(["--tree_path", str(tree_path)])

    result_path = crispr_evor.main(arguments)
    document = json.loads(result_path.read_text(encoding="utf-8"))
    assert document["tree_policy"] == ("provided_shared" if provided_tree else "estimated_separately")
    result = document["results"][0]
    assert result["group"] == group.name
    assert result["prediction"] in {"Forward", "Reverse", "ND"}
    assert result["bias_corrections"]["forward"]["rho"] == "disabled"
    assert result["bias_corrections"]["reverse"]["rho"] == "disabled"
    trees = json.loads((output / "additional_data/orientation_trees.json").read_text(encoding="utf-8"))
    selected = "reverse" if result["prediction"] == "Reverse" else "forward"
    assert trees["selected"][group.name] == trees[selected][group.name]
    # Legacy filenames remain available to existing standalone callers.
    assert (output / "0_protocol_w_orientation.csv").is_file()
    assert (output / "0_final_oriented_protocol.csv").is_file()


def test_likelihood_alias_resolution_retains_custom_functions_and_resolves_ode(monkeypatch):
    functions = [lambda branch, loss, alpha: 1.0]
    monkeypatch.setattr(run_experiments, "load_ode_likelihood_functions", lambda: functions)
    for likelihood in (None, "simple", "trusted-cache.pickle", functions):
        assert run_experiments._resolve_likelihood_function(likelihood) is likelihood
    for alias in ("ode_based", "precise"):
        assert run_experiments._resolve_likelihood_function(alias) is functions
