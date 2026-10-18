import csv
import json
from types import SimpleNamespace

import pytest

import crispr_evor


def test_prepare_arguments_forces_orientation_once():
    assert crispr_evor.prepare_arguments(["input", "output"])[-1] == "--determine_orientation"
    prepared = crispr_evor.prepare_arguments(["input", "output", "--determine_orientation"])
    assert prepared.count("--determine_orientation") == 1


@pytest.mark.parametrize(
    ("raw", "prediction", "reverse"),
    [("True", "Reverse", True), ("False", "Forward", False), ("nd", "ND", None)],
)
def test_prediction_mapping(raw, prediction, reverse):
    assert crispr_evor._prediction(raw) == (prediction, reverse)


def test_prediction_mapping_rejects_unknown_value():
    with pytest.raises(ValueError, match="Unknown orientation"):
        crispr_evor._prediction("maybe")


def test_write_orientation_json_has_unambiguous_schema(tmp_path):
    detailed = tmp_path / "detailed_results"
    detailed.mkdir()
    comparison = detailed / "0_detailed_results_w_orientation.csv"
    fields = [
        "name",
        "ln_lh_bdm",
        "reversed_ln_lh_bdm",
        "ln_lh_bdm - reversed_ln_lh_bdm",
        "recommend reversing array",
    ]
    with comparison.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "name": "group",
                "ln_lh_bdm": "-10",
                "reversed_ln_lh_bdm": "-4",
                "ln_lh_bdm - reversed_ln_lh_bdm": "-6",
                "recommend reversing array": "True",
            }
        )
    additional = tmp_path / "additional_data"
    additional.mkdir()
    (additional / "orientation_trees.json").write_text(
        json.dumps({policy: {"group": "(a,b);"} for policy in ("forward", "reverse", "selected")}),
        encoding="utf-8",
    )
    result_path = crispr_evor.write_orientation_json(tmp_path, 5)
    document = json.loads(result_path.read_text(encoding="utf-8"))
    assert document["tool"] == "CRISPR-evOr"
    assert document["status"] == "completed"
    assert document["tree_policy"] == "estimated_separately"
    assert document["orientation_trees"] == "additional_data/orientation_trees.json"
    assert document["bias_correction_policy"] == {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fail",
        "fallback_scope": "per_group_per_orientation",
        "orientation_likelihood_effect": (
            "none; rho fallback changes rate estimates, not optimizer likelihoods or the IDM/BDM LRT; "
            "orientation uses the BDM likelihood"
        ),
    }
    assert document["results"] == [
        {
            "group": "group",
            "prediction": "Reverse",
            "recommended_reverse": True,
            "decisive": True,
            "confidence_threshold": 5,
            "forward_ln_likelihood_bdm": -10.0,
            "reverse_ln_likelihood_bdm": -4.0,
            "forward_minus_reverse_ln_likelihood_bdm": -6.0,
            "bias_corrections": {
                direction: {
                    "alpha": "applied",
                    "rho": "applied",
                    "warning_code": None,
                    "warning": None,
                }
                for direction in ("forward", "reverse")
            },
        }
    ]
    assert not list(tmp_path.glob(".crispr_evor_results.*.tmp"))


def _write_inputs(tmp_path, rows, fieldnames=None, tree_groups=None):
    fields = fieldnames or [
        "name",
        "ln_lh_bdm",
        "reversed_ln_lh_bdm",
        "ln_lh_bdm - reversed_ln_lh_bdm",
        "recommend reversing array",
    ]
    detailed = tmp_path / "detailed_results"
    detailed.mkdir(exist_ok=True)
    with (detailed / "0_detailed_results_w_orientation.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    groups = tree_groups if tree_groups is not None else [row.get("name") for row in rows if row.get("name")]
    additional = tmp_path / "additional_data"
    additional.mkdir(exist_ok=True)
    (additional / "orientation_trees.json").write_text(
        json.dumps({policy: {group: "(a,b);" for group in groups} for policy in ("forward", "reverse", "selected")}),
        encoding="utf-8",
    )


def test_prepare_arguments_inserts_required_flag_before_option_terminator():
    assert crispr_evor.prepare_arguments(["--", "-input", "-output"]) == [
        "--determine_orientation",
        "--",
        "-input",
        "-output",
    ]


def test_writer_validates_threshold_schema_duplicates_and_scientific_consistency(tmp_path):
    valid = {
        "name": "group",
        "ln_lh_bdm": "-10",
        "reversed_ln_lh_bdm": "-4",
        "ln_lh_bdm - reversed_ln_lh_bdm": "-6",
        "recommend reversing array": "True",
    }
    for invalid_threshold in (-1, float("nan"), float("inf"), True, None):
        with pytest.raises(ValueError, match="confidence_threshold"):
            crispr_evor.write_orientation_json(tmp_path, invalid_threshold)

    _write_inputs(tmp_path, [valid], fieldnames=["name"])
    with pytest.raises(ValueError, match="missing fields"):
        crispr_evor.write_orientation_json(tmp_path, 5)

    _write_inputs(tmp_path, [valid, valid])
    with pytest.raises(ValueError, match="Duplicate orientation"):
        crispr_evor.write_orientation_json(tmp_path, 5)

    inconsistent = dict(valid, **{"ln_lh_bdm - reversed_ln_lh_bdm": "-7"})
    _write_inputs(tmp_path, [inconsistent])
    with pytest.raises(ValueError, match="Inconsistent likelihood"):
        crispr_evor.write_orientation_json(tmp_path, 5)

    disagreeing = dict(valid, **{"recommend reversing array": "False"})
    _write_inputs(tmp_path, [disagreeing])
    with pytest.raises(ValueError, match="disagrees"):
        crispr_evor.write_orientation_json(tmp_path, 5)


def test_writer_validates_tree_maps_and_cleans_temp_file_on_failure(tmp_path, monkeypatch):
    row = {
        "name": "group",
        "ln_lh_bdm": "-10",
        "reversed_ln_lh_bdm": "-4",
        "ln_lh_bdm - reversed_ln_lh_bdm": "-6",
        "recommend reversing array": "True",
    }
    _write_inputs(tmp_path, [row], tree_groups=[])
    with pytest.raises(ValueError, match="missing groups"):
        crispr_evor.write_orientation_json(tmp_path, 5)

    _write_inputs(tmp_path, [row])
    tree_output = tmp_path / "additional_data" / "orientation_trees.json"
    tree_document = json.loads(tree_output.read_text(encoding="utf-8"))
    tree_document["selected"]["group"] = ""
    tree_output.write_text(json.dumps(tree_document), encoding="utf-8")
    with pytest.raises(ValueError, match="invalid trees"):
        crispr_evor.write_orientation_json(tmp_path, 5)

    _write_inputs(tmp_path, [row])

    def fail_dump(*_args, **_kwargs):
        raise RuntimeError("serialization failed")

    monkeypatch.setattr(crispr_evor.json, "dump", fail_dump)
    with pytest.raises(RuntimeError, match="serialization failed"):
        crispr_evor.write_orientation_json(tmp_path, 5)
    assert not list(tmp_path.glob(".crispr_evor_results.*.tmp"))


@pytest.mark.parametrize('prediction,difference', [('True', '-6'), ('False', '6'), ('nd', '0')])
def test_writer_rejects_selected_tree_disagreeing_with_orientation(tmp_path, prediction, difference):
    row = {
        'name': 'group',
        'ln_lh_bdm': str(-10 + float(difference)),
        'reversed_ln_lh_bdm': '-10',
        'ln_lh_bdm - reversed_ln_lh_bdm': difference,
        'recommend reversing array': prediction,
    }
    _write_inputs(tmp_path, [row])
    tree_path = tmp_path / 'additional_data' / 'orientation_trees.json'
    tree_maps = json.loads(tree_path.read_text(encoding='utf-8'))
    tree_maps['selected']['group'] = '(different,leaves);'
    tree_path.write_text(json.dumps(tree_maps), encoding='utf-8')
    with pytest.raises(ValueError, match='Selected tree disagrees'):
        crispr_evor.write_orientation_json(tmp_path, 5)


def test_writer_rejects_stale_tree_groups_and_unknown_fallback_provenance(tmp_path):
    row = {
        'name': 'group', 'ln_lh_bdm': '-10', 'reversed_ln_lh_bdm': '-4',
        'ln_lh_bdm - reversed_ln_lh_bdm': '-6', 'recommend reversing array': 'True',
    }
    _write_inputs(tmp_path, [row], tree_groups=['group', 'stale'])
    with pytest.raises(ValueError, match='unexpected groups'):
        crispr_evor.write_orientation_json(tmp_path, 5)
    _write_inputs(tmp_path, [row])
    with pytest.raises(ValueError, match='explicit effective correction provenance'):
        crispr_evor.write_orientation_json(tmp_path, 5, rho_bias_failure_policy='fallback_uncorrected')


def test_main_runs_core_once_and_returns_result_path(tmp_path, monkeypatch):
    parsed = SimpleNamespace(
        output_path=str(tmp_path),
        orientation_decision_boundary=5,
        tree_path=None,
        no_alpha_bias_correction=False,
        no_rho_bias_correction=False,
        rho_bias_failure_policy="fallback_uncorrected",
    )
    calls = []

    def fake_core(arguments, **kwargs):
        calls.append((arguments, kwargs))
        return parsed

    expected = tmp_path / "crispr_evor_results.json"
    monkeypatch.setattr(crispr_evor.spacerplacer, "main", fake_core)
    monkeypatch.setattr(
        crispr_evor.spacerplacer,
        "parse_args",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("parsed twice")),
    )
    monkeypatch.setattr(crispr_evor, "write_orientation_json", lambda *_args, **_kwargs: expected)

    assert crispr_evor.main(["input", "output"]) == expected
    assert len(calls) == 1
    assert calls[0][0].count("--determine_orientation") == 1


def test_writer_records_per_orientation_rho_fallback(tmp_path):
    warning = (
        "Group 'group' (primary): rho bias correction was unavailable and uncorrected rho estimates were used. "
        "The BDM likelihood used for orientation is unchanged."
    )
    row = {
        "name": "group",
        "ln_lh_bdm": "-10",
        "reversed_ln_lh_bdm": "-4",
        "ln_lh_bdm - reversed_ln_lh_bdm": "-6",
        "recommend reversing array": "True",
        "alpha_bias_correction_status": "applied",
        "rho_bias_correction_status": "fallback_uncorrected",
        "rho_bias_failure_policy": "fallback_uncorrected",
        "rho_bias_correction_warning_code": "rho_bias_fit_unavailable_used_uncorrected",
        "rho_bias_correction_warning": warning,
        "reversed_alpha_bias_correction_status": "applied",
        "reversed_rho_bias_correction_status": "applied",
        "reversed_rho_bias_failure_policy": "fallback_uncorrected",
        "reversed_rho_bias_correction_warning_code": "",
        "reversed_rho_bias_correction_warning": "",
    }
    _write_inputs(tmp_path, [row], fieldnames=list(row))

    result_path = crispr_evor.write_orientation_json(
        tmp_path,
        5,
        rho_bias_failure_policy="fallback_uncorrected",
    )
    result = json.loads(result_path.read_text(encoding="utf-8"))["results"][0]

    assert result["bias_corrections"]["forward"] == {
        "alpha": "applied",
        "rho": "fallback_uncorrected",
        "warning_code": "rho_bias_fit_unavailable_used_uncorrected",
        "warning": warning,
    }
    assert result["bias_corrections"]["reverse"]["rho"] == "applied"


def test_writer_emits_explicit_no_results_status(tmp_path):
    _write_inputs(tmp_path, [], tree_groups=[])
    result_path = crispr_evor.write_orientation_json(tmp_path, 5)
    document = json.loads(result_path.read_text(encoding="utf-8"))
    assert document["status"] == "completed_no_results"
    assert document["results"] == []


def test_console_entry_point_discards_programmatic_path(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(crispr_evor, "main", lambda: calls.append(tmp_path / "result.json"))

    assert crispr_evor.cli() is None


def test_help_preserves_threshold_equality_and_plot_independent_events(capsys):
    import pytest
    import spacerplacer
    with pytest.raises(SystemExit) as caught:
        spacerplacer.main(["--help"])
    assert caught.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split())
    assert "less than or equal to" in help_text
    assert "Event export also works with reconstruction plotting disabled" in help_text
