import json
from pathlib import Path
import sys

import pytest

from crispr_integration.runner import (
    StageError,
    _adapter_group_leaf_names,
    _downstream_arguments,
    load_runner_config,
    run_pipeline,
)


FIXTURES = Path(__file__).parent / "fixtures"


def test_downstream_arguments_export_events_without_plotting(tmp_path):
    arguments = _downstream_arguments(
        tmp_path / "ccdb",
        tmp_path / "output",
        spacer_distance=1,
    )

    no_plot_index = arguments.index("--no_plot_reconstruction")
    assert arguments[no_plot_index + 1] == "--save_reconstructed_events"


FAKE_WORKER = '''
import json
from pathlib import Path
import shutil
import sys

kind = sys.argv[1]
args = sys.argv[2:]


def ccdb_tree(ccdb_path, group_name):
    leaves = sorted(path.stem for path in (Path(ccdb_path) / group_name).rglob("*.fasta"))
    return "(" + ",".join(f"{leaf}:1" for leaf in leaves) + ");"


if kind == "identify":
    fixture = Path(args.pop(0))
    output = Path(args[args.index("--result_folder") + 1]) / "genome_A"
    output.mkdir(parents=True)
    shutil.copyfile(fixture, output / "report.json")
elif kind == "spacerplacer":
    group = next(entry.name for entry in sorted(Path(args[0]).iterdir()) if entry.is_dir())
    output = Path(args[1])
    tree_dir = output / "additional_data"
    tree_dir.mkdir(parents=True)
    tree = ccdb_tree(args[0], group)
    alpha_requested = "--no_alpha_bias_correction" not in args
    rho_requested = "--no_rho_bias_correction" not in args
    rho_failure = args[args.index("--rho_bias_failure_policy") + 1]
    (tree_dir / "dict_nwk_trees.json").write_text(json.dumps({group: tree}))
    (tree_dir / "aligned_groups.json").write_text(json.dumps({"cache": True}))
    (output / "0_results.csv").write_text(
        "name,deletion_rate_bdm,insertion_rate_bdm,ln_lh_bdm,nb of spacers in alignment,"
        "nb of reconstructed insertions,nb of reconstructed deletions,"
        "Deletion model preferred by LRT,tree length,run_time,rho_bias_failure_policy,"
        "alpha_bias_correction_status,rho_bias_correction_status,"
        "rho_bias_correction_warning_code,rho_bias_correction_warning\\n"
        f"{group},1,2,-4,3,4,4,BDM,{tree.count(':')},6,{rho_failure},"
        f"{'applied' if alpha_requested else 'disabled'},"
        f"{'applied' if rho_requested else 'disabled'},,\\n"
    )
elif kind == "crispr_evor":
    assert "--tree_path" not in args
    assert "--aligned_groups_path" not in args
    group = next(entry.name for entry in sorted(Path(args[0]).iterdir()) if entry.is_dir())
    output = Path(args[1])
    output.mkdir(parents=True)
    alpha_requested = "--no_alpha_bias_correction" not in args
    rho_requested = "--no_rho_bias_correction" not in args
    rho_failure = args[args.index("--rho_bias_failure_policy") + 1]
    bias_corrections = {
        direction: {
            "alpha": "applied" if alpha_requested else "disabled",
            "rho": "applied" if rho_requested else "disabled",
            "warning_code": None,
            "warning": None,
        }
        for direction in ("forward", "reverse")
    }
    tree = ccdb_tree(args[0], group)
    (output / "0_results.csv").write_text(
        "name,predicted orientation,recommend reversing array,confidence threshold,ln_lh_bdm,"
        "reversed_ln_lh_bdm,ln_lh_bdm - reversed_ln_lh_bdm,"
        "nb of reconstructed insertions,Deletion model preferred by LRT,tree length,deletion_rate_bdm,insertion_rate_bdm,nb of reconstructed deletions,rho_bias_failure_policy,alpha_bias_correction_status,rho_bias_correction_status,rho_bias_correction_warning_code,rho_bias_correction_warning\\n"
        f"{group},Forward,False,5,-4,-10,6,4,BDM,{tree.count(':')},1,2,4,{rho_failure},{bias_corrections['forward']['alpha']},{bias_corrections['forward']['rho']},,\\n"
    )
    tree_dir = output / "additional_data"
    tree_dir.mkdir()
    (tree_dir / "orientation_trees.json").write_text(json.dumps({
        policy: {group: tree}
        for policy in ("forward", "reverse", "selected")
    }))
    (output / "crispr_evor_results.json").write_text(json.dumps({
        "tool": "CRISPR-evOr",
        "version": "test",
        "status": "completed",
        "orientation_semantics": {
            "Forward": "input order",
            "Reverse": "reverse input order",
            "ND": "not determined",
        },
        "bias_correction_policy": {
            "alpha_requested": alpha_requested,
            "rho_requested": rho_requested,
            "rho_failure": rho_failure,
            "fallback_scope": "per_group_per_orientation",
            "orientation_likelihood_effect": "none for orientation BDM likelihood",
        },
        "tree_policy": "estimated_separately",
        "orientation_trees": "additional_data/orientation_trees.json",
        "results": [{
            "group": group,
            "prediction": "Forward",
            "recommended_reverse": False,
            "decisive": True,
            "confidence_threshold": 5,
            "forward_ln_likelihood_bdm": -4.0,
            "reverse_ln_likelihood_bdm": -10.0,
            "forward_minus_reverse_ln_likelihood_bdm": 6.0,
            "bias_corrections": bias_corrections,
        }],
    }))
else:
    raise SystemExit(2)
'''


OUTPUT_CONTRACT_WORKER = """
import json
from pathlib import Path
import sys

scenario = sys.argv[1]
args = sys.argv[2:]
output = Path(args[1])
output.mkdir(parents=True)
group = next(entry.name for entry in sorted(Path(args[0]).iterdir()) if entry.is_dir())
leaf_names = sorted(path.stem for path in (Path(args[0]) / group).rglob("*.fasta"))


def tree_from_leaves(leaves):
    if len(leaves) == 1:
        return f"{leaves[0]}:1;"
    return "(" + ",".join(f"{leaf}:1" for leaf in leaves) + ");"


valid_tree = tree_from_leaves(leaf_names)
reverse_tree = tree_from_leaves(list(reversed(leaf_names)))


def identity_variant(kind):
    leaves = list(leaf_names)
    if kind == "wrong":
        leaves[0] = "unexpected_leaf"
    elif kind == "missing":
        leaves.pop()
    elif kind == "extra":
        leaves.append("unexpected_leaf")
    elif kind == "duplicate":
        leaves.append(leaves[0])
    else:
        raise ValueError(kind)
    return tree_from_leaves(leaves)


alpha_requested = "--no_alpha_bias_correction" not in args
rho_requested = "--no_rho_bias_correction" not in args
rho_failure = (
    args[args.index("--rho_bias_failure_policy") + 1]
    if "--rho_bias_failure_policy" in args else "fail"
)


def result_entry(
    prediction="Forward",
    result_group=None,
    *,
    forward=-4.0,
    reverse=-10.0,
    difference=6.0,
):
    return {
        "group": result_group or group,
        "prediction": prediction,
        "recommended_reverse": {"Forward": False, "Reverse": True, "ND": None}[prediction],
        "decisive": prediction != "ND",
        "confidence_threshold": 5,
        "forward_ln_likelihood_bdm": forward,
        "reverse_ln_likelihood_bdm": reverse,
        "forward_minus_reverse_ln_likelihood_bdm": difference,
        "bias_corrections": {
            direction: {
                "alpha": "applied" if alpha_requested else "disabled",
                "rho": "applied" if rho_requested else "disabled",
                "warning_code": None,
                "warning": None,
            }
            for direction in ("forward", "reverse")
        },
    }


def write_trees(groups=None, *, tree=None, policy_overrides=None):
    if groups is None:
        groups = (group,)
    policy_overrides = policy_overrides or {}
    default_trees = (
        {policy: tree for policy in ("forward", "reverse", "selected")}
        if tree is not None
        else {"forward": valid_tree, "reverse": reverse_tree, "selected": valid_tree}
    )
    additional = output / "additional_data"
    additional.mkdir()
    (additional / "orientation_trees.json").write_text(json.dumps({
        policy: {group: policy_overrides.get(policy, default_trees[policy]) for group in groups}
        for policy in ("forward", "reverse", "selected")
    }))


def write_spacerplacer_csv(groups=None, values=None):
    if groups is None:
        groups = (group,)
    selected_values = {
        "deletion_rate_bdm": 1,
        "insertion_rate_bdm": 2,
        "ln_lh_bdm": -4,
        "nb of spacers in alignment": 3,
        "nb of reconstructed insertions": 4,
        "nb of reconstructed deletions": 4,
        "Deletion model preferred by LRT": "BDM",
        "tree length": len(leaf_names),
        "run_time": 6,
        "rho_bias_failure_policy": rho_failure,
        "alpha_bias_correction_status": "applied" if alpha_requested else "disabled",
        "rho_bias_correction_status": "applied" if rho_requested else "disabled",
        "rho_bias_correction_warning_code": "",
        "rho_bias_correction_warning": "",
    }
    selected_values.update(values or {})
    columns = tuple(selected_values)
    header = "name," + ",".join(columns) + "\\n"
    rows = "".join(
        f"{name},{','.join(str(selected_values[column]) for column in columns)}\\n"
        for name in groups
    )
    (output / "0_results.csv").write_text(header + rows)


def write_spacerplacer_trees(groups=None, *, tree=None):
    if groups is None:
        groups = (group,)
    if tree is None:
        tree = valid_tree
    additional = output / "additional_data"
    additional.mkdir()
    (additional / "dict_nwk_trees.json").write_text(json.dumps({
        name: tree for name in groups
    }))


def write_evor_csv(
    groups=None,
    *,
    prediction="Forward",
    recommendation="False",
    forward=-4,
    reverse=-10,
    difference=6,
    insertions=4,
    preferred_model="BDM",
    tree_length=None,
):
    if groups is None:
        groups = (group,)
    if tree_length is None:
        tree_length = len(leaf_names)
    header = (
        "name,predicted orientation,recommend reversing array,confidence threshold,ln_lh_bdm,"
        "reversed_ln_lh_bdm,ln_lh_bdm - reversed_ln_lh_bdm,"
        "nb of reconstructed insertions,Deletion model preferred by LRT,tree length,deletion_rate_bdm,insertion_rate_bdm,nb of reconstructed deletions,rho_bias_failure_policy,alpha_bias_correction_status,rho_bias_correction_status,rho_bias_correction_warning_code,rho_bias_correction_warning\\n"
    )
    rows = "".join(
        f"{name},{prediction},{recommendation},5,{forward},{reverse},{difference},"
        f"{insertions},{preferred_model},{tree_length},1,2,4,{rho_failure},{'applied' if alpha_requested else 'disabled'},{'applied' if rho_requested else 'disabled'},,\\n"
        for name in groups
    )
    (output / "0_results.csv").write_text(header + rows)


def write_metadata(results=None, tree_policy="estimated_separately"):
    if results is None:
        results = [result_entry()]
    (output / "crispr_evor_results.json").write_text(json.dumps({
        "tool": "CRISPR-evOr",
        "version": "test",
        "status": "completed" if results else "completed_no_results",
        "orientation_semantics": {
            "Forward": "input order",
            "Reverse": "reverse input order",
            "ND": "not determined",
        },
        "bias_correction_policy": {
            "alpha_requested": alpha_requested,
            "rho_requested": rho_requested,
            "rho_failure": rho_failure,
            "fallback_scope": "per_group_per_orientation",
            "orientation_likelihood_effect": "none for orientation BDM likelihood",
        },
        "tree_policy": tree_policy,
        "orientation_trees": "additional_data/orientation_trees.json",
        "results": results,
    }))


if scenario == "spacerplacer_missing_csv":
    pass
elif scenario == "spacerplacer_empty_csv":
    (output / "0_results.csv").touch()
elif scenario == "spacerplacer_header_only":
    write_spacerplacer_csv(groups=())
elif scenario == "spacerplacer_missing_trees":
    write_spacerplacer_csv()
elif scenario == "spacerplacer_invalid_tree":
    write_spacerplacer_csv()
    write_spacerplacer_trees(tree="/nonexistent/tree.nwk")
elif scenario == "spacerplacer_negative_branch":
    write_spacerplacer_csv()
    write_spacerplacer_trees(tree=valid_tree.replace(":1", ":-1", 1))
elif scenario.startswith("spacerplacer_tree_identity_"):
    write_spacerplacer_csv()
    write_spacerplacer_trees(
        tree=identity_variant(scenario.removeprefix("spacerplacer_tree_identity_"))
    )
elif scenario == "spacerplacer_blank_rate":
    write_spacerplacer_csv(values={"deletion_rate_bdm": ""})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_text_rate":
    write_spacerplacer_csv(values={"deletion_rate_bdm": "not-a-number"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_nan_rate":
    write_spacerplacer_csv(values={"deletion_rate_bdm": "NaN"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_infinite_rate":
    write_spacerplacer_csv(values={"deletion_rate_bdm": "Inf"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_negative_rate":
    write_spacerplacer_csv(values={"deletion_rate_bdm": -1})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_fractional_count":
    write_spacerplacer_csv(values={"nb of spacers in alignment": 3.5})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_negative_count":
    write_spacerplacer_csv(values={"nb of reconstructed deletions": -1})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_negative_tree_length":
    write_spacerplacer_csv(values={"tree length": -1})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_negative_runtime":
    write_spacerplacer_csv(values={"run_time": -1})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_blank_likelihood":
    write_spacerplacer_csv(values={"ln_lh_bdm": ""})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_nan_likelihood":
    write_spacerplacer_csv(values={"ln_lh_bdm": "NaN"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_infinite_likelihood":
    write_spacerplacer_csv(values={"ln_lh_bdm": "Inf"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_negative_insertions":
    write_spacerplacer_csv(values={"nb of reconstructed insertions": -1})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_fractional_insertions":
    write_spacerplacer_csv(values={"nb of reconstructed insertions": 1.5})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_invalid_model":
    write_spacerplacer_csv(values={"Deletion model preferred by LRT": "OTHER"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_tree_length_mismatch":
    write_spacerplacer_csv(values={"tree length": 999})
    write_spacerplacer_trees()
elif scenario in {"spacerplacer_valid_fallback", "spacerplacer_fallback_missing_warning"}:
    write_spacerplacer_csv(values={
        "rho_bias_correction_status": "fallback_uncorrected",
        "rho_bias_correction_warning_code": "rho_bias_fit_unavailable_used_uncorrected",
        "rho_bias_correction_warning": (
            "rho correction unavailable; input-order BDM likelihood unchanged"
            if scenario == "spacerplacer_valid_fallback" else ""
        ),
    })
    write_spacerplacer_trees()
elif scenario == "spacerplacer_policy_mismatch":
    write_spacerplacer_csv(values={"rho_bias_failure_policy": "fail"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_invalid_alpha_status":
    write_spacerplacer_csv(values={"alpha_bias_correction_status": "disabled"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_invalid_rho_status":
    write_spacerplacer_csv(values={"rho_bias_correction_status": "disabled"})
    write_spacerplacer_trees()
elif scenario == "spacerplacer_spurious_warning":
    write_spacerplacer_csv(values={
        "rho_bias_correction_warning_code": "rho_bias_fit_unavailable_used_uncorrected",
        "rho_bias_correction_warning": "spurious",
    })
    write_spacerplacer_trees()
elif scenario == "evor_missing_csv":
    write_trees()
    write_metadata()
elif scenario == "evor_malformed_json":
    write_evor_csv()
    (output / "crispr_evor_results.json").write_text("{")
elif scenario == "evor_wrong_tree_policy":
    write_evor_csv()
    write_trees()
    write_metadata(tree_policy="provided_shared")
elif scenario == "evor_missing_orientation_trees":
    write_evor_csv()
    write_metadata()
elif scenario in {"evor_valid_fallback", "evor_fallback_missing_warning"}:
    write_evor_csv()
    write_trees()
    result = result_entry()
    result["bias_corrections"]["forward"].update({
        "rho": "fallback_uncorrected",
        "warning_code": "rho_bias_fit_unavailable_used_uncorrected",
        "warning": (
            "rho correction unavailable; orientation likelihood unchanged"
            if scenario == "evor_valid_fallback" else None
        ),
    })
    if scenario == "evor_valid_fallback":
        path = output / "0_results.csv"
        path.write_text(path.read_text().replace(
            ",applied,applied,,", ",applied,fallback_uncorrected,rho_bias_fit_unavailable_used_uncorrected,rho correction unavailable; orientation likelihood unchanged"))
    write_metadata(results=[result])
elif scenario == "evor_no_results":
    write_evor_csv(groups=())
    write_trees(groups=())
    write_metadata(results=[])
elif scenario == "evor_csv_metadata_mismatch":
    write_evor_csv(groups=("unexpected_group",))
    write_trees()
    write_metadata()
elif scenario == "evor_csv_value_mismatch":
    write_evor_csv(difference=7)
    write_trees()
    write_metadata()
elif scenario == "evor_blank_likelihood":
    write_evor_csv(forward="")
    write_trees()
    write_metadata()
elif scenario == "evor_nan_likelihood":
    write_evor_csv(forward="NaN")
    write_trees()
    write_metadata()
elif scenario == "evor_infinite_likelihood":
    write_evor_csv(forward="Inf")
    write_trees()
    write_metadata()
elif scenario == "evor_negative_insertions":
    write_evor_csv(insertions=-1)
    write_trees()
    write_metadata()
elif scenario == "evor_fractional_insertions":
    write_evor_csv(insertions=1.5)
    write_trees()
    write_metadata()
elif scenario == "evor_invalid_model":
    write_evor_csv(preferred_model="OTHER")
    write_trees()
    write_metadata()
elif scenario == "evor_tree_length_mismatch":
    write_evor_csv(tree_length=999)
    write_trees()
    write_metadata()
elif scenario == "evor_manifest_group_mismatch":
    write_evor_csv(groups=("unexpected_group",))
    write_trees(groups=("unexpected_group",))
    write_metadata(results=[result_entry(result_group="unexpected_group")])
elif scenario == "evor_inconsistent_likelihood_delta":
    write_evor_csv(difference=7)
    write_trees()
    write_metadata(results=[result_entry(difference=7)])
elif scenario == "evor_impossible_prediction":
    write_evor_csv(forward=-6, difference=4)
    write_trees()
    write_metadata(results=[result_entry(forward=-6, difference=4)])
elif scenario == "evor_invalid_newick":
    write_evor_csv()
    write_trees(tree="/nonexistent/tree.nwk")
    write_metadata()
elif scenario == "evor_negative_branch":
    write_evor_csv()
    write_trees(policy_overrides={"reverse": reverse_tree.replace(":1", ":-1", 1)})
    write_metadata()
elif scenario == "evor_valid_reverse":
    write_evor_csv(
        prediction="Reverse", recommendation="True", forward=-4, reverse=-10, difference=-6
    )
    write_trees(policy_overrides={"selected": reverse_tree})
    write_metadata(results=[result_entry(
        prediction="Reverse", forward=-10, reverse=-4, difference=-6
    )])
elif scenario == "evor_wrong_selected_policy_forward":
    write_evor_csv()
    write_trees(policy_overrides={"selected": reverse_tree})
    write_metadata()
elif scenario == "evor_wrong_selected_policy_reverse":
    write_evor_csv(
        prediction="Reverse", recommendation="True", forward=-4, reverse=-10, difference=-6
    )
    write_trees(policy_overrides={"selected": valid_tree})
    write_metadata(results=[result_entry(
        prediction="Reverse", forward=-10, reverse=-4, difference=-6
    )])
elif scenario == "evor_wrong_selected_policy_nd":
    write_evor_csv(prediction="ND", recommendation="ND", forward=-6, difference=4)
    write_trees(policy_overrides={"selected": reverse_tree})
    write_metadata(results=[result_entry(prediction="ND", forward=-6, difference=4)])
elif scenario.startswith("evor_tree_identity_"):
    policy, kind = scenario.removeprefix("evor_tree_identity_").split("_", 1)
    write_evor_csv()
    write_trees(policy_overrides={policy: identity_variant(kind)})
    write_metadata()
else:
    raise SystemExit(2)
"""


def test_complete_pipeline_uses_evor_selected_reconstruction_without_inferred_tree(tmp_path):
    worker = tmp_path / "fake_worker.py"
    worker.write_text(FAKE_WORKER, encoding="utf-8")
    config_path = tmp_path / "runners.json"
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable, str(worker), "identify", str(FIXTURES / "current_report.json")],
        "spacerplacer": [sys.executable, str(worker), "spacerplacer"],
        "crispr_evor": [sys.executable, str(worker), "crispr_evor"],
    }), encoding="utf-8")
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">genome_A\nACGTACGTACGT\n", encoding="utf-8")

    result = run_pipeline(
        fasta,
        tmp_path / "job",
        runner_config=load_runner_config(config_path),
        timeout_seconds=10,
        bias_corrections=False,
        cpu=3,
        identify_folder_jobs=2,
    )

    assert result["status"] == "completed"
    assert result["stages"]["crispridentify"]["status"] == "completed"
    identify_command = result["stages"]["crispridentify"]["command"]
    assert identify_command[identify_command.index("--accepted_only") + 1] == "True"
    assert result["stages"]["adapter"]["emitted_group_count"] == 1
    assert result["stages"]["adapter"]["source_binding"]["status"] == "verified"
    assert result["stages"]["spacerplacer"]["status"] == "fulfilled_by_crispr_evor"
    assert result["stages"]["spacerplacer"]["results"].endswith("crispr_evor/0_results.csv")
    assert result["stages"]["crispr_evor"]["tree_source"] == "crispr_evor_orientation_specific_estimates"
    assert result["stages"]["crispr_evor"]["tree_policy"] == "estimated_separately"
    assert result["stages"]["crispr_evor"]["result_status"] == "completed"
    assert result["stages"]["crispr_evor"]["result_count"] == 1
    assert result["stages"]["crispr_evor"]["result_json"].endswith("crispr_evor_results.json")
    assert result["stages"]["crispr_evor"]["orientation_trees"].endswith("orientation_trees.json")
    assert "--no_rho_bias_correction" in result["stages"]["crispr_evor"]["command"]
    evor_command = result["stages"]["crispr_evor"]["command"]
    assert evor_command[evor_command.index("--rho_bias_failure_policy") + 1] == "fallback_uncorrected"
    assert evor_command[evor_command.index("--alignment_workers") + 1] == "3"
    assert "--tree_path" not in evor_command
    assert "--aligned_groups_path" not in evor_command
    assert result["stages"]["crispr_evor"]["bias_fallback_count"] == 0
    assert result["stages"]["crispr_evor"]["bias_correction_policy"]["rho_requested"] is False
    assert result["stages"]["crispridentify"]["command"][-4:] == [
        "--cpu", "3", "--folder_jobs", "2",
    ]
    assert (tmp_path / "job" / "integration_result.json").is_file()
    assert (tmp_path / "job" / "crispr_evor" / "crispr_evor_results.json").is_file()
    assert not (tmp_path / "job" / "spacerplacer").exists()


def test_spacerplacer_only_still_runs_the_standalone_worker(tmp_path):
    worker = tmp_path / "fake_worker.py"
    worker.write_text(FAKE_WORKER, encoding="utf-8")
    config_path = tmp_path / "runners.json"
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable],
        "spacerplacer": [sys.executable, str(worker), "spacerplacer"],
        "crispr_evor": [sys.executable],
    }), encoding="utf-8")

    result = run_pipeline(
        FIXTURES / "current_report.json",
        tmp_path / "job",
        input_kind="identify-results",
        runner_config=load_runner_config(config_path),
        tools=("spacerplacer",),
        cpu=2,
    )

    assert result["status"] == "completed"
    assert result["stages"]["spacerplacer"]["status"] == "completed"
    command = result["stages"]["spacerplacer"]["command"]
    assert command[command.index("--alignment_workers") + 1] == "2"
    assert command[command.index("--rho_bias_failure_policy") + 1] == "fallback_uncorrected"
    assert result["stages"]["spacerplacer"]["bias_fallback_count"] == 0
    assert result["stages"]["spacerplacer"]["bias_correction_policy"] == {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group",
        "orientation_likelihood_effect": (
            "none; standalone SpacerPlacer preserves the supplied input orientation"
        ),
    }
    assert "crispr_evor" not in result["stages"]
    assert (tmp_path / "job" / "spacerplacer" / "0_results.csv").is_file()
    assert result["stages"]["spacerplacer"]["results"].endswith("spacerplacer/0_results.csv")


def test_no_eligible_groups_skips_workers(tmp_path):
    config_path = tmp_path / "runners.json"
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable],
        "spacerplacer": [sys.executable],
        "crispr_evor": [sys.executable],
    }), encoding="utf-8")

    result = run_pipeline(
        FIXTURES / "current_report.json",
        tmp_path / "job",
        input_kind="identify-results",
        runner_config=load_runner_config(config_path),
        categories=("Low score",),
    )

    assert result["status"] == "completed_no_eligible_groups"
    assert result["stages"]["spacerplacer"]["status"] == "skipped"
    assert result["stages"]["crispr_evor"]["status"] == "skipped"


def test_identify_only_is_complete_even_without_comparative_groups(tmp_path):
    config_path = tmp_path / "runners.json"
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable],
        "spacerplacer": [sys.executable],
        "crispr_evor": [sys.executable],
    }), encoding="utf-8")

    result = run_pipeline(
        FIXTURES / "current_report.json",
        tmp_path / "job",
        input_kind="identify-results",
        runner_config=load_runner_config(config_path),
        tools=(),
        categories=("Low score",),
    )

    assert result["status"] == "completed"
    assert result["stages"]["adapter"]["emitted_group_count"] == 0
    assert set(result["stages"]) == {"crispridentify", "adapter"}


def _adapter_manifest_with_files(files, *, array_count=None):
    arrays = [{"ccdb_file": value} for value in files]
    return {
        "emitted_group_count": 1,
        "groups": [{
            "name": "group_1",
            "array_count": len(arrays) if array_count is None else array_count,
            "arrays": arrays,
        }],
    }


def test_adapter_group_leaf_contract_accepts_exact_ccdb_paths():
    manifest = _adapter_manifest_with_files([
        "ccdb/group_1/pos_strand/array_a.fasta",
        "ccdb/group_1/neg_strand/array_b.fasta",
    ])

    assert _adapter_group_leaf_names(manifest) == {"group_1": {"array_a", "array_b"}}


@pytest.mark.parametrize("invalid_path", [
    "/ccdb/group_1/pos_strand/array_a.fasta",
    "ccdb/other_group/pos_strand/array_a.fasta",
    "ccdb/group_1/other/array_a.fasta",
    "ccdb/group_1/pos_strand/nested/array_a.fasta",
    "ccdb/group_1/pos_strand/array_a.txt",
])
def test_adapter_group_leaf_contract_rejects_noncanonical_paths(invalid_path):
    manifest = _adapter_manifest_with_files([
        invalid_path,
        "ccdb/group_1/neg_strand/array_b.fasta",
    ])

    with pytest.raises(StageError, match="ccdb_file is unsafe"):
        _adapter_group_leaf_names(manifest)


def test_adapter_group_leaf_contract_requires_two_arrays():
    manifest = _adapter_manifest_with_files(
        ["ccdb/group_1/pos_strand/array_a.fasta"],
        array_count=1,
    )

    with pytest.raises(StageError, match="invalid array metadata"):
        _adapter_group_leaf_names(manifest)


def test_adapter_group_leaf_contract_rejects_duplicate_leaf_names():
    manifest = _adapter_manifest_with_files([
        "ccdb/group_1/pos_strand/array_a.fasta",
        "ccdb/group_1/neg_strand/array_a.fasta",
    ])

    with pytest.raises(StageError, match="duplicate tree leaf"):
        _adapter_group_leaf_names(manifest)


@pytest.mark.parametrize("folder_jobs", [0, -1, True])
def test_identify_folder_jobs_must_be_a_positive_integer(tmp_path, folder_jobs):
    with pytest.raises(ValueError, match="identify_folder_jobs must be a positive integer"):
        run_pipeline(
            FIXTURES / "current_report.json",
            tmp_path / "job",
            input_kind="identify-results",
            identify_folder_jobs=folder_jobs,
        )


def _output_contract_config(tmp_path, scenario, selected_tool):
    worker = tmp_path / "output_contract_worker.py"
    worker.write_text(OUTPUT_CONTRACT_WORKER, encoding="utf-8")
    config_path = tmp_path / "runners.json"
    command = [sys.executable, str(worker), scenario]
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable],
        "spacerplacer": command if selected_tool == "spacerplacer" else [sys.executable],
        "crispr_evor": command if selected_tool == "crispr_evor" else [sys.executable],
    }), encoding="utf-8")
    return load_runner_config(config_path)


@pytest.mark.parametrize(
    ("scenario", "selected_tool", "message"),
    [
        ("spacerplacer_missing_csv", "spacerplacer", "SpacerPlacer result CSV is missing"),
        ("spacerplacer_empty_csv", "spacerplacer", "SpacerPlacer result CSV is empty"),
        (
            "spacerplacer_header_only",
            "spacerplacer",
            "SpacerPlacer result CSV groups disagree with the adapter manifest",
        ),
        ("spacerplacer_missing_trees", "spacerplacer", "SpacerPlacer tree output is missing"),
        (
            "spacerplacer_invalid_tree",
            "spacerplacer",
            "SpacerPlacer tree for group .* has a missing or unsafe node label",
        ),
        (
            "spacerplacer_negative_branch",
            "spacerplacer",
            "SpacerPlacer tree for group .* has a negative branch length",
        ),
        *[
            (
                f"spacerplacer_tree_identity_{kind}",
                "spacerplacer",
                (
                    "SpacerPlacer tree for group .* has duplicate leaf labels"
                    if kind == "duplicate"
                    else "SpacerPlacer tree for group .* leaf identities disagree with the adapter manifest"
                ),
            )
            for kind in ("wrong", "missing", "extra", "duplicate")
        ],
        (
            "spacerplacer_blank_rate",
            "spacerplacer",
            "field 'deletion_rate_bdm' is invalid",
        ),
        (
            "spacerplacer_text_rate",
            "spacerplacer",
            "field 'deletion_rate_bdm' is invalid",
        ),
        (
            "spacerplacer_nan_rate",
            "spacerplacer",
            "field 'deletion_rate_bdm' is non-finite",
        ),
        (
            "spacerplacer_infinite_rate",
            "spacerplacer",
            "field 'deletion_rate_bdm' is non-finite",
        ),
        (
            "spacerplacer_negative_rate",
            "spacerplacer",
            "field 'deletion_rate_bdm' must be non-negative",
        ),
        (
            "spacerplacer_fractional_count",
            "spacerplacer",
            "field 'nb of spacers in alignment' must be a non-negative integer",
        ),
        (
            "spacerplacer_negative_count",
            "spacerplacer",
            "field 'nb of reconstructed deletions' must be a non-negative integer",
        ),
        (
            "spacerplacer_negative_tree_length",
            "spacerplacer",
            "field 'tree length' must be non-negative",
        ),
        (
            "spacerplacer_negative_runtime",
            "spacerplacer",
            "field 'run_time' must be non-negative",
        ),
        (
            "spacerplacer_blank_likelihood",
            "spacerplacer",
            "field 'ln_lh_bdm' is invalid",
        ),
        (
            "spacerplacer_nan_likelihood",
            "spacerplacer",
            "field 'ln_lh_bdm' is non-finite",
        ),
        (
            "spacerplacer_infinite_likelihood",
            "spacerplacer",
            "field 'ln_lh_bdm' is non-finite",
        ),
        (
            "spacerplacer_negative_insertions",
            "spacerplacer",
            "field 'nb of reconstructed insertions' must be a non-negative integer",
        ),
        (
            "spacerplacer_fractional_insertions",
            "spacerplacer",
            "field 'nb of reconstructed insertions' must be a non-negative integer",
        ),
        (
            "spacerplacer_invalid_model",
            "spacerplacer",
            "field 'Deletion model preferred by LRT' is invalid",
        ),
        (
            "spacerplacer_tree_length_mismatch",
            "spacerplacer",
            "field 'tree length' disagrees with the Newick branch sum",
        ),
        (
            "spacerplacer_fallback_missing_warning",
            "spacerplacer",
            "rho fallback for group .* lacks warning provenance",
        ),
        (
            "spacerplacer_policy_mismatch",
            "spacerplacer",
            "rho failure policy for group .* disagrees with the requested policy",
        ),
        (
            "spacerplacer_invalid_alpha_status",
            "spacerplacer",
            "alpha correction status for group .* disagrees with the request",
        ),
        (
            "spacerplacer_invalid_rho_status",
            "spacerplacer",
            "rho correction status for group .* disagrees with the request",
        ),
        (
            "spacerplacer_spurious_warning",
            "spacerplacer",
            "has a spurious rho fallback warning",
        ),
        ("evor_missing_csv", "crispr_evor", "CRISPR-evOr result CSV is missing"),
        ("evor_malformed_json", "crispr_evor", "result metadata is not valid JSON"),
        ("evor_wrong_tree_policy", "crispr_evor", "reported tree_policy 'provided_shared'"),
        (
            "evor_missing_orientation_trees",
            "crispr_evor",
            "cannot resolve CRISPR-evOr orientation tree output",
        ),
        (
            "evor_fallback_missing_warning",
            "crispr_evor",
            "fallback lacks warning provenance",
        ),
        (
            "evor_no_results",
            "crispr_evor",
            "result metadata groups disagree with the adapter manifest",
        ),
        (
            "evor_csv_metadata_mismatch",
            "crispr_evor",
            "result CSV groups disagree with result metadata",
        ),
        (
            "evor_csv_value_mismatch",
            "crispr_evor",
            "CSV field 'ln_lh_bdm - reversed_ln_lh_bdm' disagrees with metadata",
        ),
        (
            "evor_blank_likelihood",
            "crispr_evor",
            "field 'ln_lh_bdm' is invalid",
        ),
        (
            "evor_nan_likelihood",
            "crispr_evor",
            "field 'ln_lh_bdm' is non-finite",
        ),
        (
            "evor_infinite_likelihood",
            "crispr_evor",
            "field 'ln_lh_bdm' is non-finite",
        ),
        (
            "evor_negative_insertions",
            "crispr_evor",
            "field 'nb of reconstructed insertions' must be a non-negative integer",
        ),
        (
            "evor_fractional_insertions",
            "crispr_evor",
            "field 'nb of reconstructed insertions' must be a non-negative integer",
        ),
        (
            "evor_invalid_model",
            "crispr_evor",
            "field 'Deletion model preferred by LRT' is invalid",
        ),
        (
            "evor_tree_length_mismatch",
            "crispr_evor",
            "field 'tree length' disagrees with the Newick branch sum",
        ),
        (
            "evor_manifest_group_mismatch",
            "crispr_evor",
            "result metadata groups disagree with the adapter manifest",
        ),
        (
            "evor_inconsistent_likelihood_delta",
            "crispr_evor",
            "has an inconsistent likelihood difference",
        ),
        (
            "evor_impossible_prediction",
            "crispr_evor",
            "prediction disagrees with its likelihood difference and threshold",
        ),
        (
            "evor_invalid_newick",
            "crispr_evor",
            "tree for group .* has a missing or unsafe node label",
        ),
        (
            "evor_negative_branch",
            "crispr_evor",
            "CRISPR-evOr 'reverse' tree for group .* has a negative branch length",
        ),
        *[
            (
                f"evor_wrong_selected_policy_{prediction}",
                "crispr_evor",
                f"selected tree for group .* disagrees with prediction '{prediction.title() if prediction != 'nd' else 'ND'}'",
            )
            for prediction in ("forward", "reverse", "nd")
        ],
        *[
            (
                f"evor_tree_identity_{policy}_{kind}",
                "crispr_evor",
                (
                    f"CRISPR-evOr '{policy}' tree for group .* has duplicate leaf labels"
                    if kind == "duplicate"
                    else (
                        f"CRISPR-evOr '{policy}' tree for group .* leaf identities disagree "
                        "with the adapter manifest"
                    )
                ),
            )
            for policy in ("forward", "reverse", "selected")
            for kind in ("wrong", "missing", "extra", "duplicate")
        ],
    ],
)
def test_successful_worker_process_requires_valid_output_contract(
    tmp_path,
    scenario,
    selected_tool,
    message,
):
    job = tmp_path / "job"
    with pytest.raises(StageError, match=message):
        run_pipeline(
            FIXTURES / "current_report.json",
            job,
            input_kind="identify-results",
            runner_config=_output_contract_config(tmp_path, scenario, selected_tool),
            tools=(selected_tool,),
        )

    result = json.loads((job / "integration_result.json").read_text(encoding="utf-8"))
    assert result["status"] == "failed"
    assert result["error"]["stage"] == selected_tool
    assert result["stages"][selected_tool]["status"] == "failed_output_validation"


def test_crispr_evor_surfaces_valid_per_direction_bias_fallback(tmp_path):
    result = run_pipeline(
        FIXTURES / "current_report.json",
        tmp_path / "job",
        input_kind="identify-results",
        runner_config=_output_contract_config(tmp_path, "evor_valid_fallback", "crispr_evor"),
        tools=("crispr_evor",),
    )

    stage = result["stages"]["crispr_evor"]
    assert result["status"] == "completed"
    assert stage["bias_fallback_count"] == 1
    assert stage["bias_warnings"] == [{
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "forward",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable; orientation likelihood unchanged",
    }]


def test_spacerplacer_surfaces_valid_input_order_bias_fallback(tmp_path):
    result = run_pipeline(
        FIXTURES / "current_report.json",
        tmp_path / "job",
        input_kind="identify-results",
        runner_config=_output_contract_config(
            tmp_path, "spacerplacer_valid_fallback", "spacerplacer"
        ),
        tools=("spacerplacer",),
    )

    stage = result["stages"]["spacerplacer"]
    assert result["status"] == "completed"
    assert stage["bias_fallback_count"] == 1
    assert stage["bias_warnings"] == [{
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "input_order",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable; input-order BDM likelihood unchanged",
    }]


def test_spacerplacer_rejects_reported_fallback_under_strict_policy(tmp_path):
    with pytest.raises(StageError, match="reported a rho fallback under a strict policy"):
        run_pipeline(
            FIXTURES / "current_report.json",
            tmp_path / "job",
            input_kind="identify-results",
            runner_config=_output_contract_config(
                tmp_path, "spacerplacer_valid_fallback", "spacerplacer"
            ),
            tools=("spacerplacer",),
            rho_bias_failure_policy="fail",
        )


def test_crispr_evor_accepts_selected_reverse_csv_semantics(tmp_path):
    result = run_pipeline(
        FIXTURES / "current_report.json",
        tmp_path / "job",
        input_kind="identify-results",
        runner_config=_output_contract_config(tmp_path, "evor_valid_reverse", "crispr_evor"),
        tools=("crispr_evor",),
    )

    assert result["status"] == "completed"
    assert result["stages"]["crispr_evor"]["result_count"] == 1


@pytest.mark.parametrize("tool", ["spacerplacer", "crispr_evor"])
def test_scientific_options_reach_native_process_argv(tmp_path, tool):
    from dataclasses import replace
    worker = tmp_path / "fake_worker.py"
    worker.write_text(FAKE_WORKER, encoding="utf-8")
    config_path = tmp_path / "runners.json"
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable, str(worker), "identify", str(FIXTURES / "current_report.json")],
        "spacerplacer": [sys.executable, str(worker), "spacerplacer"],
        "crispr_evor": [sys.executable, str(worker), "crispr_evor"],
    }), encoding="utf-8")
    config = replace(load_runner_config(config_path), tool_options={
        "crispridentify": {"--fast_run": False, "--min_repeats": 4},
        tool: {"--tree_construction_method": "nj", "--seed": 71, "--insertion_rate": 0.0002},
    })
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">genome_A\nACGTACGTACGT\n", encoding="utf-8")
    result = run_pipeline(fasta, tmp_path / "job", runner_config=config, tools=(tool,))
    assert result["status"] == "completed"
    identify = result["stages"]["crispridentify"]["command"]
    assert identify[identify.index("--fast_run") + 1] == "False"
    assert identify[identify.index("--min_repeats") + 1] == "4"
    command = result["stages"][tool]["command"]
    assert command[command.index("--tree_construction_method") + 1] == "nj"
    assert command[command.index("--seed") + 1] == "71"
    assert command[command.index("--insertion_rate") + 1] == "0.0002"


def test_evor_rejects_result_that_ignored_orientation_override(tmp_path):
    from dataclasses import replace
    worker = tmp_path / "fake_worker.py"
    worker.write_text(FAKE_WORKER, encoding="utf-8")
    config_path = tmp_path / "runners.json"
    config_path.write_text(json.dumps({
        "crispridentify": [sys.executable],
        "spacerplacer": [sys.executable],
        "crispr_evor": [sys.executable, str(worker), "crispr_evor"],
    }), encoding="utf-8")
    config = replace(load_runner_config(config_path), tool_options={
        "crispr_evor": {"--orientation_decision_boundary": 3},
    })
    with pytest.raises(StageError, match="requested orientation boundary"):
        run_pipeline(FIXTURES / "current_report.json", tmp_path / "job", input_kind="identify-results", runner_config=config, tools=("crispr_evor",))


def test_optional_context_runners_use_their_installed_entrypoints(tmp_path):
    path = tmp_path / "runners.json"
    path.write_text(json.dumps({"crispridentify": ["identify"], "spacerplacer": ["spacerplacer"], "crispr_evor": ["crispr-evor"]}))
    config = load_runner_config(path)
    assert config.phigaro == ("phigaro-v2", "predict")
    assert config.defensefinder == ("defense-context",)
    assert config.rnaalifold == ("rnaalifold-v2",)
