#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def option(args, name):
    return args[args.index(name) + 1]


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def has_case(records, name):
    return any(name in path.name for path in records)


def render_tree(leaves, mutation=None):
    rendered = list(leaves)
    if mutation == "wrong":
        rendered[0] = "unrelated_leaf"
    elif mutation == "missing":
        rendered.pop()
    elif mutation == "extra":
        rendered.append("unexpected_leaf")
    elif mutation == "duplicate":
        rendered[-1] = rendered[0]
    branches = [f"{leaf}:1" for leaf in rendered]
    if mutation == "missing_branch":
        branches[0] = rendered[0]
    return "(" + ",".join(branches) + ");"


def doublefork_helper(pid_path):
    if os.fork():
        os._exit(0)
    os.setsid()
    if os.fork():
        os._exit(0)
    Path(pid_path).write_text(str(os.getpid()), encoding="ascii")
    time.sleep(30)


def main():

    args = sys.argv[1:]
    if not args or args[0] != "run":
        return 2
    input_root = Path(args[1])
    output = Path(args[2])
    tool = option(args, "--tools")
    output.mkdir(parents=True)
    (output / "logs").mkdir()
    (output / "logs" / "integration.log").write_text(
        "safe integration log\n", encoding="utf-8"
    )
    records = sorted(input_root.glob("*.fasta"))

    if has_case(records, "doublefork"):
        pid_path = output / "doublefork.pid"
        subprocess.Popen(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "doublefork-helper",
                str(pid_path),
            ],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            close_fds=True,
        )
        deadline = time.monotonic() + 3
        while not pid_path.is_file() and time.monotonic() < deadline:
            time.sleep(0.01)
        if not pid_path.is_file():
            return 3
    if has_case(records, "slow"):
        (output / "runner.pid").write_text(str(os.getpid()), encoding="ascii")
        time.sleep(30)
    if has_case(records, "detach"):
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
            close_fds=True,
        )
        (output / "detached.pid").write_text(str(child.pid), encoding="ascii")
    if has_case(records, "logflood"):
        sys.stdout.buffer.write(b"x" * 262_144)
        sys.stdout.buffer.flush()
        time.sleep(30)
    if has_case(records, "diskflood"):
        with (output / "flood.bin").open("wb") as stream:
            stream.truncate(12_000_000)
        time.sleep(30)
    if has_case(records, "envcheck"):
        write_json(
            output / "environment.json",
            {
                name: os.environ.get(name)
                for name in ("TMPDIR", "TMP", "TEMP", "XDG_CACHE_HOME")
            },
        )

    no_eligible = has_case(records, "noeligible")
    invalid_contract = has_case(records, "invalidcontract")
    invalid_adapter = has_case(records, "invalidadapter")
    invalid_source_binding = has_case(records, "stalesource")
    invalid_result_groups = has_case(records, "invalidgroups")
    path_leak = has_case(records, "pathleak")
    header_only = has_case(records, "headeronly")
    invalid_tree = has_case(records, "invalidtree")
    invalid_orientation = has_case(records, "invalidorientation")
    compatible_minor = has_case(records, "compatibleminor")
    unknown_schema_major = has_case(records, "unknownschema")
    malformed_schema = has_case(records, "malformedschema")
    if compatible_minor:
        identify_schema_version = "1.8.4"
        adapter_schema_version = "1.9.0"
    elif unknown_schema_major:
        identify_schema_version = "2.0.0"
        adapter_schema_version = "2.0.0"
    elif malformed_schema:
        identify_schema_version = "1.01.0"
        adapter_schema_version = "1.1.0"
    else:
        identify_schema_version = "1.1.0"
        adapter_schema_version = "1.1.0"

    for record in records:
        source_id = record.stem.split("_", 1)[1]
        if "missingreport" in record.name:
            continue
        record_sequence = "".join(record.read_text(encoding="ascii").splitlines()[1:])
        source_sha256 = hashlib.sha256(
            record_sequence.upper().encode("ascii")
        ).hexdigest()
        report = {
            "schema": {
                "name": "CRISPRidentify-report",
                "version": identify_schema_version,
            },
            "additive_extension": {"accepted": True} if compatible_minor else None,
            "source": {
                "id": source_id,
                "length": len(record_sequence),
                "sha256": "0" * 64 if invalid_source_binding else source_sha256,
                "sha256_scope": "uppercase_sequence",
            },
            "arrays": [
                {
                    "id": "/tmp/private-array" if path_leak else f"{source_id}-bf",
                    "category": "Bona-fide",
                    "certainty_score": 0.82,
                    "repeat_count": 3,
                    "spacer_count": 2,
                    "orientation": {"strand": "+", "label": "Forward"},
                    "consensus_repeat": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
                    "spacers": [{"sequence": "ACGT"}, {"sequence": "TGCA"}],
                },
                {
                    "id": f"{source_id}-possible",
                    "category": "Possible",
                    "certainty_score": 0.61,
                    "repeat_count": 3,
                    "spacer_count": 2,
                    "orientation": {"strand": "-", "label": "Reverse"},
                    "consensus_repeat": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
                    "spacers": [{"sequence": "ACGT"}, {"sequence": "AAAA"}],
                },
                {
                    "id": f"{source_id}-unknown",
                    "category": "Possible",
                    "certainty_score": 0.55,
                    "orientation": {"strand": None, "label": "Not determined"},
                    "consensus_repeat": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
                    "spacers": [{"sequence": "ACGT"}],
                },
            ],
        }
        write_json(output / "crispridentify" / source_id / "report.json", report)
        (output / "crispridentify" / source_id / "Summary.csv").write_text(
            "ID,Category\narray,Bona-fide\n", encoding="utf-8"
        )

    group_arrays = []
    for record in records:
        source_id = record.stem.split("_", 1)[1]
        for suffix, strand in (("bf", "pos_strand"), ("possible", "neg_strand")):
            leaf = f"{source_id}_{suffix}"
            group_arrays.append(
                {
                    "array_id": f"{source_id}-{suffix}",
                    "source_id": source_id,
                    "ccdb_file": f"ccdb/group_1/{strand}/{leaf}.fasta",
                    "report_path": f"crispridentify/{source_id}/report.json",
                }
            )
    expected_leaves = [Path(entry["ccdb_file"]).stem for entry in group_arrays]
    global_tree_mutation = next(
        (
            mutation
            for marker, mutation in (
                ("treemissingbranch", "missing_branch"),
                ("treewrong", "wrong"),
                ("treemissing", "missing"),
                ("treeextra", "extra"),
                ("treeduplicate", "duplicate"),
            )
            if has_case(records, marker)
        ),
        None,
    )
    skipped = [
        {
            "array_id": "unknown",
            "source_id": "sample",
            "reason": (
                "/etc/private-skip-reason"
                if path_leak
                else "orientation_not_determined"
            ),
        }
    ]
    manifest = {
        "schema": {
            "name": "crispr-tool-integration-adapter",
            "version": "0.0.0" if invalid_adapter else adapter_schema_version,
        },
        "input": "/srv/crispr/jobs/private/crispridentify",
        "runner_config": "/etc/crispr-tools/runners.json",
        "input_array_count": len(records) * 3,
        "selected_array_count": len(records) * 2,
        "emitted_array_count": (
            0
            if no_eligible
            else len(group_arrays) + int(has_case(records, "invalidarraycount"))
        ),
        "emitted_group_count": 0 if no_eligible else 1,
        "groups": (
            []
            if no_eligible
            else [
                {
                    "name": "group_1",
                    "array_count": len(group_arrays),
                    "repeat_key": "/opt/private-repeat" if path_leak else "ACGT",
                    "arrays": group_arrays,
                }
            ]
        ),
        "skipped_arrays": skipped,
        "additive_extension": {"accepted": True} if compatible_minor else None,
    }
    write_json(output / "adapter" / "manifest.json", manifest)

    stages = {
        "crispridentify": {
            "status": "completed",
            "duration_seconds": 0.01,
            "command": ["/operator/crispridentify", "--cpu", option(args, "--cpu")],
            "return_code": 0,
        },
        "adapter": {
            "status": "completed",
            "duration_seconds": 0.01,
            "emitted_group_count": manifest["emitted_group_count"],
            "emitted_array_count": manifest["emitted_array_count"],
        },
    }
    status = "completed"
    if no_eligible and tool != "identify":
        status = "completed_no_eligible_groups"
        stage_name = "crispr_evor" if tool == "crispr-evor" else "spacerplacer"
        stages[stage_name] = {
            "status": "skipped",
            "reason": "adapter produced no group with enough related arrays",
        }
    elif tool == "crispr-evor":
        evor = output / "crispr_evor"
        (evor / "detailed_results").mkdir(parents=True)
        if has_case(records, "forwardprediction"):
            prediction = "Forward"
            recommended_reverse = False
            csv_recommendation = "False"
            decisive = True
            forward_likelihood = -4.0
            reverse_likelihood = -10.0
        elif has_case(records, "ndprediction"):
            prediction = "ND"
            recommended_reverse = None
            csv_recommendation = "ND"
            decisive = False
            forward_likelihood = -5.0
            reverse_likelihood = -4.0
        else:
            prediction = "Reverse"
            recommended_reverse = True
            csv_recommendation = "True"
            decisive = True
            forward_likelihood = -10.0
            reverse_likelihood = -4.0
        likelihood_difference = forward_likelihood - reverse_likelihood
        selected_likelihood = (
            reverse_likelihood if prediction == "Reverse" else forward_likelihood
        )
        alternate_likelihood = (
            forward_likelihood if prediction == "Reverse" else reverse_likelihood
        )
        reported_tree_length = float(len(expected_leaves))
        if has_case(records, "roundedtreelength"):
            reported_tree_length += len(expected_leaves) * 4e-6
        if has_case(records, "treelengthmismatch"):
            reported_tree_length += 1
        model = "arbitrary" if has_case(records, "invalidmodel") else "BDM"
        insertion_count = (
            "-1"
            if has_case(records, "negativeinsertions")
            else "1.5"
            if has_case(records, "fractionalinsertions")
            else "2"
        )
        csv_selected_likelihood = (
            "NaN" if has_case(records, "nanlikelihood") else str(selected_likelihood)
        )
        fields = [
            "name",
            "predicted orientation",
            "recommend reversing array",
            "confidence threshold",
            "ln_lh_bdm",
            "reversed_ln_lh_bdm",
            "ln_lh_bdm - reversed_ln_lh_bdm",
            "Deletion model preferred by LRT",
            "nb of reconstructed insertions",
            "tree length",
        ]
        with (evor / "0_results.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            if not header_only:
                writer.writerow(
                    {
                        "name": "group_1",
                        "predicted orientation": prediction,
                        "recommend reversing array": csv_recommendation,
                        "confidence threshold": "5",
                        "ln_lh_bdm": csv_selected_likelihood,
                        "reversed_ln_lh_bdm": str(alternate_likelihood),
                        "ln_lh_bdm - reversed_ln_lh_bdm": str(likelihood_difference),
                        "Deletion model preferred by LRT": model,
                        "nb of reconstructed insertions": insertion_count,
                        "tree length": str(reported_tree_length),
                    }
                )
        (evor / "detailed_results" / "events.csv").write_text(
            "branch,event\nroot,insertion\n", encoding="utf-8"
        )
        difference = (
            likelihood_difference + 0.5
            if invalid_orientation
            else likelihood_difference
        )
        write_json(
            evor / "crispr_evor_results.json",
            {
                "tool": "CRISPR-evOr",
                "version": "1.0.1.post1",
                "status": "completed",
                "orientation_semantics": {
                    "Forward": "input order is the predicted acquisition orientation",
                    "Reverse": "reverse input order is predicted",
                    "ND": "likelihood difference is below the decision threshold",
                },
                "bias_correction_policy": {
                    "alpha_requested": "--disable-bias-corrections" not in args,
                    "rho_requested": "--disable-bias-corrections" not in args,
                    "rho_failure": "fail",
                    "fallback_scope": "per_group_per_orientation",
                    "orientation_likelihood_effect": "none for orientation BDM likelihood",
                },
                "tree_policy": "estimated_separately",
                "orientation_trees": "additional_data/orientation_trees.json",
                "results": [
                    {
                        "group": "group_1",
                        "prediction": prediction,
                        "recommended_reverse": recommended_reverse,
                        "decisive": decisive,
                        "confidence_threshold": 5,
                        "forward_ln_likelihood_bdm": forward_likelihood,
                        "reverse_ln_likelihood_bdm": reverse_likelihood,
                        "forward_minus_reverse_ln_likelihood_bdm": difference,
                    }
                ],
            },
        )
        selected_policy = "reverse" if prediction == "Reverse" else "forward"
        policy_trees = {}
        for policy in ("forward", "reverse", "selected"):
            mutation = global_tree_mutation
            if has_case(records, f"{policy}bad"):
                mutation = "wrong"
            policy_leaves = (
                list(reversed(expected_leaves))
                if policy == "reverse"
                or (policy == "selected" and selected_policy == "reverse")
                else expected_leaves
            )
            policy_trees[policy] = {
                "group_1": (
                    "unsafe" if invalid_tree else render_tree(policy_leaves, mutation)
                )
            }
        if has_case(records, "wrongselectedpolicy"):
            wrong_policy = "forward" if selected_policy == "reverse" else "reverse"
            policy_trees["selected"] = dict(policy_trees[wrong_policy])
        write_json(
            evor / "additional_data" / "orientation_trees.json",
            policy_trees,
        )
        (evor / "unsafe.pkl").write_bytes(b"pickle")
        (evor / "report.html").write_text("<html>unsafe</html>", encoding="utf-8")
        try:
            (evor / "external.log").symlink_to("/etc/passwd")
        except OSError:
            pass
        stages["crispr_evor"] = {
            "status": "completed",
            "duration_seconds": 0.02,
            "command": [
                "/operator/crispr-evor",
                "--alignment_workers",
                option(args, "--cpu"),
            ],
            "return_code": 0,
            "tree_policy": "estimated_separately",
            "result_count": 1,
            "result_groups": "group_1" if invalid_result_groups else ["group_1"],
        }
    elif tool == "spacerplacer":
        placer = output / "spacerplacer"
        additional = placer / "additional_data"
        additional.mkdir(parents=True)
        reported_tree_length = float(len(expected_leaves))
        if has_case(records, "roundedtreelength"):
            reported_tree_length += len(expected_leaves) * 4e-6
        if has_case(records, "treelengthmismatch"):
            reported_tree_length += 1
        header = (
            "name,deletion_rate_bdm,insertion_rate_bdm,ln_lh_bdm,"
            "Deletion model preferred by LRT,nb of spacers in alignment,"
            "nb of reconstructed insertions,nb of reconstructed deletions,"
            "tree length,run_time\n"
        )
        values = [
            "1",
            "2",
            "-10",
            "BDM",
            "3",
            "2",
            "4",
            str(reported_tree_length),
            "6",
        ]
        malformed_values = {
            "numericblank": (0, ""),
            "numerictext": (1, "garbage"),
            "numericnan": (7, "NaN"),
            "numericinf": (8, "Inf"),
            "negativerate": (0, "-1"),
            "negativetree": (7, "-1"),
            "negativeruntime": (8, "-1"),
            "negativecount": (4, "-1"),
            "fractionalcount": (6, "4.5"),
            "nanlikelihood": (2, "NaN"),
            "invalidmodel": (3, "arbitrary"),
            "negativeinsertions": (5, "-1"),
            "fractionalinsertions": (5, "1.5"),
        }
        for marker, (index, value) in malformed_values.items():
            if has_case(records, marker):
                values[index] = value
        rows = "" if header_only else f"group_1,{','.join(values)}\n"
        (placer / "0_results.csv").write_text(header + rows, encoding="utf-8")
        tree = (
            "unsafe"
            if invalid_tree
            else render_tree(expected_leaves, global_tree_mutation)
        )
        write_json(additional / "dict_nwk_trees.json", {"group_1": tree})
        (placer / "details.log").write_text(
            "reconstruction details\n", encoding="utf-8"
        )
        stages["spacerplacer"] = {
            "status": "completed",
            "duration_seconds": 0.02,
            "command": ["/operator/spacerplacer"],
            "return_code": 0,
            "result_count": 1,
            "result_groups": "group_1" if invalid_result_groups else ["group_1"],
        }

    # Native-shaped correction columns and metadata are part of the independent gate.
    if status == "completed" and tool in {"spacerplacer", "crispr-evor"}:
        stage_name = "crispr_evor" if tool == "crispr-evor" else "spacerplacer"
        active = "--disable-bias-corrections" not in args
        effective = {"alpha": "applied" if active else "disabled", "rho": "applied" if active else "disabled", "warning_code": None, "warning": None}
        policy = {"alpha_requested": active, "rho_requested": active, "rho_failure": "fallback_uncorrected", "fallback_scope": "per_group_per_orientation" if tool == "crispr-evor" else "per_group", "orientation_likelihood_effect": "none for orientation likelihood"}
        table = output / stage_name / "0_results.csv"
        with table.open(newline="") as stream:
            reader = csv.DictReader(stream)
            columns, data = list(reader.fieldnames), list(reader)
        additions = {"deletion_rate_bdm": "1", "insertion_rate_bdm": "2", "nb of reconstructed deletions": "4", "rho_bias_failure_policy": "fallback_uncorrected", "alpha_bias_correction_status": effective["alpha"], "rho_bias_correction_status": effective["rho"], "rho_bias_correction_warning_code": "", "rho_bias_correction_warning": ""}
        extra = {key: value for key, value in additions.items() if key not in columns}
        with table.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns + list(extra))
            writer.writeheader()
            writer.writerows([{**row, **extra} for row in data])
        if tool == "crispr-evor":
            path = output / stage_name / "crispr_evor_results.json"
            doc = json.loads(path.read_text())
            doc["bias_correction_policy"] = policy
            for entry in doc["results"]:
                entry["bias_corrections"] = {direction: dict(effective) for direction in ("forward", "reverse")}
            write_json(path, doc)
        stages[stage_name].update(bias_correction_policy=policy, bias_fallback_count=0, bias_warnings=[])

    integration_result = {
        "tool_options": json.loads(os.environ.get("CRISPRLOCI_TOOL_OPTIONS", "{}")),
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "status": status,
        "stages": stages,
        "fake_launcher_args": args,
    }
    if invalid_contract:
        integration_result = {"status": "completed"}
    write_json(output / "integration_result.json", integration_result)

    if has_case(records, "tampernormalized"):
        (input_root.parent / "normalized.fasta").chmod(0o600)
        (input_root.parent / "normalized.fasta").write_text(
            ">tampered\nAAAA\n", encoding="ascii"
        )
    if has_case(records, "tamperrecord"):
        records[0].chmod(0o600)
        records[0].write_text(">tampered\nAAAA\n", encoding="ascii")
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "doublefork-helper":
        doublefork_helper(sys.argv[2])
        raise SystemExit(0)
    raise SystemExit(main())
