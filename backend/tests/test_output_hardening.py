from __future__ import annotations

import csv
import json
import os
import zipfile

import pytest

from crispr_web.artifacts import (
    _ArchiveLimitExceeded,
    _record,
    _write_zip,
    collect_artifacts,
)
from crispr_web.summary import build_summary, sanitize_public_summary


def _summary_for_stages(tmp_path, stages, *, max_rows=10):
    output_root = tmp_path / "output"
    output_root.mkdir(parents=True)
    (output_root / "integration_result.json").write_text(
        json.dumps({"status": "completed", "stages": stages}),
        encoding="utf-8",
    )
    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=max_rows,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )
    return sanitize_public_summary(summary, private_roots=(tmp_path,))


def _write_reconstruction_bundle(
    output_root,
    *,
    orientation,
    group,
    identities,
    top_order,
    states,
    gains,
    losses,
    other_events=None,
):
    prefix = f"{group}_rec"
    root = output_root / "crispr_evor" / orientation / "reconstructed_events" / prefix
    root.mkdir(parents=True)
    (root / f"{prefix}_tree.nwk").write_text("(leaf:0.1)root:0.0;", encoding="utf-8")
    (root / f"{prefix}_spacer_names_to_numbers.json").write_text(
        json.dumps(identities), encoding="utf-8"
    )
    (root / f"{prefix}_top_order.json").write_text(
        json.dumps(top_order), encoding="utf-8"
    )
    (root / f"{prefix}_rec_spacers.json").write_text(
        json.dumps({"rec_spacers": states}), encoding="utf-8"
    )
    (root / f"{prefix}_rec_gains_losses.json").write_text(
        json.dumps({"rec_gains": gains, "rec_losses": losses}),
        encoding="utf-8",
    )
    (root / f"{prefix}_other_events.json").write_text(
        json.dumps(other_events or {}), encoding="utf-8"
    )
    return root


def test_empty_artifact_registration_fails_closed(tmp_path):
    job_root = tmp_path / "job"
    export_root = job_root / "exports"
    export_root.mkdir(parents=True)
    empty = export_root / "empty.txt"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="empty artifact"):
        _record(job_root, empty, "empty.txt", checkpoint=None)


def test_collect_artifacts_omits_zero_byte_bundle_members(tmp_path):
    job_root = tmp_path / "job"
    output_root = job_root / "attempt-001"
    identify_root = output_root / "crispridentify" / "sample"
    log_root = output_root / "logs"
    worker_log_root = job_root / "worker-logs"
    for directory in (identify_root, log_root, worker_log_root, job_root / "input"):
        directory.mkdir(parents=True)

    (job_root / "input" / "normalized.fasta").write_text(
        ">sample\nACGT\n", encoding="ascii"
    )
    (identify_root / "Summary.csv").write_text("ID,Category\n", encoding="utf-8")
    (identify_root / "Low_Score_Candidates.txt").write_text("", encoding="utf-8")
    (log_root / "empty.log").write_text("", encoding="utf-8")
    (worker_log_root / "attempt-001.stderr.log").write_text("", encoding="utf-8")

    records, warnings = collect_artifacts(
        job_root,
        output_root,
        provenance={"input": {}, "software": {}, "execution": {}, "citations": []},
        max_files=100,
        max_file_bytes=100_000,
        max_archive_bytes=1_000_000,
        private_roots=(),
    )

    assert all(record["size_bytes"] > 0 for record in records)
    assert any("Low_Score_Candidates.txt" in warning for warning in warnings)
    assert any("empty.log" in warning for warning in warnings)
    assert any("attempt-001.stderr.log" in warning for warning in warnings)
    assert "crispridentify/sample/Summary.csv" in {record["name"] for record in records}

    bundle = next(record for record in records if record["name"] == "results.zip")
    with zipfile.ZipFile(job_root / bundle["relative_path"]) as archive:
        members = [info for info in archive.infolist() if not info.is_dir()]
        assert members
        assert all(info.file_size > 0 for info in members)
        names = {info.filename for info in members}
        assert "results/crispridentify/sample/Summary.csv" in names
        assert "results/crispridentify/sample/Low_Score_Candidates.txt" not in names
        assert "results/logs/empty.log" not in names
        assert "logs/attempt-001.stderr.log" not in names


def test_adapter_summary_exposes_only_safe_exact_membership(tmp_path):
    output_root = tmp_path / "output"
    adapter_root = output_root / "adapter"
    adapter_root.mkdir(parents=True)
    (output_root / "integration_result.json").write_text(
        json.dumps({"status": "completed", "stages": {}}), encoding="utf-8"
    )
    (adapter_root / "manifest.json").write_text(
        json.dumps(
            {
                "input_array_count": 2,
                "selected_array_count": 2,
                "emitted_array_count": 2,
                "emitted_group_count": 1,
                "skipped_arrays": [],
                "groups": [
                    {
                        "name": "repeat_group_001",
                        "array_count": 2,
                        "repeat_key": "ACGTACGT",
                        "arrays": [
                            {
                                "source_id": "example_record_01",
                                "array_id": "array_01",
                                "category": "Bona-fide",
                                "spacer_count": 12,
                                "strand": "-",
                                "input_sequence_orientation": "reverse",
                                "ccdb_strand": "Reverse",
                                "report_path": "/srv/private/report.json",
                                "ccdb_file": "/srv/private/group.ccdb",
                            },
                            {
                                "source_id": "example_record_02",
                                "array_id": "array_02",
                                "category": "Possible",
                                "spacer_count": 10,
                                "strand": "+",
                                "input_sequence_orientation": "forward",
                                "ccdb_strand": "Forward",
                            },
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=10,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )

    group = summary["adapter"]["groups"][0]
    assert group["arrays"] == [
        {
            "source_id": "example_record_01",
            "array_id": "array_01",
            "category": "Bona-fide",
            "spacer_count": 12,
            "strand": "-",
            "input_sequence_orientation": "reverse",
            "ccdb_strand": "Reverse",
        },
        {
            "source_id": "example_record_02",
            "array_id": "array_02",
            "category": "Possible",
            "spacer_count": 10,
            "strand": "+",
            "input_sequence_orientation": "forward",
            "ccdb_strand": "Forward",
        },
    ]
    assert group["arrays_truncated"] is False
    assert "report_path" not in json.dumps(group)
    assert "ccdb_file" not in json.dumps(group)
    assert "/srv/private" not in json.dumps(group)


def test_wide_result_csv_exposes_only_scientific_allowlist(tmp_path):
    output_root = tmp_path / "output"
    result_root = output_root / "spacerplacer"
    result_root.mkdir(parents=True)
    junk_fields = [f"private_{index}" for index in range(2_000)]
    fields = ["name", *junk_fields, "tree length"]
    with (result_root / "0_results.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {"name": "group_1", "tree length": "2.5", **dict.fromkeys(junk_fields, "x")}
        )

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=10,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )
    row = summary["reconstruction"]["results"][0]
    assert row == {"name": "group_1", "tree length": 2.5}
    assert "private_0" not in json.dumps(summary)


def test_summary_parser_never_follows_output_symlink(tmp_path):
    output_root = tmp_path / "output"
    result_root = output_root / "spacerplacer"
    result_root.mkdir(parents=True)
    outside = tmp_path / "outside.csv"
    outside.write_text("name,tree length\nDO_NOT_LEAK,9\n", encoding="utf-8")
    try:
        (result_root / "0_results.csv").symlink_to(outside)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=10,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )
    assert summary["reconstruction"] is None
    assert "DO_NOT_LEAK" not in json.dumps(summary)


def test_orientation_reconstructions_use_private_cross_hypothesis_spacer_ids(
    tmp_path,
):
    output_root = tmp_path / "output"
    spacer_a = "ACGTACGTACGTACGTACGTACGTACGTACGA"
    spacer_b = "TTGCTTGCTTGCTTGCTTGCTTGCTTGCTTGC"
    spacer_c = "GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC"
    _write_reconstruction_bundle(
        output_root,
        orientation="0_forward",
        group="group_alpha",
        identities={spacer_a: 1, spacer_b: 2, spacer_c: 3},
        top_order=[spacer_a, spacer_b, spacer_c],
        states={"root": [1, 0, 1], "leaf": [1, 1, 0]},
        gains={"root": [1], "leaf": [2]},
        losses={"root": [], "leaf": [[3]]},
        other_events={
            "rec_contra_dict": {"root": [3]},
            "rec_duplications_dict": {"leaf": [2]},
            "rec_rearrangements_dict": {},
            "rec_reacquisition_dict": {"root": [1]},
            "rec_indep_gain_dict": {"leaf": [2]},
            "rec_other_dup_events_dict": {"leaf": [3]},
        },
    )
    _write_reconstruction_bundle(
        output_root,
        orientation="0_reversed",
        group="group_alpha",
        identities={spacer_c: 1, spacer_b: 2, spacer_a: 3},
        top_order=[spacer_c, spacer_b, spacer_a],
        states={"root": [1, 0, 1], "leaf": [0, 1, 1]},
        gains={"root": [1], "leaf": [3]},
        losses={"root": [], "leaf": [[2, 3]]},
        other_events={
            "rec_contra_dict": {"leaf": [1]},
            "rec_duplications_dict": {},
            "rec_rearrangements_dict": {"leaf": [2]},
            "rec_reacquisition_dict": {"root": [3]},
            "rec_indep_gain_dict": {},
            "rec_other_dup_events_dict": {},
        },
    )

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=10,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )

    reconstructions = summary["orientation"]["reconstructions"]
    assert summary["orientation"]["reconstructions_truncated"] is False
    assert [row["hypothesis"] for row in reconstructions] == ["input", "reverse"]
    input_row, reverse_row = reconstructions
    assert input_row["spacer_order"] == [1, 2, 3]
    assert reverse_row["spacer_order"] == [3, 2, 1]
    assert input_row["newick"] == reverse_row["newick"] == "(leaf:0.1)root:0.0;"
    assert input_row["nodes"] == [
        {
            "name": "root",
            "spacers": [1, 3],
            "gains": [1],
            "loss_blocks": [],
            "contradictions": [3],
            "duplications": [],
            "rearrangements": [],
            "reacquisitions": [1],
            "independent_gains": [],
            "other_duplication_events": [],
        },
        {
            "name": "leaf",
            "spacers": [1, 2],
            "gains": [2],
            "loss_blocks": [[3]],
            "contradictions": [],
            "duplications": [2],
            "rearrangements": [],
            "reacquisitions": [],
            "independent_gains": [2],
            "other_duplication_events": [3],
        },
    ]
    assert reverse_row["nodes"][0]["spacers"] == [3, 1]
    assert reverse_row["nodes"][1]["gains"] == [1]
    assert reverse_row["nodes"][1]["loss_blocks"] == [[2, 1]]
    assert input_row["acquisition_count"] == reverse_row["acquisition_count"] == 2
    assert input_row["deletion_count"] == 1
    assert reverse_row["deletion_count"] == 2
    rendered = json.dumps(summary)
    assert spacer_a not in rendered
    assert spacer_b not in rendered
    assert spacer_c not in rendered
    assert not any("reconstruction" in warning for warning in summary["warnings"])


def test_malformed_optional_reconstruction_events_are_warned_and_omitted(tmp_path):
    output_root = tmp_path / "output"
    private_spacer = "ACACACACACACACACACACACACACACACAC"
    root = _write_reconstruction_bundle(
        output_root,
        orientation="0_forward",
        group="group_optional",
        identities={private_spacer: 1},
        top_order=[private_spacer],
        states={"root": [1], "leaf": [1]},
        gains={"root": [1], "leaf": []},
        losses={"root": [], "leaf": []},
        other_events={"rec_contra_dict": {"root": [private_spacer]}},
    )
    prefix = "group_optional_rec"
    (root / f"{prefix}_rec_gains_losses.json").write_text(
        json.dumps(
            {
                "rec_gains": {"root": [1], "leaf": []},
                "rec_losses": {"root": [[private_spacer]], "leaf": []},
            }
        ),
        encoding="utf-8",
    )

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=10,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )

    row = summary["orientation"]["reconstructions"][0]
    assert row["nodes"][0]["spacers"] == [1]
    assert row["nodes"][0]["gains"] == [1]
    assert row["nodes"][0]["contradictions"] == []
    assert row["acquisition_count"] == 1
    assert row["deletion_count"] == 0
    assert any("branch losses were malformed" in item for item in summary["warnings"])
    assert any(
        "contradictions events were malformed" in item for item in summary["warnings"]
    )
    assert private_spacer not in json.dumps(summary)


def test_reconstruction_parser_never_follows_private_identity_symlink(tmp_path):
    output_root = tmp_path / "output"
    public_spacer = "GATCGATCGATCGATCGATCGATCGATCGATC"
    input_root = _write_reconstruction_bundle(
        output_root,
        orientation="0_forward",
        group="group_symlink",
        identities={public_spacer: 1},
        top_order=[public_spacer],
        states={"root": [1], "leaf": [1]},
        gains={"root": [1], "leaf": []},
        losses={"root": [], "leaf": []},
    )
    reverse_root = _write_reconstruction_bundle(
        output_root,
        orientation="0_reversed",
        group="group_symlink",
        identities={public_spacer: 1},
        top_order=[public_spacer],
        states={"root": [1], "leaf": [1]},
        gains={"root": [1], "leaf": []},
        losses={"root": [], "leaf": []},
    )
    assert input_root.is_dir()
    private_identity = "DO_NOT_LEAK_PRIVATE_SPACER_SEQUENCE"
    outside = tmp_path / "private-map.json"
    outside.write_text(json.dumps({private_identity: 1}), encoding="utf-8")
    mapping_path = reverse_root / "group_symlink_rec_spacer_names_to_numbers.json"
    mapping_path.unlink()
    try:
        mapping_path.symlink_to(outside)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=10,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )

    reconstructions = summary["orientation"]["reconstructions"]
    assert [row["hypothesis"] for row in reconstructions] == ["input"]
    assert any("group_symlink (reverse)" in item for item in summary["warnings"])
    assert private_identity not in json.dumps(summary)
    assert public_spacer not in json.dumps(summary)


def test_orientation_reconstruction_group_count_is_bounded(tmp_path):
    output_root = tmp_path / "output"
    for group, spacer in (("group_a", "A"), ("group_b", "C")):
        for orientation in ("0_forward", "0_reversed"):
            _write_reconstruction_bundle(
                output_root,
                orientation=orientation,
                group=group,
                identities={spacer: 1},
                top_order=[spacer],
                states={"root": [1], "leaf": [1]},
                gains={"root": [1], "leaf": []},
                losses={"root": [], "leaf": []},
            )

    summary = build_summary(
        output_root,
        max_arrays=10,
        max_rows=1,
        max_files=100,
        accepted_categories=("Bona-fide", "Possible"),
    )

    assert [
        (row["group"], row["hypothesis"])
        for row in summary["orientation"]["reconstructions"]
    ] == [("group_a", "input"), ("group_a", "reverse")]
    assert summary["orientation"]["reconstructions_truncated"] is True


def test_public_summary_has_final_serialized_byte_cap():
    summary = {
        "schema_version": "1.0.0",
        "pipeline_status": "completed",
        "detection": {"arrays": ["x" * 10_000 for _ in range(600)]},
        "warnings": [],
    }
    sanitized = sanitize_public_summary(summary, private_roots=())
    encoded = json.dumps(sanitized, separators=(",", ":"), allow_nan=False).encode(
        "utf-8"
    )
    assert len(encoded) <= 5_000_000
    assert sanitized["detection"]["arrays"] == []
    assert sanitized["detection"]["arrays_truncated"] is True
    assert any("byte limit" in warning for warning in sanitized["warnings"])


def test_pipeline_summary_surfaces_valid_bias_fallback_provenance(tmp_path):
    evor_warning = {
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "forward",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable; orientation likelihood unchanged",
    }
    spacerplacer_warning = {
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "input_order",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable; input-order BDM likelihood unchanged",
    }
    evor_policy = {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group_per_orientation",
        "orientation_likelihood_effect": "none for orientation BDM likelihood",
    }
    spacerplacer_policy = {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group",
        "orientation_likelihood_effect": (
            "none; standalone SpacerPlacer preserves the supplied input orientation"
        ),
    }
    summary = _summary_for_stages(
        tmp_path,
        {
            "crispr_evor": {
                "status": "completed",
                "bias_correction_policy": {**evor_policy, "unknown": "/srv/private"},
                "bias_fallback_count": 1,
                "bias_warnings": [{**evor_warning, "internal_path": "/srv/private"}],
            },
            "spacerplacer": {
                "status": "completed",
                "bias_correction_policy": spacerplacer_policy,
                "bias_fallback_count": 1,
                "bias_warnings": [spacerplacer_warning],
            },
        },
    )

    stages = summary["pipeline"]["stages"]
    assert stages["crispr_evor"]["bias_correction_policy"] == evor_policy
    assert stages["crispr_evor"]["bias_fallback_count"] == 1
    assert stages["crispr_evor"]["bias_warnings"] == [evor_warning]
    assert stages["spacerplacer"]["bias_correction_policy"] == spacerplacer_policy
    assert stages["spacerplacer"]["bias_fallback_count"] == 1
    assert stages["spacerplacer"]["bias_warnings"] == [spacerplacer_warning]
    assert summary["warnings"] == [
        {
            "stage": "crispr_evor",
            **evor_warning,
            "title": (
                "CRISPR-evOr bias fallback for repeat_03e58e5dd347_group_001 (forward)"
            ),
        },
        {
            "stage": "spacerplacer",
            **spacerplacer_warning,
            "title": (
                "SpacerPlacer bias fallback for "
                "repeat_03e58e5dd347_group_001 (input order)"
            ),
        },
    ]
    assert "/srv/private" not in json.dumps(summary)


def test_pipeline_summary_rejects_malformed_bias_fallback_provenance(tmp_path):
    valid_warning = {
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "forward",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable",
    }
    evor_policy = {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group_per_orientation",
        "orientation_likelihood_effect": "none for orientation BDM likelihood",
    }
    malformed = [
        {"bias_fallback_count": True, "bias_warnings": [valid_warning]},
        {"bias_fallback_count": 2, "bias_warnings": [valid_warning]},
        {
            "bias_fallback_count": 1,
            "bias_warnings": [{**valid_warning, "group": "   "}],
        },
        {
            "bias_fallback_count": 1,
            "bias_warnings": [{**valid_warning, "message": ""}],
        },
        {
            "bias_fallback_count": 1,
            "bias_warnings": [{**valid_warning, "orientation": "sideways"}],
        },
        {
            "bias_fallback_count": 1,
            "bias_warnings": [{**valid_warning, "orientation": ["forward"]}],
        },
        {
            "bias_fallback_count": 1,
            "bias_warnings": [{**valid_warning, "orientation": "input_order"}],
        },
        {
            "bias_fallback_count": 1,
            "bias_warnings": [
                {
                    key: value
                    for key, value in valid_warning.items()
                    if key != "orientation"
                }
            ],
        },
        {
            "bias_fallback_count": 1,
            "bias_warnings": [{**valid_warning, "code": "unknown_warning"}],
        },
    ]

    for index, fallback_fields in enumerate(malformed):
        summary = _summary_for_stages(
            tmp_path / str(index),
            {
                "crispr_evor": {
                    "status": "completed",
                    "bias_correction_policy": evor_policy,
                    **fallback_fields,
                }
            },
        )
        stage = summary["pipeline"]["stages"]["crispr_evor"]
        assert "bias_fallback_count" not in stage
        assert "bias_warnings" not in stage
        assert summary["warnings"] == [
            "Malformed crispr_evor bias-fallback provenance was omitted."
        ]


def test_pipeline_summary_requires_coherent_policy_for_bias_fallbacks(tmp_path):
    warning = {
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "forward",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable",
    }
    fallback_policy = {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group_per_orientation",
        "orientation_likelihood_effect": "none for orientation BDM likelihood",
    }
    incoherent_stages = [
        {
            "bias_fallback_count": 1,
            "bias_warnings": [warning],
        },
        {
            "bias_correction_policy": {**fallback_policy, "rho_failure": "fail"},
            "bias_fallback_count": 1,
            "bias_warnings": [warning],
        },
        {
            "bias_correction_policy": {
                **fallback_policy,
                "rho_requested": False,
            },
            "bias_fallback_count": 1,
            "bias_warnings": [warning],
        },
    ]

    for index, raw_stage in enumerate(incoherent_stages):
        summary = _summary_for_stages(
            tmp_path / str(index),
            {"crispr_evor": {"status": "completed", **raw_stage}},
        )
        stage = summary["pipeline"]["stages"]["crispr_evor"]
        assert "bias_fallback_count" not in stage
        assert "bias_warnings" not in stage
        assert summary["warnings"] == [
            "Malformed crispr_evor bias-fallback provenance was omitted."
        ]

    strict_zero = _summary_for_stages(
        tmp_path / "strict-zero",
        {
            "crispr_evor": {
                "status": "completed",
                "bias_correction_policy": {
                    **fallback_policy,
                    "rho_failure": "fail",
                },
                "bias_fallback_count": 0,
                "bias_warnings": [],
            }
        },
    )
    strict_stage = strict_zero["pipeline"]["stages"]["crispr_evor"]
    assert strict_stage["bias_fallback_count"] == 0
    assert strict_stage["bias_warnings"] == []
    assert strict_zero["warnings"] == []


def test_pipeline_summary_requires_exact_bias_policy_enums(tmp_path):
    valid_policy = {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group",
        "orientation_likelihood_effect": (
            "none; standalone SpacerPlacer preserves the supplied input orientation"
        ),
    }
    invalid_policies = [
        {**valid_policy, "rho_failure": "best_effort"},
        {**valid_policy, "rho_failure": ["fallback_uncorrected"]},
        {**valid_policy, "fallback_scope": "per_group_per_orientation"},
        {**valid_policy, "alpha_requested": "true"},
        {**valid_policy, "orientation_likelihood_effect": ""},
    ]

    for index, policy in enumerate(invalid_policies):
        summary = _summary_for_stages(
            tmp_path / str(index),
            {
                "spacerplacer": {
                    "status": "completed",
                    "bias_correction_policy": policy,
                }
            },
        )
        stage = summary["pipeline"]["stages"]["spacerplacer"]
        assert "bias_correction_policy" not in stage
        assert summary["warnings"] == [
            "Malformed spacerplacer bias-correction policy was omitted."
        ]


def test_pipeline_summary_bounds_and_redacts_untrusted_bias_warnings(tmp_path):
    valid_warning = {
        "group": "repeat_03e58e5dd347_group_001",
        "orientation": "input_order",
        "code": "rho_bias_fit_unavailable_used_uncorrected",
        "message": "rho correction unavailable",
    }
    policy = {
        "alpha_requested": True,
        "rho_requested": True,
        "rho_failure": "fallback_uncorrected",
        "fallback_scope": "per_group",
        "orientation_likelihood_effect": (
            "none; standalone SpacerPlacer preserves the supplied input orientation"
        ),
    }
    unsafe_cases = [
        [{**valid_warning, "group": "g" * 241}],
        [{**valid_warning, "message": "m" * 2_001}],
        [{**valid_warning, "message": "details at /srv/crispr/private.log"}],
        [{**valid_warning, "message": "Authorization: Bearer super-secret-value"}],
        [{**valid_warning, "message": f"unexpected {'T' * 48}"}],
        [valid_warning, valid_warning, valid_warning],
    ]

    for index, raw_warnings in enumerate(unsafe_cases):
        summary = _summary_for_stages(
            tmp_path / str(index),
            {
                "spacerplacer": {
                    "status": "completed",
                    "bias_correction_policy": policy,
                    "bias_fallback_count": len(raw_warnings),
                    "bias_warnings": raw_warnings,
                }
            },
            max_rows=1,
        )
        stage = summary["pipeline"]["stages"]["spacerplacer"]
        assert "bias_fallback_count" not in stage
        assert "bias_warnings" not in stage
        rendered = json.dumps(summary)
        assert "/srv/crispr/private.log" not in rendered
        assert "super-secret-value" not in rendered
        assert "T" * 48 not in rendered
        assert summary["warnings"] == [
            "Malformed spacerplacer bias-fallback provenance was omitted."
        ]


def test_zip_writer_cannot_cross_exact_archive_cap(tmp_path):
    source = tmp_path / "random.bin"
    source.write_bytes(os.urandom(64_000))
    temporary = tmp_path / "results.zip.tmp"
    cap = 128

    with pytest.raises(_ArchiveLimitExceeded):
        _write_zip(
            temporary,
            [("random.bin", source)],
            max_bytes=cap,
            checkpoint=None,
        )
    assert temporary.is_file()
    assert temporary.stat().st_size <= cap
