import json
from pathlib import Path
import shutil

from crispr_integration.adapter import adapt_reports


FIXTURES = Path(__file__).parent / "fixtures"


def test_adapter_writes_related_arrays_and_preserves_source_orientation(tmp_path):
    output = tmp_path / "adapter"
    manifest = adapt_reports(FIXTURES / "current_report.json", output)

    assert manifest["input_array_count"] == 3
    assert manifest["selected_array_count"] == 2
    assert manifest["emitted_group_count"] == 1
    assert manifest["emitted_array_count"] == 2
    assert manifest["skipped_arrays"][0]["reason"] == "category_not_selected"
    assert manifest["schema"]["version"] == "1.1.0"
    assert manifest["input"] == "current_report.json"
    assert manifest["input_kind"] == "file"
    assert manifest["source_binding"]["status"] == "not_requested"
    assert manifest["path_contract"]["base"] == "adapter_input_root"
    assert {
        array["report_path"]
        for array in manifest["groups"][0]["arrays"]
    } == {"current_report.json"}
    assert manifest["skipped_arrays"][0]["report_path"] == "current_report.json"

    group = output / "ccdb" / manifest["groups"][0]["name"]
    positive = list((group / "pos_strand").glob("*.fasta"))
    negative = list((group / "neg_strand").glob("*.fasta"))
    assert len(positive) == 1
    assert len(negative) == 1
    assert "GGGGAAAACCCCTTTTGGGG" in negative[0].read_text(encoding="utf-8")

    on_disk = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert on_disk == manifest
    assert str(tmp_path) not in json.dumps(on_disk)


def test_adapter_records_verified_source_binding(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">genome_A example\nacgtacgtacgt\n", encoding="utf-8")

    manifest = adapt_reports(
        FIXTURES / "current_report.json",
        tmp_path / "adapter",
        expected_fasta_path=fasta,
    )

    assert manifest["source_binding"] == {
        "status": "verified",
        "identity": ["source.id", "source.length", "source.sha256"],
        "sha256_scope": "uppercase_sequence",
    }


def test_future_biological_sequences_are_always_stored_positive(tmp_path):
    output = tmp_path / "adapter"
    manifest = adapt_reports(FIXTURES / "future_report.json", output)

    group = output / "ccdb" / manifest["groups"][0]["name"]
    assert len(list((group / "pos_strand").glob("*.fasta"))) == 2
    assert not (group / "neg_strand").exists()


def test_singletons_are_a_successful_no_result(tmp_path):
    output = tmp_path / "adapter"
    manifest = adapt_reports(
        FIXTURES / "current_report.json",
        output,
        categories=("Low score",),
    )

    assert manifest["emitted_group_count"] == 0
    assert {entry["reason"] for entry in manifest["skipped_arrays"]} == {
        "category_not_selected", "orientation_not_determined"
    }


def test_unknown_source_orientation_is_excluded_instead_of_assumed_forward(tmp_path):
    output = tmp_path / "adapter"
    manifest = adapt_reports(
        FIXTURES / "current_report.json",
        output,
        categories=("Low score",),
    )

    unknown = [
        entry for entry in manifest["skipped_arrays"]
        if entry["reason"] == "orientation_not_determined"
    ]
    assert len(unknown) == 1
    assert unknown[0]["strand"] is None
    assert manifest["emitted_array_count"] == 0
    assert not list((output / "ccdb").rglob("*.fasta"))


def test_realistic_workflow_fixture_forms_one_five_array_group(tmp_path):
    output = tmp_path / "adapter"
    manifest = adapt_reports(FIXTURES / "workflow_complete.json", output)

    assert manifest["emitted_group_count"] == 1
    assert manifest["emitted_array_count"] == 5
    group = output / "ccdb" / manifest["groups"][0]["name"]
    assert len(list((group / "pos_strand").glob("*.fasta"))) == 4
    assert len(list((group / "neg_strand").glob("*.fasta"))) == 1


def test_directory_manifest_paths_are_relative_and_location_independent(tmp_path):
    manifests = []
    ccdb_names = []
    for root_name in ("first-host-root", "different-host-root"):
        input_root = tmp_path / root_name / "identify-results"
        report = input_root / "nested" / "sequence-a" / "report.json"
        report.parent.mkdir(parents=True)
        shutil.copyfile(FIXTURES / "workflow_complete.json", report)

        output = tmp_path / f"output-{root_name}"
        manifest = adapt_reports(input_root, output)
        manifests.append(manifest)
        ccdb_names.append(sorted(path.name for path in (output / "ccdb").rglob("*.fasta")))

        serialized = (output / "manifest.json").read_text(encoding="utf-8")
        assert str(tmp_path) not in serialized
        assert "/srv/" not in serialized
        assert "/opt/" not in serialized
        assert manifest["input"] == "."
        assert manifest["input_kind"] == "directory"
        assert {
            array["report_path"]
            for group in manifest["groups"]
            for array in group["arrays"]
        } == {"nested/sequence-a/report.json"}

    assert manifests[0] == manifests[1]
    assert ccdb_names[0] == ccdb_names[1]
