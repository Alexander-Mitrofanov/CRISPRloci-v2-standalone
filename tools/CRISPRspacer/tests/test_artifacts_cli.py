import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from crisprspacer.artifacts import build_reference, load_reference, validate_result, write_result
from crisprspacer.cli import main
from crisprspacer.inputs import from_spacers
from crisprspacer.io import InputError, LimitError, digest, verify_manifest, write_manifest
from crisprspacer.search import Policy
from conftest import record


def fixture():
    return from_spacers([record("ACGT", "=query")], digest(b"fixture")), [record("ACGTACGT")]


def test_result_source_metrics_tables_and_replay(tmp_path):
    inventory, targets = fixture()
    out = tmp_path / "result"
    write_result(out, inventory, targets, Policy())
    result = validate_result(out)
    assert result["status"] == "valid"
    assert result["biological_validation"] is False
    assert "'=query" in (out / "queries.tsv").read_text()


@pytest.mark.parametrize("field,value", [("substitutions", 1), ("relative_strand", "-"),
                                       ("target_start", 2), ("excluded", True), ("pam_status", "compatible")])
def test_rehashed_alignment_tampering_is_rejected(tmp_path, field, value):
    inventory, targets = fixture()
    out = tmp_path / "result"
    write_result(out, inventory, targets, Policy())
    report = json.loads((out / "analysis.json").read_text())
    report["matches"][0][field] = value
    (out / "analysis.json").write_text(json.dumps(report))
    (out / "manifest.json").unlink()
    write_manifest(out, "crisprspacer.result")
    with pytest.raises(InputError):
        validate_result(out)


def test_rehashed_removed_hit_fails_completeness_replay(tmp_path):
    inventory, targets = fixture()
    out = tmp_path / "result"
    write_result(out, inventory, targets, Policy())
    report = json.loads((out / "analysis.json").read_text())
    report["matches"].pop()
    (out / "analysis.json").write_text(json.dumps(report))
    (out / "manifest.json").unlink()
    write_manifest(out, "crisprspacer.result")
    with pytest.raises(InputError, match="replay"):
        validate_result(out)


def test_rehashed_table_corruption_rejected(tmp_path):
    inventory, targets = fixture()
    out = tmp_path / "result"
    write_result(out, inventory, targets, Policy())
    (out / "matches.tsv").write_text("incorrect\n")
    (out / "manifest.json").unlink()
    write_manifest(out, "crisprspacer.result")
    with pytest.raises(InputError, match="Tabular"):
        validate_result(out)


def test_atomic_failure_and_overwrite_refusal(tmp_path):
    inventory, targets = fixture()
    out = tmp_path / "result"
    with pytest.raises(LimitError):
        write_result(out, inventory, targets, Policy(max_hits=1))
    assert not out.exists() and not out.with_name("result.partial").exists()
    write_result(out, inventory, targets, Policy())
    before = (out / "manifest.json").read_bytes()
    with pytest.raises(InputError, match="already exists"):
        write_result(out, inventory, targets, Policy())
    assert (out / "manifest.json").read_bytes() == before


def test_reference_roundtrip_and_tampering(tmp_path):
    _, targets = fixture()
    out = tmp_path / "reference"
    build_reference(out, targets, "synthetic", digest(b"fixture"), "synthetic", "v1")
    loaded, provenance = load_reference(out)
    assert loaded == targets and provenance["name"] == "synthetic"
    with pytest.raises(InputError):
        build_reference(out, targets, "synthetic", digest(b"fixture"), "synthetic", "v1")
    content = json.loads((out / "reference.json").read_text())
    content["records"][0]["sequence"] = "AAAAAAAA"
    (out / "reference.json").write_text(json.dumps(content))
    (out / "manifest.json").unlink()
    write_manifest(out, "crisprspacer.reference")
    with pytest.raises(InputError, match="hash"):
        load_reference(out)


@pytest.mark.parametrize("field", ["record_count", "base_count"])
@pytest.mark.parametrize("value", [True, 1.0])
def test_reference_rejects_noninteger_inventory_counts(tmp_path, field, value):
    out = tmp_path / "reference"
    targets = [record("A")]
    build_reference(out, targets, "synthetic", digest(b"fixture"), "synthetic", "v1")
    content = json.loads((out / "reference.json").read_text())
    content[field] = value
    (out / "reference.json").write_text(json.dumps(content))
    (out / "manifest.json").unlink()
    write_manifest(out, "crisprspacer.reference")
    with pytest.raises(InputError, match="integer"):
        load_reference(out)


@pytest.mark.parametrize("value", [True, 1.0])
def test_manifest_rejects_noninteger_artifact_byte_counts(tmp_path, value):
    # A one-byte artifact distinguishes bool acceptance from a mere size mismatch.
    (tmp_path / "artifact").write_bytes(b"x")
    write_manifest(tmp_path, "synthetic")
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    manifest["files"]["artifact"]["bytes"] = value
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(InputError, match="integer"):
        verify_manifest(tmp_path, "synthetic", {"artifact"})


def test_manifest_roundtrip_preserves_verified_bytes(tmp_path):
    (tmp_path / "artifact").write_bytes(b"x")
    write_manifest(tmp_path, "synthetic")
    manifest, payloads = verify_manifest(tmp_path, "synthetic", {"artifact"})
    assert manifest["files"]["artifact"] == {"sha256": digest(b"x"), "bytes": 1}
    assert payloads == {"artifact": b"x"}


def test_manifest_cannot_escape_bundle(tmp_path):
    _, targets = fixture()
    out = tmp_path / "reference"
    build_reference(out, targets, "synthetic", digest(b"fixture"), "synthetic", "v1")
    manifest = json.loads((out / "manifest.json").read_text())
    manifest["files"]["../outside"] = manifest["files"].pop("reference.json")
    (out / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(InputError, match="exact required"):
        load_reference(out)


def test_cli_self_prepare_and_search(native_fixture, tmp_path, capsys):
    _, _, report, fasta = native_fixture
    prepared, self_result, searched = tmp_path / "prepared", tmp_path / "self", tmp_path / "searched"
    assert main(["prepare", "--source", str(fasta), "--identify", str(report), "--output", str(prepared)]) == 0
    assert main(["self", "--source", str(fasta), "--identify", str(report), "--output", str(self_result)]) == 0
    assert main(["search", "--inventory", str(prepared), "--targets", str(fasta), "--output", str(searched)]) == 0
    assert main(["validate", str(self_result)]) == 0
    assert main(["validate", str(searched)]) == 0
    lines = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert lines[1]["counts"]["excluded_array_matches"] == 1
    assert lines[2]["counts"]["excluded_array_matches"] == 0


def test_cli_index_and_reference_search(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">q\nACGT\n")
    reference = tmp_path / "reference"
    assert main(["index", "--targets", str(fasta), "--name", "fixture", "--output", str(reference)]) == 0
    assert main(["validate-reference", str(reference)]) == 0
    out = tmp_path / "out"
    assert main(["search", "--spacers", str(fasta), "--reference", str(reference), "--output", str(out)]) == 0
    assert validate_result(out)["status"] == "valid"


def test_cli_errors_have_nonzero_status_and_no_result(tmp_path, capsys):
    path = tmp_path / "input.fa"
    path.write_text(">q\nAAA\n")
    out = tmp_path / "out"
    assert main(["search", "--spacers", str(path), "--targets", str(path), "--max-work", "1", "--output", str(out)]) == 3
    assert not out.exists()
    assert json.loads(capsys.readouterr().err)["search_complete"] is False


def test_standalone_module_from_unrelated_working_directory(tmp_path):
    src = Path(__file__).resolve().parents[1] / "src"
    env = {**os.environ, "PYTHONPATH": str(src)}
    result = subprocess.run([sys.executable, "-m", "crisprspacer", "--version"], cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0 and "CRISPRspacer 0.3.1" in result.stdout


def test_forged_reference_inventory_rejected(tmp_path):
    inventory, targets = fixture()
    out = tmp_path / "result"
    write_result(out, inventory, targets, Policy())
    fake = {"name": "invented", "source_label": "invented", "source_release": "invented",
            "inventory_sha256": "0"*64, "reference_sha256": "1"*64}
    for name in ("inputs.json", "analysis.json"):
        data = json.loads((out / name).read_text())
        data["reference"] = fake
        (out / name).write_text(json.dumps(data))
    (out / "manifest.json").unlink()
    write_manifest(out, "crisprspacer.result")
    with pytest.raises(InputError, match="bind retained targets"):
        validate_result(out)


def test_reference_parse_uses_verified_snapshot(tmp_path, monkeypatch):
    import crisprspacer.artifacts as artifacts
    _, targets = fixture()
    out = tmp_path / "reference"
    build_reference(out, targets, "synthetic", digest(b"fixture"), "synthetic", "v1")
    verify = artifacts.verify_manifest
    def replace_after_verify(*args):
        checked = verify(*args)
        (out / "reference.json").write_text("{}")
        return checked
    monkeypatch.setattr(artifacts, "verify_manifest", replace_after_verify)
    loaded, _ = load_reference(out)
    assert loaded == targets


def test_huge_json_integer_has_structured_cli_error(tmp_path, capsys):
    path = tmp_path / "inventory.json"
    path.write_text('{"a":' + '9'*5000 + '}')
    assert main(["search", "--inventory", str(path), "--targets", str(path), "--output", str(tmp_path / "out")]) == 2
    assert json.loads(capsys.readouterr().err)["status"] == "error"
