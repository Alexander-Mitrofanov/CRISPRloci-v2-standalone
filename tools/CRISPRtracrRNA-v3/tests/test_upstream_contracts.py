import copy
import json

import pytest

from crisprtracrrna.adapters import parse_casandra, parse_identify, sha256_file
from crisprtracrrna.errors import ContractError
from crisprtracrrna.models import Sequence
from crisprtracrrna.sequence import write_fasta


def identify_report(source):
    return {
        "schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
        "coordinate_system": {
            "indexing": 1,
            "end_inclusive": True,
            "reference": "input_forward_strand",
        },
        "source": {"id": source.id, "length": len(source.bases), "sha256": source.sha256},
        "validation": {"status": "valid"},
        "arrays": [
            {
                "id": "array1",
                "category": "Bona-fide",
                "certainty_score": 0.9,
                "source_interval": {"start": 1, "end": 4},
                "source_sequence": "ACGT",
                "consensus_repeat": "ACGT",
                "orientation": {"strand": "-"},
                "validation": {"status": "valid"},
            }
        ],
    }


def test_identify_additive_updates_and_unaccepted_arrays(tmp_path):
    source = Sequence("seq", "seq", "ACGT" * 20)
    report = identify_report(source)
    report["schema"]["version"] = "1.9.0"
    report["future_optional_field"] = {"new": True}
    alternative = copy.deepcopy(report["arrays"][0])
    alternative["category"] = "Alternative"
    report["arrays"].append(alternative)
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    arrays, preserved = parse_identify(path, source)
    assert len(arrays) == 1
    assert arrays[0].strand == "-"
    assert preserved["future_optional_field"] == {"new": True}


@pytest.mark.parametrize(
    "change", ["version", "hash", "id", "coordinates", "validation", "sequence", "nan"]
)
def test_identify_fails_closed_on_incompatible_or_stale_reports(tmp_path, change):
    source = Sequence("seq", "seq", "ACGT" * 20)
    report = identify_report(source)
    if change == "version":
        report["schema"]["version"] = "2.0.0"
    elif change == "hash":
        report["source"]["sha256"] = "0" * 64
    elif change == "id":
        report["source"]["id"] = "other"
    elif change == "coordinates":
        report["coordinate_system"]["indexing"] = 0
    elif change == "validation":
        report["validation"]["status"] = "invalid"
    elif change == "sequence":
        report["arrays"][0]["source_sequence"] = "AAAA"
    elif change == "nan":
        report["arrays"][0]["certainty_score"] = float("nan")
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    with pytest.raises(ContractError):
        parse_identify(path, source)


def casandra_fixture(directory, fasta, source):
    gene = {
        "protein_id": "p1",
        "contig_id": source.id,
        "start_1based": 1,
        "end_1based_inclusive": 9,
        "prediction": {"is_cas": True},
    }
    cassette = {
        "cassette_id": "c1",
        "contig_id": source.id,
        "start_1based": 1,
        "end_1based_inclusive": 12,
        "final_classification": {"type": "II"},
    }
    (directory / "proteins.jsonl").write_text(json.dumps(gene) + "\n")
    (directory / "cassettes.jsonl").write_text(json.dumps(cassette) + "\n")
    run = {
        "schema_version": 5,
        "program": "CasAndra",
        "genes": 1,
        "cas_cassettes": 1,
        "inputs": [{"kind": "genome_fasta", "sha256": sha256_file(fasta)}],
    }
    (directory / "run.json").write_text(json.dumps(run))
    manifest = {
        "schema_version": 5,
        "files": {
            name: {"sha256": sha256_file(directory / name)}
            for name in ["proteins.jsonl", "cassettes.jsonl"]
        },
    }
    (directory / "manifest.json").write_text(json.dumps(manifest))


def test_standalone_casandra_contract_and_tampering(tmp_path):
    source = Sequence("seq", "seq", "ACGT" * 20)
    fasta = tmp_path / "input.fa"
    write_fasta(fasta, [source])
    casandra_fixture(tmp_path, fasta, source)
    result = parse_casandra(tmp_path, fasta, source)
    assert result["status"] == "computed"
    assert result["cassettes"][0]["cassette_id"] == "c1"
    (tmp_path / "cassettes.jsonl").write_text("")
    with pytest.raises(ContractError, match="checksum"):
        parse_casandra(tmp_path, fasta, source)


def test_casandra_future_additions_are_preserved_but_unknown_schema_is_rejected(tmp_path):
    source = Sequence("seq", "seq", "ACGT" * 20)
    fasta = tmp_path / "input.fa"
    write_fasta(fasta, [source])
    casandra_fixture(tmp_path, fasta, source)
    run_path = tmp_path / "run.json"
    run = json.loads(run_path.read_text())
    run["new_optional_field"] = 42
    run_path.write_text(json.dumps(run))
    assert parse_casandra(tmp_path, fasta, source)["run"]["new_optional_field"] == 42
    run["schema_version"] = 6
    run_path.write_text(json.dumps(run))
    with pytest.raises(ContractError, match="schema"):
        parse_casandra(tmp_path, fasta, source)
