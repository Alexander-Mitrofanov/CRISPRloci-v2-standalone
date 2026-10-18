"""External-report failures and frozen experimental gate regression coverage."""

import copy
import json
from dataclasses import replace

import pytest

from crisprtracrrna.adapters import load_json, parse_casandra, parse_identify, sha256_file
from crisprtracrrna.config import Settings
from crisprtracrrna.errors import ContractError
from crisprtracrrna.evidence import short_rescue_admission
from crisprtracrrna.models import AntiRepeat, Array, Candidate, Hit, Interval, Sequence
from crisprtracrrna.sequence import write_fasta


@pytest.fixture
def source():
    return Sequence("seq", "seq", "ACGT" * 20)


@pytest.fixture
def identify(source):
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
                "consensus_repeat": "ACGT",
                "orientation": {"strand": "-"},
                "validation": {"status": "valid"},
            }
        ],
    }


@pytest.fixture
def casandra(tmp_path, source):
    fasta = tmp_path / "input.fa"
    write_fasta(fasta, [source])
    records = {
        "proteins.jsonl": {
            "contig_id": source.id,
            "start_1based": 1,
            "end_1based_inclusive": 9,
            "prediction": {"is_cas": True},
        },
        "cassettes.jsonl": {
            "cassette_id": "c1",
            "contig_id": source.id,
            "start_1based": 1,
            "end_1based_inclusive": 12,
            "final_classification": {"type": "II"},
        },
    }
    for name, record in records.items():
        (tmp_path / name).write_text(json.dumps(record) + "\n", encoding="utf-8")
    run = {
        "schema_version": 5,
        "program": "CasAndra",
        "genes": 1,
        "cas_cassettes": 1,
        "inputs": [{"kind": "genome_fasta", "sha256": sha256_file(fasta)}],
    }
    (tmp_path / "run.json").write_text(json.dumps(run), encoding="utf-8")
    manifest = {
        "schema_version": 5,
        "files": {name: {"sha256": sha256_file(tmp_path / name)} for name in records},
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return tmp_path, fasta


def rewrite_checked_record(directory, filename, text):
    """Keep the checksum valid so tests reach the structured-record boundary."""
    path = directory / filename
    path.write_text(text, encoding="utf-8")
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"][filename]["sha256"] = sha256_file(path)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def update_run(directory, **updates):
    path = directory / "run.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    report.update(updates)
    path.write_text(json.dumps(report), encoding="utf-8")


@pytest.mark.parametrize("number", ["NaN", "Infinity", "-Infinity", "1e999", "-1e999"])
def test_structured_json_rejects_nonfinite_additive_fields(tmp_path, number):
    path = tmp_path / "report.json"
    path.write_text('{"optional": [1, {"number": ' + number + "}]}", encoding="utf-8")
    with pytest.raises(ContractError, match="Non-finite") as error:
        load_json(path)
    assert str(path) in str(error.value)


@pytest.mark.parametrize("text", ["[]", "null", "5"])
def test_structured_json_requires_an_object(tmp_path, text):
    path = tmp_path / "report.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ContractError, match="Expected JSON object"):
        load_json(path)


@pytest.mark.parametrize(
    "text",
    [
        '{"status": "failed", "status": "complete"}',
        '{"optional": {"score": 0.1, "score": 1.0}}',
        '{"optional": [{"valid": false, "valid": true}]}',
    ],
)
def test_structured_json_rejects_duplicate_keys_at_every_depth(tmp_path, text):
    path = tmp_path / "report.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ContractError, match="Duplicate JSON key") as error:
        load_json(path)
    assert str(path) in str(error.value)


@pytest.mark.parametrize("filename", ["proteins.jsonl", "cassettes.jsonl"])
def test_casandra_duplicate_json_key_has_record_context(casandra, source, filename):
    directory, fasta = casandra
    original = (directory / filename).read_text(encoding="utf-8")
    malformed = original.rstrip()[:-1] + ', "optional": 1, "optional": 2}\n'
    rewrite_checked_record(directory, filename, original + "\n" + malformed)
    with pytest.raises(ContractError, match="Duplicate JSON key") as error:
        parse_casandra(directory, fasta, source)
    assert f"{directory / filename}:3:" in str(error.value)


def test_casandra_duplicate_cassette_ids_are_not_independent_evidence(casandra, source):
    directory, fasta = casandra
    filename = "cassettes.jsonl"
    original = (directory / filename).read_text(encoding="utf-8")
    rewrite_checked_record(directory, filename, original + original)
    update_run(directory, cas_cassettes=2)
    with pytest.raises(ContractError, match="Duplicate CasAndra cassette id") as error:
        parse_casandra(directory, fasta, source)
    assert f"{directory / filename}:2:" in str(error.value)


@pytest.mark.parametrize("score", [True, False, "0.9", "1", None, [], {}])
def test_identify_rejects_coercible_and_nonnumeric_certainty(tmp_path, identify, source, score):
    identify["arrays"][0]["certainty_score"] = score
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    with pytest.raises(ContractError, match="Invalid array score"):
        parse_identify(path, source)


@pytest.mark.parametrize("score", [0, 1, 0.9])
def test_identify_retains_valid_numeric_certainty(tmp_path, identify, source, score):
    identify["arrays"][0]["certainty_score"] = score
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    arrays, preserved = parse_identify(path, source)
    assert arrays[0].score == score
    assert preserved == identify


@pytest.mark.parametrize("filename", ["proteins.jsonl", "cassettes.jsonl"])
@pytest.mark.parametrize("number", ["NaN", "Infinity", "-Infinity", "1e999", "-1e999"])
def test_casandra_rejects_nonfinite_additive_fields(casandra, source, filename, number):
    directory, fasta = casandra
    original = (directory / filename).read_text(encoding="utf-8")
    extra = json.loads(original)
    extra["optional"] = {"nonwinning": [0, "NONFINITE"]}
    malformed = json.dumps(extra).replace('"NONFINITE"', number)
    rewrite_checked_record(directory, filename, original + "\n" + malformed + "\n")
    with pytest.raises(ContractError, match="Non-finite") as error:
        parse_casandra(directory, fasta, source)
    assert f"{directory / filename}:3:" in str(error.value)


@pytest.mark.parametrize("text", ["[]", "null", "5", "{broken"])
def test_casandra_malformed_later_line_has_file_and_physical_line(casandra, source, text):
    directory, fasta = casandra
    filename = "cassettes.jsonl"
    original = (directory / filename).read_text(encoding="utf-8")
    rewrite_checked_record(directory, filename, original + "\n" + text + "\n")
    with pytest.raises(ContractError) as error:
        parse_casandra(directory, fasta, source)
    assert f"{directory / filename}:3:" in str(error.value)


def test_casandra_preserves_utf8_finite_additions_blank_lines_and_row_order(casandra, source):
    directory, fasta = casandra
    expected = {}
    for filename, key in (("proteins.jsonl", "genes"), ("cassettes.jsonl", "cassettes")):
        first = json.loads((directory / filename).read_text(encoding="utf-8"))
        first["optional"] = {"label": "première", "finite": [1e100, -1e-100, 0]}
        second = copy.deepcopy(first)
        second["optional"]["label"] = "第二"
        if key == "cassettes":
            second["cassette_id"] = "c2"
        expected[key] = [first, second]
        lines = "\n \n" + "\n\n".join(json.dumps(r, ensure_ascii=False) for r in expected[key])
        rewrite_checked_record(directory, filename, lines + "\n")
    update_run(directory, genes=2, cas_cassettes=2, optional={"finite": 1e100})
    result = parse_casandra(directory, fasta, source)
    assert result["genes"] == expected["genes"]
    assert result["cassettes"] == expected["cassettes"]
    assert result["run"]["optional"] == {"finite": 1e100}


@pytest.mark.parametrize("identifier", [None, "", 42, False, []])
def test_identify_rejects_missing_or_invalid_consumed_ids(tmp_path, identify, source, identifier):
    row = identify["arrays"][0]
    if identifier is None:
        del row["id"]
    else:
        row["id"] = identifier
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    with pytest.raises(ContractError) as error:
        parse_identify(path, source)
    assert str(path) in str(error.value)


@pytest.mark.parametrize("identifier", [None, "", 42, False, []])
def test_casandra_rejects_missing_or_invalid_consumed_ids(casandra, source, identifier):
    directory, fasta = casandra
    filename = "cassettes.jsonl"
    row = json.loads((directory / filename).read_text(encoding="utf-8"))
    if identifier is None:
        del row["cassette_id"]
    else:
        row["cassette_id"] = identifier
    rewrite_checked_record(directory, filename, json.dumps(row) + "\n")
    with pytest.raises(ContractError) as error:
        parse_casandra(directory, fasta, source)
    assert f"{directory / filename}:1:" in str(error.value)


@pytest.mark.parametrize("key", ["genes", "cas_cassettes"])
@pytest.mark.parametrize("count", [True, False, -1, 1.0, 1.5, "1"])
def test_casandra_consumed_counts_are_nonnegative_integers(casandra, source, key, count):
    directory, fasta = casandra
    update_run(directory, **{key: count})
    with pytest.raises(ContractError, match="nonnegative integer"):
        parse_casandra(directory, fasta, source)


def test_casandra_accepts_zero_counts(casandra, source):
    directory, fasta = casandra
    for filename in ("proteins.jsonl", "cassettes.jsonl"):
        rewrite_checked_record(directory, filename, "\n")
    update_run(directory, genes=0, cas_cassettes=0)
    result = parse_casandra(directory, fasta, source)
    assert result["genes"] == [] and result["cassettes"] == []


@pytest.mark.parametrize("repeat", ["", "ACGT!", None, 42])
def test_identify_malformed_repeat_uses_path_scoped_contract_error(
    tmp_path, identify, source, repeat
):
    identify["arrays"][0]["consensus_repeat"] = repeat
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    with pytest.raises(ContractError) as error:
        parse_identify(path, source)
    assert str(path) in str(error.value)


def test_identify_preserves_valid_unaccepted_arrays_and_additive_fields(tmp_path, identify, source):
    alternative = copy.deepcopy(identify["arrays"][0])
    alternative.update(category="Alternative", future_optional_field={"new": True})
    identify["arrays"].append(alternative)
    identify["arrays"][0]["consensus_repeat"] = "acgu"
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    arrays, preserved = parse_identify(path, source)
    assert len(arrays) == 1 and arrays[0].repeat == "ACGT"
    assert preserved == identify


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("source_interval", {"start": -10, "end": 99999}),
        ("source_interval", {"start": 1, "end": 99999}),
        ("source_interval", {"start": True, "end": 4}),
        ("source_sequence", "AAAA"),
        ("validation", {"status": "invalid"}),
        ("certainty_score", -1),
        ("certainty_score", 2),
        ("certainty_score", "0.9"),
        ("orientation", {"strand": "unknown"}),
        ("consensus_repeat", "!"),
        ("consensus_repeat", ""),
        ("id", []),
    ],
)
def test_identify_validates_unaccepted_array_rows(tmp_path, identify, source, field, value):
    alternative = copy.deepcopy(identify["arrays"][0])
    alternative.update(category="Alternative")
    alternative[field] = value
    identify["arrays"].append(alternative)
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    with pytest.raises(ContractError) as error:
        parse_identify(path, source)
    assert str(path) in str(error.value)


@pytest.mark.parametrize("coordinate", [True, 1.0, "1", None])
def test_identify_invalid_coordinate_types_remain_contract_errors(
    tmp_path, identify, source, coordinate
):
    identify["arrays"][0]["source_interval"]["start"] = coordinate
    path = tmp_path / "report.json"
    path.write_text(json.dumps(identify), encoding="utf-8")
    with pytest.raises(ContractError, match="coordinates must be integers"):
        parse_identify(path, source)


@pytest.mark.parametrize("coordinate", [True, 1.0, "1", None])
def test_casandra_invalid_coordinate_types_remain_contract_errors(casandra, source, coordinate):
    directory, fasta = casandra
    filename = "proteins.jsonl"
    row = json.loads((directory / filename).read_text(encoding="utf-8"))
    row["start_1based"] = coordinate
    rewrite_checked_record(directory, filename, json.dumps(row))
    with pytest.raises(ContractError, match="coordinates must be integers"):
        parse_casandra(directory, fasta, source)


def test_casandra_checksum_failure_precedes_json_validation(casandra, source):
    directory, fasta = casandra
    (directory / "proteins.jsonl").write_text("{broken", encoding="utf-8")
    with pytest.raises(ContractError, match="checksum mismatch"):
        parse_casandra(directory, fasta, source)


def test_casandra_input_digest_requires_matching_kind(casandra, source):
    directory, fasta = casandra
    digest = sha256_file(fasta)
    update_run(
        directory,
        inputs=[
            {"kind": "protein_fasta", "sha256": digest},
            {"kind": "genome_fasta", "sha256": "x"},
        ],
    )
    with pytest.raises(ContractError, match="input digest"):
        parse_casandra(directory, fasta, source)
    update_run(
        directory,
        inputs=[
            {"kind": "protein_fasta", "sha256": "x"},
            {"kind": "genome_fasta", "sha256": digest},
        ],
    )
    assert parse_casandra(directory, fasta, source)["status"] == "computed"


@pytest.fixture
def short_candidate():
    array = Array("array", Interval(250, 300), "T" * 20, "+", "Possible", 0.8)
    anti = AntiRepeat(array, array.repeat, "+", Interval(100, 120), "+", 0.9, 0.5, True)
    return Candidate("c", anti, Interval(95, 180), "A" * 85, Interval(5, 25))


def short_duplex(energy=-4, paired=20):
    return {
        "target_interval": Interval(5, 25).public(),
        "energy_kcal_mol": energy,
        "target_paired_positions": list(range(6, 6 + paired)),
        "repeat_paired_positions": list(range(1, 1 + paired)),
    }


@pytest.mark.parametrize("energy, accepted", [(-0.0001, True), (0, False), (0.0001, False)])
def test_short_admission_preserves_strict_energy_boundary(short_candidate, energy, accepted):
    model = Hit("c", Interval(25, 60), "+", 5, 0.01)
    result = short_rescue_admission(short_candidate, [model], [], short_duplex(energy), Settings())
    assert result["accepted"] is accepted
    assert result["rejection_reasons"] == ([] if accepted else ["nonnegative_energy"])


def test_short_admission_preserves_literal_thresholds_and_rejection_order(short_candidate):
    result = short_rescue_admission(
        replace(short_candidate, sequence="C" * 85), [], [], short_duplex(0, 9), Settings()
    )
    assert result["rejection_reasons"] == [
        "nonnegative_energy",
        "noncanonical_or_ambiguous_pair",
        "fewer_than_10_anchored_paired_bases",
        "less_than_half_anchor_paired",
        "no_independent_regional_support",
    ]
    # Compare serialized output to the literal pre-change schema, retaining key
    # order and integer/float representation as well as numeric values.
    assert json.dumps(result["thresholds"]) == json.dumps(
        {
            "min_query_span": 12,
            "min_identity": 0.8,
            "min_anchored_paired_bases": 10,
            "min_anchored_paired_fraction": 0.5,
            "energy_strictly_below": 0,
            "cm_evalue_max": 0.01,
            "cm_duplex_overlap_strictly_below": 10,
            "terminator_evalue_max": 0.001,
        }
    )
