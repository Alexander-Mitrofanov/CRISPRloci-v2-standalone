import copy
import gzip
import json

import pytest

from crisprspacer.inputs import from_identify, validate_inventory
from crisprspacer.io import InputError, metadata, read_fasta, read_json
from conftest import record, report


def test_deletion_preserves_native_ordinal_gap(tmp_path):
    source = record("ACACGATTAC")
    data = report(source, ((1, 2), (3, 4), (9, 10)))
    path = tmp_path / "report.json"
    path.write_text(json.dumps(data))
    inventory = from_identify([source], [path])
    observed = next(q for q in inventory["queries"] if q["kind"] == "observed")
    deletion = next(q for q in inventory["queries"] if q["kind"] == "deletion")
    assert observed["ordinal"] == 2 and observed["sequence"] == "GATT"
    assert deletion["ordinal"] == 1 and deletion["sequence"] is None and deletion["boundary"] == 2


@pytest.mark.parametrize("mutation", [
    lambda r: r["source"].update(sha256="0"*64),
    lambda r: r["coordinate_system"].update(indexing=0),
    lambda r: r["arrays"][0]["spacers"][0].update(sequence="AAAA"),
    lambda r: r["arrays"][0]["spacers"][0]["source_interval"].update(start=6),
    lambda r: r["arrays"][0].update(spacer_count=2),
    lambda r: r["arrays"][0].update(source_sequence="AAAA"),
    lambda r: r["arrays"][0]["repeats"][0].update(sequence="AA"),
    lambda r: r["arrays"][0]["orientation"].update(strand="unknown"),
    lambda r: r["schema"].update(version="2.0.0"),
])
def test_report_claims_are_not_trusted(native_fixture, mutation):
    source, data, path, _ = native_fixture
    mutation(data)
    path.write_text(json.dumps(data))
    with pytest.raises(InputError):
        from_identify([source], [path])


def test_duplicate_source_report_rejected(native_fixture):
    source, _, path, _ = native_fixture
    with pytest.raises(InputError, match="duplicate"):
        from_identify([source], [path, path])


def test_native_category_and_unknown_orientation_preserved(native_fixture):
    source, _, path, _ = native_fixture
    inventory = from_identify([source], [path])
    assert inventory["queries"][0]["category"] == "Bona-fide"
    assert inventory["queries"][0]["source_orientation"] is None


def test_inventory_coordinate_tampering_rejected(native_fixture):
    source, _, path, _ = native_fixture
    inventory = from_identify([source], [path])
    tampered = copy.deepcopy(inventory)
    tampered["queries"][0]["start"] += 1
    with pytest.raises(InputError):
        validate_inventory(tampered)


@pytest.mark.parametrize("content", [">x\nACG\n>x\nTTT\n", "ACGT\n", ">x\n", ">x\nACGU\n", ">x\nAC GT\n", ">x\nA-C\n"])
def test_invalid_fasta(content, tmp_path):
    path = tmp_path / "bad.fasta"
    path.write_text(content)
    with pytest.raises(InputError):
        read_fasta(path)


def test_gzip_fasta_normalization(tmp_path):
    path = tmp_path / "sequences.fa.gz"
    path.write_bytes(gzip.compress(b">one description\nacgtn\n>two\nTTAA\n"))
    result = read_fasta(path)
    assert [r["id"] for r in result] == ["one", "two"]
    assert result[0]["sequence"] == "ACGTN"


def test_symlink_input_is_rejected(tmp_path):
    actual, link = tmp_path / "actual", tmp_path / "link"
    actual.write_text(">x\nACG\n")
    link.symlink_to(actual)
    with pytest.raises(InputError, match="Symlinks"):
        read_fasta(link)


@pytest.mark.parametrize("text", ['{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":1e9999}',
                                  '{"a":'+'9'*5000+'}', '{'])
def test_json_is_strict(tmp_path, text):
    path = tmp_path / "input.json"
    path.write_text(text)
    with pytest.raises(InputError):
        read_json(path)


def test_unknown_metadata_record_is_rejected(tmp_path):
    path = tmp_path / "meta.json"
    path.write_text(json.dumps({"schema_version": "1.0.0", "records": {"missing": {"assembly_id": "A"}}}))
    with pytest.raises(InputError):
        metadata([record("ACGT")], path)


def test_fasta_hash_comes_from_parsed_snapshot(tmp_path):
    import hashlib
    from crisprspacer.io import read_fasta_snapshot
    path = tmp_path / "input.fa"
    data = b">record\nacgt\n"
    path.write_bytes(data)
    rows, raw_hash = read_fasta_snapshot(path)
    assert raw_hash == hashlib.sha256(data).hexdigest()
    assert rows[0]["sequence"] == "ACGT"
