from copy import deepcopy

import pytest

from crisprmap.io import InputError, canonical, digest
from crisprmap.reference import build_reference, load_reference, validate_reference


def example():
    return build_reference([{"id": "a", "sequence": "aauu", "description": "historical repeat"}],
                           {"a": {"family": "F1", "motif": None, "accession": "NC_000001"}},
                           {"namespace": "example", "source_url": "https://example.org/repeats",
                            "source_sha256": "0" * 64, "license": "operator must assess"})


def test_reference_roundtrip_and_digest(tmp_path):
    reference = example()
    assert reference["entries"][0]["sequence"] == "AATT"
    assert reference["entries"][0]["sequence_sha256"] == digest("AATT")
    assert validate_reference(reference) == reference
    path = tmp_path / "reference.json"
    path.write_bytes(canonical(reference))
    loaded, raw_sha = load_reference(path)
    assert loaded == reference
    assert raw_sha == digest(path.read_bytes())


@pytest.mark.parametrize("field,value", [("sequence", "AAAA"), ("family", "F2"), ("id", "b"),
                                         ("sequence_sha256", "1" * 64)])
def test_modified_reference_binding_rejected(field, value):
    reference = example()
    reference["entries"][0][field] = value
    with pytest.raises(InputError):
        validate_reference(reference)


@pytest.mark.parametrize("field,value", [("id", []), ("sequence", []), ("family", True),
                                         ("motif", "bad\nlabel"), ("description", None)])
def test_malformed_reference_fields_raise_contract_error(field, value):
    reference = example()
    reference["entries"][0][field] = value
    with pytest.raises(InputError):
        validate_reference(reference)


def test_missing_labels_are_not_invented():
    ref = build_reference([{"id": "a", "sequence": "ACGT"}], {}, example()["provenance"])
    assert ref["entries"][0]["family"] is None
    assert ref["entries"][0]["motif"] is None


@pytest.mark.parametrize("mutation", ["duplicate", "extra", "schema", "lowercase", "digest"])
def test_reference_contract_mutations(mutation):
    ref = deepcopy(example())
    if mutation == "duplicate":
        ref["entries"].append(deepcopy(ref["entries"][0]))
    elif mutation == "extra":
        ref["unexpected"] = True
    elif mutation == "schema":
        ref["schema"]["version"] = "2.0.0"
    elif mutation == "lowercase":
        ref["entries"][0]["sequence"] = "aatt"
    else:
        ref["content_sha256"] = "1" * 64
    with pytest.raises(InputError):
        validate_reference(ref)


def test_unknown_metadata_and_unprovenanced_reference_rejected():
    with pytest.raises(InputError):
        build_reference([{"id": "a", "sequence": "ACGT"}], {"unknown": {}}, example()["provenance"])
    with pytest.raises(InputError):
        build_reference([{"id": "a", "sequence": "ACGT"}], {"a": {"subtype": "I"}}, example()["provenance"])
    for key in ("namespace", "source_url", "source_sha256"):
        provenance = dict(example()["provenance"])
        provenance.pop(key)
        with pytest.raises(InputError):
            build_reference([{"id": "a", "sequence": "ACGT"}], {}, provenance)
