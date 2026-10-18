"""Historical identifier and original-result semantics are independent of core."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location("fetch_legacy", Path(__file__).parents[1] / "scripts/fetch_legacy.py")
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


def test_historical_unknowns_and_accession_are_not_inferred():
    records = [
        {"id": "Taxon_A_G_1_M6_F8", "description": "NC_123456_1_4_20"},
        {"id": "Other_A_P_2_X_F10", "description": "CP001398_1_5_22"},
        {"id": "Unknown_X_X", "description": "undocumented"},
    ]
    rows = legacy.metadata_from_legacy(records)
    assert rows[records[0]["id"]] == {"family": "F8", "motif": "M6", "accession": "NC_123456"}
    assert rows[records[1]["id"]] == {"family": "F10", "motif": None, "accession": "CP001398"}
    assert rows[records[2]["id"]] == {"family": None, "motif": None, "accession": None}


def test_invalid_historical_identifier_is_rejected():
    with pytest.raises(ValueError, match="Unrecognized"):
        legacy.metadata_from_legacy([{"id": "no_labels"}])


def test_original_results_preserve_missing_and_original_row_order():
    rows = legacy.parse_original_results("28\t-\tE\tID1 UserId=example\tGAUC\n-\t4\t-\tID2 UserId=other CRISPRmap=Taxon_X_F4\tGTTA\n")
    assert rows[0]["motif"] == "M28"
    assert rows[0]["family"] is None
    assert rows[0]["sequence"] == "GATC"
    assert rows[1]["superclass"] is None
    assert rows[1]["original_reference_id"] == "Taxon_X_F4"


@pytest.mark.parametrize("text", ["1\t2\tC\tID2 UserId=q\tACGT", "1\t2\tC\tID1 UserId=q\tACNT", "x\t2\tC\tID1 UserId=q\tACGT"])
def test_original_results_reject_unexpected_shapes(text):
    with pytest.raises(ValueError):
        legacy.parse_original_results(text)


def test_pin_checks_both_length_and_digest():
    import hashlib
    valid = b"reference"
    pin = {"bytes": len(valid), "sha256": hashlib.sha256(valid).hexdigest()}
    assert legacy.checked_bytes(valid, pin) == valid
    with pytest.raises(ValueError, match="pinned"):
        legacy.checked_bytes(b"Reference", pin)
