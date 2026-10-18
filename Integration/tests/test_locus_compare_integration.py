import copy
import json
from pathlib import Path
import pytest
from crispr_integration.locus_compare_validation import validate_report, sha, translate
from crispr_integration.locus_compare import feature_document

FIXTURE = Path(__file__).resolve().parents[2] / "tests/fixtures/locus_compare"


def fixture():
    raw = (FIXTURE / "input.fasta").read_bytes()
    features = (FIXTURE / "features.json").read_bytes()
    records = {row.splitlines()[0]: "".join(row.splitlines()[1:]) for row in raw.decode().split(">")[1:]}
    return (
        json.loads((FIXTURE / "comparison.json").read_text()),
        json.loads(features),
        records,
        sha(raw),
        sha(features),
        0.3,
    )


def test_native_fixture_source_translation_and_alignment():
    assert validate_report(*fixture()) == {"loci": 2, "genes": 2, "links": 1, "evaluated_protein_pairs": 1}


@pytest.mark.parametrize(
    "key,value",
    [("identity", 0.9), ("similarity", 0.4), ("score", -100), ("aligned_query", "A" * 26), ("target_gene", "unknown")],
)
def test_reject_alignment_tampering(key, value):
    args = list(fixture())
    args[0]["links"][0][key] = value
    with pytest.raises(ValueError):
        validate_report(*args)


def test_reject_rewritten_gene_translation():
    args = list(fixture())
    args[0]["loci"][0]["genes"][0]["translation"] = "A" * 26
    args[1]["loci"][0]["genes"][0]["translation"] = "A" * 26
    with pytest.raises(ValueError, match="translation"):
        validate_report(*args)


def test_reject_duplicate_gene_links():
    args = list(fixture())
    args[0]["links"] *= 2
    with pytest.raises(ValueError):
        validate_report(*args)


def test_genetic_codes_and_ambiguous_residues():
    assert translate("ATGTGATAA", 4) == "MW"
    assert translate("RAY", 11) == "B"
    assert translate("SAR", 11) == "Z"
