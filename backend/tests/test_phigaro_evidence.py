"""Native positive fixture, exact-match replay and hostile evidence checks."""

from copy import deepcopy
import json
from pathlib import Path
import pytest
from crispr_web.execution.validation_prophage_evidence import (
    validate_report as backend_report,
)
from crispr_web.execution.validation_source_spacer_evidence import (
    validate_search as backend_search,
)
from crispr_integration.prophage_evidence import validate_report as integration_report
from crispr_integration.source_spacer_evidence import (
    validate_search as integration_search,
)

FIXTURES = Path(__file__).parent / "fixtures"


def native_case():
    folder = FIXTURES / "phigaro"
    payload = (folder / "input.fasta").read_bytes()
    lines = payload.decode().splitlines()
    return (
        json.loads((folder / "report.json").read_text()),
        {lines[0][1:]: "".join(lines[1:])},
        payload,
        (folder / "proteins.faa").read_bytes(),
        (folder / "hits.tsv").read_bytes(),
    )


@pytest.mark.parametrize("validate", [backend_report, integration_report])
def test_native_phigaro_positive(validate):
    assert validate(*native_case()) == {
        "sources": 1,
        "assessed_sources": 1,
        "genes": 61,
        "regions": 1,
    }


@pytest.mark.parametrize(
    "mutation",
    [
        "protein",
        "end",
        "omitted_region",
        "false_gc",
        "model",
        "nan_score",
        "native_hit",
        "source",
        "terminal_gene",
    ],
)
@pytest.mark.parametrize("validate", [backend_report, integration_report])
def test_phigaro_evidence_tampering(validate, mutation):
    report, records, payload, proteins, hits = native_case()
    report = deepcopy(report)
    if mutation == "protein":
        report["genes"][0]["translation"] = "AAAA"
    if mutation == "end":
        report["regions"][0]["end"] -= 1
    if mutation == "omitted_region":
        report["regions"] = []
    if mutation == "false_gc":
        report["genes"][0]["source_gc"] = False
    if mutation == "model":
        report["reference"]["sha256"] = "0" * 64
    if mutation == "nan_score":
        report["regions"][0]["peak_score"] = float("nan")
    if mutation == "native_hit":
        hits = b""
    if mutation == "source":
        records = {key: "C" + value[1:] for key, value in records.items()}
    if mutation == "terminal_gene":
        report["regions"][0]["gene_ids"].pop()
    with pytest.raises((ValueError, KeyError)):
        validate(report, records, payload, proteins, hits)


def search_case():
    return json.loads((FIXTURES / "prophage-spacer/fixture.json").read_text())


@pytest.mark.parametrize("validate", [backend_search, integration_search])
def test_complete_exact_search_and_interval_join(validate):
    case = search_case()
    matches, counts = validate(
        case["search"], case["reports"], case["records"], case["regions"]
    )
    assert counts["excluded_array_matches"] == 1 and counts["candidate_matches"] == 2
    assert {m["strand"] for m in matches} == {"+", "-"}
    assert sum(bool(m["region_ids"]) for m in matches) == 1
    assert all(m["match_class"] == "sequence_match" for m in matches)


@pytest.mark.parametrize(
    "mutation",
    [
        "omitted",
        "exclude",
        "query",
        "native_source",
        "native_spacer",
        "coordinate",
        "policy",
    ],
)
@pytest.mark.parametrize("validate", [backend_search, integration_search])
def test_spacer_match_tampering(validate, mutation):
    case = search_case()
    report = case["search"]
    if mutation == "omitted":
        report["matches"].pop()
    if mutation == "exclude":
        report["matches"][0]["excluded"] = False
    if mutation == "query":
        report["queries"][0]["sequence"] = "A" * 16
    if mutation == "native_source":
        case["reports"][0]["source"]["sha256"] = "0" * 64
    if mutation == "native_spacer":
        case["reports"][0]["arrays"][0]["spacers"][0]["source_interval"]["start"] += 1
    if mutation == "coordinate":
        report["matches"][1]["target_start"] += 1
    if mutation == "policy":
        report["policy"]["substitutions"] = 1
    with pytest.raises((ValueError, KeyError)):
        validate(report, case["reports"], case["records"], case["regions"])
