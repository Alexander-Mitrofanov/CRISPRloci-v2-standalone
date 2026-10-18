import json
from pathlib import Path
import pytest
from crispr_web.execution.validation_locus_compare import (
    validate_report,
    sha,
    translate,
)
from crispr_web.summaries.locus_compare import locus_compare_summary
from crispr_web.publication.locus_compare import write_locus_compare_exports

FIXTURE = Path(__file__).resolve().parents[2] / "tests/fixtures/locus_compare"


def fixture():
    raw = (FIXTURE / "input.fasta").read_bytes()
    features = (FIXTURE / "features.json").read_bytes()
    records = {
        row.splitlines()[0]: "".join(row.splitlines()[1:])
        for row in raw.decode().split(">")[1:]
    }
    return (
        json.loads((FIXTURE / "comparison.json").read_text()),
        json.loads(features),
        records,
        sha(raw),
        sha(features),
        0.3,
    )


def test_independent_native_report_gate():
    assert validate_report(*fixture()) == {
        "loci": 2,
        "genes": 2,
        "links": 1,
        "evaluated_protein_pairs": 1,
    }


@pytest.mark.parametrize(
    "key,value",
    [
        ("identity", 0.9),
        ("similarity", 0.4),
        ("score", -100),
        ("aligned_query", "A" * 26),
        ("target_gene", "unknown"),
    ],
)
def test_reject_forged_alignment(key, value):
    args = list(fixture())
    args[0]["links"][0][key] = value
    with pytest.raises(ValueError):
        validate_report(*args)


def test_reject_forged_source_translation():
    args = list(fixture())
    args[0]["loci"][0]["genes"][0]["translation"] = "A" * 26
    args[1]["loci"][0]["genes"][0]["translation"] = "A" * 26
    with pytest.raises(ValueError):
        validate_report(*args)


def test_bounded_preview_and_recreated_exports(tmp_path):
    native = tmp_path / "locus-compare"
    native.mkdir()
    (native / "comparison.json").write_bytes((FIXTURE / "comparison.json").read_bytes())
    summary = locus_compare_summary(tmp_path, max_rows=0)
    assert summary["links"] == [] and summary["links_truncated"]
    assert "translation" not in summary["loci"][0]["genes"][0]
    write_locus_compare_exports(tmp_path, {"locus_compare": summary})
    report = json.loads((tmp_path / "locus-comparison.json").read_text())
    assert len(report["links"]) == 1
    assert (tmp_path / "locus-proteins.faa").read_text().count(">") == 2
    assert "aligned_query" in (tmp_path / "locus-alignments.tsv").read_text()


def test_independent_genetic_code():
    assert translate("ATGTGATAA", 4) == "MW"
    assert translate("RAY", 11) == "B"
