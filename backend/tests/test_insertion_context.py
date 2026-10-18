import csv
import io
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import pytest
from crispr_web.execution.validation_insertion_context import (
    validate_insertion_context_pipeline,
    MODELS,
)
from crispr_web.execution.validation_insertion_evidence import document, sha, FIELDS
from crispr_web.summaries.insertion_context import insertion_context_summary
from crispr_web.publication.insertion_context import write_insertion_context_exports

SOURCE = b">f000001.s000001\nAAAACCCCTTTT\n"


def case(tmp_path):
    root = tmp_path / "result"
    root.mkdir()
    inputs = tmp_path / "input"
    inputs.mkdir()
    (inputs / "record.fasta").write_bytes(SOURCE)
    (root / "insertion-input").mkdir()
    (root / "insertion-input/source.fasta").write_bytes(SOURCE)
    (root / "context-input").mkdir()
    (root / "context-input/record-000001.fasta").write_bytes(SOURCE)
    shutil.copytree(Path(__file__).parent / "fixtures/insertion_cas_context", root / "casandra/record-000001")
    (root / "isescan").mkdir()
    stream = io.StringIO()
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(FIELDS)
    writer.writerow(
        [
            "seq000001",
            "IS1",
            "IS1_1",
            "1",
            "12",
            "12",
            "2",
            "1",
            "4",
            "9",
            "12",
            "8",
            "4",
            "4",
            "0",
            "0",
            "0",
            "",
            "1",
            "1e-20",
            "-1",
            "p",
            "0",
            "AAAA:AAAA",
        ]
    )
    native = stream.getvalue().encode()
    (root / "isescan/native.tsv").write_bytes(native)
    report = document(SOURCE, native, {"threads": 1, "include_partial": True}, MODELS)
    (root / "isescan/isescan.json").write_text(json.dumps(report))
    integration = {
        "workflow": "insertion_context",
        "stages": {
            "casandra": {"status": "completed", "gene_count": 0, "cas_gene_count": 0, "cassette_count": 0},
            "isescan": {
                "status": "completed",
                "counts": {
                    "sources": 1,
                    "elements": 1,
                    "complete": 0,
                    "partial": 1,
                    "with_terminal_inverted_repeats": 1,
                },
            }
        },
    }
    args = (
        SimpleNamespace(molecule="DNA"),
        root,
        {"f000001.s000001": (12, sha("AAAACCCCTTTT"))},
        integration,
        inputs,
    )
    return args, report


def test_independent_validation_and_safe_exports(tmp_path):
    args, report = case(tmp_path)
    assert validate_insertion_context_pipeline(*args)
    summary = {"insertion_context": insertion_context_summary(args[1], max_rows=100)}
    assert summary["insertion_context"]["elements"][0]["evalue_this_copy"] is None
    write_insertion_context_exports(args[1], summary)
    exported = json.loads((args[1] / "insertion-context.json").read_text())
    assert exported["elements"][0]["target_site_duplication_status"] == "not_assessed"
    assert (args[1] / "insertion-elements.tsv").is_file()


@pytest.mark.parametrize("change", ["source", "native", "report", "symlink", "counts"])
def test_rejects_mismatched_scientific_evidence(tmp_path, change):
    args, report = case(tmp_path)
    root = args[1]
    if change == "source":
        (args[-1] / "record.fasta").write_bytes(SOURCE.replace(b"CCCC", b"GGGG"))
    if change == "native":
        (root / "isescan/native.tsv").write_bytes(
            (root / "isescan/native.tsv")
            .read_bytes()
            .replace(b"AAAA:AAAA", b"AAAT:AAAA")
        )
    if change == "report":
        report["elements"][0]["transposase_orf"] = {
            "start": 1,
            "end": 12,
            "strand": "+",
        }
        (root / "isescan/isescan.json").write_text(json.dumps(report))
    if change == "symlink":
        target = root / "isescan/native.tsv"
        data = target.read_bytes()
        target.unlink()
        outside = tmp_path / "outside.tsv"
        outside.write_bytes(data)
        target.symlink_to(outside)
    if change == "counts":
        args[3]["stages"]["isescan"]["counts"]["elements"] = 0
    assert not validate_insertion_context_pipeline(*args)
