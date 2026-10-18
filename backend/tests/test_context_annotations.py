"""Cas-context completeness, independent binding and bounded coordinate joins."""
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from crispr_web.execution.validation_context_annotations import validate_cas_context
from crispr_web.summaries.context_annotations import cas_neighbors


@pytest.fixture
def native(tmp_path):
    fixture = Path(__file__).parent / "fixtures/defense_context/input.fasta"
    payload = fixture.read_bytes()
    sequence = "".join(payload.decode().splitlines()[1:])
    (tmp_path / "context-input").mkdir()
    (tmp_path / "context-input/record-000001.fasta").write_bytes(payload)
    shutil.copytree(Path(__file__).parent / "fixtures/defense_cas_context", tmp_path / "casandra/record-000001")
    return tmp_path, {"f000001.s000001": sequence}, {"status": "completed", "gene_count": 42, "cas_gene_count": 0, "cassette_count": 0}


def test_actual_native_context(native):
    validate_cas_context(*native)


@pytest.mark.parametrize("mutation", ["source", "extra_directory", "extra_input", "sequence", "digest", "counts"])
def test_cas_context_mutation_rejected(native, mutation):
    root, records, stage = native
    folder = root / "casandra/record-000001"
    if mutation == "source":
        (root / "context-input/record-000001.fasta").write_text(">other\nACGT\n")
    elif mutation == "extra_directory":
        (root / "casandra/record-000002").mkdir()
    elif mutation == "extra_input":
        (root / "context-input/extra.fasta").write_text(">extra\nACGT\n")
    elif mutation == "sequence":
        path = folder / "proteins.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        rows[0]["nucleotide_sequence"] = "A" + rows[0]["nucleotide_sequence"][1:]
        if rows[0]["nucleotide_sequence"] == json.loads(path.read_text().splitlines()[0])["nucleotide_sequence"]:
            rows[0]["nucleotide_sequence"] = "C" + rows[0]["nucleotide_sequence"][1:]
        path.write_text("".join(json.dumps(row)+"\n" for row in rows))
        doc = json.loads((folder / "manifest.json").read_text())
        doc["files"]["proteins.jsonl"]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        (folder / "manifest.json").write_text(json.dumps(doc))
    elif mutation == "digest":
        with (folder / "proteins.jsonl").open("a") as handle:
            handle.write("\n")
    else:
        stage["gene_count"] -= 1
    with pytest.raises(ValueError):
        validate_cas_context(root, records, stage)


def test_context_joins_ties_overlap_adjacency_and_source_scope(tmp_path):
    folder = tmp_path / "casandra/record-000001"
    folder.mkdir(parents=True)
    cassettes = [
        {"cassette_id": "left", "contig_id": "a", "start_1based": 1, "end_1based_inclusive": 10},
        {"cassette_id": "right", "contig_id": "a", "start_1based": 21, "end_1based_inclusive": 30},
        {"cassette_id": "other", "contig_id": "b", "start_1based": 10, "end_1based_inclusive": 20},
    ]
    (folder / "cassettes.jsonl").write_text("".join(json.dumps(row)+"\n" for row in cassettes))
    records = cas_neighbors(tmp_path, [{"id": "interval", "source_id": "a", "start": 11, "end": 20}])
    assert {row["cassette_id"] for row in records[0]["nearest_cassettes"]} == {"left", "right"}
    assert all(row["overlap_bases"] == 0 and row["intervening_bases"] == 0 for row in records[0]["nearest_cassettes"])
    records = cas_neighbors(tmp_path, [{"id": "interval", "source_id": "a", "start": 9, "end": 20}])
    assert records[0]["nearest_cassettes"][0]["overlap_bases"] == 2


def test_context_tie_preview_is_bounded_and_explicit(tmp_path):
    folder = tmp_path / "casandra/record-000001"
    folder.mkdir(parents=True)
    (folder / "cassettes.jsonl").write_text("".join(json.dumps({"cassette_id": str(index), "contig_id": "a", "start_1based": 1, "end_1based_inclusive": 10})+"\n" for index in range(100)))
    row = cas_neighbors(tmp_path, [{"id": "interval", "source_id": "a", "start": 5, "end": 10}])[0]
    assert row["nearest_cassette_count"] == 100 and row["nearest_cassettes_truncated"]
    assert len(row["nearest_cassettes"]) == 32


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("caller", ["Pyrodigal", "unsupported-caller"])
def test_caller_normalized_ambiguity_remains_source_bound(tmp_path, reverse, caller):
    fixture = Path(__file__).parent / "fixtures/defense_ambiguity"
    sequence = "".join((fixture / "input.fasta").read_text().splitlines()[1:])
    folder = tmp_path / "casandra/record-000001"
    shutil.copytree(fixture / "casandra", folder)
    rows = [json.loads(line) for line in (folder / "proteins.jsonl").read_text().splitlines()]
    if reverse:
        sequence = sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
        for row in rows:
            row["start_1based"], row["end_1based_inclusive"] = len(sequence)-row["end_1based_inclusive"]+1, len(sequence)-row["start_1based"]+1
            row["strand"] = "+" if row["strand"] == "-" else "-"
    for row in rows:
        row["caller"]["program"] = caller
    raw = "".join(json.dumps(row)+"\n" for row in rows).encode()
    (folder / "proteins.jsonl").write_bytes(raw)
    manifest = json.loads((folder / "manifest.json").read_text())
    manifest["files"]["proteins.jsonl"]["sha256"] = hashlib.sha256(raw).hexdigest()
    (folder / "manifest.json").write_text(json.dumps(manifest))
    payload = f">f000001.s000001\n{sequence}\n".encode()
    (tmp_path / "context-input").mkdir()
    (tmp_path / "context-input/record-000001.fasta").write_bytes(payload)
    run = json.loads((folder / "run.json").read_text())
    for item in run["inputs"]:
        if item["kind"] == "genome_fasta":
            item["sha256"] = hashlib.sha256(payload).hexdigest()
    (folder / "run.json").write_text(json.dumps(run))
    args = (tmp_path, {"f000001.s000001": sequence}, {"status": "completed", "gene_count": run["genes"], "cas_gene_count": run["cas_proteins"], "cassette_count": run["cas_cassettes"]})
    if caller == "Pyrodigal":
        validate_cas_context(*args)
    else:
        with pytest.raises(ValueError, match="gene sequence"):
            validate_cas_context(*args)
