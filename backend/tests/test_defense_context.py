"""Real-native fixture validation, hostile evidence and complete export checks."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from crispr_integration.defense_context import validate_report, MODEL
from crispr_web.execution.validation_defense_context import validate_defense_context_pipeline
from crispr_web.summaries.defense_context import defense_context_summary
from crispr_web.publication.defense_context import write_defense_context_exports

FIXTURE = Path(__file__).parent / "fixtures/defense_context"


@pytest.fixture
def evidence(tmp_path):
    root = tmp_path / "output"
    shutil.copytree(FIXTURE, root / "defense-context")
    source = tmp_path / "input"
    source.mkdir()
    shutil.copy(FIXTURE / "input.fasta", source / "source.fasta")
    seq = "".join((source / "source.fasta").read_text().splitlines()[1:])
    (root / "context-input").mkdir()
    shutil.copyfile(source / "source.fasta", root / "context-input/record-000001.fasta")
    shutil.copytree(Path(__file__).parent / "fixtures/defense_cas_context", root / "casandra/record-000001")
    doc = json.loads((root / "defense-context/report.json").read_text())
    integration = {"workflow": "defense_context", "stages": {"casandra": {"status": "completed", "gene_count": 42, "cas_gene_count": 0, "cassette_count": 0}, "defensefinder": {
        "status": "completed", "counts": doc["counts"], "model_id": MODEL["id"],
        "model_sha256": MODEL["manifest_sha256"]}}}
    expected = {"f000001.s000001": (len(seq), hashlib.sha256(seq.encode()).hexdigest())}
    return root, source, seq, doc, integration, expected


def valid(data):
    root, source, _, _, integration, expected = data
    return validate_defense_context_pipeline(SimpleNamespace(molecule="DNA"), root, expected, integration, source)


def test_real_fixture_and_independent_gates(evidence):
    _, _, seq, doc, _, _ = evidence
    validate_report(doc, {"f000001.s000001": seq}, doc["input"]["sha256"])
    assert valid(evidence)
    assert doc["counts"] == {"sources": 1, "genes": 42, "systems": 1, "components": 2, "hits": 5}
    assert doc["systems"][0]["subtype"] == "RM_Type_II"


@pytest.mark.parametrize("mutation", ["translation", "source", "order", "bounds", "model", "score", "membership", "native", "symlink"])
def test_tampered_evidence_fails(evidence, mutation):
    root, _, _, doc, _, _ = evidence
    doc = deepcopy(doc)
    if mutation == "translation":
        doc["genes"][0]["protein"] = "A" + doc["genes"][0]["protein"][1:]
    elif mutation == "source":
        doc["genes"][0]["source_id"] = "other"
    elif mutation == "order":
        doc["genes"][0]["position"] = True
    elif mutation == "bounds":
        doc["systems"][0]["start"] -= 1
    elif mutation == "model":
        doc["model"]["version"] = "3.1.0"
    elif mutation == "score":
        doc["hits"][0]["score"] = True
    elif mutation == "membership":
        doc["systems"][0]["gene_ids"].pop()
    elif mutation == "native":
        path = root / "defense-context/work/RM/best_solution.tsv"
        path.write_text("# No Systems found\n")
    elif mutation == "symlink":
        path = root / "defense-context/work/RM/best_solution.tsv"
        contents = path.read_text()
        other = root.parent / "outside.tsv"
        other.write_text(contents)
        path.unlink()
        path.symlink_to(other)
    (root / "defense-context/report.json").write_text(json.dumps(doc))
    assert not valid(evidence)


def test_preview_truncation_does_not_truncate_download(evidence):
    root = evidence[0]
    preview = defense_context_summary(root, max_rows=0)
    assert preview["systems"] == [] and preview["systems_truncated"]
    write_defense_context_exports(root, {"defense_context": preview})
    public = json.loads((root / "defense-context.json").read_text())
    assert len(public["systems"]) == 1 and len(public["hits"]) == 5 and len(public["genes"]) == 42
    assert "work/" not in (root / "defense-context.json").read_text()


@pytest.mark.parametrize("mutation", ["extra_base", "missing_stop", "false_initiator"])
def test_coherently_rehashed_cds_boundary_tampering_fails_both_gates(evidence, mutation):
    root, _, sequence, doc, _, _ = evidence
    doc = deepcopy(doc)
    gene = doc["genes"][0]  # Complete reverse-strand gene whose protein is unchanged by these edits.
    assert gene["strand"] == -1 and not gene["partial_begin"] and not gene["partial_end"]
    if mutation == "extra_base":
        gene["start"] -= 1
    elif mutation == "missing_stop":
        gene["start"] += 3
    else:
        gene["initiator_methionine"] = False
    cds = sequence[gene["start"]-1:gene["end"]].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    gene["cds_sha256"] = hashlib.sha256(cds.encode()).hexdigest()
    with pytest.raises(ValueError, match="boundary"):
        validate_report(doc, {"f000001.s000001": sequence}, doc["input"]["sha256"])
    (root / "defense-context/report.json").write_text(json.dumps(doc))
    assert not valid(evidence)


def test_real_ambiguous_cds_preserves_exact_source_hash_and_native_x(tmp_path):
    fixture = Path(__file__).parent / "fixtures/defense_ambiguity"
    root = tmp_path / "output"
    source = tmp_path / "input"
    source.mkdir()
    shutil.copy(fixture / "input.fasta", source / "source.fasta")
    shutil.copytree(fixture, root / "defense-context")
    (root / "context-input").mkdir()
    shutil.copy(source / "source.fasta", root / "context-input/record-000001.fasta")
    shutil.copytree(fixture / "casandra", root / "casandra/record-000001")
    sequence = "".join((source / "source.fasta").read_text().splitlines()[1:])
    doc = json.loads((root / "defense-context/report.json").read_text())
    run = json.loads((root / "casandra/record-000001/run.json").read_text())
    integration = {"workflow": "defense_context", "stages": {
        "casandra": {"status": "completed", "gene_count": run["genes"], "cas_gene_count": run["cas_proteins"], "cassette_count": run["cas_cassettes"]},
        "defensefinder": {"status": "completed", "counts": doc["counts"], "model_id": MODEL["id"], "model_sha256": MODEL["manifest_sha256"]}}}
    expected = {"f000001.s000001": (len(sequence), hashlib.sha256(sequence.encode()).hexdigest())}
    validate_report(doc, {"f000001.s000001": sequence}, doc["input"]["sha256"])
    assert valid((root, source, sequence, doc, integration, expected))
    affected = [g for g in doc["genes"] if g["start"] <= 301 and g["end"] >= 311]
    assert affected and all("X" in gene["protein"] for gene in affected)
    # A caller-normalized DNA digest cannot replace the exact submitted CDS hash.
    gene = affected[0]
    cds = sequence[gene["start"]-1:gene["end"]]
    if gene["strand"] == -1:
        cds = cds.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
    gene["cds_sha256"] = hashlib.sha256("".join(b if b in "ACGT" else "N" for b in cds).encode()).hexdigest()
    with pytest.raises(ValueError, match="translation"):
        validate_report(doc, {"f000001.s000001": sequence}, doc["input"]["sha256"])
    (root / "defense-context/report.json").write_text(json.dumps(doc))
    assert not valid((root, source, sequence, doc, integration, expected))


@pytest.mark.parametrize("table,column", [("best", "hit_pos"), ("best", "hit_seq_len"), ("hmm", 2), ("hmm", 3), ("fasta", None)])
def test_native_global_rank_and_protein_length_are_bound(evidence, table, column):
    import csv
    import io
    root = evidence[0]
    assert valid(evidence)
    if table == "fasta":
        path = root / "defense-context/work/proteins.faa"
        rows = path.read_text().splitlines(keepends=True)
        path.write_text("".join(rows[2:4] + rows[:2] + rows[4:]))
    elif table == "best":
        path = root / "defense-context/work/RM/best_solution.tsv"
        rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))
        rows[0][column] = str(int(rows[0][column]) + 1)
        with path.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    else:
        paths = (root / "defense-context/work").glob("*/hmmer_results/*.res_hmm_extract")
        for path in paths:
            lines = path.read_text().splitlines()
            indices = [i for i, row in enumerate(lines) if row and not row.startswith("#")]
            if indices:
                index = indices[0]
                cells = lines[index].split("\t")
                cells[column] = str(int(cells[column]) + 1)
                lines[index] = "\t".join(cells)
                path.write_text("\n".join(lines) + "\n")
                break
        else:
            raise AssertionError("positive fixture must contain native hits")
    assert not valid(evidence)
