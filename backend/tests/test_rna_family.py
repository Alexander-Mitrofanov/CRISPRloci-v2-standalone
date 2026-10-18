"""Native synthetic family fixture and independent source/publication gates."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from crispr_web.execution.validation_rna_family import validate_rna_family_pipeline
from crispr_web.summaries.rna_family import rna_family_summary
from crispr_web.publication.rna_family import write_rna_family_exports

FIXTURE = Path(__file__).parent / "fixtures/rna_family"


@pytest.fixture
def evidence(tmp_path):
    root = tmp_path / "output"
    shutil.copytree(FIXTURE, root / "rna-family")
    source = tmp_path / "input"
    source.mkdir()
    records = {}
    for part in (FIXTURE / "input.fasta").read_text().split(">")[1:]:
        key, *parts = part.splitlines()
        records[key] = "".join(parts)
    for index, (key, sequence) in enumerate(records.items()):
        (source / f"{index:06d}.fa").write_text(f">{key}\n{sequence}\n")
    doc = json.loads((root / "rna-family/report.json").read_text())
    integration = {"workflow": "rna_family", "stages": {"rnaalifold": {"status": "completed", "counts": doc["counts"]}}}
    expected = {key: (len(seq), hashlib.sha256(seq.encode()).hexdigest()) for key, seq in records.items()}
    return root, source, doc, integration, expected


def valid(data):
    root, source, _, integration, expected = data
    return validate_rna_family_pipeline(SimpleNamespace(molecule="DNA"), root, expected, integration, source)


def test_native_fixture_valid_and_descriptive_covariation(evidence):
    assert valid(evidence)
    assert evidence[2]["pairs"][0]["compensatory_sequence_pairs"] == 2
    assert evidence[2]["counts"]["unique_sequences"] == 3


@pytest.mark.parametrize("mutation", ["sequence", "members", "consensus", "structure", "probability", "mass", "covariation",
                                     "native_alignment", "native_fold", "symlink", "policy", "completion", "boolean"])
def test_modified_evidence_rejected(evidence, mutation):
    root, _, doc, integration, _ = evidence
    doc = deepcopy(doc)
    if mutation == "sequence":
        doc["alignment"][0]["sequence"] = "A" + doc["alignment"][0]["sequence"][1:]
    elif mutation == "members":
        doc["alignment"][0]["member_ids"].pop()
    elif mutation == "consensus":
        doc["fold"]["consensus"] = "A" + doc["fold"]["consensus"][1:]
    elif mutation == "structure":
        doc["fold"]["structure"] = ")" + doc["fold"]["structure"][1:]
    elif mutation == "probability":
        doc["probabilities"][0]["probability"] = -0.1
    elif mutation == "mass":
        doc["probabilities"][0]["probability"] = 1.0
        doc["probabilities"][1]["probability"] = 1.0
    elif mutation == "covariation":
        doc["pairs"][0]["compensatory_sequence_pairs"] += 1
    elif mutation == "native_alignment":
        (root / "rna-family/work/aligned.fasta").write_text(">family000001\nAAAA\n")
    elif mutation == "native_fold":
        p = root / "rna-family/work/fold.json"
        native = json.loads(p.read_text())
        native["mfe_score_kcal_mol"] = 0
        p.write_text(json.dumps(native))
    elif mutation == "symlink":
        p = root / "rna-family/work/fold.json"
        other = root.parent / "elsewhere.json"
        other.write_bytes(p.read_bytes())
        p.unlink()
        p.symlink_to(other)
    elif mutation == "policy":
        doc["policy"]["homology"] = "proven"
    elif mutation == "completion":
        integration["stages"]["rnaalifold"]["status"] = "failed"
    elif mutation == "boolean":
        doc["probabilities"][0]["left"] = True
    (root / "rna-family/report.json").write_text(json.dumps(doc))
    assert not valid(evidence)


def test_preview_bound_does_not_remove_complete_export(evidence):
    root = evidence[0]
    preview = rna_family_summary(root, max_rows=1)
    assert len(preview["pairs"]) == 1 and preview["pairs_truncated"]
    write_rna_family_exports(root, {"rna_family": preview})
    public = json.loads((root / "rna-family.json").read_text())
    assert len(public["pairs"]) == 12 and len(public["probabilities"]) == 14
    assert (root / "rna-family-consensus.dbn").exists()
    assert "work/" not in (root / "rna-family.json").read_text()


def test_absent_family_summary_is_unavailable(tmp_path):
    assert rna_family_summary(tmp_path, max_rows=10)["available"] is False
