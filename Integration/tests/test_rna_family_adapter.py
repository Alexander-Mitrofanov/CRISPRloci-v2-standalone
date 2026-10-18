"""Independent adapter gate over a real native synthetic family fixture."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from crispr_integration.rna_family import records_from_input, validate_report, sha

FIXTURE = Path(__file__).resolve().parents[2] / "backend/tests/fixtures/rna_family"


def test_native_family_contract():
    records = records_from_input(FIXTURE / "input.fasta", "DNA")
    doc = json.loads((FIXTURE / "report.json").read_text())
    validate_report(doc, records, "DNA", sha((FIXTURE / "input.fasta").read_bytes()))


@pytest.mark.parametrize("field", ["sources", "alignment", "pairs", "probabilities", "columns"])
def test_missing_evidence_rejected(field):
    records = records_from_input(FIXTURE / "input.fasta", "DNA")
    doc = deepcopy(json.loads((FIXTURE / "report.json").read_text()))
    doc[field].pop(0)
    with pytest.raises(ValueError):
        validate_report(doc, records, "DNA", sha((FIXTURE / "input.fasta").read_bytes()))


def test_family_rejects_ambiguous_and_identical_sequences(tmp_path):
    path = tmp_path / "family.fa"
    path.write_text(">a\n" + "N"*40 + "\n>b\n" + "A"*40 + "\n")
    with pytest.raises(ValueError):
        records_from_input(path, "DNA")
    path.write_text(">a\n" + "A"*40 + "\n>b\n" + "A"*40 + "\n")
    with pytest.raises(ValueError):
        records_from_input(path, "RNA")
