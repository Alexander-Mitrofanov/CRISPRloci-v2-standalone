from pathlib import Path

import pytest

from casandra.runtime import parse_fasta


def test_protein_fasta_rejects_duplicate_ids(tmp_path: Path):
    path = tmp_path / "duplicate.faa"
    path.write_text(">same\nMKT\n>same\nAAA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate"):
        list(parse_fasta(path))


def test_protein_fasta_rejects_empty_and_invalid_records(tmp_path: Path):
    empty = tmp_path / "empty.faa"
    empty.write_text(">empty\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no sequence"):
        list(parse_fasta(empty))
    invalid = tmp_path / "invalid.faa"
    invalid.write_text(">bad\nMK?\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Invalid amino-acid"):
        list(parse_fasta(invalid))

