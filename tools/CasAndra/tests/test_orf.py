"""Tests for reproducible prokaryotic protein extraction."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from casandra.orf import Contig, choose_gene_mode, read_nucleotide_fasta


class OrfTests(unittest.TestCase):
    def test_fasta_preserves_description_and_normalizes_rna(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fna"
            path.write_text(">ctg1 full description\nAUGCNN\n", encoding="utf-8")
            records = read_nucleotide_fasta(path)
            self.assertEqual(records[0].description, "ctg1 full description")
            self.assertEqual(records[0].sequence, "ATGCNN")

    def test_duplicate_contig_ids_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fna"
            path.write_text(">ctg first\nATG\n>ctg second\nATG\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                read_nucleotide_fasta(path)

    def test_auto_mode_uses_documented_ideal_training_boundary(self) -> None:
        small = [Contig("x", "x", "A" * 99_999, 0)]
        large = [Contig("x", "x", "A" * 100_000, 0)]
        self.assertEqual(choose_gene_mode(small, "auto"), "meta")
        self.assertEqual(choose_gene_mode(large, "auto"), "single")


if __name__ == "__main__":
    unittest.main()
