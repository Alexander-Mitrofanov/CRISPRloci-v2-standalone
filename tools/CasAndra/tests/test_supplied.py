from pathlib import Path

from casandra.supplied import extract_gff3_proteins


def test_gff3_adapter_translates_and_preserves_coordinates(tmp_path: Path):
    fasta = tmp_path / "genome.fna"
    gff = tmp_path / "genes.gff3"
    fasta.write_text(">ctg description\nATGAAATAA\n", encoding="utf-8")
    gff.write_text("##gff-version 3\nctg\ttest\tCDS\t1\t9\t.\t+\t0\tID=p1\n", encoding="utf-8")
    output = tmp_path / "out"
    extract_gff3_proteins(fasta, gff, output)
    row = __import__("json").loads((output / "genes.jsonl").read_text())
    assert row["protein_id"] == "p1"
    assert row["protein_sequence"] == "MK"
    assert (row["start_1based"], row["end_1based_inclusive"]) == (1, 9)
    manifest = __import__("json").loads((output / "manifest.json").read_text())
    assert set(manifest["inputs"]) == {"annotation", "genome"}


def test_gff3_adapter_rejects_duplicate_protein_ids(tmp_path: Path):
    fasta = tmp_path / "genome.fna"
    gff = tmp_path / "genes.gff3"
    fasta.write_text(">ctg\nATGAAATAAATGAAATAA\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\n"
        "ctg\ttest\tCDS\t1\t9\t.\t+\t0\tID=p1\n"
        "ctg\ttest\tCDS\t10\t18\t.\t+\t0\tID=p1\n",
        encoding="utf-8",
    )
    with __import__("pytest").raises(ValueError, match="Duplicate"):
        extract_gff3_proteins(fasta, gff, tmp_path / "out")
