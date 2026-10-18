import json
from casandra.orf_fast import extract_genome_proteins_fast


def test_fast_extraction_hashes_once_and_writes_manifest(tmp_path):
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\n" + "ATG" * 200 + "TAA\n")
    out = tmp_path / "out"
    extract_genome_proteins_fast(fasta, out, mode="meta")
    audit = json.loads((out / "audit.json").read_text())
    assert audit["schema_version"] == 3
    assert "computed once" in audit["performance_fix"]
    assert (out / "manifest.json").exists()
