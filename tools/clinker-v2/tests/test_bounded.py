import json
import pytest
from Bio.Seq import Seq
from clinker.bounded import compare, SCHEMA, main


def inputs():
    dna = "ATG" + "GCT" * 25 + "TAA"
    gene = {"id": "cds1", "start": 1, "end": len(dna), "strand": "+", "phase": 0,
            "translation_table": 11, "partial_5prime": False, "partial_3prime": False,
            "translation": "M" + "A" * 25}
    document = {"schema": SCHEMA, "coordinate_system": "1-based-inclusive-source-forward", "loci": [
        {"id": "one", "source_id": "one", "start": 1, "end": len(dna), "genes": [gene]},
        {"id": "two", "source_id": "two", "start": 1, "end": len(dna), "genes": [{**gene, "strand": "-"}]}]}
    fasta = f">one\n{dna}\n>two\n{Seq(dna).reverse_complement()}\n".encode()
    return fasta, document


def test_true_alignment_and_reverse_source_binding(tmp_path):
    fasta, doc = inputs()
    report, _ = compare(fasta, json.dumps(doc).encode())
    assert report["evaluated_protein_pairs"] == 1
    assert report["links"][0]["identity"] == 1
    assert report["links"][0]["score"] == 105.0
    assert report["loci"][0]["genes"][0]["nucleotide_sha256"] == report["loci"][1]["genes"][0]["nucleotide_sha256"]
    (tmp_path / "input.fa").write_bytes(fasta)
    (tmp_path / "features.json").write_text(json.dumps(doc))
    args = ["compare", "--input", str(tmp_path / "input.fa"), "--features", str(tmp_path / "features.json"), "--outdir", str(tmp_path / "results")]
    main(args)
    assert sorted(p.name for p in (tmp_path / "results").iterdir()) == ["alignments.tsv", "comparison.json", "loci.gbk", "proteins.faa"]
    with pytest.raises(SystemExit):
        main(args)


@pytest.mark.parametrize("field,value", [("translation", "A"*26), ("translation_table", 99), ("phase", True), ("end", 100000), ("partial_5prime", None)])
def test_reject_forged_or_unbounded_features(field, value):
    fasta, doc = inputs()
    doc["loci"][0]["genes"][0][field] = value
    with pytest.raises(ValueError):
        compare(fasta, json.dumps(doc).encode())


def test_alternative_genetic_code_and_partial_phase():
    gene = {"id": "g", "start": 1, "end": 11, "strand": "+", "phase": 1,
            "translation_table": 4, "partial_5prime": True, "partial_3prime": True,
            "translation": "MWK"}
    doc = {"schema": SCHEMA, "coordinate_system": "1-based-inclusive-source-forward", "loci": [
        {"id": n, "source_id": n, "start": 1, "end": 11, "genes": [gene]} for n in ["a", "b"]]}
    report, _ = compare(b">a\nAATGTGAAAAC\n>b\nAATGTGAAAAC\n", json.dumps(doc).encode())
    assert report["links"][0]["identity"] == 1


def test_no_genes_is_valid_negative_result():
    fasta, doc = inputs()
    for locus in doc["loci"]:
        locus["genes"] = []
    report, _ = compare(fasta, json.dumps(doc).encode())
    assert report["links"] == []
    assert report["evaluated_protein_pairs"] == 0


def test_work_limit_checked_before_alignment():
    fasta, doc = inputs()
    for locus in doc["loci"]:
        locus["genes"] = [{**locus["genes"][0], "id": f"g{i:03}"} for i in range(70)]
    with pytest.raises(ValueError, match="Total CDS"):
        compare(fasta, json.dumps(doc).encode())


def test_duplicate_json_rejected():
    with pytest.raises(ValueError, match="Duplicate JSON"):
        compare(b">a\nAAA\n", b'{"schema":1,"schema":2}')


def test_insertion_produces_real_gap_evidence():
    _,doc=inputs()
    a="ATG"+"GCT"*25+"TAA"
    b="ATG"+"GCT"*12+"GGT"+"GCT"*13+"TAA"
    doc["loci"][1]["end"]=len(b)
    doc["loci"][1]["genes"][0].update(end=len(b),strand="+",translation="M"+"A"*12+"G"+"A"*13)
    report,_=compare(f">one\n{a}\n>two\n{b}\n".encode(),json.dumps(doc).encode())
    link=report["links"][0]
    assert "-" in link["aligned_query"]
    assert link["identity"]==26/27
    assert link["score"]==95


def test_export_preserves_leading_gap_residues_and_shields_identifiers(tmp_path):
    import csv
    from clinker.bounded import write_exports
    fasta, features = inputs()
    report, records = compare(fasta, json.dumps(features).encode())
    link = report["links"][0]
    link.update(query_locus="=1+1", aligned_query="--M" + "A" * 25, aligned_target="M" + "A" * 25 + "--")
    write_exports(tmp_path, report, records)
    with (tmp_path / "alignments.tsv").open() as handle:
        exported = next(csv.DictReader(handle, delimiter="\t"))
    assert exported["aligned_query"] == link["aligned_query"]
    assert exported["aligned_target"] == link["aligned_target"]
    assert exported["query_locus"] == "'=1+1"


def test_export_preserves_negative_scientific_score_and_shields_numeric_identifiers(tmp_path):
    import csv
    from clinker.bounded import write_exports
    gene = {"start": 1, "end": 12, "strand": "+", "phase": 0,
            "translation_table": 11, "partial_5prime": False, "partial_3prime": False}
    features = {"schema": SCHEMA, "coordinate_system": "1-based-inclusive-source-forward", "loci": [
        {"id": "-1.0", "source_id": "a", "start": 1, "end": 12,
         "genes": [{**gene, "id": "=1+1", "translation": "MCC"}]},
        {"id": "+1", "source_id": "b", "start": 1, "end": 12,
         "genes": [{**gene, "id": "-1.0", "translation": "MDD"}]}]}
    report, records = compare(b">a\nATGTGTTGTTAA\n>b\nATGGATGATTAA\n", json.dumps(features).encode())
    link = report["links"][0]
    assert link["identity"] == 1 / 3
    assert link["score"] == -1.0
    write_exports(tmp_path, report, records)
    with (tmp_path / "alignments.tsv").open() as handle:
        exported = next(csv.DictReader(handle, delimiter="\t"))
    exported_json = json.loads((tmp_path / "comparison.json").read_text())["links"][0]
    for key in ("identity", "similarity", "score"):
        assert float(exported[key]) == exported_json[key] == link[key]
    for key in ("query_locus", "query_gene", "target_locus", "target_gene"):
        assert exported[key] == "'" + link[key]


def test_report_metadata_mutation_cannot_change_later_comparisons(monkeypatch):
    from clinker import bounded
    # Keep this regression isolated even against the old aliasing behavior.
    monkeypatch.setattr(bounded, "ALIGNER", bounded.ALIGNER.copy())
    monkeypatch.setattr(bounded, "LIMITS", bounded.LIMITS.copy())
    settings, limits = bounded.ALIGNER.copy(), bounded.LIMITS.copy()
    gene = {"id": "g", "start": 1, "strand": "+", "phase": 0,
            "translation_table": 11, "partial_5prime": False, "partial_3prime": False}
    features = {"schema": SCHEMA, "coordinate_system": "1-based-inclusive-source-forward", "loci": [
        {"id": "a", "source_id": "a", "start": 1, "end": 9,
         "genes": [{**gene, "end": 9, "translation": "MA"}]},
        {"id": "b", "source_id": "b", "start": 1, "end": 12,
         "genes": [{**gene, "end": 12, "translation": "MAA"}]}]}
    fasta, raw = b">a\nATGGCTTAA\n>b\nATGGCTGCTTAA\n", json.dumps(features).encode()
    report, _ = compare(fasta, raw)
    original_link = report["links"][0].copy()
    assert "-" in original_link["aligned_query"]
    assert original_link["score"] == -1.0
    report["options"]["aligner"]["open_gap_score"] = 100
    report["limits"]["loci"] = 1
    assert bounded.ALIGNER == settings
    assert bounded.LIMITS == limits
    subsequent, _ = compare(fasta, raw)
    assert subsequent["options"]["aligner"] == settings
    assert subsequent["limits"] == limits
    assert subsequent["links"][0] == original_link
