"""Structured JSONL, TSV, GFF3, and HTML report writers."""

from __future__ import annotations

import csv
import html
import json
from pathlib import Path
from urllib.parse import quote


def read_jsonl(path: Path):
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def write_jsonl(path: Path, rows) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        for row in rows:
            json.dump(row, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")


def _gff_value(value: object, *, sequence_id: bool = False) -> str:
    safe = "._:-" if sequence_id else "._-|="
    return quote(str(value), safe=safe)


def write_outputs(staging: Path, genes, predictions, cassettes) -> None:
    prediction_by_id = {row["sequence_id"]: row for row in predictions}
    protein_rows = [
        {**gene, "prediction": prediction_by_id[gene["protein_id"]]} for gene in genes
    ]
    write_jsonl(staging / "proteins.jsonl", protein_rows)
    write_jsonl(staging / "cassettes.jsonl", cassettes)

    with (staging / "cas_proteins.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["protein_id", "contig", "start", "end", "strand", "type", "subtype", "margin", "profile"]
        )
        for gene in genes:
            prediction = prediction_by_id[gene["protein_id"]]
            if not prediction["is_cas"]:
                continue
            writer.writerow(
                [
                    gene["protein_id"],
                    gene["contig_id"],
                    gene["start_1based"],
                    gene["end_1based_inclusive"],
                    gene["strand"],
                    prediction["classification"]["type"],
                    prediction["classification"]["subtype"],
                    f"{prediction['score_margin']:.6g}",
                    prediction["best_positive_profile"],
                ]
            )

    with (staging / "cassettes.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["cassette_id", "contig", "start", "end", "cas_genes", "class", "type", "subtype", "method", "confidence"]
        )
        for cassette in cassettes:
            final = cassette["final_classification"]
            writer.writerow(
                [
                    cassette["cassette_id"],
                    cassette["contig_id"],
                    cassette["start_1based"],
                    cassette["end_1based_inclusive"],
                    cassette["cas_gene_count"],
                    final["class"],
                    final["type"],
                    final["subtype"],
                    final["method"],
                    f"{final['confidence']:.6g}",
                ]
            )

    with (staging / "casandra.gff3").open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for cassette in cassettes:
            final = cassette["final_classification"]
            attributes = (
                f"ID={_gff_value(cassette['cassette_id'])};"
                f"cas_type={_gff_value(final['type'] or '')};"
                f"cas_subtype={_gff_value(final['subtype'] or '')};"
                f"method={_gff_value(final['method'])};confidence={final['confidence']:.6g}"
            )
            handle.write(
                f"{_gff_value(cassette['contig_id'], sequence_id=True)}\tCasAndra\tCas_cassette\t"
                f"{cassette['start_1based']}\t{cassette['end_1based_inclusive']}\t.\t.\t.\t"
                f"{attributes}\n"
            )
        for gene in genes:
            prediction = prediction_by_id[gene["protein_id"]]
            if not prediction["is_cas"]:
                continue
            attributes = (
                f"ID={_gff_value(gene['protein_id'])};"
                f"cas_type={_gff_value(prediction['classification']['type'] or '')};"
                f"cas_subtype={_gff_value(prediction['classification']['subtype'] or '')};"
                f"profile={_gff_value(prediction['best_positive_profile'] or '')};"
                f"margin={prediction['score_margin']:.6g}"
            )
            handle.write(
                f"{_gff_value(gene['contig_id'], sequence_id=True)}\tCasAndra\tCas_protein\t"
                f"{gene['start_1based']}\t{gene['end_1based_inclusive']}\t"
                f"{prediction['positive_profile_score']:.6g}\t{gene['strand']}\t0\t{attributes}\n"
            )

    rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(value if value is not None else ''))}</td>"
            for value in (
                cassette["contig_id"],
                cassette["start_1based"],
                cassette["end_1based_inclusive"],
                cassette["cas_gene_count"],
                cassette["final_classification"]["type"],
                cassette["final_classification"]["subtype"],
                cassette["final_classification"]["method"],
            )
        )
        + "</tr>"
        for cassette in cassettes
    )
    (staging / "report.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>CasAndra report</title>"
        "<style>body{font-family:system-ui;margin:2rem}table{border-collapse:collapse}"
        "th,td{border:1px solid #bbb;padding:.35rem}</style><h1>CasAndra Cas-only report</h1>"
        "<p>CRISPR arrays are outside this tool's scope.</p><table><thead><tr>"
        "<th>Contig</th><th>Start</th><th>End</th><th>Cas genes</th><th>Type</th>"
        "<th>Subtype</th><th>Method</th></tr></thead><tbody>"
        + rows
        + "</tbody></table>",
        encoding="utf-8",
    )

