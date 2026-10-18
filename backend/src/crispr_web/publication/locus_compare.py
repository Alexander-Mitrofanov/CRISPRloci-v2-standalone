"""Protected sanitized comparison evidence; upstream HTML stays private."""

import csv
from ..summaries.locus_compare import public_report
from .sequence_ports import write_document
from .context_text import atomic_text
from .tabular import scientific_cell


def write_locus_compare_exports(root, summary, checkpoint=None):
    if summary.get("locus_compare", {}).get("available") is not True:
        return
    report = public_report(root, checkpoint)
    write_document(
        root,
        "locus-comparison.json",
        report,
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    with atomic_text(root, "locus-alignments.tsv", checkpoint=checkpoint) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        keys = (
            "query_locus",
            "query_gene",
            "target_locus",
            "target_gene",
            "identity",
            "similarity",
            "score",
            "aligned_query",
            "aligned_target",
        )
        writer.writerow(keys)
        for link in report["links"]:
            if checkpoint:
                checkpoint()
            row = []
            for key in keys:
                value = link[key]
                value = scientific_cell(value, protein_alignment=key in {"aligned_query", "aligned_target"})
                row.append(value)
            writer.writerow(row)
    # Recreate FASTA from checked JSON rather than copying unchecked native exports.
    with atomic_text(
        root, "locus-proteins.faa", encoding="ascii", checkpoint=checkpoint
    ) as handle:
        for i, locus in enumerate(report["loci"], 1):
            for j, gene in enumerate(locus["genes"], 1):
                if checkpoint:
                    checkpoint()
                handle.write(
                    f">locus{i:03d}_gene{j:03d} locus={locus['id']} gene={gene['id']} source={locus['source_id']}:{gene['start']}-{gene['end']} strand={gene['strand']}\n"
                )
                for pos in range(0, len(gene["translation"]), 80):
                    handle.write(gene["translation"][pos : pos + 80] + "\n")
