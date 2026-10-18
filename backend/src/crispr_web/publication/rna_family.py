"""Allowlisted complete family evidence and portable alignment downloads."""

import csv

from ..summaries.common import _json
from ..summaries.rna_family import public_report
from .sequence_ports import write_document
from .context_text import atomic_text


def write_rna_family_exports(root, summary, checkpoint=None):
    if summary.get("rna_family", {}).get("available") is not True:
        return
    doc = public_report(
        _json(
            root,
            root / "rna-family/report.json",
            max_bytes=20_000_000,
            checkpoint=checkpoint,
        )
    )
    write_document(
        root, "rna-family.json", doc, max_bytes=20_000_000, checkpoint=checkpoint
    )
    with atomic_text(
        root, "rna-family-alignment.fasta", encoding="ascii", checkpoint=checkpoint
    ) as stream:
        for row in doc["alignment"]:
            stream.write(f">{row['id']}\n{row['sequence']}\n")
    with atomic_text(
        root, "rna-family-consensus.dbn", encoding="ascii", checkpoint=checkpoint
    ) as stream:
        stream.write(
            ">RNA_family_consensus\n"
            + doc["fold"]["consensus"]
            + "\n"
            + doc["fold"]["structure"]
            + "\n"
        )
    with atomic_text(root, "rna-family-pairs.tsv", checkpoint=checkpoint) as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        fields = (
            "left",
            "right",
            "probability",
            "canonical_sequences",
            "noncanonical_sequences",
            "gapped_sequences",
            "compensatory_sequence_pairs",
        )
        writer.writerow(fields)
        for row in doc["pairs"]:
            if checkpoint:
                checkpoint()
            writer.writerow([row[key] for key in fields])
