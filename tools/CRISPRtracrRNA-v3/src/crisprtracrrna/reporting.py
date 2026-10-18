"""Publication of retained JSON evidence and deterministic prediction exports."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any
from urllib.parse import quote
from uuid import uuid4

LOG = logging.getLogger(__name__)


def atomic_json(path: Path, value: Any) -> None:
    """Replace a JSON file atomically, cleaning up an abandoned publication attempt."""
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    created = False
    try:
        # Exclusive creation gives each attempt its own sibling and respects umask.
        with temporary.open("x", encoding="utf-8") as handle:
            created = True
            handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
        temporary.replace(path)
    finally:
        if created:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                # Cleanup is diagnostic; it must not hide the publication failure.
                LOG.warning("Could not remove temporary JSON file %s", temporary, exc_info=True)


def _assign_export_ids(record: dict[str, Any]) -> None:
    """Give predictions the same identities in retained records and final exports."""
    for ordinal, prediction in enumerate(record["predictions"], 1):
        prediction["export_id"] = f"{record['key']}.p{ordinal:06d}"


def write_exports(output: Path, records: list[dict[str, Any]]) -> None:
    fields = [
        "source_file",
        "sequence_id",
        "prediction_id",
        "candidate_origin",
        "admission_branch",
        "evidence_class",
        "boundary_kind",
        "start",
        "end",
        "strand",
        "sequence",
        "score",
        "evalue",
        "array_id",
        "array_category",
        "anti_repeat_identity",
        "anti_repeat_coverage",
        "interaction_energy_kcal_mol",
        "model_score",
        "terminator_evalue",
        "cas_cassette_id",
        "cas_distance_bases",
    ]
    with (
        (output / "summary.csv").open("w", newline="", encoding="utf-8") as handle,
        (output / "predictions.fasta").open("w", encoding="utf-8") as fasta,
        (output / "predictions.gff3").open("w", encoding="utf-8") as gff,
    ):
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        gff.write("##gff-version 3\n")
        for record in records:
            _assign_export_ids(record)
            source = record["source"]
            # Prefix with record key to avoid collisions across files with reused contig IDs.
            seqid = quote(record["key"] + ":" + source["id"], safe="._:-")
            gff.write(f"##sequence-region {seqid} 1 {source['length']}\n")
            for prediction in record["predictions"]:
                identity = prediction["export_id"]
                interval = prediction["interval"]
                ranking = prediction.get("ranking", {})
                interaction_row = prediction.get("interaction") or {}
                model = prediction.get("model_hit") or {}
                term = prediction.get("terminator") or {}
                cas = prediction.get("cas_context") or {}
                anti = prediction.get("anti_repeat") or {}
                array = prediction.get("array") or {}
                writer.writerow(
                    {
                        "source_file": source["file"],
                        "sequence_id": source["id"],
                        "prediction_id": identity,
                        "candidate_origin": prediction.get("candidate_origin", "covariance_model"),
                        "admission_branch": prediction.get("admission_branch", "model_mode"),
                        "evidence_class": prediction.get("evidence_class", "cm_only"),
                        "boundary_kind": prediction.get("boundary_kind", "merged_model_region"),
                        "start": interval["start"],
                        "end": interval["end"],
                        "strand": prediction["strand"],
                        "sequence": prediction["sequence"],
                        "score": ranking.get("score", prediction.get("best_score")),
                        "evalue": prediction.get("best_evalue"),
                        "array_id": array.get("id"),
                        "array_category": array.get("category"),
                        "anti_repeat_identity": anti.get("identity"),
                        "anti_repeat_coverage": anti.get("coverage"),
                        "interaction_energy_kcal_mol": interaction_row.get("energy_kcal_mol"),
                        "model_score": model.get("score"),
                        "terminator_evalue": term.get("evalue"),
                        "cas_cassette_id": cas.get("id"),
                        "cas_distance_bases": cas.get("distance_bases"),
                    }
                )
                fasta.write(
                    f">{identity} source={quote(source['id'], safe='._-')} "
                    f"start={interval['start']} end={interval['end']} "
                    f"strand={prediction['strand']}\n"
                    f"{prediction['sequence']}\n"
                )
                gff.write(
                    f"{seqid}\tCRISPRtracrRNA\tncRNA\t{interval['start']}\t{interval['end']}\t.\t"
                    f"{prediction['strand']}\t.\tID={identity};inference=ab%20initio%20prediction;"
                    "product=tracrRNA%20candidate\n"
                )
