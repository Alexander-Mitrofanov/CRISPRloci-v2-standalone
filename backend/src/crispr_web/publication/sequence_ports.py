"""Complete sequence handoffs, bounded independently of browser previews."""

from __future__ import annotations

import hashlib
import json
import os

from ..exchange import INPUT_KINDS, MAX_EXCHANGE_BYTES, sequence_set, strict_json
from ..fasta import normalize_fasta
from ..summaries.bounded_io import _read_bounded
from .array_data import array_documents


def write_document(root, name, document, *, max_bytes, checkpoint=None):
    """Bound serialized bytes while writing; never publish partial JSON."""
    temporary = root / f".{name}.tmp"
    size = 0
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            for part in json.JSONEncoder(indent=2, allow_nan=False).iterencode(
                document
            ):
                if checkpoint:
                    checkpoint()
                size += len(part.encode("utf-8"))
                if size > min(max_bytes, MAX_EXCHANGE_BYTES):
                    raise ValueError("sequence export exceeds byte limit")
                stream.write(part)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary, 0o600)
        os.replace(temporary, root / name)
    finally:
        temporary.unlink(missing_ok=True)


def write_sequence_exports(root, job_root, claimed, settings, checkpoint=None):
    warnings = []
    limit = min(settings.max_artifact_file_bytes, MAX_EXCHANGE_BYTES)
    try:
        payload = _read_bounded(
            job_root,
            job_root / "input/normalized.fasta",
            max_bytes=limit,
            checkpoint=checkpoint,
        )
        normalized = normalize_fasta(
            payload.decode("ascii"),
            max_total_bases=settings.max_total_bases,
            max_record_bases=settings.max_record_bases,
            max_records=settings.max_records,
            max_header_characters=settings.max_header_characters,
            molecule=claimed.molecule,
        )
        origins = {}
        metadata = job_root / "handoff.json"
        if metadata.exists():
            origins = strict_json(
                _read_bounded(
                    job_root, metadata, max_bytes=limit, checkpoint=checkpoint
                ).decode("utf-8")
            )["origins"]
        if claimed.mode == "array_compare":
            from ..array_exchange import native_input

            mapping = strict_json(_read_bounded(
                job_root, job_root / "array-map.json", max_bytes=limit,
                checkpoint=checkpoint,
            ).decode("utf-8"))
            document = mapping["document"]
            _, expected = native_input(document)
            if expected != {row.identifier: row.sequence for row in normalized.records}:
                raise ValueError("comparison sequence map differs from input")
            sources_by_id = {row["id"]: row for row in document["sources"]}
            for index, array in enumerate(document["arrays"]):
                for unit in array["spacers"]:
                    if checkpoint:
                        checkpoint()
                    if unit.get("status") == "deletion":
                        continue
                    origins[f"occ{index:06d}_{unit['ordinal']:06d}"] = {
                        "source_id": array["source_id"],
                        "source_sha256": sources_by_id[array["source_id"]]["sha256"],
                        "array_id": array["id"],
                        "ordinal": unit["ordinal"],
                        "start": unit["start"],
                        "end": unit["end"],
                    }
        records = [
            {
                "id": row.identifier,
                "sequence": row.sequence,
                **(
                    {"origin": origins[row.identifier]}
                    if row.identifier in origins
                    else {}
                ),
            }
            for row in normalized.records
        ]
        packet = sequence_set("spacer" if claimed.mode == "array_compare" else INPUT_KINDS[claimed.mode], claimed.molecule, records)
        write_document(
            root, "input-sequences.json", packet, max_bytes=limit, checkpoint=checkpoint
        )
        if claimed.mode == "array_compare":
            write_document(root, "spacers.json", packet, max_bytes=limit, checkpoint=checkpoint)
        sources = {
            row.identifier: (
                row.sequence,
                hashlib.sha256(row.sequence.encode()).hexdigest(),
            )
            for row in normalized.records
        }
        if claimed.mode in {
            "detection",
            "orientation",
            "reconstruction",
            "leader",
            "repeat_context",
            "prophage_context",
            "loci",
        }:
            from ..service import CATEGORY_POLICY_CATEGORIES

            documents = array_documents(
                root,
                sources,
                categories=CATEGORY_POLICY_CATEGORIES[claimed.category_policy],
                mode=claimed.mode,
                max_files=settings.max_bundle_files,
                max_bytes=limit,
                checkpoint=checkpoint,
            )
            for name, document in documents.items():
                write_document(
                    root, name, document, max_bytes=limit, checkpoint=checkpoint
                )
            if any(
                not documents[name]["complete"]
                for name in ("repeats.json", "spacers.json")
            ):
                warnings.append(
                    "Some array sequences lack verified source intervals; incomplete handoffs cannot be submitted."
                )
    except (ValueError, OSError, UnicodeError, RecursionError):
        warnings.append(
            "Reusable sequence exports could not be completed within the validated data and size limits."
        )
    return warnings
