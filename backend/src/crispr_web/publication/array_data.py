"""Narrow array data and source-verified repeat/spacer ports.

No summary limits, native feature dumps, duplicated architecture, paths or logs.
Coordinates are 1-based inclusive in the source-forward frame, irrespective of
the predicted biological orientation. Deletion slots are not DNA sequences.
"""

from __future__ import annotations

import hashlib
import json
import re

from ..exchange import sequence_set
from ..summaries.bounded_io import _json, _safe_file
from ..summaries.detection import _finite_score, _first, _report_arrays


def _interval(value):
    if isinstance(value, dict):
        start, end = value.get("start"), value.get("end")
        if type(start) is int and type(end) is int and 1 <= start <= end:
            return {"start": start, "end": end}
    return None


def _units(raw, *, kind, source_id, source, digest, array_id, port, complete):
    units = []
    if not isinstance(raw, list):
        complete[kind] = False
        return None
    for index, item in enumerate(raw, 1):
        item = {"sequence": item} if isinstance(item, str) else item
        if not isinstance(item, dict):
            complete[kind] = False
            units.append({"ordinal": index, "sequence": None, "status": "not_reported"})
            continue
        interval = _interval(item.get("source_interval"))
        sequence = item.get("sequence")
        deletion = item.get("kind") == "deletion"
        if deletion and (sequence not in (None, "") or interval is not None):
            raise ValueError("deletion slot contains a nonempty sequence or interval")
        verified = (
            not deletion
            and isinstance(sequence, str)
            and bool(sequence)
            and interval
            and sequence == source[interval["start"] - 1 : interval["end"]]
        )
        row = {"ordinal": index, "sequence": sequence if verified else None}
        if interval:
            row.update(interval)
        if deletion:
            row["status"] = "deletion"
            position = item.get("boundary_position")
            if type(position) is int and 0 <= position <= len(source):
                row["boundary_position"] = position
        elif not verified:
            row["status"] = "not_source_verified"
            complete[kind] = False
        else:
            identity = json.dumps(
                [source_id, digest, array_id, kind, index], separators=(",", ":")
            )
            row["id"] = f"{kind}_" + hashlib.sha256(identity.encode()).hexdigest()[:24]
            port.append(
                {
                    "id": row["id"],
                    "sequence": sequence,
                    "origin": {
                        "source_id": source_id,
                        "source_sha256": digest,
                        "array_id": array_id,
                        "ordinal": index,
                        **interval,
                    },
                }
            )
        units.append(row)
    return units


def array_documents(
    root, sources, *, categories, mode, max_files, max_bytes, checkpoint=None
):
    search_root = root / ("tracrrna/work" if mode == "loci" else "crispridentify")
    arrays, source_rows, seen = [], [], set()
    ports = {"repeat": [], "spacer": []}
    complete = {"repeat": True, "spacer": True}
    budget = max_bytes
    source_aliases = {}
    if mode == "loci":
        report = root / "tracrrna/report.json"
        doc = _json(root, report, max_bytes=budget, checkpoint=checkpoint)
        budget -= report.stat().st_size
        records = doc.get("records")
        if doc.get("status") != "complete" or not isinstance(records, list):
            raise ValueError("incomplete locus source mapping")
        mapped = set()
        for record in records:
            if checkpoint:
                checkpoint()
            if not isinstance(record, dict):
                raise ValueError("invalid locus source mapping")
            key, binding = record.get("key"), record.get("source")
            if not isinstance(binding, dict):
                raise ValueError("invalid locus source mapping")
            source_id = binding.get("id")
            if (
                not isinstance(key, str)
                or not re.fullmatch(r"f[0-9]{6}\.s[0-9]{6}", key)
                or key in source_aliases
                or not isinstance(source_id, str)
                or source_id not in sources
                or source_id in mapped
                or binding.get("length") != len(sources[source_id][0])
                or binding.get("sha256") != sources[source_id][1]
            ):
                raise ValueError("locus source mapping differs from input")
            source_aliases[key] = source_id
            mapped.add(source_id)
        if mapped != set(sources):
            raise ValueError("locus source mapping does not cover input")
    scanned = 0
    for path in search_root.rglob("*"):
        scanned += 1
        if checkpoint:
            checkpoint()
        if scanned > max_files:
            raise ValueError("array discovery exceeds file limit")
        if path.name != "report.json" or (
            mode == "loci" and "identify" not in path.relative_to(search_root).parts
        ):
            continue
        if _safe_file(root, path) is None:
            raise ValueError("unsafe array report")
        doc = _json(root, path, max_bytes=budget, checkpoint=checkpoint)
        budget -= path.stat().st_size
        schema, binding = doc.get("schema", {}), doc.get("source", {})
        source_id = binding.get("id")
        if mode == "loci":
            # The complete workflow renames each retained input to its record
            # key. Bind that key and directory before restoring the public ID.
            key = path.relative_to(search_root).parts[0]
            if source_id != key or key not in source_aliases:
                raise ValueError("array source differs from locus record")
            source_id = source_aliases[key]
        if (
            schema.get("name") != "CRISPRidentify-report"
            or source_id not in sources
            or source_id in seen
        ):
            raise ValueError("invalid array source")
        source, digest = sources[source_id]
        if (
            binding.get("sha256") != digest
            or binding.get("length") != len(source)
            or binding.get("sha256_scope") != "uppercase_sequence"
        ):
            raise ValueError("array source differs from input")
        seen.add(source_id)
        source_rows.append({"id": source_id, "length": len(source), "sha256": digest})
        array_ids = set()
        for raw in _report_arrays(doc):
            if checkpoint:
                checkpoint()
            category = _first(raw, ("category", "classification", "category_name"))
            if category not in categories:
                continue
            array_id = _first(raw, ("id", "array_id", "name"))
            if not isinstance(array_id, str) or not array_id or array_id in array_ids:
                raise ValueError("array identifier missing or duplicated")
            array_ids.add(array_id)
            orientation = _first(raw, ("orientation", "direction", "strand"))
            if isinstance(orientation, dict):
                orientation = _first(
                    orientation, ("strand", "label", "direction", "status")
                )
            # Keep unknown orientation explicit; never default it to forward.
            strand = (
                {"+": "+", "-": "-", "Forward": "+", "Reverse": "-"}.get(orientation)
                if isinstance(orientation, str)
                else None
            )
            row = {
                "id": array_id,
                "source_id": source_id,
                "category": category,
                "model_score": _finite_score(
                    _first(raw, ("model_score", "certainty_score", "score"))
                ),
                "strand": strand,
                "interval": _interval(raw.get("source_interval")),
            }
            boundary = raw.get("boundary_assessment", {})
            row["boundary"] = {
                "exact": boundary.get("reported_geometry", {}).get(
                    "exact_boundary_claim"
                ),
                **{
                    name: value.get("status")
                    for name, value in boundary.get("dimensions", {}).items()
                    if name in {"repeat_phase", "occurrence_path"}
                    and isinstance(value, dict)
                },
            }
            for kind, aliases in (
                ("repeat", ("repeats", "repeat_sequences")),
                ("spacer", ("spacers", "spacer_sequences")),
            ):
                row[kind + "s"] = _units(
                    _first(raw, aliases),
                    kind=kind,
                    source_id=source_id,
                    source=source,
                    digest=digest,
                    array_id=array_id,
                    port=ports[kind],
                    complete=complete,
                )
            arrays.append(row)
    if seen != set(sources):
        raise ValueError("array reports do not cover input sources")
    return {
        "arrays.json": {
            "schema": {"name": "crisprloci.arrays", "version": "1.0.0"},
            "coordinate_system": "1-based-inclusive-source-forward",
            "score_semantics": "uncalibrated-model-score",
            "categories": list(categories),
            "sources": sorted(source_rows, key=lambda row: row["id"]),
            "arrays": sorted(arrays, key=lambda row: (row["source_id"], row["id"])),
        },
        **{
            kind + "s.json": sequence_set(
                kind,
                "DNA",
                sorted(rows, key=lambda row: row["id"]),
                complete=complete[kind],
            )
            for kind, rows in ports.items()
        },
    }
