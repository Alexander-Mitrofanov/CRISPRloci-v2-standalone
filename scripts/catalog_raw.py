"""Normalize exceptional committed JSON outcomes without a predecessor runtime.

May be appended to catalog_snapshot.py for execution over SSH. Only the raw
documents selected by the verified import plan are read; no remote file is
created or changed.
"""

import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sqlite3

if "COLUMNS" not in globals():
    from catalog_snapshot import COLUMNS, create_tables, normalize_database


RAW_PROFILE = "5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540"
RAW_LIMIT = 256 * 1024**2


def raw_dict(value):
    return value if isinstance(value, dict) else {}


def raw_text(value):
    return str(value) if isinstance(value, (str, int, float)) and not isinstance(value, bool) else None


def raw_integer(value):
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def raw_number(value):
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) else None


def raw_boolean(value):
    return int(value) if isinstance(value, bool) else None


def raw_interval(value):
    interval = raw_dict(value)
    start, end, length = (raw_integer(interval.get(key)) for key in ("start", "end", "length"))
    if length is None and start is not None and end is not None and end >= start:
        length = end - start + 1
    return {"start": start, "end": end, "length": length}


def raw_classification(item):
    native = raw_dict(item.get("native"))
    classification = raw_dict(item.get("classification")) or raw_dict(native.get("classification"))
    return {column: raw_text(classification.get(key)) or raw_text(native.get(key)) or raw_text(native.get(key.capitalize()))
            for column, key in (("cas_class", "class"), ("cas_type", "type"), ("subtype", "subtype"))}


def raw_list(item, field):
    values = item.get(field, [])
    if values is None:
        return []
    if not isinstance(values, list) or any(not isinstance(value, dict) for value in values):
        raise ValueError(f"Malformed raw scientific list: {field}")
    return values


def raw_insert(connection, table, values):
    columns = COLUMNS[table][1:]
    cursor = connection.execute(
        f"INSERT INTO {table} ({','.join(columns)}) VALUES ({','.join('?' for _ in columns)})",
        [values.get(column) for column in columns],
    )
    return cursor.lastrowid


def raw_tool(connection, sequence_pk, tool):
    name, version = raw_text(tool.get("name")), raw_text(tool.get("version"))
    if not name:
        raise ValueError("Raw scientific tool name is absent")

    def common(item, index, prefix):
        return {"sequence_pk": sequence_pk, "tool": raw_text(item.get("tool")) or name,
                "tool_version": version, "call_id": raw_text(item.get("call_id")) or f"{prefix}-{index}",
                "reported_sequence_id": raw_text(item.get("sequence_id")),
                **raw_interval(item.get("interval"))}

    arrays, systems, genes = (raw_list(tool, field) for field in ("arrays", "cas_systems", "cas_genes"))
    for index, item in enumerate(arrays, 1):
        confidence = raw_dict(item.get("confidence"))
        array_pk = raw_insert(connection, "arrays", {
            **common(item, index, "array"), "strand": raw_text(item.get("strand")),
            "repeat_consensus": raw_text(item.get("repeat_consensus")),
            "repeat_count": raw_integer(item.get("repeat_count")),
            "spacer_count": raw_integer(item.get("spacer_count")),
            "evidence_level": raw_text(confidence.get("evidence_level")),
            "category": raw_text(confidence.get("category")), "certainty_score": raw_number(confidence.get("score")),
        })
        for field, table, kind in (("repeats", "repeats", None), ("spacers", "spacers", "observed"), ("spacer_deletions", "spacers", "deletion")):
            for ordinal, child in enumerate(raw_list(item, field), 1):
                values = {"array_pk": array_pk, "ordinal": raw_integer(child.get("index")) or ordinal,
                          **raw_interval(child.get("interval")), "sequence": raw_text(child.get("sequence")),
                          "reported_sequence_id": raw_text(child.get("sequence_id"))}
                if table == "spacers":
                    values.update(kind=kind, boundary_position=raw_integer(raw_dict(child.get("source_boundary")).get("position")))
                if kind == "deletion":
                    values["length"] = raw_integer(child.get("length"))
                raw_insert(connection, table, values)
    gene_to_system = {}
    for index, item in enumerate(systems, 1):
        native = raw_dict(item.get("native"))
        classification = raw_dict(item.get("classification")) or raw_dict(native.get("classification"))
        confidence = raw_number(classification.get("confidence"))
        if confidence is None:
            confidence = raw_number(native.get("confidence"))
        system_pk = raw_insert(connection, "cas_systems", {
            **common(item, index, "cas-system"), **raw_classification(item),
            "confidence": confidence,
            "confidence_kind": raw_text(classification.get("confidence_kind")) or raw_text(native.get("confidence_kind")),
            "experimental": raw_boolean(native.get("experimental")),
            "missing_gene_count": raw_integer(native.get("missing_gene_count")),
        })
        for gene in raw_list(native, "genes"):
            gene_id = raw_text(gene.get("protein_id")) or raw_text(gene.get("call_id"))
            if gene_id is not None:
                gene_to_system.setdefault(gene_id, system_pk)
    for index, item in enumerate(genes, 1):
        native = raw_dict(item.get("native"))
        values = common(item, index, "cas-gene")
        profile = native.get("profile")
        profile_text = raw_text(profile)
        if profile is not None and profile_text is None:
            profile_text = json.dumps(profile, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        raw_insert(connection, "cas_genes", {
            **values, **raw_classification(item), "cas_system_pk": gene_to_system.get(values["call_id"]),
            "strand": raw_text(item.get("strand")) or raw_text(native.get("strand")),
            "annotation": raw_text(native.get("annotation")), "profile": profile_text,
            "score": raw_number(native.get("score")), "score_kind": raw_text(native.get("score_kind")),
            "protein_length": raw_integer(native.get("protein_length")),
            "translation_table": raw_integer(native.get("translation_table")), "is_cas": raw_boolean(native.get("is_cas")),
        })
    return len(arrays), len(systems), len(genes)


def raw_assembly(connection, payload, entry):
    accession = payload.get("assembly", {}).get("accession")
    status = payload.get("status")
    if accession != entry["accession"] or status != entry["status"] or status not in {"complete", "failed"}:
        raise ValueError("Raw outcome disagrees with the selected source inventory")
    if raw_dict(payload.get("provenance")).get("analysis_profile_sha256") != RAW_PROFILE:
        raise ValueError("Unexpected raw scientific profile")
    validation_status = raw_text(raw_dict(payload.get("validation")).get("status"))
    if status == "complete" and validation_status != "valid":
        raise ValueError("Raw successful result has not passed source validation")
    assembly_pk = raw_insert(connection, "assemblies", {
        "accession": accession, "status": status, "validation_status": validation_status,
        "sequence_count": 0, "array_count": 0, "cas_system_count": 0, "cas_gene_count": 0,
    })
    if status == "failed":
        return
    totals = [0, 0, 0, 0]
    for ordinal, sequence in enumerate(raw_list(payload, "sequences"), 1):
        sequence_id = raw_text(sequence.get("id"))
        if not sequence_id:
            raise ValueError("Raw sequence identifier is absent")
        sequence_pk = raw_insert(connection, "sequences", {
            "assembly_pk": assembly_pk, "sequence_id": sequence_id, "sequence_ordinal": ordinal,
            "header": raw_text(sequence.get("header")), "length": raw_integer(sequence.get("length")),
            "sha256": raw_text(sequence.get("sha256")),
        })
        totals[0] += 1
        for tool in raw_list(sequence, "tools"):
            counts = raw_tool(connection, sequence_pk, tool)
            for index, count in enumerate(counts, 1):
                totals[index] += count
    connection.execute("UPDATE assemblies SET sequence_count=?,array_count=?,cas_system_count=?,cas_gene_count=? WHERE assembly_pk=?", [*totals, assembly_pk])


def transform_raw(items, shard_id):
    """Return ``(receipt, gzip_sqlite)`` for selected raw and unavailable records."""
    source = sqlite3.connect(":memory:")
    documents = []
    selected = set()
    try:
        create_tables(source)
        for item in items:
            path = Path(item["path"])
            if path.is_symlink() or not path.is_file() or path.stat().st_size > RAW_LIMIT:
                raise ValueError("Invalid or oversized raw source")
            with path.open("rb") as stream:
                data = stream.read(RAW_LIMIT + 1)
            if len(data) > RAW_LIMIT or hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise ValueError("Raw source checksum mismatch")
            decoded = data
            if path.suffix == ".gz":
                with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
                    decoded = stream.read(RAW_LIMIT + 1)
            if len(decoded) > RAW_LIMIT:
                raise ValueError("Expanded raw source exceeds bound")
            payload = json.loads(decoded)
            documents.append({"sha256": item["sha256"], "kind": item["kind"],
                              "accessions": [entry["accession"] for entry in item["entries"]]})
            for entry in item["entries"]:
                if entry["accession"] in selected or entry["sha256"] != item["sha256"]:
                    raise ValueError("Duplicated or inconsistent selected raw outcome")
                selected.add(entry["accession"])
                if item["kind"] == "unavailable":
                    if entry["status"] != "unavailable" or entry["accession"] not in payload.get("unavailable", {}):
                        raise ValueError("Unavailable outcome is not present in its receipt")
                    raw_insert(source, "assemblies", {
                        "accession": entry["accession"], "status": "unavailable", "validation_status": None,
                        "sequence_count": 0, "array_count": 0, "cas_system_count": 0, "cas_gene_count": 0,
                    })
                elif item["kind"] == "raw":
                    raw_assembly(source, payload, entry)
                else:
                    raise ValueError("Unexpected raw source kind")
        source.commit()
        compressed, registry, counts, size = normalize_database(source.serialize(), shard_id)
        source_hash = hashlib.sha256(json.dumps(documents, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        receipt = {
            "id": shard_id, "sha256": hashlib.sha256(compressed).hexdigest(), "bytes": len(compressed),
            "uncompressed_bytes": size, "counts": counts, "assemblies": registry,
            "source_sha256": source_hash, "source_documents": documents,
            "analysis_profile_sha256": RAW_PROFILE,
        }
        return receipt, compressed
    finally:
        source.close()
