"""Source-bound query inventories and a strict CRISPRidentify 1.2.0 adapter."""

from __future__ import annotations

from pathlib import Path

from .io import DNA, MAX_BASES, MAX_RECORDS, InputError, LimitError, decode_json, digest, identifier, integer, read_bytes

SCHEMA = {"name": "crisprspacer.inventory", "version": "1.0.0"}
CATEGORIES = {"Bona-fide", "Possible", "Alternative", "Possible Discarded", "Low score"}
RECORD_KEYS = {"id", "sequence", "sha256", "length", "assembly_id", "topology"}
QUERY_KEYS = {"id", "sequence", "sha256", "kind", "source_id", "array_id", "ordinal", "start", "end",
              "boundary", "category", "source_orientation"}
ARRAY_KEYS = {"id", "native_id", "source_id", "start", "end", "category", "orientation"}


def validate_records(records: object) -> dict[str, dict]:
    if not isinstance(records, list) or len(records) > MAX_RECORDS:
        raise InputError("Invalid or oversized record inventory")
    result, total = {}, 0
    for row in records:
        if not isinstance(row, dict) or set(row) != RECORD_KEYS:
            raise InputError("Invalid normalized record fields")
        key = identifier(row["id"])
        if key in result:
            raise InputError("Duplicate normalized record ID")
        sequence = row["sequence"]
        if not isinstance(sequence, str) or not sequence or not set(sequence) <= DNA:
            raise InputError("Invalid normalized DNA sequence")
        integer(row["length"], "record length", 1, MAX_BASES)
        if row["length"] != len(sequence) or row["sha256"] != digest(sequence):
            raise InputError("Normalized record length/hash mismatch")
        if row["assembly_id"] is not None:
            identifier(row["assembly_id"], "assembly_id")
        if row["topology"] not in {"linear", "circular"}:
            raise InputError("Invalid topology")
        total += len(sequence)
        if total > MAX_BASES:
            raise LimitError("Record inventory exceeds base budget")
        result[key] = row
    return result


def interval(value: object, size: int) -> tuple[int, int]:
    if not isinstance(value, dict):
        raise InputError("Missing source interval")
    start = integer(value.get("start"), "interval start", 1, size)
    end = integer(value.get("end"), "interval end", start, size)
    if "length" in value and (type(value["length"]) is not int or value["length"] != end - start + 1):
        raise InputError("Interval length mismatch")
    return start, end


def validate_inventory(data: object) -> dict:
    if not isinstance(data, dict) or set(data) != {"schema", "sources", "arrays", "queries", "annotated_records", "provenance"}:
        raise InputError("Invalid inventory fields")
    if data["schema"] != SCHEMA:
        raise InputError("Unsupported inventory schema")
    sources = validate_records(data["sources"])
    annotated = data["annotated_records"]
    if not isinstance(annotated, list) or any(not isinstance(x, str) for x in annotated):
        raise InputError("Invalid annotation coverage")
    if len(set(annotated)) != len(annotated) or not set(annotated) <= set(sources):
        raise InputError("Unknown/duplicate annotation source")
    if not isinstance(data["arrays"], list) or not isinstance(data["queries"], list):
        raise InputError("Array/query inventory must be lists")
    if len(data["queries"]) > MAX_RECORDS or len(data["arrays"]) > MAX_RECORDS:
        raise LimitError("Too many array/query occurrences")
    arrays = {}
    native_ids = set()
    for row in data["arrays"]:
        if not isinstance(row, dict) or set(row) != ARRAY_KEYS:
            raise InputError("Invalid array fields")
        key = identifier(row["id"])
        native = identifier(row["native_id"])
        source = row["source_id"]
        if source not in sources or source not in annotated or key in arrays or (source, native) in native_ids:
            raise InputError("Invalid/duplicate array source or identity")
        interval(row, sources[source]["length"])
        if row["category"] not in CATEGORIES or row["orientation"] not in {"+", "-", None}:
            raise InputError("Invalid category/orientation")
        arrays[key] = row
        native_ids.add((source, native))
    query_ids, slots = set(), set()
    total_bases = 0
    for row in data["queries"]:
        if not isinstance(row, dict) or set(row) != QUERY_KEYS:
            raise InputError("Invalid query fields")
        key = identifier(row["id"])
        if key in query_ids:
            raise InputError("Duplicate query occurrence ID")
        query_ids.add(key)
        if row["kind"] not in {"observed", "deletion"}:
            raise InputError("Invalid spacer kind")
        seq = row["sequence"]
        if row["kind"] == "observed":
            if not isinstance(seq, str) or not seq or not set(seq) <= DNA or len(seq) > 10_000:
                raise InputError("Observed spacers require 1–10,000 IUPAC DNA bases")
            if row["sha256"] != digest(seq) or row["boundary"] is not None:
                raise InputError("Query hash/boundary mismatch")
            total_bases += len(seq)
        elif seq is not None or row["sha256"] is not None or row["start"] is not None or row["end"] is not None:
            raise InputError("Deletions have no sequence or inclusive interval")
        if row["source_id"] is None:
            if row["kind"] != "observed" or any(row[k] is not None for k in
                ("array_id", "ordinal", "start", "end", "boundary", "category", "source_orientation")):
                raise InputError("Unbound spacer FASTA cannot claim array context")
            continue
        source = sources.get(row["source_id"])
        array = arrays.get(row["array_id"])
        if source is None or array is None or array["source_id"] != row["source_id"]:
            raise InputError("Query source/array missing or mismatched")
        ordinal = integer(row["ordinal"], "spacer ordinal", 1, MAX_RECORDS)
        if (row["array_id"], ordinal) in slots:
            raise InputError("Duplicate source spacer ordinal")
        slots.add((row["array_id"], ordinal))
        if row["category"] != array["category"] or row["source_orientation"] != array["orientation"]:
            raise InputError("Query category/orientation disagrees with source array")
        if row["kind"] == "observed":
            start, end = interval(row, source["length"])
            if start < array["start"] or end > array["end"] or seq != source["sequence"][start-1:end]:
                raise InputError("Spacer does not reconstruct from its source interval")
        else:
            integer(row["boundary"], "deletion boundary", array["start"] - 1, array["end"])
    if total_bases > MAX_BASES:
        raise LimitError("Query inventory exceeds base budget")
    if not isinstance(data["provenance"], list):
        raise InputError("Invalid inventory provenance")
    for row in data["provenance"]:
        if not isinstance(row, dict) or set(row) != {"kind", "sha256", "source_id", "schema_version"}:
            raise InputError("Invalid provenance fields")
        if row["kind"] not in {"spacer_fasta", "identify_report"}:
            raise InputError("Invalid provenance kind")
        if not isinstance(row["sha256"], str) or len(row["sha256"]) != 64 or set(row["sha256"]) - set("0123456789abcdef"):
            raise InputError("Invalid provenance digest")
        if row["source_id"] is not None and row["source_id"] not in sources:
            raise InputError("Unknown provenance source")
        if row["schema_version"] not in {None, "1.2.0"}:
            raise InputError("Unsupported source schema")
    return data


def from_spacers(records: list[dict], raw_hash: str) -> dict:
    data = {"schema": SCHEMA, "sources": [], "arrays": [], "queries": [], "annotated_records": [],
            "provenance": [{"kind": "spacer_fasta", "sha256": raw_hash, "source_id": None, "schema_version": None}]}
    for record in records:
        query = dict.fromkeys(QUERY_KEYS)
        query.update(id=record["id"], sequence=record["sequence"], sha256=record["sha256"], kind="observed")
        data["queries"].append(query)
    return validate_inventory(data)


def from_identify(sources: list[dict], reports: list[Path]) -> dict:
    source_map = validate_records(sources)
    data = {"schema": SCHEMA, "sources": sources, "arrays": [], "queries": [], "annotated_records": [], "provenance": []}
    for path in reports:
        report_bytes = read_bytes(path)
        report = decode_json(report_bytes)
        if not isinstance(report, dict) or report.get("schema") != {"name": "CRISPRidentify-report", "version": "1.2.0"}:
            raise InputError("prepare supports canonical CRISPRidentify-report 1.2.0 only")
        coords = report.get("coordinate_system", {})
        if not isinstance(coords, dict) or type(coords.get("indexing")) is not int or (coords.get("indexing"), coords.get("end_inclusive"), coords.get("reference"),
            coords.get("orientation_is_separate")) != (1, True, "input_forward_strand", True):
            raise InputError("Unsupported report coordinate frame")
        if coords["end_inclusive"] is not True or coords["orientation_is_separate"] is not True:
            raise InputError("Coordinate flags must be boolean")
        claim = report.get("source", {})
        if not isinstance(claim, dict):
            raise InputError("Invalid report source")
        source = source_map.get(claim.get("id"))
        if source is None or source["id"] in data["annotated_records"]:
            raise InputError("Unknown or duplicate source report")
        if type(claim.get("length")) is not int or claim.get("length") != source["length"] or claim.get("sha256") != source["sha256"] or claim.get("sha256_scope") != "uppercase_sequence":
            raise InputError("CRISPRidentify source binding mismatch")
        rows = report.get("arrays")
        if not isinstance(rows, list) or len(rows) > MAX_RECORDS:
            raise InputError("Invalid array list")
        data["annotated_records"].append(source["id"])
        data["provenance"].append({"kind": "identify_report", "source_id": source["id"],
                                   "schema_version": "1.2.0", "sha256": digest(report_bytes)})
        for native in rows:
            if len(data["arrays"]) >= MAX_RECORDS:
                raise LimitError("Too many source array occurrences")
            if not isinstance(native, dict):
                raise InputError("Invalid native array")
            start, end = interval(native.get("source_interval"), source["length"])
            if native.get("source_sequence", source["sequence"][start-1:end]) != source["sequence"][start-1:end]:
                raise InputError("Array sequence is not source-bound")
            orientation = native.get("orientation")
            if not isinstance(orientation, dict):
                raise InputError("Missing orientation object")
            array_id = f"array-{len(data['arrays'])+1:06d}"
            array = {"id": array_id, "native_id": native.get("id"), "source_id": source["id"],
                     "start": start, "end": end, "category": native.get("category"), "orientation": orientation.get("strand")}
            data["arrays"].append(array)
            repeats = native.get("repeats")
            spacers = native.get("spacers")
            deletions = native.get("spacer_deletions", [])
            if not all(isinstance(x, list) for x in (repeats, spacers, deletions)) or len(repeats) < 2:
                raise InputError("Missing canonical repeat/spacer architecture")
            counts = {"repeat_count": len(repeats), "spacer_count": len(spacers)+len(deletions)}
            for key, count in counts.items():
                if type(native.get(key)) is not int or native[key] != count:
                    raise InputError("Native array count mismatch")
            if len(repeats) != len(spacers) + len(deletions) + 1:
                raise InputError("Repeat/spacer architecture count mismatch")
            for key, count in (("non_deleted_spacer_count", len(spacers)), ("spacer_deletion_count", len(deletions))):
                if key in native and (type(native[key]) is not int or native[key] != count):
                    raise InputError("Native spacer count mismatch")
            repeat_intervals = []
            for ordinal, repeat in enumerate(repeats, 1):
                if not isinstance(repeat, dict) or type(repeat.get("index")) is not int or repeat["index"] != ordinal:
                    raise InputError("Invalid repeat ordinal")
                left, right = interval(repeat.get("source_interval"), source["length"])
                if repeat.get("sequence") != source["sequence"][left-1:right] or left < start or right > end:
                    raise InputError("Repeat source reconstruction mismatch")
                if repeat_intervals and left <= repeat_intervals[-1][1]:
                    raise InputError("Repeat intervals overlap or are unordered")
                repeat_intervals.append((left, right))
            if repeat_intervals[0][0] != start or repeat_intervals[-1][1] != end:
                raise InputError("Array interval does not match outer repeats")
            seen = set()
            for kind, units in (("observed", spacers), ("deletion", deletions)):
                for unit in units:
                    if len(data["queries"]) >= MAX_RECORDS:
                        raise LimitError("Too many source spacer occurrences")
                    if not isinstance(unit, dict):
                        raise InputError("Invalid spacer entry")
                    ordinal = integer(unit.get("index"), "native spacer ordinal", 1, len(repeats)-1)
                    if ordinal in seen:
                        raise InputError("Duplicate spacer ordinal")
                    seen.add(ordinal)
                    left, right = repeat_intervals[ordinal-1][1]+1, repeat_intervals[ordinal][0]-1
                    query = {"id": f"spacer-{len(data['queries'])+1:06d}", "kind": kind, "source_id": source["id"],
                             "array_id": array_id, "ordinal": ordinal, "category": array["category"],
                             "source_orientation": array["orientation"], "boundary": None,
                             "start": None, "end": None, "sequence": None, "sha256": None}
                    if kind == "observed":
                        if interval(unit.get("source_interval"), source["length"]) != (left, right):
                            raise InputError("Spacer does not occupy its inter-repeat gap")
                        seq = unit.get("sequence")
                        if seq != source["sequence"][left-1:right]:
                            raise InputError("Spacer source reconstruction mismatch")
                        query.update(sequence=seq, sha256=digest(seq), start=left, end=right)
                    else:
                        boundary = unit.get("source_boundary", {})
                        if right != left-1 or unit.get("sequence") != "" or boundary != {
                            "position": right, "coordinate_system": "0-based interbase", "reference": "input_forward_strand"}:
                            raise InputError("Invalid deletion boundary")
                        query["boundary"] = right
                    data["queries"].append(query)
    data["annotated_records"].sort()
    return validate_inventory(data)
