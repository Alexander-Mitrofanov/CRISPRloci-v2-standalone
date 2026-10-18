"""Read checksum-pinned campaign sources without importing an upstream runtime.

Private source paths are inputs only. Public metadata contains content identities,
never a copy of a campaign configuration, native report, or connection details.
"""
from __future__ import annotations

from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tarfile
import tempfile

from .io import InputError, decode_json

DIGEST = re.compile(r"[0-9a-f]{64}\Z")
ACCESSION = re.compile(r"GC[AF]_\d+\.\d+\Z")
MAX_METADATA = 64 * 1024**2
MAX_DATABASE = 1024**3
MAX_RAW = 256 * 1024**2
MAX_SQLITE_VALUE = 8 * 1024**2
SOURCE_TABLES = {
    "assemblies": {"assembly_pk", "accession", "status", "validation_status", "sequence_count", "array_count", "metadata_json"},
    "sequences": {"sequence_pk", "assembly_pk", "sequence_id", "length", "sha256"},
    "arrays": {"array_pk", "sequence_pk", "tool", "tool_version", "call_id", "start", "end", "strand", "repeat_consensus", "repeat_count", "spacer_count", "category", "evidence_level", "certainty_score"},
    "repeats": {"repeat_pk", "array_pk", "ordinal", "start", "end", "sequence"},
}


def _preflight_sqlite_fields(connection, table, columns):
    """Bound projected row payload in SQLite before Python materializes values.

    Identifiers come only from this module's allowlists. CAST AS BLOB measures
    encoded bytes, including embedded NULs and multibyte text. Numeric lengths
    are the only values returned to Python by this preflight.
    """
    payload = "+".join(f'coalesce(length(CAST("{column}" AS BLOB)),0)' for column in sorted(columns))
    if connection.execute(f'SELECT 1 FROM "{table}" WHERE ({payload})>? LIMIT 1',
                          (MAX_SQLITE_VALUE,)).fetchone() is not None:
        raise InputError("Source SQLite projected row exceeds its byte limit")


def _configure_sqlite_length_limit(connection):
    setter = getattr(connection, "setlimit", None)
    if callable(setter):
        setter(sqlite3.SQLITE_LIMIT_LENGTH, MAX_SQLITE_VALUE)
        return False
    # Python3.10 lacks Connection.setlimit. The immutable, read-only database
    # cannot change between these SQL-only checks and subsequent value fetches.
    _preflight_sqlite_fields(connection, "sqlite_master", {"type", "name", "tbl_name", "sql"})
    return True


def sha256_file(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""):
            value.update(block)
    return value.hexdigest()


def required_digest(value, label="checksum"):
    if not isinstance(value, str) or DIGEST.fullmatch(value) is None:
        raise InputError(f"Invalid {label}")
    return value


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     allow_nan=False, separators=(",", ":")).encode()).hexdigest()


def checked_path(path):
    result = Path(path).absolute()
    if any(parent.is_symlink() for parent in (result, *result.parents)) or not result.is_file():
        raise InputError("Source must be a regular file without symlink traversal")
    return result


def _source_json(data):
    result = decode_json(data)
    if not isinstance(result, dict):
        raise InputError("Source JSON must be an object")
    return result


def bounded_json(path, maximum=MAX_RAW):
    path = checked_path(path)
    if path.stat().st_size > maximum:
        raise InputError("Compressed source exceeds its byte limit")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as stream:
        data = stream.read(maximum + 1)
    if len(data) > maximum:
        raise InputError("Expanded source exceeds its byte limit")
    return _source_json(data)


def validate_free_reserve(value):
    if type(value) is not int or not 0 <= value <= 2**63 - 1:
        raise InputError("Free-space byte budgets must be nonnegative integers below 2**63")
    return value


def check_free(directory, minimum, extra=0):
    validate_free_reserve(minimum)
    validate_free_reserve(extra)
    if shutil.disk_usage(directory).free < minimum + extra:
        raise InputError("Import stopped before exhausting the configured free-space reserve")


def source_selection(source):
    """Return selected accessions, optionally bound to source document digests."""
    entries = source.get("entries")
    if entries is not None:
        if not isinstance(entries, list):
            raise InputError("Selected source entries must be a list")
        result = {}
        for entry in entries:
            if not isinstance(entry, dict):
                raise InputError("Malformed selected source entry")
            accession = entry.get("accession")
            if not isinstance(accession, str) or not ACCESSION.fullmatch(accession) or accession in result:
                raise InputError("Invalid or duplicate selected accession")
            if entry.get("status") not in {"complete", "failed", "unavailable"}:
                raise InputError("Selected source outcome is not terminal")
            result[accession] = {"accession": accession, "status": entry["status"],
                                 "sha256": required_digest(entry.get("sha256"), "source document checksum")}
        return result
    values = source.get("selected")
    if values is None:
        return None
    if not isinstance(values, list):
        raise InputError("Selected accessions must be a unique list")
    if any(not isinstance(value, str) or not ACCESSION.fullmatch(value) for value in values):
        raise InputError("Invalid selected accession")
    if len(values) != len(set(values)):
        raise InputError("Selected accessions must be a unique list")
    return {value: {"accession": value} for value in values}


def source_identity(source):
    if not isinstance(source, dict):
        raise InputError("Malformed source plan entry")
    kind = source.get("kind", source.get("source_kind"))
    if kind not in {"bundle", "raw", "unavailable"}:
        raise InputError("Unsupported source kind")
    selected = source_selection(source)
    selection = None if selected is None else [selected[key] for key in sorted(selected)]
    selection_sha = canonical_digest(selection)
    sha = required_digest(source.get("sha256"), "source checksum")
    return {"source_key": canonical_digest({"kind": kind, "sha256": sha, "selection_sha256": selection_sha}),
            "kind": kind, "source_sha256": sha, "selection_sha256": selection_sha}


def plan_sources(plan):
    if not isinstance(plan, dict):
        raise InputError("Source plan must be an object")
    if "sources" in plan:
        values = plan["sources"]
    elif "jobs" in plan:
        values = plan["jobs"]
    elif "nodes" in plan:
        if (not isinstance(plan["nodes"], dict)
                or any(not isinstance(group, list) for group in plan["nodes"].values())):
            raise InputError("Malformed source plan nodes")
        values = [item for group in plan["nodes"].values() for item in group]
    else:
        raise InputError("Source plan has no sources")
    if not isinstance(values, list):
        raise InputError("Source plan entries must be a list")
    for value in values:
        if not isinstance(value, dict):
            raise InputError("Malformed source plan entry")
        if value.get("kind") == "raw_group":
            items = value.get("items")
            if not isinstance(items, list):
                raise InputError("Malformed raw source group")
            for item in items:
                source_identity(item)
                yield item
        else:
            source_identity(value)
            yield value


def taxonomy(metadata):
    """Allowlist source taxon identifiers and scientific organism names only."""
    if isinstance(metadata, str):
        if len(metadata) > MAX_METADATA:
            raise InputError("Source assembly metadata exceeds its byte limit")
        try:
            metadata = decode_json(metadata.encode("utf-8"))
        except (InputError, UnicodeError) as error:
            raise InputError("Malformed source assembly metadata") from error
    if not isinstance(metadata, dict):
        return None, None
    organism = metadata.get("organism")
    organism = organism if isinstance(organism, dict) else {}
    taxon = organism.get("taxId", organism.get("tax_id", metadata.get("taxId", metadata.get("tax_id"))))
    if isinstance(taxon, str) and taxon.isdecimal():
        # Bound conversion even when metadata contains thousands of digits;
        # an out-of-range identifier is missing taxonomy, not an import error.
        significant = taxon.lstrip("0")
        taxon = int(significant or "0") if len(significant) <= 19 else None
    if type(taxon) is not int or not 0 < taxon < 2**63:
        taxon = None
    name = organism.get("organismName", organism.get("scientificName", metadata.get("organism_name")))
    if name is None and isinstance(metadata.get("organism"), str):
        name = metadata["organism"]
    if not isinstance(name, str) or not 1 <= len(name) <= 500 or any(ord(c) < 32 for c in name):
        name = None
    # Names are fields, not instructions or arbitrary native metadata.
    return taxon, name


class BundleSource:
    def __init__(self, connection, selected, metadata, expected_counts):
        self.connection = connection
        self.selected = selected
        self.metadata = metadata
        self.expected_counts = expected_counts
        self.selected_pks = set()
        self.expected_arrays = 0
        actual = set()
        for row in connection.execute("SELECT assembly_pk,accession,status,validation_status,array_count FROM assemblies"):
            pk, accession, status, validation, arrays = row
            if accession not in selected:
                continue
            entry = selected[accession]
            if status != entry["status"] or (status == "complete" and validation != "valid"):
                raise InputError("Source assembly outcome or validation differs from its collection")
            actual.add(accession)
            self.selected_pks.add(pk)
            self.expected_arrays += arrays or 0
        if actual != set(selected):
            raise InputError("Source database does not contain every selected accession")

    def assemblies(self):
        for row in self.connection.execute("SELECT assembly_pk,accession,status,validation_status,sequence_count,metadata_json FROM assemblies ORDER BY assembly_pk"):
            if row[0] in self.selected_pks:
                yield {"native_pk": row[0], "accession": row[1], "status": row[2],
                       "validation_status": row[3], "sequence_count": row[4],
                       "source_document_sha256": self.selected[row[1]]["sha256"], "taxonomy": taxonomy(row[5])}

    def contigs(self):
        sql = "SELECT sequence_pk,assembly_pk,sequence_id,sha256,length FROM sequences WHERE sequence_pk IN (SELECT sequence_pk FROM arrays) ORDER BY sequence_pk"
        for row in self.connection.execute(sql):
            if row[1] in self.selected_pks:
                yield {"native_pk": row[0], "assembly_pk": row[1], "sequence_id": row[2], "sequence_sha256": row[3], "length": row[4]}

    def arrays(self):
        fields = ["array_pk", "sequence_pk", "call_id", "tool", "tool_version", "start", "end", "strand", "repeat_consensus", "repeat_count", "spacer_count", "category", "evidence_level", "certainty_score"]
        sql = "SELECT " + ",".join("r." + field for field in fields) + ",s.assembly_pk FROM arrays r JOIN sequences s ON s.sequence_pk=r.sequence_pk ORDER BY r.array_pk"
        for row in self.connection.execute(sql):
            if row[-1] in self.selected_pks:
                result = dict(zip(fields, row[:-1]))
                result["native_pk"] = result.pop("array_pk")
                yield result

    def repeats(self):
        sql = "SELECT r.repeat_pk,r.array_pk,r.ordinal,r.start,r.end,r.sequence,s.assembly_pk FROM repeats r JOIN arrays a ON a.array_pk=r.array_pk JOIN sequences s ON s.sequence_pk=a.sequence_pk ORDER BY r.array_pk,r.ordinal,r.repeat_pk"
        for row in self.connection.execute(sql):
            if row[-1] in self.selected_pks:
                yield dict(zip(("native_pk", "array_pk", "ordinal", "start", "end", "sequence"), row[:-1]))


def _raw_object(value, label):
    if not isinstance(value, dict):
        raise InputError(f"Raw source {label} must be an object")
    return value


def _raw_list(value, label):
    if not isinstance(value, list):
        raise InputError(f"Raw source {label} must be a list")
    return value


class RawSource:
    def __init__(self, payload, source, expected_profile):
        self.payload = payload
        self.kind = source_identity(source)["kind"]
        self.metadata = {**source_identity(source), "database_sha256": None,
                         "collection_index_sha256": None, "profile_sha256": expected_profile,
                         "profile_file_sha256": None}
        selected = source_selection(source)
        if not selected or any(not {"status", "sha256"} <= set(entry) for entry in selected.values()):
            raise InputError("Raw sources require explicit selected source entries")
        self.selected = selected
        if self.kind == "unavailable":
            unavailable = payload.get("unavailable")
            if not isinstance(unavailable, dict) or not set(selected) <= set(unavailable):
                raise InputError("Unavailable receipt does not contain every selected accession")
            if any(entry.get("status") != "unavailable" or entry["sha256"] != source["sha256"] for entry in selected.values()):
                raise InputError("Unavailable outcome identity differs from its receipt")
            self.expected_arrays = 0
            return
        assembly = _raw_object(payload.get("assembly", {}), "assembly")
        accession = assembly.get("accession")
        if set(selected) != {accession}:
            raise InputError("Raw assembly differs from selected accession")
        entry = selected[accession]
        if entry["status"] != payload.get("status") or entry["sha256"] != source["sha256"]:
            raise InputError("Raw assembly identity differs from its selected outcome")
        provenance = _raw_object(payload.get("provenance", {}), "provenance")
        if provenance.get("analysis_profile_sha256") != expected_profile:
            raise InputError("Raw source scientific profile differs")
        if "analysis_profile" in provenance and canonical_digest(provenance["analysis_profile"]) != expected_profile:
            raise InputError("Raw source scientific profile content differs from its digest")
        validation = _raw_object(payload.get("validation", {}), "validation")
        if payload["status"] == "complete":
            _raw_list(payload.get("sequences", []), "contigs")
        if payload["status"] == "complete" and validation.get("status") != "valid":
            raise InputError("Raw successful assembly was not source validated")
        summary = _raw_object(payload.get("summary", {}), "summary")
        # Different preserved assembly schema versions use either key.
        self.expected_arrays = summary.get("array_count", summary.get("arrays"))
        if self.expected_arrays is None and isinstance(summary.get("tools"), dict):
            values = [value.get("arrays") for value in summary["tools"].values() if isinstance(value, dict)]
            if values and all(type(value) is int and value >= 0 for value in values):
                self.expected_arrays = sum(values)

    def assemblies(self):
        if self.kind == "unavailable":
            for index, accession in enumerate(sorted(self.selected), 1):
                yield {"native_pk": index, "accession": accession, "status": "unavailable", "validation_status": None,
                       "sequence_count": None, "source_document_sha256": self.selected[accession]["sha256"], "taxonomy": (None, None)}
            return
        assembly = self.payload["assembly"]
        yield {"native_pk": 1, "accession": assembly["accession"], "status": self.payload["status"],
               "validation_status": self.payload.get("validation", {}).get("status"),
               "sequence_count": len(self.payload.get("sequences", [])) if self.payload["status"] == "complete" else None,
               "source_document_sha256": self.selected[assembly["accession"]]["sha256"],
               "taxonomy": taxonomy(assembly.get("ncbi_metadata"))}

    def _sequences(self):
        if self.kind == "unavailable" or self.payload.get("status") != "complete":
            return
        for index, sequence in enumerate(self.payload.get("sequences", []), 1):
            yield index, _raw_object(sequence, "contig")

    def _tools(self, sequence):
        for tool in _raw_list(sequence.get("tools", []), "tools"):
            tool = _raw_object(tool, "tool")
            yield tool, _raw_list(tool.get("arrays", []), "arrays")

    def contigs(self):
        for index, sequence in self._sequences():
            if any(arrays for _, arrays in self._tools(sequence)):
                yield {"native_pk": index, "assembly_pk": 1, "sequence_id": sequence.get("id"),
                       "sequence_sha256": sequence.get("sha256"), "length": sequence.get("length")}

    def _arrays(self):
        count = 0
        for sequence_pk, sequence in self._sequences():
            for tool, arrays in self._tools(sequence):
                for array in arrays:
                    count += 1
                    yield count, sequence_pk, tool, _raw_object(array, "array")

    def arrays(self):
        for index, sequence_pk, tool, array in self._arrays():
            interval = _raw_object(array.get("interval", {}), "array interval")
            confidence = _raw_object(array.get("confidence", {}), "array confidence")
            yield {"native_pk": index, "sequence_pk": sequence_pk, "call_id": array.get("call_id"),
                   "tool": array.get("tool", tool.get("name")), "tool_version": tool.get("version"),
                   "start": interval.get("start"), "end": interval.get("end"), "strand": array.get("strand"),
                   "repeat_consensus": array.get("repeat_consensus"), "repeat_count": array.get("repeat_count"),
                   "spacer_count": array.get("spacer_count"), "category": confidence.get("category"),
                   "evidence_level": confidence.get("evidence_level"), "certainty_score": confidence.get("score")}

    def repeats(self):
        index = 0
        for array_pk, _, _, array in self._arrays():
            for ordinal, repeat in enumerate(_raw_list(array.get("repeats", []), "repeats"), 1):
                index += 1
                repeat = _raw_object(repeat, "repeat")
                interval = _raw_object(repeat.get("interval", {}), "repeat interval")
                yield {"native_pk": index, "array_pk": array_pk, "ordinal": repeat.get("index", ordinal),
                       "start": interval.get("start"), "end": interval.get("end"), "sequence": repeat.get("sequence")}


@contextmanager
def open_source(source, expected_profile, scratch_dir, min_free_bytes=2 * 1024**3):
    identity = source_identity(source)
    required_digest(expected_profile, "scientific profile checksum")
    path = checked_path(source.get("path"))
    if sha256_file(path) != identity["source_sha256"]:
        raise InputError("Source file checksum mismatch")
    if source.get("bytes") is not None and path.stat().st_size != source["bytes"]:
        raise InputError("Source file size mismatch")
    if identity["kind"] != "bundle":
        yield RawSource(bounded_json(path), source, expected_profile)
        return
    scratch = Path(scratch_dir)
    check_free(scratch, min_free_bytes)
    descriptor, filename = tempfile.mkstemp(prefix=".crisprmap-source-", suffix=".sqlite3", dir=scratch)
    os.close(descriptor)
    temporary = Path(filename)
    connection = None
    try:
        captured, database_sha, database_size = {}, None, None
        limits = {"results.sqlite3": MAX_DATABASE, "export_manifest.json": MAX_METADATA,
                  "collection_index.json": MAX_METADATA, "analysis-profile.json": MAX_METADATA}
        with tarfile.open(path, "r|gz") as archive:
            for member in archive:
                if member.name not in limits:
                    continue
                if member.name in captured or not member.isfile() or not 0 <= member.size <= limits[member.name]:
                    raise InputError("Duplicate, unsafe or oversized source archive member")
                handle = archive.extractfile(member)
                if member.name == "results.sqlite3":
                    check_free(scratch, min_free_bytes, member.size)
                    value = hashlib.sha256()
                    remaining = member.size
                    with temporary.open("wb") as output:
                        while remaining:
                            block = handle.read(min(1024**2, remaining))
                            if not block:
                                raise InputError("Truncated source database")
                            output.write(block)
                            value.update(block)
                            remaining -= len(block)
                    database_sha, database_size = value.hexdigest(), member.size
                    captured[member.name] = True
                else:
                    raw = handle.read(limits[member.name] + 1)
                    if len(raw) != member.size:
                        raise InputError("Truncated source metadata")
                    captured[member.name] = raw
        if set(captured) != set(limits):
            raise InputError("Source bundle lacks required committed artifacts")
        manifest = _source_json(captured["export_manifest.json"])
        if (manifest.get("database_schema_version") != 1 or manifest.get("status") != "complete"
                or manifest.get("commit_marker") != "export_manifest.json"
                or manifest.get("counts", {}).get("export_issues") != 0):
            raise InputError("Source export is uncommitted, incompatible or contains validation issues")
        database = manifest.get("artifacts", {}).get("results.sqlite3", {})
        if database.get("sha256") != database_sha or database.get("bytes") != database_size:
            raise InputError("Source database does not match its export manifest")
        profile = _source_json(captured["analysis-profile.json"])
        if profile.get("analysis_profile_sha256") != expected_profile or canonical_digest(profile.get("analysis_profile")) != expected_profile:
            raise InputError("Source profile receipt or profile content differs from the pinned profile")
        collection = _source_json(captured["collection_index.json"])
        collection_sha = hashlib.sha256(captured["collection_index.json"]).hexdigest()
        if collection.get("analysis_profile_sha256") != expected_profile:
            raise InputError("Source collection profile differs")
        if source.get("collection_index_sha256") is not None and source["collection_index_sha256"] != collection_sha:
            raise InputError("Source collection checksum differs from the selected plan")
        entries = {}
        for entry in collection.get("entries", []):
            accession = entry.get("accession")
            if not isinstance(accession, str) or not ACCESSION.fullmatch(accession) or accession in entries:
                raise InputError("Invalid or duplicate source collection accession")
            if entry.get("status") not in {"complete", "failed"}:
                raise InputError("Nonterminal source collection outcome")
            entries[accession] = {"accession": accession, "status": entry["status"],
                                  "sha256": required_digest(entry.get("sha256"))}
        selection = source_selection(source)
        if selection is not None:
            if not set(selection) <= set(entries):
                raise InputError("Selected accessions absent from source collection")
            for accession, selected in selection.items():
                if any(key in selected and selected[key] != entries[accession][key] for key in ("status", "sha256")):
                    raise InputError("Selected source document differs from committed collection")
            entries = {accession: entries[accession] for accession in selection}
        connection = sqlite3.connect(temporary.as_uri() + "?mode=ro&immutable=1", uri=True)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        needs_field_preflight = _configure_sqlite_length_limit(connection)
        if connection.execute("PRAGMA user_version").fetchone()[0] != 1:
            raise InputError("Source database schema version differs from its manifest")
        for table, columns in SOURCE_TABLES.items():
            row = connection.execute("SELECT type,sql FROM sqlite_master WHERE name=?", (table,)).fetchone()
            if row is None or row[0] != "table" or "VIRTUAL" in row[1].upper():
                raise InputError("Source scientific table is missing or substituted")
            actual_columns = {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}
            if not columns <= actual_columns:
                raise InputError("Source scientific table schema is incompatible")
            if needs_field_preflight:
                _preflight_sqlite_fields(connection, table, columns)
        metadata = {**identity, "database_sha256": database_sha, "collection_index_sha256": collection_sha,
                    "profile_sha256": expected_profile,
                    "profile_file_sha256": hashlib.sha256(captured["analysis-profile.json"]).hexdigest()}
        yield BundleSource(connection, entries, metadata, manifest["counts"])
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)
