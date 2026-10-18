#!/usr/bin/env python3
"""Offline, dependency-free normalization of checksum-pinned CRISPRcasDb exports.

The same source is sent over SSH to transform exports next to their storage.
It never modifies scientific sources or starts services. Compressed SQLite shards
retain the published scalar fields; private paths and native JSON are omitted.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import tarfile

STRIDE = 1_000_000_000
CAMPAIGN_PROFILE = "5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540"
MAX_DATABASE = 512 * 1024**2
MAX_SHARD = 256 * 1024**2
TABLES = {
    "assemblies": "assembly_pk INTEGER PRIMARY KEY, accession TEXT NOT NULL UNIQUE, status TEXT NOT NULL, validation_status TEXT, sequence_count INTEGER, array_count INTEGER, cas_system_count INTEGER, cas_gene_count INTEGER",
    "sequences": "sequence_pk INTEGER PRIMARY KEY, assembly_pk INTEGER NOT NULL REFERENCES assemblies(assembly_pk), sequence_id TEXT NOT NULL, header TEXT, length INTEGER, sha256 TEXT, sequence_ordinal INTEGER",
    "arrays": "array_pk INTEGER PRIMARY KEY, sequence_pk INTEGER NOT NULL REFERENCES sequences(sequence_pk), tool TEXT, tool_version TEXT, call_id TEXT, reported_sequence_id TEXT, start INTEGER, end INTEGER, length INTEGER, strand TEXT, repeat_consensus TEXT, repeat_count INTEGER, spacer_count INTEGER, evidence_level TEXT, category TEXT, certainty_score REAL",
    "repeats": "repeat_pk INTEGER PRIMARY KEY, array_pk INTEGER NOT NULL REFERENCES arrays(array_pk), ordinal INTEGER, start INTEGER, end INTEGER, length INTEGER, sequence TEXT, reported_sequence_id TEXT",
    "spacers": "spacer_pk INTEGER PRIMARY KEY, array_pk INTEGER NOT NULL REFERENCES arrays(array_pk), ordinal INTEGER, kind TEXT, start INTEGER, end INTEGER, length INTEGER, boundary_position INTEGER, sequence TEXT, reported_sequence_id TEXT",
    "cas_systems": "cas_system_pk INTEGER PRIMARY KEY, sequence_pk INTEGER NOT NULL REFERENCES sequences(sequence_pk), tool TEXT, tool_version TEXT, call_id TEXT, reported_sequence_id TEXT, start INTEGER, end INTEGER, length INTEGER, cas_class TEXT, cas_type TEXT, subtype TEXT, confidence REAL, confidence_kind TEXT, experimental INTEGER, missing_gene_count INTEGER",
    "cas_genes": "cas_gene_pk INTEGER PRIMARY KEY, sequence_pk INTEGER NOT NULL REFERENCES sequences(sequence_pk), cas_system_pk INTEGER REFERENCES cas_systems(cas_system_pk), tool TEXT, tool_version TEXT, call_id TEXT, reported_sequence_id TEXT, start INTEGER, end INTEGER, length INTEGER, strand TEXT, annotation TEXT, profile TEXT, cas_class TEXT, cas_type TEXT, subtype TEXT, score REAL, score_kind TEXT, protein_length INTEGER, translation_table INTEGER, is_cas INTEGER",
}
COLUMNS = {table: [field.strip().split()[0] for field in definition.split(",")] for table, definition in TABLES.items()}
INDEXES = {
    "sequences": [("assembly_pk",)],
    # Scientific searches are genome-scoped: parent-key indexes bound each
    # search without storing redundant global sequence/profile indexes.
    "arrays": [("sequence_pk",)],
    "repeats": [("array_pk", "ordinal")],
    "spacers": [("array_pk", "ordinal")],
    "cas_systems": [("sequence_pk",)],
    "cas_genes": [("sequence_pk",), ("cas_system_pk",)],
}


def create_tables(connection):
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA user_version=1")
    connection.execute("CREATE TABLE export_metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)")
    for table, definition in TABLES.items():
        connection.execute(f"CREATE TABLE {table} ({definition})")


def add_indexes(connection):
    for table, groups in INDEXES.items():
        for columns in groups:
            name = "idx_" + table + "_" + "_".join(columns)
            connection.execute(f"CREATE INDEX {name} ON {table} ({','.join(columns)})")


def checked_archive(path, expected_sha256):
    """Read named members only, never extract archive paths to disk."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("source archive is not a regular file")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
    if digest.hexdigest() != expected_sha256:
        raise ValueError("source bundle checksum mismatch")
    needed = {"results.sqlite3": MAX_DATABASE, "export_manifest.json": 1024**2,
              "collection_index.json": 8 * 1024**2, "analysis-profile.json": 1024**2}
    members = {}
    with tarfile.open(path, "r|gz") as archive:
        for member in archive:
            if member.name not in needed:
                continue
            if member.name in members or not member.isfile() or member.size > needed[member.name]:
                raise ValueError("invalid or oversized export member")
            handle = archive.extractfile(member)
            data = handle.read(needed[member.name] + 1)
            if len(data) != member.size:
                raise ValueError("truncated export member")
            members[member.name] = data
    if members.keys() != needed.keys():
        raise ValueError("incomplete source publication")
    manifest = json.loads(members["export_manifest.json"])
    record = manifest["artifacts"]["results.sqlite3"]
    data = members["results.sqlite3"]
    if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise ValueError("source database checksum mismatch")
    if manifest.get("database_schema_version") != 1 or manifest["counts"]["export_issues"]:
        raise ValueError("unsupported or scientifically invalid export")
    profile = json.loads(members["analysis-profile.json"])
    profile_digest = hashlib.sha256(json.dumps(profile.get("analysis_profile"), sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if profile.get("analysis_profile_sha256") != CAMPAIGN_PROFILE or profile_digest != CAMPAIGN_PROFILE:
        raise ValueError("source publication has a different scientific profile")
    collection = json.loads(members["collection_index.json"])
    if collection.get("analysis_profile_sha256") != CAMPAIGN_PROFILE:
        raise ValueError("source collection has a different scientific profile")
    return data, manifest, members


def normalize_database(data, shard_id, expected_counts=None, selected_accessions=None):
    """Return a self-contained compressed shard, assembly registry rows and counts."""
    if not 0 < shard_id < 1_000_000:
        raise ValueError("invalid shard id")
    source = sqlite3.connect(":memory:")
    target = sqlite3.connect(":memory:")
    try:
        source.deserialize(data)
        source.execute("PRAGMA query_only=ON")
        source.execute("PRAGMA trusted_schema=OFF")
        source.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 8 * 1024**2)
        kinds = dict(source.execute("SELECT name,type FROM sqlite_master"))
        if any(kinds.get(table) != "table" for table in TABLES):
            raise ValueError("source tables are missing or substituted")
        if expected_counts is not None:
            for table in TABLES:
                if source.execute(f"SELECT count(*) FROM {table}").fetchone()[0] != expected_counts[table]:
                    raise ValueError("source row counts disagree with manifest")
        create_tables(target)
        offset = shard_id * STRIDE
        allowed_assemblies = None
        allowed_sequences = set()
        allowed_arrays = set()
        if selected_accessions is not None:
            allowed_assemblies = {row[0] for row in source.execute("SELECT assembly_pk,accession FROM assemblies") if row[1] in selected_accessions}
        counts = {}
        registry = []
        for table, columns in COLUMNS.items():
            indices = [i for i, col in enumerate(columns) if col.endswith("_pk")]
            count = 0
            insert = f"INSERT INTO {table} ({','.join(columns)}) VALUES ({','.join('?' for _ in columns)})"
            selection = ""
            if table == "sequences":
                # Preserve every feature's genomic context. Feature-free contigs
                # have no browsable records; their total remains on the assembly.
                selection = " WHERE sequence_pk IN (SELECT sequence_pk FROM arrays UNION SELECT sequence_pk FROM cas_systems UNION SELECT sequence_pk FROM cas_genes)"
            for original in source.execute(f"SELECT {','.join(columns)} FROM {table}{selection} ORDER BY {columns[0]}"):
                row = list(original)
                if allowed_assemblies is not None:
                    if table == "assemblies" and row[0] not in allowed_assemblies:
                        continue
                    if table == "sequences" and row[1] not in allowed_assemblies:
                        continue
                    if table in {"arrays", "cas_systems", "cas_genes"} and row[1] not in allowed_sequences:
                        continue
                    if table in {"repeats", "spacers"} and row[1] not in allowed_arrays:
                        continue
                    if table == "sequences":
                        allowed_sequences.add(row[0])
                    if table == "arrays":
                        allowed_arrays.add(row[0])
                for i in indices:
                    if row[i] is not None:
                        if not isinstance(row[i], int) or not 0 < row[i] < STRIDE:
                            raise ValueError("invalid source key")
                        row[i] += offset
                for value in row:
                    if isinstance(value, str) and (len(value) > 65536 or "\x00" in value):
                        raise ValueError("oversized or invalid scalar")
                if table == "assemblies" and row[2] != "complete":
                    # Failure/unavailability is not a scientifically measured zero.
                    row[4:] = [None] * 4
                target.execute(insert, row)
                if table == "assemblies":
                    registry.append(row)
                count += 1
            counts[table] = count
        if target.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise ValueError("source relations are broken")
        add_indexes(target)
        target.commit()
        normalized = target.serialize()
        if len(normalized) > MAX_SHARD:
            raise ValueError("normalized shard exceeds memory budget")
        compressed = gzip.compress(normalized, compresslevel=6, mtime=0)
        return compressed, registry, counts, len(normalized)
    finally:
        source.close()
        target.close()


def transform_bundle(path, digest, shard_id, selected_accessions=None):
    data, manifest, members = checked_archive(path, digest)
    compressed, registry, counts, size = normalize_database(data, shard_id, manifest["counts"], selected_accessions)
    receipt = {
        "id": shard_id, "sha256": hashlib.sha256(compressed).hexdigest(),
        "bytes": len(compressed), "uncompressed_bytes": size, "counts": counts,
        "assemblies": registry, "source_sha256": digest,
        "source_database_sha256": manifest["artifacts"]["results.sqlite3"]["sha256"],
        "collection_index_sha256": hashlib.sha256(members["collection_index.json"]).hexdigest(),
        "analysis_profile_file_sha256": hashlib.sha256(members["analysis-profile.json"]).hexdigest(),
    }
    return receipt, compressed


def emit_frame(output, receipt, data):
    header = json.dumps(receipt, separators=(",", ":")).encode() + b"\n"
    output.write(header)
    output.write(data)
    output.flush()


def read_frame(stream):
    header = stream.readline(8 * 1024**2 + 1)
    if not header:
        return None
    if not header.endswith(b"\n") or len(header) > 8 * 1024**2:
        raise ValueError("invalid transfer header")
    receipt = json.loads(header)
    size = receipt["bytes"]
    if not isinstance(size, int) or not 0 < size <= 64 * 1024**2:
        raise ValueError("invalid transfer size")
    content = io.BytesIO()
    while content.tell() < size:
        block = stream.read(min(1024**2, size - content.tell()))
        if not block:
            raise ValueError("truncated transfer")
        content.write(block)
    data = content.getvalue()
    if hashlib.sha256(data).hexdigest() != receipt["sha256"]:
        raise ValueError("transfer checksum mismatch")
    return receipt, data
