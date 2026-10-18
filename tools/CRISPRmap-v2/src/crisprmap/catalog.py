"""Resumable, provenance-preserving repeat catalogues built beside their sources.

A source is committed in one SQLite transaction. ``manifest.json`` is published
last and is the only completion marker. Build files are never valid references.
"""
from __future__ import annotations

from collections import OrderedDict
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile
import time
from datetime import datetime, timezone

from .core import reverse_complement
from .io import InputError, atomic_json
from .source_import import (ACCESSION, bounded_json, canonical_digest, check_free, checked_path,
                            open_source, plan_sources, required_digest, sha256_file, source_identity, validate_free_reserve)

SCHEMA = {"name": "crisprmap.catalog", "version": "1.0.0"}
STRIDE = 1_000_000_000
NORMALIZATION = "uppercase_U_to_T_canonical_reverse_complement"
NODE_COUNTS = ("bona_fide_consensus_count", "possible_consensus_count", "bona_fide_repeat_count", "possible_repeat_count")
SCHEMA_SQL = """
PRAGMA user_version=1;
CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
CREATE TABLE sources(
 source_id INTEGER PRIMARY KEY,source_key TEXT NOT NULL UNIQUE,kind TEXT NOT NULL,
 source_sha256 TEXT NOT NULL,database_sha256 TEXT,collection_index_sha256 TEXT,
 profile_sha256 TEXT NOT NULL,profile_file_sha256 TEXT,selection_sha256 TEXT NOT NULL,
 counts_json TEXT NOT NULL);
CREATE TABLE assemblies(
 assembly_id INTEGER PRIMARY KEY,source_id INTEGER NOT NULL REFERENCES sources(source_id),
 source_native_pk INTEGER NOT NULL,accession TEXT NOT NULL UNIQUE,accession_stem TEXT NOT NULL,
 status TEXT NOT NULL,validation_status TEXT,source_document_sha256 TEXT NOT NULL,
 tax_id INTEGER,organism_name TEXT,source_sequence_count INTEGER);
CREATE TABLE contigs(
 contig_id INTEGER PRIMARY KEY,assembly_id INTEGER NOT NULL REFERENCES assemblies(assembly_id),
 source_native_pk INTEGER NOT NULL,sequence_id TEXT NOT NULL,sequence_sha256 TEXT,length INTEGER);
CREATE TABLE nodes(
 node_id INTEGER PRIMARY KEY,sequence TEXT NOT NULL UNIQUE,sequence_sha256 TEXT NOT NULL UNIQUE,
 length INTEGER NOT NULL,eligible INTEGER NOT NULL,
 bona_fide_consensus_count INTEGER NOT NULL DEFAULT 0,possible_consensus_count INTEGER NOT NULL DEFAULT 0,
 bona_fide_repeat_count INTEGER NOT NULL DEFAULT 0,possible_repeat_count INTEGER NOT NULL DEFAULT 0);
CREATE TABLE arrays(
 array_id INTEGER PRIMARY KEY,contig_id INTEGER NOT NULL REFERENCES contigs(contig_id),
 source_native_pk INTEGER NOT NULL,call_id TEXT NOT NULL,tool TEXT NOT NULL,tool_version TEXT,
 start INTEGER,end INTEGER,strand TEXT,consensus_node_id INTEGER REFERENCES nodes(node_id),
 consensus_orientation INTEGER,consensus_raw TEXT,repeat_count INTEGER,spacer_count INTEGER,
 category TEXT NOT NULL,evidence_level TEXT,certainty_score REAL);
CREATE TABLE repeat_occurrences(
 repeat_id INTEGER PRIMARY KEY,array_id INTEGER NOT NULL REFERENCES arrays(array_id),
 source_native_pk INTEGER NOT NULL,ordinal INTEGER NOT NULL,start INTEGER,end INTEGER,
 node_id INTEGER REFERENCES nodes(node_id),orientation INTEGER,raw_sequence TEXT);
"""
INDEX_SQL = """
CREATE INDEX IF NOT EXISTS idx_assemblies_source ON assemblies(source_id);
CREATE INDEX IF NOT EXISTS idx_assemblies_stem ON assemblies(accession_stem);
CREATE INDEX IF NOT EXISTS idx_assemblies_taxon ON assemblies(tax_id);
CREATE INDEX IF NOT EXISTS idx_contigs_assembly ON contigs(assembly_id);
CREATE INDEX IF NOT EXISTS idx_contigs_sha ON contigs(sequence_sha256);
CREATE INDEX IF NOT EXISTS idx_arrays_contig ON arrays(contig_id);
CREATE INDEX IF NOT EXISTS idx_arrays_consensus ON arrays(consensus_node_id);
CREATE INDEX IF NOT EXISTS idx_repeats_array ON repeat_occurrences(array_id,ordinal);
CREATE INDEX IF NOT EXISTS idx_repeats_node ON repeat_occurrences(node_id);
CREATE INDEX IF NOT EXISTS idx_nodes_length ON nodes(length) WHERE eligible=1;
"""


def _directory(path):
    path = Path(path).absolute()
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise InputError("Catalogue directory must not traverse symlinks")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _connect(path, readonly=False):
    if readonly:
        checked_path(path)
        connection = sqlite3.connect(Path(path).absolute().as_uri() + "?mode=ro&immutable=1", uri=True)
        connection.execute("PRAGMA query_only=ON")
    else:
        if Path(path).is_symlink():
            raise InputError("Catalogue file must not be a symlink")
        # ATTACH must honor mode=ro&immutable=1 even on SQLite builds without
        # USE_URI. Enable URI parsing on the main connection explicitly.
        connection = sqlite3.connect(Path(path).absolute().as_uri(), uri=True, timeout=30)
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA trusted_schema=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-32768")
    return connection


def _text(value, label, *, nullable=False, maximum=1000):
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not 1 <= len(value) <= maximum or any(ord(c) < 32 for c in value):
        raise InputError(f"Invalid source {label}")
    return value


def _integer(value, label, *, nullable=False, low=0):
    if value is None and nullable:
        return None
    if type(value) is not int or not low <= value < 2**63:
        raise InputError(f"Invalid source {label}")
    return value


def _key(source_id, native_pk):
    if type(native_pk) is not int or not 0 < native_pk < STRIDE:
        raise InputError("Source native primary key is outside the supported range")
    return source_id * STRIDE + native_pk


def normalized_sequence(value):
    if value is None:
        return None, None, None
    _text(value, "nucleotide sequence", maximum=10000)
    normalized = value.upper().replace("U", "T")
    if not set(normalized) <= set("ACGTRYSWKMBDHVN"):
        raise InputError("Source repeat contains unsupported non-nucleotide symbols")
    reverse = reverse_complement(normalized)
    canonical = min(normalized, reverse)
    orientation = 0 if normalized == reverse else (1 if normalized == canonical else -1)
    return canonical, orientation, value if value != normalized else None


def reconstruct_sequence(canonical, orientation, raw=None):
    if raw is not None:
        return raw
    if canonical is None:
        return None
    if orientation not in {-1, 0, 1}:
        raise InputError("Invalid stored sequence orientation")
    return reverse_complement(canonical) if orientation == -1 else canonical


class CatalogWriter:
    def __init__(self, directory, expected_profile, input_digest, build_kind, min_free_bytes):
        validate_free_reserve(min_free_bytes)
        self.directory = _directory(directory)
        self.profile = required_digest(expected_profile)
        self.minimum = min_free_bytes
        self.path = self.directory / "catalog.building.sqlite3"
        published = self.directory / "catalog.sqlite3"
        if (self.directory / "manifest.json").exists():
            raise InputError("Catalogue is already published")
        if published.exists():
            # A crash after the database rename but before the final marker may
            # resume only when the unfinished publication belongs to this build.
            if self.path.exists():
                raise InputError("Both build and unpublished final databases exist")
            prior = _connect(published, readonly=True)
            try:
                prior_metadata = dict(prior.execute("SELECT key,value FROM metadata"))
            finally:
                prior.close()
            if (prior_metadata.get("input_manifest_sha256") != input_digest
                    or prior_metadata.get("profile_sha256") != self.profile
                    or prior_metadata.get("build_kind") != build_kind):
                raise InputError("Unpublished final database belongs to a different build")
            os.replace(published, self.path)
        check_free(self.directory, self.minimum)
        existed = self.path.exists()
        self.db = _connect(self.path)
        self.cache = OrderedDict()
        if not existed:
            self.db.executescript(SCHEMA_SQL)
            self.db.executemany("INSERT INTO metadata VALUES (?,?)", [
                ("schema", json.dumps(SCHEMA, sort_keys=True)), ("profile_sha256", self.profile),
                ("input_manifest_sha256", input_digest), ("build_kind", build_kind),
                ("normalization", NORMALIZATION),
            ])
            self.db.commit()
        metadata = dict(self.db.execute("SELECT key,value FROM metadata"))
        if (metadata.get("schema") != json.dumps(SCHEMA, sort_keys=True)
                or metadata.get("profile_sha256") != self.profile
                or metadata.get("input_manifest_sha256") != input_digest
                or metadata.get("build_kind") != build_kind):
            self.db.close()
            raise InputError("Resume requires the same plan, profile and catalogue schema")

    def close(self):
        self.db.close()

    def _node(self, sequence, category, kind):
        canonical, orientation, raw = normalized_sequence(sequence)
        if canonical is None:
            return None, orientation, raw
        field = ("bona_fide_" if category == "Bona-fide" else "possible_") + kind + "_count"
        ident = self.cache.get(canonical)
        if ident is None:
            row = self.db.execute("SELECT node_id FROM nodes WHERE sequence=?", (canonical,)).fetchone()
            if row is None:
                cursor = self.db.execute("INSERT INTO nodes(sequence,sequence_sha256,length,eligible) VALUES (?,?,?,?)",
                                         (canonical, hashlib.sha256(canonical.encode()).hexdigest(), len(canonical),
                                          int(set(canonical) <= set("ACGT") and len(canonical) <= 200)))
                ident = cursor.lastrowid
            else:
                ident = row[0]
            self.cache[canonical] = ident
            if len(self.cache) > 100_000:
                self.cache.popitem(last=False)
        else:
            self.cache.move_to_end(canonical)
        self.db.execute(f"UPDATE nodes SET {field}={field}+1 WHERE node_id=?", (ident,))
        return ident, orientation, raw

    def ingest(self, source):
        identity = source_identity(source)
        if self.db.execute("SELECT 1 FROM sources WHERE source_key=?", (identity["source_key"],)).fetchone():
            return False
        check_free(self.directory, self.minimum)
        with open_source(source, self.profile, self.directory, self.minimum) as opened:
            metadata = opened.metadata
            counts = {"assemblies": 0, "contigs": 0, "arrays": 0, "repeat_occurrences": 0,
                      "complete_assemblies": 0, "failed_assemblies": 0, "unavailable_assemblies": 0,
                      "bona_fide_arrays": 0, "possible_arrays": 0,
                      "missing_consensus_sequences": 0, "missing_repeat_sequences": 0}
            self.db.execute("BEGIN IMMEDIATE")
            try:
                source_id = self.db.execute("INSERT INTO sources(source_key,kind,source_sha256,database_sha256,collection_index_sha256,profile_sha256,profile_file_sha256,selection_sha256,counts_json) VALUES (?,?,?,?,?,?,?,?,?)",
                    tuple(metadata.get(key) for key in ("source_key", "kind", "source_sha256", "database_sha256", "collection_index_sha256", "profile_sha256", "profile_file_sha256", "selection_sha256")) + ("{}",)).lastrowid
                for row in opened.assemblies():
                    accession = row["accession"]
                    if not isinstance(accession, str) or not ACCESSION.fullmatch(accession):
                        raise InputError("Invalid source assembly accession")
                    status = row["status"]
                    if status not in {"complete", "failed", "unavailable"}:
                        raise InputError("Source assembly outcome is nonterminal")
                    tax_id, organism = row["taxonomy"]
                    self.db.execute("INSERT INTO assemblies VALUES (?,?,?,?,?,?,?,?,?,?,?)", (
                        _key(source_id, row["native_pk"]), source_id, row["native_pk"], accession, accession.rsplit(".", 1)[0],
                        status, row["validation_status"], required_digest(row["source_document_sha256"]), tax_id, organism,
                        _integer(row["sequence_count"], "source sequence count", nullable=True)))
                    counts["assemblies"] += 1
                    counts[status + "_assemblies"] += 1
                for row in opened.contigs():
                    self.db.execute("INSERT INTO contigs VALUES (?,?,?,?,?,?)", (
                        _key(source_id, row["native_pk"]), _key(source_id, row["assembly_pk"]), row["native_pk"],
                        _text(row["sequence_id"], "contig identifier"),
                        required_digest(row["sequence_sha256"]) if row["sequence_sha256"] is not None else None,
                        _integer(row["length"], "contig length", nullable=True, low=1)))
                    counts["contigs"] += 1
                arrays = {}
                for row in opened.arrays():
                    category = row["category"]
                    if category not in {"Bona-fide", "Possible"}:
                        raise InputError("Source array tier differs from the accepted-only campaign")
                    if row["strand"] is not None:
                        raise InputError("Source profile did not compute biological strand; expected null")
                    array_id = _key(source_id, row["native_pk"])
                    if array_id in arrays or len(arrays) >= 1_000_000:
                        raise InputError("Duplicate source array or per-source array budget exceeded")
                    repeat_count = _integer(row["repeat_count"], "repeat count", nullable=True)
                    arrays[array_id] = [category, repeat_count, 0]
                    node_id, orientation, raw = self._node(row["repeat_consensus"], category, "consensus")
                    start, end = _integer(row["start"], "array start", nullable=True, low=1), _integer(row["end"], "array end", nullable=True, low=1)
                    if start is not None and end is not None and start > end:
                        raise InputError("Source array interval is reversed")
                    score = row["certainty_score"]
                    if score is not None and (type(score) not in (int, float) or not math.isfinite(score)):
                        raise InputError("Source certainty score is invalid")
                    self.db.execute("INSERT INTO arrays VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                        array_id, _key(source_id, row["sequence_pk"]), row["native_pk"], _text(row["call_id"], "array call identifier"),
                        _text(row["tool"], "tool"), _text(row["tool_version"], "tool version", nullable=True), start, end, None,
                        node_id, orientation, raw, repeat_count, _integer(row["spacer_count"], "spacer count", nullable=True),
                        category, _text(str(row["evidence_level"]), "evidence level") if row["evidence_level"] is not None else None, score))
                    counts["arrays"] += 1
                    counts["bona_fide_arrays" if category == "Bona-fide" else "possible_arrays"] += 1
                    counts["missing_consensus_sequences"] += node_id is None
                # Source-local primary-key range avoids rescanning the catalogue
                # on each import. Failed receipts remain valid outcomes, but
                # their partial features must never contribute to core evidence.
                lower = source_id * STRIDE
                invalid = self.db.execute(
                    "SELECT 1 FROM arrays a CROSS JOIN contigs c ON c.contig_id=a.contig_id "
                    "CROSS JOIN assemblies s ON s.assembly_id=c.assembly_id "
                    "WHERE a.array_id>? AND a.array_id<? "
                    "AND (s.status<>'complete' OR s.validation_status IS NOT 'valid') LIMIT 1",
                    (lower, lower + STRIDE)).fetchone()
                if invalid is not None:
                    raise InputError("Source arrays require a complete, validated owning assembly")
                for row in opened.repeats():
                    array_id = _key(source_id, row["array_pk"])
                    if array_id not in arrays:
                        raise InputError("Source repeat parent is absent from selected arrays")
                    info = arrays[array_id]
                    ordinal = _integer(row["ordinal"], "repeat ordinal", low=1)
                    if ordinal != info[2] + 1:
                        raise InputError("Source repeat ordinals are not contiguous within their array")
                    info[2] += 1
                    node_id, orientation, raw = self._node(row["sequence"], info[0], "repeat")
                    start, end = _integer(row["start"], "repeat start", nullable=True, low=1), _integer(row["end"], "repeat end", nullable=True, low=1)
                    if start is not None and end is not None:
                        if start > end or (row["sequence"] is not None and len(row["sequence"]) != end - start + 1):
                            raise InputError("Source repeat length and interval disagree")
                    self.db.execute("INSERT INTO repeat_occurrences VALUES (?,?,?,?,?,?,?,?,?)", (
                        _key(source_id, row["native_pk"]), array_id, row["native_pk"], ordinal, start, end, node_id, orientation, raw))
                    counts["repeat_occurrences"] += 1
                    counts["missing_repeat_sequences"] += node_id is None
                    if counts["repeat_occurrences"] % 100_000 == 0:
                        check_free(self.directory, self.minimum)
                if opened.expected_arrays is not None and counts["arrays"] != opened.expected_arrays:
                    raise InputError("Selected source array count does not reconcile")
                if any(expected is not None and expected != observed for _, expected, observed in arrays.values()):
                    raise InputError("Source repeat occurrences do not reconcile with declared array counts")
                self.db.execute("UPDATE sources SET counts_json=? WHERE source_id=?", (json.dumps(counts, sort_keys=True), source_id))
                self.db.commit()
            except BaseException:
                self.db.rollback()
                self.cache.clear()
                raise
        self.db.execute("PRAGMA wal_checkpoint(PASSIVE)")
        return True

    def execute_guarded(self, statement, parameters=()):
        """Monitor reserve during large SQLite statements; never a disk quota."""
        failure = None
        last_check = 0.0
        directories = (self.directory, Path(os.environ.get("SQLITE_TMPDIR") or tempfile.gettempdir()))

        def progress():
            nonlocal failure, last_check
            now = time.monotonic()
            if now - last_check < 0.25:
                return 0
            try:
                for directory in directories:
                    check_free(directory, self.minimum)
            except (InputError, OSError) as error:
                failure = error
                return 1
            last_check = now
            return 0

        if progress():
            raise failure
        self.db.set_progress_handler(progress, 10000)
        try:
            return self.db.execute(statement, parameters)
        except sqlite3.Error as error:
            if failure is not None:
                raise failure from error
            raise
        finally:
            self.db.set_progress_handler(None, 0)

    def publish(self, extra=None):
        check_free(self.directory, self.minimum)
        # Each complete index is resumable; truncate its WAL before the next
        # sort/build rather than retaining a multi-index publication WAL.
        for statement in INDEX_SQL.split(";"):
            if not statement.strip():
                continue
            self.execute_guarded(statement)
            self.db.commit()
            self.execute_guarded("PRAGMA wal_checkpoint(TRUNCATE)").fetchone()
            check_free(self.directory, self.minimum)
        if self.execute_guarded("PRAGMA foreign_key_check").fetchone() is not None:
            raise InputError("Catalogue contains broken occurrence relations")
        if self.execute_guarded("PRAGMA quick_check").fetchone()[0] != "ok":
            raise InputError("Catalogue SQLite integrity check failed")
        counts = {table: self.execute_guarded(f"SELECT count(*) FROM {table}").fetchone()[0]
                  for table in ("sources", "assemblies", "contigs", "arrays", "repeat_occurrences", "nodes")}
        counts["eligible_nodes"] = self.execute_guarded("SELECT count(*) FROM nodes WHERE eligible=1").fetchone()[0]
        counts.update({status + "_assemblies": count for status, count in self.execute_guarded("SELECT status,count(*) FROM assemblies GROUP BY status")})
        expected = {}
        selection = hashlib.sha256()
        for source_key, raw in self.execute_guarded("SELECT source_key,counts_json FROM sources ORDER BY source_key"):
            selection.update((source_key + "\n").encode())
            for key, count in json.loads(raw).items():
                expected[key] = expected.get(key, 0) + count
        for key, count in expected.items():
            if key in counts and counts[key] != count:
                raise InputError("Catalogue row totals differ from committed source receipts")
            counts[key] = count
        node_totals = self.execute_guarded("SELECT " + ",".join(f"coalesce(sum({key}),0)" for key in NODE_COUNTS) + " FROM nodes").fetchone()
        if sum(node_totals[:2]) != counts["arrays"] - counts.get("missing_consensus_sequences", 0) or sum(node_totals[2:]) != counts["repeat_occurrences"] - counts.get("missing_repeat_sequences", 0):
            raise InputError("Catalogue node evidence counts do not reconcile with occurrences")
        metadata = dict(self.execute_guarded("SELECT key,value FROM metadata"))
        self.execute_guarded("PRAGMA wal_checkpoint(TRUNCATE)")
        self.execute_guarded("PRAGMA journal_mode=DELETE")
        self.db.close()
        published = self.directory / "catalog.sqlite3"
        with self.path.open("rb") as handle:
            os.fsync(handle.fileno())
        check_free(self.directory, self.minimum)
        os.replace(self.path, published)
        manifest = {"schema": dict(SCHEMA), "status": "complete", "profile_sha256": self.profile,
                    "build_kind": metadata["build_kind"], "input_manifest_sha256": metadata["input_manifest_sha256"],
                    "normalization": NORMALIZATION, "coordinate_system": "1-based-inclusive-forward-source",
                    "biological_orientation": "not_computed", "selection_sha256": selection.hexdigest(),
                    "source_native_pk_scope": "original SQLite primary key for bundle sources; deterministic per-source traversal ordinal for raw JSON sources",
                    "occurrence_orientation_scope": "relative to canonical sequence; biological strand remains null",
                    "database": {"file": "catalog.sqlite3", "sha256": sha256_file(published), "bytes": published.stat().st_size},
                    "counts": counts, "node_evidence_counts": dict(zip(NODE_COUNTS, node_totals)),
                    "generated_utc": datetime.now(timezone.utc).isoformat(),
                    "interpretation": "Predicted source repeat observations; no family, structural motif, Cas subtype or host assignment."}
        if extra:
            manifest.update(extra)
        atomic_json(self.directory / "manifest.json", manifest)
        return manifest


def load_catalog(directory, verify_hash=True):
    directory = Path(directory).absolute()
    manifest = bounded_json(directory / "manifest.json", 8 * 1024**2)
    record = manifest.get("database", {})
    if (manifest.get("schema") != SCHEMA or manifest.get("status") != "complete"
            or record.get("file") != "catalog.sqlite3" or manifest.get("normalization") != NORMALIZATION):
        raise InputError("Unsupported or incomplete catalogue manifest")
    required_digest(record.get("sha256"))
    required_digest(manifest.get("profile_sha256"))
    path = checked_path(directory / "catalog.sqlite3")
    if path.stat().st_size != record.get("bytes") or (verify_hash and sha256_file(path) != record["sha256"]):
        raise InputError("Catalogue file differs from its published manifest")
    connection = _connect(path, readonly=True)
    metadata = dict(connection.execute("SELECT key,value FROM metadata"))
    if metadata.get("profile_sha256") != manifest["profile_sha256"] or metadata.get("schema") != json.dumps(SCHEMA, sort_keys=True):
        connection.close()
        raise InputError("Catalogue metadata differs from its manifest")
    return connection, manifest


def build_partition(plan_path, output_dir, expected_profile, workers=1, *, min_free_bytes=2 * 1024**3):
    validate_free_reserve(min_free_bytes)
    if workers != 1:
        raise InputError("Each partition has one writer; parallelize independent partition plans")
    plan = bounded_json(plan_path, 256 * 1024**2)
    if plan.get("profile_sha256", plan.get("analysis_profile_sha256")) != expected_profile:
        raise InputError("Source plan differs from the requested scientific profile")
    sources = list(plan_sources(plan))
    keys = [source_identity(source)["source_key"] for source in sources]
    if len(keys) != len(set(keys)):
        raise InputError("Partition plan contains duplicate source selections")
    # Bind only scientific identities, not private paths, hostnames or gzip metadata.
    input_digest = canonical_digest({"profile_sha256": expected_profile, "source_keys": keys})
    output = Path(output_dir)
    if (output / "manifest.json").exists():
        connection, manifest = load_catalog(output)
        connection.close()
        if manifest.get("input_manifest_sha256") != input_digest:
            raise InputError("Published partition differs from the requested source selection")
        return manifest
    writer = CatalogWriter(output, expected_profile, input_digest, "partition", min_free_bytes)
    try:
        committed = {row[0] for row in writer.db.execute("SELECT source_key FROM sources")}
        if not committed <= set(keys):
            raise InputError("Build file contains sources outside its requested selection")
        for index, source in enumerate(sources, 1):
            changed = writer.ingest(source)
            if changed and (index % 25 == 0 or index == len(sources)):
                print(json.dumps({"event": "catalog_source_progress", "completed": index, "total": len(sources)}, sort_keys=True), flush=True)
        if writer.db.execute("SELECT count(*) FROM sources").fetchone()[0] != len(keys):
            raise InputError("Partition does not contain every requested source")
        return writer.publish()
    except BaseException:
        try:
            writer.close()
        except sqlite3.Error:
            pass
        raise


def merge_catalog(partition_dirs, destination, *, min_free_bytes=2 * 1024**3):
    """Merge immutable partitions using disk-backed ID maps and SQL table copies."""
    validate_free_reserve(min_free_bytes)
    # Only one small manifest per partition is retained; no occurrence lists.
    parts = []
    profile = None
    for directory in partition_dirs:
        connection, manifest = load_catalog(directory)
        connection.close()
        if profile is not None and profile != manifest["profile_sha256"]:
            raise InputError("Cannot merge different scientific profiles")
        profile = manifest["profile_sha256"]
        parts.append((Path(directory).absolute(), manifest))
    if not parts:
        raise InputError("At least one partition is required")
    parts.sort(key=lambda item: item[1]["database"]["sha256"])
    hashes = [manifest["database"]["sha256"] for _, manifest in parts]
    if len(hashes) != len(set(hashes)):
        raise InputError("Duplicate partition database")
    input_digest = canonical_digest({"profile_sha256": profile, "partition_database_sha256": hashes})
    destination = Path(destination)
    if (destination / "manifest.json").exists():
        connection, manifest = load_catalog(destination)
        connection.close()
        if manifest.get("input_manifest_sha256") != input_digest:
            raise InputError("Published catalogue differs from the requested partitions")
        return manifest
    writer = CatalogWriter(destination, profile, input_digest, "merged", min_free_bytes)
    db = writer.db
    execute = writer.execute_guarded
    try:
        execute("CREATE TABLE IF NOT EXISTS merged_parts(database_sha256 TEXT PRIMARY KEY)")
        db.commit()
        for directory, manifest in parts:
            part_sha = manifest["database"]["sha256"]
            if execute("SELECT 1 FROM merged_parts WHERE database_sha256=?", (part_sha,)).fetchone():
                continue
            check_free(writer.directory, min_free_bytes, 2 * manifest["database"]["bytes"])
            uri = (directory / "catalog.sqlite3").as_uri() + "?mode=ro&immutable=1"
            execute("ATTACH DATABASE ? AS incoming", (uri,))
            try:
                execute("BEGIN IMMEDIATE")
                if execute("SELECT 1 FROM incoming.sources i JOIN sources s USING(source_key) LIMIT 1").fetchone():
                    raise InputError("Partition selections overlap in source identity")
                if execute("SELECT 1 FROM incoming.assemblies i JOIN assemblies a USING(accession) LIMIT 1").fetchone():
                    raise InputError("Partition selections overlap in assembly accession")
                source_offset = execute("SELECT coalesce(max(source_id),0) FROM sources").fetchone()[0]
                offset = source_offset * STRIDE
                execute("INSERT INTO sources SELECT source_id+?,source_key,kind,source_sha256,database_sha256,collection_index_sha256,profile_sha256,profile_file_sha256,selection_sha256,counts_json FROM incoming.sources", (source_offset,))
                execute("INSERT INTO assemblies SELECT assembly_id+?,source_id+?,source_native_pk,accession,accession_stem,status,validation_status,source_document_sha256,tax_id,organism_name,source_sequence_count FROM incoming.assemblies", (offset, source_offset))
                execute("INSERT INTO contigs SELECT contig_id+?,assembly_id+?,source_native_pk,sequence_id,sequence_sha256,length FROM incoming.contigs", (offset, offset))
                fields = ",".join(NODE_COUNTS)
                update = ",".join(f"{field}={field}+excluded.{field}" for field in NODE_COUNTS)
                execute(f"INSERT INTO nodes(sequence,sequence_sha256,length,eligible,{fields}) SELECT sequence,sequence_sha256,length,eligible,{fields} FROM incoming.nodes WHERE 1 ON CONFLICT(sequence) DO UPDATE SET {update}")
                execute("CREATE TEMP TABLE node_map(old_id INTEGER PRIMARY KEY,new_id INTEGER NOT NULL)")
                execute("INSERT INTO node_map SELECT i.node_id,n.node_id FROM incoming.nodes i JOIN nodes n ON n.sequence=i.sequence")
                execute("INSERT INTO arrays SELECT a.array_id+?,a.contig_id+?,a.source_native_pk,a.call_id,a.tool,a.tool_version,a.start,a.end,a.strand,m.new_id,a.consensus_orientation,a.consensus_raw,a.repeat_count,a.spacer_count,a.category,a.evidence_level,a.certainty_score FROM incoming.arrays a LEFT JOIN node_map m ON m.old_id=a.consensus_node_id", (offset, offset))
                execute("INSERT INTO repeat_occurrences SELECT r.repeat_id+?,r.array_id+?,r.source_native_pk,r.ordinal,r.start,r.end,m.new_id,r.orientation,r.raw_sequence FROM incoming.repeat_occurrences r LEFT JOIN node_map m ON m.old_id=r.node_id", (offset, offset))
                execute("DROP TABLE node_map")
                execute("INSERT INTO merged_parts VALUES (?)", (part_sha,))
                db.commit()
            except BaseException:
                db.rollback()
                raise
            finally:
                # Cleanup must remain possible after the reserve watchdog fires.
                db.execute("DETACH DATABASE incoming")
            execute("PRAGMA wal_checkpoint(TRUNCATE)")
            print(json.dumps({"event": "catalog_partition_merged", "completed": execute("SELECT count(*) FROM merged_parts").fetchone()[0], "total": len(parts)}, sort_keys=True), flush=True)
        return writer.publish({"partition_database_sha256": hashes})
    except BaseException:
        try:
            writer.close()
        except sqlite3.Error:
            pass
        raise
