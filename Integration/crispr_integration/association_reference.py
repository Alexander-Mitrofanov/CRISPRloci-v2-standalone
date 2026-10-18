"""Immutable broad-reference verification and bounded target sequence reads."""

from pathlib import Path
from collections import OrderedDict
import hashlib
import re
import sqlite3
import stat

from .protospacer import _payload, _json as _parse_json, _records


def _read_manifest(root):
    path = root / "reference.json"
    if path.stat().st_size > 2 * 1024 * 1024:
        raise ValueError("Association manifest exceeds budget")
    return _payload(root, path)


def _read_genome(path):
    return _records(path)


SUPPORTED_VERSIONS = {"2.1.0", "2.1.1"}
SHA = re.compile(r"[0-9a-f]{64}\Z")
ID = re.compile(r"[A-Za-z0-9_.:-]{1,160}\Z")
GENOME = re.compile(r"genomes/[A-Za-z0-9_.-]+\Z")
INDEX = re.compile(r"db/[A-Za-z0-9_.-]+\Z")
MAX_GENOMES = 100_000
MAX_BASES = 10_000_000_000
MAX_GENOME = 5_000_000
MAX_CATALOGUE = 128 * 1024 * 1024
MAX_INDEX = 16 * 1024 * 1024 * 1024
MAX_INDEX_TOTAL = 64 * 1024 * 1024 * 1024
PUBLIC_FIELDS = ("name", "engine_version", "control_method", "genome_count", "total_bases", "scope", "source")


def _check_path(path):
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise ValueError("Unsafe association reference path")


def _identity(path):
    info = path.stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def _verified_file(path, identity, limit, checkpoint):
    if (
        not isinstance(identity, dict)
        or set(identity) != {"sha256", "size_bytes"}
        or not isinstance(identity["sha256"], str)
        or not SHA.fullmatch(identity["sha256"])
        or type(identity["size_bytes"]) is not int
        or not 0 <= identity["size_bytes"] <= limit
    ):
        raise ValueError("Invalid association reference file identity")
    _check_path(path)
    info = path.stat()
    if not stat.S_ISREG(info.st_mode) or info.st_size != identity["size_bytes"]:
        raise ValueError("Association reference file size/type changed")
    before = _identity(path)
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            if checkpoint:
                checkpoint()
            size += len(chunk)
            if size > identity["size_bytes"]:
                raise ValueError("Association reference file grew")
            digest.update(chunk)
    if size != identity["size_bytes"] or digest.hexdigest() != identity["sha256"] or before != _identity(path):
        raise ValueError("Association reference checksum or identity changed")


class TargetSequences:
    """Read-only SQLite identity lookup with at most eight loaded genome sequences."""

    def __init__(self, root, count, total, checkpoint=None):
        self.root, self.checkpoint = root, checkpoint
        self.cache = OrderedDict()
        self.connection = None
        self.connection = sqlite3.connect((root / "catalogue.sqlite3").as_uri() + "?mode=ro&immutable=1", uri=True)
        try:
            self.connection.execute("PRAGMA query_only = ON")
            self.connection.execute("PRAGMA trusted_schema = OFF")
            self.connection.execute("PRAGMA cache_size = -4096")
            schema = self.connection.execute(
                "SELECT type, name FROM sqlite_master WHERE name NOT LIKE 'sqlite_%'"
            ).fetchall()
            columns = self.connection.execute("PRAGMA table_info(genomes)").fetchall()
            if schema != [("table", "genomes")] or [(row[1], row[2], row[5]) for row in columns] != [
                ("id", "TEXT", 1),
                ("length", "INTEGER", 0),
                ("sha256", "TEXT", 0),
                ("file", "TEXT", 0),
            ]:
                raise ValueError("Invalid association catalogue schema")
            found, bases = 0, 0
            for ident, length, digest, filename in self.connection.execute("SELECT id,length,sha256,file FROM genomes"):
                if checkpoint:
                    checkpoint()
                found += 1
                if (
                    found > count
                    or not isinstance(ident, str)
                    or not ID.fullmatch(ident)
                    or type(length) is not int
                    or not 1 <= length <= MAX_GENOME
                    or not isinstance(digest, str)
                    or not SHA.fullmatch(digest)
                    or not isinstance(filename, str)
                    or not GENOME.fullmatch(filename)
                    or filename.split("/")[1] in {".", ".."}
                ):
                    raise ValueError("Invalid association target metadata")
                genome_path = root / filename
                _check_path(genome_path)
                if not genome_path.is_file() or genome_path.stat().st_size > 6_000_000:
                    raise ValueError("Missing or oversized association target file")
                bases += length
                if bases > total:
                    raise ValueError("Association catalogue base count exceeded")
            if found != count or bases != total:
                raise ValueError("Association catalogue counts changed")
            if self.connection.execute("SELECT COUNT(DISTINCT file) FROM genomes").fetchone()[0] != count:
                raise ValueError("Duplicate association target source files")
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        self.cache.clear()

    def __del__(self):
        self.close()

    def __contains__(self, ident):
        return (
            isinstance(ident, str)
            and self.connection.execute("SELECT 1 FROM genomes WHERE id=?", (ident,)).fetchone() is not None
        )

    def __getitem__(self, ident):
        if ident in self.cache:
            self.cache.move_to_end(ident)
            return self.cache[ident]
        if self.checkpoint:
            self.checkpoint()
        row = self.connection.execute("SELECT length,sha256,file FROM genomes WHERE id=?", (ident,)).fetchone()
        if row is None:
            raise KeyError(ident)
        length, digest, filename = row
        path = self.root / filename
        _check_path(path)
        if not path.is_file() or path.stat().st_size > 6_000_000:
            raise ValueError("Unsafe association target file")
        data = _read_genome(path)
        if set(data) != {ident}:
            raise ValueError("Association target identifier changed")
        sequence = data[ident]
        if (
            len(sequence) != length
            or set(sequence) - set("ACGTN")
            or hashlib.sha256(sequence.encode("ascii")).hexdigest() != digest
        ):
            raise ValueError("Association target source checksum changed")
        self.cache[ident] = sequence
        if len(self.cache) > 8:
            self.cache.popitem(last=False)
        return sequence


def reference_data(root, expected, checkpoint=None):
    root = Path(root).absolute()
    _check_path(root)
    raw = _read_manifest(root)
    if len(raw) > 2 * 1024 * 1024 or hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("Association reference digest changed")
    doc = _parse_json(raw)
    fields = {"schema_version", *PUBLIC_FIELDS, "catalogue", "files"}
    if (
        not isinstance(doc, dict)
        or set(doc) != fields
        or doc["schema_version"] != "spacepharer.reference/2.0.0"
        or doc["engine_version"] not in SUPPORTED_VERSIONS
        or doc["control_method"] != "reversed_translated_orfs"
        or not isinstance(doc["name"], str)
        or not ID.fullmatch(doc["name"])
        or type(doc["genome_count"]) is not int
        or not 1 <= doc["genome_count"] <= MAX_GENOMES
        or type(doc["total_bases"]) is not int
        or not doc["genome_count"] <= doc["total_bases"] <= MAX_BASES
        or not isinstance(doc["scope"], str)
        or doc["scope"] not in {"broad_phage_reference", "targeted_panel"}
        or not isinstance(doc["source"], dict)
        or set(doc["source"]) != {"collection", "release", "retrieved_at", "selection_policy"}
        or any(
            not isinstance(value, str) or not 1 <= len(value) <= 1000 or any(ord(c) < 32 for c in value)
            for value in doc["source"].values()
        )
        or not isinstance(doc["files"], dict)
        or not 1 <= len(doc["files"]) <= 3000
    ):
        raise ValueError("Invalid immutable association reference")
    db = root / "db"
    _check_path(db)
    if not db.is_dir():
        raise ValueError("Missing association index directory")
    inventory = set()
    for path in db.iterdir():
        inventory.add("db/" + path.name)
        if len(inventory) > 3000:
            raise ValueError("Association index inventory exceeded")
    required = {f"db/{prefix}{suffix}" for prefix in ("target", "control") for suffix in ("", ".index", ".dbtype")}
    if inventory != set(doc["files"]) or not required <= inventory:
        raise ValueError("Incomplete association index inventory")
    total_bytes = 0
    for name, identity in doc["files"].items():
        if not isinstance(name, str) or not INDEX.fullmatch(name) or name.split("/")[1] in {".", ".."}:
            raise ValueError("Invalid association index name")
        _verified_file(root / name, identity, MAX_INDEX, checkpoint)
        total_bytes += identity["size_bytes"]
        if total_bytes > MAX_INDEX_TOTAL:
            raise ValueError("Association index bytes exceed reference budget")
    catalogue = doc["catalogue"]
    if (
        not isinstance(catalogue, dict)
        or set(catalogue) != {"file", "sha256", "size_bytes"}
        or catalogue["file"] != "catalogue.sqlite3"
    ):
        raise ValueError("Invalid association catalogue identity")
    _verified_file(
        root / catalogue["file"], {key: catalogue[key] for key in ("sha256", "size_bytes")}, MAX_CATALOGUE, checkpoint
    )
    try:
        targets = TargetSequences(root, doc["genome_count"], doc["total_bases"], checkpoint)
    except sqlite3.Error as error:
        raise ValueError("Invalid association catalogue") from error
    return {key: doc[key] for key in PUBLIC_FIELDS} | {"manifest_sha256": expected}, targets
