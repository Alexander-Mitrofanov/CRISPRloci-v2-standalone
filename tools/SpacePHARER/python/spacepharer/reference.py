"""Operator-built, immutable references with bounded lazy sequence access."""

from __future__ import annotations

import io
import re
import shutil
import sqlite3
import tarfile
from collections import OrderedDict
from collections.abc import Iterator, Mapping
from pathlib import Path

from . import __version__
from .engine import run
from .io import (
    ID,
    InputError,
    digest,
    fasta,
    file_digest,
    json_value,
    read,
    reject_symlinks,
    transaction,
    write_json,
)

SCHEMA = "spacepharer.reference/2.0.0"
# Patch validation changes do not alter the target/control database format.
REFERENCE_ENGINE_VERSIONS = {"2.1.0", "2.1.1"}
MAX_GENOMES = 100_000
MAX_BASES = 10_000_000_000
MAX_GENOME_BASES = 5_000_000
MAX_CATALOGUE_BYTES = 128 * 1024 * 1024
MAX_INDEX_BYTES = 16 * 1024**3
MAX_INDEX_TOTAL_BYTES = 64 * 1024**3
SHA = re.compile(r"[0-9a-f]{64}\Z")
SOURCE_FIELDS = {"collection", "release", "retrieved_at", "selection_policy"}
PUBLIC_FIELDS = ("name", "engine_version", "control_method", "genome_count", "total_bases", "scope", "source")


def source_metadata(value: dict | None) -> dict:
    if value is None:
        value = {
            "collection": "Operator-provided FASTAs",
            "release": "unspecified",
            "retrieved_at": "unspecified",
            "selection_policy": "All explicitly supplied single-genome DNA FASTAs; no automatic selection.",
        }
    if (
        not isinstance(value, dict)
        or set(value) != SOURCE_FIELDS
        or any(
            not isinstance(text, str) or not 1 <= len(text) <= 1000 or any(ord(c) < 32 for c in text)
            for text in value.values()
        )
    ):
        raise InputError("Reference source metadata requires four bounded plain-text fields")
    return dict(value)


def genome_list(path: Path) -> list[Path]:
    """A bounded newline list avoids OS argument-length limits for operator imports."""
    text = read(path, 32 * 1024 * 1024).decode("utf-8")
    rows = text.splitlines()
    if not 1 <= len(rows) <= MAX_GENOMES or any(not row or len(row) > 4096 or "\0" in row for row in rows):
        raise InputError("Genome list must contain 1–100000 nonempty paths")
    return [Path(row) if Path(row).is_absolute() else path.parent / row for row in rows]


def prepare(
    inputs: list[Path],
    output: Path,
    name: str,
    threads: int = 1,
    *,
    source: dict | None = None,
    scope: str = "targeted_panel",
    build_timeout: int = 14_400,
) -> dict:
    if not ID.fullmatch(name) or not 1 <= len(inputs) <= MAX_GENOMES:
        raise InputError("Reference needs a bounded name and 1–100000 single-genome FASTAs")
    if (
        type(threads) is not int
        or not 1 <= threads <= 8
        or type(build_timeout) is not int
        or not 60 <= build_timeout <= 14_400
    ):
        raise InputError("Invalid reference build resource limits")
    if scope not in {"targeted_panel", "broad_phage_reference"}:
        raise InputError("Invalid reference scope")
    source = source_metadata(source)
    with transaction(output) as temp:
        (temp / "genomes").mkdir()
        (temp / "db").mkdir()
        (temp / "tmp").mkdir()
        catalogue_path = temp / "catalogue.sqlite3"
        catalogue = sqlite3.connect(catalogue_path)
        catalogue.execute(
            "CREATE TABLE genomes(id TEXT PRIMARY KEY, length INTEGER NOT NULL, sha256 TEXT NOT NULL, file TEXT NOT NULL UNIQUE)"
        )
        seen: set[str] = set()
        total_bases = 0
        archive_path = temp / "tmp/genomes.tar"
        try:
            # Native tar2db preserves one FASTA per genome while avoiding ARG_MAX.
            # USTAR names are generated, bounded and free of operator path details.
            with tarfile.open(archive_path, "w", format=tarfile.USTAR_FORMAT) as archive:
                for i, input_path in enumerate(inputs):
                    row = fasta(input_path, maximum=1, length=MAX_GENOME_BASES, limit=6_000_000)[0]
                    if row["id"] in seen or set(row["sequence"]) - set("ACGTN"):
                        raise InputError("Reference IDs must be unique and DNA containing only ACGTN")
                    seen.add(row["id"])
                    total_bases += len(row["sequence"])
                    if total_bases > MAX_BASES:
                        raise InputError("Reference exceeds 10 billion bases")
                    relative = f"genomes/target{i:06d}.fasta"
                    data = f">{row['id']}\n{row['sequence']}\n".encode("ascii")
                    (temp / relative).write_bytes(data)
                    catalogue.execute(
                        "INSERT INTO genomes VALUES(?,?,?,?)",
                        (row["id"], len(row["sequence"]), digest(row["sequence"].encode()), relative),
                    )
                    member = tarfile.TarInfo(f"target{i:06d}.fasta")
                    member.size = len(data)
                    member.mode = 0o600
                    archive.addfile(member, io.BytesIO(data))
            catalogue.commit()
        finally:
            catalogue.close()
        if catalogue_path.stat().st_size > MAX_CATALOGUE_BYTES:
            raise InputError("Reference catalogue exceeds 128 MiB")
        log = temp / "prepare.log"
        tar_db = temp / "tmp/tardb"
        sequence_db = temp / "tmp/seqdb"
        run(["tar2db", str(archive_path), str(tar_db), "--threads", str(threads)], log, timeout=build_timeout)
        run(["createdb", str(tar_db), str(sequence_db)], log, timeout=build_timeout)
        for prefix, reverse in (("target", 0), ("control", 1)):
            run(
                [
                    "createsetdb",
                    str(sequence_db),
                    str(temp / "db" / prefix),
                    str(temp / "tmp" / prefix),
                    "--reverse-fragments",
                    str(reverse),
                    "--threads",
                    str(threads),
                ],
                log,
                timeout=build_timeout,
            )
        # Native aliases can point within the transaction. Stream them to ordinary
        # files so the atomic rename leaves a self-contained, symlink-free tree.
        for path in (temp / "db").iterdir():
            if path.is_symlink():
                resolved = path.resolve()
                if (
                    temp not in resolved.parents
                    or not resolved.is_file()
                    or resolved.stat().st_size > MAX_INDEX_BYTES
                ):
                    raise InputError("Native reference produced an unsafe alias")
                replacement = path.with_name(path.name + ".materializing")
                with resolved.open("rb") as src, replacement.open("xb") as dest:
                    shutil.copyfileobj(src, dest, 1024 * 1024)
                path.unlink()
                replacement.rename(path)
        shutil.rmtree(temp / "tmp")
        inventory: dict[str, dict] = {}
        total_bytes = 0
        for path in sorted((temp / "db").iterdir()):
            size = path.stat().st_size
            if not path.is_file() or path.is_symlink() or size > MAX_INDEX_BYTES:
                raise InputError("Native reference file is invalid or too large")
            total_bytes += size
            if total_bytes > MAX_INDEX_TOTAL_BYTES or len(inventory) >= 3000:
                raise InputError("Native reference index inventory exceeds budget")
            inventory[str(path.relative_to(temp))] = {
                "sha256": file_digest(path, MAX_INDEX_BYTES),
                "size_bytes": size,
            }
        manifest = {
            "schema_version": SCHEMA,
            "name": name,
            "engine_version": __version__,
            "control_method": "reversed_translated_orfs",
            "genome_count": len(seen),
            "total_bases": total_bases,
            "scope": scope,
            "source": source,
            "catalogue": {
                "file": "catalogue.sqlite3",
                "sha256": file_digest(catalogue_path, MAX_CATALOGUE_BYTES),
                "size_bytes": catalogue_path.stat().st_size,
            },
            "files": inventory,
        }
        write_json(temp / "reference.json", manifest)
    return manifest


def _catalogue(root: Path, manifest: dict) -> sqlite3.Connection:
    receipt = manifest.get("catalogue")
    if (
        not isinstance(receipt, dict)
        or set(receipt) != {"file", "sha256", "size_bytes"}
        or receipt["file"] != "catalogue.sqlite3"
        or not isinstance(receipt["sha256"], str)
        or not SHA.fullmatch(receipt["sha256"])
        or type(receipt["size_bytes"]) is not int
        or not 1 <= receipt["size_bytes"] <= MAX_CATALOGUE_BYTES
    ):
        raise InputError("Invalid reference catalogue receipt")
    path = root / receipt["file"]
    reject_symlinks(path)
    if (
        not path.is_file()
        or path.stat().st_size != receipt["size_bytes"]
        or file_digest(path, MAX_CATALOGUE_BYTES) != receipt["sha256"]
    ):
        raise InputError("Reference catalogue checksum mismatch")
    connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro&immutable=1", uri=True)
    connection.execute("PRAGMA trusted_schema=OFF")
    connection.execute("PRAGMA query_only=ON")
    connection.execute("PRAGMA cache_size=-4096")
    try:
        objects = connection.execute(
            "SELECT type,name,tbl_name FROM sqlite_master WHERE name NOT LIKE 'sqlite_%' LIMIT 2"
        ).fetchall()
        if objects != [("table", "genomes", "genomes")]:
            raise InputError("Invalid reference catalogue schema")
        columns = connection.execute("PRAGMA table_info(genomes)").fetchall()
        if [(row[1], row[2], row[5]) for row in columns] != [
            ("id", "TEXT", 1),
            ("length", "INTEGER", 0),
            ("sha256", "TEXT", 0),
            ("file", "TEXT", 0),
        ]:
            raise InputError("Invalid reference catalogue columns")
        count, bases = 0, 0
        seen_files: set[str] = set()
        for ident, length, sha, file in connection.execute("SELECT id,length,sha256,file FROM genomes"):
            if (
                not isinstance(ident, str)
                or not ID.fullmatch(ident)
                or type(length) is not int
                or not 1 <= length <= MAX_GENOME_BASES
                or not isinstance(sha, str)
                or not SHA.fullmatch(sha)
                or not isinstance(file, str)
                or not re.fullmatch(r"genomes/[A-Za-z0-9_.-]+", file)
                or file in seen_files
            ):
                raise InputError("Invalid reference genome catalogue row")
            count += 1
            bases += length
            seen_files.add(file)
            if count > MAX_GENOMES or bases > MAX_BASES:
                raise InputError("Reference catalogue exceeds bounds")
            genome_path = root / file
            reject_symlinks(genome_path)
            if not genome_path.is_file() or genome_path.stat().st_size > 6_000_000:
                raise InputError("Missing or oversized reference genome")
        if count != manifest["genome_count"] or bases != manifest["total_bases"]:
            raise InputError("Reference catalogue count mismatch")
        if connection.execute("SELECT COUNT(DISTINCT id) FROM genomes").fetchone()[0] != count:
            raise InputError("Duplicate reference genome identifiers")
    except BaseException:
        connection.close()
        raise
    return connection


class ReferenceSequences(Mapping[str, str]):
    """Only matching target genomes are read; at most eight sequences stay cached."""

    def __init__(self, root: Path, connection: sqlite3.Connection, count: int):
        self.root = root
        self.connection = connection
        self.count = count
        self.cache: OrderedDict[str, str] = OrderedDict()

    def __len__(self) -> int:
        return self.count

    def __iter__(self) -> Iterator[str]:
        for row in self.connection.execute("SELECT id FROM genomes ORDER BY id"):
            yield row[0]

    def __contains__(self, key: object) -> bool:
        return (
            isinstance(key, str)
            and self.connection.execute("SELECT 1 FROM genomes WHERE id=?", (key,)).fetchone() is not None
        )

    def __getitem__(self, key: str) -> str:
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        row = self.connection.execute("SELECT length,sha256,file FROM genomes WHERE id=?", (key,)).fetchone()
        if row is None:
            raise KeyError(key)
        sequence = fasta(self.root / row[2], maximum=1, length=MAX_GENOME_BASES, limit=6_000_000)[0]
        if (
            sequence["id"] != key
            or len(sequence["sequence"]) != row[0]
            or digest(sequence["sequence"].encode()) != row[1]
            or set(sequence["sequence"]) - set("ACGTN")
        ):
            raise InputError("Reference genome binding mismatch")
        self.cache[key] = sequence["sequence"]
        if len(self.cache) > 8:
            self.cache.popitem(last=False)
        return sequence["sequence"]

    def close(self) -> None:
        self.connection.close()
        self.cache.clear()

    def __del__(self) -> None:
        self.close()


def load(path: Path) -> tuple[dict, ReferenceSequences]:
    reject_symlinks(path)
    if not path.is_dir():
        raise InputError("Reference must be an operator-owned directory")
    raw = read(path / "reference.json", 2 * 1024 * 1024)
    manifest = json_value(raw)
    if (
        set(manifest) != {"schema_version", *PUBLIC_FIELDS, "catalogue", "files"}
        or manifest.get("schema_version") != SCHEMA
        or not isinstance(manifest.get("engine_version"), str)
        or manifest.get("engine_version") not in REFERENCE_ENGINE_VERSIONS
        or manifest.get("control_method") != "reversed_translated_orfs"
        or not isinstance(manifest.get("name"), str)
        or not ID.fullmatch(manifest["name"])
        or manifest.get("scope") not in {"targeted_panel", "broad_phage_reference"}
        or type(manifest.get("genome_count")) is not int
        or not 1 <= manifest["genome_count"] <= MAX_GENOMES
        or type(manifest.get("total_bases")) is not int
        or not 1 <= manifest["total_bases"] <= MAX_BASES
    ):
        raise InputError("Invalid reference schema or engine identity")
    source_metadata(manifest["source"])
    if not isinstance(manifest.get("files"), dict) or not 1 <= len(manifest["files"]) <= 3000:
        raise InputError("Invalid reference inventory")
    reject_symlinks(path / "db")
    inventory = set()
    for member in (path / "db").iterdir():
        inventory.add(str(member.relative_to(path)))
        if len(inventory) > 3000:
            raise InputError("Reference index inventory exceeds budget")
    required = {
        f"db/{prefix}{suffix}" for prefix in ("target", "control") for suffix in ("", ".index", ".dbtype")
    }
    if inventory != set(manifest["files"]) or not required <= inventory:
        raise InputError("Incomplete reference index inventory")
    total_bytes = 0
    for name, receipt in manifest["files"].items():
        if (
            not isinstance(name, str)
            or not re.fullmatch(r"db/[A-Za-z0-9_.-]+", name)
            or not isinstance(receipt, dict)
            or set(receipt) != {"sha256", "size_bytes"}
            or not isinstance(receipt["sha256"], str)
            or not SHA.fullmatch(receipt["sha256"])
            or type(receipt["size_bytes"]) is not int
            or not 0 <= receipt["size_bytes"] <= MAX_INDEX_BYTES
        ):
            raise InputError("Invalid reference member")
        member = path / name
        reject_symlinks(member)
        if not member.is_file() or member.stat().st_size != receipt["size_bytes"]:
            raise InputError("Reference index size/checksum mismatch")
        total_bytes += receipt["size_bytes"]
        if total_bytes > MAX_INDEX_TOTAL_BYTES:
            raise InputError("Reference indexes exceed 64 GiB")
        if file_digest(member, MAX_INDEX_BYTES) != receipt["sha256"]:
            raise InputError("Reference index checksum mismatch")
    connection = _catalogue(path, manifest)
    manifest["manifest_sha256"] = digest(raw)
    return manifest, ReferenceSequences(path, connection, manifest["genome_count"])


def verify(path: Path) -> dict:
    """Operator preactivation check, including every source sequence."""
    manifest, sequences = load(path)
    try:
        for ident in sequences:
            sequences[ident]
    finally:
        sequences.close()
    return {key: manifest[key] for key in (*PUBLIC_FIELDS, "manifest_sha256")}
