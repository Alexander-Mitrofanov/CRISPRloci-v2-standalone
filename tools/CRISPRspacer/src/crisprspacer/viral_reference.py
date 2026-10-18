"""Operator-built, source-bound local viral BLAST references; no network access."""
from __future__ import annotations

from contextlib import contextmanager
import gzip
import hashlib
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
import time
from urllib.parse import urlsplit
import zlib

from .blast_process import executable_identity, run_blast
from .io import (DNA, InputError, LimitError, atomic_directory, canonical, decode_json,
                 identifier, integer, read_bytes, regular, write_json)

SCHEMA = {"name": "crisprspacer.viral-reference", "version": "1.0.0"}
MAX_CORPUS_BYTES = 20 * 1024**3
MAX_BUNDLE_BYTES = 80 * 1024**3
MAX_RECORDS = 1_000_000
MAX_SEQUENCE_BASES = 20_000_000
MAX_SOURCES = 1000
MAX_FILES = 4096
MIN_FREE_BYTES = 512 * 1024**2
_NATIVE = re.compile(r"db(?:\.[0-9]+)?\.(?:ndb|nhr|nin|nog|nos|not|nsq|ntf|nto|njs|nal)\Z")
_HEX = re.compile(r"[0-9a-f]{64}\Z")


def _sha(value):
    return isinstance(value, str) and _HEX.fullmatch(value) is not None


def _https(value):
    if not isinstance(value, str) or len(value) > 2000:
        raise InputError("Invalid reference source URL")
    parsed = urlsplit(value)
    if parsed.scheme != "https" or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise InputError("Reference source URL must be HTTPS without credentials/query/fragment")
    return value


def _file_identity(path, checkpoint=lambda: None, limit=MAX_BUNDLE_BYTES):
    path = regular(path)
    sha, total = hashlib.sha256(), 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            total += len(chunk)
            if total > limit:
                raise LimitError("Viral reference file exceeds byte budget")
            sha.update(chunk)
            checkpoint()
    return {"sha256": sha.hexdigest(), "bytes": total}


def _deadline(timeout):
    integer(timeout, "timeout", 1, 86400)
    deadline = time.monotonic() + timeout

    def checkpoint():
        if time.monotonic() >= deadline:
            raise LimitError("Viral reference operation exceeded its deadline")

    def remaining():
        checkpoint()
        return max(1, int(deadline - time.monotonic()))

    return checkpoint, remaining


def _source_snapshot(source, snapshot, checkpoint, remaining_bytes):
    sha, md5, total = hashlib.sha256(), hashlib.md5(), 0
    with regular(source).open("rb") as src, snapshot.open("xb") as dst:
        while chunk := src.read(1024 * 1024):
            total += len(chunk)
            if total > remaining_bytes:
                raise LimitError("Viral source corpus exceeds compressed byte budget")
            if shutil.disk_usage(snapshot.parent).free < MIN_FREE_BYTES + len(chunk):
                raise LimitError("Insufficient free disk space for viral reference import")
            sha.update(chunk)
            md5.update(chunk)
            dst.write(chunk)
            checkpoint()
    return {"name": source.name, "bytes": total, "sha256": sha.hexdigest(), "md5": md5.hexdigest()}


def _records(path, byte_counter, checkpoint):
    """Keep at most one bounded genome in memory; corpus limits are aggregate."""
    current, pieces, bases = None, [], 0
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rb") as handle:
            while raw := handle.readline(1_000_001):
                byte_counter[0] += len(raw)
                if len(raw) > 1_000_000 or byte_counter[0] > MAX_CORPUS_BYTES:
                    raise LimitError("Viral FASTA exceeds line/decompressed corpus byte budget")
                checkpoint()
                line = raw.decode("ascii").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if len(line) > 2000 or any(ord(c) < 32 or ord(c) > 126 for c in line):
                        raise InputError("Viral FASTA header exceeds printable 2000-character limit")
                    if current is not None:
                        if not bases:
                            raise InputError("Empty viral FASTA record")
                        yield (*current, "".join(pieces))
                    parts = line[1:].split(maxsplit=1)
                    accession = identifier(parts[0] if parts else None, "viral accession")
                    current = accession, parts[1] if len(parts) == 2 else ""
                    pieces, bases = [], 0
                else:
                    if current is None:
                        raise InputError("Viral FASTA sequence precedes its header")
                    line = line.upper()
                    if not set(line) <= DNA:
                        raise InputError("Viral FASTA accepts DNA IUPAC letters only")
                    bases += len(line)
                    if bases > MAX_SEQUENCE_BASES:
                        raise LimitError("Viral FASTA record exceeds sequence length budget")
                    pieces.append(line)
            if current is None or not bases:
                raise InputError("Empty viral FASTA/source record")
            yield (*current, "".join(pieces))
    except (UnicodeError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise InputError("Malformed or corrupt viral FASTA") from exc


def _normalized_chunks(key, sequence):
    yield f">{key}\n".encode("ascii")
    for start in range(0, len(sequence), 80):
        yield (sequence[start:start + 80] + "\n").encode("ascii")


def _content_row(row):
    return {key: row[key] for key in ("id", "accession", "title", "length", "sha256")}


def _source_receipt(path, source_files, source_release, source_url):
    if path is None:
        return None
    raw = read_bytes(path, 2 * 1024**2)
    receipt = decode_json(raw)
    if (not isinstance(receipt, dict)
            or receipt.get("schema") != {"name": "crisprspacer.ncbi-viral-download", "version": "1.0.0"}
            or receipt.get("source_release") != source_release or receipt.get("source_url") != source_url
            or not isinstance(receipt.get("files"), list) or len(receipt["files"]) != len(source_files)):
        raise InputError("Source download receipt does not match this viral reference")
    expected = {}
    for item in receipt["files"]:
        if not isinstance(item, dict) or set(item) != {"name", "url", "bytes", "md5", "sha256"} or item["name"] in expected:
            raise InputError("Malformed source receipt file inventory")
        expected[item["name"]] = item
    for item in source_files:
        upstream = expected.get(item["name"])
        if upstream is None or any(upstream[k] != item[k] for k in ("name", "bytes", "md5", "sha256")):
            raise InputError("Source download receipt byte/hash mismatch")
        if "url" in item and item["url"] != upstream["url"]:
            raise InputError("Source download receipt URL mismatch")
        item["url"] = _https(upstream["url"])
    return raw


def _check_aliases(directory, files):
    for name in files:
        if not name.endswith(".nal"):
            continue
        try:
            text = read_bytes(directory / name, 1024 * 1024).decode("ascii")
        except UnicodeError as exc:
            raise InputError("Malformed BLAST database alias") from exc
        seen = False
        for line in text.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            key, _, value = line.partition(" ")
            if key == "DBLIST":
                volumes = value.split()
                if not volumes or any(not re.fullmatch(r"db\.[0-9]+", v) or v + ".nin" not in files for v in volumes):
                    raise InputError("BLAST database alias escapes the reference bundle")
                seen = True
            elif key not in {"TITLE", "NSEQ", "LENGTH"}:
                raise InputError("Unsupported BLAST database alias directive")
        if not seen:
            raise InputError("BLAST database alias has no local volume inventory")


def build_viral_reference(sources: list[Path], destination: Path, *, name, source_release, source_url,
                          makeblastdb="makeblastdb", blastdbcmd="blastdbcmd", source_manifest=None, timeout=3600):
    """Build a new local snapshot; callers provision/download sources separately."""
    identifier(name)
    identifier(source_release)
    _https(source_url)
    if not isinstance(sources, (list, tuple)) or not 1 <= len(sources) <= MAX_SOURCES:
        raise InputError("Viral reference requires 1–1000 local FASTA sources")
    sources = [Path(s) for s in sources]
    if len({s.name for s in sources}) != len(sources):
        raise InputError("Viral source basenames must be unique")
    checkpoint, remaining = _deadline(timeout)
    identities = {"makeblastdb": executable_identity(str(makeblastdb)), "blastdbcmd": executable_identity(str(blastdbcmd))}
    logical = ["makeblastdb", "-in", "subjects.fasta", "-dbtype", "nucl", "-parse_seqids",
               "-blastdb_version", "5", "-max_file_sz", "1GB", "-title", name, "-out", "db"]
    source_files, compressed, uncompressed, count, bases = [], 0, [0], 0, 0
    content_hash = hashlib.sha256()
    with atomic_directory(destination) as work:
        connection = sqlite3.connect(work / "subjects.sqlite3")
        try:
            connection.execute("PRAGMA journal_mode=DELETE")
            connection.execute("CREATE TABLE subjects (id TEXT PRIMARY KEY, accession TEXT UNIQUE NOT NULL, title TEXT NOT NULL, sequence TEXT NOT NULL, length INTEGER NOT NULL, sha256 TEXT NOT NULL)")
            with (work / "subjects.fasta").open("xb") as fasta:
                for number, source in enumerate(sources):
                    snapshot = work / (f"source-{number}.gz" if source.suffix == ".gz" else f"source-{number}.fasta")
                    info = _source_snapshot(source, snapshot, checkpoint, MAX_CORPUS_BYTES - compressed)
                    compressed += info["bytes"]
                    source_files.append(info)
                    for accession, title, sequence in _records(snapshot, uncompressed, checkpoint):
                        count += 1
                        if count > MAX_RECORDS:
                            raise LimitError("Viral reference exceeds record budget")
                        if shutil.disk_usage(work).free < MIN_FREE_BYTES + 3 * len(sequence):
                            raise LimitError("Insufficient free disk space for viral reference import")
                        bases += len(sequence)
                        row = {"id": f"v{count:09d}", "accession": accession, "title": title,
                               "sequence": sequence, "length": len(sequence), "sha256": hashlib.sha256(sequence.encode()).hexdigest()}
                        try:
                            connection.execute("INSERT INTO subjects VALUES (:id,:accession,:title,:sequence,:length,:sha256)", row)
                        except sqlite3.IntegrityError as exc:
                            raise InputError("Duplicate viral accession in source collection") from exc
                        content_hash.update(canonical(_content_row(row)))
                        for chunk in _normalized_chunks(row["id"], sequence):
                            fasta.write(chunk)
                    snapshot.unlink()
            connection.commit()
        finally:
            connection.close()
        receipt_raw = _source_receipt(source_manifest, source_files, source_release, source_url)
        if receipt_raw is not None:
            (work / "source-manifest.json").write_bytes(receipt_raw)
        if shutil.disk_usage(work).free < MIN_FREE_BYTES + 2 * bases + count * 1024:
            raise LimitError("Insufficient free disk space for native viral reference index")
        run_blast([identities["makeblastdb"]["path"], *logical[1:]], cwd=work,
                  stdout=work / "build.stdout", stderr=work / "build.stderr", timeout=remaining(), max_output_bytes=1024**2)
        (work / "build.stdout").unlink()
        (work / "build.stderr").unlink()
        files = {}
        total_bytes = 0
        for path in sorted(work.iterdir()):
            if path.name not in {"subjects.fasta", "subjects.sqlite3", "source-manifest.json"} and not _NATIVE.fullmatch(path.name):
                raise InputError("Unexpected native BLAST build output")
            files[path.name] = _file_identity(path, checkpoint)
            total_bytes += files[path.name]["bytes"]
            if total_bytes > MAX_BUNDLE_BYTES or len(files) > MAX_FILES:
                raise LimitError("Viral reference exceeds aggregate output budget")
        if "db.nin" not in files and "db.nal" not in files:
            raise InputError("makeblastdb did not create a nucleotide database")
        _check_aliases(work, files)
        reference = {"schema": SCHEMA, "name": name, "source_release": source_release, "source_url": source_url,
                     "record_count": count, "total_bases": bases, "source_files": source_files,
                     "source_manifest": "source-manifest.json" if receipt_raw is not None else None,
                     "executables": identities, "build": {"command": logical, "database_format": 5,
                     "normalization": "uppercase DNA IUPAC; input file/record order; generated IDs; 80-base lines",
                     "topology": "linear", "taxonomy": "not_supplied"},
                     "content_sha256": content_hash.hexdigest(), "files": files}
        write_json(work / "reference.json", reference)
        validate_viral_reference(work, blastdbcmd=identities["blastdbcmd"]["path"], timeout=remaining())
        checkpoint()
    return reference


def load_viral_reference(directory: Path, expected=None, checkpoint=lambda: None):
    """Check a complete content-addressed snapshot before any search subprocess."""
    directory = Path(directory).absolute()
    raw = read_bytes(directory / "reference.json", 4 * 1024**2)
    reference_sha = hashlib.sha256(raw).hexdigest()
    if expected is not None and (not _sha(expected) or expected != reference_sha):
        raise InputError("Viral reference does not match the pinned digest")
    ref = decode_json(raw)
    required = {"schema", "name", "source_release", "source_url", "record_count", "total_bases", "source_files",
                "source_manifest", "executables", "build", "content_sha256", "files"}
    if not isinstance(ref, dict) or set(ref) != required or ref["schema"] != SCHEMA or not _sha(ref["content_sha256"]):
        raise InputError("Unsupported viral reference schema")
    identifier(ref["name"])
    identifier(ref["source_release"])
    _https(ref["source_url"])
    integer(ref["record_count"], "viral reference record count", 1, MAX_RECORDS)
    integer(ref["total_bases"], "viral reference base count", 1, MAX_CORPUS_BYTES)
    if not isinstance(ref["source_files"], list) or not 1 <= len(ref["source_files"]) <= MAX_SOURCES:
        raise InputError("Invalid viral source inventory")
    seen, source_bytes = set(), 0
    for info in ref["source_files"]:
        if not isinstance(info, dict) or set(info) not in ({"name", "bytes", "sha256", "md5"}, {"name", "bytes", "sha256", "md5", "url"}):
            raise InputError("Malformed viral source identity")
        key = info["name"]
        if not isinstance(key, str) or Path(key).name != key or key in seen or not key or not _sha(info["sha256"]) or not isinstance(info["md5"], str) or not re.fullmatch(r"[0-9a-f]{32}", info["md5"]):
            raise InputError("Invalid or duplicate viral source identity")
        seen.add(key)
        integer(info["bytes"], "source bytes", 1, MAX_CORPUS_BYTES)
        source_bytes += info["bytes"]
        if "url" in info:
            _https(info["url"])
    if source_bytes > MAX_CORPUS_BYTES:
        raise LimitError("Viral reference source byte budget exceeded")
    if not isinstance(ref["executables"], dict) or set(ref["executables"]) != {"makeblastdb", "blastdbcmd"}:
        raise InputError("Invalid viral reference executable receipt")
    for identity in ref["executables"].values():
        if (not isinstance(identity, dict) or set(identity) != {"path", "sha256", "version"}
                or not isinstance(identity["path"], str) or not Path(identity["path"]).is_absolute()
                or not _sha(identity["sha256"]) or not isinstance(identity["version"], str)
                or not 1 <= len(identity["version"]) <= 8192):
            raise InputError("Malformed viral reference executable identity")
    if (not isinstance(ref["build"], dict) or ref["build"].get("database_format") != 5
            or ref["build"].get("topology") != "linear" or ref["build"].get("taxonomy") != "not_supplied"):
        raise InputError("Unsupported viral reference build policy")
    files = ref["files"]
    if not isinstance(files, dict) or not 3 <= len(files) <= MAX_FILES or not {"subjects.fasta", "subjects.sqlite3"} <= set(files):
        raise InputError("Invalid viral reference file inventory")
    if ref["source_manifest"] not in (None, "source-manifest.json") or (ref["source_manifest"] is not None) != ("source-manifest.json" in files):
        raise InputError("Invalid retained source receipt identity")
    if {p.name for p in directory.iterdir()} != set(files) | {"reference.json"}:
        raise InputError("Unexpected or missing files in viral reference")
    total = 0
    for name, identity in files.items():
        if name not in {"subjects.fasta", "subjects.sqlite3", "source-manifest.json"} and not _NATIVE.fullmatch(name):
            raise InputError("Unsafe or unexpected viral reference file")
        if not isinstance(identity, dict) or set(identity) != {"sha256", "bytes"} or not _sha(identity["sha256"]):
            raise InputError("Invalid viral reference file identity")
        integer(identity["bytes"], "reference file bytes", 0, MAX_BUNDLE_BYTES)
        total += identity["bytes"]
        if total > MAX_BUNDLE_BYTES:
            raise LimitError("Viral reference aggregate byte budget exceeded")
        if _file_identity(directory / name, checkpoint, identity["bytes"]) != identity:
            raise InputError("Viral reference file integrity check failed")
    if "db.nin" not in files and "db.nal" not in files:
        raise InputError("Viral reference has no nucleotide BLAST database")
    _check_aliases(directory, files)
    if ref["source_manifest"]:
        _source_receipt(directory / ref["source_manifest"], [dict(x) for x in ref["source_files"]], ref["source_release"], ref["source_url"])
    checkpoint()
    return ref, reference_sha


@contextmanager
def bounded_catalog(connection, checkpoint):
    """Keep SQL scans inside the caller's deadline and cancellation policy.

    SQLite swallows exceptions raised by progress callbacks. Retain the original
    exception so a resource limit or cancellation is reported with its real cause.
    """
    if connection is None:
        yield
        return
    interrupted = []

    def progress():
        try:
            checkpoint()
        except BaseException as exc:
            interrupted.append(exc)
            return 1
        return 0

    checkpoint()
    connection.set_progress_handler(progress, 10_000)
    try:
        yield
    except sqlite3.OperationalError as exc:
        if interrupted:
            raise interrupted[0] from exc
        raise
    finally:
        connection.set_progress_handler(None, 0)


def open_subjects(directory: Path, checkpoint=lambda: None):
    path = regular(Path(directory) / "subjects.sqlite3")
    connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        with bounded_catalog(connection, checkpoint):
            connection.execute("PRAGMA trusted_schema=OFF")
            connection.execute("PRAGMA query_only=ON")
            objects = {(r[0], r[1]) for r in connection.execute("SELECT type,name FROM sqlite_master WHERE name NOT LIKE 'sqlite_%'")}
            columns = [r[1] for r in connection.execute("PRAGMA table_info(subjects)")]
            if objects != {("table", "subjects")} or columns != ["id", "accession", "title", "sequence", "length", "sha256"]:
                raise InputError("Unexpected viral catalog schema")
            if connection.execute("SELECT 1 FROM subjects WHERE typeof(id)!='text' OR length(id)!=10 OR typeof(accession)!='text' OR length(accession)>500 OR typeof(title)!='text' OR length(title)>2000 OR typeof(sequence)!='text' OR length(sequence)>? OR typeof(length)!='integer' OR length<1 OR length!=length(sequence) OR typeof(sha256)!='text' OR length(sha256)!=64 LIMIT 1", (MAX_SEQUENCE_BASES,)).fetchone():
                raise InputError("Malformed or oversized viral catalog row")
    except BaseException:
        connection.close()
        raise
    return connection


def validate_viral_reference(directory: Path, expected=None, *, blastdbcmd=None, timeout=3600):
    """Deeply bind the catalog, normalized FASTA and every BLAST subject sequence."""
    directory = Path(directory).absolute()
    checkpoint, remaining = _deadline(timeout)
    ref, ref_sha = load_viral_reference(directory, expected, checkpoint)
    executable = executable_identity(str(blastdbcmd or ref["executables"]["blastdbcmd"]["path"]))
    bound = ref["executables"]["blastdbcmd"]
    if executable["sha256"] != bound["sha256"] or executable["version"] != bound["version"]:
        raise InputError("blastdbcmd does not match the reference build receipt")
    connection = open_subjects(directory, checkpoint)
    try:
        content_hash, fasta_hash, count, bases = hashlib.sha256(), hashlib.sha256(), 0, 0
        with bounded_catalog(connection, checkpoint):
            for row in connection.execute("SELECT * FROM subjects ORDER BY id"):
                checkpoint()
                count += 1
                if count > MAX_RECORDS or row["id"] != f"v{count:09d}":
                    raise InputError("Invalid viral catalog record order")
                identifier(row["accession"], "viral accession")
                if any(ord(c) < 32 or ord(c) > 126 for c in row["title"]):
                    raise InputError("Invalid viral catalog title")
                sequence = row["sequence"]
                if set(sequence) - DNA or hashlib.sha256(sequence.encode()).hexdigest() != row["sha256"]:
                    raise InputError("Viral catalog sequence/hash mismatch")
                bases += row["length"]
                content_hash.update(canonical(_content_row(row)))
                for chunk in _normalized_chunks(row["id"], sequence):
                    fasta_hash.update(chunk)
        if (count != ref["record_count"] or bases != ref["total_bases"] or content_hash.hexdigest() != ref["content_sha256"]
                or fasta_hash.hexdigest() != ref["files"]["subjects.fasta"]["sha256"]):
            raise InputError("Viral reference normalized content/counts mismatch")
        if shutil.disk_usage(directory.parent).free < MIN_FREE_BYTES + bases + count * 32:
            raise LimitError("Insufficient free disk space for BLAST reference roundtrip")
        with tempfile.TemporaryDirectory(prefix=".crisprspacer-roundtrip-", dir=directory.parent) as tmp:
            work = Path(tmp)
            output, stderr = work / "subjects.tsv", work / "stderr"
            run_blast([executable["path"], "-db", str(directory / "db"), "-entry", "all", "-outfmt", "%i\t%s"],
                      cwd=work, stdout=output, stderr=stderr, timeout=remaining(), max_output_bytes=MAX_CORPUS_BYTES)
            with bounded_catalog(connection, checkpoint), output.open("rb") as handle:
                expected_rows = iter(connection.execute("SELECT id,sequence FROM subjects ORDER BY id"))
                while raw := handle.readline(MAX_SEQUENCE_BASES + 100):
                    checkpoint()
                    if len(raw) > MAX_SEQUENCE_BASES + 64 or not raw.endswith(b"\n"):
                        raise InputError("Oversized/malformed BLAST subject roundtrip")
                    parts = raw.rstrip(b"\r\n").decode("ascii").split("\t")
                    row = next(expected_rows, None)
                    if len(parts) != 2 or row is None or parts[0] not in (row["id"], "lcl|" + row["id"]) or parts[1] != row["sequence"]:
                        raise InputError("BLAST subject roundtrip differs from normalized catalog")
            if next(expected_rows, None) is not None:
                raise InputError("BLAST subject roundtrip omitted catalog records")
    except (sqlite3.Error, UnicodeError) as exc:
        raise InputError("Invalid viral subject catalog") from exc
    finally:
        connection.close()
    checkpoint()
    return {"status": "valid", "reference": ref, "reference_sha256": ref_sha}
