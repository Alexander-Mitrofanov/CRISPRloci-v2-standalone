"""Immutable, disk-indexed spacer collections and bounded exact target searches.

The fixed prefix nominates candidates only; the entire spacer must match. Each
orientation is indexed separately, and palindromes once. No network operations.
"""
from __future__ import annotations

from collections import defaultdict
import gzip
import hashlib
from pathlib import Path
import sqlite3
import time
from urllib.parse import urlsplit
import zlib

from . import __version__
from .artifacts import write_table
from .inputs import validate_records
from .io import (DNA, InputError, LimitError, atomic_directory, canonical, decode_json, digest,
                 identifier, integer, read_bytes, read_fasta_snapshot, regular, verify_manifest,
                 write_json, write_manifest)
from .search import reverse_complement

SEED = 16
MAX_LENGTH = 1000
MAX_DATABASE_BYTES = 4 * 1024**3
SCHEMA = {"name": "crisprspacer.database-search", "version": "1.0.0"}
POLICY = {"method": "full-length-exact", "substitutions": 0, "gaps": 0,
          "strands": "both", "topology": "linear", "array_overlap": "unknown", "pam": "not_evaluated"}
FIELDS = ["spacer_id", "spacer_sequence", "spacer_sha256", "target_id", "target_sha256",
          "target_start", "target_end", "relative_strand", "identities", "query_length",
          "substitutions", "gaps", "query_coverage", "array_overlap", "pam_status", "confidence", "evalue"]
RESULT_FILES = {"analysis.json", "inputs.json", "matches.tsv"}


def file_digest(path, limit=MAX_DATABASE_BYTES, checkpoint=lambda: None):
    path = regular(path)
    total = 0
    result = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            total += len(chunk)
            if total > limit:
                raise LimitError("Database/source exceeds its byte budget")
            result.update(chunk)
            checkpoint()
    return result.hexdigest()


def _records(path):
    """Operator import: streaming bounds larger than interactive FASTA limits."""
    header, pieces, bases, count, total = None, [], 0, 0, 0
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rb") as handle:
            while raw := handle.readline(100_001):
                total += len(raw)
                if len(raw) > 100_000 or total > 2 * 1024**3:
                    raise LimitError("Spacer FASTA exceeds import byte/line budget")
                line = raw.decode("ascii").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        if not bases:
                            raise InputError("Empty spacer record")
                        yield header, "".join(pieces)
                    header = identifier(line[1:].split()[0] if line[1:].split() else None)
                    pieces, bases = [], 0
                    count += 1
                    if count > 5_000_000:
                        raise LimitError("Spacer collection exceeds five million records")
                else:
                    if header is None:
                        raise InputError("Sequence precedes FASTA header")
                    line = line.upper()
                    if set(line) - DNA:
                        raise InputError("Spacer database accepts DNA IUPAC only")
                    bases += len(line)
                    if bases > 100_000:
                        raise LimitError("Spacer source record exceeds import length budget")
                    pieces.append(line)
            if header is None or not bases:
                raise InputError("Empty spacer FASTA/record")
            yield header, "".join(pieces)
    except (UnicodeError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise InputError("Malformed spacer FASTA") from exc


def build_database(source: Path, destination: Path, name: str, source_url: str, source_release: str):
    identifier(name)
    identifier(source_release)
    parsed = urlsplit(source_url)
    if parsed.scheme != "https" or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment or len(source_url) > 1000:
        raise InputError("Source URL must be a public HTTPS URL without credentials/query/fragment")
    source = regular(source)
    with atomic_directory(destination) as work:
        # Import a byte snapshot; do not hash one source version and parse another.
        snapshot = work / ("source.gz" if source.suffix == ".gz" else "source.fasta")
        total = 0
        with source.open("rb") as src, snapshot.open("xb") as dst:
            while chunk := src.read(1024 * 1024):
                total += len(chunk)
                if total > 512 * 1024**2:
                    raise LimitError("Compressed/source import exceeds 512 MiB")
                dst.write(chunk)
        ref = {"name": name, "source_url": source_url, "source_release": source_release,
               "source_file_sha256": file_digest(snapshot), "record_count": 0,
               "eligible_record_count": 0, "skipped_ambiguous": 0, "skipped_length": 0,
               "seed_length": SEED, "min_length": SEED, "max_length": MAX_LENGTH}
        path = work / "spacers.sqlite3"
        connection = sqlite3.connect(path)
        try:
            connection.executescript("""
                PRAGMA journal_mode=DELETE;
                CREATE TABLE metadata (value TEXT NOT NULL);
                CREATE TABLE spacers (id TEXT PRIMARY KEY, sequence TEXT NOT NULL, status TEXT NOT NULL);
                CREATE TABLE seeds (prefix TEXT NOT NULL, spacer_id TEXT NOT NULL,
                                    strand TEXT NOT NULL, oriented TEXT NOT NULL);
            """)
            for key, seq in _records(snapshot):
                ref["record_count"] += 1
                status = "eligible"
                if not SEED <= len(seq) <= MAX_LENGTH:
                    status = "skipped_length"
                elif set(seq) - set("ACGT"):
                    status = "skipped_ambiguous"
                ref["eligible_record_count" if status == "eligible" else status] += 1
                try:
                    connection.execute("INSERT INTO spacers VALUES (?,?,?)", (key, seq, status))
                except sqlite3.IntegrityError as exc:
                    raise InputError("Duplicate spacer FASTA ID") from exc
                if status == "eligible":
                    reverse = reverse_complement(seq)
                    rows = [(seq[:SEED], key, "both" if seq == reverse else "+", seq)]
                    if seq != reverse:
                        rows.append((reverse[:SEED], key, "-", reverse))
                    connection.executemany("INSERT INTO seeds VALUES (?,?,?,?)", rows)
                if ref["record_count"] % 10_000 == 0 and path.stat().st_size > MAX_DATABASE_BYTES:
                    raise LimitError("Index exceeds disk budget")
            connection.execute("CREATE INDEX prefix_lookup ON seeds(prefix)")
            connection.execute("CREATE UNIQUE INDEX spacer_orientation ON seeds(spacer_id,strand)")
            connection.execute("INSERT INTO metadata VALUES (?)", (canonical(ref).decode(),))
            connection.commit()
            _validate_contents(connection, ref)
        finally:
            connection.close()
        snapshot.unlink()
        ref["database_sha256"] = file_digest(path)
        write_json(work / "reference.json", ref)
    return ref


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and not set(value) - set("0123456789abcdef")


def _validate_contents(connection, reference):
    """Stream both ordered inventories to prove no omitted/extra seed rows."""
    objects = set(connection.execute("SELECT type,name FROM sqlite_master WHERE name NOT LIKE 'sqlite_%'"))
    if objects != {("table", "metadata"), ("table", "spacers"), ("table", "seeds"),
                   ("index", "prefix_lookup"), ("index", "spacer_orientation")}:
        raise InputError("Unexpected database schema/indexes")
    expected_columns = {"metadata": ["value"], "spacers": ["id", "sequence", "status"],
                        "seeds": ["prefix", "spacer_id", "strand", "oriented"]}
    for table, columns in expected_columns.items():
        if [r[1] for r in connection.execute("PRAGMA table_info(" + table + ")")] != columns:
            raise InputError("Unexpected database columns")
    # Check SQL lengths before copying any row into Python, including corrupt imports.
    if connection.execute("SELECT 1 FROM spacers WHERE typeof(id)!='text' OR typeof(sequence)!='text' OR typeof(status)!='text' OR length(id)>500 OR length(sequence)>100000 OR length(status)>30 LIMIT 1").fetchone():
        raise InputError("Malformed or oversized spacer row")
    if connection.execute("SELECT 1 FROM seeds WHERE typeof(prefix)!='text' OR typeof(spacer_id)!='text' OR typeof(strand)!='text' OR typeof(oriented)!='text' OR length(prefix)!=16 OR length(spacer_id)>500 OR length(strand)>4 OR length(oriented)>1000 LIMIT 1").fetchone():
        raise InputError("Malformed or oversized seed row")
    counts = {"record_count": 0, "eligible_record_count": 0, "skipped_length": 0, "skipped_ambiguous": 0}
    seed_rows = iter(connection.execute("SELECT spacer_id,strand,prefix,oriented FROM seeds ORDER BY spacer_id,strand"))
    sentinel = object()
    for key, sequence, status in connection.execute("SELECT id,sequence,status FROM spacers ORDER BY id"):
        counts["record_count"] += 1
        if counts["record_count"] > 5_000_000:
            raise LimitError("Database record validation limit exceeded")
        identifier(key)
        if not sequence or set(sequence) - DNA:
            raise InputError("Invalid database spacer sequence")
        expected = "skipped_length" if not SEED <= len(sequence) <= MAX_LENGTH else "skipped_ambiguous" if set(sequence) - set("ACGT") else "eligible"
        if status != expected:
            raise InputError("Spacer eligibility mismatch")
        counts["eligible_record_count" if status == "eligible" else status] += 1
        if status == "eligible":
            reverse = reverse_complement(sequence)
            rows = [(key, "both" if reverse == sequence else "+", sequence[:SEED], sequence)]
            if reverse != sequence:
                rows.append((key, "-", reverse[:SEED], reverse))
            for row in rows:
                if next(seed_rows, sentinel) != row:
                    raise InputError("Spacer seed inventory is incomplete or inconsistent")
    if next(seed_rows, sentinel) is not sentinel or any(reference[k] != v for k, v in counts.items()):
        raise InputError("Spacer database has extra seeds or inconsistent counts")


def validate_database(directory, expected=None):
    connection, reference = open_database(directory, expected)
    try:
        _validate_contents(connection, reference)
    finally:
        connection.close()
    return reference


def open_database(directory: Path, expected=None, checkpoint=lambda: None):
    ref = decode_json(read_bytes(directory / "reference.json", 16_384))
    expected_keys = {"name", "source_url", "source_release", "source_file_sha256", "database_sha256",
                     "record_count", "eligible_record_count", "skipped_ambiguous", "skipped_length",
                     "seed_length", "min_length", "max_length"}
    if not isinstance(ref, dict) or set(ref) != expected_keys or not _sha(ref.get("database_sha256")) or not _sha(ref.get("source_file_sha256")):
        raise InputError("Unsupported spacer database reference")
    for key in ("name", "source_release"):
        identifier(ref[key])
    if not isinstance(ref["source_url"], str) or not ref["source_url"].startswith("https://"):
        raise InputError("Invalid source URL")
    for key in ("record_count", "eligible_record_count", "skipped_ambiguous", "skipped_length"):
        integer(ref[key], key, 0, 5_000_000)
    for key, value in (("seed_length", SEED), ("min_length", SEED), ("max_length", MAX_LENGTH)):
        integer(ref[key], key, value, value)
    if ref["record_count"] != sum(ref[k] for k in ("eligible_record_count", "skipped_ambiguous", "skipped_length")):
        raise InputError("Invalid spacer database eligibility counts/policy")
    if expected is not None and (not _sha(expected) or expected != ref["database_sha256"]):
        raise InputError("Spacer database does not match the pinned digest")
    path = regular(directory / "spacers.sqlite3")
    if file_digest(path, checkpoint=checkpoint) != ref["database_sha256"]:
        raise InputError("Spacer database integrity check failed")
    connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
    try:
        connection.execute("PRAGMA trusted_schema=OFF")
        connection.execute("PRAGMA query_only=ON")
        if connection.execute("SELECT type FROM sqlite_master WHERE name='metadata'").fetchone() != ("table",):
            raise InputError("Database metadata must be a regular table")
        sizes = connection.execute("SELECT typeof(value),length(value) FROM metadata LIMIT 2").fetchall()
        if len(sizes) != 1 or sizes[0][0] != "text" or not 1 <= sizes[0][1] <= 16_384:
            raise InputError("Malformed or oversized database metadata")
        rows = connection.execute("SELECT value FROM metadata LIMIT 2").fetchall()
        bound_reference = {k: v for k, v in ref.items() if k != "database_sha256"}
        if (len(rows) != 1
                or canonical(decode_json(rows[0][0].encode())) != canonical(bound_reference)):
            raise InputError("Database metadata does not bind reference identity")
    except Exception:
        connection.close()
        raise
    return connection, ref


def analyze_database(directory, targets, expected=None, *, timeout=300, max_hits=10_000, max_candidates=2_000_000):
    integer(timeout, "timeout", 1, 3600)
    integer(max_hits, "max_hits", 1, 10_000)
    integer(max_candidates, "max_candidates", 1, 2_000_000)
    if not validate_records(targets):
        raise InputError("At least one target record is required")
    if any(t["topology"] != "linear" for t in targets):
        raise InputError("Database search supports linear targets only")
    if sum(t["length"] for t in targets) > 5_000_000 or len(targets) > 1000:
        raise LimitError("Database search allows at most five million target bases / 1000 records")
    deadline = time.monotonic() + timeout

    def checkpoint():
        if time.monotonic() > deadline:
            raise LimitError("Database search deadline exceeded; no completed result")

    connection, reference = open_database(directory, expected, checkpoint)
    connection.set_progress_handler(lambda: int(time.monotonic() > deadline), 10_000)
    counts = {"matches": 0, "matched_spacers": 0, "targets": len(targets),
              "target_bases": sum(t["length"] for t in targets), "seed_windows": 0,
              "ambiguous_seed_windows": 0, "candidates_checked": 0}
    matches, ids, evidence_bytes = [], set(), 0
    try:
        for target in targets:
            seq = target["sequence"]
            windows = max(0, len(seq) - SEED + 1)
            counts["seed_windows"] += windows
            for begin in range(0, windows, 400):
                checkpoint()
                positions = defaultdict(list)
                for pos in range(begin, min(windows, begin + 400)):
                    prefix = seq[pos:pos+SEED]
                    if set(prefix) - set("ACGT"):
                        counts["ambiguous_seed_windows"] += 1
                    else:
                        positions[prefix].append(pos)
                if not positions:
                    continue
                placeholders = ",".join("?" for _ in positions)
                # Values are bound parameters; only placeholder count is formatted.
                cursor = connection.execute(
                    "SELECT prefix,spacer_id,strand,oriented FROM seeds WHERE prefix IN (" + placeholders + ")",
                    tuple(positions))
                for prefix, key, strand, oriented in cursor:
                    for pos in positions[prefix]:
                        counts["candidates_checked"] += 1
                        if counts["candidates_checked"] > max_candidates:
                            raise LimitError("Database candidate budget exceeded; no completed result")
                        checkpoint()
                        if seq[pos:pos+len(oriented)] != oriented:
                            continue
                        if len(matches) >= max_hits:
                            raise LimitError("Database hit budget exceeded; no completed result")
                        spacer = reverse_complement(oriented) if strand == "-" else oriented
                        row = {"spacer_id": key, "spacer_sequence": spacer, "spacer_sha256": digest(spacer),
                               "target_id": target["id"], "target_sha256": target["sha256"],
                               "target_start": pos+1, "target_end": pos+len(spacer), "relative_strand": strand,
                               "identities": len(spacer), "query_length": len(spacer), "substitutions": 0,
                               "gaps": 0, "query_coverage": 1.0, "array_overlap": "unknown",
                               "pam_status": "not_evaluated", "confidence": None, "evalue": None}
                        evidence_bytes += len(canonical(row))
                        if evidence_bytes > 16 * 1024**2:
                            raise LimitError("Database evidence byte budget exceeded; no completed result")
                        matches.append(row)
                        ids.add(key)
    except sqlite3.OperationalError as exc:
        if time.monotonic() > deadline:
            raise LimitError("Database SQL deadline exceeded; no completed result") from exc
        raise
    finally:
        connection.close()
    checkpoint()
    matches.sort(key=lambda r: (r["target_id"], r["target_start"], r["target_end"], r["spacer_id"], r["relative_strand"]))
    counts.update(matches=len(matches), matched_spacers=len(ids))
    assessable = counts["seed_windows"] > counts["ambiguous_seed_windows"] and reference["eligible_record_count"] > 0
    return {"schema": SCHEMA, "tool": {"name": "CRISPRspacer", "version": __version__},
            "status": "completed", "search_complete": True,
            "outcome": "matches" if matches else "no_matches" if assessable else "no_assessable_windows",
            "reference": reference, "policy": POLICY,
            "targets": [{k: t[k] for k in ("id", "length", "sha256")} for t in targets],
            "counts": counts, "matches": matches}


def search_database(database, target_path, destination, expected=None, **limits):
    targets, source_hash = read_fasta_snapshot(target_path)
    deadline = time.monotonic() + limits.get("timeout", 300)
    with atomic_directory(destination) as work:
        report = analyze_database(database, targets, expected, **limits)
        write_json(work / "analysis.json", report)
        write_json(work / "inputs.json", {"schema_version": "1.0.0", "targets": targets,
                                         "target_file_sha256": source_hash, "reference": report["reference"]})
        write_table(work / "matches.tsv", report["matches"], FIELDS)
        write_manifest(work, "crisprspacer.database-search")
        if time.monotonic() > deadline:
            raise LimitError("Database publication deadline exceeded")
    return report


def validate_database_result(result, database):
    _, payloads = verify_manifest(result, "crisprspacer.database-search", RESULT_FILES)
    inputs = decode_json(payloads["inputs.json"])
    report = decode_json(payloads["analysis.json"])
    if not isinstance(inputs, dict) or set(inputs) != {"schema_version", "targets", "target_file_sha256", "reference"} or inputs["schema_version"] != "1.0.0" or not _sha(inputs["target_file_sha256"]):
        raise InputError("Invalid retained database-search inputs")
    validate_records(inputs["targets"])
    replay = analyze_database(database, inputs["targets"], inputs["reference"]["database_sha256"])
    if not isinstance(report, dict) or report.get("tool") not in [
            {"name": "CRISPRspacer", "version": version} for version in ("0.2.0", "0.3.0", __version__)]:
        raise InputError("Unsupported database-search tool identity")
    replay["tool"]["version"] = report["tool"]["version"]
    if canonical(replay) != canonical(report) or canonical(report["reference"]) != canonical(inputs["reference"]):
        raise InputError("Database search replay differs from retained evidence")
    import tempfile
    with tempfile.TemporaryDirectory() as scratch:
        table = Path(scratch) / "matches.tsv"
        write_table(table, report["matches"], FIELDS)
        if table.read_bytes() != payloads["matches.tsv"]:
            raise InputError("Database TSV differs from authoritative JSON")
    return {"status": "valid", "matches": report["counts"]["matches"], "database_sha256": report["reference"]["database_sha256"]}
