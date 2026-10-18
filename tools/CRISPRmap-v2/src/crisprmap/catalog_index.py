"""Disk-backed, lossless partition candidates for short global edit distance.

For a reference split into p nonempty segments, fewer than p unit edits leave
at least one segment intact. Prefix and suffix length differences lower-bound
the edits outside that segment. Both orientations and every feasible offset
under this necessary bound are enumerated.
Candidates are evidence to verify, not matches or a heuristic top-hit list.
"""

from __future__ import annotations

from contextlib import contextmanager
from fractions import Fraction
import hashlib
import math
from pathlib import Path
import shutil
import sqlite3
import time

from .core import reverse_complement
from .io import InputError, LimitError, atomic_json, integer, read_json

SCHEMA = {"name": "crisprmap.catalog-index", "version": "2.0.0"}
MAX_LENGTH = 200
MAX_DATABASE_BYTES = 256 * 1024**3
CORE_PREDICATE = "eligible_and_bona_fide_consensus_count_positive"
COUNT_FIELDS = ("bona_fide_consensus_count", "possible_consensus_count",
                "bona_fide_repeat_count", "possible_repeat_count")
NODE_COLUMNS = "node_id,sequence,sequence_sha256,length,eligible," + ",".join(COUNT_FIELDS)


def file_sha256(path: Path) -> str:
    path = Path(path)
    if not path.is_file() or any(p.is_symlink() for p in (path, *path.parents)):
        raise InputError("Expected a regular file without symlink traversal")
    result = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            result.update(chunk)
    return result.hexdigest()


def similarity_fraction(value) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise InputError("Similarity must be a number between 0.5 and 1")
    try:
        result = Fraction(str(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise InputError("Invalid minimum similarity") from exc
    if not Fraction(1, 2) <= result <= 1 or result.denominator > 10000:
        raise InputError("Similarity must be in [0.5, 1] with at most four decimal places")
    return result


def distance_threshold(a: int, b: int, similarity) -> int:
    similarity = similarity if isinstance(similarity, Fraction) else similarity_fraction(similarity)
    return (max(a, b) * (similarity.denominator - similarity.numerator)) // similarity.denominator


def partition(length: int, distance: int) -> list[tuple[int, int, int]]:
    """Return ordinal, start and size of balanced, nonempty reference segments."""
    count = min(length, distance + 1)
    return [(ordinal, ordinal * length // count,
             (ordinal + 1) * length // count - ordinal * length // count) for ordinal in range(count)]


def _profiles(length: int, maximum: int, similarity: Fraction) -> list[int]:
    # A qualifying longer query M obeys M <= L/similarity. Therefore this
    # threshold upper-bounds every supported relative-distance comparison to L.
    similarities = {similarity, *(value for value in (Fraction(9, 10), Fraction(19, 20)) if value >= similarity)}
    bounds = [distance_threshold(length, min(MAX_LENGTH, length * value.denominator // value.numerator), value)
              for value in similarities]
    return sorted({min(length - 1, value) for value in (min(3, maximum), maximum, *bounds)})


def _policy(manifest):
    try:
        policy = manifest["policy"]
        maximum = integer(policy["max_distance"], "indexed max_distance", 0, MAX_LENGTH)
        fraction = policy["min_similarity_fraction"]
        if (not isinstance(fraction, list) or len(fraction) != 2
                or any(type(n) is not int for n in fraction)):
            raise ValueError
        similarity = Fraction(*fraction)
        if not Fraction(1, 2) <= similarity <= 1 or policy["max_sequence_length"] != MAX_LENGTH:
            raise ValueError
        profiles = [[value.numerator, value.denominator] for value in sorted(
            {similarity, *(value for value in (Fraction(9, 10), Fraction(19, 20)) if value >= similarity)})]
        if (policy.get("relative_profile_fractions") != profiles
                or policy.get("core_predicate") != CORE_PREDICATE
                or policy.get("method") != "balanced_partition_all_offsets_both_query_orientations"):
            raise ValueError
        return maximum, similarity
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise InputError("Unsupported candidate-index policy") from exc


def _readonly(path):
    connection = sqlite3.connect(Path(path).absolute().as_uri() + "?mode=ro&immutable=1", uri=True, timeout=5)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA trusted_schema=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-32768")
    return connection


def _check_free(directory, minimum):
    if shutil.disk_usage(directory).free < minimum:
        raise LimitError("Index build stopped before exhausting the configured free-space reserve")


def build_index(catalog_dir, destination, max_distance=6, min_similarity=0.85,
                *, max_database_bytes=MAX_DATABASE_BYTES, max_nodes=100_000_000,
                deadline_seconds=172800, min_free_bytes=2 * 1024**3) -> dict:
    """Build a new immutable candidate index; publish its manifest last."""
    from .catalog import load_catalog

    integer(max_distance, "max_distance", 0, MAX_LENGTH)
    integer(max_database_bytes, "max_database_bytes", 1024 * 1024, MAX_DATABASE_BYTES)
    integer(max_nodes, "max_nodes", 1, 100_000_000)
    integer(min_free_bytes, "min_free_bytes", 0, MAX_DATABASE_BYTES)
    similarity = similarity_fraction(min_similarity)
    if type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds) or deadline_seconds <= 0:
        raise InputError("Deadline must be positive and finite")
    started = time.monotonic()
    catalog_dir, destination = Path(catalog_dir).absolute(), Path(destination).absolute()
    if any(p.is_symlink() for p in (destination, *destination.parents)):
        raise InputError("Output paths must not traverse symlinks")
    catalog, source = load_catalog(catalog_dir, verify_hash=True)
    database = None
    owned = False
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        _check_free(destination.parent, min_free_bytes)
        destination.mkdir()
        owned = True
        database_path = destination / "index.sqlite3"
        database = sqlite3.connect(database_path)
        database.execute("PRAGMA page_size=4096")
        database.execute("PRAGMA journal_mode=DELETE")
        database.execute("PRAGMA synchronous=FULL")
        database.execute("PRAGMA temp_store=FILE")
        database.execute("PRAGMA cache_size=-65536")
        database.execute(f"PRAGMA max_page_count={max_database_bytes // 4096}")
        database.executescript("""
            CREATE TABLE profiles(length INTEGER NOT NULL,k INTEGER NOT NULL,ordinal INTEGER NOT NULL,
              start INTEGER NOT NULL,size INTEGER NOT NULL,PRIMARY KEY(length,k,ordinal)) WITHOUT ROWID;
            CREATE TABLE postings(length INTEGER NOT NULL,k INTEGER NOT NULL,ordinal INTEGER NOT NULL,
              segment TEXT NOT NULL,is_core INTEGER NOT NULL CHECK(is_core IN (0,1)),node_id INTEGER NOT NULL,
              PRIMARY KEY(length,k,ordinal,segment,is_core,node_id)) WITHOUT ROWID;
            CREATE INDEX postings_node ON postings(node_id);
        """)
        counts = {"nodes": 0, "postings": 0, "profiles": 0, "core_nodes": 0, "core_postings": 0}
        profiles_by_length = {}
        for row in catalog.execute("SELECT node_id,sequence,length,bona_fide_consensus_count "
                                   "FROM nodes WHERE eligible=1 ORDER BY node_id"):
            node_id, sequence, length, core_count = row
            is_core = int(core_count > 0)
            if not 1 <= length <= MAX_LENGTH or len(sequence) != length or set(sequence) - set("ACGT"):
                raise InputError("Eligible catalog node is not a bounded unambiguous sequence")
            if time.monotonic() - started > deadline_seconds:
                raise LimitError("Index build exceeded deadline")
            counts["nodes"] += 1
            if counts["nodes"] > max_nodes:
                raise LimitError("Index build exceeded node budget")
            if length not in profiles_by_length:
                profiles = [(k, partition(length, k)) for k in _profiles(length, max_distance, similarity)]
                for k, parts in profiles:
                    database.executemany("INSERT INTO profiles VALUES (?,?,?,?,?)",
                                         [(length, k, ordinal, start, size) for ordinal, start, size in parts])
                    counts["profiles"] += len(parts)
                profiles_by_length[length] = profiles
            profiles = profiles_by_length[length]
            counts["core_nodes"] += is_core
            for k, parts in profiles:
                rows = [(length, k, ordinal, sequence[start:start + size], is_core, node_id)
                        for ordinal, start, size in parts]
                database.executemany("INSERT INTO postings VALUES (?,?,?,?,?,?)", rows)
                counts["postings"] += len(rows)
                counts["core_postings"] += is_core * len(rows)
            if counts["nodes"] % 10000 == 0:
                database.commit()
                _check_free(destination, min_free_bytes)
        database.commit()
        _check_free(destination, min_free_bytes)
        database.execute("ANALYZE")
        database.commit()
        database.close()
        database = None
        manifest = {"schema": dict(SCHEMA), "status": "complete",
                    "catalog": {"manifest_sha256": file_sha256(catalog_dir / "manifest.json"),
                                "database_sha256": source["database"]["sha256"]},
                    "database": {"file": "index.sqlite3", "sha256": file_sha256(database_path),
                                 "bytes": database_path.stat().st_size},
                    "policy": {"max_distance": max_distance,
                               "min_similarity_fraction": [similarity.numerator, similarity.denominator],
                               "relative_profile_fractions": [[value.numerator, value.denominator] for value in
                                   sorted({similarity, *(value for value in (Fraction(9, 10), Fraction(19, 20))
                                                         if value >= similarity)})],
                               "max_sequence_length": MAX_LENGTH,
                               "core_predicate": CORE_PREDICATE,
                               "method": "balanced_partition_all_offsets_both_query_orientations"},
                    "counts": counts,
                    "build_limits": {"max_database_bytes": max_database_bytes, "max_nodes": max_nodes,
                                     "deadline_seconds": deadline_seconds, "min_free_bytes": min_free_bytes}}
        if time.monotonic() - started > deadline_seconds:
            raise LimitError("Index build exceeded deadline before publication")
        _check_free(destination, min_free_bytes)
        atomic_json(destination / "manifest.json", manifest)
        return manifest
    except (sqlite3.Error, OSError) as exc:
        raise InputError(f"Index build failed: {type(exc).__name__}") from exc
    finally:
        catalog.close()
        if database is not None:
            database.close()
        if owned and not (destination / "manifest.json").exists():
            shutil.rmtree(destination)


class CandidateIndex:
    """Verified read-only reference/index pair with disk-backed query scratch."""

    def __init__(self, catalog_dir, index_dir, verify_hash=True):
        from .catalog import load_catalog

        self.catalog = None
        self.connection = None
        self.catalog_dir, self.index_dir = Path(catalog_dir).absolute(), Path(index_dir).absolute()
        try:
            self.catalog, self.catalog_manifest = load_catalog(self.catalog_dir, verify_hash=verify_hash)
            self.manifest = read_json(self.index_dir / "manifest.json")
            if not isinstance(self.manifest, dict) or self.manifest.get("schema") != SCHEMA:
                raise InputError("Unsupported candidate-index schema")
            if self.manifest.get("status") != "complete":
                raise InputError("Candidate index is not complete")
            self.max_distance, self.min_similarity = _policy(self.manifest)
            source = self.manifest.get("catalog", {})
            if (source.get("manifest_sha256") != file_sha256(self.catalog_dir / "manifest.json")
                    or source.get("database_sha256") != self.catalog_manifest["database"]["sha256"]):
                raise InputError("Index was built against a different catalog")
            descriptor = self.manifest.get("database", {})
            path = self.index_dir / "index.sqlite3"
            if (descriptor.get("file") != path.name or not path.is_file()
                    or any(p.is_symlink() for p in (path, *path.parents))
                    or descriptor.get("bytes") != path.stat().st_size
                    or not 0 < path.stat().st_size <= MAX_DATABASE_BYTES):
                raise InputError("Invalid candidate-index database file")
            if verify_hash and file_sha256(path) != descriptor.get("sha256"):
                raise InputError("Candidate-index database checksum mismatch")
            self.connection = _readonly(path)
            self.connection.execute("ATTACH DATABASE ? AS catalog",
                                    ((self.catalog_dir / "catalog.sqlite3").as_uri() + "?mode=ro&immutable=1",))
            self.profiles = {}
            for row in self.connection.execute("SELECT length,k,ordinal,start,size FROM profiles ORDER BY length,k,ordinal"):
                length, k, ordinal, start, size = row
                self.profiles.setdefault(length, {}).setdefault(k, []).append((ordinal, start, size))
            for length, profiles in self.profiles.items():
                if not 1 <= length <= MAX_LENGTH or sorted(profiles) != _profiles(
                    length, self.max_distance, self.min_similarity
                ) or any(parts != partition(length, k) for k, parts in profiles.items()):
                    raise InputError("Invalid reference partition inventory")
            lengths = {row[0] for row in self.catalog.execute("SELECT DISTINCT length FROM nodes WHERE eligible=1")}
            if set(self.profiles) != lengths:
                raise InputError("Index length inventory differs from the catalog")
            self.connection.executescript("""
                CREATE TEMP TABLE query_keys(length INTEGER,k INTEGER,ordinal INTEGER,segment TEXT,
                  PRIMARY KEY(length,k,ordinal,segment)) WITHOUT ROWID;
                CREATE TEMP TABLE candidates(node_id INTEGER PRIMARY KEY);
            """)
            self.last_stats = {}
        except (sqlite3.Error, OSError, KeyError, TypeError) as exc:
            self.close()
            raise InputError("Malformed or unreadable candidate-index artifact") from exc
        except Exception:
            self.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if self.catalog is not None:
            self.catalog.close()
            self.catalog = None

    @contextmanager
    def _deadline(self, seconds):
        if type(seconds) not in (int, float) or not math.isfinite(seconds) or seconds <= 0:
            raise InputError("Deadline must be finite and positive")
        until = time.monotonic() + seconds
        self.connection.set_progress_handler(lambda: int(time.monotonic() > until), 1000)
        try:
            yield until
        except sqlite3.OperationalError as exc:
            if time.monotonic() > until:
                raise LimitError("Candidate retrieval exceeded deadline") from exc
            raise InputError("Candidate SQLite operation failed") from exc
        finally:
            self.connection.set_progress_handler(None, 0)

    def iter_candidates(self, sequence, *, max_distance=3, min_similarity=None, core_only=False,
                        min_node_id=None, max_candidates=2_000_000, max_postings=10_000_000,
                        deadline_seconds=300):
        """Yield complete potential neighbors; the caller must verify edit distance.

        max_distance=None selects relative-similarity mode. A limit error means
        the candidate inventory was not completed. Calls on one handle are serial.
        """
        if not isinstance(sequence, str) or not 1 <= len(sequence) <= MAX_LENGTH or set(sequence) - set("ACGT"):
            raise InputError("Indexed query must be uppercase unambiguous DNA of length 1–200")
        integer(max_candidates, "max_candidates", 1, 100_000_000)
        integer(max_postings, "max_postings", 1, 1_000_000_000)
        if type(core_only) is not bool:
            raise InputError("core_only must be boolean")
        if min_node_id is not None:
            integer(min_node_id, "min_node_id", 0, 2**63 - 1)
        if min_similarity is not None:
            if max_distance is not None:
                raise InputError("Choose fixed distance or relative similarity, not both")
            similarity = similarity_fraction(min_similarity)
            if similarity < self.min_similarity:
                raise InputError("Requested similarity is below the index coverage policy")
        else:
            integer(max_distance, "max_distance", 0, self.max_distance)
            similarity = None
        connection = self.connection
        with self._deadline(deadline_seconds) as until:
            connection.execute("DELETE FROM temp.query_keys")
            connection.execute("DELETE FROM temp.candidates")
            fallback_lengths = []
            for length, profiles in self.profiles.items():
                threshold = max_distance if similarity is None else distance_threshold(len(sequence), length, similarity)
                if abs(length - len(sequence)) > threshold:
                    continue
                available = [k for k in profiles if k >= threshold]
                if not available:
                    if threshold < length:
                        raise InputError("Index has no complete partition profile for requested threshold")
                    fallback_lengths.append(length)
                    continue
                k = min(available)
                # If an intact segment shifts by delta and the total length
                # difference is D, its prefix needs >=abs(delta) edits and its
                # suffix >=abs(D-delta). Thus every qualifying alignment has
                # ceil((D-threshold)/2) <= delta <= floor((D+threshold)/2).
                # Use the actual edit threshold, not the stored profile's k.
                difference = len(sequence) - length
                minimum_shift = -((threshold - difference) // 2)
                maximum_shift = (threshold + difference) // 2
                for query in {sequence, reverse_complement(sequence)}:
                    for ordinal, start, size in profiles[k]:
                        for offset in range(max(0, start + minimum_shift),
                                            min(len(query) - size, start + maximum_shift) + 1):
                            connection.execute("INSERT OR IGNORE INTO temp.query_keys VALUES (?,?,?,?)",
                                               (length, k, ordinal, query[offset:offset + size]))
            # Put core/id predicates on the posting key itself. A core query
            # must not visit millions of noncore postings and look up every
            # source node only to discard them. Eligibility and core flags are
            # bound to the immutable catalogue and checked by deep validation.
            posting_restriction = " AND p.is_core=1" if core_only else ""
            fallback_restriction = " AND n.bona_fide_consensus_count>0" if core_only else ""
            parameters = []
            if min_node_id is not None:
                posting_restriction += " AND p.node_id>?"
                fallback_restriction += " AND n.node_id>?"
                parameters.append(min_node_id)
            sql = ("SELECT p.node_id FROM temp.query_keys q CROSS JOIN postings p ON "
                   "p.length=q.length AND p.k=q.k AND p.ordinal=q.ordinal AND p.segment=q.segment "
                   "WHERE 1" + posting_restriction)
            posting_parameters = list(parameters)
            if fallback_lengths:
                placeholders = ",".join("?" for _ in fallback_lengths)
                sql += (" UNION ALL SELECT n.node_id FROM catalog.nodes n WHERE n.eligible=1 AND n.length IN ("
                        + placeholders + ")" + fallback_restriction)
                posting_parameters.extend([*fallback_lengths, *parameters])
            # Preserve duplicate posting accounting, including one posting per
            # fallback node, without a Python callback/INSERT for every posting.
            # LIMIT is only an over-budget sentinel: it can never produce a
            # completed truncated inventory. Successful counts exhaust the join.
            postings = connection.execute("SELECT count(*) FROM (" + sql + " LIMIT ?) AS bounded_postings",
                                          [*posting_parameters, max_postings + 1]).fetchone()[0]
            if postings > max_postings:
                raise LimitError("Candidate retrieval exceeded posting budget")
            if time.monotonic() > until:
                raise LimitError("Candidate retrieval exceeded deadline")
            # Likewise insert at most max_candidates+1 distinct node IDs. This
            # bounds scratch growth before rejection rather than materializing
            # an arbitrarily large set and discovering the overflow afterwards.
            connection.execute("INSERT OR IGNORE INTO temp.candidates SELECT DISTINCT node_id FROM ("
                               + sql + ") AS nominated_nodes LIMIT ?",
                               [*posting_parameters, max_candidates + 1])
            candidates = connection.execute("SELECT count(*) FROM temp.candidates").fetchone()[0]
            if candidates > max_candidates:
                raise LimitError("Candidate retrieval exceeded candidate budget")
            if time.monotonic() > until:
                raise LimitError("Candidate retrieval exceeded deadline")
            self.last_stats = {"postings": postings, "candidates": candidates,
                               "query_keys": connection.execute("SELECT count(*) FROM temp.query_keys").fetchone()[0],
                               "fallback_length_count": len(fallback_lengths)}
            columns = ",".join("n." + name for name in NODE_COLUMNS.split(","))
            # Force candidate-first primary-key lookups. An unconstrained JOIN
            # ordered by n.node_id can scan the entire catalogue for each query
            # merely to avoid sorting, even when only one candidate was found.
            for row in connection.execute("SELECT " + columns + " FROM temp.candidates c CROSS JOIN catalog.nodes n "
                                          "ON n.node_id=c.node_id ORDER BY c.node_id"):
                if time.monotonic() > until:
                    raise LimitError("Candidate retrieval exceeded deadline")
                yield dict(row)


def validate_index(catalog_dir, index_dir, *, deep=True, deadline_seconds=172800) -> dict:
    """Check hashes; optionally prove every posting against every catalog node."""
    with CandidateIndex(catalog_dir, index_dir) as reader:
        with reader._deadline(deadline_seconds) as until:
            actual = reader.connection.execute("SELECT count(*) FROM postings").fetchone()[0]
            if actual != reader.manifest["counts"]["postings"]:
                raise InputError("Candidate posting count mismatch")
            if deep:
                expected_count = nodes = core_nodes = core_postings = 0
                for node_id, sequence, length, core_count in reader.catalog.execute(
                    "SELECT node_id,sequence,length,bona_fide_consensus_count "
                    "FROM nodes WHERE eligible=1 ORDER BY node_id"
                ):
                    if time.monotonic() > until:
                        raise LimitError("Deep index validation exceeded deadline")
                    is_core = int(core_count > 0)
                    expected = {(length, k, ordinal, sequence[start:start + size], is_core)
                                for k, parts in reader.profiles[length].items() for ordinal, start, size in parts}
                    actual_rows = {tuple(row) for row in reader.connection.execute(
                        "SELECT length,k,ordinal,segment,is_core FROM postings WHERE node_id=?", (node_id,))}
                    if expected != actual_rows:
                        raise InputError("Candidate posting inventory differs from source sequences")
                    expected_count += len(expected)
                    nodes += 1
                    core_nodes += is_core
                    core_postings += is_core * len(expected)
                if (expected_count != actual or nodes != reader.manifest["counts"]["nodes"]
                        or core_nodes != reader.manifest["counts"]["core_nodes"]
                        or core_postings != reader.manifest["counts"]["core_postings"]):
                    raise InputError("Candidate index has omitted or unexpected nodes")
        return reader.manifest
