"""Regression controls for bounded search inputs and typed evidence replay."""

import gzip
import json
import sqlite3

import pytest

from crisprspacer import __version__
from crisprspacer.artifacts import validate_result, write_result
from crisprspacer.cli import main
from crisprspacer.database import (_records, analyze_database, build_database, file_digest,
                                   open_database, search_database, validate_database_result)
from crisprspacer.inputs import from_spacers
from crisprspacer.io import InputError, LimitError, digest, read_fasta, write_manifest
from crisprspacer.search import Policy
from crisprspacer.viral_reference import bounded_catalog, open_subjects
from crisprspacer.viral_search import _queries, validate_viral_search
from conftest import record

SEQ = "AACCGGTTAACCGGTTAA"


def rehash(directory, kind):
    (directory / "manifest.json").unlink()
    write_manifest(directory, kind)


def exact_result(tmp_path):
    source = tmp_path / "spacers.fasta"
    source.write_text(f">spacer\n{SEQ}\n")
    database = tmp_path / "database"
    build_database(source, database, "synthetic", "https://example.org/synthetic", "fixture")
    output = tmp_path / "result"
    search_database(database, source, output)
    return database, source, output


def corrupt_gzip(tmp_path):
    raw = gzip.compress(f">query\n{SEQ}\n".encode(), mtime=0)
    path = tmp_path / "corrupt.fasta.gz"
    # Invalid deflate block after a valid gzip header: zlib.error, not BadGzipFile.
    path.write_bytes(raw[:10] + bytes([raw[10] ^ 255]) + raw[11:])
    return path


@pytest.mark.parametrize("reader", [read_fasta, _queries, lambda path: list(_records(path))])
def test_invalid_deflate_payload_is_a_documented_input_error(tmp_path, reader):
    with pytest.raises(InputError, match="[Mm]alformed"):
        reader(corrupt_gzip(tmp_path))


@pytest.mark.parametrize("command", ["search", "database-index", "viral-search"])
def test_corrupt_gzip_cli_does_not_traceback_or_publish(tmp_path, capsys, command):
    source = corrupt_gzip(tmp_path)
    output = tmp_path / "result"
    if command == "search":
        args = ["--spacers", str(source), "--targets", str(source)]
    elif command == "database-index":
        args = ["--spacers", str(source), "--name", "synthetic", "--source-release", "fixture",
                "--source-url", "https://example.org/synthetic"]
    else:
        args = ["--spacers", str(source), "--reference", str(tmp_path / "missing-reference")]
    assert main([command, *args, "--output", str(output)]) == 2
    assert json.loads(capsys.readouterr().err)["status"] == "error"
    assert not output.exists() and not output.with_suffix(".partial").exists()


@pytest.mark.parametrize("field,value", [("sha256", "0" * 64), ("length", True),
                                         ("topology", "circular"), ("sequence", "ACGU")])
def test_exact_database_api_rejects_false_target_identity_or_traversal(tmp_path, field, value):
    database, _, _ = exact_result(tmp_path)
    target = record(SEQ)
    target[field] = value
    with pytest.raises(InputError):
        analyze_database(database, [target])


@pytest.mark.parametrize("targets", [[], [record(SEQ), record(SEQ)]])
def test_exact_database_api_requires_nonempty_unique_target_records(tmp_path, targets):
    database, _, _ = exact_result(tmp_path)
    with pytest.raises(InputError):
        analyze_database(database, targets)


@pytest.mark.parametrize("section,field,value", [
    ("counts", "matches", True), ("counts", "ambiguous_seed_windows", False),
    ("policy", "substitutions", False), (None, "search_complete", 1),
])
def test_exact_replay_rejects_booleans_as_numbers(tmp_path, section, field, value):
    database, _, output = exact_result(tmp_path)
    report = json.loads((output / "analysis.json").read_text())
    (report if section is None else report[section])[field] = value
    (output / "analysis.json").write_text(json.dumps(report))
    rehash(output, "crisprspacer.database-search")
    with pytest.raises(InputError, match="replay"):
        validate_database_result(output, database)


@pytest.mark.parametrize("value", ["x" * 20_000, b"{}"])
def test_metadata_is_size_and_type_checked_before_materializing(tmp_path, value):
    database, _, _ = exact_result(tmp_path)
    with sqlite3.connect(database / "spacers.sqlite3") as connection:
        connection.execute("UPDATE metadata SET value=?", (value,))
    reference = json.loads((database / "reference.json").read_text())
    reference["database_sha256"] = file_digest(database / "spacers.sqlite3")
    (database / "reference.json").write_text(json.dumps(reference))
    with pytest.raises(InputError, match="metadata"):
        open_database(database)


def test_metadata_view_cannot_execute_an_unbounded_query_before_search_deadline(tmp_path):
    database, _, _ = exact_result(tmp_path)
    with sqlite3.connect(database / "spacers.sqlite3") as connection:
        connection.execute("DROP TABLE metadata")
        connection.execute("""CREATE VIEW metadata AS
            WITH RECURSIVE slow(x) AS (SELECT 1 UNION ALL SELECT x+1 FROM slow WHERE x<100000000)
            SELECT sum(x) AS value FROM slow""")
    reference = json.loads((database / "reference.json").read_text())
    reference["database_sha256"] = file_digest(database / "spacers.sqlite3")
    (database / "reference.json").write_text(json.dumps(reference))
    with pytest.raises(InputError, match="regular table"):
        analyze_database(database, [record(SEQ)], timeout=1)


def test_prior_exact_and_core_reports_remain_replayable(tmp_path):
    database, _, output = exact_result(tmp_path)
    report = json.loads((output / "analysis.json").read_text())
    assert report["tool"]["version"] == __version__
    for version in ("0.2.0", "0.3.0"):
        report["tool"]["version"] = version
        (output / "analysis.json").write_text(json.dumps(report))
        rehash(output, "crisprspacer.database-search")
        assert validate_database_result(output, database)["status"] == "valid"
    core = tmp_path / "core"
    inventory = from_spacers([record(SEQ)], digest(b"fixture"))
    report = write_result(core, inventory, [record(SEQ)], Policy())
    report["tool_version"] = "0.3.0"
    (core / "analysis.json").write_text(json.dumps(report))
    rehash(core, "crisprspacer.result")
    assert validate_result(core)["status"] == "valid"


@pytest.mark.parametrize("error", [LimitError("catalog deadline"), KeyboardInterrupt()])
def test_catalog_scan_preserves_checkpoint_failure_and_cancellation(tmp_path, error):
    directory = tmp_path / "catalog"
    directory.mkdir()
    with sqlite3.connect(directory / "subjects.sqlite3") as connection:
        connection.execute("CREATE TABLE subjects(id TEXT PRIMARY KEY, accession TEXT, title TEXT, "
                           "sequence TEXT, length INTEGER, sha256 TEXT)")
        connection.executemany("INSERT INTO subjects VALUES (?,?,?,?,?,?)", [
            (f"v{i:09d}", f"SYNTH{i}.1", "generated", SEQ, len(SEQ), digest(SEQ)) for i in range(1, 1001)])
    calls = 0

    def checkpoint():
        nonlocal calls
        calls += 1
        if calls == 2:
            raise error

    with pytest.raises(type(error)) as raised:
        open_subjects(directory, checkpoint)
    assert raised.value is error and calls == 2


def test_catalog_evidence_queries_are_also_interruptible():
    connection = sqlite3.connect(":memory:")
    calls = 0

    def checkpoint():
        nonlocal calls
        calls += 1
        if calls == 2:
            raise LimitError("evidence deadline")

    try:
        with pytest.raises(LimitError, match="evidence deadline"), bounded_catalog(connection, checkpoint):
            connection.execute("""WITH RECURSIVE slow(x) AS
                (SELECT 1 UNION ALL SELECT x+1 FROM slow WHERE x<100000000) SELECT sum(x) FROM slow""").fetchone()
        # Handler removal also leaves the connection usable after an interruption.
        assert connection.execute("SELECT 1").fetchone() == (1,)
    finally:
        connection.close()


def test_prior_viral_report_version_is_accepted(tmp_path):
    # Reuse a real native output through the existing module's fixture-building helpers.
    from test_viral_search import SEQ as viral_sequence, spacers
    from crisprspacer.viral_reference import build_viral_reference
    from crisprspacer.viral_search import search_viruses
    import shutil
    if not all(shutil.which(name) for name in ("blastn", "makeblastdb", "blastdbcmd")):
        pytest.skip("Local BLAST+ required")
    source = spacers(tmp_path, [("query", viral_sequence)])
    reference = tmp_path / "reference"
    build_viral_reference([source], reference, name="synthetic", source_release="fixture",
                          source_url="https://example.org/synthetic")
    output = tmp_path / "result"
    report = search_viruses(source, reference, output)
    assert report["tool_version"] == __version__
    report["tool_version"] = "0.3.0"
    (output / "result.json").write_text(json.dumps(report))
    rehash(output, "crisprspacer.viral-search")
    assert validate_viral_search(output, reference) == report


@pytest.mark.parametrize("section,field,value", [
    ("counts", "query_occurrences", True), ("queries", "eligible", 1),
    ("candidate_viruses", "query_occurrence_count", True),
])
def test_viral_replay_rejects_booleans_as_numbers(tmp_path, section, field, value):
    from test_viral_search import SEQ as viral_sequence, spacers
    from crisprspacer.viral_reference import build_viral_reference
    from crisprspacer.viral_search import search_viruses
    import shutil
    if not all(shutil.which(name) for name in ("blastn", "makeblastdb", "blastdbcmd")):
        pytest.skip("Local BLAST+ required")
    source = spacers(tmp_path, [("query", viral_sequence)])
    reference = tmp_path / "reference"
    build_viral_reference([source], reference, name="synthetic", source_release="fixture",
                          source_url="https://example.org/synthetic")
    output = tmp_path / "result"
    report = search_viruses(source, reference, output)
    row = report[section] if section == "counts" else report[section][0]
    row[field] = value
    (output / "result.json").write_text(json.dumps(report))
    rehash(output, "crisprspacer.viral-search")
    with pytest.raises(InputError):
        validate_viral_search(output, reference)
