from __future__ import annotations

import json
import os
import random
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from spacepharer.engine import EngineError, run
from spacepharer.io import InputError, digest, fasta, fasta_bytes, groups, json_value, read, transaction
from spacepharer.pipeline import associate
from spacepharer.reference import load, prepare
from spacepharer.report import oriented, parse

ROOT = Path(__file__).resolve().parents[1]
NATIVE = os.environ.get("SPACEPHARER_NATIVE")
ORACLE = os.environ.get("SPACEPHARER_ORACLE")
BASELINE = os.environ.get("SPACEPHARER_CALIBRATION_BASELINE")


def child_is_running(pid: int) -> bool:
    """A terminated child may disappear between /proc lookup and reading."""
    try:
        return Path(f"/proc/{pid}/stat").read_text().split()[2] != "Z"
    except (FileNotFoundError, ProcessLookupError):
        return False


@pytest.mark.parametrize("error", [FileNotFoundError, ProcessLookupError])
def test_child_is_running_after_disappearance(monkeypatch, error):
    def missing(path):
        assert path == Path("/proc/123/stat")
        raise error()

    monkeypatch.setattr(Path, "read_text", missing)
    assert not child_is_running(123)


def test_child_is_running_propagates_permission_error(monkeypatch):
    def denied(path):
        assert path == Path("/proc/123/stat")
        raise PermissionError()

    monkeypatch.setattr(Path, "read_text", denied)
    with pytest.raises(PermissionError):
        child_is_running(123)


def test_json_duplicate_and_nonfinite():
    for data in (b'{"x":1,"x":2}', b'{"x":NaN}', b"[]"):
        with pytest.raises(InputError):
            json_value(data)


@pytest.mark.parametrize(
    "text", [">x\nACGT\n>x\nACGT\n", "AAAA\n", ">bad/id\nACGT\n", ">x\n", ">x\n" + "A" * 201]
)
def test_invalid_fasta(tmp_path, text):
    path = tmp_path / "input.fasta"
    path.write_text(text)
    with pytest.raises(InputError):
        fasta(path)


def test_group_partition(tmp_path):
    path = tmp_path / "groups.json"
    records = [{"id": "a", "sequence": "ACGT"}, {"id": "b", "sequence": "GGGG"}]
    for members in (["a"], ["a", "a", "b"], ["a", "b", "c"]):
        path.write_text(
            json.dumps(
                {"schema_version": "spacepharer.groups/1.0.0", "groups": [{"id": "host", "members": members}]}
            )
        )
        with pytest.raises(InputError):
            groups(path, records)


def test_symlink_ancestor_and_atomic_failure(tmp_path):
    directory = tmp_path / "dir"
    directory.mkdir()
    (directory / "file").write_text("hi")
    (tmp_path / "link").symlink_to(directory, target_is_directory=True)
    with pytest.raises(InputError):
        read(tmp_path / "link/file")
    output = tmp_path / "output"
    with pytest.raises(RuntimeError), transaction(output) as staging:
        (staging / "file").write_text("partial")
        raise RuntimeError("failure")
    assert not output.exists()
    assert not list(tmp_path.glob(".spacepharer-*"))


def test_reverse_source_binding():
    assert oriented("AACCGGTT", 6, 3) == "CCGG"
    query = {"q": {"sequence": "AAGCCC", "occurrence_ids": ["s"], "group_id": "g"}}
    text = "#group\tt\t10\t-1\t1\n>q\tt\t0.001\t6\t1\t1\t6\nGGGCTT\nGGGCTT\n"
    result = parse(text, query, {"t": "GGGCTT"}, {"group": "g"}, False)
    assert result[0]["native_fdr"] is None
    assert result[0]["status"] == "uncalibrated_candidate"
    with pytest.raises(InputError):
        parse(text.replace("\t6\t1\t1\t6", "\t6\t2\t1\t6"), query, {"t": "GGGCTT"}, {"group": "g"}, False)
    with pytest.raises(InputError):
        parse(text, query, {"t": "GGGCTT"}, {"group": "unrelated"}, False)


@pytest.mark.parametrize("mutation", ["association", "support", "double_gap", "count"])
def test_native_report_rejects_inflated_evidence(mutation):
    query = {"q": {"sequence": "AAAA", "occurrence_ids": ["s"], "group_id": "g"}}
    header = "#group\tt\t10\t-1\t1\n"
    hit = ">q\tt\t0.001\t1\t4\t1\t4\nAAAA\nAAAA\n"
    text = header + hit
    if mutation == "association":
        text += header + hit
    elif mutation == "support":
        text = header.replace("\t1\n", "\t2\n") + hit + hit
    elif mutation == "double_gap":
        text = text.replace("AAAA\n", "AA-AA\n")
    else:
        text = header.replace("\t1\n", "\t10001\n") + hit
    with pytest.raises(InputError):
        parse(text, query, {"t": "AAAA"}, {"group": "g"}, False)


def test_distinct_target_intervals_remain_independent_support():
    query = {"q": {"sequence": "AAAA", "occurrence_ids": ["s", "copy"], "group_id": "g"}}
    text = (
        "#group\tt\t10\t-1\t2\n>q\tt\t0.001\t1\t4\t1\t4\nAAAA\nAAAA\n>q\tt\t0.002\t1\t4\t5\t8\nAAAA\nAAAA\n"
    )
    rows = parse(text, query, {"t": "AAAAAAAA"}, {"group": "g"}, False)
    assert [hit["occurrence_ids"] for hit in rows[0]["hits"]] == [["s", "copy"], ["s", "copy"]]


def test_multiline_genome_fasta_preserves_records_and_limits():
    sequence = "ACGTN" * 20000
    data = (">first\n" + "\n".join(sequence) + "\n>second\nacgt\n").encode()
    assert fasta_bytes(data, maximum=2, length=len(sequence)) == [
        {"id": "first", "sequence": sequence},
        {"id": "second", "sequence": "ACGT"},
    ]
    with pytest.raises(InputError, match="length limit"):
        fasta_bytes(data, maximum=2, length=len(sequence) - 1)


def database(path: Path, scores: list[float]):
    data = "".join(f"{i}\t{score:.8f}\t0.75\n" for i, score in enumerate(scores)).encode() + b"\0"
    path.write_bytes(data)
    path.with_suffix(".index").write_text(f"0\t0\t{len(data)}\n")
    path.with_suffix(".dbtype").write_bytes((5).to_bytes(4, "little"))


@pytest.mark.skipif(not NATIVE or not ORACLE, reason="native and original executables required")
def test_calibrated_original_native_parity(tmp_path):
    rng = random.Random(491)
    database(
        tmp_path / "positive",
        [rng.uniform(10, 30) for _ in range(200)] + [rng.expovariate(1) for _ in range(100)],
    )
    database(tmp_path / "negative", [rng.expovariate(1) for _ in range(300)])
    for name, executable in (("original", ORACLE), ("modern", NATIVE)):
        subprocess.run(
            [
                executable,
                "filtermatchbyfdr",
                str(tmp_path / "positive"),
                str(tmp_path / "negative"),
                str(tmp_path / name),
                "--report-fdr",
                "1",
                "--threads",
                "1",
                "-v",
                "0",
            ],
            check=True,
        )
    assert (tmp_path / "original").read_bytes() == (tmp_path / "modern").read_bytes()
    receipt = json.loads((tmp_path / "modern.calibration.json").read_text())
    assert receipt["status"] == "estimated"
    assert receipt["control_score_count"] == 300
    assert 0 < len((tmp_path / "modern").read_text().splitlines()) < 300


@pytest.mark.skipif(not NATIVE or not BASELINE, reason="native and prepatch calibration oracle required")
@pytest.mark.parametrize("population", ["calibrated", "ties", "unavailable"])
def test_calibration_lookup_exact_prepatch_parity(tmp_path, population):
    rng = random.Random(491)
    targets = [rng.uniform(10, 30) for _ in range(200)] + [rng.expovariate(1) for _ in range(100)]
    controls = [rng.expovariate(1) for _ in range(300)]
    if population == "ties":
        targets = [round(score, 1) for score in targets] * 3
        controls = [round(score, 1) for score in controls] * 3
    elif population == "unavailable":
        controls = []
    database(tmp_path / "positive", targets)
    database(tmp_path / "negative", controls)
    for name, executable in (("baseline", BASELINE), ("modern", NATIVE)):
        subprocess.run(
            [
                executable,
                "filtermatchbyfdr",
                str(tmp_path / "positive"),
                str(tmp_path / "negative"),
                str(tmp_path / name),
                "--report-fdr",
                "1",
                "--threads",
                "1",
                "-v",
                "0",
            ],
            check=True,
            timeout=30,
        )
    for suffix in ("", ".index", ".dbtype", ".calibration.json"):
        assert (tmp_path / ("baseline" + suffix)).read_bytes() == (
            tmp_path / ("modern" + suffix)
        ).read_bytes()
    receipt = json.loads((tmp_path / "modern.calibration.json").read_text())
    assert receipt["status"] == ("unavailable" if population == "unavailable" else "estimated")


def test_reference_accepts_existing_210_indexes(reference, tmp_path):
    import shutil

    panel = tmp_path / "legacy-panel"
    shutil.copytree(reference, panel)
    path = panel / "reference.json"
    manifest = json.loads(path.read_text())
    manifest["engine_version"] = "2.1.0"
    path.write_text(json.dumps(manifest))
    identity, sequences = load(panel)
    try:
        assert identity["engine_version"] == "2.1.0"
        assert len(sequences) == identity["genome_count"]
    finally:
        sequences.close()


@pytest.mark.skipif(not NATIVE, reason="native executable required")
def test_short_calibration_is_explicit_at_quiet_verbosity(tmp_path):
    database(tmp_path / "positive", [10, 5])
    database(tmp_path / "negative", [1])
    subprocess.run(
        [
            NATIVE,
            "filtermatchbyfdr",
            str(tmp_path / "positive"),
            str(tmp_path / "negative"),
            str(tmp_path / "result"),
            "--report-fdr",
            "1",
            "--threads",
            "1",
            "-v",
            "0",
        ],
        check=True,
    )
    receipt = json.loads((tmp_path / "result.calibration.json").read_text())
    assert receipt["status"] == "unavailable"
    assert all(
        float(row.split()[3]) == -1 for row in (tmp_path / "result").read_text().strip("\0\n").splitlines()
    )


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    if not NATIVE:
        pytest.skip("native executable required")
    import gzip

    directory = tmp_path_factory.mktemp("reference")
    files = []
    for path in sorted((ROOT / "examples").glob("*.gz")):
        target = directory / path.stem
        target.write_bytes(gzip.decompress(path.read_bytes()))
        files.append(target)
    prepare(files, directory / "reference", "upstream-three-phages")
    return directory / "reference"


def test_reference_hash_tamper(reference, tmp_path):
    import shutil

    copy = tmp_path / "reference"
    shutil.copytree(reference, copy)
    (copy / "db/target.index").write_text("0\t0\t1\n")
    with pytest.raises(InputError, match="checksum"):
        load(copy)


def test_real_biological_pipeline_original_parity(reference, tmp_path):
    fixture = ROOT / "tests/fixtures"
    report = associate(fixture / "spacers.fasta", fixture / "groups.json", reference, tmp_path / "report")
    assert report["calibration_status"] == "unavailable"
    assert [(x["target_id"], x["combined_score"], x["native_hit_count"]) for x in report["associations"]] == [
        ("AY369265.2", 12.24, 1),
        ("AJ242593.1", 74.56, 5),
    ]
    assert all(
        x["native_fdr"] is None and x["status"] == "uncalibrated_candidate" for x in report["associations"]
    )
    original_lines = (fixture / "original-results.tsv").read_text().splitlines()
    original_hits = []
    for i, line in enumerate(original_lines):
        if line.startswith(">"):
            parts = line.split("\t")
            original_hits.append(
                (parts[0][1:], float(parts[2]), list(map(int, parts[3:])), original_lines[i + 1 : i + 3])
            )
    modern_hits = [
        (
            hit["occurrence_ids"][0],
            hit["p_best_hit"],
            [hit[k] for k in ("query_start", "query_end", "target_start", "target_end")],
            [hit["query_alignment"], hit["target_alignment"]],
        )
        for row in report["associations"]
        for hit in row["hits"]
    ]
    assert sorted(modern_hits) == sorted(original_hits)


def test_unsupported_and_duplicate_occurrences(reference, tmp_path):
    sequence = "CGGCATACCAAGTTTTTTCGGGCCTTTGGGAAAGAA"
    fasta_path = tmp_path / "input.fasta"
    fasta_path.write_text(f">a\n{sequence}\n>copy\n{sequence}\n>ambiguous\nACGTNN\n")
    group_path = tmp_path / "groups.json"
    group_path.write_text(
        json.dumps(
            {
                "schema_version": "spacepharer.groups/1.0.0",
                "groups": [{"id": "sample", "members": ["a", "copy", "ambiguous"]}],
            }
        )
    )
    report = associate(fasta_path, group_path, reference, tmp_path / "report")
    assert report["groups"][0]["unique_sequence_count"] == 1
    assert report["groups"][0]["unsupported_occurrence_ids"] == ["ambiguous"]
    assert len(report["records"]) == 3


def test_missing_native_version(monkeypatch, tmp_path):
    executable = tmp_path / "wrong"
    executable.write_text("#!/bin/sh\necho 1.0.0\n")
    executable.chmod(0o755)
    monkeypatch.setenv("SPACEPHARER_NATIVE", str(executable))
    with pytest.raises(EngineError, match="version"):
        run([], tmp_path / "log")


def test_native_descendant_cleanup_after_leader_exit(monkeypatch, tmp_path):
    executable = tmp_path / "native"
    pidfile = tmp_path / "child.pid"
    executable.write_text(
        f"#!{sys.executable}\nimport sys,subprocess\nif sys.argv[1:] == ['version']:\n print('2.1.1')\nelse:\n p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'])\n open({str(pidfile)!r},'w').write(str(p.pid))\n"
    )
    executable.chmod(0o755)
    monkeypatch.setenv("SPACEPHARER_NATIVE", str(executable))
    run(["test"], tmp_path / "log")
    pid = int(pidfile.read_text())
    for _ in range(30):
        if not child_is_running(pid):
            break
        time.sleep(0.01)
    else:
        os.kill(pid, signal.SIGKILL)
        pytest.fail("native descendant survived completed leader")


def test_missing_calibration_does_not_imply_estimation(reference, tmp_path, monkeypatch):
    from spacepharer import pipeline

    input_path = tmp_path / "input.fasta"
    input_path.write_text(">a\nACGTTTCATCGGATCGTGGCTAACGTTAGC\n")
    group_path = tmp_path / "groups.json"
    group_path.write_text(
        json.dumps({"schema_version": "spacepharer.groups/1.0.0", "groups": [{"id": "g", "members": ["a"]}]})
    )
    monkeypatch.setattr(pipeline, "run", lambda *args, **kwargs: None)
    with pytest.raises(InputError, match="receipt missing"):
        pipeline.associate(input_path, group_path, reference, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("target_count", [0, 1])
def test_estimated_calibration_requires_target_population(reference, tmp_path, monkeypatch, target_count):
    from spacepharer import pipeline

    input_path = tmp_path / "input.fasta"
    input_path.write_text(">a\nACGTTTCATCGGATCGTGGCTAACGTTAGC\n")
    group_path = tmp_path / "groups.json"
    group_path.write_text(
        json.dumps({"schema_version": "spacepharer.groups/1.0.0", "groups": [{"id": "g", "members": ["a"]}]})
    )

    def fake_run(args, log):
        if args[0] == "predictmatch":
            receipt = Path(args[5]) / "receipt"
            receipt.mkdir(parents=True)
            (receipt / "match.calibration.json").write_text(
                json.dumps(
                    {
                        "schema_version": "spacepharer.calibration/1.0.0",
                        "status": "estimated",
                        "target_score_count": target_count,
                        "control_score_count": 300,
                    }
                )
            )

    monkeypatch.setattr(pipeline, "run", fake_run)
    with pytest.raises(InputError, match="calibration support"):
        pipeline.associate(input_path, group_path, reference, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_report_hashes_bind_parsed_snapshots_and_close_catalogue(reference, tmp_path, monkeypatch):
    from spacepharer import pipeline

    input_path = tmp_path / "input.fasta"
    input_raw = b">a\nACGTTTCATCGGATCGTGGCTAACGTTAGC\n"
    input_path.write_bytes(input_raw)
    group_path = tmp_path / "groups.json"
    group_raw = json.dumps(
        {"schema_version": "spacepharer.groups/1.0.0", "groups": [{"id": "g", "members": ["a"]}]}
    ).encode()
    group_path.write_bytes(group_raw)
    manifest, targets = load(reference)
    monkeypatch.setattr(pipeline, "load", lambda path: (manifest, targets))

    def fake_run(args, log):
        input_path.write_text(">changed\nAAAA\n")
        group_path.write_text("{}")
        if args[0] == "predictmatch":
            receipt = Path(args[5]) / "receipt"
            receipt.mkdir(parents=True)
            (receipt / "match.calibration.json").write_text(
                json.dumps(
                    {
                        "schema_version": "spacepharer.calibration/1.0.0",
                        "status": "unavailable",
                        "target_score_count": 0,
                        "control_score_count": 0,
                    }
                )
            )
            Path(args[4]).write_text("")

    monkeypatch.setattr(pipeline, "run", fake_run)
    report = pipeline.associate(input_path, group_path, reference, tmp_path / "out")
    assert report["input_sha256"] == digest(input_raw)
    assert report["groups_sha256"] == digest(group_raw)
    assert report["records"] == [{"id": "a", "sequence": input_raw.decode().splitlines()[1]}]
    import sqlite3

    with pytest.raises(sqlite3.ProgrammingError, match="closed database"):
        targets.connection.execute("SELECT 1")


@pytest.mark.parametrize("options", [{"threads": True}, {"threads": 1.5}, {"memory_mb": 2048.5}])
def test_resource_limits_require_integers(tmp_path, options):
    with pytest.raises(InputError, match="resource limits"):
        associate(
            tmp_path / "input", tmp_path / "groups", tmp_path / "reference", tmp_path / "out", **options
        )


def test_all_unsupported_never_becomes_negative(reference, tmp_path):
    input_path = tmp_path / "input.fasta"
    input_path.write_text(">a\nACGTNNNN\n>b\nACGT\n")
    group_path = tmp_path / "groups.json"
    group_path.write_text(
        json.dumps(
            {"schema_version": "spacepharer.groups/1.0.0", "groups": [{"id": "g", "members": ["a", "b"]}]}
        )
    )
    report = associate(input_path, group_path, reference, tmp_path / "out")
    assert report["calibration_status"] == "not_run"
    assert report["calibration"] is None
    assert report["groups"][0]["unsupported_occurrence_ids"] == ["a", "b"]
    assert report["associations"] == []


@pytest.mark.parametrize("interrupt", ["timeout", "sigterm"])
def test_native_group_killed_on_interruption(tmp_path, interrupt):
    executable = tmp_path / "native"
    pidfile = tmp_path / "child.pid"
    executable.write_text(
        f"#!{sys.executable}\nimport sys,subprocess,time\nif sys.argv[1:] == ['version']:\n print('2.1.1')\nelse:\n p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'])\n open({str(pidfile)!r},'w').write(str(p.pid))\n time.sleep(60)\n"
    )
    executable.chmod(0o755)
    code = f"from spacepharer.engine import run; from pathlib import Path; run(['test'],Path({str(tmp_path / 'log')!r}),timeout=1)"
    runner = subprocess.Popen(
        [sys.executable, "-c", code],
        env={**os.environ, "SPACEPHARER_NATIVE": str(executable)},
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    for _ in range(100):
        if pidfile.exists():
            break
        time.sleep(0.01)
    assert pidfile.exists()
    if interrupt == "sigterm":
        runner.terminate()
    assert runner.wait(timeout=5) != 0
    child = int(pidfile.read_text())
    for _ in range(100):
        if not child_is_running(child):
            break
        time.sleep(0.01)
    else:
        os.kill(child, signal.SIGKILL)
        pytest.fail("native child survived interrupted runner")


def test_native_threshold_boundary_remains_candidate():
    query = {"q": {"sequence": "AAAA", "occurrence_ids": ["s"], "group_id": "g"}}
    text = "#group\tt\t10\t0.05586\t1\n>q\tt\t0.001\t1\t4\t1\t4\nAAAA\nAAAA\n"
    rows = parse(text, query, {"t": "AAAA"}, {"group": "g"}, True, 0.05)
    assert rows[0]["status"] == "above_cutoff_candidate"
    assert rows[0]["native_fdr"] == 0.05586
    with pytest.raises(InputError, match="calibration receipt"):
        parse(text, query, {"t": "AAAA"}, {"group": "g"}, False)


def test_source_archive_contains_required_upstream_files():
    # Applies both in checkouts and unpacked source archives. Scientific models
    # and vendored build recipes must never disappear behind ignore rules.
    inventory = json.loads((ROOT / "docs/upstream-regular-files.json").read_text())
    assert len(inventory) > 1200
    for name in inventory:
        assert (ROOT / name).is_file(), name
    assert (ROOT / "lib/mmseqs/data/VTML40.out").stat().st_size > 0
    assert (ROOT / "lib/mmseqs/lib/zstd/build/cmake/lib/CMakeLists.txt").is_file()


@pytest.mark.parametrize("mutation", ["extra", "missing", "symlink"])
def test_reference_inventory_boundary(reference, tmp_path, mutation):
    import shutil

    panel = tmp_path / "panel"
    shutil.copytree(reference, panel)
    if mutation == "extra":
        (panel / "db/unlisted").write_text("unexpected")
    elif mutation == "missing":
        document = json.loads((panel / "reference.json").read_text())
        del document["files"]["db/target.index"]
        (panel / "db/target.index").unlink()
        (panel / "reference.json").write_text(json.dumps(document))
    else:
        (panel / "db").rename(tmp_path / "indexes")
        (panel / "db").symlink_to(tmp_path / "indexes", target_is_directory=True)
    with pytest.raises(InputError):
        load(panel)


@pytest.mark.skipif(not NATIVE, reason="native executable required")
def test_createdb_preserves_more_than_65535_source_groups(tmp_path):
    """The former uint16 source index silently merged genome 65536 into 0."""
    import io
    import tarfile

    archive_path = tmp_path / "genomes.tar"
    count = 65_537
    with tarfile.open(archive_path, "w", format=tarfile.USTAR_FORMAT) as archive:
        for index in range(count):
            data = f">seq{index}\nACGTACGTACGTACGTACGTACGTACGTACGT\n".encode()
            member = tarfile.TarInfo(f"source{index}.fasta")
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
    for args in (
        ["tar2db", str(archive_path), str(tmp_path / "tardb"), "--threads", "2", "-v", "0"],
        ["createdb", str(tmp_path / "tardb"), str(tmp_path / "seqdb"), "-v", "0"],
    ):
        subprocess.run([NATIVE, *args], check=True, timeout=120, stdout=subprocess.DEVNULL)
    sources = dict(line.split("\t") for line in (tmp_path / "seqdb.source").read_text().splitlines())
    rows = [line.split("\t") for line in (tmp_path / "seqdb.lookup").read_text().splitlines()]
    assert len(rows) == len(sources) == count
    for _, name, source in rows:
        assert sources[source] == f"source{name[3:]}.fasta"
    assert len({row[2] for row in rows}) == count


@pytest.mark.skipif(not NATIVE, reason="native executable required")
def test_broad_reference_over_100_and_lazy_sequence_binding(tmp_path):
    from spacepharer.reference import verify

    rng = random.Random(613)
    files = []
    for index in range(101):
        path = tmp_path / f"genome{index}.fasta"
        sequence = "ATG" + "".join(rng.choices("ACGT", k=297)) + "NNN"
        path.write_text(f">target{index}\n{sequence}\n")
        files.append(path)
    panel = tmp_path / "panel"
    prepare(files, panel, "broad-test", scope="broad_phage_reference")
    manifest, sequences = load(panel)
    assert manifest["genome_count"] == 101
    assert manifest["total_bases"] == 101 * 303
    assert "genomes" not in manifest
    assert not sequences.cache
    assert "target100" in sequences
    assert not sequences.cache
    assert sequences["target100"] == files[100].read_text().splitlines()[1]
    for index in range(12):
        sequences[f"target{index}"]
    assert len(sequences.cache) == 8
    assert verify(panel)["genome_count"] == 101
    sequences.close()
    # A genome is loaded and hashed only when needed. Operator verification must
    # additionally catch corrupt genomes that no particular submitted job hits.
    target = panel / "genomes/target000100.fasta"
    target.write_text(target.read_text().replace("ATG", "ACG", 1))
    _, lazy = load(panel)
    assert "target100" in lazy
    with pytest.raises(InputError, match="binding mismatch"):
        lazy["target100"]
    with pytest.raises(InputError, match="binding mismatch"):
        verify(panel)
    lazy.close()


def test_operator_genome_list_and_metadata_limits(tmp_path):
    from spacepharer.reference import genome_list, source_metadata

    path = tmp_path / "paths.txt"
    path.write_text("first.fasta\nsecond.fasta\n")
    assert genome_list(path) == [tmp_path / "first.fasta", tmp_path / "second.fasta"]
    path.write_text("first.fasta\n\nsecond.fasta\n")
    with pytest.raises(InputError):
        genome_list(path)
    with pytest.raises(InputError):
        source_metadata({"collection": "missing required provenance"})
    value = source_metadata(None)
    value["selection_policy"] = "control\ncharacter"
    with pytest.raises(InputError):
        source_metadata(value)


@pytest.mark.parametrize("mutation", ["catalogue", "count", "path", "view", "unindexed"])
def test_catalogue_integrity_and_schema(reference, tmp_path, mutation):
    import hashlib
    import shutil
    import sqlite3

    panel = tmp_path / "panel"
    shutil.copytree(reference, panel)
    manifest_path = panel / "reference.json"
    manifest = json.loads(manifest_path.read_text())
    catalogue = panel / "catalogue.sqlite3"
    if mutation == "count":
        manifest["genome_count"] += 1
    elif mutation == "catalogue":
        with catalogue.open("ab") as handle:
            handle.write(b"tampered")
    else:
        with sqlite3.connect(catalogue) as connection:
            if mutation == "path":
                connection.execute("UPDATE genomes SET file='../outside.fasta' WHERE rowid=1")
            elif mutation == "unindexed":
                connection.execute("ALTER TABLE genomes RENAME TO original")
                connection.execute("CREATE TABLE genomes(id TEXT, length INTEGER, sha256 TEXT, file TEXT)")
                connection.execute("INSERT INTO genomes SELECT * FROM original")
                connection.execute("DROP TABLE original")
            else:
                connection.execute("CREATE VIEW malicious AS SELECT 1")
        manifest["catalogue"]["sha256"] = hashlib.sha256(catalogue.read_bytes()).hexdigest()
        manifest["catalogue"]["size_bytes"] = catalogue.stat().st_size
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(InputError):
        load(panel)


@pytest.mark.parametrize("score_count", [1_000_001, 10_000_000, 10_000_001])
def test_broad_reference_calibration_population_bound(reference, tmp_path, monkeypatch, score_count):
    from spacepharer import pipeline

    query = tmp_path / "query.fasta"
    query.write_text(">s\nACGTTTCATCGGATCGTGGCTAACGTTAGC\n")
    assignment = tmp_path / "groups.json"
    assignment.write_text(
        json.dumps({"schema_version": "spacepharer.groups/1.0.0", "groups": [{"id": "g", "members": ["s"]}]})
    )

    def native_receipt(args, *unused, **kwargs):
        if args[0] == "predictmatch":
            directory = Path(args[5]) / "native-run"
            directory.mkdir(parents=True)
            (directory / "match.calibration.json").write_text(
                json.dumps(
                    {
                        "schema_version": "spacepharer.calibration/1.0.0",
                        "status": "estimated",
                        "target_score_count": score_count,
                        "control_score_count": score_count,
                    }
                )
            )
            Path(args[4]).write_text("")

    monkeypatch.setattr(pipeline, "run", native_receipt)
    if score_count <= 10_000_000:
        report = pipeline.associate(query, assignment, reference, tmp_path / "output")
        assert report["calibration"]["target_score_count"] == score_count
        assert report["calibration_status"] == "estimated"
    else:
        with pytest.raises(InputError, match="calibration count"):
            pipeline.associate(query, assignment, reference, tmp_path / "output")
        assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("version,exit_code", [("2.1.1", 0), ("1.0.0", 0), ("2.1.1", 1)])
@pytest.mark.parametrize("worker_thread", [False, True])
def test_version_probe_cleans_descendants(monkeypatch, tmp_path, version, exit_code, worker_thread):
    from spacepharer.engine import executable

    native = tmp_path / "native"
    pidfile = tmp_path / "version-child.pid"
    native.write_text(
        f"#!{sys.executable}\nimport subprocess,sys\n"
        "p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'],"
        "stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)\n"
        f"open({str(pidfile)!r},'w').write(str(p.pid))\n"
        f"print({version!r})\nsys.exit({exit_code})\n"
    )
    native.chmod(0o755)
    monkeypatch.setenv("SPACEPHARER_NATIVE", str(native))
    handlers = [signal.getsignal(signum) for signum in (signal.SIGTERM, signal.SIGINT)]

    def probe():
        if worker_thread:
            with ThreadPoolExecutor(max_workers=1) as pool:
                return pool.submit(executable).result(timeout=15)
        return executable()

    if version == "2.1.1" and exit_code == 0:
        assert probe() == str(native.resolve())
    else:
        with pytest.raises(EngineError):
            probe()
    assert handlers == [signal.getsignal(signum) for signum in (signal.SIGTERM, signal.SIGINT)]
    child = int(pidfile.read_text())
    for _ in range(30):
        if not child_is_running(child):
            break
        time.sleep(0.01)
    else:
        os.kill(child, signal.SIGKILL)
        pytest.fail("native version descendant survived identity verification")


@pytest.mark.parametrize("interrupt", ["timeout", "sigterm"])
def test_version_probe_killed_on_interruption(tmp_path, interrupt):
    native = tmp_path / "native"
    pidfile = tmp_path / "child.pid"
    native.write_text(
        f"#!{sys.executable}\nimport subprocess,sys,time\n"
        "p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'],"
        "stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)\n"
        f"open({str(pidfile)!r},'w').write(str(p.pid))\n"
        "print('2.1.1',flush=True)\ntime.sleep(60)\n"
    )
    native.chmod(0o755)
    runner = subprocess.Popen(
        [sys.executable, "-c", "from spacepharer.engine import executable; executable()"],
        env={**os.environ, "SPACEPHARER_NATIVE": str(native)},
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    for _ in range(100):
        if pidfile.exists():
            break
        time.sleep(0.01)
    assert pidfile.exists()
    if interrupt == "sigterm":
        runner.terminate()
    assert runner.wait(timeout=15) != 0
    child = int(pidfile.read_text())
    for _ in range(30):
        if not child_is_running(child):
            break
        time.sleep(0.01)
    else:
        os.kill(child, signal.SIGKILL)
        pytest.fail("native version child survived interrupted verification")
