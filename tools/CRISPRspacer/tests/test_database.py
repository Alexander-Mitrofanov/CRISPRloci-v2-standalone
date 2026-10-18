import gzip
import json
from pathlib import Path
import random
import sqlite3

import pytest

from crisprspacer.cli import main
from crisprspacer.database import analyze_database, build_database, file_digest, open_database, search_database, validate_database, validate_database_result
from crisprspacer.io import InputError, LimitError, digest
from crisprspacer.search import reverse_complement


def target(seq, key="target"):
    return {"id": key, "sequence": seq, "length": len(seq), "sha256": digest(seq), "topology": "linear", "assembly_id": None}


def database(tmp_path, rows):
    source = tmp_path / "spacers.fasta"
    source.write_text("".join(f">{key}\n{seq}\n" for key, seq in rows))
    out = tmp_path / "db"
    ref = build_database(source, out, "fixture", "https://example.org/fixture.fasta", "test")
    return out, ref


def test_three_record_ncbi_sample_cli(tmp_path):
    source = Path(__file__).parents[1] / "examples/ncbi-spacerome-three.fasta"
    db = tmp_path / "db"
    assert main(["database-index", "--spacers", str(source), "--name", "sample-three-only",
                 "--source-url", "https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz",
                 "--source-release", "2017-04-25-sample", "--output", str(db)]) == 0
    ref = json.loads((db / "reference.json").read_text())
    assert ref["record_count"] == ref["eligible_record_count"] == 3
    seq = source.read_text().splitlines()[1]
    targets = tmp_path / "targets.fasta"
    targets.write_text(">forward\nNN" + seq + "NN\n>reverse\n" + reverse_complement(seq) + "\n")
    result = tmp_path / "result"
    assert main(["database-search", "--database", str(db), "--targets", str(targets),
                 "--expected-database-sha256", ref["database_sha256"], "--output", str(result)]) == 0
    report = json.loads((result / "analysis.json").read_text())
    assert [(r["target_id"], r["target_start"], r["relative_strand"]) for r in report["matches"]] == [("forward", 3, "+"), ("reverse", 1, "-")]
    assert main(["validate-database", str(db)]) == 0
    assert main(["validate-database-result", str(result), "--database", str(db)]) == 0


def test_eligibility_duplicate_occurrences_palindrome_and_ambiguity(tmp_path):
    seq = "ACGTTGCA" * 2
    pal = "ACGT" * 4
    db, ref = database(tmp_path, [("a", seq), ("duplicate", seq), ("pal", pal), ("short", "ACG"),
                                ("ambiguous", "N" * 20), ("long", "A" * 1001)])
    assert (ref["record_count"], ref["eligible_record_count"], ref["skipped_length"], ref["skipped_ambiguous"]) == (6, 3, 2, 1)
    report = analyze_database(db, [target(seq + "NN" + pal)])
    assert report["counts"]["matches"] == 3
    assert report["matches"][-1]["relative_strand"] == "both"
    assert analyze_database(db, [target("N"*25)])["outcome"] == "no_assessable_windows"
    assert analyze_database(db, [target("AAA")])["outcome"] == "no_assessable_windows"
    assert analyze_database(db, [target("T"*25)])["outcome"] == "no_matches"


def test_exact_search_does_not_cross_records_or_match_ambiguous_tails(tmp_path):
    seq = "AACCGGTT" * 3
    db, _ = database(tmp_path, [("spacer", seq)])
    report = analyze_database(db, [target(seq[:16], "left"), target(seq[16:], "right"), target(seq[:16] + "N"*8, "ambiguous")])
    assert report["matches"] == []
    assert report["counts"]["candidates_checked"] > 0


def test_randomized_exhaustive_oracle(tmp_path):
    rng = random.Random(401)
    rows = [(f"s{i}", "".join(rng.choices("ACGT", k=rng.randrange(16, 41)))) for i in range(90)]
    rows += [("duplicate", rows[0][1]), ("palindrome", "ACGT"*4)]
    db, _ = database(tmp_path, rows)
    targets = []
    for i in range(160):
        chosen = rng.choice(rows)[1]
        if i % 2:
            chosen = reverse_complement(chosen)
        if i % 3 == 0:
            chosen = chosen[:17] + "N" + chosen[18:]
        targets.append(target("".join(rng.choices("ACGTN", k=5)) + chosen + "".join(rng.choices("ACGTN", k=10)), f"t{i}"))
    actual = analyze_database(db, targets)
    observed = {(h["spacer_id"], h["target_id"], h["target_start"], h["relative_strand"]) for h in actual["matches"]}
    expected = set()
    for key, seq in rows:
        rc = reverse_complement(seq)
        orientations = [(seq, "both")] if rc == seq else [(seq, "+"), (rc, "-")]
        for t in targets:
            for oriented, strand in orientations:
                for pos in range(len(t["sequence"]) - len(seq) + 1):
                    if t["sequence"][pos:pos+len(seq)] == oriented:
                        expected.add((key, t["id"], pos+1, strand))
    assert expected == observed
    assert len(actual["matches"]) == len(expected)


def test_gzip_and_duplicate_rejection_atomicity(tmp_path):
    source = tmp_path / "source.gz"
    source.write_bytes(gzip.compress(b">same\nACGTACGTACGTACGT\n>same\nACGTACGTACGTACGT\n"))
    dest = tmp_path / "db"
    with pytest.raises(InputError, match="Duplicate"):
        build_database(source, dest, "fixture", "https://example.org/source", "test")
    assert not dest.exists()
    assert not dest.with_name("db.partial").exists()


def test_fail_closed_limits_identity_and_symlinks(tmp_path):
    db, ref = database(tmp_path, [("spacer", "A"*16)])
    with pytest.raises(InputError, match="pinned"):
        open_database(db, "0" * 64)
    with pytest.raises(LimitError, match="candidate"):
        analyze_database(db, [target("A"*30)], max_candidates=1)
    fasta = tmp_path / "targets.fa"
    fasta.write_text(">a\n" + "A"*30 + "\n")
    output = tmp_path / "result"
    with pytest.raises(LimitError, match="hit budget"):
        search_database(db, fasta, output, max_hits=1)
    assert not output.exists()
    with (db / "spacers.sqlite3").open("ab") as f:
        f.write(b"tampered")
    with pytest.raises(InputError, match="integrity"):
        open_database(db, ref["database_sha256"])


def test_no_eligible_and_replay_rejects_tampering(tmp_path):
    db, _ = database(tmp_path, [("short", "ACGT")])
    assert analyze_database(db, [target("A"*50)])["outcome"] == "no_assessable_windows"
    fasta = tmp_path / "target.fa"
    fasta.write_text(">t\n" + "A"*50 + "\n")
    output = tmp_path / "result"
    search_database(db, fasta, output)
    assert validate_database_result(output, db)["status"] == "valid"
    (output / "matches.tsv").write_text("tampered\n")
    with pytest.raises(InputError, match="integrity"):
        validate_database_result(output, db)


@pytest.mark.parametrize("mutation", ["missing", "extra", "sequence", "status", "prefix"])
def test_deep_validation_rejects_semantically_broken_new_snapshot(tmp_path, mutation):
    db, ref = database(tmp_path, [("s", "AACCGGTTAACCGGTTAA")])
    sql = {"missing": "DELETE FROM seeds", "extra": "INSERT INTO seeds VALUES ('AAAAAAAAAAAAAAAA','extra','+','AAAAAAAAAAAAAAAA')",
           "sequence": "UPDATE spacers SET sequence='AAAAAAAAAAAAAAAAAA'",
           "status": "UPDATE spacers SET status='skipped_ambiguous'",
           "prefix": "UPDATE seeds SET prefix='AAAAAAAAAAAAAAAA'"}[mutation]
    with sqlite3.connect(db / "spacers.sqlite3") as connection:
        connection.execute(sql)
    ref["database_sha256"] = file_digest(db / "spacers.sqlite3")
    (db / "reference.json").write_text(json.dumps(ref))
    with pytest.raises(InputError):
        validate_database(db, ref["database_sha256"])


def test_sqlite_query_is_interrupted_even_without_candidates(tmp_path):
    db, ref = database(tmp_path, [("s", "A"*16)])
    # Simulate a defective operator-pinned query plan; this must not hang a worker.
    with sqlite3.connect(db / "spacers.sqlite3") as connection:
        connection.execute("DROP TABLE seeds")
        connection.execute("""CREATE VIEW seeds AS
            WITH RECURSIVE slow(x) AS (SELECT 1 UNION ALL SELECT x+1 FROM slow WHERE x<100000000)
            SELECT printf('%016x',x) AS prefix, 's' AS spacer_id, '+' AS strand,
                   'AAAAAAAAAAAAAAAA' AS oriented FROM slow""")
    ref["database_sha256"] = file_digest(db / "spacers.sqlite3")
    (db / "reference.json").write_text(json.dumps(ref))
    with pytest.raises(LimitError, match="deadline"):
        analyze_database(db, [target("A"*16)], timeout=1)


def test_database_metadata_binds_json_types(tmp_path):
    directory, ref = database(tmp_path, [("query", "ACGTTGCAACGTTGCAAC")])
    with sqlite3.connect(directory / "spacers.sqlite3") as connection:
        metadata = json.loads(connection.execute("SELECT value FROM metadata").fetchone()[0])
        metadata["eligible_record_count"] = True
        connection.execute("UPDATE metadata SET value=?", (json.dumps(metadata),))
    ref["database_sha256"] = file_digest(directory / "spacers.sqlite3")
    (directory / "reference.json").write_text(json.dumps(ref))
    with pytest.raises(InputError, match="metadata"):
        open_database(directory)


@pytest.mark.parametrize("field", ["seed_length", "min_length", "max_length"])
def test_database_policy_requires_integer_receipt_fields(tmp_path, field):
    directory, ref = database(tmp_path, [("query", "ACGTTGCAACGTTGCAAC")])
    ref[field] = float(ref[field])
    with sqlite3.connect(directory / "spacers.sqlite3") as connection:
        connection.execute("UPDATE metadata SET value=?", (json.dumps({k: v for k, v in ref.items()
                            if k != "database_sha256"}),))
    ref["database_sha256"] = file_digest(directory / "spacers.sqlite3")
    (directory / "reference.json").write_text(json.dumps(ref))
    with pytest.raises(InputError):
        open_database(directory)
