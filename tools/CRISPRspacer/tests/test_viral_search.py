"""Synthetic references only: scientific evidence, not a virus corpus test."""

import gzip
import json
from pathlib import Path
import shutil
import sqlite3

import pytest

from crisprspacer.io import InputError, LimitError, digest, write_manifest
from crisprspacer.search import reverse_complement
from crisprspacer.viral_reference import build_viral_reference
from crisprspacer.viral_search import (_evidence, _queries, search_viruses, validate_viral_search)

SEQ = "ACGTTGCCATGATCGTACCGATTCGGCATAGTCAAGTCCTA"


def spacers(tmp_path, entries):
    path = tmp_path / "spacers.fasta"
    path.write_text("".join(f">{key}\n{seq}\n" for key, seq in entries))
    return path


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    if not shutil.which("blastn") or not shutil.which("makeblastdb") or not shutil.which("blastdbcmd"):
        pytest.skip("Local BLAST+ binaries required for synthetic integration test")
    root = tmp_path_factory.mktemp("synthetic-viral-reference")
    source = root / "synthetic.fasta"
    one_mismatch = SEQ[:20] + ("A" if SEQ[20] != "A" else "C") + SEQ[21:]
    terminal_mismatch = ("A" if SEQ[0] != "A" else "C") + SEQ[1:]
    source.write_text(f">SYNTH001.1 synthetic forward\nNN{SEQ}NN{SEQ}NN\n"
                      f">SYNTH002.1 synthetic reverse\nNN{reverse_complement(SEQ)}NN\n"
                      f">SYNTH003.1 synthetic internal substitution\n{one_mismatch}\n"
                      f">SYNTH004.1 synthetic terminal substitution\n{terminal_mismatch}\n")
    destination = root / "reference"
    build_viral_reference([source], destination, name="synthetic-fixture", source_release="fixture-1",
                          source_url="https://example.org/synthetic-fixture.fasta")
    return destination


def test_real_blast_source_binding_duplicates_strands_and_counts(tmp_path, reference):
    path = spacers(tmp_path, [("repeat-id", SEQ.lower()), ("repeat-id", SEQ), ("ambiguity", "N" * 30),
                              ("too-short", "ACGT")])
    destination = tmp_path / "result"
    report = search_viruses(path, reference, destination)
    assert report["execution_completed"] is True
    assert report["retrieval_completeness"] == "heuristic"
    assert report["counts"]["query_occurrences"] == 4
    assert report["counts"]["eligible_query_occurrences"] == 2
    assert report["counts"]["unique_eligible_query_sequences"] == 1
    assert {hit["relative_strand"] for hit in report["matches"]} == {"+", "-"}
    candidate = next(row for row in report["candidate_viruses"] if row["accession"] == "SYNTH001.1")
    assert candidate["query_occurrence_count"] == 2
    assert candidate["unique_query_sequence_count"] == 1
    assert candidate["location_count"] == 2
    assert candidate["accepted_hsps"] == 4
    assert all(hit["query_id"] == "repeat-id" and hit["taxid"] is None for hit in report["matches"])
    assert all(hit["subject_accession"] != "SYNTH004.1" for hit in report["matches"])
    assert report["counts"]["rejection_reasons"]["partial_query"] > 0
    assert validate_viral_search(destination, reference) == report
    assert validate_viral_search(destination) == report
    assert not (destination / "queries.fasta").exists()
    assert not (destination / "subjects.fasta").exists()


def test_exact_profile_excludes_internal_substitutions(tmp_path, reference):
    path = spacers(tmp_path, [("q", SEQ)])
    report = search_viruses(path, reference, tmp_path / "result", max_mismatches=0)
    assert {h["subject_accession"] for h in report["matches"]} == {"SYNTH001.1", "SYNTH002.1"}
    assert report["counts"]["rejection_reasons"]["too_many_substitutions"] > 0


def test_all_ineligible_skips_blast_but_retains_explicit_outcome(tmp_path, reference, monkeypatch):
    path = spacers(tmp_path, [("q", "N" * 18), ("short", "AT"), ("long", "A" * 81)])
    def forbidden(*args, **kwargs):
        pytest.fail("BLAST must not run with no eligible queries")
    monkeypatch.setattr("crisprspacer.viral_search.run_blast", forbidden)
    report = search_viruses(path, reference, tmp_path / "result")
    assert report["outcome"] == "no_eligible_queries"
    assert report["counts"]["raw_hsps"] == 0
    assert [q["skip_reasons"] for q in report["queries"]] == [["ambiguous_query"], ["unsupported_length"], ["unsupported_length"]]
    validate_viral_search(tmp_path / "result", reference)


@pytest.mark.parametrize("option,value", [("max_mismatches", -1), ("max_mismatches", 3),
                                         ("threads", 0), ("max_hits", 0), ("max_hsps", 0),
                                         ("timeout", 0), ("timeout", 86401), ("timeout", float("nan"))])
def test_invalid_policy_rejected_before_output(tmp_path, option, value):
    with pytest.raises(InputError):
        search_viruses(Path("missing"), Path("missing"), tmp_path / "out", **{option: value})
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("option", ["max_hits", "max_hsps"])
def test_limits_fail_atomically(tmp_path, reference, option):
    path = spacers(tmp_path, [("q", SEQ)])
    destination = tmp_path / "result"
    with pytest.raises(LimitError):
        search_viruses(path, reference, destination, **{option: 1})
    assert not destination.exists()
    assert not destination.with_name("result.partial").exists()


def test_reference_pin_mismatch(tmp_path, reference):
    path = spacers(tmp_path, [("q", SEQ)])
    with pytest.raises(InputError):
        search_viruses(path, reference, tmp_path / "result", expected_reference_sha256="0" * 64)
    assert not (tmp_path / "result").exists()


@pytest.mark.parametrize("content", [">empty\n", "ACGT\n", ">rna\nACGU\n", ">gap\nAC-GT\n"])
def test_malformed_queries_rejected(tmp_path, content):
    source = tmp_path / "queries.fasta"
    source.write_text(content)
    with pytest.raises(InputError):
        _queries(source)


def test_gzip_lowercase_and_duplicate_original_ids(tmp_path):
    source = tmp_path / "spacers.gz"
    raw = (f">same description\n{SEQ.lower()}\n>same duplicate description\n{SEQ}\n").encode()
    source.write_bytes(gzip.compress(raw))
    queries, sha = _queries(source)
    assert sha == digest(source.read_bytes())
    assert [q["original_id"] for q in queries] == ["same", "same"]
    assert [q["occurrence_id"] for q in queries] == ["q000000001", "q000000002"]
    assert all(q["sequence"] == SEQ for q in queries)


def catalog(tmp_path, sequence):
    reference = tmp_path / "catalog"
    reference.mkdir()
    connection = sqlite3.connect(reference / "subjects.sqlite3")
    connection.execute("CREATE TABLE subjects(id TEXT PRIMARY KEY, accession TEXT, title TEXT, sequence TEXT, length INTEGER, sha256 TEXT)")
    connection.execute("INSERT INTO subjects VALUES(?,?,?,?,?,?)", ("v000000001", "SYNTH.1", "synthetic", sequence, len(sequence), digest(sequence)))
    connection.commit()
    connection.close()
    return reference


def raw_line(query, subject, *, qstart=1, qend=None, sstart=1, send=None, qseq=None, sseq=None):
    qend = len(query) if qend is None else qend
    send = len(subject) if send is None else send
    qseq = query if qseq is None else qseq
    sseq = subject if sseq is None else sseq
    identities = sum(a == b and a != "-" for a, b in zip(qseq, sseq))
    mismatch = sum(a != b and "-" not in (a, b) for a, b in zip(qseq, sseq))
    gaps = qseq.count("-") + sseq.count("-")
    return "\t".join(map(str, ["q000000001", "v000000001", len(query), len(subject), qstart, qend,
                              sstart, send, len(qseq), identities, mismatch, gaps, "1e-5", 40, 20, qseq, sseq])) + "\n"


def evidence(tmp_path, query, subject, lines, *, mismatches=2):
    reference = catalog(tmp_path, subject)
    queries, _ = _queries(spacers(tmp_path, [("q", query)]))
    raw = tmp_path / "raw.tsv"
    raw.write_text(lines)
    return _evidence(raw, queries, reference, mismatches, 100, 100, 1)


def test_reverse_query_coordinates_preserve_relative_frame(tmp_path):
    line = raw_line(SEQ, SEQ, qstart=len(SEQ), qend=1, sstart=len(SEQ), send=1,
                    qseq=reverse_complement(SEQ), sseq=reverse_complement(SEQ))
    report = evidence(tmp_path, SEQ, SEQ, line)
    assert report["matches"][0]["relative_strand"] == "+"
    assert report["matches"][0]["query_start"] == len(SEQ)
    assert report["matches"][0]["subject_sequence_query_frame"] == SEQ
    assert report["queries"][0]["subject_limit_saturated"] is True


def test_disjoint_hsps_never_sum_to_full_query_coverage(tmp_path):
    cut = len(SEQ) // 2
    lines = raw_line(SEQ, SEQ, qend=cut, send=cut, qseq=SEQ[:cut], sseq=SEQ[:cut])
    lines += raw_line(SEQ, SEQ, qstart=cut + 1, sstart=cut + 1, qseq=SEQ[cut:], sseq=SEQ[cut:])
    report = evidence(tmp_path, SEQ, SEQ, lines)
    assert report["outcome"] == "all_hsps_filtered"
    assert report["matches"] == []
    assert report["counts"]["rejection_reasons"]["partial_query"] == 2


def test_subject_ambiguity_and_gaps_are_explicit_rejections(tmp_path):
    subject = SEQ[:20] + "N" + SEQ[21:]
    report = evidence(tmp_path, SEQ, subject, raw_line(SEQ, subject))
    assert report["counts"]["rejection_reasons"]["ambiguous_subject"] == 1
    other = tmp_path / "other"
    other.mkdir()
    subject = SEQ[:20] + "A" + SEQ[20:]
    line = raw_line(SEQ, subject, qseq=SEQ[:20] + "-" + SEQ[20:])
    report = evidence(other, SEQ, subject, line)
    assert report["counts"]["rejection_reasons"]["gapped"] == 1


def test_no_raw_hsps_distinct_from_filtered(tmp_path):
    report = evidence(tmp_path, SEQ, SEQ, "")
    assert report["outcome"] == "no_raw_hsps"
    assert report["queries"][0]["outcome"] == "no_raw_hsps"


def test_local_identifier_prefixes_are_bound_to_safe_inventory(tmp_path):
    line = raw_line(SEQ, SEQ).replace("q000000001\tv000000001", "lcl|q000000001\tlcl|v000000001")
    report = evidence(tmp_path, SEQ, SEQ, line)
    assert report["matches"][0]["subject_id"] == "v000000001"
    assert report["matches"][0]["query_occurrence_id"] == "q000000001"


@pytest.mark.parametrize("count", [0, 1, 2, 3])
def test_full_query_substitution_threshold_is_literal(tmp_path, count):
    bases = list(SEQ)
    for position in [12, 20, 28][:count]:
        bases[position] = "A" if bases[position] != "A" else "C"
    subject = "".join(bases)
    report = evidence(tmp_path, SEQ, subject, raw_line(SEQ, subject))
    assert bool(report["matches"]) is (count <= 2)
    if count <= 2:
        assert report["matches"][0]["substitutions"] == count


@pytest.mark.parametrize("field,value", [("qseqid", "unknown"), ("sseqid", "unknown"), ("qlen", "1"),
                                        ("slen", "1"), ("qstart", "0"), ("sstart", "0"),
                                        ("nident", "0"), ("evalue", "nan"), ("bitscore", "-1"),
                                        ("qseq", "A" * len(SEQ)), ("sseq", "A" * len(SEQ))])
def test_raw_hsps_are_verified_not_trusted(tmp_path, field, value):
    from crisprspacer.viral_search import OUTFMT_FIELDS
    fields = raw_line(SEQ, SEQ).rstrip().split("\t")
    fields[OUTFMT_FIELDS.index(field)] = value
    with pytest.raises(InputError):
        evidence(tmp_path, SEQ, SEQ, "\t".join(fields) + "\n")


@pytest.mark.parametrize("field,value", [("substitutions", 2), ("subject_start", 1),
                                        ("query_sequence", "A" * len(SEQ)), ("relative_strand", "-"),
                                        ("confidence", 0.99)])
def test_rehashed_result_tampering_fails_evidence_validation(tmp_path, reference, field, value):
    source = spacers(tmp_path, [("q", SEQ)])
    output = tmp_path / "result"
    search_viruses(source, reference, output)
    report = json.loads((output / "result.json").read_text())
    report["matches"][0][field] = ("+" if report["matches"][0][field] == "-" else "-") if field == "relative_strand" else value
    (output / "result.json").write_text(json.dumps(report))
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    with pytest.raises(InputError):
        validate_viral_search(output, reference)


@pytest.mark.parametrize("section,field,value", [
    ("settings", "task", "megablast"), ("settings", "strand", "plus"),
    ("settings", "max_target_seqs", 1), ("settings", "max_hsps", 1),
    ("settings", "num_threads", 0), ("settings", "outfmt", "6"),
    ("executable", "sha256", "invalid"), ("executable", "version", None),
    ("limits", "max_blast_output_bytes", 1), ("limits", "on_limit", "truncate"),
    ("limits", "timeout_seconds", 86401), ("limits", "max_raw_hsps", 0),
    ("policy", "min_query_length", 18.0), ("policy", "max_query_length", 80.0),
    ("counts", "accepted_hsps", 900), ("counts", "matched_unique_query_sequences", 900),
])
def test_rehashed_provenance_and_counts_rejected_without_reference(tmp_path, reference, section, field, value):
    source = spacers(tmp_path, [("q", SEQ)])
    output = tmp_path / "result"
    search_viruses(source, reference, output)
    report = json.loads((output / "result.json").read_text())
    target = report["blast"][section] if section in {"settings", "executable"} else report[section]
    target[field] = value
    (output / "result.json").write_text(json.dumps(report))
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    with pytest.raises(InputError):
        validate_viral_search(output)


def test_rehashed_raw_hsp_removal_fails_internal_replay(tmp_path, reference):
    source = spacers(tmp_path, [("q", SEQ)])
    output = tmp_path / "result"
    search_viruses(source, reference, output)
    raw = output / "blast-hsps.tsv"
    raw.write_text("".join(raw.read_text().splitlines(keepends=True)[1:]))
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    with pytest.raises(InputError):
        validate_viral_search(output)


def test_no_reference_validator_checks_combined_output_budget(tmp_path, reference, monkeypatch):
    source = spacers(tmp_path, [("q", SEQ)])
    output = tmp_path / "result"
    report = search_viruses(source, reference, output)
    raw_size = (output / "blast-hsps.tsv").stat().st_size
    monkeypatch.setattr("crisprspacer.viral_search.MAX_OUTPUT_BYTES", raw_size + 4)
    report["limits"]["max_blast_output_bytes"] = raw_size + 4
    (output / "blast.stderr.txt").write_text("12345")
    (output / "result.json").write_text(json.dumps(report))
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    with pytest.raises(LimitError):
        validate_viral_search(output)


@pytest.mark.parametrize("field,value", [("reward", True), ("ungapped", 1), ("soft_masking", 0)])
def test_blast_provenance_rejects_boolean_numeric_coercion(tmp_path, reference, field, value):
    path = spacers(tmp_path, [("query", SEQ)])
    output = tmp_path / "result"
    report = search_viruses(path, reference, output)
    report["blast"]["settings"][field] = value
    (output / "result.json").write_text(json.dumps(report))
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    with pytest.raises(InputError, match="settings"):
        validate_viral_search(output)


def test_ineligible_query_length_cannot_be_boolean(tmp_path, reference):
    path = spacers(tmp_path, [("short", "A")])
    output = tmp_path / "result"
    report = search_viruses(path, reference, output)
    inputs = json.loads((output / "inputs.json").read_text())
    inputs["queries"][0]["length"] = True
    report["queries"][0]["length"] = True
    (output / "inputs.json").write_text(json.dumps(inputs))
    (output / "result.json").write_text(json.dumps(report))
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    with pytest.raises(InputError, match="length"):
        validate_viral_search(output)


def test_validator_replays_the_verified_raw_hsp_snapshot(tmp_path, reference, monkeypatch):
    from crisprspacer import viral_search

    path = spacers(tmp_path, [("query", SEQ)])
    output = tmp_path / "result"
    search_viruses(path, reference, output)
    raw = output / "blast-hsps.tsv"
    original = raw.read_bytes()
    assert original
    # The verified bundle omits every HSP while its report still claims matches.
    raw.write_bytes(b"")
    (output / "manifest.json").unlink()
    write_manifest(output, "crisprspacer.viral-search")
    verify = viral_search.verify_manifest

    def replace_after_verification(*args):
        snapshot = verify(*args)
        raw.write_bytes(original)
        return snapshot

    monkeypatch.setattr(viral_search, "verify_manifest", replace_after_verification)
    with pytest.raises(InputError, match="replay"):
        validate_viral_search(output)
