"""Observed short-spacer BLAST HSPs against an immutable viral reference.

This is a heuristic local-alignment retrieval workflow.  Even an exact, full-query
acceptance profile is not an exhaustive search, and a match is not host evidence.
"""

from __future__ import annotations

from collections import Counter
import csv
import gzip
import io
import math
from pathlib import Path
import time
import zlib

from . import __version__
from .artifacts import write_table
from .blast_process import executable_identity, run_blast
from .io import (DNA, MAX_BASES, MAX_BYTES, MAX_RECORDS, InputError, LimitError,
                 atomic_directory, canonical, decode_json, digest, identifier, integer,
                 read_bytes, verify_manifest, write_json, write_manifest)
from .search import reverse_complement
from .viral_reference import bounded_catalog, load_viral_reference, open_subjects

SCHEMA = {"name": "crisprspacer.viral-search", "version": "1.0.0"}
FILES = {"result.json", "inputs.json", "matches.tsv", "blast-hsps.tsv", "blast.stderr.txt"}
OUTFMT_FIELDS = ["qseqid", "sseqid", "qlen", "slen", "qstart", "qend", "sstart", "send",
                 "length", "nident", "mismatch", "gaps", "evalue", "bitscore", "score", "qseq", "sseq"]
OUTFMT = "6 " + " ".join(OUTFMT_FIELDS)
MAX_OUTPUT_BYTES = 64 * 1024 * 1024
REASONS = ("partial_query", "gapped", "ambiguous_subject", "too_many_substitutions")
MATCH_FIELDS = ["id", "query_occurrence_id", "query_id", "query_sha256", "query_sequence",
                "subject_id", "subject_accession", "subject_title", "subject_sha256", "subject_length",
                "query_start", "query_end", "subject_start", "subject_end", "relative_strand",
                "query_length", "alignment_length", "identities", "substitutions", "gaps",
                "query_coverage", "query_alignment", "subject_alignment", "subject_sequence_forward",
                "subject_sequence_query_frame", "evalue", "bitscore", "raw_score", "taxid",
                "pam_status", "array_overlap", "confidence", "raw_hsp_row"]


def _queries(path: Path) -> tuple[list[dict], str]:
    """Keep every occurrence, including repeated original FASTA identifiers."""
    raw = read_bytes(path)
    stream = gzip.GzipFile(fileobj=io.BytesIO(raw)) if Path(path).suffix == ".gz" else io.BytesIO(raw)
    rows, chunks, current = [], [], None
    total = consumed = 0

    def finish():
        if current is None:
            return
        sequence = "".join(chunks)
        if not sequence:
            raise InputError("Empty spacer FASTA records are not supported")
        reasons = []
        if not 18 <= len(sequence) <= 80:
            reasons.append("unsupported_length")
        if set(sequence) - set("ACGT"):
            reasons.append("ambiguous_query")
        rows.append({"occurrence_id": f"q{len(rows) + 1:09d}", "original_id": current,
                     "sequence": sequence, "length": len(sequence), "sha256": digest(sequence),
                     "eligible": not reasons, "skip_reasons": reasons})
        if len(rows) > MAX_RECORDS:
            raise LimitError("Too many spacer occurrences")

    try:
        with stream as handle:
            while raw_line := handle.readline(1_000_001):
                consumed += len(raw_line)
                if len(raw_line) > 1_000_000 or consumed > MAX_BYTES:
                    raise LimitError("Spacer FASTA exceeds decompressed byte or line budget")
                line = raw_line.decode("ascii").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish()
                    header = line[1:].split()
                    current = identifier(header[0] if header else None, "FASTA ID")
                    chunks = []
                else:
                    if current is None:
                        raise InputError("FASTA sequence precedes its header")
                    sequence = line.upper()
                    if not set(sequence) <= DNA:
                        raise InputError("Only DNA IUPAC letters are supported; gaps/RNA are not accepted")
                    total += len(sequence)
                    if total > MAX_BASES:
                        raise LimitError("Spacer FASTA exceeds base budget")
                    chunks.append(sequence)
        finish()
    except (UnicodeError, EOFError, gzip.BadGzipFile, OSError, zlib.error) as exc:
        raise InputError("Malformed or corrupt spacer FASTA") from exc
    if not rows:
        raise InputError("Spacer FASTA contains no records")
    return rows, digest(raw)


def _policy(mismatches: int) -> dict:
    return {"profile": "full-query-ungapped-hsp", "max_substitutions": mismatches,
            "min_query_length": 18, "max_query_length": 80, "query_alphabet": "ACGT",
            "strands": "both", "topology": "linear", "pam": "not_evaluated",
            "array_overlap": "unknown", "host_inference": "not_performed"}


def _identity(reference: dict, reference_sha256: str) -> dict:
    return {key: reference[key] for key in
            ("name", "source_release", "source_url", "record_count", "total_bases", "content_sha256")} | {
                "reference_sha256": reference_sha256}


def _blast_settings(threads: int, record_count: int) -> dict:
    return {"task": "blastn-short", "word_size": 7, "strand": "both", "dust": "no",
            "soft_masking": False, "evalue": 1000, "ungapped": True, "reward": 1,
            "penalty": -3, "num_threads": threads, "max_target_seqs": record_count,
            "max_hsps": None, "outfmt": OUTFMT}


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and not set(value) - set("0123456789abcdef")


def _number(value: str, name: str) -> float:
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise InputError(f"Invalid BLAST {name}") from exc
    if not math.isfinite(number) or number < 0:
        raise InputError(f"Invalid BLAST {name}")
    return number


def _parse_hsp(line: str, queries: dict, connection, subject_metadata=None) -> tuple[dict, dict, str]:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != len(OUTFMT_FIELDS):
        raise InputError("BLAST output does not match the declared tabular format")
    raw = dict(zip(OUTFMT_FIELDS, fields))
    # Parsed local sequence identifiers may acquire NCBI's explicit local prefix.
    for key in ("qseqid", "sseqid"):
        if raw[key].startswith("lcl|"):
            raw[key] = raw[key][4:]
    query = queries.get(raw["qseqid"])
    if query is None or not query["eligible"]:
        raise InputError("BLAST returned an unknown or ineligible query")
    hsp = {"query_occurrence_id": raw["qseqid"], "subject_id": raw["sseqid"]}
    for field in OUTFMT_FIELDS[2:12]:
        try:
            hsp[field] = int(raw[field])
        except ValueError as exc:
            raise InputError("Invalid integer in BLAST output") from exc
    for field in ("evalue", "bitscore", "score"):
        hsp[field] = _number(raw[field], field)
    qstart, qend, sstart, send = (hsp[key] for key in ("qstart", "qend", "sstart", "send"))
    if hsp["qlen"] != query["length"] or not (1 <= qstart <= hsp["qlen"] and 1 <= qend <= hsp["qlen"]):
        raise InputError("BLAST query coordinates or length are inconsistent")
    if hsp["slen"] < 1 or not (1 <= sstart <= hsp["slen"] and 1 <= send <= hsp["slen"]):
        raise InputError("BLAST subject coordinates are invalid")
    start, end = sorted((sstart, send))
    qseq, sseq = raw["qseq"].upper(), raw["sseq"].upper()
    if (not qseq or len(qseq) != len(sseq) or len(qseq) != hsp["length"]
            or set(qseq + sseq) - (DNA | {"-"})
            or any(q == s == "-" for q, s in zip(qseq, sseq))):
        raise InputError("Invalid BLAST aligned sequence columns")
    if connection is not None:
        subject = connection.execute(
            "SELECT accession,title,length,sha256,substr(sequence, ?, ?) AS window FROM subjects WHERE id=?",
            (start, end - start + 1, raw["sseqid"])).fetchone()
    else:
        # Without a reference, only raw-HSP/internal evidence can be replayed.
        # A missing metadata row cannot become an accepted result on replay.
        subject = dict((subject_metadata or {}).get(raw["sseqid"], {
            "accession": None, "title": None, "sha256": None, "length": hsp["slen"]}))
        window = sseq.replace("-", "")
        subject["window"] = window if sstart <= send else reverse_complement(window)
    if subject is None or subject["length"] != hsp["slen"]:
        raise InputError("BLAST subject is absent from or inconsistent with the source inventory")
    if len(subject["window"]) != end - start + 1:
        raise InputError("BLAST subject span disagrees with aligned bases")
    query_window = query["sequence"][min(qstart, qend) - 1:max(qstart, qend)]
    if qstart > qend:
        query_window = reverse_complement(query_window)
    subject_window = subject["window"]
    if sstart > send:
        subject_window = reverse_complement(subject_window)
    if qseq.replace("-", "") != query_window or sseq.replace("-", "") != subject_window:
        raise InputError("BLAST aligned sequences do not reconstruct from their source intervals")
    identities = sum(q == s and q != "-" for q, s in zip(qseq, sseq))
    substitutions = sum(q != s and "-" not in (q, s) for q, s in zip(qseq, sseq))
    gaps = qseq.count("-") + sseq.count("-")
    if (hsp["nident"], hsp["mismatch"], hsp["gaps"]) != (identities, substitutions, gaps):
        raise InputError("BLAST alignment statistics do not agree with aligned bases")
    hsp.update(query_alignment=qseq, subject_alignment=sseq,
               relative_strand="+" if (qstart <= qend) == (sstart <= send) else "-",
               subject_sequence_forward=subject["window"])
    return hsp, dict(subject), query["occurrence_id"]


def _filter(hsp: dict, mismatches: int) -> list[str]:
    reasons = []
    if min(hsp["qstart"], hsp["qend"]) != 1 or max(hsp["qstart"], hsp["qend"]) != hsp["qlen"]:
        reasons.append("partial_query")
    if hsp["gaps"]:
        reasons.append("gapped")
    if set(hsp["subject_alignment"].replace("-", "")) - set("ACGT"):
        reasons.append("ambiguous_subject")
    if hsp["mismatch"] > mismatches:
        reasons.append("too_many_substitutions")
    return reasons


def _match(hsp: dict, subject: dict, query: dict, row: int, index: int) -> dict:
    forward = hsp["subject_sequence_forward"]
    return {"id": f"h{index:09d}", "query_occurrence_id": query["occurrence_id"],
            "query_id": query["original_id"], "query_sha256": query["sha256"],
            "query_sequence": query["sequence"], "subject_id": hsp["subject_id"],
            "subject_accession": subject["accession"], "subject_title": subject["title"],
            "subject_sha256": subject["sha256"], "subject_length": subject["length"],
            "query_start": hsp["qstart"], "query_end": hsp["qend"],
            "subject_start": min(hsp["sstart"], hsp["send"]),
            "subject_end": max(hsp["sstart"], hsp["send"]),
            "blast_subject_start": hsp["sstart"], "blast_subject_end": hsp["send"],
            "relative_strand": hsp["relative_strand"], "query_length": query["length"],
            "alignment_length": hsp["length"], "identities": hsp["nident"],
            "substitutions": hsp["mismatch"], "gaps": hsp["gaps"], "query_coverage": 1.0,
            "query_alignment": hsp["query_alignment"], "subject_alignment": hsp["subject_alignment"],
            "subject_sequence_forward": forward,
            "subject_sequence_query_frame": forward if hsp["relative_strand"] == "+" else reverse_complement(forward),
            "evalue": hsp["evalue"], "bitscore": hsp["bitscore"], "raw_score": hsp["score"],
            "taxid": None, "pam_status": "not_evaluated", "array_overlap": "unknown",
            "confidence": None, "raw_hsp_row": row}


def _aggregate(matches: list[dict]) -> list[dict]:
    rows = {}
    for hit in matches:
        key = hit["subject_accession"]
        row = rows.setdefault(key, {"accession": key, "title": hit["subject_title"],
                                   "subject_sha256": hit["subject_sha256"], "taxid": None,
                                   "query_occurrences": set(), "query_sequences": set(),
                                   "locations": set(), "accepted_hsps": 0, "min_substitutions": 2})
        if (row["title"], row["subject_sha256"]) != (hit["subject_title"], hit["subject_sha256"]):
            raise InputError("One accession maps to inconsistent source records")
        row["query_occurrences"].add(hit["query_occurrence_id"])
        row["query_sequences"].add(hit["query_sha256"])
        row["locations"].add((hit["subject_start"], hit["subject_end"], hit["relative_strand"]))
        row["accepted_hsps"] += 1
        row["min_substitutions"] = min(row["min_substitutions"], hit["substitutions"])
    result = []
    for key in sorted(rows):
        row = rows[key]
        result.append({"accession": row["accession"], "title": row["title"],
                       "subject_sha256": row["subject_sha256"], "taxid": None,
                       "query_occurrence_count": len(row["query_occurrences"]),
                       "unique_query_sequence_count": len(row["query_sequences"]),
                       "location_count": len(row["locations"]), "accepted_hsps": row["accepted_hsps"],
                       "min_substitutions": row["min_substitutions"]})
    return result


def _evidence(path: Path | bytes, queries: list[dict], reference: Path | None, mismatches: int,
              max_hsps: int, max_hits: int, subject_limit: int, checkpoint=lambda: None,
              subject_metadata=None) -> dict:
    mapping = {q["occurrence_id"]: q for q in queries}
    statistics = {key: {"raw_hsps": 0, "accepted_hsps": 0, "rejected_hsps": 0,
                        "rejection_reasons": dict.fromkeys(REASONS, 0), "subjects": set()}
                  for key in mapping}
    matches, raw_count, rejected_count = [], 0, 0
    connection = open_subjects(reference, checkpoint) if reference is not None else None
    try:
        # Validation passes the manifest-verified bytes; searches stream their new file.
        if isinstance(path, bytes):
            stream = io.TextIOWrapper(io.BytesIO(path), encoding="ascii")
        else:
            stream = path.open("r", encoding="ascii")
        with bounded_catalog(connection, checkpoint), stream as handle:
            for line in handle:
                checkpoint()
                raw_count += 1
                if raw_count > max_hsps:
                    raise LimitError("BLAST raw HSP row budget exceeded; no completed result exists")
                hsp, subject, key = _parse_hsp(line, mapping, connection, subject_metadata)
                stats = statistics[key]
                stats["raw_hsps"] += 1
                stats["subjects"].add(hsp["subject_id"])
                reasons = _filter(hsp, mismatches)
                if reasons:
                    rejected_count += 1
                    stats["rejected_hsps"] += 1
                    for reason in reasons:
                        stats["rejection_reasons"][reason] += 1
                else:
                    if subject["accession"] is None:
                        raise InputError("Accepted raw HSP has no reported subject metadata")
                    if len(matches) >= max_hits:
                        raise LimitError("Accepted HSP budget exceeded; no completed result exists")
                    stats["accepted_hsps"] += 1
                    matches.append(_match(hsp, subject, mapping[key], raw_count, len(matches) + 1))
    except UnicodeError as exc:
        raise InputError("BLAST output is not ASCII") from exc
    finally:
        if connection is not None:
            connection.close()
    summaries = []
    for query in queries:
        stats = statistics[query["occurrence_id"]]
        returned_subjects = len(stats.pop("subjects"))
        outcome = ("ineligible" if not query["eligible"] else "matches" if stats["accepted_hsps"]
                   else "no_raw_hsps" if not stats["raw_hsps"] else "all_hsps_filtered")
        summaries.append({"occurrence_id": query["occurrence_id"], "original_id": query["original_id"],
                          "sha256": query["sha256"], "length": query["length"],
                          "eligible": query["eligible"], "skip_reasons": query["skip_reasons"],
                          "outcome": outcome, "returned_subjects": returned_subjects,
                          "subject_limit_saturated": returned_subjects >= subject_limit, **stats})
    counts = {"query_occurrences": len(queries), "eligible_query_occurrences": sum(q["eligible"] for q in queries),
              "unique_eligible_query_sequences": len({q["sha256"] for q in queries if q["eligible"]}),
              "raw_hsps": raw_count, "accepted_hsps": len(matches), "rejected_hsps": rejected_count,
              "matched_query_occurrences": len({h["query_occurrence_id"] for h in matches}),
              "matched_unique_query_sequences": len({h["query_sha256"] for h in matches}),
              "matched_subject_accessions": len({h["subject_accession"] for h in matches}),
              "rejection_reasons": dict(Counter({r: sum(q["rejection_reasons"][r] for q in summaries)
                                                 for r in REASONS}))}
    outcome = ("no_eligible_queries" if not counts["eligible_query_occurrences"] else "matches" if matches
               else "no_raw_hsps" if not raw_count else "all_hsps_filtered")
    return {"queries": summaries, "matches": matches, "candidate_viruses": _aggregate(matches),
            "counts": counts, "outcome": outcome}


def search_viruses(spacers: Path, reference: Path, destination: Path, *, blastn="blastn",
                   expected_reference_sha256=None, max_mismatches=2, timeout=300,
                   threads=2, max_hits=10000, max_hsps=200000) -> dict:
    integer(max_mismatches, "max_mismatches", 0, 2)
    integer(threads, "threads", 1, 128)
    integer(max_hits, "max_hits", 1, 100000)
    integer(max_hsps, "max_hsps", 1, 2000000)
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or not 0 < timeout <= 86400:
        raise InputError("timeout must be a positive finite number at most 86400 seconds")
    deadline = time.monotonic() + timeout

    def checkpoint():
        if time.monotonic() >= deadline:
            raise LimitError("Viral search deadline exceeded; no completed result exists")

    queries, source_sha = _queries(spacers)
    ref, ref_sha = load_viral_reference(reference, expected_reference_sha256, checkpoint=checkpoint)
    reference_id = _identity(ref, ref_sha)
    executable = executable_identity(blastn)
    checkpoint()
    settings = _blast_settings(threads, ref["record_count"])
    limits = {"timeout_seconds": timeout, "max_accepted_hsps": max_hits, "max_raw_hsps": max_hsps,
              "max_blast_output_bytes": MAX_OUTPUT_BYTES, "on_limit": "fail_without_completed_result"}
    inputs = {"schema": SCHEMA, "spacer_fasta_sha256": source_sha, "queries": queries,
              "reference": reference_id, "policy": _policy(max_mismatches)}
    with atomic_directory(destination) as work:
        query_path = work / "queries.fasta"
        query_path.write_text("".join(f">{q['occurrence_id']}\n{q['sequence']}\n" for q in queries if q["eligible"]),
                              encoding="ascii")
        raw_path, stderr_path = work / "blast-hsps.tsv", work / "blast.stderr.txt"
        if any(q["eligible"] for q in queries):
            command = [executable["path"], "-task", "blastn-short", "-query", str(query_path),
                       "-db", str(Path(reference).absolute() / "db"), "-word_size", "7", "-strand", "both",
                       "-dust", "no", "-soft_masking", "false", "-evalue", "1000", "-ungapped",
                       "-reward", "1", "-penalty", "-3", "-num_threads", str(threads),
                       "-max_target_seqs", str(ref["record_count"]), "-outfmt", OUTFMT]
            run_blast(command, cwd=work, stdout=raw_path, stderr=stderr_path,
                      timeout=max(1, math.ceil(deadline - time.monotonic())), max_output_bytes=MAX_OUTPUT_BYTES)
        else:
            raw_path.touch()
            stderr_path.touch()
        query_path.unlink()
        checkpoint()
        if raw_path.stat().st_size + stderr_path.stat().st_size > MAX_OUTPUT_BYTES:
            raise LimitError("BLAST combined output exceeds byte budget")
        evidence = _evidence(raw_path, queries, reference, max_mismatches, max_hsps, max_hits,
                             ref["record_count"], checkpoint)
        report = {"schema": SCHEMA, "tool_version": __version__, "execution_completed": True,
                  "retrieval_completeness": "heuristic", "reference": reference_id,
                  "policy": _policy(max_mismatches), "blast": {"executable": executable, "settings": settings},
                  "limits": limits, "coordinate_system": "1-based-inclusive",
                  "limitations": ["BLAST retrieval is heuristic; no-hit results do not establish absence.",
                                  "Local alignments may clip terminal mismatches and fail the full-query filter.",
                                  "Sequence matches do not establish infection, host range or functional targeting.",
                                  "E-values and bit scores describe the observed HSP in this database snapshot.",
                                  "Taxonomy, PAM compatibility and CRISPR-array overlap are not evaluated."],
                  **evidence}
        write_json(work / "inputs.json", inputs)
        write_json(work / "result.json", report)
        write_table(work / "matches.tsv", report["matches"], MATCH_FIELDS)
        write_manifest(work, "crisprspacer.viral-search")
        checkpoint()
    return report


def validate_viral_search(directory: Path, reference: Path | None = None,
                          expected_reference_sha256=None) -> dict:
    """Validate integrity and source-bound evidence; reference enables full HSP replay."""
    _, payloads = verify_manifest(directory, "crisprspacer.viral-search", FILES)
    report, inputs = (decode_json(payloads[name]) for name in ("result.json", "inputs.json"))
    if not isinstance(report, dict) or not isinstance(inputs, dict) or report.get("schema") != SCHEMA or inputs.get("schema") != SCHEMA:
        raise InputError("Unsupported viral-search schema")
    if set(inputs) != {"schema", "spacer_fasta_sha256", "queries", "reference", "policy"} or not _sha(inputs["spacer_fasta_sha256"]):
        raise InputError("Invalid viral-search input provenance")
    if set(report) != {"schema", "tool_version", "execution_completed", "retrieval_completeness", "reference",
                       "policy", "blast", "limits", "coordinate_system", "limitations", "queries",
                       "matches", "candidate_viruses", "counts", "outcome"} or report["tool_version"] not in {"0.3.0", __version__}:
        raise InputError("Unsupported viral-search result fields/version")
    if (report.get("execution_completed") is not True or report.get("retrieval_completeness") != "heuristic"
            or report.get("coordinate_system") != "1-based-inclusive"
            or report.get("policy") != inputs.get("policy") or report.get("reference") != inputs.get("reference")):
        raise InputError("Viral-search evidence contract mismatch")
    policy = report["policy"]
    if not isinstance(policy, dict):
        raise InputError("Missing viral-search policy")
    mismatches = integer(policy.get("max_substitutions"), "max_substitutions", 0, 2)
    if canonical(policy) != canonical(_policy(mismatches)):
        raise InputError("Unsupported viral-search policy")
    reference_id = report["reference"]
    if (not isinstance(reference_id, dict) or set(reference_id) != {
            "name", "source_release", "source_url", "record_count", "total_bases", "content_sha256", "reference_sha256"}
            or not _sha(reference_id["content_sha256"]) or not _sha(reference_id["reference_sha256"])):
        raise InputError("Invalid viral-search reference identity")
    identifier(reference_id["name"])
    identifier(reference_id["source_release"])
    if not isinstance(reference_id["source_url"], str) or not reference_id["source_url"].startswith("https://"):
        raise InputError("Invalid reference source URL")
    subject_limit = integer(reference_id["record_count"], "reference record count", 1, 10_000_000)
    integer(reference_id["total_bases"], "reference base count", 1, 2**63 - 1)
    blast = report.get("blast")
    if not isinstance(blast, dict) or set(blast) != {"executable", "settings"} or not isinstance(blast["settings"], dict):
        raise InputError("Missing BLAST provenance")
    threads = integer(blast["settings"].get("num_threads"), "threads", 1, 128)
    if canonical(blast["settings"]) != canonical(_blast_settings(threads, subject_limit)):
        raise InputError("BLAST retrieval settings mismatch")
    executable = blast["executable"]
    if (not isinstance(executable, dict) or set(executable) != {"path", "sha256", "version"}
            or not _sha(executable["sha256"]) or not isinstance(executable["path"], str)
            or not isinstance(executable["version"], str) or not executable["version"]):
        raise InputError("Invalid BLAST executable identity")
    limits = report.get("limits")
    if (not isinstance(limits, dict) or set(limits) != {"timeout_seconds", "max_accepted_hsps", "max_raw_hsps",
                                                     "max_blast_output_bytes", "on_limit"}
            or limits["max_blast_output_bytes"] != MAX_OUTPUT_BYTES
            or limits["on_limit"] != "fail_without_completed_result"):
        raise InputError("Invalid viral-search resource policy")
    timeout = limits["timeout_seconds"]
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or not 0 < timeout <= 86400:
        raise InputError("Invalid viral-search timeout")
    max_hsps = integer(limits.get("max_raw_hsps"), "max_raw_hsps", 1, 2000000)
    max_hits = integer(limits.get("max_accepted_hsps"), "max_accepted_hsps", 1, 100000)
    if len(payloads["blast-hsps.tsv"]) + len(payloads["blast.stderr.txt"]) > MAX_OUTPUT_BYTES:
        raise LimitError("BLAST combined output exceeds byte budget")
    queries = inputs.get("queries")
    if not isinstance(queries, list) or not 1 <= len(queries) <= MAX_RECORDS:
        raise InputError("Invalid viral-search query inventory")
    mapping = {}
    total_bases = 0
    for index, query in enumerate(queries, 1):
        if (not isinstance(query, dict) or set(query) != {"occurrence_id", "original_id", "sequence", "length",
                                                         "sha256", "eligible", "skip_reasons"}
                or query.get("occurrence_id") != f"q{index:09d}"):
            raise InputError("Invalid query occurrence identity")
        identifier(query.get("original_id"))
        seq = query.get("sequence")
        if not isinstance(seq, str) or not seq or set(seq) - DNA or query.get("sha256") != digest(seq) or query.get("length") != len(seq):
            raise InputError("Invalid query sequence/hash/length")
        integer(query["length"], "query length", 1, MAX_BASES)
        reasons = ([] if 18 <= len(seq) <= 80 else ["unsupported_length"]) + ([] if set(seq) <= set("ACGT") else ["ambiguous_query"])
        if query.get("skip_reasons") != reasons or query.get("eligible") is not (not reasons):
            raise InputError("Query eligibility mismatch")
        mapping[query["occurrence_id"]] = query
        total_bases += len(seq)
        if total_bases > MAX_BASES:
            raise LimitError("Query inventory exceeds base budget")
    matches = report.get("matches")
    if not isinstance(matches, list) or len(matches) > 100000:
        raise InputError("Invalid viral-search match list")
    subject_metadata = {}
    for index, hit in enumerate(matches, 1):
        if (not isinstance(hit, dict) or set(hit) != set(MATCH_FIELDS) | {"blast_subject_start", "blast_subject_end"}
                or hit.get("id") != f"h{index:09d}"):
            raise InputError("Invalid match identity")
        query = mapping.get(hit.get("query_occurrence_id"))
        if query is None or not query["eligible"]:
            raise InputError("Match has no eligible query")
        if (hit.get("query_id"), hit.get("query_sequence"), hit.get("query_sha256"), hit.get("query_length")) != (
                query["original_id"], query["sequence"], query["sha256"], query["length"]):
            raise InputError("Match query binding mismatch")
        forward = hit.get("subject_sequence_forward")
        if not isinstance(forward, str) or len(forward) != query["length"] or set(forward) - set("ACGT"):
            raise InputError("Invalid match subject window")
        strand = hit.get("relative_strand")
        if strand not in {"+", "-"}:
            raise InputError("Invalid match relative strand")
        oriented = forward if strand == "+" else reverse_complement(forward)
        substitutions = sum(a != b for a, b in zip(query["sequence"], oriented))
        if (hit.get("subject_sequence_query_frame") != oriented or hit.get("substitutions") != substitutions
                or substitutions > mismatches or hit.get("identities") != len(forward) - substitutions
                or hit.get("gaps") != 0 or hit.get("query_coverage") != 1.0 or hit.get("alignment_length") != len(forward)):
            raise InputError("Match alignment evidence mismatch")
        integer(hit.get("subject_length"), "subject length", len(forward), 2**63 - 1)
        start = integer(hit.get("subject_start"), "subject start", 1, hit["subject_length"])
        end = integer(hit.get("subject_end"), "subject end", start, hit["subject_length"])
        if end - start + 1 != len(forward):
            raise InputError("Match coordinate span mismatch")
        if (hit.get("taxid") is not None or hit.get("confidence") is not None or
                hit.get("pam_status") != "not_evaluated" or hit.get("array_overlap") != "unknown"):
            raise InputError("Unsupported inferred match annotations")
        for name in ("evalue", "bitscore", "raw_score"):
            value = hit.get(name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise InputError("Invalid match alignment score")
        identifier(hit["subject_id"])
        identifier(hit["subject_accession"])
        if not _sha(hit["subject_sha256"]) or not isinstance(hit["subject_title"], str):
            raise InputError("Invalid match subject provenance")
        subject = {"accession": hit["subject_accession"], "title": hit["subject_title"],
                   "sha256": hit["subject_sha256"], "length": hit["subject_length"]}
        if hit["subject_id"] in subject_metadata and subject_metadata[hit["subject_id"]] != subject:
            raise InputError("Inconsistent subject metadata")
        subject_metadata[hit["subject_id"]] = subject
    if canonical(report.get("candidate_viruses")) != canonical(_aggregate(matches)):
        raise InputError("Candidate-virus aggregation mismatch")
    expected_rows = io.StringIO(newline="")
    from .artifacts import cell
    writer = csv.writer(expected_rows, delimiter="\t", lineterminator="\n")
    writer.writerow(MATCH_FIELDS)
    for hit in matches:
        writer.writerow([cell(hit.get(key)) for key in MATCH_FIELDS])
    if expected_rows.getvalue().encode() != payloads["matches.tsv"]:
        raise InputError("Match table does not match JSON evidence")
    if expected_reference_sha256 is not None and report["reference"].get("reference_sha256") != expected_reference_sha256:
        raise InputError("Viral-search reference digest mismatch")
    if reference is not None:
        ref, sha = load_viral_reference(reference, report["reference"].get("reference_sha256"))
        if canonical(_identity(ref, sha)) != canonical(report["reference"]):
            raise InputError("Viral-search reference identity mismatch")
    replay = _evidence(payloads["blast-hsps.tsv"], queries, reference, mismatches,
                       max_hsps, max_hits, subject_limit, subject_metadata=subject_metadata)
    if any(canonical(report.get(key)) != canonical(value) for key, value in replay.items()):
        raise InputError("Viral-search HSP evidence replay mismatch")
    return report
