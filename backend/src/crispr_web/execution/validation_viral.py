"""Independent source-bound viral evidence gate; no scientific-package imports."""

from collections import Counter, defaultdict
import hashlib
import json
import math
import sqlite3
from pathlib import Path
from urllib.parse import urlsplit

from .validation_annotations import _payload, _strict_json, _fasta

DNA = set("ACGTRYSWKMBDHVN")
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
REASONS = ("partial_query", "gapped", "ambiguous_subject", "too_many_substitutions")


def sha(value):
    return hashlib.sha256(
        value.encode("ascii") if isinstance(value, str) else value
    ).hexdigest()


def same(left, right):
    """Compare structured evidence without accepting booleans as integer counts."""
    return json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)


def rc(value):
    return value.translate(COMP)[::-1]


def _catalog(reference, expected_sha, checkpoint):
    raw = _payload(reference, reference / "reference.json", checkpoint)
    if sha(raw) != expected_sha:
        raise ValueError("Unpinned viral reference")
    meta = _strict_json(raw)
    path = reference / "subjects.sqlite3"
    if any(p.is_symlink() for p in (path, *path.parents)) or not path.is_file():
        raise ValueError("Unsafe viral catalog")
    wanted = meta["files"]["subjects.sqlite3"]
    if (
        path.stat().st_size != wanted["bytes"]
        or not 0 < path.stat().st_size <= 80 * 1024**3
    ):
        raise ValueError("Viral catalog size mismatch")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024**2), b""):
            if checkpoint:
                checkpoint()
            digest.update(chunk)
    if digest.hexdigest() != wanted["sha256"]:
        raise ValueError("Viral catalog checksum mismatch")
    connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA trusted_schema=OFF")
    connection.execute("PRAGMA query_only=ON")
    return meta, connection


def validate_viral_pipeline(
    claimed, root, expected, integration, input_root, reference, checkpoint=None
):
    connection = None
    try:
        if (
            input_root is None
            or reference is None
            or claimed.molecule != "DNA"
            or integration.get("workflow") != "viral_search"
        ):
            return False
        reference = Path(reference).absolute()
        source = {}
        for path in sorted(input_root.iterdir()):
            key, sequence = _fasta(_payload(input_root, path, checkpoint))
            if (
                key in source
                or expected.get(key) != (len(sequence), sha(sequence))
                or len(sequence) > 80
            ):
                return False
            source[key] = sequence
        if set(source) != set(expected) or not source:
            return False
        native = root / "viral"
        manifest = _strict_json(_payload(root, native / "manifest.json", checkpoint))
        filenames = {
            "result.json",
            "inputs.json",
            "matches.tsv",
            "blast-hsps.tsv",
            "blast.stderr.txt",
        }
        if (
            manifest.get("schema_version") != "1.0.0"
            or manifest.get("kind") != "crisprspacer.viral-search"
            or set(manifest["files"]) != filenames
        ):
            return False
        payloads = {
            name: _payload(root, native / name, checkpoint) for name in filenames
        }
        if any(
            manifest["files"][n] != {"sha256": sha(b), "bytes": len(b)}
            for n, b in payloads.items()
        ):
            return False
        doc = _strict_json(payloads["result.json"])
        inputs = _strict_json(payloads["inputs.json"])
        schema = {"name": "crisprspacer.viral-search", "version": "1.0.0"}
        if (
            doc.get("schema") != schema
            or inputs.get("schema") != schema
            or doc.get("tool_version") not in {"0.3.0", "0.3.1"}
            or doc.get("execution_completed") is not True
            or doc.get("retrieval_completeness") != "heuristic"
            or doc.get("coordinate_system") != "1-based-inclusive"
        ):
            return False
        policy = {
            "profile": "full-query-ungapped-hsp",
            "max_substitutions": claimed.viral_max_mismatches,
            "min_query_length": 18,
            "max_query_length": 80,
            "query_alphabet": "ACGT",
            "strands": "both",
            "topology": "linear",
            "pam": "not_evaluated",
            "array_overlap": "unknown",
            "host_inference": "not_performed",
        }
        if not same(doc["policy"], policy) or not same(inputs["policy"], policy):
            return False
        meta, connection = _catalog(
            reference, claimed.viral_reference_sha256, checkpoint
        )
        identity = {
            k: meta[k]
            for k in (
                "name",
                "source_release",
                "source_url",
                "record_count",
                "total_bases",
                "content_sha256",
            )
        }
        identity["reference_sha256"] = claimed.viral_reference_sha256
        if not same(doc["reference"], identity) or not same(
            inputs["reference"], identity
        ):
            return False
        url = urlsplit(identity["source_url"])
        if url.scheme != "https" or not url.netloc or url.username or url.password:
            return False
        combined = "".join(f">{k}\n{v}\n" for k, v in source.items())
        queries = []
        for i, (key, sequence) in enumerate(source.items(), 1):
            reasons = ([] if 18 <= len(sequence) <= 80 else ["unsupported_length"]) + (
                [] if set(sequence) <= set("ACGT") else ["ambiguous_query"]
            )
            queries.append(
                {
                    "occurrence_id": f"q{i:09d}",
                    "original_id": key,
                    "sequence": sequence,
                    "length": len(sequence),
                    "sha256": sha(sequence),
                    "eligible": not reasons,
                    "skip_reasons": reasons,
                }
            )
        if not same(inputs["queries"], queries) or inputs["spacer_fasta_sha256"] != sha(
            combined
        ):
            return False
        settings = doc["blast"]["settings"]
        threads = settings.get("num_threads")
        outfmt = "6 qseqid sseqid qlen slen qstart qend sstart send length nident mismatch gaps evalue bitscore score qseq sseq"
        if (
            type(threads) is not int
            or not 1 <= threads <= 8
            or settings
            != {
                "task": "blastn-short",
                "word_size": 7,
                "strand": "both",
                "dust": "no",
                "soft_masking": False,
                "evalue": 1000,
                "ungapped": True,
                "reward": 1,
                "penalty": -3,
                "num_threads": threads,
                "max_target_seqs": meta["record_count"],
                "max_hsps": None,
                "outfmt": outfmt,
            }
        ):
            return False
        limits = doc["limits"]
        if (
            set(limits)
            != {
                "timeout_seconds",
                "max_accepted_hsps",
                "max_raw_hsps",
                "max_blast_output_bytes",
                "on_limit",
            }
            or type(limits["timeout_seconds"]) not in (int, float)
            or not math.isfinite(limits["timeout_seconds"])
            or not 0 < limits["timeout_seconds"] <= 300
            or limits["max_blast_output_bytes"] != 64 * 1024**2
            or limits["max_accepted_hsps"] != 10000
            or limits["max_raw_hsps"] != 200000
            or limits["on_limit"] != "fail_without_completed_result"
        ):
            return False
        lookup = {q["occurrence_id"]: q for q in queries}
        stats = {
            q["occurrence_id"]: {
                "raw_hsps": 0,
                "accepted_hsps": 0,
                "rejected_hsps": 0,
                "rejection_reasons": dict.fromkeys(REASONS, 0),
                "subjects": set(),
            }
            for q in queries
        }
        matches = doc["matches"]
        accepted = []
        if not isinstance(matches, list) or len(matches) > 10000:
            return False
        raw_lines = payloads["blast-hsps.tsv"].decode("ascii").splitlines()
        if len(raw_lines) > 200000:
            return False
        for index, line in enumerate(raw_lines, 1):
            if checkpoint:
                checkpoint()
            f = line.split("\t")
            if len(f) != 17:
                return False
            q = lookup[f[0]]
            subject_id = f[1]
            qlen, slen, qs, qe, ss, se, length, nident, mismatch, gaps = map(
                int, f[2:12]
            )
            ev, bit, score = map(float, f[12:15])
            qseq, sseq = f[15:17]
            if (
                not q["eligible"]
                or qlen != q["length"]
                or not 1 <= qs <= qlen
                or not 1 <= qe <= qlen
                or not 1 <= ss <= slen
                or not 1 <= se <= slen
                or not 1 <= length <= 160
            ):
                return False
            if (
                len(qseq) != length
                or len(sseq) != length
                or set(qseq + sseq) - (DNA | {"-"})
                or any(not math.isfinite(x) or x < 0 for x in (ev, bit, score))
            ):
                return False
            start, end = sorted((ss, se))
            subject = connection.execute(
                "SELECT accession,title,length,sha256,substr(sequence,?,?) AS window FROM subjects WHERE id=?",
                (start, end - start + 1, subject_id),
            ).fetchone()
            if subject is None or subject["length"] != slen:
                return False
            qwindow = q["sequence"][min(qs, qe) - 1 : max(qs, qe)]
            if qseq.replace("-", "") != (
                qwindow if qs <= qe else rc(qwindow)
            ) or sseq.replace("-", "") != (
                subject["window"] if ss <= se else rc(subject["window"])
            ):
                return False
            identities = sum(a == b and a != "-" for a, b in zip(qseq, sseq))
            subs = sum(a != b and "-" not in (a, b) for a, b in zip(qseq, sseq))
            gap = qseq.count("-") + sseq.count("-")
            if (nident, mismatch, gaps) != (identities, subs, gap):
                return False
            reasons = []
            if min(qs, qe) != 1 or max(qs, qe) != qlen:
                reasons.append("partial_query")
            if gaps:
                reasons.append("gapped")
            if set(sseq.replace("-", "")) - set("ACGT"):
                reasons.append("ambiguous_subject")
            if mismatch > claimed.viral_max_mismatches:
                reasons.append("too_many_substitutions")
            stat = stats[q["occurrence_id"]]
            stat["raw_hsps"] += 1
            stat["subjects"].add(subject_id)
            if reasons:
                stat["rejected_hsps"] += 1
                for reason in reasons:
                    stat["rejection_reasons"][reason] += 1
                continue
            stat["accepted_hsps"] += 1
            strand = "+" if (qs <= qe) == (ss <= se) else "-"
            hit = {
                "id": f"h{len(accepted) + 1:09d}",
                "query_occurrence_id": q["occurrence_id"],
                "query_id": q["original_id"],
                "query_sha256": q["sha256"],
                "query_sequence": q["sequence"],
                "subject_id": subject_id,
                "subject_accession": subject["accession"],
                "subject_title": subject["title"],
                "subject_sha256": subject["sha256"],
                "subject_length": slen,
                "query_start": qs,
                "query_end": qe,
                "subject_start": start,
                "subject_end": end,
                "blast_subject_start": ss,
                "blast_subject_end": se,
                "relative_strand": strand,
                "query_length": qlen,
                "alignment_length": length,
                "identities": nident,
                "substitutions": subs,
                "gaps": 0,
                "query_coverage": 1.0,
                "query_alignment": qseq,
                "subject_alignment": sseq,
                "subject_sequence_forward": subject["window"],
                "subject_sequence_query_frame": subject["window"]
                if strand == "+"
                else rc(subject["window"]),
                "evalue": ev,
                "bitscore": bit,
                "raw_score": score,
                "taxid": None,
                "pam_status": "not_evaluated",
                "array_overlap": "unknown",
                "confidence": None,
                "raw_hsp_row": index,
            }
            accepted.append(hit)
        if json.dumps(accepted, sort_keys=True) != json.dumps(matches, sort_keys=True):
            return False
        summaries = []
        for q in queries:
            s = stats[q["occurrence_id"]]
            returned = len(s.pop("subjects"))
            summaries.append(
                {
                    k: q[k]
                    for k in (
                        "occurrence_id",
                        "original_id",
                        "sha256",
                        "length",
                        "eligible",
                        "skip_reasons",
                    )
                }
                | {
                    "outcome": "ineligible"
                    if not q["eligible"]
                    else "matches"
                    if s["accepted_hsps"]
                    else "no_raw_hsps"
                    if not s["raw_hsps"]
                    else "all_hsps_filtered",
                    "returned_subjects": returned,
                    "subject_limit_saturated": returned >= meta["record_count"],
                    **s,
                }
            )
        if not same(doc["queries"], summaries):
            return False
        counts = {
            "query_occurrences": len(queries),
            "eligible_query_occurrences": sum(q["eligible"] for q in queries),
            "unique_eligible_query_sequences": len(
                {q["sha256"] for q in queries if q["eligible"]}
            ),
            "raw_hsps": len(raw_lines),
            "accepted_hsps": len(matches),
            "rejected_hsps": len(raw_lines) - len(matches),
            "matched_query_occurrences": len(
                {h["query_occurrence_id"] for h in matches}
            ),
            "matched_unique_query_sequences": len({h["query_sha256"] for h in matches}),
            "matched_subject_accessions": len(
                {h["subject_accession"] for h in matches}
            ),
            "rejection_reasons": dict(
                Counter(
                    {
                        r: sum(q["rejection_reasons"][r] for q in summaries)
                        for r in REASONS
                    }
                )
            ),
        }
        if not same(doc["counts"], counts):
            return False
        candidates = []
        grouped = defaultdict(list)
        for hit in matches:
            grouped[hit["subject_accession"]].append(hit)
        for accession, rows in sorted(grouped.items()):
            if checkpoint:
                checkpoint()
            first = rows[0]
            candidates.append(
                {
                    "accession": accession,
                    "title": first["subject_title"],
                    "subject_sha256": first["subject_sha256"],
                    "taxid": None,
                    "query_occurrence_count": len(
                        {h["query_occurrence_id"] for h in rows}
                    ),
                    "unique_query_sequence_count": len(
                        {h["query_sha256"] for h in rows}
                    ),
                    "location_count": len(
                        {
                            (h["subject_start"], h["subject_end"], h["relative_strand"])
                            for h in rows
                        }
                    ),
                    "accepted_hsps": len(rows),
                    "min_substitutions": min(h["substitutions"] for h in rows),
                }
            )
        if not same(doc["candidate_viruses"], candidates):
            return False
        outcome = (
            "no_eligible_queries"
            if not counts["eligible_query_occurrences"]
            else "matches"
            if matches
            else "no_raw_hsps"
            if not raw_lines
            else "all_hsps_filtered"
        )
        stage = integration["stages"]["crisprspacer"]
        return (
            doc["outcome"] == outcome
            and stage.get("status") == "completed"
            and same(stage.get("counts"), counts)
            and stage.get("reference_sha256") == claimed.viral_reference_sha256
        )
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        AttributeError,
        IndexError,
        RecursionError,
        sqlite3.Error,
    ):
        return False
    finally:
        if connection is not None:
            connection.close()
