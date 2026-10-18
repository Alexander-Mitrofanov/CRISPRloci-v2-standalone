"""Independent source, alignment and reference validation for CRISPRspacer.

This boundary deliberately imports neither Integration nor CRISPRspacer.
"""

from __future__ import annotations

import hashlib
import re
from urllib.parse import urlsplit

from .validation_annotations import _fasta, _payload, _strict_json

POLICY = {
    "method": "full-length-exact", "substitutions": 0, "gaps": 0,
    "strands": "both", "topology": "linear", "array_overlap": "unknown", "pam": "not_evaluated",
}
REFERENCE_FIELDS = {
    "name", "source_url", "source_release", "source_file_sha256", "database_sha256",
    "record_count", "eligible_record_count", "skipped_ambiguous", "skipped_length",
    "seed_length", "min_length", "max_length",
}
MATCH_FIELDS = {
    "spacer_id", "spacer_sequence", "spacer_sha256", "target_id", "target_sha256",
    "target_start", "target_end", "relative_strand", "identities", "query_length",
    "substitutions", "gaps", "query_coverage", "array_overlap", "pam_status", "confidence", "evalue",
}
COUNT_FIELDS = {
    "matches", "matched_spacers", "targets", "target_bases", "seed_windows",
    "ambiguous_seed_windows", "candidates_checked",
}


def _digest(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def _sha(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _text(value, limit=500):
    return isinstance(value, str) and 0 < len(value) <= limit and all(32 <= ord(c) < 127 for c in value)


def validate_protospacer_pipeline(claimed, root, expected, integration, input_root, checkpoint=None, *, database_path=None):
    if input_root is None or not _sha(claimed.spacer_database_sha256) or claimed.molecule != "DNA":
        return False
    try:
        if integration.get("workflow") != "protospacer" or integration.get("status") != "completed":
            return False
        source = {}
        for path in sorted(input_root.iterdir()):
            key, sequence = _fasta(_payload(input_root, path, checkpoint))
            if key in source or key not in expected or expected[key] != (len(sequence), _digest(sequence)):
                return False
            source[key] = sequence
        if not source or set(source) != set(expected):
            return False
        native = root / "protospacer"
        manifest = _strict_json(_payload(root, native / "manifest.json", checkpoint))
        if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != "crisprspacer.database-search" or set(manifest.get("files", {})) != {"analysis.json", "inputs.json", "matches.tsv"}:
            return False
        payloads = {}
        for name in ("analysis.json", "inputs.json", "matches.tsv"):
            data = _payload(root, native / name, checkpoint)
            if manifest["files"][name] != {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}:
                return False
            payloads[name] = data
        doc = _strict_json(payloads["analysis.json"])
        retained = _strict_json(payloads["inputs.json"])
        if (
            doc.get("schema") != {"name": "crisprspacer.database-search", "version": "1.0.0"}
            or doc.get("tool") not in ({"name": "CRISPRspacer", "version": "0.2.0"}, {"name": "CRISPRspacer", "version": "0.3.0"},
                                      {"name": "CRISPRspacer", "version": "0.3.1"})
            or doc.get("status") != "completed" or doc.get("search_complete") is not True
            or doc.get("policy") != POLICY
        ):
            return False
        reference = doc["reference"]
        if not isinstance(reference, dict) or set(reference) != REFERENCE_FIELDS:
            return False
        if reference["database_sha256"] != claimed.spacer_database_sha256 or not _sha(reference["source_file_sha256"]):
            return False
        if not all(_text(reference[name], 2000 if name == "source_url" else 500) for name in ("name", "source_url", "source_release")):
            return False
        url = urlsplit(reference["source_url"])
        if url.scheme != "https" or not url.netloc or url.username is not None or url.password is not None:
            return False
        for name in ("record_count", "eligible_record_count", "skipped_ambiguous", "skipped_length", "seed_length", "min_length", "max_length"):
            if type(reference[name]) is not int or reference[name] < 0:
                return False
        if (reference["seed_length"], reference["min_length"], reference["max_length"]) != (16, 16, 1000):
            return False
        if reference["record_count"] != reference["eligible_record_count"] + reference["skipped_ambiguous"] + reference["skipped_length"]:
            return False
        targets = doc["targets"]
        if not isinstance(targets, list) or len(targets) != len(source):
            return False
        target_ids = set()
        for item in targets:
            key = item["id"]
            if key not in source or key in target_ids or item != {"id": key, "length": len(source[key]), "sha256": _digest(source[key])}:
                return False
            target_ids.add(key)
        retained_targets = [{"id": key, "sequence": value, "sha256": _digest(value), "length": len(value), "assembly_id": None, "topology": "linear"} for key, value in source.items()]
        combined = "".join(f">{key}\n" + "".join(value[offset:offset+80] + "\n" for offset in range(0, len(value), 80))
                           for key, value in source.items())
        if retained != {"schema_version": "1.0.0", "targets": retained_targets, "target_file_sha256": _digest(combined), "reference": reference}:
            return False
        rows = doc["matches"]
        if not isinstance(rows, list) or len(rows) > 10_000:
            return False
        unique, spacers = set(), {}
        hashes = {key: _digest(value) for key, value in source.items()}
        for row in rows:
            if checkpoint:
                checkpoint()
            if not isinstance(row, dict) or set(row) != MATCH_FIELDS:
                return False
            sid, query, key = row["spacer_id"], row["spacer_sequence"], row["target_id"]
            if not _text(sid) or any(c.isspace() for c in sid) or not isinstance(query, str) or not 16 <= len(query) <= 1000 or set(query) - set("ACGT"):
                return False
            if key not in source or row["target_sha256"] != hashes[key] or row["spacer_sha256"] != _digest(query):
                return False
            if sid in spacers and spacers[sid] != query:
                return False
            spacers[sid] = query
            start, end, strand = row["target_start"], row["target_end"], row["relative_strand"]
            if type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(source[key]) or end-start+1 != len(query):
                return False
            reverse = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            if strand not in {"+", "-", "both"} or (strand == "both") != (query == reverse):
                return False
            if source[key][start-1:end] != (reverse if strand == "-" else query):
                return False
            identity = sid, key, start, end, strand
            if identity in unique:
                return False
            unique.add(identity)
            if any(type(row[name]) is not int or row[name] != len(query) for name in ("identities", "query_length")):
                return False
            if any(type(row[name]) is not int or row[name] != 0 for name in ("substitutions", "gaps")):
                return False
            if type(row["query_coverage"]) not in {int, float} or row["query_coverage"] != 1 or row["array_overlap"] != "unknown" or row["pam_status"] != "not_evaluated" or row["confidence"] is not None or row["evalue"] is not None:
                return False
        counts = doc["counts"]
        if not isinstance(counts, dict) or set(counts) != COUNT_FIELDS or any(type(value) is not int or value < 0 for value in counts.values()):
            return False
        windows = ambiguous = 0
        for sequence in source.values():
            run = 0
            windows += max(0, len(sequence)-15)
            for index, letter in enumerate(sequence):
                if checkpoint and index % 65536 == 0:
                    checkpoint()
                run = run + 1 if letter in "ACGT" else 0
                if index >= 15 and run < 16:
                    ambiguous += 1
        if counts["matches"] != len(rows) or counts["matched_spacers"] != len(spacers) or counts["targets"] != len(source) or counts["target_bases"] != sum(map(len, source.values())) or counts["seed_windows"] != windows or counts["ambiguous_seed_windows"] != ambiguous or counts["candidates_checked"] < len(rows) or len(spacers) > reference["eligible_record_count"]:
            return False
        outcome = "no_assessable_windows" if windows == ambiguous or reference["eligible_record_count"] == 0 else "matches" if rows else "no_matches"
        if doc["outcome"] != outcome:
            return False
        from .validation_protospacer_reference import validate_reference
        validate_reference(database_path, claimed.spacer_database_sha256, reference, spacers, checkpoint)
        stage = integration["stages"]["crisprspacer"]
        return stage.get("status") == "completed" and stage.get("counts") == counts and stage.get("database_sha256") == claimed.spacer_database_sha256
    except (OSError, UnicodeError, ValueError, TypeError, KeyError, AttributeError, IndexError, RecursionError):
        return False
