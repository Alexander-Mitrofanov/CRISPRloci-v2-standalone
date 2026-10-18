"""Independent repeat mapping verification against pinned reference bytes.

No Integration or scientific-tool imports. Recompute nearest ties with banded
edit distance, including undisplayed conflicting/missing annotation evidence.
"""

import hashlib
import json
from pathlib import Path

from .validation_annotations import _payload, _strict_json, _fasta

DNA = set("ACGT")
IUPAC = set("ACGTRYSWKMBDHVN")
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def same(a, b):
    return json.dumps(a, sort_keys=True, allow_nan=False) == json.dumps(b, sort_keys=True, allow_nan=False)


def rc(s):
    return s.translate(COMP)[::-1]


def distance(a, b, limit=3):
    """Banded unit-cost DP; returns limit+1 for distances outside the budget."""
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    previous = {j: j for j in range(min(len(b), limit) + 1)}
    for i, base in enumerate(a, 1):
        current = {0: i} if i <= limit else {}
        for j in range(max(1, i - limit), min(len(b), i + limit) + 1):
            current[j] = min(previous.get(j, limit + 1) + 1,
                             current.get(j - 1, limit + 1) + 1,
                             previous.get(j - 1, limit + 1) + (base != b[j - 1]))
        if min(current.values(), default=limit + 1) > limit:
            return limit + 1
        previous = current
    return min(previous.get(len(b), limit + 1), limit + 1)


def annotation(entries, field, exact):
    values = sorted({e[field] for e in entries if e[field] is not None})
    missing = sum(e[field] is None for e in entries)
    status = ("unannotated" if not values else "conflicting" if len(values) > 1
              else "incomplete" if missing else "historical_association" if exact
              else "candidate" if field == "family" else "reference_evidence")
    return {"status": status, "values": values,
            "value": values[0] if status == "historical_association" else None,
            "missing_reference_count": missing}


def validate_repeat_map_pipeline(claimed, root, expected, integration, input_root, reference, checkpoint=None):
    try:
        return _validate(claimed, root, expected, integration, input_root, reference, checkpoint)
    except (OSError, ValueError, TypeError, KeyError, AttributeError, IndexError, RecursionError, UnicodeError):
        return False


def _validate(claimed, root, expected, integration, input_root, reference, checkpoint):
    options = claimed.tool_options.get("crisprmap", {})
    max_distance = options.get("--max-distance", 3)
    top_k = options.get("--top-k", 20)
    if input_root is None or reference is None or claimed.molecule != "DNA" or integration.get("workflow") != "repeat_map":
        return False
    reference = Path(reference).absolute()
    if any(p.is_symlink() for p in (reference, *reference.parents)):
        return False
    raw = _payload(reference.parent, reference, checkpoint)
    if sha(raw) != claimed.repeat_reference_sha256:
        return False
    ref = _strict_json(raw)
    content = dict(ref)
    checksum = content.pop("content_sha256")
    if (sha((json.dumps(content, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()) != checksum
            or ref["schema"] != {"name": "crisprmap.reference", "version": "1.0.0"}
            or not 0 < len(ref["entries"]) <= 25000):
        return False
    entries = {}
    groups = {}
    for e in ref["entries"]:
        if e["id"] in entries or not 1 <= len(e["sequence"]) <= 200 or set(e["sequence"]) - IUPAC or sha(e["sequence"]) != e["sequence_sha256"]:
            return False
        entries[e["id"]] = e
        if set(e["sequence"]) <= DNA:
            groups.setdefault(min(e["sequence"], rc(e["sequence"])), []).append(e)
    source = {}
    for path in sorted(input_root.iterdir()):
        key, sequence = _fasta(_payload(input_root, path, checkpoint))
        if key in source or expected.get(key) != (len(sequence), sha(sequence)) or not 1 <= len(sequence) <= 200:
            return False
        source[key] = sequence
    if set(source) != set(expected) or not 1 <= len(source) <= 100:
        return False
    doc = _strict_json(_payload(root, root / "repeat-map/report.json", checkpoint))
    policy = doc["policy"]
    deadline = policy["deadline_seconds"]
    if (type(deadline) not in (int, float) or not 0 < deadline <= 300
            or not same(policy, {"max_distance": max_distance, "top_k": top_k, "max_comparisons": 100000,
                "deadline_seconds": deadline, "normalization": "uppercase_U_to_T",
                "max_report_bytes": 32 * 1024**2, "retrieval": "all_closest_ties_within_edit_budget",
                "ambiguous_bases": "unsupported", "similarity": "1-distance/max(query_length,reference_length)",
                "annotation_scope": "all_best_reference_occurrences_before_display_truncation",
                "orientation_scope": "relative_to_stored_reference_not_transcription_strand"})
            or doc["schema"] != {"name": "crisprmap.search", "version": "1.0.0"}
            or doc["execution_completed"] is not True
            or doc["method"] != "exhaustive_nearest_global_unit_edit_distance_both_orientations"
            or doc["software"]["name"] != "CRISPRmap-v2"
            or doc["software"]["version"] not in {"2.0.0a1", "2.1.0a1", "2.1.0a2"}
            or not same(doc["input"], {"sha256": sha("".join(f">{k}\n{v}\n" for k, v in source.items())), "record_count": len(source)})
            or not same(doc["reference"], {"namespace": ref["provenance"]["namespace"], "content_sha256": checksum,
                "file_sha256": claimed.repeat_reference_sha256, "entry_count": len(entries),
                "eligible_entry_count": sum(map(len, groups.values())), "canonical_group_count": len(groups),
                "provenance": ref["provenance"]}) or len(doc["results"]) != len(source)):
        return False
    cache = {}
    comparisons = 0
    upper = sum(sum(1 if t == rc(t) else 2 for t in groups if abs(len(s) - len(t)) <= max_distance)
                for s in set(source.values()) if set(s) <= DNA and min(s, rc(s)) not in groups)
    if upper > 100000:
        return False
    for row, (key, sequence) in zip(doc["results"], source.items()):
        if checkpoint:
            checkpoint()
        if (row["id"] != key or row["description"] != "" or row["sequence"] != sequence
                or row["sequence_sha256"] != sha(sequence)):
            return False
        if sequence not in cache:
            best, ties = max_distance + 1, []
            canonical = min(sequence, rc(sequence))
            if set(sequence) <= DNA:
                if canonical in groups:
                    best, ties = 0, groups[canonical]
                else:
                    for target in sorted(groups):
                        if checkpoint:
                            checkpoint()
                        if abs(len(sequence) - len(target)) > min(max_distance, best):
                            continue
                        comparisons += 1 if target == rc(target) else 2
                        d = min(distance(sequence, target, max_distance), distance(sequence, rc(target), max_distance))
                        if d <= max_distance and d <= best:
                            ties = groups[target][:] if d < best else ties + groups[target]
                            best = d
            cache[sequence] = best, sorted(ties, key=lambda e: e["id"])
        best, ties = cache[sequence]
        status = "unsupported" if set(sequence) - DNA else "exact" if best == 0 else "near" if ties else "no_hit"
        reason = "ambiguous_bases" if status == "unsupported" else None if ties else "no_reference_within_edit_budget" if groups else "no_unambiguous_reference_sequences"
        if not same({k: row[k] for k in ("status", "reason", "best_distance", "best_hit_count", "reported_hit_count", "hits_truncated")},
                    {"status": status, "reason": reason, "best_distance": best if ties else None,
                     "best_hit_count": len(ties), "reported_hit_count": min(top_k, len(ties)), "hits_truncated": len(ties) > top_k}):
            return False
        annotations = {f: annotation(ties, f, best == 0) if ties else
                       {"status": "not_assessed", "values": [], "value": None, "missing_reference_count": 0} for f in ("family", "motif")}
        if not same(row["annotations"], annotations) or len(row["hits"]) != min(top_k, len(ties)):
            return False
        for hit, e in zip(row["hits"], ties[:top_k]):
            forward, reverse = distance(sequence, e["sequence"], max_distance), distance(sequence, rc(e["sequence"]), max_distance)
            orientation = "both" if forward == reverse == best else "forward" if forward == best else "reverse"
            target = rc(e["sequence"]) if orientation == "reverse" else e["sequence"]
            a, b = hit["query_aligned"], hit["reference_aligned"]
            if (hit["reference_id"] != e["id"] or hit["orientation"] != orientation
                    or hit["aligned_orientation"] != ("reverse" if orientation == "reverse" else "forward")
                    or type(hit["distance"]) is not int or hit["distance"] != best
                    or type(hit["similarity"]) not in (int, float) or hit["similarity"] != 1 - best / max(len(sequence), len(target))
                    or any(hit[f] != e[f] for f in ("family", "motif", "accession"))
                    or hit["alignment_tie_policy"] != "diagonal_deletion_insertion"
                    or a.replace("-", "") != sequence or b.replace("-", "") != target
                    or len(a) != len(b) or len(a) > 400 or any(x == y == "-" for x, y in zip(a, b))
                    or sum(x != y for x, y in zip(a, b)) != best):
                return False
    counts = {"query_count": len(source), "unique_sequence_count": len(cache),
              "comparison_upper_bound": upper, "distance_comparisons": comparisons,
              **{f"{s}_count": sum(r["status"] == s for r in doc["results"]) for s in ("exact", "near", "no_hit", "unsupported")}}
    stage = integration["stages"]["crisprmap"]
    return (same(doc["counts"], counts) and same(stage["counts"], counts)
            and stage["status"] == "completed" and stage["reference_sha256"] == claimed.repeat_reference_sha256)
