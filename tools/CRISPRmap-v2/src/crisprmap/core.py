"""Complete nearest-reference retrieval; similarity does not classify a repeat."""

from __future__ import annotations

from copy import deepcopy
import math
import time

from .io import InputError, LimitError, canonical, digest, integer, normalize_records
from .reference import validate_reference

DNA = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
MAX_REPORT_BYTES = 32 * 1024 * 1024


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def canonical_sequence(sequence: str) -> str:
    return min(sequence, reverse_complement(sequence))


def levenshtein(query: str, target: str) -> int:
    """Exact unit-cost global Levenshtein distance using Myers bit vectors."""
    if not query:
        return len(target)
    if not target:
        return len(query)
    # The smaller pattern bounds Python's arbitrary-precision integer work.
    if len(query) > len(target):
        query, target = target, query
    masks = {}
    for index, base in enumerate(query):
        masks[base] = masks.get(base, 0) | (1 << index)
    width = len(query)
    mask = (1 << width) - 1
    high = 1 << (width - 1)
    positive, negative, score = mask, 0, width
    for base in target:
        equal = masks.get(base, 0)
        vertical = equal | negative
        horizontal = (((equal & positive) + positive) ^ positive) | equal
        ph = negative | ~(horizontal | positive)
        mh = positive & horizontal
        score += bool(ph & high) - bool(mh & high)
        ph = (ph << 1) | 1
        mh <<= 1
        positive = (mh | ~(vertical | ph)) & mask
        negative = (ph & vertical) & mask
    return score


def global_alignment(query: str, target: str) -> dict:
    """One optimal alignment, with diagonal/deletion/insertion tie priority."""
    if query == target:
        return {"query_aligned": query, "reference_aligned": target,
                "alignment_tie_policy": "diagonal_deletion_insertion"}
    matrix = [list(range(len(target) + 1))]
    for i, q in enumerate(query, 1):
        row = [i]
        for j, t in enumerate(target, 1):
            row.append(min(matrix[i - 1][j - 1] + (q != t), matrix[i - 1][j] + 1, row[j - 1] + 1))
        matrix.append(row)
    a, b = [], []
    i, j = len(query), len(target)
    while i or j:
        if i and j and matrix[i][j] == matrix[i - 1][j - 1] + (query[i - 1] != target[j - 1]):
            a.append(query[i - 1])
            b.append(target[j - 1])
            i -= 1
            j -= 1
        elif i and matrix[i][j] == matrix[i - 1][j] + 1:
            a.append(query[i - 1])
            b.append("-")
            i -= 1
        else:
            a.append("-")
            b.append(target[j - 1])
            j -= 1
    return {"query_aligned": "".join(reversed(a)), "reference_aligned": "".join(reversed(b)),
            "alignment_tie_policy": "diagonal_deletion_insertion"}


def _annotation(entries: list[dict], field: str, exact: bool) -> dict:
    values = sorted({entry[field] for entry in entries if entry[field] is not None})
    missing = sum(entry[field] is None for entry in entries)
    if not values:
        status = "unannotated"
    elif len(values) > 1:
        status = "conflicting"
    elif missing:
        status = "incomplete"
    else:
        status = "historical_association" if exact else ("candidate" if field == "family" else "reference_evidence")
    return {"status": status, "values": values,
            "value": values[0] if status == "historical_association" else None,
            "missing_reference_count": missing}


def _empty(status: str, reason: str) -> dict:
    return {"status": status, "reason": reason, "best_distance": None, "best_hit_count": 0,
            "reported_hit_count": 0, "hits_truncated": False, "hits": [],
            "annotations": {field: {"status": "not_assessed", "values": [], "value": None,
                                    "missing_reference_count": 0} for field in ("family", "motif")}}


def search(records: list[dict], reference: dict, *, max_distance: int = 3, top_k: int = 20,
           max_comparisons: int = 5_000_000, deadline_seconds: float = 300) -> dict:
    """Return all closest ties within the edit budget; display at most top_k.

    Exact canonical matches use a complete lookup index. Other searches examine
    every eligible canonical group unless the length lower bound proves that it
    cannot tie or improve the incumbent distance. Both orientations are assessed.
    Annotation consensus always covers all closest occurrences, before display
    truncation. Reference-relative orientation is not a transcription-strand call.
    """
    started = time.monotonic()
    integer(max_distance, "max_distance", 0, 200)
    integer(top_k, "top_k", 1, 1000)
    integer(max_comparisons, "max_comparisons", 1, 100_000_000)
    if (type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds)
            or not 0 < deadline_seconds <= 86400):
        raise InputError("deadline_seconds must be finite and in (0, 86400]")
    records = normalize_records(records)
    if len(records) > 10_000:
        raise LimitError("Search supports at most 10000 query occurrences")
    reference = validate_reference(reference)
    groups = {}
    eligible = 0
    for entry in reference["entries"]:
        if set(entry["sequence"]) <= DNA:
            groups.setdefault(canonical_sequence(entry["sequence"]), []).append(entry)
            eligible += 1
    ordered_groups = sorted(groups)
    reverse_groups = {target: reverse_complement(target) for target in ordered_groups}
    comparisons = 0
    # Reserve the worst-case work before starting any expensive query. Exact
    # index hits and repeated normalized queries do not consume distance calls.
    group_lengths = {}
    for target in ordered_groups:
        cost = 1 if target == reverse_groups[target] else 2
        group_lengths[len(target)] = group_lengths.get(len(target), 0) + cost
    comparison_upper_bound = sum(
        sum(count for length, count in group_lengths.items() if abs(len(sequence) - length) <= max_distance)
        for sequence in {row["sequence"] for row in records}
        if set(sequence) <= DNA and canonical_sequence(sequence) not in groups
    )
    if comparison_upper_bound > max_comparisons:
        raise LimitError("Search exceeds conservative comparison budget before execution")

    def check_deadline():
        if time.monotonic() - started > deadline_seconds:
            raise LimitError("Search exceeded deadline; no complete report exists")

    def distance(query, target):
        nonlocal comparisons
        comparisons += 1
        if comparisons > max_comparisons:
            raise LimitError("Search exceeded comparison budget; no complete report exists")
        check_deadline()
        return levenshtein(query, target)

    def analyze(sequence):
        if not set(sequence) <= DNA:
            return _empty("unsupported", "ambiguous_bases")
        if not groups:
            return _empty("no_hit", "no_unambiguous_reference_sequences")
        key = canonical_sequence(sequence)
        closest = {}
        if key in groups:
            # A zero-distance reference cannot be beaten, and the index includes
            # every exact/RC occurrence; no other group can tie distance zero.
            closest[key] = (0 if key == sequence else None, 0 if reverse_groups[key] == sequence else None)
            best = 0
        else:
            best = max_distance + 1
            for target in ordered_groups:
                check_deadline()
                if abs(len(sequence) - len(target)) > min(max_distance, best):
                    continue
                forward = distance(sequence, target)
                rc = reverse_groups[target]
                reverse = forward if rc == target else distance(sequence, rc)
                current = min(forward, reverse)
                if current > max_distance or current > best:
                    continue
                if current < best:
                    closest.clear()
                    best = current
                closest[target] = (forward if forward == best else None, reverse if reverse == best else None)
        if not closest:
            return _empty("no_hit", "no_reference_within_edit_budget")
        entries = sorted((entry for key in closest for entry in groups[key]), key=lambda row: row["id"])
        hits, alignments = [], {}
        for entry in entries[:top_k]:
            check_deadline()
            key = canonical_sequence(entry["sequence"])
            forward, reverse = closest[key]
            if entry["sequence"] != key:
                forward, reverse = reverse, forward
            orientation = "both" if forward is not None and reverse is not None else (
                "forward" if forward is not None else "reverse")
            aligned_orientation = "reverse" if orientation == "reverse" else "forward"
            target = reverse_complement(entry["sequence"]) if aligned_orientation == "reverse" else entry["sequence"]
            if target not in alignments:
                alignments[target] = global_alignment(sequence, target)
            hits.append({"reference_id": entry["id"], "orientation": orientation,
                         "aligned_orientation": aligned_orientation, "distance": best,
                         "similarity": 1 - best / max(len(sequence), len(target)),
                         **{field: entry[field] for field in ("family", "motif", "accession")},
                         **alignments[target]})
        return {"status": "exact" if best == 0 else "near", "reason": None,
                "best_distance": best, "best_hit_count": len(entries), "reported_hit_count": len(hits),
                "hits_truncated": len(entries) > len(hits), "hits": hits,
                "annotations": {field: _annotation(entries, field, best == 0) for field in ("family", "motif")}}

    results, cache, result_bytes = [], {}, 4096
    for record in records:
        check_deadline()
        sequence = record["sequence"]
        if sequence not in cache:
            cache[sequence] = analyze(sequence)
        row = {**record, "sequence_sha256": digest(sequence), **cache[sequence]}
        result_bytes += len(canonical(row))
        if result_bytes > MAX_REPORT_BYTES:
            raise LimitError("Search exceeds report byte budget; no complete report exists")
        results.append(deepcopy(row))
    check_deadline()
    report = {"schema": {"name": "crisprmap.search", "version": "1.0.0"},
            "method": "exhaustive_nearest_global_unit_edit_distance_both_orientations",
            "execution_completed": True,
            "reference": {"namespace": reference["provenance"]["namespace"],
                          "content_sha256": reference["content_sha256"],
                          "entry_count": len(reference["entries"]), "eligible_entry_count": eligible,
                          "canonical_group_count": len(groups)},
            "policy": {"max_distance": max_distance, "top_k": top_k, "max_comparisons": max_comparisons,
                       "deadline_seconds": deadline_seconds, "normalization": "uppercase_U_to_T",
                       "max_report_bytes": MAX_REPORT_BYTES,
                       "retrieval": "all_closest_ties_within_edit_budget", "ambiguous_bases": "unsupported",
                       "similarity": "1-distance/max(query_length,reference_length)",
                       "annotation_scope": "all_best_reference_occurrences_before_display_truncation",
                       "orientation_scope": "relative_to_stored_reference_not_transcription_strand"},
            "counts": {"query_count": len(records), "unique_sequence_count": len(cache),
                       "distance_comparisons": comparisons,
                       "comparison_upper_bound": comparison_upper_bound,
                       **{f"{status}_count": sum(r["status"] == status for r in results)
                          for status in ("exact", "near", "no_hit", "unsupported")}},
            "results": results}
    if len(canonical(report)) > MAX_REPORT_BYTES:
        raise LimitError("Search exceeds report byte budget; no complete report exists")
    return report
