"""Source-bound, exact same-record spacer replay for prophage context.

The bounded replay verifies all matches, including excluded array matches. It
uses original canonical detector reports, never adapter-derived query claims.
"""

import hashlib

COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
CATEGORIES = {"Bona-fide", "Possible", "Alternative", "Possible Discarded", "Low score"}


def sha(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def interval(unit, sequence):
    pos = unit["source_interval"]
    a, b = pos["start"], pos["end"]
    if type(a) is not int or type(b) is not int or not 1 <= a <= b <= len(sequence):
        raise ValueError("Invalid native interval")
    if (
        unit.get("sequence", unit.get("source_sequence", sequence[a - 1 : b]))
        != sequence[a - 1 : b]
    ):
        raise ValueError("Native sequence differs")
    return a, b


def inventory(reports, records):
    arrays = []
    queries = []
    seen = set()
    native_ids = set()
    for report in reports:
        source = report["source"]
        key = source["id"]
        sequence = records[key]
        if (
            key in seen
            or report["schema"] != {"name": "CRISPRidentify-report", "version": "1.2.0"}
            or type(source["length"]) is not int
            or source["length"] != len(sequence)
            or source["sha256"] != sha(sequence)
            or source["sha256_scope"] != "uppercase_sequence"
        ):
            raise ValueError("Identify source differs")
        seen.add(key)
        coords = report["coordinate_system"]
        if (
            type(coords["indexing"]) is not int
            or coords["indexing"] != 1
            or coords["end_inclusive"] is not True
            or coords["reference"] != "input_forward_strand"
            or coords["orientation_is_separate"] is not True
        ):
            raise ValueError("Coordinate frame differs")
        for native in report["arrays"]:
            native_id = native["id"]
            if (
                not isinstance(native_id, str)
                or not 1 <= len(native_id) <= 500
                or any(not 33 <= ord(char) <= 126 for char in native_id)
                or (key, native_id) in native_ids
            ):
                raise ValueError("Invalid or duplicate native array identity")
            native_ids.add((key, native_id))
            a, b = interval(native, sequence)
            orientation = native["orientation"]["strand"]
            category = native["category"]
            if category not in CATEGORIES or orientation not in ("+", "-", None):
                raise ValueError("Native category/orientation")
            aid = f"array-{len(arrays) + 1:06d}"
            array = {
                "id": aid,
                "native_id": native["id"],
                "source_id": key,
                "start": a,
                "end": b,
                "category": category,
                "orientation": orientation,
            }
            arrays.append(array)
            repeats = native["repeats"]
            spacers = native["spacers"]
            deleted = native.get("spacer_deletions", [])
            if (
                len(repeats) < 2
                or len(repeats) != len(spacers) + len(deleted) + 1
                or native["repeat_count"] != len(repeats)
                or native["spacer_count"] != len(spacers) + len(deleted)
            ):
                raise ValueError("Native array architecture")
            boundaries = []
            for ordinal, repeat in enumerate(repeats, 1):
                left, right = interval(repeat, sequence)
                if (
                    type(repeat["index"]) is not int
                    or repeat["index"] != ordinal
                    or (boundaries and left <= boundaries[-1][1])
                ):
                    raise ValueError("Repeat order")
                boundaries.append((left, right))
            if boundaries[0][0] != a or boundaries[-1][1] != b:
                raise ValueError("Array extent")
            ordinals = set()
            for kind, units in (("observed", spacers), ("deletion", deleted)):
                for unit in units:
                    ordinal = unit["index"]
                    if (
                        type(ordinal) is not int
                        or not 1 <= ordinal < len(repeats)
                        or ordinal in ordinals
                    ):
                        raise ValueError("Spacer ordinal")
                    ordinals.add(ordinal)
                    left = boundaries[ordinal - 1][1] + 1
                    right = boundaries[ordinal][0] - 1
                    query = {
                        "id": f"spacer-{len(queries) + 1:06d}",
                        "kind": kind,
                        "source_id": key,
                        "array_id": aid,
                        "ordinal": ordinal,
                        "category": category,
                        "source_orientation": orientation,
                        "boundary": None,
                        "start": None,
                        "end": None,
                        "sequence": None,
                        "sha256": None,
                    }
                    if kind == "deletion":
                        if (
                            right != left - 1
                            or unit["sequence"] != ""
                            or unit["source_boundary"]
                            != {
                                "position": right,
                                "coordinate_system": "0-based interbase",
                                "reference": "input_forward_strand",
                            }
                        ):
                            raise ValueError("Deletion differs")
                        query["boundary"] = right
                        status = "deletion"
                    else:
                        if interval(unit, sequence) != (left, right):
                            raise ValueError("Spacer differs")
                        seq = sequence[left - 1 : right]
                        query.update(
                            sequence=seq, sha256=sha(seq), start=left, end=right
                        )
                        status = (
                            "category_excluded"
                            if category not in {"Bona-fide", "Possible"}
                            else "unsupported_length"
                            if not 16 <= len(seq) <= 1000
                            else "ambiguous_query"
                            if set(seq) - set("ACGT")
                            else "eligible"
                        )
                    query["status"] = status
                    queries.append(query)
                    if len(queries) > 10000:
                        raise ValueError("Query limit")
            if len(arrays) > 10000:
                raise ValueError("Array limit")
    if seen != set(records):
        raise ValueError("Missing identify source")
    return arrays, queries


def validate_search(report, reports, records, regions):
    arrays, queries = inventory(reports, records)
    policy = report["policy"]
    fixed = {
        "substitutions": 0,
        "min_length": 16,
        "max_length": 1000,
        "categories": ["Bona-fide", "Possible"],
        "scope": "same-record",
        "flank": 10,
        "array_padding": 0,
        "max_hits": 10000,
        "max_candidates": 2000000,
        "max_pairs": 100000,
        "max_evidence_bytes": 20000000,
        "max_work": 500000000,
    }
    if (
        report["schema"] != {"name": "crisprspacer.analysis", "version": "1.0.0"}
        or report["execution_status"] != "completed"
        or report["search_complete"] is not True
        or report["results_truncated"] is not False
        or report["source_targets"] is not True
        or report["retrieval"] != "exhaustive_full_length_hamming"
        or {k: policy[k] for k in fixed} != fixed
        or type(policy["timeout"]) is not int
        or not 1 <= policy["timeout"] <= 86400
        or report["queries"] != queries
    ):
        raise ValueError("Spacer query or policy differs")
    expected = []
    work = 0
    cache = {}
    for query in queries:
        if query["status"] != "eligible":
            continue
        source = query["source_id"]
        seq = records[source]
        q = query["sequence"]
        reverse = q.translate(COMP)[::-1]
        key = (source, q)
        if key not in cache:
            hits = []
            for pattern, strand in (
                [(q, "both")] if q == reverse else [(q, "+"), (reverse, "-")]
            ):
                work += len(seq)
                if work > 500000000:
                    raise ValueError("Spacer replay work limit")
                pos = seq.find(pattern)
                while pos != -1:
                    hits.append((pos + 1, pos + len(q), strand))
                    if len(hits) > 10000:
                        raise ValueError("Spacer replay hit limit")
                    pos = seq.find(pattern, pos + 1)
            cache[key] = hits
        for start, end, strand in cache[key]:
            overlap = [
                a
                for a in arrays
                if a["source_id"] == source and a["start"] <= end and start <= a["end"]
            ]
            classification = (
                "source_array"
                if any(a["id"] == query["array_id"] for a in overlap)
                else "other_array"
                if overlap
                else "outside_reported_arrays"
            )
            expected.append(
                (
                    query["id"],
                    source,
                    start,
                    end,
                    strand,
                    classification,
                    tuple(a["id"] for a in overlap),
                )
            )
            if len(expected) > 10000:
                raise ValueError("Spacer match limit")
    expected.sort(key=lambda x: (x[0], x[1], x[2], x[4]))
    if len(report["matches"]) != len(expected):
        raise ValueError("Omitted or added spacer match")
    querymap = {q["id"]: q for q in queries}
    matches = []
    for item, row in zip(expected, report["matches"], strict=True):
        qid, source, start, end, strand, classification, overlap = item
        query = querymap[qid]
        length = len(query["sequence"])
        if (
            (
                row["query_id"],
                row["target_id"],
                row["target_start"],
                row["target_end"],
                row["relative_strand"],
                row["array_overlap"],
                tuple(row["overlapping_array_ids"]),
            )
            != item
            or row["query_sha256"] != query["sha256"]
            or row["target_sha256"] != sha(records[source])
            or row["target_segments"] != [{"start": start, "end": end}]
            or row["wraps_origin"] is not False
            or row["query_start"] != 1
            or row["query_end"] != length
            or row["query_length"] != length
            or row["identities"] != length
            or row["query_coverage"] != 1.0
            or row["substitutions"] != 0
            or row["gaps"] != 0
            or row["excluded"] is not bool(overlap)
            or row["same_source_record"] is not True
            or row["pam_status"] != "not_evaluated"
            or row["evalue"] is not None
            or row["confidence"] is not None
            or row["target_sequence_forward"] != records[source][start - 1 : end]
            or row["target_sequence_query_frame"] != query["sequence"]
        ):
            raise ValueError("Spacer alignment differs")
        if not overlap:
            matches.append(
                {
                    "query_id": qid,
                    "source_id": source,
                    "array_id": query["array_id"],
                    "spacer_index": query["ordinal"],
                    "target_source_id": source,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "region_ids": [
                        r["id"]
                        for r in regions
                        if r["source_id"] == source
                        and start <= r["end"]
                        and r["start"] <= end
                    ],
                    "match_class": "sequence_match",
                    "identity": 1.0,
                }
            )
    counts = {
        "query_occurrences": len(queries),
        "eligible_queries": sum(q["status"] == "eligible" for q in queries),
        "targets": len(records),
        "alignment_occurrences": len(expected),
        "candidate_matches": len(matches),
        "excluded_array_matches": len(expected) - len(matches),
    }
    if any(
        type(report["counts"][k]) is not int or report["counts"][k] != v
        for k, v in counts.items()
    ):
        raise ValueError("Spacer count differs")
    return matches, counts


def context_document(report, reports, records, regions):
    """Reconstruct the complete source crosswalk independently of the adapter."""
    matches, counts = validate_search(report, reports, records, regions)
    arrays, queries = inventory(reports, records)
    return {
        "search_policy": {'scope': 'same-record', 'query_categories': ['Bona-fide', 'Possible'], 'query_alphabet': 'ACGT', 'min_query_length': 16, 'max_query_length': 1000, 'method': 'exhaustive_full_length_exact', 'substitutions': 0, 'gaps': 0, 'strands': 'both', 'palindromic_strand': 'both', 'array_exclusion': 'overlap_with_all_reported_arrays', 'array_padding': 0, 'search_complete': True, 'results_truncated': False, 'pam': 'not_evaluated', 'targeting_activity': 'not_assessed', 'short_source_matching': 'included_independently_of_prophage_eligibility'},
        "source_inventory_complete": True,
        "source_arrays": arrays,
        "source_spacers": queries,
        "spacer_matches": matches,
        "search_counts": counts,
    }
