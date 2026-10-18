"""Optional frozen local reference retrieval; similarity never assigns a subtype."""

from __future__ import annotations

import hashlib

from .inputs import DNA, MAX_REPORT_BYTES, bounded_bytes, parse_json, require_object, reverse_complement, sha


def within_distance(query, target, maximum):
    """Exact threshold decision with O(length * maximum) work and O(maximum) memory.

    Any alignment with at most maximum edits remains inside this diagonal band.
    Full traceback is only needed for references which pass this decision.
    """
    if abs(len(query) - len(target)) > maximum:
        return False
    previous = {j: j for j in range(min(len(target), maximum) + 1)}
    for i, q in enumerate(query, 1):
        row = {}
        for j in range(max(0, i - maximum), min(len(target), i + maximum) + 1):
            row[j] = min(
                previous.get(j - 1, maximum + 1) + (q != target[j - 1]) if j else maximum + 1,
                previous.get(j, maximum + 1) + 1,
                row.get(j - 1, maximum + 1) + 1,
            )
        if not row or min(row.values()) > maximum:
            return False
        previous = row
    return previous.get(len(target), maximum + 1) <= maximum


def align(query, target):
    """Full-length unit-cost edit alignment; diagonal, deletion, insertion tie order."""
    costs = [list(range(len(target) + 1))]
    for i, q in enumerate(query, 1):
        row = [i]
        for j, t in enumerate(target, 1):
            row.append(min(costs[i - 1][j - 1] + (q != t), costs[i - 1][j] + 1, row[j - 1] + 1))
        costs.append(row)
    i, j = len(query), len(target)
    a, b = [], []
    while i or j:
        if i and j and costs[i][j] == costs[i - 1][j - 1] + (query[i - 1] != target[j - 1]):
            a.append(query[i - 1])
            b.append(target[j - 1])
            i -= 1
            j -= 1
        elif i and costs[i][j] == costs[i - 1][j] + 1:
            a.append(query[i - 1])
            b.append("-")
            i -= 1
        else:
            a.append("-")
            b.append(target[j - 1])
            j -= 1
    a, b = "".join(reversed(a)), "".join(reversed(b))
    return {
        "edit_distance": costs[-1][-1],
        "query_aligned": a,
        "reference_aligned": b,
        "query_coverage": 1.0,
        "reference_coverage": 1.0,
        "identity": sum(x == y for x, y in zip(a, b)) / len(a),
        "mismatches": sum(x != y and x != "-" and y != "-" for x, y in zip(a, b)),
        "gap_columns": sum(x == "-" or y == "-" for x, y in zip(a, b)),
    }


class References:
    def __init__(self, path, max_distance=3):
        if type(max_distance) is not int or not 0 <= max_distance <= 3:
            raise ValueError("Reference maximum distance must be an integer between 0 and 3")
        payload = bounded_bytes(path, MAX_REPORT_BYTES)
        self.document = parse_json(payload.decode("utf-8"))
        self.sha256 = hashlib.sha256(payload).hexdigest()
        if self.document.get("schema") != {"name": "CRISPRrepeat-reference", "version": "1.0.0"}:
            raise ValueError("Unsupported reference schema")
        for key in ("release", "license", "source_url", "build_date"):
            if not isinstance(self.document.get(key), str) or not self.document[key]:
                raise ValueError("Reference provenance is incomplete")
        entries = self.document.get("entries")
        if not isinstance(entries, list) or not 1 <= len(entries) <= 1000:
            raise ValueError("Reference requires 1–1000 entries")
        seen = set()
        for entry in entries:
            require_object(entry, "reference entry")
            if (not isinstance(entry.get("id"), str) or not entry["id"]
                    or len(entry["id"]) > 200 or entry["id"] in seen):
                raise ValueError("Invalid reference identity")
            seen.add(entry["id"])
            sequence = entry.get("sequence")
            if (
                not isinstance(sequence, str)
                or not 1 <= len(sequence) <= 200
                or set(sequence) - (DNA - set("RYSWKMBDHVN"))
            ):
                raise ValueError("Reference requires unambiguous DNA")
            if (
                entry.get("sha256") != sha(sequence)
                or not isinstance(entry.get("occurrences"), list)
                or not entry["occurrences"]
            ):
                raise ValueError("Missing reference sequence binding/occurrences")
            for occurrence in entry["occurrences"]:
                require_object(occurrence, "reference occurrence")
                if not all(
                    isinstance(occurrence.get(k), str) and occurrence[k]
                    for k in ("accession", "source_url", "evidence")
                ):
                    raise ValueError("Reference occurrence provenance is incomplete")
        self.entries = entries
        # Index unique oriented targets once, retaining every reference origin.
        self.targets = {}
        for entry in entries:
            for orientation, target in (("forward", entry["sequence"]),
                                        ("reverse", reverse_complement(entry["sequence"]))):
                self.targets.setdefault(target, []).append((entry, orientation))
        self.max_distance = max_distance
        self.cache = {}

    def search(self, sequence):
        if not isinstance(sequence, str) or not 1 <= len(sequence) <= 200:
            raise ValueError("Reference query must contain 1–200 nucleotides")
        sequence = sequence.replace("U", "T")
        if sequence in self.cache:
            return self.cache[sequence]
        if set(sequence) - set("ACGT"):
            return {"status": "unsupported", "reason": "ambiguous_bases", "hits": []}
        hits = []
        for target, origins in self.targets.items():
            if not within_distance(sequence, target, self.max_distance):
                continue
            result = align(sequence, target)
            for entry, orientation in origins:
                hits.append(
                    {
                        "reference_id": entry["id"],
                        "orientation": orientation,
                        "match": "exact" if result["edit_distance"] == 0 else "near",
                        "occurrences": entry["occurrences"],
                        **result,
                    }
                )
        hits.sort(key=lambda hit: (hit["edit_distance"], hit["reference_id"], hit["orientation"]))
        result = {
            "status": "truncated" if len(hits) > 50 else ("completed" if hits else "no_hit"),
            "release": self.document["release"],
            "reference_sha256": self.sha256,
            "method": "global_unit_edit_distance_both_orientations",
            "max_distance": self.max_distance,
            "total_hit_count": len(hits),
            "reported_hit_count": min(len(hits), 50),
            "hits": hits[:50],
        }
        self.cache[sequence] = result
        return result
