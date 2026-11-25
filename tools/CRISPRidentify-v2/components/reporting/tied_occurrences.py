"""Source-bound reporting of exact occurrences shared by unselected reconstruction peers."""

import hashlib
import json
import math

from components.reconstruction.reconstruction_boundary_family import geometry_id
from components.reconstruction.reconstruction_transaction import EVIDENCE_FILENAME

EXTENSION = "tied_reconstruction_occurrences"
COORDINATES = "zero_based_half_open"
AUDIT_SCHEMA = {"name": "CRISPRidentify-reconstruction-evidence", "version": "1.0.0"}


def audit_text(audit):
    """The native sidecar writer and its digest must cover exactly the same bytes."""
    return (
        json.dumps(audit, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    )


def _require(condition, message):
    if not condition:
        raise ValueError("Tied reconstruction occurrences: " + message)


def _geometry(location, category, source_length):
    _require(location["category_index"] == category, "invalid native category")
    intervals = []
    for interval in location["repeat_intervals"]:
        _require(
            isinstance(interval, (tuple, list)) and len(interval) == 2,
            "invalid repeat interval",
        )
        start, end = interval
        _require(
            type(start) is int
            and type(end) is int
            and 0 <= start < end <= source_length,
            "repeat outside source",
        )
        _require(
            not intervals or start >= intervals[-1][1],
            "unordered or overlapping repeats",
        )
        intervals.append((start, end))
    _require(
        len(intervals) >= 2 and len(intervals) == location["repeat_count"],
        "repeat count mismatch",
    )
    _require(
        list(location["source_span"]) == [intervals[0][0], intervals[-1][1]],
        "invalid source span",
    )
    _require(
        geometry_id(intervals) == location["geometry_id"],
        "geometry ID disagrees with source intervals",
    )
    return intervals


def _score(value):
    _require(
        type(value) in (int, float) and math.isfinite(value) and 0 <= value <= 1,
        "invalid score",
    )
    return value


def _accepted_anchor(arrays, anchor, core, source):
    matches = []
    for array in arrays:
        if (
            array["category"] != "Bona-fide"
            or array["legacy_interval_key"] != anchor["locus_key"]
        ):
            continue
        if array["certainty_score"] != anchor["score"]:
            continue
        path = [
            (r["source_interval"]["start"] - 1, r["source_interval"]["end"])
            for r in array["repeats"]
        ]
        if path == core:
            _require(
                all(
                    r["sequence"] == source[a:b]
                    for r, (a, b) in zip(array["repeats"], core)
                ),
                "published anchor repeats differ from source",
            )
            matches.append(array["id"])
    _require(len(matches) <= 1, "ambiguous canonical anchor binding")
    return matches[0] if matches else None


def _observations(index, decision, audit, report, source):
    ev = decision["evidence"]
    policy = audit["policy"]
    _require(
        decision["rule"] == "conservative-core-reconstruction"
        and decision["coordinate_system"] == COORDINATES
        and decision["selected"] is None
        and ev["policy"] == policy,
        "invalid native decision",
    )
    _require(
        ev["direct_peer_validation_complete"] is True
        and ev["direct_peer_rows"] == ev["freshly_scored_peer_rows"]
        and ev["direct_peer_rows"] >= len(ev["same_rank_peers"]),
        "incomplete direct-peer validation",
    )
    _require(
        all(
            ev["anchor_validation"][k] is True
            for k in (
                "source_valid",
                "structural_valid",
                "negatives_passed",
                "score_fresh",
                "qualified",
            )
        )
        and _score(ev["anchor_score"]) >= policy["minimum_anchor_score"],
        "anchor failed validation",
    )
    anchor = decision["anchor"]
    core = _geometry(anchor, 0, len(source))
    peer_ids = decision["tied_geometry_ids"]
    _require(
        len(peer_ids) >= 2 and peer_ids == sorted(set(peer_ids)),
        "invalid tied geometry set",
    )
    best = _score(ev["best_reconstruction_score"])
    _require(
        best >= policy["minimum_candidate_score"]
        and ev["best_repeat_count"] > len(core),
        "ranked peers are ineligible",
    )
    paths = {}
    scores = []
    for peer in ev["same_rank_peers"]:
        _require(peer["vetoes"] == [], "vetoed peer in unresolved best-geometry set")
        score = _score(peer["fresh_score"])
        _require(
            score >= policy["minimum_candidate_score"]
            and abs(best - score) <= policy["score_tie_tolerance"],
            "peer score is outside tied rank",
        )
        scores.append(score)
        location = peer["candidate"]
        path = _geometry(location, 1, len(source))
        _require(
            len(path) == ev["best_repeat_count"],
            "peer repeat count is outside tied rank",
        )
        paths[location["geometry_id"]] = set(path)
    _require(
        sorted(paths) == peer_ids
        and len(paths) == ev["same_rank_peer_geometries"]
        and best == max(scores),
        "incomplete or inconsistent native tie set",
    )

    common = sorted(set.intersection(*paths.values()))
    outside = [
        interval
        for interval in common
        if all(interval[1] <= a or interval[0] >= b for a, b in core)
    ]
    if not outside:
        return None
    published_id = _accepted_anchor(report["arrays"], anchor, core, source)
    if published_id is None:
        # Later filtering or projection may have removed this previous anchor.
        return None
    return {
        "anchor_array_id": published_id,
        "anchor_geometry_id": anchor["geometry_id"],
        "decision_index": index,
        "tied_geometry_ids": list(peer_ids),
        "common_repeat_intervals": [list(x) for x in common],
        "additional_repeat_intervals": [list(x) for x in outside],
    }


def tied_reconstruction_occurrences(report, audit, source):
    """Build or replay the optional diagnostic from the published anchor, native audit and FASTA."""
    if audit is None:
        return None
    ties = [
        (i, d)
        for i, d in enumerate(audit["decisions"])
        if d["action"] == "abstained" and d["reason"] == "ambiguous_best_geometry"
    ]
    if not ties:
        return None
    meta = report["source"]
    _require(
        audit["schema"] == AUDIT_SCHEMA and audit["coordinate_system"] == COORDINATES,
        "unknown audit format",
    )
    _require(
        len(source) == meta["length"]
        and audit["source_sha256"] == meta["sha256"]
        and hashlib.sha256(source.encode("ascii")).hexdigest() == meta["sha256"],
        "source binding differs",
    )
    anchors = []
    seen = set()
    for index, decision in ties:
        row = _observations(index, decision, audit, report, source)
        if row is not None:
            _require(row["anchor_array_id"] not in seen, "repeated canonical anchor")
            seen.add(row["anchor_array_id"])
            anchors.append(row)
    if not anchors:
        return None
    return {
        "schema": {
            "name": "CRISPRidentify-tied-reconstruction-occurrences",
            "version": "1.0.0",
        },
        "coordinate_system": COORDINATES,
        "source": {key: meta[key] for key in ("id", "length", "sha256")},
        "reconstruction_evidence": {
            "filename": EVIDENCE_FILENAME,
            "sha256": hashlib.sha256(audit_text(audit).encode("utf-8")).hexdigest(),
        },
        "status": "unaccepted_shared_repeat_observations",
        "anchors": anchors,
    }


def validate_tied_reconstruction_occurrences(report, audit, source):
    expected = tied_reconstruction_occurrences(report, audit, source)
    _require(
        report.get("extensions", {}).get(EXTENSION) == expected,
        "published extension failed replay",
    )
