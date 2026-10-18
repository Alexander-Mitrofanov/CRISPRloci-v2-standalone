"""Coupled actual-member views of a boundary-uncertain reconstruction.

This module knows neither model scoring nor report construction. A descriptor
summarizes admitted *whole* geometries, never combinations of their end sets.
Full source/retained-membership replay belongs to the transaction validator.
"""
import hashlib
import json
import math

from components.occurrence_protection import checked_geometry, occurrence_relation


METHOD = "source_backed_boundary_uncertain_reconstruction"
RULE = "shared-start-coupled-boundary-family-v1"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def geometry_id(geometry):
    return hashlib.sha256(json.dumps(geometry, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def negative_signature(evidence):
    """Applicability/outcome agreement, not equality of numeric measurements."""
    names = ("likely_tandem", "coding_tandem", "known_negative_annotations")
    if not all(isinstance(evidence.get(name), dict) for name in names):
        raise ValueError("boundary family requires explicit negative attestations")
    statuses = tuple(evidence[name].get("status") for name in names)
    if any(status not in ("passed", "inapplicable") for status in statuses):
        raise ValueError("boundary family negative assessment is unresolved")
    return statuses


def family_occurrence_evidence(references, geometries):
    """Check each actual geometry directly against all original S1 witnesses.

    Callers validate the original source-bound references before passing them.
    Aliases do not add votes or duplicate checks. No witness is remapped through
    a displayed member, and only compact witness identities enter this ledger.
    """
    references = {(ref["representation_sha256"], ref["source_span_sha256"], ref["score"]): ref
                  for ref in references if ref is not None}
    geometries = {geometry_id(geometry): geometry for geometry in geometries}
    witnesses = [{"representation_sha256": key[0], "source_span_sha256": key[1], "score": key[2]}
                 for key in sorted(references)]
    conflicts = []
    for gid in sorted(geometries):
        for key in sorted(references):
            status = occurrence_relation(references[key]["repeat_intervals"], geometries[gid])["status"]
            if status != "unique":
                conflicts.append({"geometry_id": gid, "witness": {
                    "representation_sha256": key[0], "source_span_sha256": key[1], "score": key[2]},
                    "status": status})
    return {"status": "conflict" if conflicts else "preserved",
            "ordinary_witnesses": witnesses, "geometry_ids": sorted(geometries),
            "direct_checks": len(references) * len(geometries), "conflicts": conflicts}


def describe_geometry_set(geometries):
    """Neutral coupled geometry arithmetic, without accepted-core authority."""
    geometries = [checked_geometry(geometry) for geometry in geometries]
    if not geometries:
        raise ValueError("coupled geometry set is empty")
    ids = [geometry_id(geometry) for geometry in geometries]
    if len(set(ids)) != len(ids) or ids != sorted(ids):
        raise ValueError("boundary family member geometries must be unique and sorted")
    starts = tuple(a for a, _b in geometries[0])
    if any(tuple(a for a, _b in geometry) != starts for geometry in geometries):
        raise ValueError("boundary family occurrence starts differ")
    ends = [sorted({geometry[index][1] for geometry in geometries})
            for index in range(len(starts))]
    return {"shared_starts": list(starts), "observed_end_sets": ends,
            "admissible_geometry_ids": ids, "supported_repeat_count": len(starts),
            "union_span": [starts[0], max(geometry[-1][1] for geometry in geometries)],
            "variable_repeat_indexes": [index for index, values in enumerate(ends) if len(values) > 1]}


def describe_family(source_sha256, anchor, members, display_representation_id):
    """Derive a descriptor from exact snapshots; callers verify admission."""
    if (not isinstance(source_sha256, str) or len(source_sha256) != 64
            or any(c not in "0123456789abcdef" for c in source_sha256)):
        raise ValueError("boundary family source binding is malformed")
    if not isinstance(members, list) or len(members) < 2:
        raise ValueError("boundary family requires distinct full geometries")
    geometries = [checked_geometry(m["fresh"]["repeat_intervals"]) for m in members]
    neutral = describe_geometry_set(geometries)
    ids = neutral["admissible_geometry_ids"]
    starts = neutral["shared_starts"]
    mappings = [m["relationship"]["correspondence"]["anchor_to_candidate_indexes"] for m in members]
    if any(mapping != mappings[0] for mapping in mappings):
        raise ValueError("boundary family core mappings differ")
    scores = [m["fresh"]["score"] for m in members]
    if any(type(s) not in (int, float) or not math.isfinite(s) or not 0 <= s <= 1 for s in scores):
        raise ValueError("boundary family member score is invalid")
    display = [m for m in members if m["location"]["representation_id"] == display_representation_id]
    if len(display) != 1:
        raise ValueError("boundary family display member is absent or ambiguous")
    display = display[0]
    observed_ends = neutral["observed_end_sets"]
    return {
        "rule": RULE, "status": "boundary_uncertain",
        "coordinate_semantics": "actual_admissible_member_view",
        "coordinate_system": "zero_based_half_open",
        "family_id": digest([RULE, source_sha256, anchor["repeat_intervals"], ids, mappings[0]]),
        "supported_repeat_count": len(starts), "core_repeat_count": len(anchor["repeat_intervals"]),
        "member_count": len(members), "admissible_geometry_ids": ids,
        "display_geometry_id": geometry_id(display["fresh"]["repeat_intervals"]),
        "display_representation_id": display_representation_id,
        "display_own_score": display["fresh"]["score"],
        "member_score_min": min(scores), "member_score_max": max(scores),
        "union_span": [starts[0], max(g[-1][1] for g in geometries)],
        "variable_repeat_indexes": [i for i, ends in enumerate(observed_ends) if len(ends) > 1],
        "shared_starts": list(starts), "observed_end_sets": observed_ends,
        "core_mapping": list(mappings[0]),
    }


def _basis_summary(basis):
    if not isinstance(basis, dict) or basis.get("method") != METHOD:
        return None
    # An orthogonal coding-tandem review is replayed by the report validator.
    # It must not erase family labels in any accepted-row export.
    if "accepted_family_tandem_review" in basis:
        if not isinstance(basis["accepted_family_tandem_review"], dict):
            raise ValueError("boundary family tandem wrapper is malformed")
        basis = {k: v for k, v in basis.items() if k != "accepted_family_tandem_review"}
    if set(basis) != {"method", "boundary_family", "evidence_sha256"}:
        raise ValueError("boundary family basis keys differ")
    proof = basis["boundary_family"]
    if digest(proof) != basis["evidence_sha256"]:
        raise ValueError("boundary family evidence digest differs")
    descriptor = describe_family(proof["source_sha256"], proof["anchor_original"],
                                 proof["members"], proof["display_representation_id"])
    if proof["descriptor"] != descriptor:
        raise ValueError("boundary family descriptor differs from actual members")
    summary = {k: v for k, v in descriptor.items()
               if k not in ("rule", "shared_starts", "observed_end_sets", "core_mapping")}
    summary["proof_sha256"] = basis["evidence_sha256"]
    return summary


def boundary_family_summary(record):
    """Compact row annotation after full publication validation; None if exact."""
    basis = record[3].get("classification_basis") if len(record) > 3 else None
    result = _basis_summary(basis)
    if result is not None:
        candidate = record[1]
        intervals = [[int(a), int(a) + len(r)] for a, r in
                     zip(candidate.list_repeat_starts, candidate.list_repeats)]
        if (geometry_id(intervals) != result["display_geometry_id"]
                or float(record[0]) != result["display_own_score"]):
            raise ValueError("boundary family displayed row differs")
    return result


def report_boundary_family_summary(array):
    """Validate compact canonical-array view without importing report_model."""
    result = _basis_summary(array.get("classification_basis"))
    if result is None:
        if "boundary_family" in array:
            raise ValueError("boundary family summary has no family proof")
        return None
    if array.get("boundary_family") != result:
        raise ValueError("boundary family canonical summary differs")
    intervals = [[r["source_interval"]["start"] - 1, r["source_interval"]["end"]]
                 for r in array["repeats"]]
    if (geometry_id(intervals) != result["display_geometry_id"]
            or array["certainty_score"] != result["display_own_score"]
            or array["repeat_count"] != result["supported_repeat_count"]
            or array["source_interval"]["start"] != intervals[0][0] + 1
            or array["source_interval"]["end"] != intervals[-1][1]):
        raise ValueError("boundary family canonical display differs")
    return result


def validate_boundary_family_audit(report, audit):
    """Bind published family selections to the separately emitted transaction.

    This detects disagreement between the two records. It is not an
    authenticity signature or independent replay of unrecorded donor ranking;
    coherently rewriting both complete ledgers is outside this check's scope.
    Ordinary-only callers retain their existing publication contract.
    """
    proofs = [array["classification_basis"]["boundary_family"]
              for array in report["arrays"]
              if array.get("classification_basis", {}).get("method") == METHOD]
    decisions = ([] if audit is None else [decision for decision in audit["decisions"]
                 if decision["action"] == "selected" and decision.get("boundary_uncertain")])
    if not proofs and not decisions:
        return
    if (audit is None or audit["source_sha256"] != report["source"]["sha256"]
            or audit["coordinate_system"] != "zero_based_half_open"):
        raise ValueError("boundary family transaction source is missing or differs")
    if (sorted(digest(proof["selection"]) for proof in proofs)
            != sorted(digest(decision) for decision in decisions)
            or audit["stats"]["boundary_uncertain_selected"] != len(proofs)):
        raise ValueError("boundary family transaction decisions differ from the report")
    if any(proof["source_sha256"] != audit["source_sha256"] for proof in proofs):
        raise ValueError("boundary family proof source differs from transaction")
