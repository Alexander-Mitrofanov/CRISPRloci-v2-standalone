"""Bound, nonrecursive candidate birth origins for conservative competition.

Constructors create ordinary candidates. Factories deriving a candidate must
inherit the parent origins before emission; the terminal repair factory replaces
the provisional constructor birth with a source-bound repair birth. Copies never
create an independent ordinary witness. Legacy classification does not inspect
these observations. Invalid inheritance is recorded and rejected only when the
conservative competition policy consumes it.
"""

import hashlib
import json
from typing import NamedTuple


class RepairBirth(NamedTuple):
    source_sha256: str
    attestation_sha256: str
    geometry: tuple
    representation: tuple


class CandidateOrigins(NamedTuple):
    representation: tuple
    ordinary_birth: object = None
    repair_births: tuple = ()
    error: object = None


def _representation(candidate):
    # Tuple snapshots share immutable strings. There is no constructor-time
    # hashing, recursive history, or copying of sequence bytes.
    return (tuple(candidate.list_repeat_starts), tuple(candidate.list_repeats),
            tuple(candidate.list_repeats_gaped), tuple(candidate.list_spacers))


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def mark_ordinary_birth(candidate):
    representation = _representation(candidate)
    candidate.competition_origins = CandidateOrigins(representation, representation)


def _origins(candidate):
    origins = getattr(candidate, "competition_origins", None)
    if not isinstance(origins, CandidateOrigins):
        raise ValueError("candidate competition origin is missing or malformed")
    if origins.error:
        raise ValueError("candidate competition origin is unresolved: " + origins.error)
    if origins.representation != _representation(candidate):
        raise ValueError("candidate competition origin representation binding is stale")
    if origins.ordinary_birth is None and not origins.repair_births:
        raise ValueError("candidate competition origin has no birth witness")
    return origins


def inherit_candidate_origins(target, parents):
    """Replace a provisional birth with the independent parent-origin union.

    An ordinary witness needs only one immutable birth representation; retaining
    more equal-authority witnesses would grow history without changing ownership.
    Seed repair attestations stay distinct from the current representation.
    """
    parents = tuple(parents)
    representation = _representation(target)
    try:
        if not parents:
            raise ValueError("derived candidate has no parent")
        origins = [_origins(parent) for parent in parents]
        ordinary = [origin.ordinary_birth for origin in origins
                    if origin.ordinary_birth is not None]
        repairs = tuple(sorted({birth for origin in origins for birth in origin.repair_births},
                               key=lambda birth: (birth.source_sha256, birth.attestation_sha256)))
        target.competition_origins = CandidateOrigins(
            representation, min(ordinary) if ordinary else None, repairs)
    except (TypeError, ValueError) as error:
        target.competition_origins = CandidateOrigins(representation, error=str(error))


def mark_repair_origin(candidate, attestation):
    geometry = (tuple(candidate.list_repeat_starts), tuple(
        start + len(repeat) for start, repeat in zip(
            candidate.list_repeat_starts, candidate.list_repeats)))
    source_hash = attestation.get("source_sha256")
    if (attestation.get("status") != "generated"
            or attestation.get("policy") != "terminal-only-trim-v1"
            or attestation.get("coordinate_system") != "zero_based_half_open"
            or tuple(tuple(x) for x in attestation.get("child_geometry", ())) != geometry
            or not isinstance(source_hash, str) or len(source_hash) != 64
            or any(char not in "0123456789abcdef" for char in source_hash)):
        raise ValueError("repair birth is not bound to a generated source geometry")
    representation = _representation(candidate)
    candidate.competition_origins = CandidateOrigins(
        representation, repair_births=(RepairBirth(
            source_hash, _digest(attestation), geometry, representation),))


def candidate_origin_kind(candidate):
    origins = _origins(candidate)
    return "ordinary" if origins.ordinary_birth is not None else "repair_only"


def candidate_origin_evidence(candidate, source_sha256=None):
    origins = _origins(candidate)
    if source_sha256 is not None and any(
            birth.source_sha256 != source_sha256 for birth in origins.repair_births):
        raise ValueError("candidate repair birth source differs from the current source")
    return {
        "kind": candidate_origin_kind(candidate),
        "current_representation_sha256": _digest(origins.representation),
        "ordinary_birth_sha256": (_digest(origins.ordinary_birth)
                                    if origins.ordinary_birth is not None else None),
        "repair_births": [
            {"source_sha256": birth.source_sha256,
             "attestation_sha256": birth.attestation_sha256,
             "seed_geometry": [list(x) for x in birth.geometry],
             "seed_representation_sha256": _digest(birth.representation)}
            for birth in origins.repair_births],
    }
