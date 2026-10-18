"""Pure reconstruction selection anchored to a verified accepted core.

Coordinates are zero-based, half-open. Search keys never establish same-array
membership. Every proposal directly maps all core repeats, and every newly
included repeat has independent local sequence support. This planner does not
mutate categories, inherit a core's score, synthesize repeat paths, or perform
external scoring. The caller validates/scores candidates in a batch and commits
the returned decisions transactionally.
"""

from bisect import bisect_left
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
import hashlib
import json
import math
from numbers import Integral, Real

from components.initial_competition import protected_ordinary_reference, competition_evidence
from components.occurrence_protection import occurrence_relation
from components.reconstruction_boundary_family import (
    negative_signature, family_occurrence_evidence, boundary_family_summary,
)
from components.reconstruction_reference import (
    MAPPED_CORE_REFERENCE, assess_added_repeat_reference, check_reference_rule,
)


RULE_NAME = "conservative-core-reconstruction"
RULE_VERSION = "1.0.0"


@dataclass(frozen=True)
class ReconstructionPolicy:
    minimum_anchor_score: float = 0.75
    minimum_candidate_score: float = 0.60
    boundary_tolerance: int = 3
    minimum_reciprocal_overlap: float = 0.80
    maximum_spacer_length: int = 140
    minimum_added_repeat_identity: float = 0.80
    score_tie_tolerance: float = 1e-12

    def __post_init__(self):
        for name in ("minimum_anchor_score", "minimum_candidate_score",
                     "minimum_reciprocal_overlap", "minimum_added_repeat_identity"):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not math.isfinite(value) or not 0 < value <= 1):
                raise ValueError("{} must be in (0, 1]".format(name))
        if self.minimum_candidate_score > self.minimum_anchor_score:
            raise ValueError("candidate floor cannot exceed anchor floor")
        for name in ("boundary_tolerance", "maximum_spacer_length"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError("{} must be a non-negative integer".format(name))
        if (isinstance(self.score_tie_tolerance, bool)
                or not isinstance(self.score_tie_tolerance, Real)
                or not math.isfinite(self.score_tie_tolerance)
                or self.score_tie_tolerance < 0):
            raise ValueError("score tie tolerance must be finite and non-negative")

    def as_dict(self):
        return dict(vars(self))


def _score(row):
    if not isinstance(row, (list, tuple)) or len(row) not in (3, 4):
        raise ValueError("reconstruction rows need score, candidate and features")
    value = row[0]
    if (isinstance(value, bool) or not isinstance(value, Real)
            or not math.isfinite(value) or not 0 <= value <= 1):
        raise ValueError("reconstruction scores must be finite in [0, 1]")
    return float(value)


def candidate_geometry(candidate):
    starts = tuple(candidate.list_repeat_starts)
    repeats = tuple(candidate.list_repeats)
    if not starts or len(starts) != len(repeats):
        raise ValueError("repeat geometry must have matching starts and sequences")
    geometry = []
    for start, repeat in zip(starts, repeats):
        if (isinstance(start, bool) or not isinstance(start, Integral) or start < 0
                or not isinstance(repeat, str) or not repeat):
            raise ValueError("repeat geometry contains invalid coordinates/sequence")
        end = int(start) + len(repeat)
        if geometry and start < geometry[-1][1]:
            raise ValueError("repeat intervals must be ordered and nonoverlapping")
        geometry.append((int(start), end))
    return tuple(geometry)


def _digest(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def candidate_geometry_id(candidate):
    """Component-local geometry ID; caller binds it to the source DNA hash."""
    return _digest(candidate_geometry(candidate))


def _candidate_signature(candidate):
    return (tuple(candidate.list_repeat_starts), tuple(candidate.list_repeats),
            tuple(candidate.list_repeats_gaped), tuple(candidate.list_spacers),
            getattr(candidate, "consensus", None),
            getattr(candidate, "consensus_gaped", None))


@dataclass(frozen=True)
class RowLocation:
    category_index: int
    locus_key: object
    row_index: int
    row: object = field(repr=False, compare=False)
    geometry: tuple = field(init=False)
    geometry_id: str = field(init=False)
    representation_id: str = field(init=False)
    signature: tuple = field(init=False, repr=False)

    def __post_init__(self):
        score = _score(self.row)
        geometry = candidate_geometry(self.row[1])
        signature = _candidate_signature(self.row[1])
        object.__setattr__(self, "geometry", geometry)
        object.__setattr__(self, "geometry_id", _digest(geometry))
        object.__setattr__(self, "signature", signature)
        object.__setattr__(self, "representation_id", _digest((signature, score)))

    @property
    def span(self):
        return self.geometry[0][0], self.geometry[-1][1]

    def as_dict(self):
        return {"category_index": self.category_index, "locus_key": self.locus_key,
                "row_index": self.row_index, "score": _score(self.row),
                "repeat_count": len(self.geometry), "source_span": list(self.span),
                "repeat_intervals": [list(interval) for interval in self.geometry],
                "geometry_id": self.geometry_id,
                "representation_id": self.representation_id}


@dataclass(frozen=True)
class CandidateValidation:
    """Fresh row plus explicit validation attestations for that exact candidate.

    ``None`` is permitted for an unknown gate, and always fails qualification.
    ``negatives_passed`` attests all applicable checks, not merely absence of a
    negative tag. Evidence must describe applicability and each check's outcome.
    """

    row: object
    source_valid: object
    structural_valid: object
    negatives_passed: object
    score_fresh: object
    evidence: object = field(default_factory=dict)

    def __post_init__(self):
        _score(self.row)
        for name in ("source_valid", "structural_valid", "negatives_passed", "score_fresh"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, bool):
                raise TypeError("{} must be boolean or None".format(name))
        if not isinstance(self.evidence, Mapping):
            raise TypeError("validation evidence must be a mapping")
        # Validate serializability at the callback boundary, before any commit.
        json.dumps(dict(self.evidence), allow_nan=False)

    @property
    def qualified(self):
        return all(value is True for value in (
            self.source_valid, self.structural_valid, self.negatives_passed,
            self.score_fresh))

    def as_dict(self):
        return {"source_valid": self.source_valid,
                "structural_valid": self.structural_valid,
                "negatives_passed": self.negatives_passed,
                "score_fresh": self.score_fresh,
                "qualified": self.qualified,
                "evidence": json.loads(json.dumps(dict(self.evidence), allow_nan=False))}


@dataclass(frozen=True)
class RepeatCorrespondence:
    status: str
    mapping: tuple = ()
    maximum_boundary_shift: int = 0
    minimum_reciprocal_overlap: float = 0.0

    def as_dict(self):
        return {"status": self.status, "anchor_to_candidate_indexes": list(self.mapping),
                "maximum_boundary_shift": self.maximum_boundary_shift,
                "minimum_reciprocal_overlap": self.minimum_reciprocal_overlap}


def _overlap(first, second):
    return max(0, min(first[1], second[1]) - max(first[0], second[0]))


def monotone_repeat_mapping(anchor_geometry, candidate_geometry_, policy=None):
    """Find a unique ordered injective mapping, allowing internal insertions.

    Counts are capped at two; dynamic programming avoids exponential enumeration
    and retains ambiguity even when several earlier mappings converge. The
    function accepts already validated ordered, nonoverlapping interval tuples.
    """
    policy = policy or ReconstructionPolicy()
    if not anchor_geometry or len(candidate_geometry_) < len(anchor_geometry):
        return RepeatCorrespondence("no_mapping")
    starts = [item[0] for item in candidate_geometry_]
    layers = []
    for anchor in anchor_geometry:
        begin = bisect_left(starts, anchor[0] - policy.boundary_tolerance)
        end = bisect_left(starts, anchor[0] + policy.boundary_tolerance + 1)
        options = []
        for index in range(begin, end):
            candidate = candidate_geometry_[index]
            shift = max(abs(anchor[0] - candidate[0]), abs(anchor[1] - candidate[1]))
            overlap = _overlap(anchor, candidate)
            reciprocal = min(overlap / (anchor[1] - anchor[0]),
                             overlap / (candidate[1] - candidate[0]))
            if shift <= policy.boundary_tolerance and reciprocal >= policy.minimum_reciprocal_overlap:
                options.append(index)
        if not layers:
            current = {index: (1, None) for index in options}
        else:
            previous = layers[-1]
            previous_indexes = sorted(previous)
            position = 0
            count = 0
            predecessor = None
            current = {}
            for index in options:
                while position < len(previous_indexes) and previous_indexes[position] < index:
                    previous_index = previous_indexes[position]
                    ways = previous[previous_index][0]
                    predecessor = previous_index if count == 0 and ways == 1 else None
                    count = min(2, count + ways)
                    position += 1
                if count:
                    current[index] = (count, predecessor if count == 1 else None)
        if not current:
            return RepeatCorrespondence("no_mapping")
        layers.append(current)
    if sum(value[0] for value in layers[-1].values()) != 1:
        return RepeatCorrespondence("ambiguous_mapping")
    index = next(iter(layers[-1]))
    mapping = []
    for layer in reversed(layers):
        mapping.append(index)
        index = layer[index][1]
    mapping.reverse()
    shifts, overlaps = [], []
    for anchor, index in zip(anchor_geometry, mapping):
        candidate = candidate_geometry_[index]
        shifts.append(max(abs(anchor[0] - candidate[0]), abs(anchor[1] - candidate[1])))
        overlap = _overlap(anchor, candidate)
        overlaps.append(min(overlap / (anchor[1] - anchor[0]),
                            overlap / (candidate[1] - candidate[0])))
    return RepeatCorrespondence("unique", tuple(mapping), max(shifts), min(overlaps))


def row_locations(categories, *, _possible=False):
    """Snapshot Bona-fide, Alternative and Possible row locations only."""
    if not isinstance(categories, (tuple, list)) or len(categories) != 5:
        raise ValueError("reconstruction selection requires five categories")
    result = []
    for category_index in ((0, 1, 2, 3) if _possible else (0, 1, 2)):
        category = categories[category_index]
        if not isinstance(category, dict):
            raise ValueError("each category must be a dictionary")
        for key, rows in category.items():
            for row_index, row in enumerate(rows):
                result.append(RowLocation(category_index, key, row_index, row))
    return result


def _relationship(anchor, donor, policy, reference_rule):
    if len(donor.geometry) <= len(anchor.geometry) or not _overlap(anchor.span, donor.span):
        return None, "not_longer_direct_overlap"
    mapping = monotone_repeat_mapping(anchor.geometry, donor.geometry, policy)
    if mapping.status != "unique":
        return None, mapping.status
    if any(right[0] - left[1] > policy.maximum_spacer_length
           for left, right in zip(donor.geometry, donor.geometry[1:])):
        return None, "spacer_above_continuity_limit"
    mapped = set(mapping.mapping)
    added = [index for index in range(len(donor.geometry)) if index not in mapped]
    support = assess_added_repeat_reference(
        anchor.row[1].list_repeats, getattr(anchor.row[1], "consensus", None), anchor.geometry,
        donor.row[1].list_repeats, donor.geometry, mapping.mapping, added,
        policy.minimum_added_repeat_identity, reference_rule=reference_rule)
    if not support["qualified"]:
        return None, support["reason"]
    identities = support["identities"]
    internal = [index for index in added if mapping.mapping[0] < index < mapping.mapping[-1]]
    flank = [index for index in added if index not in internal]
    kind = ("internal_and_flank_extension" if internal and flank else
            "internal_repeat_insertion" if internal else "flank_extension")
    relation = {"kind": kind, "correspondence": mapping.as_dict(),
            "added_repeat_indexes": added, "internal_added_repeat_indexes": internal,
            "flank_added_repeat_indexes": flank,
            "minimum_added_repeat_identity": min(identities),
            "added_repeat_identities": identities}
    if support["reference_evidence"] is not None:
        relation["added_repeat_reference"] = support["reference_evidence"]
    return relation, None


def relationship(anchor_candidate, donor_candidate, policy=None, *,
                 reference_rule=MAPPED_CORE_REFERENCE):
    """Replay direct geometry/local-sequence evidence without trusting scores.

    The result is a JSON-safe dictionary with ``qualified`` and ``reason``;
    qualified relationships additionally include kind, mapping and local evidence.
    Source reconstruction, structural filters and negatives remain caller gates.
    """
    check_reference_rule(reference_rule)
    anchor = RowLocation(0, None, 0, [0.0, anchor_candidate, []])
    donor = RowLocation(1, None, 0, [0.0, donor_candidate, []])
    relation, reason = _relationship(anchor, donor, policy or ReconstructionPolicy(), reference_rule)
    return {"qualified": relation is not None, "reason": reason, **(relation or {})}


def _prepare(categories, policy, reference_rule, *, _possible=False):
    check_reference_rule(reference_rule)
    locations = row_locations(categories, _possible=_possible)
    accepted = [location for location in locations if location.category_index in (0, 2)]
    anchor_index, donor_index = (2, 3) if _possible else (0, 1)
    donors = [location for location in locations if location.category_index == donor_index]
    # The Possible adapter follows the completed BF transaction. Its barriers
    # include every coupled family member, not just the displayed geometry.
    spans = {}
    for location in accepted:
        family = boundary_family_summary(location.row) if _possible else None
        spans[id(location)] = tuple(family["union_span"]) if family else location.span
    # A complete direct map must include an edge from the first core repeat.
    # Index those possible edges, without a source-window or candidate cap.
    starts = sorted((start, ordinal) for ordinal, donor in enumerate(donors)
                    if _score(donor.row) >= policy.minimum_candidate_score
                    for start, _end in donor.geometry) if _possible else []
    start_values = [start for start, _ordinal in starts]
    # Bind once to original scored rows. Fresh-score callbacks and category
    # movement must not turn the selected geometry into its own witness.
    references = {id(location.row): protected_ordinary_reference(location.row)
                  for location in locations}
    guarded = {}
    plans = []
    for anchor in accepted:
        if anchor.category_index != anchor_index:
            continue
        proposals = []
        rejected = Counter()
        if _score(anchor.row) >= policy.minimum_anchor_score:
            if _possible:
                start = anchor.geometry[0][0]
                left = bisect_left(start_values, start - policy.boundary_tolerance)
                right = bisect_left(start_values, start + policy.boundary_tolerance + 1)
                scoped_donors = [donors[index] for index in sorted({
                    ordinal for _start, ordinal in starts[left:right]})]
            else:
                scoped_donors = donors
            for donor in scoped_donors:
                if _score(donor.row) < policy.minimum_candidate_score:
                    continue
                relation, reason = _relationship(anchor, donor, policy, reference_rule)
                if reason:
                    if _overlap(anchor.span, donor.span):
                        rejected[reason] += 1
                    continue
                preserves = True
                for original in (anchor, donor):
                    reference = references[id(original.row)]
                    if reference is None:
                        continue
                    key = (reference["representation_sha256"], donor.geometry)
                    if key not in guarded:
                        guarded[key] = occurrence_relation(reference["repeat_intervals"],
                                                           donor.geometry)["status"] == "unique"
                    preserves &= guarded[key]
                vetoes = []
                if _possible:
                    origin_scope = competition_evidence(donor.row)
                    if origin_scope is not None and origin_scope["decision"] == "repair_scope_conflict":
                        rejected["repair_scope_conflict"] += 1
                        vetoes.append("repair_scope_conflict")
                if not preserves:
                    rejected["protected_ordinary_occurrence_conflict"] += 1
                    vetoes.append("protected_ordinary_occurrence_conflict")
                # Even a compatible second accepted core is a fusion risk;
                # shared search keys or transitive relations grant no exception.
                if any(other is not anchor and _overlap(donor.span, spans[id(other)])
                       for other in accepted):
                    rejected["overlaps_other_accepted_core"] += 1
                    vetoes.append("overlaps_other_accepted_core")
                proposals.append((donor, relation, tuple(vetoes)))
        # Without any geometrically admissible survivor there is no possible
        # selection or family. Otherwise keep vetoed direct peers for fresh
        # same-rank auditing; filtering them could manufacture uniqueness.
        if not any(not vetoes for _donor, _relation, vetoes in proposals):
            proposals = []
        plans.append((anchor, proposals, dict(rejected)))
    return locations, plans


def reconstruction_candidate_locations(categories, policy=None, *,
                                       reference_rule=MAPPED_CORE_REFERENCE, _possible=False,
                                       _prepared=None):
    """Return geometry-qualified donor locations for one batched fresh scoring.

    Anchors are not included. The validation callback must also explicitly
    certify anchor rows; callers can batch them with these donor rows if needed.
    Aliases of an identical scored row are returned only once, by object identity.
    """
    policy = policy or ReconstructionPolicy()
    _locations_, plans = (_prepare(categories, policy, reference_rule, _possible=_possible)
                          if _prepared is None else _prepared)
    seen = set()
    donors = []
    for _anchor, proposals, _rejected in plans:
        for location, _relation, _vetoes in proposals:
            if id(location.row) not in seen:
                donors.append(location)
                seen.add(id(location.row))
    return tuple(donors)


@dataclass(frozen=True)
class ReconstructionDecision:
    anchor: RowLocation
    action: str
    reason: str
    selected: object = None
    validated_row: object = field(default=None, repr=False)
    selected_locations: tuple = ()
    relation: object = None
    tied_geometry_ids: tuple = ()
    evidence: object = field(default_factory=dict)
    family_members: tuple = field(default=(), repr=False)

    def as_dict(self):
        return {"rule": RULE_NAME, "rule_version": RULE_VERSION,
                "action": self.action, "reason": self.reason,
                "coordinate_system": "zero_based_half_open",
                "anchor": self.anchor.as_dict(),
                "selected": self.selected.as_dict() if self.selected else None,
                "reconstruction_score": _score(self.validated_row) if self.validated_row else None,
                "selected_aliases": [location.as_dict() for location in self.selected_locations],
                "relationship": self.relation,
                "tied_geometry_ids": list(self.tied_geometry_ids),
                "boundary_uncertain": bool(self.family_members),
                "evidence": self.evidence}

    @property
    def collision_span(self):
        if self.family_members:
            return (min(item[0].span[0] for item in self.family_members),
                    max(item[0].span[1] for item in self.family_members))
        return self.selected.span if self.selected else None


@dataclass(frozen=True)
class ReconstructionSelectionResult:
    decisions: tuple
    stats: dict


def select_reconstructions(categories, validate_candidate, policy=None, *,
                           reference_rule=MAPPED_CORE_REFERENCE, _possible=False, _prepared=None):
    """Choose directly supported extent; leave commitment to the caller.

    Only final Alternative rows donate geometry. All other accepted arrays are
    protected. Lower-ranked and non-direct failures do not veto supported
    insertions. Direct, equally ranked peers must all pass validation; distinct
    geometries may form a coupled shared-start boundary family. Other ties
    abstain. Alignment aliases do not vote; conflicting proposals abstain
    symmetrically.
    """
    policy = policy or ReconstructionPolicy()
    if not callable(validate_candidate):
        raise TypeError("validate_candidate must be callable")
    locations, plans = (_prepare(categories, policy, reference_rule, _possible=_possible)
                        if _prepared is None else _prepared)
    validations = {}

    def validate(location):
        identity = id(location.row)
        if identity not in validations:
            result = validate_candidate(location.row)
            if not isinstance(result, CandidateValidation):
                raise TypeError("validator must return CandidateValidation")
            if result.row[1] is not location.row[1]:
                raise ValueError("validation must retain the exact candidate object")
            if _candidate_signature(location.row[1]) != location.signature:
                raise ValueError("validation mutated the candidate representation")
            validations[identity] = result
        return validations[identity]

    decisions = []
    for anchor, proposals, rejected in plans:
        evidence = {"policy": policy.as_dict(), "geometry_rejections": rejected}
        if not proposals:
            decisions.append(ReconstructionDecision(anchor, "retained", "no_supported_extension", evidence=evidence))
            continue
        core = validate(anchor)
        evidence["anchor_validation"] = core.as_dict()
        evidence["anchor_score"] = _score(core.row)
        if not core.qualified or _score(core.row) < policy.minimum_anchor_score:
            decisions.append(ReconstructionDecision(anchor, "abstained", "anchor_not_qualified", evidence=evidence))
            continue
        qualified = []
        observed = []
        validation_rejections = []
        for donor, relation, vetoes in proposals:
            validation = validate(donor)
            reasons = list(vetoes)
            if not validation.qualified or _score(validation.row) < policy.minimum_candidate_score:
                reasons.append("candidate_validation_or_score_failed")
                validation_rejections.append({"candidate": donor.as_dict(),
                                              "fresh_score": _score(validation.row),
                                              "validation": validation.as_dict()})
            if not reasons:
                qualified.append((donor, relation, validation))
            observed.append((donor, relation, validation, tuple(reasons)))
        evidence["validation_rejections"] = validation_rejections
        evidence["direct_peer_rows"] = len(observed)
        evidence["direct_peer_geometries"] = len({item[0].geometry_id for item in observed})
        evidence["direct_peer_validation_complete"] = True
        evidence["freshly_scored_peer_rows"] = len(observed)
        if not qualified:
            decisions.append(ReconstructionDecision(anchor, "abstained", "no_validated_extension", evidence=evidence))
            continue
        maximum_count = max(len(item[0].geometry) for item in qualified)
        longest = [item for item in qualified if len(item[0].geometry) == maximum_count]
        maximum_score = max(_score(item[2].row) for item in longest)
        tied = [item for item in longest if maximum_score - _score(item[2].row) <= policy.score_tie_tolerance]
        tied_geometries = tuple(sorted({item[0].geometry_id for item in tied}))
        evidence["qualified_extension_representations"] = len(qualified)
        evidence["qualified_extension_geometries"] = len({item[0].geometry_id for item in qualified})
        evidence["best_repeat_count"] = maximum_count
        evidence["best_reconstruction_score"] = maximum_score
        peers = [item for item in observed if len(item[0].geometry) == maximum_count
                 and abs(_score(item[2].row) - maximum_score) <= policy.score_tie_tolerance]
        evidence["same_rank_peers"] = [{"candidate": item[0].as_dict(),
            "fresh_score": _score(item[2].row), "vetoes": list(item[3])} for item in peers]
        evidence["same_rank_peer_geometries"] = len({item[0].geometry_id for item in peers})
        if any(item[3] for item in peers):
            decisions.append(ReconstructionDecision(anchor, "abstained", "same_rank_peer_veto",
                tied_geometry_ids=tuple(sorted({item[0].geometry_id for item in peers})), evidence=evidence))
            continue
        family = ()
        if len(tied_geometries) > 1:
            if _possible:
                decisions.append(ReconstructionDecision(anchor, "abstained", "ambiguous_possible_geometry",
                    tied_geometry_ids=tied_geometries, evidence=evidence))
                continue
            starts = {tuple(a for a, _b in item[0].geometry) for item in tied}
            mappings = {tuple(item[1]["correspondence"]["anchor_to_candidate_indexes"]) for item in tied}
            if len(starts) != 1 or len(mappings) != 1:
                decisions.append(ReconstructionDecision(anchor, "abstained", "ambiguous_best_geometry",
                                                        tied_geometry_ids=tied_geometries, evidence=evidence))
                continue
            family_guard = family_occurrence_evidence(
                [protected_ordinary_reference(original.row)
                 for original in [anchor] + [item[0] for item in tied]],
                [item[0].geometry for item in tied])
            evidence["family_occurrence_guard"] = family_guard
            if family_guard["status"] != "preserved":
                decisions.append(ReconstructionDecision(anchor, "abstained", "boundary_family_occurrence_conflict",
                                                        tied_geometry_ids=tied_geometries, evidence=evidence))
                continue
            try:
                negative_signatures = {negative_signature(item[2].evidence) for item in tied}
            except ValueError:
                negative_signatures = set()
            if len(negative_signatures) != 1:
                decisions.append(ReconstructionDecision(anchor, "abstained", "boundary_family_negative_disagreement",
                                                        tied_geometry_ids=tied_geometries, evidence=evidence))
                continue
            # Aliases do not multiply geometric support. Keep one real, best
            # scored representation per geometry; the full peer list remains
            # in evidence and all original rows are retained by the adapter.
            family = tuple(min((item for item in tied if item[0].geometry_id == gid),
                               key=lambda item: (-_score(item[2].row), item[0].representation_id))
                           for gid in tied_geometries)
        # Exact source geometry, best own score; deterministic alignment choice.
        donor, relation, validation = min(family or tied, key=lambda item:
            (-_score(item[2].row), item[0].geometry_id, item[0].representation_id)
            if family else (item[0].representation_id,))
        aliases = tuple(location for location in locations if location.category_index == (3 if _possible else 1)
                        and location.geometry_id == donor.geometry_id)
        evidence["selected_validation"] = validation.as_dict()
        evidence["selected_alignment_tie_count"] = len(tied)
        reason = "validated_shared_occurrence_extent" if family else "validated_direct_reconstruction"
        decisions.append(ReconstructionDecision(anchor, "selected", reason,
                                                donor, validation.row, aliases, relation,
                                                tied_geometries, evidence, family_members=family))

    collisions = {}
    for index, first in enumerate(decisions):
        if first.action != "selected":
            continue
        for other_index in range(index + 1, len(decisions)):
            second = decisions[other_index]
            if second.action == "selected" and _overlap(first.collision_span, second.collision_span):
                collisions.setdefault(index, []).append(second.anchor.geometry_id)
                collisions.setdefault(other_index, []).append(first.anchor.geometry_id)
    for index, conflicting in collisions.items():
        original = decisions[index]
        evidence = dict(original.evidence,
                        conflicting_anchor_geometry_ids=sorted(conflicting),
                        abstained_proposal=original.selected.as_dict(),
                        abstained_relationship=original.relation)
        decisions[index] = replace(original, action="abstained", reason="proposal_collision",
                                   selected=None, validated_row=None, selected_locations=(),
                                   relation=None, evidence=evidence, family_members=())
    stats = {"anchors": len(plans), "validation_requests": len(validations),
             "selected": sum(d.action == "selected" for d in decisions),
             "retained": sum(d.action == "retained" for d in decisions),
             "abstained": sum(d.action == "abstained" for d in decisions),
             "boundary_uncertain_selected": sum(bool(d.family_members) for d in decisions),
             "added_repeats": sum(len(d.selected.geometry) - len(d.anchor.geometry)
                                  for d in decisions if d.action == "selected"),
             "reason_counts": dict(Counter(d.reason for d in decisions))}
    return ReconstructionSelectionResult(tuple(decisions), stats)
