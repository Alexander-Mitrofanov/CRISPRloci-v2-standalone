"""Glue between scored candidates and source-backed architecture policy.

This module contains no pipeline category mutation.  It validates the scored
row contract, derives the explicitly named model features used by the rescue
policy, and separates a permissive *augmentation seed* decision from the
strict final automatic-rescue decision.  Seed eligibility can authorize only
source-backed proposal generation; it can never authorize an output category.
"""

from dataclasses import replace
import math
from numbers import Real

from components.architecture_rescue_policy import (
    AUDITED_ARCHITECTURE_RESCUE_POLICY,
    ELIGIBLE,
    INELIGIBLE,
    REVIEW,
    evaluate_architecture_rescue_eligibility,
)
from components.components_evaluated_arrays_enhancement import (
    compare_source_backed_phase_candidates,
)


SOURCE_AUGMENTATION_SATELLITE_POLICY = replace(
    AUDITED_ARCHITECTURE_RESCUE_POLICY,
    name="frozen300-source-augmentation-satellite",
    version="1.0.0",
    min_repeat_count=3,
    min_classifier_score=0.35,
    max_spacer_length_cv=0.20,
)

SOURCE_CROSS_FAMILY_PHASE_ARBITRATION_RULE = {
    "name": "source-cross-family-phase-arbitration",
    "version": "1.0.0",
    "minimum_repeat_count": 5,
    "maximum_boundary_shift": 12,
}


def _json_scalar(value):
    item = getattr(value, "item", None)
    if callable(item):
        try:
            scalar = item()
        except (TypeError, ValueError):
            pass
        else:
            if scalar is not value:
                return _json_scalar(scalar)
    return value


def _flatten_feature_vector(value):
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        value = tolist()
    while (
        isinstance(value, (list, tuple))
        and len(value) == 1
        and isinstance(value[0], (list, tuple))
    ):
        value = value[0]
    if not isinstance(value, (list, tuple)):
        raise ValueError("a scored feature vector must be a sequence")
    flattened = tuple(_json_scalar(item) for item in value)
    if any(isinstance(item, (list, tuple, dict)) for item in flattened):
        raise ValueError("a scored feature vector must be one-dimensional")
    return flattened


def named_features_from_scored_row(row, feature_name_sets):
    """Return each available, validated named feature in a score row.

    The compact runtime supplies one vector. Duplicate declarations, if a
    caller provides them, must contain numerically identical values.
    """

    if not isinstance(row, (list, tuple)) or len(row) < 3:
        raise ValueError("a scored row must contain score, candidate, features")
    try:
        raw_sets = tuple(row[2])
        name_sets = tuple(tuple(names) for names in feature_name_sets)
    except TypeError as error:
        raise ValueError("scored feature sets must be iterable") from error
    if len(raw_sets) != len(name_sets):
        raise ValueError(
            "scored feature vectors and feature-name sets have different counts"
        )

    named = {}
    for names, raw_vector in zip(name_sets, raw_sets):
        values = _flatten_feature_vector(raw_vector)
        if len(values) != len(names):
            raise ValueError(
                "a scored feature vector has the wrong number of values"
            )
        for name, raw_value in zip(names, values):
            if not isinstance(name, str) or not name:
                raise ValueError("feature names must be non-empty strings")
            value = _json_scalar(raw_value)
            if value is None:
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(float(value))
            ):
                raise ValueError(
                    "feature {!r} must be a finite number".format(name)
                )
            value = float(value)
            previous = named.get(name)
            if previous is not None and not math.isclose(
                previous,
                value,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise ValueError(
                    "feature {!r} differs between model vectors".format(name)
                )
            named[name] = value
    return named


def evaluate_scored_architecture_candidate(
    genome_sequence,
    row,
    feature_name_sets,
    structural_filter,
    *,
    policy=AUDITED_ARCHITECTURE_RESCUE_POLICY,
    coding_evidence=None
):
    """Evaluate one scored row with explicit optional negative evidence."""

    if not isinstance(row, (list, tuple)) or len(row) < 2:
        raise ValueError("architecture rescue received a malformed score row")
    score, candidate = row[:2]
    named = named_features_from_scored_row(row, feature_name_sets)
    required = ("repeat_similarity", "spacer_similarity")
    missing = tuple(name for name in required if name not in named)
    if missing:
        raise ValueError(
            "architecture rescue features are incomplete: {}".format(
                ", ".join(missing)
            )
        )

    optional = {}
    if coding_evidence is not None:
        optional.update({
            "coding_overlap_fraction": coding_evidence[0],
            "context_cds_score": coding_evidence[1],
        })

    return evaluate_architecture_rescue_eligibility(
        genome_sequence,
        candidate,
        classifier_score=score,
        repeat_similarity=named["repeat_similarity"],
        spacer_similarity=named["spacer_similarity"],
        afsf_pass=bool(structural_filter(candidate)),
        policy=policy,
        **optional,
    )


def select_augmentation_seed_ids(
    genome_sequence,
    scored_rows,
    feature_name_sets,
    structural_filter,
):
    """Return source-augmentation seed IDs and their audited decisions.

    Five-or-more-repeat primary paths must already pass the strict structural
    policy at score >=0.40.  Three- or four-repeat satellites may use the
    audited score >=0.35 / spacer-CV <=0.20 construction gate.  Satellites can
    only contribute source anchors; every final output is freshly scored and
    judged again by the strict policy.  Review and rejected decisions do not
    seed because their phase has not been independently affirmed here.
    """

    selected = set()
    decisions = {}
    roles = {}
    for rows in scored_rows.values():
        for row in rows:
            repeat_count = len(row[1].list_repeats)
            if repeat_count >= 5:
                policy = AUDITED_ARCHITECTURE_RESCUE_POLICY
                role = "primary"
            elif repeat_count in (3, 4):
                policy = SOURCE_AUGMENTATION_SATELLITE_POLICY
                role = "satellite"
            else:
                continue
            decision = evaluate_scored_architecture_candidate(
                genome_sequence,
                row,
                feature_name_sets,
                structural_filter,
                policy=policy,

            )
            candidate_id = id(row[1])
            decisions[candidate_id] = decision
            if decision.decision == ELIGIBLE:
                selected.add(candidate_id)
                roles[candidate_id] = role
            elif decision.decision not in (INELIGIBLE, REVIEW):
                raise RuntimeError("unknown architecture seed decision")
    return frozenset(selected), decisions, roles


def candidate_repeat_path(candidate):
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts,
            candidate.list_repeats,
        )
    )


def path_is_contiguous_subpath(shorter, longer):
    """Return whether one exact source-repeat path occurs within another."""

    shorter = tuple(shorter)
    longer = tuple(longer)
    if not shorter or len(shorter) > len(longer):
        return False
    return any(
        longer[offset:offset + len(shorter)] == shorter
        for offset in range(len(longer) - len(shorter) + 1)
    )


def _path_interval(path):
    return path[0][0], path[-1][1]


def _intervals_overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def lineage_candidate_families(lineage):
    """Return deterministic connected ancestry families from path lineage."""

    entries = tuple(lineage)
    by_path_id = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("lineage entries must be mappings")
        path_id = entry.get("path_id")
        if not isinstance(path_id, str) or not path_id:
            raise ValueError("every lineage entry needs a non-empty path ID")
        if path_id in by_path_id:
            raise ValueError("lineage path IDs must be unique")
        by_path_id[path_id] = entry

    adjacency = {path_id: set() for path_id in by_path_id}
    for path_id, entry in by_path_id.items():
        try:
            parent_ids = tuple(entry.get("parent_path_ids", ()))
        except TypeError as error:
            raise ValueError("lineage parent IDs must be iterable") from error
        for parent_id in parent_ids:
            parent = by_path_id.get(parent_id)
            if parent is None:
                raise ValueError("lineage refers to an unknown parent path")
            if parent.get("superlocus_id") != entry.get("superlocus_id"):
                raise ValueError("lineage cannot cross source superloci")
            adjacency[path_id].add(parent_id)
            adjacency[parent_id].add(path_id)

    unseen = set(adjacency)
    families = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        stack = [first]
        family_ids = []
        while stack:
            current = stack.pop()
            family_ids.append(current)
            neighbours = sorted(
                adjacency[current].intersection(unseen), reverse=True
            )
            for neighbour in neighbours:
                unseen.remove(neighbour)
                stack.append(neighbour)
        families.append(tuple(
            by_path_id[path_id] for path_id in sorted(family_ids)
        ))
    return tuple(families)


def _proposal_rank(proposal):
    row = proposal["row"]
    path = candidate_repeat_path(row[1])
    interval = _path_interval(path)
    return (
        len(path),
        interval[1] - interval[0],
        float(row[0]),
        tuple((-start, -end) for start, end in path),
    )


def _proposal_phase_identity(proposal):
    """Return validated source-family identity for one rescue proposal."""

    if not isinstance(proposal, dict):
        raise ValueError("architecture proposals must be mappings")
    row = proposal.get("row")
    if not isinstance(row, (tuple, list)) or len(row) < 2:
        raise ValueError("architecture proposal contains a malformed row")
    lineage = proposal.get("lineage")
    if not isinstance(lineage, dict):
        raise ValueError("architecture proposal lineage must be a mapping")
    path_id = lineage.get("path_id")
    if not isinstance(path_id, str) or not path_id:
        raise ValueError("architecture proposal requires a lineage path ID")
    if "superlocus_id" not in lineage:
        raise ValueError("architecture proposal requires a source superlocus")
    try:
        family_path_ids = tuple(proposal["family_path_ids"])
    except (KeyError, TypeError) as error:
        raise ValueError(
            "architecture proposal requires iterable family path IDs"
        ) from error
    if (
            not family_path_ids
            or any(
                not isinstance(member, str) or not member
                for member in family_path_ids
            )
            or len(set(family_path_ids)) != len(family_path_ids)):
        raise ValueError(
            "architecture proposal family path IDs must be unique strings"
        )
    return {
        "candidate": row[1],
        "geometry": candidate_repeat_path(row[1]),
        "path_id": path_id,
        "superlocus_id": lineage["superlocus_id"],
        "family_path_ids": tuple(sorted(family_path_ids)),
    }


def _cross_family_phase_equivalent(first, second, rule):
    """Return whether proposals are alternate phases of the same occurrences."""

    if (
            first["superlocus_id"] != second["superlocus_id"]
            or first["family_path_ids"] == second["family_path_ids"]):
        return False
    first_path = first["geometry"]
    second_path = second["geometry"]
    if (
            len(first_path) != len(second_path)
            or len(first_path) < rule["minimum_repeat_count"]):
        return False
    maximum_shift = rule["maximum_boundary_shift"]
    return all(
        abs(first_start - second_start) <= maximum_shift
        and abs(first_end - second_end) <= maximum_shift
        for (
            (first_start, first_end),
            (second_start, second_end),
        ) in zip(first_path, second_path)
    )


def _source_phase_arbitration_evidence(result):
    """Drop runtime candidate identities from source-comparison evidence."""

    return {
        "classifier_used": result["classifier_used"],
        "winner_margin": result["winner_margin"],
        "background_window": list(result["background_window"])
        if result["background_window"] is not None else None,
        "background_probabilities": dict(
            result["background_probabilities"] or {}
        ),
        "rule": dict(result["rule"]),
        "candidates": [
            {
                "geometry": [
                    [start, end]
                    for start, end in zip(
                        evidence["geometry"][0],
                        evidence["geometry"][1],
                    )
                ],
                "source_score": evidence["source_score"],
                "boundary_evidence": evidence["boundary_evidence"],
            }
            for evidence in result["candidate_evidence"]
        ],
    }


def arbitrate_cross_family_phase_proposals(full_dna, proposals):
    """Resolve only strict-eligible, cross-family boundary-phase conflicts.

    Proposal selection intentionally happens before this function, so source
    evidence can choose only among candidates that independently passed the
    final rescue policy.  It cannot create an array, extend a path, or promote
    a model-rejected candidate.  Comparisons are limited to equal occurrence
    counts with corresponding source boundaries no more than 12 bases apart.

    A clear source-only winner geometry is retained.  If source evidence
    cannot distinguish a bounded conflict, locus existence is preserved by
    retaining the lexicographically first source-valid geometry as an
    explicitly ambiguous representative.  The representative is never chosen
    by span or classifier score and callers must restrict it to ``Possible``.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    try:
        proposals = tuple(proposals)
    except TypeError as error:
        raise ValueError("architecture proposals must be iterable") from error
    if len({id(proposal) for proposal in proposals}) != len(proposals):
        raise ValueError("architecture proposal identities must be unique")

    rule = dict(SOURCE_CROSS_FAMILY_PHASE_ARBITRATION_RULE)
    identities = {
        id(proposal): _proposal_phase_identity(proposal)
        for proposal in proposals
    }
    ordered = sorted(
        proposals,
        key=lambda proposal: (
            repr(identities[id(proposal)]["superlocus_id"]),
            identities[id(proposal)]["geometry"],
            identities[id(proposal)]["family_path_ids"],
            identities[id(proposal)]["path_id"],
        ),
    )

    # Complete-link clusters prevent a transitive chain of individually
    # bounded phases from silently widening the audited 12-base search domain.
    clusters = []
    for proposal in ordered:
        identity = identities[id(proposal)]
        for cluster in clusters:
            if all(
                    _cross_family_phase_equivalent(
                        identity, identities[id(member)], rule
                    )
                    for member in cluster):
                cluster.append(proposal)
                break
        else:
            clusters.append([proposal])

    retained_ids = {id(proposal) for proposal in proposals}
    arbitration_by_proposal_id = {}
    records = []
    for cluster in clusters:
        if len(cluster) < 2:
            continue
        by_geometry = {}
        for proposal in cluster:
            geometry = identities[id(proposal)]["geometry"]
            by_geometry.setdefault(geometry, []).append(proposal)
        if len(by_geometry) < 2:
            # Duplicate source geometry needs no coordinate arbitration.
            continue

        geometries = tuple(sorted(by_geometry))
        representatives = [
            min(
                by_geometry[geometry],
                key=lambda proposal: identities[id(proposal)]["path_id"],
            )
            for geometry in geometries
        ]
        result = compare_source_backed_phase_candidates(
            full_dna,
            [
                identities[id(proposal)]["candidate"]
                for proposal in representatives
            ],
        )
        cluster_path_ids = tuple(sorted(
            identities[id(proposal)]["path_id"]
            for proposal in cluster
        ))
        if result["status"] == "selected":
            winner_representative = representatives[
                result["winner_candidate_index"]
            ]
            winner_geometry = identities[
                id(winner_representative)
            ]["geometry"]
            winner_proposals = tuple(by_geometry[winner_geometry])
            winner_ids = {id(proposal) for proposal in winner_proposals}
            removed = tuple(
                proposal for proposal in cluster
                if id(proposal) not in winner_ids
            )
            retained_ids.difference_update(
                id(proposal) for proposal in removed
            )
            status = "selected"
            winner_path_id = min(
                identities[id(proposal)]["path_id"]
                for proposal in winner_proposals
            )
            retained_path_ids = tuple(sorted(
                identities[id(proposal)]["path_id"]
                for proposal in winner_proposals
            ))
            boundary_status = "policy_selected"
        else:
            # Every proposal in this cluster independently passed the strict
            # existence policy.  A phase tie is therefore uncertainty about
            # coordinates, not evidence that the locus is absent.  Choose one
            # deterministic source-valid geometry for representation only.
            representative = representatives[0]
            representative_geometry = identities[
                id(representative)
            ]["geometry"]
            representative_proposals = tuple(
                by_geometry[representative_geometry]
            )
            representative_ids = {
                id(proposal) for proposal in representative_proposals
            }
            removed = tuple(
                proposal for proposal in cluster
                if id(proposal) not in representative_ids
            )
            retained_ids.difference_update(
                id(proposal) for proposal in removed
            )
            status = "ambiguous_representative"
            winner_path_id = None
            retained_path_ids = tuple(sorted(
                identities[id(proposal)]["path_id"]
                for proposal in representative_proposals
            ))
            boundary_status = "ambiguous"

        removed_path_ids = tuple(sorted(
            identities[id(proposal)]["path_id"]
            for proposal in removed
        ))
        record = {
            "policy": dict(rule),
            "coordinate_system": "zero_based_half_open",
            "superlocus_id": identities[id(cluster[0])]["superlocus_id"],
            "status": status,
            "reason": result["reason"],
            "winner_path_id": winner_path_id,
            "representative_path_id": (
                min(retained_path_ids) if retained_path_ids else None
            ),
            "boundary_status": boundary_status,
            "reported_geometry_role": "representative",
            "retained_path_ids": retained_path_ids,
            "removed_path_ids": removed_path_ids,
            "competitor_path_ids": cluster_path_ids,
            "source_evidence": _source_phase_arbitration_evidence(result),
        }
        records.append(record)
        for proposal in cluster:
            if id(proposal) in retained_ids:
                arbitration_by_proposal_id[id(proposal)] = record

    retained = tuple(sorted(
        (
            {
                **proposal,
                **(
                    {
                        "source_phase_arbitration": (
                            arbitration_by_proposal_id[id(proposal)]
                        )
                    }
                    if id(proposal) in arbitration_by_proposal_id
                    else {}
                ),
            }
            for proposal in proposals
            if id(proposal) in retained_ids
        ),
        key=lambda proposal: (
            _proposal_rank(proposal),
            tuple(sorted(proposal["family_path_ids"])),
            proposal["lineage"]["path_id"],
        ),
        reverse=True,
    ))
    return retained, tuple(records)


def select_lineage_family_proposals(
    lineage,
    rows_by_candidate_id,
    decisions,
    seed_roles,
):
    """Select a unique maximal eligible path from each supported family.

    A family must contain either a strict primary seed or a four-repeat
    satellite that can be conservatively terminal-completed.  A three-repeat
    satellite is corroborating graph evidence only and can never originate an
    accepted lineage by itself.  Among eligible final paths, a candidate can
    suppress another only through exact contiguous path containment; multiple
    incomparable maxima cause abstention rather than score-based selection.
    """

    proposals = []
    for family in lineage_candidate_families(lineage):
        supported = any(
            seed_roles.get(id(entry.get("candidate"))) == "primary"
            or (
                seed_roles.get(id(entry.get("candidate"))) == "satellite"
                and len(entry["candidate"].list_repeats) == 4
            )
            for entry in family
        )
        if not supported:
            continue
        legacy_locus_keys = []
        for member in family:
            details = member.get("details", {})
            if not isinstance(details, dict):
                raise ValueError("lineage details must be a mapping")
            raw_locus_keys = details.get("legacy_locus_keys", ())
            if isinstance(raw_locus_keys, (str, bytes)):
                raise ValueError("legacy locus keys must be an iterable")
            try:
                raw_locus_keys = tuple(raw_locus_keys)
            except TypeError as error:
                raise ValueError(
                    "legacy locus keys must be an iterable"
                ) from error
            for locus_key in raw_locus_keys:
                try:
                    hash(locus_key)
                except TypeError as error:
                    raise ValueError(
                        "legacy locus keys must be hashable"
                    ) from error
                if locus_key not in legacy_locus_keys:
                    legacy_locus_keys.append(locus_key)

        eligible = []
        for entry in family:
            if entry.get("architecture_commit_eligible", True) is not True:
                continue
            candidate = entry.get("candidate")
            row = rows_by_candidate_id.get(id(candidate))
            decision = decisions.get(id(candidate))
            if row is None or decision is None:
                raise ValueError(
                    "every lineage candidate must have a scored decision"
                )
            requires_source_coding_context = entry.get(
                "requires_source_coding_context", False
            )
            coding_veto_assessed = getattr(
                getattr(decision, "stats", None),
                "coding_veto_assessed",
                False,
            )
            if (
                    decision.decision == ELIGIBLE
                    and (
                        not requires_source_coding_context
                        or coding_veto_assessed is True
                    )):
                eligible.append({
                    "row": row,
                    "decision": decision,
                    "lineage": entry,
                    "family_path_ids": tuple(
                        member["path_id"] for member in family
                    ),
                    "family_seed_roles": {
                        member["path_id"]: seed_roles[
                            id(member["candidate"])
                        ]
                        for member in family
                        if id(member["candidate"]) in seed_roles
                    },
                    "family_legacy_locus_keys": tuple(legacy_locus_keys),
                })
        maximal = []
        for proposal in eligible:
            proposal_path = candidate_repeat_path(proposal["row"][1])
            if not any(
                len(candidate_repeat_path(other["row"][1]))
                > len(proposal_path)
                and path_is_contiguous_subpath(
                    proposal_path,
                    candidate_repeat_path(other["row"][1]),
                )
                for other in eligible
            ):
                maximal.append(proposal)
        if len(maximal) == 1:
            proposals.append(maximal[0])
    return tuple(sorted(proposals, key=_proposal_rank, reverse=True))


def plan_architecture_rescue_transitions(categories, proposals):
    """Plan additions/replacements without mutating category containers.

    A selected geometry may replace accepted calls only when every overlapping
    accepted repeat path is an exact contiguous subpath.  Span overlap alone
    is treated as a conflict and causes abstention.  Planning is sequential in
    deterministic architecture rank so two proposed arrays cannot both claim
    the same source interval.
    """

    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("architecture transitions require five categories")
    row_locations = {}
    accepted = []
    for category_index, category in enumerate(categories):
        for locus_key, rows in category.items():
            for row in rows:
                if not isinstance(row, (list, tuple)) or len(row) < 2:
                    raise ValueError("category contains a malformed row")
                previous = row_locations.get(id(row))
                location = (category_index, locus_key)
                if previous is not None and previous != location:
                    raise ValueError(
                        "one scored row occurs in multiple category regions"
                    )
                row_locations[id(row)] = location
                if category_index in (0, 2):
                    accepted.append({
                        "row": row,
                        "locus_key": locus_key,
                    })

    transitions = []
    for proposal in sorted(proposals, key=_proposal_rank, reverse=True):
        row = proposal["row"]
        candidate = row[1]
        selected_path = candidate_repeat_path(candidate)
        selected_interval = _path_interval(selected_path)
        overlapping_records = [
            record for record in accepted
            if _intervals_overlap(
                selected_interval,
                _path_interval(candidate_repeat_path(record["row"][1])),
            )
        ]
        dominated_records = [
            record for record in overlapping_records
            if path_is_contiguous_subpath(
                candidate_repeat_path(record["row"][1]), selected_path
            )
        ]
        overlapping = [record["row"] for record in overlapping_records]
        dominated = [record["row"] for record in dominated_records]
        dominated_ids = {id(accepted_row) for accepted_row in dominated}
        conflicts = [
            accepted_row for accepted_row in overlapping
            if id(accepted_row) not in dominated_ids
        ]
        already_accepted = any(
            record["row"][1] is candidate for record in accepted
        )

        if conflicts:
            transitions.append({
                **proposal,
                "action": "abstain",
                "reason": "overlapping_non_subpath_accepted_call",
                "dominated_rows": tuple(dominated),
                "conflicting_rows": tuple(conflicts),
            })
            continue
        architecture_override = float(row[0]) < 0.5
        if architecture_override and overlapping:
            transitions.append({
                **proposal,
                "action": "abstain",
                "reason": "architecture_rescue_cannot_replace_accepted_call",
                "dominated_rows": tuple(dominated),
                "conflicting_rows": (),
            })
            continue
        strict_extension = any(
            len(candidate_repeat_path(accepted_row[1])) < len(selected_path)
            for accepted_row in dominated
        )
        if already_accepted and not strict_extension:
            transitions.append({
                **proposal,
                "action": "unchanged",
                "reason": "selected_geometry_already_accepted",
                "dominated_rows": (),
                "conflicting_rows": (),
            })
            continue
        if dominated and not strict_extension:
            transitions.append({
                **proposal,
                "action": "abstain",
                "reason": "accepted_geometry_not_strictly_extended",
                "dominated_rows": tuple(dominated),
                "conflicting_rows": (),
            })
            continue

        action = "replace" if dominated else "add"
        own_location = row_locations.get(id(row))
        family_locus_keys = tuple(
            proposal.get("family_legacy_locus_keys", ())
        )
        relevant_locus_keys = []
        for locus_key in (
            [record["locus_key"] for record in dominated_records]
            + ([own_location[1]] if own_location is not None else [])
            + list(family_locus_keys)
        ):
            if locus_key not in relevant_locus_keys:
                relevant_locus_keys.append(locus_key)

        dominated_locus_keys = []
        for record in dominated_records:
            if record["locus_key"] not in dominated_locus_keys:
                dominated_locus_keys.append(record["locus_key"])
        if len(dominated_locus_keys) == 1:
            target_locus_key = dominated_locus_keys[0]
        elif own_location is not None:
            target_locus_key = own_location[1]
        elif relevant_locus_keys:
            target_locus_key = min(
                relevant_locus_keys,
                key=lambda value: (type(value).__name__, repr(value)),
            )
        else:
            # Synthetic/integration callers may lack augmentation lineage.
            # The exact source interval is an unambiguous last-resort region.
            target_locus_key = selected_interval

        transition = {
            **proposal,
            "action": action,
            "reason": "source_architecture_eligible",
            "dominated_rows": tuple(dominated),
            "conflicting_rows": (),
            "target_category_index": (
                2
                if (
                    proposal.get("source_phase_arbitration", {}).get(
                        "status"
                    ) == "ambiguous_representative"
                )
                else (0 if float(row[0]) >= 0.75 else 2)
            ),
            "source_interval": selected_interval,
            "target_locus_key": target_locus_key,
            "migrated_locus_keys": tuple(
                locus_key
                for locus_key in relevant_locus_keys
                if locus_key != target_locus_key
            ),
        }
        transitions.append(transition)
        accepted = [
            record for record in accepted
            if id(record["row"]) not in dominated_ids
        ]
        accepted.append({
            "row": row,
            "locus_key": target_locus_key,
        })
    return tuple(transitions)


__all__ = [
    "SOURCE_AUGMENTATION_SATELLITE_POLICY",
    "SOURCE_CROSS_FAMILY_PHASE_ARBITRATION_RULE",
    "arbitrate_cross_family_phase_proposals",
    "candidate_repeat_path",
    "evaluate_scored_architecture_candidate",
    "lineage_candidate_families",
    "named_features_from_scored_row",
    "path_is_contiguous_subpath",
    "plan_architecture_rescue_transitions",
    "select_augmentation_seed_ids",
    "select_lineage_family_proposals",
]
