"""Bounded, source-only repair of weak overextended repeat terminals.

This stage proposes geometry, never scores or classifies it. Coordinates in
its audit records are zero-based and half-open; normal reports keep their
existing one-based inclusive convention.
"""

from collections import Counter
import hashlib
from components.candidate_origins import mark_repair_origin

from components.components_evaluated_arrays_enhancement import (
    _source_phase_boundary_evidence,
    build_source_backed_candidate_from_boundaries,
    validate_candidate_source_reconstruction,
)


POLICY_VERSION = "terminal-only-trim-v1"
MAX_PARENTS_PER_LOCUS = 8
MAX_TRIM_COLUMNS = 3


def source_geometry(candidate):
    starts = tuple(candidate.list_repeat_starts)
    return starts, tuple(
        start + len(repeat)
        for start, repeat in zip(starts, candidate.list_repeats)
    )


def classification_restriction(candidate):
    eligible = getattr(candidate, "source_path_classification_eligible", True)
    if type(eligible) is not bool:
        return "invalid_eligibility_marker"
    if not eligible:
        return "classification_ineligible"
    if any(
            getattr(record, "derivation", None) == "local_chain_family_parent"
            for record in getattr(candidate, "source_path_provenance", ())):
        return "local_chain_family_parent"
    return None


def _weak_removed_columns(alignment, left, right, dominance):
    width = len(alignment[0])
    indexes = list(range(left)) + list(range(width - right, width))
    for index in indexes:
        column = [repeat[index] for repeat in alignment]
        # Match the existing terminal predicate's occupied-column definition.
        if column.count(" ") / len(column) >= 0.5:
            continue
        counts = Counter(base for base in column if base not in " -")
        if counts and max(counts.values()) / sum(counts.values()) >= dominance:
            return False
    return True


def propose_terminal_repairs(full_dna, candidates, structural_filter,
                             source_sha256=None):
    """Return at most one minimal supported child per bounded parent.

    Missing source is a backward-compatible no-op. Present invalid source
    geometry raises rather than silently reconstructing a different claim.
    The input objects and their score-relevant fields are never mutated.
    """
    stats = Counter(candidates_seen=len(candidates), parents_attempted=0,
                    hypotheses_considered=0, children_generated=0)
    records = []
    if full_dna is None:
        stats["no_source_dna"] = len(candidates)
        return [], records, dict(stats)
    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")

    by_geometry = {}
    for candidate in candidates:
        by_geometry.setdefault(source_geometry(candidate), []).append(candidate)
    stats["duplicate_geometries"] = len(candidates) - len(by_geometry)
    # Decide eligibility in deterministic priority order. Once eight parents
    # qualify, no lower-priority geometry can enter the bounded proposal set.
    # Avoid canonicalizing that remainder, and report it as unexamined rather
    # than pretending to know how many would have passed the biological gate.
    ordered = sorted(by_geometry.items(), key=lambda item: (
        -(item[0][1][-1] - item[0][0][0]), -len(item[0][0]), item[0],
    ))
    qualifying = []
    stats["parent_budget_skipped"] = 0
    stats["canonical_eligibility_unexamined"] = 0
    for geometry, members in ordered:
        # A duplicate must not launder a restriction from another derivation.
        restrictions = sorted({reason for member in members
                               if (reason := classification_restriction(member))})
        if restrictions:
            stats[restrictions[0]] += 1
            continue
        if len(geometry[0]) < 5:
            stats["insufficient_repeats"] += 1
            continue
        if any(getattr(member, "boundary_repair_records", ()) for member in members):
            stats["already_repaired"] += 1
            continue
        # Every contributing duplicate must agree with the source; choosing
        # one representative must not conceal a malformed sibling claim.
        for member in members:
            validate_candidate_source_reconstruction(
                full_dna, member.list_repeats, member.list_repeat_starts,
                member.list_spacers,
            )
        # These gates read repeat/spacer sequence and lengths, never the raw
        # alignment. All alignment-dependent biological decisions below use
        # the deterministic source-only representation, including singleton
        # inputs whose legacy tail columns happen to hide weak terminals.
        parent = members[0]
        if not structural_filter._passes_non_gapped_predicates(parent):
            stats["other_structural_failure"] += 1
            continue
        if len(qualifying) == MAX_PARENTS_PER_LOCUS:
            stats["parent_budget_skipped"] += 1
            stats["canonical_eligibility_unexamined"] += 1
            continue
        stats["canonical_eligibility_examined"] += 1
        # After source verification, identical repeat strings prove the
        # source-only alignment is gap-free and its terminal columns pass.
        if len(set(parent.list_repeats)) == 1:
            stats["canonical_eligibility_trivial_exact"] += 1
            stats["terminal_columns_pass"] += 1
            continue
        canonical = build_source_backed_candidate_from_boundaries(
            full_dna, None, *geometry,
        )
        if structural_filter._filter_by_column(canonical):
            stats["terminal_columns_pass"] += 1
            continue
        evidence = structural_filter.evidence(canonical)
        failed = [name for name, status in evidence.predicate_results
                  if status == "failed"]
        if not evidence.structural_valid or failed != ["repeat_terminal_columns"]:
            stats["other_structural_failure"] += 1
            continue
        qualifying.append((geometry, canonical, evidence))

    stats["qualifying_parents"] = len(qualifying)
    children = []
    known_geometries = set(by_geometry)
    for geometry, canonical, evidence in qualifying:
        stats["parents_attempted"] += 1
        record = {
            "policy": POLICY_VERSION,
            "coordinate_system": "zero_based_half_open",
            "source_sha256": source_sha256 or hashlib.sha256(
                full_dna.encode("ascii")).hexdigest(),
            "parent_geometry": geometry,
            "parent_eligibility": True,
            "structural_rule": dict(vars(structural_filter)),
            "parent_failures": ["repeat_terminal_columns"],
            "original_checks": evidence.as_dict(),
            "status": "abstained", "reason": "no_supported_trim",
            "hypotheses_considered": 0, "variant_rejections": {},
        }
        rejected = Counter()
        supported = {}
        alignment = canonical.list_repeats_gaped
        for left in range(MAX_TRIM_COLUMNS + 1):
            for right in range(MAX_TRIM_COLUMNS + 1):
                if not left and not right:
                    continue
                stats["hypotheses_considered"] += 1
                record["hypotheses_considered"] += 1
                if not _weak_removed_columns(
                        alignment, left, right, structural_filter.column_dominance):
                    rejected["conserved_removed_column"] += 1
                    continue
                deltas = tuple((
                    sum(base not in " -" for base in repeat[:left]),
                    sum(base not in " -" for base in repeat[len(repeat)-right:])
                    if right else 0,
                ) for repeat in alignment)
                starts = tuple(start + delta[0]
                               for start, delta in zip(geometry[0], deltas))
                ends = tuple(end - delta[1]
                             for end, delta in zip(geometry[1], deltas))
                new_geometry = starts, ends
                if new_geometry == geometry or any(
                        end <= start for start, end in zip(starts, ends)):
                    rejected["empty_or_unchanged_geometry"] += 1
                    continue
                boundary = _source_phase_boundary_evidence(
                    full_dna, new_geometry, (geometry, new_geometry),
                )
                if not boundary["passes"]:
                    rejected["weak_source_boundary"] += 1
                    continue
                child = build_source_backed_candidate_from_boundaries(
                    full_dna, None, starts, ends,
                )
                checks = structural_filter.evidence(child)
                if not checks.passes:
                    rejected["structural_filter"] += 1
                    continue
                proposal = (left + right, left, right, child, deltas,
                            boundary, checks.as_dict())
                previous = supported.get(new_geometry)
                if previous is None or proposal[:3] < previous[:3]:
                    supported[new_geometry] = proposal
        record["variant_rejections"] = dict(rejected)
        if supported:
            minimum = min(proposal[0] for proposal in supported.values())
            winners = [(geometry, proposal) for geometry, proposal in supported.items()
                       if proposal[0] == minimum]
            if len(winners) != 1:
                record["reason"] = "ambiguous_minimum_trim"
            else:
                child_geometry, proposal = winners[0]
                _, left, right, child, deltas, boundary, checks = proposal
                record.update(
                    child_geometry=child_geometry, left_trim_columns=left,
                    right_trim_columns=right, consumed_bases=deltas,
                    boundary_evidence=boundary, resulting_checks=checks,
                )
                if child_geometry in known_geometries:
                    record["reason"] = "child_geometry_already_present"
                else:
                    record.update(status="generated", reason="unique_minimum_supported_trim")
                    child.boundary_repair_records = (record,)
                    mark_repair_origin(child, record)
                    children.append(child)
                    known_geometries.add(child_geometry)
                    stats["children_generated"] += 1
        records.append(record)
    return children, records, dict(stats)
