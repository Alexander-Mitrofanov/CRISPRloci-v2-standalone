"""Extend a complete boundary family only when all its members support the same terminal site."""

from components.evidence.initial_competition import (
    protected_ordinary_reference,
    with_competition_evidence,
)
from components.reconstruction.reconstruction_boundary_family import (
    METHOD,
    boundary_family_summary,
    describe_family,
    digest,
    display_order,
    DISPLAY_RULE,
    family_occurrence_evidence,
    geometry_id,
    negative_signature,
)
from components.reconstruction.reconstruction_validation import (
    needs_coding_context,
    validate_reconstruction_candidate,
)
from components.reconstruction.reconstruction_selection import RowLocation
from components.refinement.enhancement_search import (
    build_source_backed_candidate_from_boundaries,
)
from components.refinement.periodic_completion import plan_path
from components.refinement.potential_degenerate import added_member_evidence, geometry
from components.rescue.architecture_rescue_coding_context import SourceCDS
from components.scoring.evaluation import build_advanced_fuzzy_search_filter


def terminal_projection(source, members, ordinary):
    """Retain each member's internal intervals even if fitting suggests a different width."""
    projected = []
    common = None
    for parent in members:
        old = geometry(parent)
        proposal = plan_path(source, parent, ordinary)
        if len(proposal) != len(old) + 1:
            return None
        if proposal[0][1] <= old[0][0]:
            addition, fitted = proposal[0], proposal[1:]
            full = (addition,) + old
        elif proposal[-1][0] >= old[-1][1]:
            addition, fitted = proposal[-1], proposal[:-1]
            full = old + (addition,)
        else:
            return None
        if (
            any(b <= c or d <= a for (a, b), (c, d) in zip(old, fitted))
            or common is not None
            and addition != common
        ):
            return None
        common = addition
        candidate = build_source_backed_candidate_from_boundaries(
            source, parent, *zip(*full)
        )
        candidate.potential_degenerate_repeats = {
            **getattr(parent, "potential_degenerate_repeats", {}),
            **added_member_evidence(source, parent, old, full, ordinary),
        }
        projected.append(candidate)
    return common, projected


def _parent_proof(proof):
    parent = dict(proof)
    extension = parent.pop("terminal_extension")
    parent["members"] = extension["parent_members"]
    parent["display_representation_id"] = extension["parent_display_representation_id"]
    parent["descriptor"] = describe_family(
        parent["source_sha256"],
        parent["anchor_original"],
        parent["members"],
        parent["display_representation_id"],
    )
    return parent


def _build_proof(source, parent, ordinary, addition, rows, checks):
    from components.reconstruction.reconstruction_transaction import _restore, _snapshot

    if (
        not all(check.qualified for check in checks)
        or any(row[0] < parent["policy"]["minimum_anchor_score"] for row in rows)
        or len({negative_signature(check.evidence) for check in checks}) != 1
    ):
        return None
    witnesses = [(parent["anchor_original"], parent["anchor_initial_competition"])]
    witnesses.extend(
        (m["original"], m["initial_competition"])
        for m in parent["members"] + parent["additional_peer_representations"]
    )
    references = [
        protected_ordinary_reference(
            with_competition_evidence(_restore(source, snap), obs), source=source
        )
        for snap, obs in witnesses
        if obs is not None
    ]
    if (
        family_occurrence_evidence(references, [geometry(row[1]) for row in rows])[
            "status"
        ]
        != "preserved"
    ):
        return None
    members = []
    by_id = {}
    for old, row, check in zip(parent["members"], rows, checks):
        path = geometry(row[1])
        gid = geometry_id(path)
        fresh = _snapshot(row)
        rid = RowLocation(1, None, 0, row).representation_id
        old_mapping = old["relationship"]["correspondence"][
            "anchor_to_candidate_indexes"
        ]
        offset = int(path[0] == addition)
        members.append(
            {
                "parent_geometry_id": old["location"]["geometry_id"],
                "fresh": fresh,
                "validation": check.as_dict(),
                "location": {"geometry_id": gid, "representation_id": rid},
                "relationship": {
                    "correspondence": {
                        "anchor_to_candidate_indexes": [i + offset for i in old_mapping]
                    }
                },
            }
        )
        by_id[rid] = row
    members.sort(key=lambda m: m["location"]["geometry_id"])
    anchor = parent["anchor_original"]
    shared_end = len({m["fresh"]["repeat_intervals"][-1][1] for m in members}) == 1
    display = min(
        members,
        key=lambda m: display_order(
            m["fresh"]["score"],
            m["fresh"]["repeat_intervals"],
            anchor["repeat_intervals"],
            m["location"]["geometry_id"],
            m["location"]["representation_id"],
            DISPLAY_RULE,
            shared_end,
        ),
    )
    rid = display["location"]["representation_id"]
    proof = dict(
        parent,
        members=members,
        display_representation_id=rid,
        descriptor=describe_family(parent["source_sha256"], anchor, members, rid),
        terminal_extension={
            "ordinary_edit_limit": ordinary,
            "added_interval": list(addition),
            "parent_members": parent["members"],
            "parent_display_representation_id": parent["display_representation_id"],
        },
    )
    return {
        "method": METHOD,
        "boundary_family": proof,
        "evidence_sha256": digest(proof),
    }, by_id[rid]


def _overlap(span, other):
    return span[0] < other[1] and other[0] < span[1]


def _span(row):
    family = boundary_family_summary(row)
    if family:
        return family["union_span"]
    path = geometry(row[1])
    return path[0][0], path[-1][1]


def extend_selected_families(source, output, decisions, parameters, score, source_cds):
    ordinary = parameters["param_max_edit_distance"]
    structural = build_advanced_fuzzy_search_filter(parameters)
    active = {
        id(row): _span(row)
        for index in (0, 2)
        for group in output[index].values()
        for row in group
    }
    proposals = []
    extensions = []
    for decision in decisions:
        anchor = decision.anchor
        slot = output[0][anchor.locus_key]
        row = slot[anchor.row_index]
        parent = row[3]["classification_basis"]["boundary_family"]
        native = {
            location.geometry_id: checked.row[1]
            for location, _, checked in decision.family_members
        }
        result = terminal_projection(
            source,
            [native[m["location"]["geometry_id"]] for m in parent["members"]],
            ordinary,
        )
        if result:
            addition, candidates = result
            # Reserve proposed spans together so overlapping simultaneous extensions both abstain.
            lo, hi = active[id(row)]
            active[id(row)] = min(lo, addition[0]), max(hi, addition[1])
            proposals.append(
                (row, slot, anchor.row_index, parent, addition, candidates)
            )
    for row, slot, index, parent, addition, candidates in proposals:
        if any(
            _overlap(active[id(row)], span)
            for key, span in active.items()
            if key != id(row)
        ):
            continue
        key = "boundary-family-terminal-extension"
        result = score({key: candidates})
        if (
            set(result) != {key}
            or len(result[key]) != len(candidates)
            or any(r[1] is not c for r, c in zip(result[key], candidates))
        ):
            raise ValueError("boundary family terminal scoring changed the request")
        rows = result[key]
        cds = (
            tuple(source_cds()) if any(map(needs_coding_context, candidates)) else None
        )
        checks = [
            validate_reconstruction_candidate(source, r, r, structural, source_cds=cds)
            for r in rows
        ]
        built = _build_proof(source, parent, ordinary, addition, rows, checks)
        if built:
            basis, display = built
            slot[index] = list(display[:3]) + [dict(row[3], classification_basis=basis)]
            extensions.append(
                {
                    "parent_family_id": parent["descriptor"]["family_id"],
                    "family_id": basis["boundary_family"]["descriptor"]["family_id"],
                    "added_interval": list(addition),
                    "evidence_sha256": basis["evidence_sha256"],
                }
            )
    return extensions


def validate_extension(source, current, categories, ordinary_witness_keys):
    """Replay original family admission, each independent placement, and its complete validated path."""
    from components.reconstruction.reconstruction_transaction import (
        _restore,
        _snapshot,
        validate_boundary_family_basis,
    )

    basis = current[3]["classification_basis"]
    proof = basis["boundary_family"]
    parent = _parent_proof(proof)
    parent_basis = {
        "method": METHOD,
        "boundary_family": parent,
        "evidence_sha256": digest(parent),
    }
    original_rows = [_restore(source, m["fresh"]) for m in parent["members"]]
    original = next(
        r
        for r, m in zip(original_rows, parent["members"])
        if m["location"]["representation_id"] == parent["display_representation_id"]
    )
    original = original + [{"classification_basis": parent_basis}]
    trial = [
        {
            k: [original if r is current else r for r in group]
            for k, group in categories[0].items()
        }
    ] + list(categories[1:])
    validate_boundary_family_basis(
        source,
        original[1],
        parent_basis,
        original[0],
        0,
        original,
        trial,
        ordinary_witness_keys=ordinary_witness_keys,
    )
    ext = proof["terminal_extension"]
    ordinary = ext["ordinary_edit_limit"]
    if type(ordinary) is not int or ordinary < 0:
        raise ValueError("invalid boundary family terminal edit limit")
    projected = terminal_projection(source, [r[1] for r in original_rows], ordinary)
    if projected is None or list(projected[0]) != ext["added_interval"]:
        raise ValueError("boundary family terminal placement differs")
    structural = build_advanced_fuzzy_search_filter(parent["structural_parameters"])
    by_parent = {m["parent_geometry_id"]: m for m in proof["members"]}
    rows, checks = [], []
    for member, candidate in zip(parent["members"], projected[1]):
        record = by_parent[member["location"]["geometry_id"]]
        row = _restore(source, record["fresh"])
        if (
            geometry(row[1]) != geometry(candidate)
            or row[1].list_gaped_repeats != candidate.list_gaped_repeats
        ):
            raise ValueError(
                "boundary family terminal projection changed internal boundaries or alignment"
            )
        coding = record["validation"]["evidence"]["coding_tandem"]
        cds = (
            tuple(SourceCDS(**raw) for raw in coding["considered_source_cds"])
            if coding["source_cds_assessed"]
            else None
        )
        rows.append(row)
        checks.append(
            validate_reconstruction_candidate(
                source, row, row, structural, source_cds=cds
            )
        )
    rebuilt = _build_proof(source, parent, ordinary, projected[0], rows, checks)
    if (
        rebuilt is None
        or rebuilt[0] != basis
        or _snapshot(rebuilt[1]) != _snapshot(current)
        or any(
            _overlap(_span(current), _span(r))
            for index in (0, 2)
            for group in categories[index].values()
            for r in group
            if r is not current
        )
    ):
        raise ValueError("boundary family terminal evidence or accepted span differs")
    return basis
