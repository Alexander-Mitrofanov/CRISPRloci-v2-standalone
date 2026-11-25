"""Select one source-supported interpretation of an already accepted array.

The unchanged interpretation is a retained, scored witness. Selection is
replayed at publication against the locus inventory and original chromosome.
No reviewer record, reference-tool call or accession participates.
"""

from dataclasses import asdict
from typing import NamedTuple

from components.evidence.initial_competition import (
    guard_reconstructed_occurrences,
    protected_ordinary_reference,
)
from components.evidence.occurrence_protection import occurrence_relation
from components.reconstruction.reconstruction_selection import CandidateValidation
from components.reconstruction.reconstruction_transaction import (
    STRUCTURAL_PARAMETERS,
    _digest,
    _plain,
    _restore,
    _snapshot,
)
from components.reconstruction.reconstruction_validation import (
    _known_negative_annotations,
    needs_coding_context,
    validate_reconstruction_candidate,
)
from components.refinement.enhancement_search import (
    _source_repeat_alignment,
    build_source_backed_candidate_from_boundaries,
)
from components.refinement.potential_degenerate import added_member_evidence, geometry
from components.refinement.repeat_boundary_repair import classification_restriction
from components.rescue.architecture_rescue_coding_context import SourceCDS
from components.scoring.evaluation import build_advanced_fuzzy_search_filter

METHOD = "source_supported_competing_geometry"
PARENT_METHODS = {
    "model_score_band",
    "source_backed_internal_repeat_completion",
    "source_backed_candidate_reconstruction",
}
DUPLICATION = "tandem_duplication_resegmentation"
PROVISIONAL_CORE_RULES = {
    "conserved_flank_with_adjacent_exemplars",
    "periodic_internal_insertion",
    "retained_substitution_width_extension",
}


def _overlaps(first, second):
    return first[0][0] < second[-1][1] and second[0][0] < first[-1][1]


def _basis(row):
    return row[3].get("classification_basis", {}) if len(row) == 4 else {}


def _eligible_parent(cat, parent):
    basis = _basis(parent)
    if (
        cat not in (0, 2)
        or basis.get("method") not in PARENT_METHODS
        or parent[0] < 0.5
    ):
        return False
    if basis["method"] == "model_score_band":
        if basis.get("final_category_before_source_transactions") != (
            "Bona-fide" if cat == 0 else "Possible"
        ):
            return False
    return (
        not classification_restriction(parent[1])
        and _known_negative_annotations(parent)["status"] == "passed"
    )


def _retained(categories, key, old, exclude=None):
    scores = {}
    for category in categories:
        for row in category.get(key, ()):
            if row is exclude or _basis(row).get("method") == METHOD:
                continue
            path = geometry(row[1])
            if (
                3 <= len(path) <= 128
                and _overlaps(path, old)
                and not classification_restriction(row[1])
                and _known_negative_annotations(row)["status"] == "passed"
            ):
                scores[path] = max(float(row[0]), scores.get(path, 0.0))
    return scores


def _selection(source, old, retained, parent):
    from components.refinement import (
        reviewed_boundary_selection,
        reviewed_membership_selection,
        reviewed_terminal_selection,
    )

    proposals = []
    for selector in (
        reviewed_boundary_selection,
        reviewed_membership_selection,
        reviewed_terminal_selection,
    ):
        result = selector.select(source, old, retained)
        if result is not None:
            proposals.append(result)
    if not proposals:
        return reviewed_boundary_selection.select_alignment(source, old, parent)
    if len(proposals) == 1 and proposals[0][0] != old:
        return proposals[0]
    return None


def _preservation(source, parent, path, rule):
    old = geometry(parent[1])
    relation = occurrence_relation(old, path)
    if relation["status"] == "unique" and guard_reconstructed_occurrences(
        path, parent, source=source
    ):
        return {"status": "all_parent_occurrences_preserved"}
    if rule != DUPLICATION or len(path) != len(old) + 1:
        return None
    missing = relation["reference_indexes_inside_candidate_spacers"]
    if len(missing) != 1 or relation["reference_indexes_without_edge"] != missing:
        return None
    j = missing[0]
    reference = protected_ordinary_reference(parent, source=source)
    if reference is None or tuple(map(tuple, reference["repeat_intervals"])) != old:
        return None
    if occurrence_relation(old[:j] + old[j + 1 :], path)["status"] != "unique":
        return None
    return {
        "status": "competing_spacer_duplication",
        "displaced_repeat_index": j + 1,
        "displaced_interval": list(old[j]),
        "preserved_ordinary_reference": reference,
    }


def _cds_scope(dna, candidate, rule, details, source_cds):
    coding = ()
    if needs_coding_context(candidate) or rule == "periodic_internal_insertion":
        # Bind the assessed, local CDS inventory; no whole-chromosome list is published.
        start, end = geometry(candidate)[0][0], geometry(candidate)[-1][1]
        coding = tuple(c for c in source_cds() if c.start < end and start < c.end)
    if rule == "periodic_internal_insertion":
        start, end = details["split_spacer_interval"]
        if any(c.start <= start and end <= c.end for c in coding):
            return None
    return coding


def _assess(dna, row, structural, cds, old, retained, rule):
    return validate_reconstruction_candidate(
        dna,
        row,
        row,
        structural,
        source_cds=cds,
        allele_flank_context=(old, retained)
        if rule == "retained_allele_flanks"
        else None,
    )


def _score_one(dna, parent, path, key, score):
    candidate = build_source_backed_candidate_from_boundaries(
        dna, parent[1], *zip(*path)
    )
    (fresh,) = score({key: [candidate]})[key]
    if fresh[1] is not candidate:
        raise ValueError("competing geometry scorer changed candidate identity")
    return fresh


class _Assessment(NamedTuple):
    row: list
    cds: tuple
    check: CandidateValidation

    def as_alternative(self):
        return {
            "status": "unselected",
            "candidate": _snapshot(self.row),
            "qualified": self.check.qualified,
            "validation": self.check.evidence,
            "local_cds": [asdict(c) for c in self.cds],
        }


def _width_competition(choice):
    return choice[1] in {"retained_indel_allele_extension", "ambiguous_edge_width_competition"} and bool(
        choice[2].get("alternative_repeat_intervals")
    )


def _proposal_span(choice):
    path = choice[0]
    if not _width_competition(choice):
        return path
    rival = choice[2]["alternative_repeat_intervals"]
    return ((min(path[0][0], rival[0][0]), max(path[-1][1], rival[-1][1])),)


def _rank_choice(dna, parent, cat, choice, candidates):
    if not _width_competition(choice):
        return choice
    if choice[1] == "ambiguous_edge_width_competition":
        from components.refinement.reviewed_boundary_selection import rank_alignment as rank
    else:
        from components.refinement.reviewed_membership_selection import rank_extension as rank

    compact = candidates[choice[0]]
    uniform_path = tuple(map(tuple, choice[2]["alternative_repeat_intervals"]))
    uniform = candidates[uniform_path]
    eligible = (
        uniform.check.qualified
        and _category(cat, parent[0], uniform.row[0], choice[1]) is not None
        and _preservation(dna, parent, uniform_path, choice[1]) is not None
    )
    return rank(
        choice, uniform_eligible=eligible, compact_score=compact.row[0],
        uniform_score=uniform.row[0],
    )


def _category(parent_category, parent_score, fresh_score, rule):
    if parent_score < 0.5 or (
        fresh_score < 0.5
        and (parent_category != 2 or rule not in PROVISIONAL_CORE_RULES)
    ):
        return None
    # A qualified geometry below the model's ordinary acceptance threshold
    # can only refine an already accepted locus, at provisional confidence.
    return 0 if parent_category == 0 and fresh_score >= 0.75 else 2


def _uncertainty(rule, details):
    if rule == DUPLICATION:
        return "competing_membership_in_duplicated_spacer"
    if details.get("equal_edit_alternative_boundary"):
        return "terminal_boundary_provisional"
    if rule in {
        "uniform_edge_substitutions",
        "modal_insertion_alignment",
        "ambiguous_edge_width_competition",
        "retained_polymorphic_edge_trim",
        "retained_allele_flanks",
        "retained_indel_allele_extension",
    }:
        return "repeat_boundary_provisional"
    return "source_supported_proposal"


def apply(dna, categories, parameters, score, source_cds):
    anchors = [
        (cat, key, i, row, geometry(row[1]))
        for cat in (0, 2)
        for key, rows in categories[cat].items()
        for i, row in enumerate(rows)
    ]
    plans = []
    for slot in anchors:
        cat, key, i, parent, old = slot
        if not 3 <= len(old) <= 128 or not _eligible_parent(cat, parent):
            continue
        retained = _retained(categories, key, old)
        retained[old] = float(parent[0])
        choice = _selection(dna, old, retained, parent[1])
        if choice is None:
            continue
        path, rule, details = choice
        span = _proposal_span(choice)
        if any(
            other[:3] != slot[:3] and _overlaps(span, other[-1]) for other in anchors
        ):
            continue
        if _preservation(dna, parent, path, rule) is not None:
            plans.append((slot, span, choice, retained))
    if not plans:
        return categories
    structural = build_advanced_fuzzy_search_filter(parameters)
    output = [{key: list(rows) for key, rows in cat.items()} for cat in categories]
    changes = []
    for index, (slot, span, choice, retained) in enumerate(plans):
        cat, key, i, parent, old = slot
        if any(
            index != j and _overlaps(span, other[1]) for j, other in enumerate(plans)
        ):
            continue
        path, rule, details = choice
        fresh = _score_one(dna, parent, path, key, score)
        cds = _cds_scope(dna, fresh[1], rule, details, source_cds)
        if cds is None:
            continue
        check = _assess(dna, fresh, structural, cds, old, retained, rule)
        category = _category(cat, parent[0], fresh[0], rule)
        if not check.qualified or category is None:
            continue
        candidates = {path: _Assessment(fresh, cds, check)}
        rival_record = None
        if details.get("alternative_repeat_intervals"):
            alt_path = tuple(map(tuple, details["alternative_repeat_intervals"]))
            rival = _score_one(dna, parent, alt_path, key, score)
            alt_cds = _cds_scope(dna, rival[1], "", {}, source_cds)
            alt_check = validate_reconstruction_candidate(
                dna, rival, rival, structural, source_cds=alt_cds
            )
            candidates[alt_path] = _Assessment(rival, alt_cds, alt_check)
            path, rule, details = _rank_choice(dna, parent, cat, choice, candidates)
            fresh, cds, check = candidates[path]
            other_path = tuple(map(tuple, details["alternative_repeat_intervals"]))
            rival_record = candidates[other_path].as_alternative()
            category = _category(cat, parent[0], fresh[0], rule)
        preservation = _preservation(dna, parent, path, rule)
        if preservation is None or category is None:
            continue
        proposal = fresh[1]
        proof = _plain(
            dict(
                method=METHOD,
                rule=rule,
                coordinate_system="zero_based_half_open",
                geometry_status=_uncertainty(rule, details),
                score_semantics="own_fresh_array_score",
                source_sha256=check.evidence["source_sha256"],
                evidence=details,
                previous_category_index=cat,
                previous=_snapshot(parent),
                previous_provenance_sha256=_digest(parent[3]),
                selected=_snapshot(fresh),
                occurrence_preservation=preservation,
                local_cds=[asdict(c) for c in cds],
                alternative=rival_record,
                structural_parameters={k: parameters[k] for k in STRUCTURAL_PARAMETERS},
                validation=check.evidence,
            )
        )
        motif_core = proposal if rule == DUPLICATION else parent[1]
        labels = {
            p: v
            for p, v in getattr(parent[1], "potential_degenerate_repeats", {}).items()
            if p in path
            and (rule != DUPLICATION or v["core_consensus"] == motif_core.consensus)
        }
        labels.update(
            added_member_evidence(
                dna, motif_core, old, path, parameters["param_max_edit_distance"]
            )
        )
        proposal.potential_degenerate_repeats = labels
        changes.append(
            (cat, key, i, category, list(fresh[:3]) + [{"classification_basis": proof}])
        )
    # All locations refer to the original category lists, including multi-row loci.
    for cat, key, i, category, fresh in sorted(changes, key=lambda c: -c[2]):
        parent = output[cat][key][i]
        output[cat + 1].setdefault(key, []).append(parent)
        if cat == category:
            output[cat][key][i] = fresh
            continue
        output[cat][key].pop(i)
        if not output[cat][key]:
            del output[cat][key]
        output[category].setdefault(key, []).append(fresh)
    return output if changes else categories


def validate(dna, basis, current, cat, categories, source_hash):
    if (
        basis["source_sha256"] != source_hash
        or basis["coordinate_system"] != "zero_based_half_open"
    ):
        raise ValueError("competing geometry source differs")
    if _snapshot(current) != basis["selected"]:
        raise ValueError("competing geometry score or intervals differ")
    fresh = _restore(dna, basis["selected"])
    previous = _restore(dna, basis["previous"])
    old, path = geometry(previous[1]), geometry(fresh[1])
    old_cat = basis["previous_category_index"]
    if (
        old_cat not in (0, 2)
        or _category(old_cat, previous[0], fresh[0], basis["rule"]) != cat
    ):
        raise ValueError("invalid competing geometry confidence category")
    locations = [
        k
        for k, rows in categories[cat].items()
        if any(r[1] is current[1] for r in rows)
    ]
    if len(locations) != 1:
        raise ValueError("missing competing geometry locus")
    key = locations[0]
    parent = next(
        (
            r
            for r in categories[old_cat + 1].get(key, ())
            if geometry(r[1]) == old
            and _snapshot(r) == basis["previous"]
            and len(r) == 4
            and _digest(r[3]) == basis["previous_provenance_sha256"]
        ),
        None,
    )
    if parent is None or not _eligible_parent(old_cat, parent):
        raise ValueError("competing geometry lost its accepted parent witness")
    retained = _retained(categories, key, old, exclude=current)
    retained[old] = float(parent[0])
    choice = _selection(dna, old, retained, parent[1])
    if choice is None or choice[1] != basis["rule"]:
        raise ValueError("competing geometry rule or placement failed replay")
    preservation = _preservation(dna, parent, path, basis["rule"])
    if preservation is None or preservation != basis["occurrence_preservation"]:
        raise ValueError("competing geometry occurrence evidence differs")
    span = _proposal_span(choice)
    if any(
        r[1] is not current[1] and _overlaps(span, geometry(r[1]))
        for idx in (0, 2)
        for rows in categories[idx].values()
        for r in rows
    ):
        raise ValueError("competing geometry overlaps another accepted array")
    cds = tuple(SourceCDS(**c) for c in basis["local_cds"])
    if _cds_scope(dna, fresh[1], basis["rule"], basis["evidence"], lambda: cds) != cds:
        raise ValueError("competing geometry coding context differs")
    structural = build_advanced_fuzzy_search_filter(basis["structural_parameters"])
    check = _assess(dna, fresh, structural, cds, old, retained, basis["rule"])
    if not check.qualified or _plain(check.evidence) != basis["validation"]:
        raise ValueError("competing geometry negative or structural checks differ")
    candidates = {path: _Assessment(fresh, cds, check)}
    alternative = basis.get("alternative")
    raw_alt = choice[2].get("alternative_repeat_intervals")
    if raw_alt is not None:
        if not alternative:
            raise ValueError("missing competing alternative")
        rival = _restore(dna, alternative["candidate"])
        rival_path = geometry(rival[1])
        if {path, rival_path} != {choice[0], tuple(map(tuple, raw_alt))} or list(
            rival[1].list_repeats_gaped
        ) != list(_source_repeat_alignment(rival[1].list_repeats)):
            raise ValueError("competing alternative intervals differ")
        coding = tuple(SourceCDS(**c) for c in alternative["local_cds"])
        if _cds_scope(dna, rival[1], "", {}, lambda: coding) != coding:
            raise ValueError("competing alternative coding context differs")
        assessed = validate_reconstruction_candidate(
            dna, rival, rival, structural, source_cds=coding
        )
        if (
            alternative["status"] != "unselected"
            or alternative["qualified"] is not assessed.qualified
            or alternative["validation"] != _plain(assessed.evidence)
        ):
            raise ValueError("competing alternative validation differs")
        candidates[rival_path] = _Assessment(rival, coding, assessed)
        choice = _rank_choice(dna, parent, old_cat, choice, candidates)
        if tuple(map(tuple, choice[2]["alternative_repeat_intervals"])) != rival_path:
            raise ValueError("competing alternative role differs")
    elif alternative is not None:
        raise ValueError("unexpected competing alternative")
    if _plain(choice) != _plain((path, basis["rule"], basis["evidence"])):
        raise ValueError("competing geometry rule or placement failed replay")
    if (
        basis["geometry_status"] != _uncertainty(choice[1], choice[2])
        or basis["score_semantics"] != "own_fresh_array_score"
    ):
        raise ValueError("competing geometry confidence semantics differ")
    return basis
