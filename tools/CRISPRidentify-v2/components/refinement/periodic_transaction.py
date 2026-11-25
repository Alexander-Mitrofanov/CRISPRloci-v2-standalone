"""Freshly assess locally periodic extensions before committing reconstruction seeds."""

from bisect import bisect_left
from collections import Counter

from components.evidence.candidate_origins import inherit_candidate_origins
from components.evidence.initial_competition import guard_reconstructed_occurrences
from components.evidence.occurrence_protection import occurrence_relation
from components.reconstruction.reconstruction_validation import (
    needs_coding_context,
    validate_reconstruction_candidate,
)
from components.reconstruction.reconstruction_transaction import (
    STRUCTURAL_PARAMETERS,
    _restore,
    _snapshot,
)
from components.refinement.boundary_panel import resolve_phase
from components.refinement.enhancement_search import (
    build_source_backed_candidate_from_boundaries,
)
from components.refinement.periodic_completion import _fit_width, _spacing, plan_path
from components.refinement.potential_degenerate import added_member_evidence, geometry
from components.refinement.repeat_boundary_repair import classification_restriction
from components.scoring.evaluation import build_advanced_fuzzy_search_filter
from components.rescue.architecture_rescue_coding_context import SourceCDS
from components.rescue.likely_tandem_transaction import validated_row_tags


def _covers(longer, shorter, tolerance):
    """Require a distinct, same-phase source copy for each shorter-path repeat."""
    j = 0
    for left, right in shorter:
        j = bisect_left(longer, (left - tolerance,), lo=j)
        if (
            j == len(longer)
            or abs(longer[j][0] - left) > tolerance
            or abs(longer[j][1] - right) > tolerance
        ):
            return False
        j += 1
    return True


def _overlaps(a, b):
    return a[0][0] < b[-1][1] and b[0][0] < a[-1][1]


def _retained(cats, key):
    return (
        r
        for cat in cats
        for r in cat.get(key, ())
        if not classification_restriction(r[1])
    )


def _longer_phase(dna, path, candidates):
    # Sparse arrays with unequal copy widths may also support a longer
    # repeated sequence. Two exact copies and the same third occurrence
    # justify disclosing that alternative when model support is weak.
    if len(path) != 3 or len({b - a for a, b in path}) == 1:
        return None
    widths = sorted(b - a for a, b in path)
    supported = set()
    for row in candidates:
        other = geometry(row[1])
        if len(other) != len(path) or not _covers(other, path, 7):
            continue
        sizes = {b - a for a, b in other}
        if len(sizes) != 1 or not widths[1] + 4 <= min(sizes) <= widths[1] + 7:
            continue
        counts = Counter(dna[a:b] for a, b in other)
        motif, count = counts.most_common(1)[0]
        if count >= 2 and all(
            sum(x != y for x, y in zip(seq, motif)) <= len(motif) // 5 for seq in counts
        ):
            supported.add(other)
    return supported.pop() if len(supported) == 1 else None


def _embedded_fragment_donor(dna, fragment, alternatives):
    """Find a single periodic path whose spacers explain the long repeat edges."""
    old = geometry(fragment)
    if len(old) != 3 or len(set(fragment.list_repeats)) != 1 or classification_restriction(fragment):
        return None
    supported = {}
    for row in alternatives:
        donor = row[1]
        if (row[0] < 0.75 or len(donor.list_repeats) < 3 * len(old)
                or classification_restriction(donor) or "tandem" in validated_row_tags(row)):
            continue
        path = geometry(donor)
        if path[0][1] >= old[0][0] or path[-1][0] <= old[-1][1]:
            continue
        motif, count = Counter(donor.list_repeats).most_common(1)[0]
        width = len(motif)
        if (2 * count <= len(path) or not 2 * width >= len(fragment.list_repeats[0]) >= 1.5 * width
                or any(b - a != width for a, b in path[1:-1])):
            continue
        path = _fit_width(dna, motif, path, include_ends=True)
        spacing = _spacing(path)
        if spacing is None or any(b - a != width for a, b in path):
            continue
        spacers = [dna[b:c] for (_, b), (c, _) in zip(path, path[1:])]
        typical, margin = spacing
        if (any(abs(len(s) - typical) > margin for s in spacers)
                or 4 * len(set(spacers)) < 3 * len(spacers)):
            continue
        relation = occurrence_relation(old, path)
        if relation["status"] != "unique":
            continue
        # Require the same anchored edge on every fragment repeat. Containment
        # alone would allow an unrelated short repeat anywhere in the interval.
        mapped = [path[i] for i in relation["mapping"]]
        if not (all(a == x for (a, _), (x, _) in zip(old, mapped))
                or all(b == y for (_, b), (_, y) in zip(old, mapped))):
            continue
        if path not in supported or row[0] > supported[path][0]:
            supported[path] = row
    if not supported:
        return None
    # Local search can retain multiple truncations of one fully observed path.
    # Require every shorter proposal to agree at every shared occurrence.
    full = max(supported, key=len)
    for path in supported:
        start = bisect_left(full, path[0])
        if full[start:start + len(path)] != path:
            return None
    return full, supported[full]


def _compatible_anchors(path, old, core, collisions, reconstructed):
    """Every displaced core must map, with a retained bridge across disjoint cores."""
    if any(
        not _covers(path, p, min(7, max(len(core.consensus), len(r[1].consensus)) // 4))
        or classification_restriction(r[1])
        or (
            len(r) == 4
            and r[3].get("classification_basis", {}).get("method")
            == "source_backed_boundary_uncertain_reconstruction"
        )
        for _, r, p in collisions
    ):
        return False
    if any(not _overlaps(old, p) for _, _, p in collisions):
        witnesses = (
            geometry(w[1])
            for (_, locus, _), _, _ in collisions
            for w in _retained(reconstructed, locus)
        )
        if not any(
            all(_covers(witness, p, 4) for _, _, p in collisions)
            for witness in witnesses
        ):
            return False
    return True


def _retained_alignment(reconstructed, key, path):
    """Select a single retained indel alignment on the proposed copy positions."""
    donors = [
        r
        for r in _retained(reconstructed, key)
        if r[0] >= 0.5
        and len(r[1].list_repeats) == len(path)
        and _covers(geometry(r[1]), path, 2)
    ]
    if not donors:
        return None
    best = max(donors, key=lambda r: r[0])
    if any(r[0] == best[0] and geometry(r[1]) != geometry(best[1]) for r in donors):
        return None
    return best


def refine(dna, original, reconstructed, parameters, score, source_cds):
    ordinary = parameters["param_max_edit_distance"]
    anchors = {
        (cat, key, n): (row, geometry(row[1]))
        for cat in (0, 2)
        for key, rows in reconstructed[cat].items()
        for n, row in enumerate(rows)
    }
    seeds = [{key: list(rows) for key, rows in cat.items()} for cat in original]
    replacements, hypotheses, reserved = {}, {}, []
    structural = build_advanced_fuzzy_search_filter(parameters)

    parent_witnesses = []
    motifs = parameters.get("param_effective_config", {}).get(
        "boundary_consensuses", ()
    )

    def assess(candidate, *, panel_phase=False, allele_trim_parent=None):
        rows = score({"local-period": [candidate]})
        if set(rows) != {"local-period"} or len(rows["local-period"]) != 1:
            raise ValueError("periodic refinement scoring changed request cardinality")
        fresh = rows["local-period"][0]
        if fresh[1] is not candidate:
            raise ValueError("periodic refinement scoring changed candidate identity")
        cds = tuple(source_cds()) if needs_coding_context(candidate) else None
        check = validate_reconstruction_candidate(
            dna, fresh, fresh, structural, source_cds=cds, panel_phase=panel_phase,
            allele_trim_parent=allele_trim_parent,
        )
        return fresh, check

    for slot in sorted(anchors, key=lambda loc: (-len(anchors[loc][1]), loc)):
        if slot in replacements:
            continue
        row, old = anchors[slot]
        basis = row[3].get("classification_basis", {}) if len(row) == 4 else {}
        if basis.get("method") == "source_backed_boundary_uncertain_reconstruction":
            continue
        core = row[1]
        key = slot[1]
        path = plan_path(dna, core, ordinary)
        phase = resolve_phase(dna, path, motifs) if len(path) > len(old) else None
        if path == old:
            path = _longer_phase(dna, old, _retained(reconstructed, key))
            if not path:
                retained = _embedded_fragment_donor(dna, core, reconstructed[1].get(key, ()))
                if retained is None:
                    continue
                full, donor = retained
                if (any(_overlaps(full, prior) for prior in reserved)
                        or any(_overlaps(full, other) for loc, (_, other) in anchors.items() if loc != slot)
                        or not guard_reconstructed_occurrences(full, donor, row, original[slot[0]][key][slot[2]], source=dna)):
                    continue
                candidate = build_source_backed_candidate_from_boundaries(dna, donor[1], *zip(*full))
                inherit_candidate_origins(candidate, (core, donor[1]))
                fresh, check = assess(candidate)
                if check.qualified and fresh[0] >= 0.75:
                    candidate.periodic_refinement = dict(
                        method="retained_array_spanning_fragment", status="source_supported_geometry",
                        coordinate_system="zero_based_half_open", parent_repeat_intervals=[old, geometry(donor[1])],
                        parent_model_scores=[float(row[0]), float(donor[0])],
                    )
                    replacements[slot] = (0, key, fresh)
                    reserved.append(full)
                continue
        if any(_overlaps(path, prior) for prior in reserved):
            continue
        collisions = [
            (loc, r, p) for loc, (r, p) in anchors.items() if _overlaps(path, p)
        ]
        # A periodic extension can subsume independent core fragments only
        # when a retained observation spans their junction in the same phase.
        if not _compatible_anchors(path, old, core, collisions, reconstructed):
            continue
        parents = [r[1] for _, r, _ in collisions]
        pre_rows = [original[c][k][n] for (c, k, n), _, _ in collisions]
        if not guard_reconstructed_occurrences(
            path, *pre_rows, *(r for _, r, _ in collisions), source=dna
        ):
            continue
        candidate = build_source_backed_candidate_from_boundaries(
            dna, core, *zip(*path)
        )
        inherit_candidate_origins(candidate, parents)
        fresh, check = assess(candidate, allele_trim_parent=old
                              if len(collisions) == 1 and row[0] >= 0.5 else None)
        allele_trim = bool(check.qualified and "recurrent_terminal_alleles" in check.evidence)
        if allele_trim and fresh[0] < 0.5:
            continue
        provisional = (
            fresh[0] < 0.5
            and check.qualified
            and len(collisions) == 1
            and len(old) == 3
            and row[0] >= 0.5
            and (len(path) == 3 or len(path) == 4 and set(old) <= set(path))
        )
        if not provisional and (not check.qualified or fresh[0] < 0.5):
            hypotheses[old] = dict(
                status="unaccepted_alternative",
                model_score=float(fresh[0]),
                coordinate_system="zero_based_half_open",
                repeat_intervals=path,
                assessment=check.evidence["decision"],
                failed_predicates=[
                    name
                    for name, status in structural.evidence(candidate).predicate_results
                    if status == "failed"
                ],
            )
            # A longer retained alignment may have small genuine indels where
            # the fixed-width proposal failed. Keep its source geometry intact.
            if len(path) < len(old) + 2:
                continue
            best = _retained_alignment(reconstructed, key, path)
            if best is None:
                continue
            path = geometry(best[1])
            if not guard_reconstructed_occurrences(path, *pre_rows, best, source=dna):
                continue
            candidate = build_source_backed_candidate_from_boundaries(
                dna, core, *zip(*path)
            )
            inherit_candidate_origins(candidate, (*parents, best[1]))
            fresh, check = assess(candidate)
            if not check.qualified or fresh[0] < 0.5:
                continue
        curated_selected = False
        if phase and not provisional:
            projected, references = phase
            if _covers(projected, path, 3) and guard_reconstructed_occurrences(
                projected, *pre_rows, *(r for _, r, _ in collisions), source=dna
            ):
                curated = build_source_backed_candidate_from_boundaries(
                    dna, candidate, *zip(*projected)
                )
                curated_row, curated_check = assess(curated, panel_phase=True)
                if curated_check.qualified:
                    curated_selected = True
                    parent_witnesses.append((slot[1], fresh))
                    proof = {
                        "method": "curated_periodic_boundary",
                        "status": "provisional",
                        "core": _snapshot(fresh),
                        "geometry": _snapshot(curated_row),
                        "search_core": _snapshot(row),
                        "ordinary_edit_limit": ordinary,
                        "matched_consensuses": references,
                        "source_sha256": curated_check.evidence["source_sha256"],
                        "structural_parameters": {
                            k: parameters[k] for k in STRUCTURAL_PARAMETERS
                        },
                        "validation": curated_check.evidence,
                    }
                    path, candidate = projected, curated
                    fresh = list(curated_row[:3]) + [{"classification_basis": proof}]
        # A fallback alignment or curated phase can move an endpoint into
        # another locus. Recheck the final geometry before retiring any core.
        collided_slots = {loc for loc, _, _ in collisions}
        if any(_overlaps(path, prior) for prior in reserved) or any(
            _overlaps(path, other)
            for loc, (_, other) in anchors.items()
            if loc not in collided_slots
        ):
            if curated_selected:
                parent_witnesses.pop()
            continue
        reserved.append(path)
        for loc, _, _ in collisions:
            replacements[loc] = None
        category = 0 if fresh[0] >= 0.75 and not (curated_selected or allele_trim) else 2
        if provisional or allele_trim:
            if not any(r[1] is core for r in pre_rows):
                parent_witnesses.append((key, row[:3]))
            proof = {
                "method": "recurrent_allele_boundary" if allele_trim else "core_supported_periodic_geometry",
                "status": "provisional",
                "core": _snapshot(row),
                "geometry": _snapshot(fresh),
                "source_sha256": check.evidence["source_sha256"],
                "structural_parameters": {k: parameters[k] for k in STRUCTURAL_PARAMETERS},
                "ordinary_edit_limit": ordinary,
                "validation": check.evidence,
            }
            if not allele_trim:
                proof["retained_phase"] = len(path) == len(old)
            fresh = list(fresh[:3]) + [{"classification_basis": proof}]
        replacements[slot] = (category, key, fresh)
        candidate.periodic_refinement = dict(
            coordinate_system="zero_based_half_open",
            parent_repeat_intervals=[p for _, _, p in collisions],
            parent_model_scores=[float(r[0]) for _, r, _ in collisions],
            method="local_spacer_period",
            status="source_supported_geometry",
        )
        candidate.potential_degenerate_repeats = added_member_evidence(
            dna, core, old, path, ordinary
        )
    if replacements:
        for cat, key, n in sorted(replacements, reverse=True):
            prior = seeds[cat][key].pop(n)
            if not seeds[cat][key]:
                del seeds[cat][key]
            seeds[cat + 1].setdefault(key, []).append(prior)
        for key, witness in parent_witnesses:
            seeds[3].setdefault(key, []).append(witness)
        for replacement in replacements.values():
            if replacement:
                cat, key, fresh = replacement
                seeds[cat].setdefault(key, []).append(fresh)
    return seeds if replacements else None, hypotheses


def validate_provisional(
    dna, basis, current_row, category, source_hash, categories, motifs=()
):
    """Rebind the accepted core, fresh scores, independently sliced path and negative checks."""
    curated = basis["method"] == "curated_periodic_boundary"
    allele_trim = basis["method"] == "recurrent_allele_boundary"
    if (
        category != 2
        or basis["status"] != "provisional"
        or basis["source_sha256"] != source_hash
        or _snapshot(current_row) != basis["geometry"]
        or (not (curated or allele_trim) and not 0 <= current_row[0] < 0.5)
        or (allele_trim and not 0.5 <= current_row[0] <= 1)
    ):
        raise ValueError("provisional geometry or category is inconsistent")
    parent = _restore(dna, basis["core"])
    fresh = _restore(dna, basis["geometry"])
    path, old = geometry(fresh[1]), geometry(parent[1])
    if (not (curated or allele_trim) and len(old) != 3) or parent[0] < 0.5:
        raise ValueError("provisional geometry lacks an accepted core")
    rows = [r for cat in categories for group in cat.values() for r in group]
    if not any(
        r[0] == parent[0] and geometry(r[1]) == old and _snapshot(r) == basis["core"]
        for r in rows
    ):
        raise ValueError("provisional geometry lost its core witness")
    if curated:
        search_core = _restore(dna, basis["search_core"])[1]
        original = plan_path(dna, search_core, basis["ordinary_edit_limit"])
        result = resolve_phase(dna, original, motifs)
        if (
            not result
            or not _covers(path, old, 3)
            or len(old) != len(path)
            or len(original) <= len(geometry(search_core))
            or result[1] != basis["matched_consensuses"]
        ):
            raise ValueError("curated boundary lacks a unique supported panel phase")
        expected = result[0]
    elif allele_trim:
        expected = plan_path(dna, parent[1], basis["ordinary_edit_limit"])
    elif basis["retained_phase"]:
        retained = [
            r for r in rows if r[1] is not current_row[1] and geometry(r[1]) == path
        ]
        expected = _longer_phase(dna, old, retained)
    else:
        expected = plan_path(dna, parent[1], basis["ordinary_edit_limit"])
        if len(path) != 4 or not set(old) <= set(path):
            raise ValueError("provisional extension lost its anchor copies")
    if expected != path:
        raise ValueError("provisional source placement could not be reproduced")
    structural = build_advanced_fuzzy_search_filter(basis["structural_parameters"])
    coding = basis["validation"].get("coding_tandem", {})
    cds = tuple(SourceCDS(**entry) for entry in coding.get("considered_source_cds", ()))
    checked = validate_reconstruction_candidate(
        dna, fresh, fresh, structural, source_cds=cds, panel_phase=curated,
        allele_trim_parent=old if allele_trim else None,
    )
    if allele_trim and "recurrent_terminal_alleles" not in checked.evidence:
        raise ValueError("provisional trim lacks recurrent source alleles")
    if not checked.qualified or checked.evidence != basis["validation"]:
        raise ValueError("provisional source and structural checks differ")
    return basis
