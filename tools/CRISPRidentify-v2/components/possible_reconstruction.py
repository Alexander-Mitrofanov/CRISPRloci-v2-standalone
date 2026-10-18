"""Checked retained extents of Possible cores, without upgrading confidence.

This separate transaction runs after the unchanged Bona-fide transaction.
It reuses direct relationship, fresh candidate checks and selection mechanics,
but admits only PD donors, preserves the Possible tier, and abstains on distinct
geometry ties. No source discovery, classifier threshold change or BF evidence
inheritance occurs here.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import math

from components.architecture_rescue_coding_context import SourceCDS
from components.candidate_origins import candidate_origin_evidence
from components.fast_run_rescue import compact_feature_values
from components.initial_competition import (
    METADATA_KEY, competition_evidence, guard_reconstructed_occurrences,
    ordinary_reference_key, protected_ordinary_reference,
    validate_ordinary_witness_inventory,
)
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.reconstruction_reference import MAPPED_CORE_REFERENCE, recorded_reference_rule
from components.reconstruction_selection import (
    RowLocation, _prepare, reconstruction_candidate_locations, select_reconstructions,
)
from components.reconstruction_transaction import (
    STRUCTURAL_PARAMETERS, _digest, _encoded, _inventory_state, _plain,
    _restore, _snapshot,
)
from components.reconstruction_validation import needs_coding_context, validate_reconstruction_candidate


METHOD = "source_backed_possible_reconstruction"
EVIDENCE_FILENAME = "possible_reconstruction_evidence.json"
RULE = "possible-core-reconstruction-v1"


@dataclass(frozen=True)
class PossibleReconstructionPolicy:
    """A distinct tier contract; the established BF policy is not weakened."""

    minimum_anchor_score: float = 0.50
    minimum_candidate_score: float = 0.60
    boundary_tolerance: int = 3
    minimum_reciprocal_overlap: float = 0.80
    maximum_spacer_length: int = 140
    minimum_added_repeat_identity: float = 0.80
    score_tie_tolerance: float = 1e-12

    def as_dict(self):
        return asdict(self)


POLICY = PossibleReconstructionPolicy()


def _decision(decision):
    value = decision.as_dict()
    value["rule"] = RULE
    return _plain(value)


def _observation(row, source_hash):
    return {"snapshot": _snapshot(row),
            "metadata": _plain(row[3] if len(row) == 4 else {}),
            "origins": candidate_origin_evidence(row[1], source_hash)}


def _restore_typed(source, snapshot, metadata=None):
    row = _restore(source, snapshot)
    if _encoded(_snapshot(row)) != _encoded(snapshot):
        raise ValueError("Possible reconstruction snapshot types differ")
    if metadata:
        if not isinstance(metadata, dict):
            raise ValueError("Possible original metadata must be a mapping")
        row.append(deepcopy(metadata))
    return row


def reconstruct_possible_arrays(source, categories, parameters, score_candidates, get_source_cds):
    """Return copied categories and a separate, source-bound Possible audit."""
    structural_parameters = {key: parameters[key] for key in STRUCTURAL_PARAMETERS}
    structural = build_advanced_fuzzy_search_filter(structural_parameters)
    validate_ordinary_witness_inventory(categories, source)
    before = _inventory_state(categories)
    prepared = _prepare(categories, POLICY, MAPPED_CORE_REFERENCE, _possible=True)
    donors = reconstruction_candidate_locations(categories, POLICY, _possible=True, _prepared=prepared)
    rows = [location.row for location in donors]
    rows.extend(anchor.row for anchor, proposals, _rejected in prepared[1] if proposals)
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    observations = {id(row): _digest(_observation(row, source_hash)) for row in rows}
    candidates = list({id(row[1]): row[1] for row in rows}.values())
    key = "possible-reconstruction-fresh-score-batch"
    scored = score_candidates({key: candidates}) if candidates else {}
    if _inventory_state(categories) != before:
        raise ValueError("Possible fresh scoring mutated the original inventory")
    if any(_digest(_observation(row, source_hash)) != observations[id(row)] for row in rows):
        raise ValueError("Possible fresh scoring mutated original provenance")
    if candidates and (set(scored) != {key} or len(scored[key]) != len(candidates)):
        raise ValueError("Possible fresh scoring changed request cardinality")
    fresh = {}
    for candidate, row in zip(candidates, scored.get(key, ())):
        if row[1] is not candidate:
            raise ValueError("Possible fresh scoring changed identity/order")
        if not math.isfinite(float(row[0])) or not 0 <= row[0] <= 1:
            raise ValueError("Possible fresh score is invalid")
        compact_feature_values(row)
        fresh[id(candidate)] = row
    cds = tuple(get_source_cds()) if any(needs_coding_context(c) for c in candidates) else None
    validations = {}

    def validate(row):
        if id(row) not in validations:
            if id(row[1]) not in fresh:
                raise ValueError("Possible selection requested an unscored candidate")
            validations[id(row)] = validate_reconstruction_candidate(
                source, row, fresh[id(row[1])], structural, source_cds=cds)
        return validations[id(row)]

    result = select_reconstructions(categories, validate, POLICY, _possible=True, _prepared=prepared)
    if _inventory_state(categories) != before:
        raise ValueError("Possible validation mutated the original inventory")
    if any(_digest(_observation(row, source_hash)) != observations[id(row)] for row in rows):
        raise ValueError("Possible validation mutated original provenance")
    output = [{key: list(rows) for key, rows in category.items()} for category in categories]
    selected = [d for d in result.decisions if d.action == "selected"]
    remove = {(d.selected.locus_key, d.selected.row_index) for d in selected}
    for key, rows in list(output[3].items()):
        retained = [row for ordinal, row in enumerate(rows) if (key, ordinal) not in remove]
        if retained:
            output[3][key] = retained
        else:
            del output[3][key]
    selected_records = {}
    for decision in selected:
        anchor, donor = decision.anchor, decision.selected
        core = validations[id(anchor.row)]
        checked = validations[id(donor.row)]
        anchor_observation = _observation(anchor.row, source_hash)
        donor_observation = _observation(donor.row, source_hash)
        selection = dict(_decision(decision), original_observations_sha256={
            "anchor": _digest(anchor_observation), "donor": _digest(donor_observation)})
        selected_records[id(decision)] = selection
        proof = _plain({
            "rule": RULE, "source_sha256": source_hash,
            "coordinate_system": "zero_based_half_open",
            "policy": POLICY.as_dict(), "structural_parameters": structural_parameters,
            "confidence_category": "Possible",
            "anchor_original": _snapshot(anchor.row), "anchor_fresh": _snapshot(core.row),
            "donor_original": _snapshot(donor.row), "reconstruction": _snapshot(checked.row),
            "anchor_observation": anchor_observation,
            "donor_observation": donor_observation,
            "anchor_initial_competition": competition_evidence(anchor.row, len(source), source=source),
            "donor_initial_competition": competition_evidence(donor.row, len(source), source=source),
            "anchor_validation": core.evidence, "reconstruction_validation": checked.evidence,
            "selection": selection,
            "score_semantics": "own_fresh_model_score; original_possible_core_confidence_is_separate",
        })
        metadata = deepcopy(donor.row[3]) if len(donor.row) == 4 else {}
        metadata.pop(METADATA_KEY, None)
        metadata.pop("boundary_assessment", None)
        metadata["classification_basis"] = {"method": METHOD, "reconstruction": proof,
                                            "evidence_sha256": _digest(proof)}
        output[2][anchor.locus_key][anchor.row_index] = list(checked.row[:3]) + [metadata]
        output[3].setdefault(anchor.locus_key, []).append(anchor.row)
    if sum(map(len, output[2].values())) != sum(map(len, categories[2].values())):
        raise ValueError("Possible reconstruction changed accepted locus count")
    audit = {"schema": {"name": "CRISPRidentify-possible-reconstruction-evidence", "version": "1.0.0"},
             "rule": RULE, "source_sha256": source_hash, "coordinate_system": "zero_based_half_open",
             "policy": POLICY.as_dict(), "confidence_category": "Possible",
             "biological_accuracy_established": False,
             "stats": dict(result.stats, fresh_candidates_scored=len(candidates),
                           candidates_validated=len(validations), selections_committed=len(selected)),
             "actions": dict(Counter(d.action for d in result.decisions)),
             "decisions": [selected_records.get(id(d), _decision(d)) for d in result.decisions]}
    return output, _plain(audit)


_PROOF_KEYS = frozenset((
    "rule", "source_sha256", "coordinate_system", "policy", "structural_parameters",
    "confidence_category", "anchor_original", "anchor_fresh", "donor_original", "reconstruction",
    "anchor_observation", "donor_observation", "anchor_initial_competition",
    "donor_initial_competition", "anchor_validation", "reconstruction_validation",
    "selection", "score_semantics",
))


def _proof(basis):
    if not isinstance(basis, dict) or basis.get("method") != METHOD:
        return None
    if "accepted_family_tandem_review" in basis:
        if not isinstance(basis["accepted_family_tandem_review"], dict):
            raise ValueError("Possible tandem wrapper is malformed")
        basis = {key: value for key, value in basis.items() if key != "accepted_family_tandem_review"}
    if set(basis) != {"method", "reconstruction", "evidence_sha256"}:
        raise ValueError("Possible reconstruction basis keys differ")
    proof = basis["reconstruction"]
    if not isinstance(proof, dict) or set(proof) != _PROOF_KEYS:
        raise ValueError("Possible reconstruction proof keys differ")
    if (basis["evidence_sha256"] != _digest(proof) or proof["rule"] != RULE
            or _encoded(proof["policy"]) != _encoded(POLICY.as_dict())
            or proof["coordinate_system"] != "zero_based_half_open"
            or proof["confidence_category"] != "Possible"
            or proof["score_semantics"] != "own_fresh_model_score; original_possible_core_confidence_is_separate"):
        raise ValueError("Possible reconstruction rule/policy/confidence binding differs")
    expected = {name: _digest(proof[name + "_observation"]) for name in ("anchor", "donor")}
    if proof["selection"].get("original_observations_sha256") != expected:
        raise ValueError("Possible original observations differ from selection")
    return proof


def validate_possible_reconstruction_basis(source, candidate, basis, row_score, category_index,
                                           current_row, categories, *, ordinary_witness_keys=None):
    """Replay the explicit Possible proof using current retained observations."""
    proof = _proof(basis)
    if proof is None or category_index != 2:
        raise ValueError("Possible reconstruction must have a Possible core category")
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    if source_hash != proof["source_sha256"]:
        raise ValueError("Possible reconstruction source binding differs")
    if current_row is None or current_row[1] is not candidate:
        raise ValueError("Possible reconstruction current row is missing")
    if (_encoded(_snapshot(current_row)) != _encoded(proof["reconstruction"])
            or float(row_score) != proof["reconstruction"]["score"]):
        raise ValueError("Possible reconstruction current score/representation differs")
    original_rows, fresh_rows = {}, {}
    for name, fresh_key in (("anchor", "anchor_fresh"), ("donor", "reconstruction")):
        observation = proof[name + "_observation"]
        if not isinstance(observation, dict) or set(observation) != {"snapshot", "metadata", "origins"}:
            raise ValueError("Possible original observation keys differ")
        snapshot = proof[name + "_original"]
        if _encoded(observation["snapshot"]) != _encoded(snapshot):
            raise ValueError("Possible original observation snapshot differs")
        original = _restore_typed(source, snapshot, observation["metadata"])
        fresh = _restore_typed(source, proof[fresh_key])
        for key in ("repeat_intervals", "gapped_repeats", "eligible", "eligibility_reason"):
            if _encoded(snapshot[key]) != _encoded(proof[fresh_key][key]):
                raise ValueError("Possible fresh scoring changed geometry/eligibility")
        recorded = competition_evidence(original, len(source), source=source)
        if _encoded(recorded) != _encoded(proof[name + "_initial_competition"]):
            raise ValueError("Possible original competition differs")
        original_rows[name], fresh_rows[name] = original, fresh
    retained = [row for rows in categories[3].values() for row in rows
                if _encoded(_snapshot(row)) == _encoded(proof["anchor_original"])]
    if not any(_encoded(_observation(row, source_hash)) == _encoded(proof["anchor_observation"])
               for row in retained):
        raise ValueError("Possible reconstruction has no retained exact original core")
    if (_encoded(candidate_origin_evidence(candidate, source_hash))
            != _encoded(proof["donor_observation"]["origins"])):
        raise ValueError("Possible selected donor origins differ")
    if ordinary_witness_keys is None:
        ordinary_witness_keys = validate_ordinary_witness_inventory(categories, source, observations=(
            proof["anchor_initial_competition"], proof["donor_initial_competition"]))
    for row in original_rows.values():
        reference = protected_ordinary_reference(row, source=source)
        if reference is not None and ordinary_reference_key(reference) not in ordinary_witness_keys:
            raise ValueError("Possible reconstruction has no original ordinary witness")
    if not guard_reconstructed_occurrences(proof["reconstruction"]["repeat_intervals"],
                                          *original_rows.values(), source=source):
        raise ValueError("Possible reconstruction violates original occurrence protection")
    structural = build_advanced_fuzzy_search_filter(proof["structural_parameters"])
    validations = {}
    for name, evidence_key in (("anchor", "anchor_validation"), ("donor", "reconstruction_validation")):
        recorded = proof[evidence_key]
        coding = recorded["coding_tandem"]
        cds = tuple(SourceCDS(**raw) for raw in coding["considered_source_cds"]) if coding["source_cds_assessed"] else None
        # Validation requires original/fresh rows to share one candidate object.
        original, fresh = original_rows[name], fresh_rows[name]
        fresh[1] = original[1]
        check = validate_reconstruction_candidate(source, original, fresh, structural, source_cds=cds)
        if not check.qualified or _encoded(check.evidence) != _encoded(recorded):
            raise ValueError("Possible candidate-specific validation differs")
        validations[id(original)] = check
    trial = [{key: list(rows) for key, rows in category.items()} for category in categories]
    trial[2] = {key: [row for row in rows if row[1] is not candidate] for key, rows in trial[2].items()}
    # Only the recorded core donates authority during replay; other accepted
    # Possible rows are collision witnesses, not additional reconstruction jobs.
    trial[2][("possible-replay-core",)] = [original_rows["anchor"]]
    trial[3] = {("possible-replay-donor",): [original_rows["donor"]]}
    reference_rule = recorded_reference_rule(proof["selection"]["relationship"])
    result = select_reconstructions(trial, lambda row: validations[id(row)], POLICY,
                                    reference_rule=reference_rule, _possible=True)
    decision = next(d for d in result.decisions if d.anchor.row is original_rows["anchor"])
    selection = proof["selection"]
    if (decision.action != "selected" or selection["action"] != "selected"
            or selection["rule"] != RULE or selection["boundary_uncertain"] is not False
            or _encoded(decision.relation) != _encoded(selection["relationship"])
            or selection["reconstruction_score"] != row_score):
        raise ValueError("Possible relationship/anti-fusion selection replay differs")
    for name, snapshot, category in (("anchor", proof["anchor_original"], 2),
                                     ("selected", proof["donor_original"], 3)):
        location = selection[name]
        original = original_rows["anchor" if name == "anchor" else "donor"]
        actual = RowLocation(category, location["locus_key"], location["row_index"], original).as_dict()
        if (type(location["category_index"]) is not int or location["category_index"] != category
                or type(location["row_index"]) is not int or location["row_index"] < 0
                or _encoded(location) != _encoded(actual)):
            raise ValueError("Possible selection original location differs")
    return basis


def _summary(basis):
    proof = _proof(basis)
    if proof is None:
        return None
    return {"method": METHOD, "confidence_category": "Possible",
            "core_repeat_count": len(proof["anchor_original"]["repeat_intervals"]),
            "core_original_own_score": proof["anchor_original"]["score"],
            "core_fresh_own_score": proof["anchor_fresh"]["score"],
            "repeat_count": len(proof["reconstruction"]["repeat_intervals"]),
            "own_score": proof["reconstruction"]["score"],
            "added_repeat_count": len(proof["selection"]["relationship"]["added_repeat_indexes"]),
            "evidence_sha256": _digest(proof)}


def possible_reconstruction_summary(record):
    """Compact annotation for an already fully validated scored row."""
    basis = record[3].get("classification_basis") if len(record) == 4 else None
    result = _summary(basis)
    if result is not None and _encoded(_snapshot(record)) != _encoded(basis["reconstruction"]["reconstruction"]):
        raise ValueError("Possible summary row differs from its proof")
    return result


def report_possible_reconstruction_summary(array):
    result = _summary(array.get("classification_basis"))
    if result is None:
        if "possible_reconstruction" in array:
            raise ValueError("Possible summary has no matching proof")
        return None
    proof = array["classification_basis"]["reconstruction"]
    if any(type(repeat["source_interval"][key]) is not int
           for repeat in array["repeats"] for key in ("start", "end")):
        raise ValueError("Possible canonical repeat coordinates have invalid types")
    intervals = [[repeat["source_interval"]["start"] - 1, repeat["source_interval"]["end"]]
                 for repeat in array["repeats"]]
    if (_encoded(array.get("possible_reconstruction")) != _encoded(result) or array["category"] != "Possible"
            or _encoded(intervals) != _encoded(proof["reconstruction"]["repeat_intervals"])
            or _encoded(array["certainty_score"]) != _encoded(result["own_score"])
            or _encoded(array["repeat_count"]) != _encoded(result["repeat_count"])
            or _encoded(array["source_interval"]["start"]) != _encoded(intervals[0][0] + 1)
            or _encoded(array["source_interval"]["end"]) != _encoded(intervals[-1][1])):
        raise ValueError("Possible canonical summary binding differs")
    return result


def _nested_proofs(basis):
    if not isinstance(basis, dict):
        return
    proof = _proof(basis)
    if proof is not None:
        yield proof
    elif "prior_classification_basis" in basis:
        yield from _nested_proofs(basis["prior_classification_basis"])


def validate_possible_reconstruction_audit(report, audit, *, provenance_categories=None,
                                         likely_tandem_records=None):
    """Bind decisions to published or explicitly quarantined original rows.

    Call after report-extension validation when supplying quarantine context.
    Neither absent projected rows nor self-consistent provenance hashes alone
    establish an unrecorded historical ranking or legitimate demotion.
    """
    published = {}
    for array in report["arrays"]:
        for proof in _nested_proofs(array.get("classification_basis")):
            key = _digest(proof["selection"])
            if key in published:
                raise ValueError("Possible publication duplicates a selected proof")
            published[key] = proof
    if audit is None and not published:
        return
    if (not isinstance(audit, dict) or set(audit) != {
            "schema", "rule", "source_sha256", "coordinate_system", "policy", "confidence_category",
            "biological_accuracy_established", "stats", "actions", "decisions"}
            or audit["schema"] != {"name": "CRISPRidentify-possible-reconstruction-evidence", "version": "1.0.0"}
            or audit.get("rule") != RULE or audit["source_sha256"] != report["source"]["sha256"]
            or audit["coordinate_system"] != "zero_based_half_open"
            or audit["confidence_category"] != "Possible" or audit["biological_accuracy_established"] is not False
            or _encoded(audit["policy"]) != _encoded(POLICY.as_dict())):
        raise ValueError("Possible transaction source/policy differs")
    stats = audit["stats"]
    counts = {"anchors", "validation_requests", "selected", "retained", "abstained",
              "boundary_uncertain_selected", "added_repeats", "fresh_candidates_scored",
              "candidates_validated", "selections_committed"}
    if (not isinstance(stats, dict) or set(stats) != counts | {"reason_counts"}
            or any(type(stats[key]) is not int or stats[key] < 0 for key in counts)
            or not isinstance(audit["decisions"], list)):
        raise ValueError("Possible audit counters/decisions are malformed")
    actions, reasons = Counter(), Counter()
    for decision in audit["decisions"]:
        if (decision.get("rule") != RULE or decision.get("rule_version") != "1.0.0"
                or decision.get("action") not in ("selected", "retained", "abstained")
                or decision.get("boundary_uncertain") is not False
                or not isinstance(decision.get("reason"), str)
                or type(decision["anchor"]["category_index"]) is not int
                or decision["anchor"]["category_index"] != 2
                or _encoded(decision["evidence"]["policy"]) != _encoded(POLICY.as_dict())):
            raise ValueError("Possible audit decision scope differs")
        actions[decision["action"]] += 1
        reasons[decision["reason"]] += 1
    decisions = [d for d in audit["decisions"] if d["action"] == "selected"]
    if (_encoded(audit["actions"]) != _encoded(dict(actions))
            or _encoded(stats["reason_counts"]) != _encoded(dict(reasons))
            or stats["anchors"] != len(audit["decisions"])
            or any(stats[action] != actions[action] for action in ("selected", "retained", "abstained"))
            or stats["boundary_uncertain_selected"] != 0
            or stats["validation_requests"] != stats["candidates_validated"]
            or stats["added_repeats"] != sum(len(d["relationship"]["added_repeat_indexes"]) for d in decisions)):
        raise ValueError("Possible audit full decision accounting differs")
    expected = {_digest(d): d for d in decisions}
    if len(expected) != len(decisions) or set(published).difference(expected):
        raise ValueError("Possible transaction decisions differ from published proofs")
    missing = set(expected).difference(published)
    observed = {}
    if provenance_categories is not None:
        for index, category in enumerate(provenance_categories):
            for key, rows in category.items():
                for ordinal, row in enumerate(rows):
                    basis = row[3].get("classification_basis") if len(row) == 4 else None
                    for proof in _nested_proofs(basis):
                        observed[_digest(proof["selection"])] = (index, key, ordinal, row, proof)
    records = [record.as_dict() if hasattr(record, "as_dict") else record
               for record in (likely_tandem_records or ())]
    for key in missing:
        if key not in observed:
            raise ValueError("Possible selected output is missing without retained disposition")
        index, locus, ordinal, row, proof = observed[key]
        valid = [record for record in records if record["action"] == "quarantined"
                 and record["source_category_index"] == index == 2
                 and _encoded(record["locus_key"]) == _encoded(locus)
                 and record["row_ordinal"] == ordinal and record["score"] == row[0]
                 and [[r["source_interval"]["start"], r["source_interval"]["end"]]
                      for r in record["source_architecture"]["repeats"]]
                 == proof["reconstruction"]["repeat_intervals"]]
        if len(valid) != 1:
            raise ValueError("Possible projected selection has no exact quarantine disposition")
        published[key] = proof
    if (audit["stats"]["selections_committed"] != len(expected)
            or any(proof["source_sha256"] != audit["source_sha256"] for proof in published.values())):
        raise ValueError("Possible transaction proof count/source differs")
