"""Candidate-specific evidence for conservative retained-array reconstruction.

This module never transfers a core's negative-screen clearance to a proposed
geometry.  It checks the source, the normal structural conjunction, and the
existing small-array and source-CDS tandem rules on the proposal itself.  A
policy-inapplicable test is distinct from an applicable test that was not run.
The caller owns fresh classifier evaluation and geometric conflict selection.
"""

from dataclasses import asdict
import hashlib
import math
from numbers import Real

from components.accepted_family_tandem_validator import (
    _equal_period_runs,
    validate_accepted_family_tandem,
)
from components.architecture_rescue_coding_context import SourceCDS
from components.likely_tandem import evaluate_likely_tandem
from components.likely_tandem_transaction import validated_row_tags
from components.reconstruction_selection import CandidateValidation
from components.reconstruction_validator import validate_candidate_reconstruction


RULE_NAME = "retained-reconstruction-candidate-validation"
RULE_VERSION = "1.0.0"
_SOURCE_DIGEST_CACHE = (None, None)
_NEGATIVE_METHODS = frozenset((
    "source_backed_accepted_family_tandem_demotion",
    "low_spacer_diversity_tandem_demotion",
))


def _source_digest(source):
    global _SOURCE_DIGEST_CACHE
    if source is not _SOURCE_DIGEST_CACHE[0]:
        _SOURCE_DIGEST_CACHE = (
            source, hashlib.sha256(source.upper().encode("ascii")).hexdigest(),
        )
    return _SOURCE_DIGEST_CACHE[1]


def needs_coding_context(candidate):
    """Whether the existing direct CDS-tandem rule can find a period run.

    This cheap predicate only scopes the existing rule. It is not source
    validation and must not be used as a general coding/tandem clearance.
    """
    starts = tuple(candidate.list_repeat_starts)
    if any(type(start) is not int or start < 0 for start in starts):
        raise ValueError("repeat starts must be nonnegative integer offsets")
    if any(right <= left for left, right in zip(starts, starts[1:])):
        raise ValueError("repeat starts must be strictly ordered")
    return bool(_equal_period_runs(starts))


def _row(row, label):
    if not isinstance(row, (list, tuple)) or len(row) not in (2, 3, 4):
        raise ValueError("{} must be a scored candidate row".format(label))
    score = row[0]
    if (isinstance(score, bool) or not isinstance(score, Real)
            or not math.isfinite(float(score)) or not 0 <= score <= 1):
        raise ValueError("{} score must be finite and in [0, 1]".format(label))
    if len(row) == 4 and not isinstance(row[3], dict):
        raise ValueError("{} metadata must be a mapping".format(label))
    return row[1]


def _known_negative_annotations(*rows):
    reasons = set()
    for row in rows:
        if "tandem" in validated_row_tags(row):
            reasons.add("tandem_tag")
        metadata = row[3] if len(row) == 4 else {}
        basis = metadata.get("classification_basis", {})
        if not isinstance(basis, dict):
            raise ValueError("classification_basis must be a mapping")
        if basis.get("method") in _NEGATIVE_METHODS:
            reasons.add(basis["method"])
        if basis.get("disposition") == "likely_tandem":
            reasons.add("likely_tandem_disposition")
        review = basis.get("accepted_family_tandem_review")
        if review is not None and not isinstance(review, dict):
            raise ValueError("accepted-family tandem review must be a mapping")
        if (isinstance(review, dict) and "supported" in review
                and type(review["supported"]) is not bool):
            raise ValueError("accepted-family tandem support must be boolean")
        if isinstance(review, dict) and review.get("supported") is True:
            reasons.add("recorded_source_cds_tandem_supported")
    return {
        "status": "failed" if reasons else "passed",
        "reason": "known_negative_annotation" if reasons else "none_recorded",
        "recorded_failures": sorted(reasons),
        "scope": "recorded_failures_only; not blanket negative clearance",
    }


def _coding_evidence(source, candidate, source_cds, required):
    if not needs_coding_context(candidate):
        return {
            "status": "inapplicable",
            "reason": "equal_in_frame_period_run_not_found",
            "source_cds_assessed": False,
            "considered_source_cds": [],
            "scope": "existing_source_cds_tandem_rule_on_proposed_geometry",
        }
    if not required or source_cds is None:
        return {
            "status": "not_run",
            "reason": ("source_cds_unavailable" if required
                       else "required_coding_screen_disabled"),
            "source_cds_assessed": False,
            "considered_source_cds": [],
            "scope": "existing_source_cds_tandem_rule_on_proposed_geometry",
        }
    if not isinstance(source_cds, tuple):
        raise TypeError("source_cds must be an assessed tuple of SourceCDS")
    start = candidate.list_repeat_starts[0]
    end = candidate.list_repeat_starts[-1] + len(candidate.list_repeats[-1])
    relevant = []
    for cds in source_cds:
        if not isinstance(cds, SourceCDS):
            raise TypeError("source_cds must contain SourceCDS objects")
        if cds.end > len(source):
            raise ValueError("source CDS lies outside the source")
        if cds.start < end and cds.end > start:
            relevant.append(cds)
    relevant = tuple(sorted(set(relevant), key=lambda cds: (
        cds.start, cds.end, cds.is_complement, cds.confidence,
    )))
    decision = validate_accepted_family_tandem(
        source, candidate, (candidate,), relevant,
    )
    return {
        "status": "failed" if decision.supported else "passed",
        "reason": decision.reason,
        "source_cds_assessed": True,
        "considered_source_cds": [asdict(cds) for cds in relevant],
        "decision": decision.as_dict(),
        "scope": "existing_source_cds_tandem_rule_on_proposed_geometry",
    }


def validate_reconstruction_candidate(
        source_dna, original_row, fresh_row, structural_filter, *,
        source_cds=None, coding_screen_required=True):
    """Return final-geometry validation, without mutating either candidate row.

    ``fresh_row`` must come from the caller's fresh classifier evaluation of
    the SAME candidate object, with its own feature vector. This function
    checks that contract's shape and identity; it does not rerun the model or
    independently certify score freshness. ``source_cds=None`` is unavailable;
    an empty tuple means a completed CDS prediction found no records. An
    applicable but disabled/unavailable CDS screen blocks selection.

    ``negatives_passed`` means the enumerated policy checks passed or were
    explicitly inapplicable. It is not a claim of biological correctness.
    Cas or IS annotation is not an applicable gate of this policy.
    """
    original = _row(original_row, "original_row")
    candidate = _row(fresh_row, "fresh_row")
    if candidate is not original:
        raise ValueError("fresh_row must contain the same candidate object")
    if len(fresh_row) < 3 or fresh_row[2] is None:
        raise ValueError("fresh_row must contain its own computed features")
    if type(coding_screen_required) is not bool:
        raise TypeError("coding_screen_required must be boolean")
    reconstruction = validate_candidate_reconstruction(source_dna, candidate)
    evidence = {
        "rule_name": RULE_NAME,
        "rule_version": RULE_VERSION,
        "source_reconstruction": reconstruction.as_dict(),
        "score_freshness": "caller_supplied_fresh_model_evaluation_same_candidate",
        "policy_checks_passed": False,
        "biological_correctness_certified": False,
        "limitations": [
            "negative_clearance_is_limited_to_enumerated_policy_checks",
            "cas_and_is_annotations_are_not_required_gates",
        ],
    }
    if not reconstruction.valid:
        evidence["decision"] = "source_reconstruction_failed"
        return CandidateValidation(
            row=fresh_row, source_valid=False, structural_valid=False,
            negatives_passed=False, score_fresh=True, evidence=evidence,
        )
    evidence["source_sha256"] = _source_digest(source_dna)
    evidence["source_span_sha256"] = hashlib.sha256(source_dna[
        reconstruction.array_start:reconstruction.array_end
    ].upper().encode("ascii")).hexdigest()
    # The standard independent source validator checks ungapped sequences.
    # Recheck the classifier alignment too, since structural evidence uses it.
    aligned = getattr(candidate, "list_gaped_repeats", ())
    alignment_valid = (
        len(aligned) == len(candidate.list_repeats)
        and all(isinstance(gapped, str) for gapped in aligned)
        and len({len(gapped) for gapped in aligned}) == 1
        and all(gapped.replace(" ", "").replace("-", "") == repeat
                for gapped, repeat in zip(aligned, candidate.list_repeats))
    )
    eligible = getattr(candidate, "source_path_classification_eligible", True)
    evidence["candidate_eligibility"] = {
        "status": "passed" if eligible is True else "failed",
        "reason": getattr(candidate, "source_path_classification_reason", None),
    }
    evidence["alignment_reconstruction"] = {
        "status": "passed" if alignment_valid else "failed",
    }
    structural = structural_filter.evidence(candidate)
    evidence["structural_filter"] = structural.as_dict()
    structural_valid = alignment_valid and eligible is True and structural.passes
    if not structural_valid:
        evidence["decision"] = "structure_or_eligibility_failed"
        return CandidateValidation(
            row=fresh_row, source_valid=True, structural_valid=False,
            negatives_passed=False, score_fresh=True, evidence=evidence,
        )
    canonical = all(set(sequence.upper()) <= set("ACGT")
                    for sequence in (tuple(candidate.list_repeats)
                                     + tuple(candidate.list_spacers)))
    evidence["canonical_sequence"] = {
        "status": "passed" if canonical else "not_run",
        "reason": ("canonical_repeat_and_spacer_bases" if canonical
                   else "noncanonical_bases_require_separate_review"),
    }
    evidence["known_negative_annotations"] = _known_negative_annotations(
        original_row, fresh_row,
    )
    tandem = evaluate_likely_tandem(candidate.list_spacers)
    if tandem.likely_tandem:
        tandem_status = "failed"
    elif tandem.reason_code in (
            "spacer_count_below_policy_minimum",
            "spacer_count_above_policy_maximum"):
        tandem_status = "inapplicable"
    elif tandem.complete_link_threshold_satisfied is False:
        tandem_status = "passed"
    else:
        tandem_status = "not_run"
    evidence["likely_tandem"] = {
        "status": tandem_status,
        "reason": tandem.reason_code,
        "decision": tandem.as_dict(),
    }
    evidence["coding_tandem"] = _coding_evidence(
        source_dna, candidate, source_cds, coding_screen_required,
    )
    checks = ("canonical_sequence", "known_negative_annotations",
              "likely_tandem", "coding_tandem")
    negatives_passed = all(evidence[name]["status"] in ("passed", "inapplicable")
                           for name in checks)
    evidence["policy_checks_passed"] = negatives_passed
    evidence["decision"] = "passed" if negatives_passed else "negative_policy_not_cleared"
    return CandidateValidation(
        row=fresh_row, source_valid=True, structural_valid=True,
        negatives_passed=negatives_passed, score_fresh=True, evidence=evidence,
    )


__all__ = ["needs_coding_context", "validate_reconstruction_candidate"]
