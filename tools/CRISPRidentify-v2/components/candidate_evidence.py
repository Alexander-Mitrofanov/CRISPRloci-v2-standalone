"""Opt-in, source-validated research evidence; never additional accepted arrays.

All five final internal categories are retained in a deterministic gzip JSONL
sidecar. Search scopes are comparison contexts, not inferred biological loci.
Only compact provenance is retained; archived source plus repeat intervals and
gapped repeats reconstruct each representation without duplicating whole DNA.
"""

from collections import Counter, defaultdict
import gzip
import hashlib
import json

from components.initial_competition import competition_evidence, validate_ordinary_witness_inventory
from components.candidate_origins import candidate_origin_evidence
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.reconstruction_validator import validate_candidate_reconstruction
from components.report_model import (
    INPUT_CATEGORY_NAMES, _feature_sets, _json_safe, _score_value,
)
from components.likely_tandem_transaction import validated_row_tags
from components.reconstruction_boundary_family import boundary_family_summary
from components.possible_reconstruction import possible_reconstruction_summary
from components.contextual_spacer_rescue import OBSERVATION_KEY, contextual_spacer_summary


EVIDENCE_FILENAME = "candidate_evidence.jsonl.gz"
SCHEMA = {"name": "CRISPRidentify-candidate-evidence", "version": "1.0.0"}


def _encoded(value):
    return json.dumps(value, allow_nan=False, sort_keys=True,
                      separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def _identity(value):
    return hashlib.sha256(_encoded(value)).hexdigest()


def _intervals(candidate):
    return [[start, start + len(repeat)] for start, repeat in zip(
        candidate.list_repeat_starts, candidate.list_repeats)]


def _scope_key(key):
    if hasattr(key, "start") and hasattr(key, "end"):
        return [key.start, key.end]
    return _json_safe(key)


def _features(row, names):
    sets = _feature_sets(row, names)
    # Keep the general feature-set contract as well as the compact model map.
    flat = ({item["name"]: item["value"] for item in sets[0]["features"]}
            if len(sets) == 1 else {})
    return sets, flat


def _lineage(candidate):
    records = getattr(candidate, "source_path_provenance", ()) or ()
    counts = Counter(str(getattr(record, "derivation", "unknown"))
                     for record in records)
    return {"record_count": len(records), "derivations": dict(sorted(counts.items()))}


def _describe_observation(row, *, category_index, category_position, key,
                          source_sequence, source_hash, feature_names,
                          canonical_ids):
    """Shared exact final-row identity; no score/feature recomputation."""
    legacy_key = _scope_key(key)
    scope_id = _identity([source_hash, legacy_key])
    candidate = row[1]
    validity = validate_candidate_reconstruction(source_sequence, candidate)
    if not validity.valid:
        raise ValueError("candidate evidence source reconstruction failed: {}"
                         .format(validity.as_dict()))
    gapped = candidate.list_gaped_repeats
    if (len(gapped) != len(candidate.list_repeats)
            or any(a.replace(" ", "").replace("-", "") != r
                   for a, r in zip(gapped, candidate.list_repeats))):
        raise ValueError("candidate evidence alignment/source mismatch")
    eligible = getattr(candidate, "source_path_classification_eligible", True)
    if not isinstance(eligible, bool):
        raise ValueError("candidate evidence eligibility must be boolean")
    score = _score_value(row[0])
    if score is None:
        raise ValueError("candidate evidence requires a genuine finite score")
    feature_sets, feature_values = _features(row, feature_names)
    geometry_id = _identity([source_hash, _intervals(candidate)])
    alignment_id = _identity(gapped)
    cid = _identity([geometry_id, alignment_id, score, feature_sets,
                     eligible, category_index, scope_id])
    metadata = row[3] if len(row) > 3 else {}
    if not isinstance(metadata, dict):
        raise ValueError("candidate evidence metadata must be a mapping")
    basis = metadata.get("classification_basis", {})
    if not isinstance(basis, dict):
        raise ValueError("candidate evidence classification basis must be a mapping")
    descriptor = {
        "candidate_id": cid,
        "observation_id": _identity([cid, category_index, category_position]),
        "source_geometry_id": geometry_id, "alignment_id": alignment_id,
        "comparison_scope_id": scope_id, "legacy_locus_key": legacy_key,
        "category": INPUT_CATEGORY_NAMES[category_index],
        "category_index": category_index, "category_position": category_position,
        "accepted": category_index in (0, 2), "score": score,
        "score_origin": "existing_final_row",
        "feature_repeat_count_consistent": (
            feature_values.get("number_repeats") == len(candidate.list_repeats)
            if "number_repeats" in feature_values else None),
        "classification_eligible": eligible,
        "eligibility_reason": getattr(candidate, "source_path_classification_reason", None),
        "feature_sets": feature_sets, "feature_values": feature_values,
        "canonical_array_ids": sorted(canonical_ids.get(
            (INPUT_CATEGORY_NAMES[category_index], geometry_id), [])),
    }
    return descriptor, basis


def _observation_record(descriptor, row, basis, *, structural, source_sequence,
                        source_hash, parameters):
    candidate = row[1]
    evidence = structural.evidence(candidate).as_dict()
    repeat_intervals = _intervals(candidate)
    tandem_review = basis.get("accepted_family_tandem_review")
    record = dict(descriptor,
        record_type="candidate", source_valid=True,
        source_span=[repeat_intervals[0][0], repeat_intervals[-1][1]],
        repeat_count=len(repeat_intervals), repeat_intervals=repeat_intervals,
        gapped_repeats=candidate.list_gaped_repeats,
        structural_filter=evidence,
        source_path_provenance=_lineage(candidate),
        classification_basis_method=basis.get("method"),
        negative_evidence={
            "tags": list(validated_row_tags(row)),
            "accepted_family_tandem_review": _json_safe(tandem_review),
            "accepted_family_tandem_review_status": (
                "recorded" if tandem_review is not None else "not_recorded"),
            "final_classification_method": basis.get("method"),
        })
    family = boundary_family_summary(row)
    if family is not None:
        record["boundary_family"] = family
    possible_summary = possible_reconstruction_summary(row)
    if possible_summary is not None:
        record["possible_reconstruction"] = possible_summary
    contextual_summary = contextual_spacer_summary(row)
    if contextual_summary is not None:
        record["contextual_spacer_context"] = contextual_summary
    if len(row) == 4 and OBSERVATION_KEY in row[3]:
        record[OBSERVATION_KEY] = _json_safe(row[3][OBSERVATION_KEY])
    competition = competition_evidence(row, len(source_sequence), source=source_sequence)
    if competition is not None:
        record["initial_competition"] = _json_safe(competition)
    if parameters.get("param_reconstruction_policy") == "conservative":
        record["candidate_origins"] = candidate_origin_evidence(candidate, source_hash)
    return record


def candidate_observation(row, *, category_index, category_position, key,
                          source_sequence, source_hash, feature_names,
                          structural, parameters, canonical_ids=None):
    """Return the same typed current observation as the optional sidecar.

    Reusers supply a source hash and structural filter computed once per report.
    Source slicing, alignment, origins and competition are still checked for
    every materialized observation. The result may share lists with the input;
    callers storing it must detach it before returning or mutating it.
    """
    descriptor, basis = _describe_observation(
        row, category_index=category_index, category_position=category_position,
        key=key, source_sequence=source_sequence, source_hash=source_hash,
        feature_names=feature_names, canonical_ids=canonical_ids or {})
    return _observation_record(
        descriptor, row, basis, structural=structural,
        source_sequence=source_sequence, source_hash=source_hash,
        parameters=parameters)


def stage_candidate_evidence(publication, path, *, source_header, source_sequence,
                             categories, parameters, feature_names,
                             canonical_report):
    """Validate and stage the complete final inventory in the report transaction."""
    if len(categories) != len(INPUT_CATEGORY_NAMES):
        raise ValueError("candidate evidence requires all five category slots")
    validate_ordinary_witness_inventory(categories, source_sequence)
    source_hash = hashlib.sha256(source_sequence.encode("utf-8")).hexdigest()
    source_id = source_header.split()[0].lstrip(">")
    structural = build_advanced_fuzzy_search_filter(parameters)
    canonical_ids = defaultdict(list)
    for array in canonical_report["arrays"]:
        if array["category"] not in ("Bona-fide", "Possible"):
            continue
        geometry = [[repeat["source_interval"]["start"] - 1,
                     repeat["source_interval"]["end"]]
                    for repeat in array["repeats"]]
        gid = _identity([source_hash, geometry])
        canonical_ids[(array["category"], gid)].append(array["id"])

    descriptors = []
    scopes = {}
    category_counts = Counter()
    for category_index, category in enumerate(categories):
        position = 0
        for key in sorted(category):
            legacy_key = _scope_key(key)
            scope_id = _identity([source_hash, legacy_key])
            scope = scopes.setdefault(scope_id, {
                "record_type": "comparison_scope", "comparison_scope_id": scope_id,
                "legacy_locus_key": legacy_key, "is_biological_locus": False,
                "accepted_candidate_ids": set(), "candidate_observations": 0,
            })
            for row in category[key]:
                descriptor, basis = _describe_observation(
                    row, category_index=category_index, category_position=position,
                    key=key, source_sequence=source_sequence,
                    source_hash=source_hash, feature_names=feature_names,
                    canonical_ids=canonical_ids)
                descriptors.append((descriptor, row, basis))
                scope["candidate_observations"] += 1
                if category_index in (0, 2):
                    scope["accepted_candidate_ids"].add(descriptor["candidate_id"])
                category_counts[INPUT_CATEGORY_NAMES[category_index]] += 1
                position += 1

    def write(binary_stream):
        digest = hashlib.sha256()
        with gzip.GzipFile(filename="", fileobj=binary_stream, mode="wb",
                           compresslevel=6, mtime=0) as output:
            def emit(record, include_digest=True):
                data = _encoded(record) + b"\n"
                output.write(data)
                if include_digest:
                    digest.update(data)
            emit({
                "record_type": "header", "schema": SCHEMA,
                "research_only": True, "not_additional_arrays": True,
                "capture_stage": "final_before_projection",
                "coordinate_system": "zero_based_half_open",
                "source": {"id": source_id, "sha256": source_hash,
                           "length": len(source_sequence)},
                "run": {"release_version": parameters.get("param_release_version"),
                        "preset": parameters.get("param_runtime_preset"),
                        "effective_config": parameters.get("param_effective_config", {}),
                        "retention_enabled": True},
                "comparison_scope_definition": "legacy search context; not a biological locus or evidence that candidates are the same array",
                "identity_semantics": "candidate_id identifies content; observation_id distinguishes repeated inventory rows; source_geometry_id ignores alignment variants",
                "negative_evidence_semantics": "not_recorded is unknown, never a passing negative screen",
                "provenance_retention": "source-path derivations summarized by counts; no recursive family history; final classification method and recorded accepted-family tandem review retained",
                "score_semantics": "existing final row score and features retained unchanged; no classifier recomputation or independent freshness certification",
                "annotations": "orientation, Cas and IS annotations are not inferred for unaccepted candidates",
            })
            for sid in sorted(scopes):
                scope = dict(scopes[sid])
                scope["accepted_candidate_ids"] = sorted(scope["accepted_candidate_ids"])
                emit(scope)
            for descriptor, row, basis in descriptors:
                record = _observation_record(
                    descriptor, row, basis, structural=structural,
                    source_sequence=source_sequence, source_hash=source_hash,
                    parameters=parameters)
                emit(record)
            emit({"record_type": "footer", "complete": True,
                  "candidate_observations": len(descriptors),
                  "comparison_scopes": len(scopes),
                  "category_counts": {name: category_counts[name] for name in INPUT_CATEGORY_NAMES},
                  "uncompressed_records_sha256": digest.hexdigest()}, include_digest=False)
    publication.stage_binary_writer(path, write)
