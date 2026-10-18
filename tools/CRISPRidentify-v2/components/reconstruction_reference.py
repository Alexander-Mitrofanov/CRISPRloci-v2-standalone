"""Added-repeat identity evidence from an original core, never a donor majority.

This pure helper consumes candidate strings/coordinates. Source reconstruction
and the unique strict original-core mapping remain independent caller checks.
Publication restores the strings from source and replays the same evidence.
"""

from components.likely_tandem import global_edit_distance


ANCHOR_REFERENCE = "anchor-consensus-v1"
MAPPED_CORE_REFERENCE = "mapped-identical-original-core-v1"
REFERENCE_KEYS = frozenset((
    "rule", "coordinate_system", "original_consensus", "original_core_intervals",
    "mapped_core_indexes", "mapped_core_intervals", "boundary_offsets",
    "mapped_core_sequence", "added_repeat_checks",
))


def check_reference_rule(rule):
    if rule not in (ANCHOR_REFERENCE, MAPPED_CORE_REFERENCE):
        raise ValueError("unsupported added-repeat reference rule")
    return rule


def recorded_reference_rule(relation):
    """Missing discriminator means historical anchor-only replay, not default.

    Exact source-derived evidence equality is checked by the transaction after
    this dispatch. Unknown/malformed explicit recipes never fall back to legacy.
    """
    if not isinstance(relation, dict):
        raise ValueError("reconstruction relationship must be a mapping")
    if "added_repeat_reference" not in relation:
        return ANCHOR_REFERENCE
    evidence = relation["added_repeat_reference"]
    if (not isinstance(evidence, dict) or set(evidence) != REFERENCE_KEYS
            or evidence["rule"] != MAPPED_CORE_REFERENCE):
        raise ValueError("malformed or unsupported added-repeat reference evidence")
    return MAPPED_CORE_REFERENCE


def assess_added_repeat_reference(anchor_repeats, consensus, anchor_geometry,
                                 donor_repeats, donor_geometry, mapping, added,
                                 minimum_identity, *, reference_rule=MAPPED_CORE_REFERENCE):
    """Assess one reference for all added copies; return JSON-safe diagnostics.

    The original-reference success path retains the prior identity values and
    needs no extra proof field. A failed identity may use an alternate reference
    only when every actual original repeat equals its consensus and every mapped
    core copy has the same source sequence and uniform boundary offsets. No added
    copy participates in constructing the reference. The caller has already
    established a unique strict mapping; this function does not relax its limits.
    """
    check_reference_rule(reference_rule)
    result = {"qualified": False, "reason": "added_repeat_identity_below_limit",
              "reference_status": "original_identity_below_limit",
              "identities": [], "reference_evidence": None}
    if not isinstance(consensus, str) or not consensus or set(consensus) - set("ACGT"):
        result.update(reason="anchor_consensus_not_canonical", reference_status="invalid_original_consensus")
        return result
    original = []
    for index in added:
        repeat = donor_repeats[index]
        if not repeat or set(repeat) - set("ACGT"):
            result.update(reason="added_repeat_not_canonical", reference_status="noncanonical_added_repeat")
            return result
        distance = global_edit_distance(repeat, consensus)
        denominator = max(len(repeat), len(consensus))
        identity = 1.0 - distance / denominator
        original.append({"index": index, "original_edit_distance": distance,
                         "original_denominator": denominator, "original_identity": identity})
        # Preserve the original rule's early rejection behavior on replay.
        if identity < minimum_identity and reference_rule == ANCHOR_REFERENCE:
            return result
    if all(item["original_identity"] >= minimum_identity for item in original):
        result.update(qualified=True, reason=None, reference_status="original_reference_passed",
                      identities=[item["original_identity"] for item in original])
        return result
    result["original_checks"] = original
    if not anchor_repeats or any(repeat != consensus for repeat in anchor_repeats):
        result["reference_status"] = "original_core_not_unanimous"
        return result
    if (len(mapping) != len(anchor_geometry) or len(anchor_repeats) != len(anchor_geometry)
            or not mapping or len(set(mapping)) != len(mapping)):
        raise ValueError("mapped-core reference requires complete injective core indexes")
    mapped_repeats = [donor_repeats[index] for index in mapping]
    if any(not repeat or set(repeat) - set("ACGT") for repeat in mapped_repeats):
        result["reference_status"] = "mapped_core_not_canonical"
        return result
    offsets = [(donor_geometry[index][0] - interval[0],
                donor_geometry[index][1] - interval[1])
               for interval, index in zip(anchor_geometry, mapping)]
    if len(set(offsets)) != 1:
        result["reference_status"] = "mapped_core_offsets_not_uniform"
        return result
    if len(set(mapped_repeats)) != 1:
        result["reference_status"] = "mapped_core_sequence_not_unanimous"
        return result
    reference = mapped_repeats[0]
    checks = []
    for item in original:
        repeat = donor_repeats[item["index"]]
        distance = global_edit_distance(repeat, reference)
        denominator = max(len(repeat), len(reference))
        checks.append(dict(item, mapped_edit_distance=distance,
                           mapped_denominator=denominator,
                           mapped_identity=1.0 - distance / denominator))
    evidence = {"rule": MAPPED_CORE_REFERENCE, "coordinate_system": "zero_based_half_open",
                "original_consensus": consensus,
                "original_core_intervals": [list(interval) for interval in anchor_geometry],
                "mapped_core_indexes": list(mapping),
                "mapped_core_intervals": [list(donor_geometry[i]) for i in mapping],
                "boundary_offsets": list(offsets[0]), "mapped_core_sequence": reference,
                "added_repeat_checks": checks}
    result.update(reference_evidence=evidence,
                  identities=[item["mapped_identity"] for item in checks])
    if any(item["mapped_identity"] < minimum_identity for item in checks):
        result["reference_status"] = "mapped_reference_identity_below_limit"
        return result
    result.update(qualified=True, reason=None, reference_status="mapped_reference_passed")
    return result
