"""Local core authorization resolved by a complete declared source-phase family.

Witnesses authorize only the local retained-anchor map. The pending child's
own internal attestation remains the sole authority for its complete geometry.
"""

from copy import deepcopy
import hashlib

from components.components_detection import FuzzySourceOccurrence, FuzzySourcePathProvenance
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
    compare_source_backed_phase_candidates, SOURCE_CROSS_FAMILY_PHASE_RULE,
)
from components.internal_completion_authority import (
    CURRENT_INVENTORY_POLICY, authority_state, unrestricted, parent_current_allowed,
)
from components.internal_repeat_completion import _json, replay_internal_attestation
from components.repeat_boundary_repair import source_geometry
from components.terminal_repair_completion import _geometry, _overlap, _corridor, _replay_attestation


POLICY = "complete-family-local-core-composition-v1"
MAX_FAMILY_GEOMETRIES = 32


def _digest(value):
    return hashlib.sha256(_json(value).encode("ascii")).hexdigest()


def _relevant(geometry, target):
    return len(geometry[0]) == len(target[0]) and _overlap(geometry, target)


def _restore(source, item):
    candidate = build(source, None, *_geometry(item["geometry"]))
    for name, value in item["authority"].items():
        if name == "source_path_provenance" and value is not None:
            records = []
            for record in value:
                record = list(record)
                record[3] = tuple(FuzzySourceOccurrence(*node) for node in record[3])
                for index in (4, 5, 7):
                    record[index] = tuple(record[index])
                records.append(FuzzySourcePathProvenance(*record))
            value = tuple(records)
        setattr(candidate, name, value)
    return candidate


def collect_bridge_inventory(target, boundary_records, current, pending):
    """Capture the complete caller-owned family before deciding its winner."""
    witnesses = []
    for record in boundary_records:
        if record.get("reason") not in ("unique_minimum_supported_trim", "child_geometry_already_present"):
            continue
        if _relevant(_geometry(record["child_geometry"]), target):
            witnesses.append(deepcopy(record))
    pending_records = []
    for candidate in pending:
        if _relevant(source_geometry(candidate), target):
            pending_records.append({"geometry": source_geometry(candidate),
                                    "authority": authority_state(candidate),
                                    "attestation": deepcopy(getattr(candidate, "internal_repair_record", None))})
    # Include all current source-overlapping rows, regardless of count, so
    # original witness parents and the short core keep their restrictions.
    current_rows = [{"geometry": geometry, "authority": authority_state(candidate)}
                    for geometry, candidates in current.items() if _overlap(geometry, target)
                    for candidate in candidates]
    return {"complete": True, "boundary_records": sorted(witnesses, key=_json),
            "current_rows": sorted(current_rows, key=_json),
            "pending_records": sorted(pending_records, key=_json)}


def _phase_evidence(result):
    return {name: value for name, value in result.items()
            if name not in ("winner_candidate", "candidate_evidence")} | {
        "candidate_evidence": [{name: value for name, value in item.items() if name != "candidate"}
                               for item in result["candidate_evidence"]]}


def resolve_bridge_family(source, record, inventory, structural, context):
    """Replay every witness and every contributing current/pending derivation."""
    target = _geometry(record["child_geometry"])
    if inventory.get("complete") is not True:
        raise ValueError("incomplete local-core family inventory")
    family = {target: build(source, None, *target)}
    current = {}
    for item in inventory["current_rows"]:
        candidate = _restore(source, item)
        geometry = source_geometry(candidate)
        current.setdefault(geometry, []).append(candidate)
        if _relevant(geometry, target):
            if not unrestricted(source, candidate):
                raise ValueError("restricted current source-phase family member")
            family[geometry] = build(source, None, *geometry)
    if not parent_current_allowed(source, _geometry(record["parent_geometry"]), current,
                                  structural, context.source_sha256):
        raise ValueError("restricted internal parent in local-core family")
    witnesses = []
    for attestation in inventory["boundary_records"]:
        parent, child = _replay_attestation(source, attestation, structural, context.source_sha256)
        geometry = source_geometry(child)
        if not _relevant(geometry, target):
            raise ValueError("unrelated declared boundary witness")
        if any(not unrestricted(source, candidate) for g in (parent, geometry)
               for candidate in current.get(g, ())):
            raise ValueError("restricted boundary witness derivation")
        family[geometry] = child
        witnesses.append({"identity": _digest(attestation), "parent_geometry": parent,
                          "child_geometry": geometry, "attestation": attestation})
    for item in inventory["pending_records"]:
        geometry = _geometry(item["geometry"])
        if not _relevant(geometry, target):
            raise ValueError("unrelated declared pending witness")
        authority = item["authority"]
        if (authority.get("source_path_classification_eligible") is not False
                or authority.get("source_path_classification_reason") != "pending_internal_completion"
                or authority.get("source_path_restriction_reasons", ()) not in ((), [])
                or "source_path_provenance" in authority):
            raise ValueError("restricted pending source-phase family member")
        child = replay_internal_attestation(source, item["attestation"], structural, context)
        if source_geometry(child) != geometry:
            raise ValueError("stale pending source-phase geometry")
        family[geometry] = child
    if len(family) > MAX_FAMILY_GEOMETRIES:
        raise ValueError("local-core phase family exceeds 32 geometries")
    ordered = tuple(sorted(family))
    if len(ordered) < 2:
        raise ValueError("insufficient local-core source-phase family")
    context.stats["local_core_phase_comparisons"] += 1
    result = compare_source_backed_phase_candidates(source, [family[g] for g in ordered])
    if result["status"] != "selected" or source_geometry(result["winner_candidate"]) != target:
        raise ValueError("local-core source phase did not select the pending child: " + result["reason"])
    return {"geometry_family": ordered, "phase": _phase_evidence(result),
            "witnesses": witnesses, "current": current}


def compose_local_core(source, core, record, inventory, structural, context, cache=None):
    if source != context.source:
        raise ValueError("local-core context is bound to another source")
    target = _geometry(record["child_geometry"])
    identity = _digest((context.source_sha256, inventory, vars(structural), SOURCE_CROSS_FAMILY_PHASE_RULE, target))
    cache = {} if cache is None else cache
    if identity in cache:
        context.stats["local_core_phase_cache_hits"] += 1
        resolved = cache[identity]
    else:
        resolved = resolve_bridge_family(source, record, inventory, structural, context)
        cache[identity] = resolved
    if any(not unrestricted(source, candidate) for candidate in resolved["current"].get(core, ())):
        raise ValueError("restricted accepted core in local composition")
    parent = _geometry(record["parent_geometry"])
    closed = list(zip(*_geometry(record["gap_closed_geometry"])))
    final = list(zip(*target))
    retained = [final[closed.index(interval)] for interval in zip(*parent)]
    count = len(core[0])
    maps = []
    for witness in resolved["witnesses"]:
        original, corrected = witness["parent_geometry"], witness["child_geometry"]
        offset = _corridor(core, original, corrected)
        if offset is None:
            continue
        local = list(zip(corrected[0][offset:offset+count], corrected[1][offset:offset+count]))
        parent_offsets = [index for index in range(len(retained)-count+1) if retained[index:index+count] == local]
        if len(parent_offsets) != 1:
            continue
        maps.append({"witness_identity": witness["identity"], "witness_offset": offset,
                     "original_core_geometry": core,
                     "witness_parent_local": (original[0][offset:offset+count], original[1][offset:offset+count]),
                     "corrected_core_geometry": tuple(zip(*local)), "parent_offset": parent_offsets[0]})
    if not maps or len({_json((item["corrected_core_geometry"], item["parent_offset"])) for item in maps}) != 1:
        raise ValueError("no unique agreeing local retained-anchor map")
    proof = {"policy": POLICY, "source_sha256": context.source_sha256,
            "structural_rule": dict(vars(structural)), "phase_rule": dict(SOURCE_CROSS_FAMILY_PHASE_RULE),
            "local_only": True, "whole_geometry_equivalence_claimed": False,
            "inventory": deepcopy(inventory), "inventory_sha256": _digest(inventory),
            "geometry_family": resolved["geometry_family"], "phase": resolved["phase"],
            "witnesses": resolved["witnesses"], "maps": maps}
    cache = getattr(context, "local_core_replay_cache", {})
    cache[_digest((context.source_sha256, vars(structural), SOURCE_CROSS_FAMILY_PHASE_RULE, core, record, proof))] = True
    context.local_core_replay_cache = cache
    return proof


def replay_local_core(source, core, record, proof, structural, context):
    if source != context.source:
        raise ValueError("local-core replay context is bound to another source")
    if proof.get("inventory_sha256") != _digest(proof["inventory"]):
        raise ValueError("local-core declared inventory changed")
    identity = _digest((context.source_sha256, vars(structural), SOURCE_CROSS_FAMILY_PHASE_RULE, core, record, proof))
    if identity in getattr(context, "local_core_replay_cache", {}):
        context.stats["local_core_proof_cache_hits"] += 1
        return
    replay = compose_local_core(source, core, record, proof["inventory"], structural, context)
    if _json(replay) != _json(proof):
        raise ValueError("local-core full-family proof differs from replay")


def validate_declared_family_against_run(proof, target, run_evidence):
    """Bind a report proof to the independently retained complete run audits."""
    if not isinstance(run_evidence, dict):
        raise ValueError("local-core report lacks complete run-level family audits")
    def boundary_key(record):
        record = dict(record)
        record.pop("final_geometry_outcomes", None)
        return _json(record)
    boundary = [boundary_key(record) for record in run_evidence["boundary_records"]
                if record.get("reason") in ("unique_minimum_supported_trim", "child_geometry_already_present")
                and _relevant(_geometry(record["child_geometry"]), target)]
    if sorted(boundary) != sorted(boundary_key(record) for record in proof["inventory"]["boundary_records"]):
        raise ValueError("local-core proof omits or changes a declared boundary-family member")
    pending = [_json(record) for record in run_evidence["internal_records"]
               if record.get("eligibility_decision") == "pending_single_core_transaction"
               and _relevant(_geometry(record["child_geometry"]), target)]
    if sorted(pending) != sorted(_json(item["attestation"]) for item in proof["inventory"]["pending_records"]):
        raise ValueError("local-core proof omits or changes a declared internal-family member")
    current = run_evidence.get("current_authority_inventory")
    if (not isinstance(current, dict) or current.get("complete") is not True
            or current.get("policy") != CURRENT_INVENTORY_POLICY
            or current.get("source_sha256") != proof.get("source_sha256")
            or not isinstance(current.get("target_geometries"), (list, tuple))
            or not isinstance(current.get("rows"), (list, tuple))):
        raise ValueError("local-core report lacks complete current-authority inventory")
    targets = [_geometry(geometry) for geometry in current["target_geometries"]]
    if target not in targets:
        raise ValueError("local-core current-authority inventory has a stale target scope")
    rows = []
    for item in current["rows"]:
        geometry = _geometry(item["geometry"])
        if (set(item) != {"geometry", "authority"} or not isinstance(item["authority"], dict)
                or not any(_overlap(geometry, scope) for scope in targets)):
            raise ValueError("malformed local-core current-authority inventory row")
        if _overlap(geometry, target):
            rows.append(_json(item))
    if sorted(rows) != sorted(_json(item) for item in proof["inventory"]["current_rows"]):
        raise ValueError("local-core proof omits or changes a declared current-family member")
