"""Complete, bounded source tilings of anomalous internal spacers.

The geometry analyzer has no classification authority. Runtime proposals must
also carry independently replayable original lineage and remain pending until
one source-wide accepted-core transaction. All coordinates are half-open.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from statistics import median

import regex
from components.candidate_origins import mark_repair_origin

from components.components_detection import (
    FuzzySearch, FuzzySourceOccurrence, FuzzySourcePathProvenance,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
    validate_candidate_source_reconstruction,
)
from components.fuzzy_boundary_lattice import FuzzyEditLimits
from components.repeat_boundary_repair import propose_terminal_repairs, source_geometry
from components.terminal_repair_completion import _geometry, _replay_attestation


POLICY = "source-internal-gap-dag-v1"
BROKEN_REASON = "broken_local_edge_lineage_parent"


@dataclass(frozen=True)
class InternalSearchPolicy:
    max_parents_per_component: int = 8
    max_gap_bases: int = 800
    max_parent_gap_bases: int = 8192
    max_parent_interval_tests: int = 65536
    max_source_interval_tests: int = 262144
    max_gap_nodes: int = 512
    max_gap_paths: int = 256
    max_gap_edge_checks: int = 131841
    max_parent_edge_checks: int = 1048576
    max_source_edge_checks: int = 4194304
    max_gap_dp_states: int = 32768
    max_parent_dp_states: int = 131072
    max_source_dp_states: int = 1048576
    max_gap_dp_transitions: int = 4194304
    max_parent_dp_transitions: int = 16777216
    max_source_dp_transitions: int = 67108864

    def __post_init__(self):
        if any(type(value) is not int or value <= 0 for value in asdict(self).values()):
            raise ValueError("internal search ceilings must be positive integers")


class InternalSearchContext:
    def __init__(self, source, policy=None):
        if not isinstance(source, str):
            raise ValueError("source DNA must be a string")
        self.source = source
        self.source_sha256 = hashlib.sha256(source.encode("ascii")).hexdigest()
        self.policy = policy or InternalSearchPolicy()
        self.interval_cache = {}
        self.gap_cache = {}
        self.stats = Counter()
        self.parent_work = Counter()


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _source_valid(source, candidate):
    validate_candidate_source_reconstruction(
        source, candidate.list_repeats, candidate.list_repeat_starts,
        candidate.list_spacers,
    )
    starts, ends = _geometry(source_geometry(candidate))
    if set(source[starts[0]:ends[-1]]) - set("ACGT"):
        raise ValueError("internal completion requires canonical source DNA")
    return starts, ends


def _restrictions(candidate):
    value = getattr(candidate, "source_path_restriction_reasons", ())
    if not isinstance(value, (tuple, list)) or any(type(x) is not str for x in value):
        raise ValueError("invalid restriction reason collection")
    reasons = set(value)
    reason = getattr(candidate, "source_path_classification_reason", None)
    if reason is not None:
        if type(reason) is not str:
            raise ValueError("invalid classification restriction reason")
        reasons.add(reason)
    return sorted(reasons)


def _validate_lineage(source, candidate, records):
    if not isinstance(records, tuple) or not records:
        raise ValueError("missing original source-path provenance")
    for record in records:
        if not isinstance(record, FuzzySourcePathProvenance):
            raise ValueError("malformed original source-path provenance")
        proxy = FuzzySearch.__new__(FuzzySearch)
        proxy.repeat_candidate = record.repeat_candidate
        proxy.list_absolute_start = candidate.list_repeat_starts
        proxy.list_repeats = candidate.list_repeats
        proxy.list_spacers = candidate.list_spacers
        proxy._validate_source_path_record(record)
        intervals = tuple((x.start, x.end) for x in record.parent_geometry)
        root = build(source, None, *zip(*intervals))
        if tuple(root.list_repeats) != tuple(x.sequence for x in record.parent_geometry):
            raise ValueError("lineage parent sequence differs from current source")
        if tuple(root.list_spacers) != tuple(record.parent_spacers):
            raise ValueError("lineage parent spacers differ from current source")
        if (record.derivation == "local_chain_family_parent"
                and (record.source_index_span != (0, len(record.parent_geometry))
                     or not any(len(x) > 140 for x in record.parent_spacers))):
            raise ValueError("broken-edge lineage has no original complete broken path")


def parent_authorization(source, candidate, structural, source_hash):
    """Return an auditable explicit authorization or raise; never infer it."""
    geometry = _source_valid(source, candidate)
    reasons = _restrictions(candidate)
    boundary = getattr(candidate, "boundary_repair_records", ())
    lineage = getattr(candidate, "source_path_provenance", None)
    marker = getattr(candidate, "source_path_classification_eligible", None)
    if boundary and lineage is None and marker is None and not reasons:
        if len(boundary) != 1:
            raise ValueError("ambiguous boundary authorization")
        _, rebuilt = _replay_attestation(source, boundary[0], structural, source_hash)
        if source_geometry(rebuilt) != geometry:
            raise ValueError("boundary authorization geometry differs")
        return {"kind": "terminal_repair", "boundary_attestation": deepcopy(boundary[0]),
                "eligibility": True, "restriction_reasons": []}
    if type(marker) is not bool:
        raise ValueError("missing or invalid eligibility marker")
    _validate_lineage(source, candidate, lineage)
    broken = any(r.derivation == "local_chain_family_parent" for r in lineage)
    if broken:
        if marker is not False or reasons != [BROKEN_REASON]:
            raise ValueError("broken parent has inconsistent or additional restrictions")
    elif marker is not True or reasons:
        raise ValueError("unrelated parent restriction")
    return {"kind": "fuzzy_source_path", "eligibility": marker,
            "restriction_reasons": reasons,
            "lineage": [dict(record._asdict()) for record in lineage]}


def _restore_authorization(source, geometry, authorization, structural, source_hash):
    parent = build(source, None, *_geometry(geometry))
    if authorization["kind"] == "terminal_repair":
        parent.boundary_repair_records = (authorization["boundary_attestation"],)
        mark_repair_origin(parent, authorization["boundary_attestation"])
    elif authorization["kind"] == "fuzzy_source_path":
        records = []
        for raw in authorization["lineage"]:
            raw = dict(raw)
            raw["parent_geometry"] = tuple(FuzzySourceOccurrence(*x) for x in raw["parent_geometry"])
            for name in ("parent_spacers", "occurrence_ids", "source_index_span"):
                raw[name] = tuple(raw[name])
            records.append(FuzzySourcePathProvenance(**raw))
        parent.source_path_provenance = tuple(records)
        parent.source_path_classification_eligible = authorization["eligibility"]
        parent.source_path_restriction_reasons = tuple(authorization["restriction_reasons"])
        if authorization["restriction_reasons"]:
            parent.source_path_classification_reason = authorization["restriction_reasons"][0]
    else:
        raise ValueError("unknown parent authorization")
    if _json(parent_authorization(source, parent, structural, source_hash)) != _json(authorization):
        raise ValueError("parent authorization replay differs")
    return parent


class _GraphBudgetExceeded(Exception):
    pass


def analyze_gap_graph(context, input_nodes, start, end, lower, upper):
    """Exact path counts and minima on the complete bounded interval DAG.

    This pure graph entry point also permits independent exhaustive tests.
    Every interval is one geometry, regardless of equivalent alignments.
    """
    policy = context.policy
    unique = {}
    for node in input_nodes:
        interval = tuple(node["interval"])
        if interval not in unique or node["edit_cost"] < unique[interval]["edit_cost"]:
            unique[interval] = node
    nodes = [unique[key] for key in sorted(unique)]
    work = Counter(edge_checks=0, dp_states=0, dp_transitions=0)
    result = {"complete": False, "nodes": nodes, "graph_work": work,
              "diagnostic_paths_complete": False, "paths": []}
    if len(nodes) > policy.max_gap_nodes:
        return dict(result, reason="gap_node_budget")

    def take(name):
        for scope, counts in (("gap", work), ("parent", context.parent_work),
                              ("source", context.stats)):
            if counts[name] >= getattr(policy, "max_" + scope + "_" + name):
                raise _GraphBudgetExceeded(scope + "_" + name + "_budget")
        work[name] += 1
        context.parent_work[name] += 1
        context.stats[name] += 1

    # Each anchor is a source boundary of its adjacent existing repeat. Only
    # internal interval nodes contribute repeat count or edit cost.
    intervals = [(start, start)] + [tuple(node["interval"]) for node in nodes] + [(end, end)]
    sink = len(intervals) - 1
    edges = [[] for _ in intervals]
    try:
        for left, (_, left_end) in enumerate(intervals):
            for right in range(left + 1, len(intervals)):
                take("edge_checks")
                if lower <= intervals[right][0] - left_end <= upper:
                    edges[left].append(right)
        result["edge_count"] = sum(map(len, edges))
        states = [{} for _ in intervals]
        take("dp_states")
        # total paths, minimum cost, minimum multiplicity, two witnesses.
        states[0][0] = [1, 0, 1, [()]]
        for left in range(sink):
            for right in edges[left]:
                for count, (total, cost, multiplicity, witnesses) in sorted(states[left].items()):
                    take("dp_transitions")
                    internal = right != sink
                    new_count = count + int(internal)
                    new_cost = cost + (nodes[right - 1]["edit_cost"] if internal else 0)
                    new_witnesses = [path + (right - 1,) if internal else path for path in witnesses]
                    old = states[right].get(new_count)
                    if old is None:
                        take("dp_states")
                        states[right][new_count] = [total, new_cost, multiplicity, new_witnesses]
                    else:
                        old[0] += total
                        if new_cost < old[1]:
                            old[1:] = [new_cost, multiplicity, new_witnesses]
                        elif new_cost == old[1]:
                            old[2] += multiplicity
                            old[3] = sorted(set(old[3] + new_witnesses))[:2]
    except _GraphBudgetExceeded as error:
        result["reason"] = str(error)
        return result
    final = {count: state for count, state in states[sink].items() if count > 0}
    result.update(complete=True, copy_counts=sorted(final),
                  path_count=str(sum(state[0] for state in final.values())),
                  count_results=[{"copies": count, "path_count": str(state[0]),
                                  "minimum_cost": state[1], "minimum_path_count": str(state[2])}
                                 for count, state in sorted(final.items())])
    if not final:
        result["reason"] = "no_complete_tiling"
        return result
    minimum = min(state[1] for state in final.values())
    minima = [state for state in final.values() if state[1] == minimum]
    paths = [{"nodes": path, "cost": state[1]} for state in minima for path in state[3]]
    result.update(minimum_cost=minimum, minimum_path_count=str(sum(state[2] for state in minima)),
                  paths=paths[:policy.max_gap_paths],
                  diagnostic_paths_complete=sum(state[0] for state in final.values()) == len(paths[:policy.max_gap_paths]))
    if len(final) != 1:
        result["reason"] = "copy_count_ambiguity"
    elif len(minima) != 1 or minima[0][2] != 1:
        result["reason"] = "minimum_cost_geometry_ambiguity"
    else:
        result.update(reason="unique_minimum_tiling",
                      selected_intervals=tuple(nodes[index]["interval"] for index in minima[0][3][0]))
    return result


def _gap_search(context, motif, start, end, lower, upper):
    key = motif, start, end, lower, upper
    if key in context.gap_cache:
        context.stats["gap_cache_hits"] += 1
        return deepcopy(context.gap_cache[key])
    source, policy = context.source, context.policy
    context.stats["gaps_searched"] += 1
    context.stats["searched_gap_bases"] += end - start
    result = {"gap": (start, end), "complete": False, "interval_tests": 0,
              "nodes": [], "paths": [], "reason": "incomplete"}
    matcher = regex.compile("(?:" + motif + ")" + FuzzyEditLimits().pattern_suffix(), regex.ENHANCEMATCH)
    anchors = tuple(motif[i:i + 6] for i in range(len(motif) - 5))
    for left in range(start + lower, end - lower - len(motif) + 4):
        for length in range(max(1, len(motif) - 3), len(motif) + 4):
            right = left + length
            if right > end - lower:
                continue
            interval_key = motif, left, right
            result["interval_tests"] += 1
            if interval_key not in context.interval_cache:
                if len(context.interval_cache) >= policy.max_source_interval_tests:
                    result["reason"] = "source_interval_budget"
                    return result
                sequence = source[left:right]
                match = matcher.fullmatch(sequence)
                node = None
                if match and sum(match.fuzzy_counts) * 10 <= max(len(motif), length):
                    exact_anchor = next((anchor for anchor in anchors if anchor in sequence), None)
                    if exact_anchor is not None:
                        node = {"interval": (left, right), "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                                "edits": tuple(match.fuzzy_counts), "edit_cost": sum(match.fuzzy_counts),
                                "exact_anchor": exact_anchor}
                context.interval_cache[interval_key] = node
                context.stats["distinct_interval_tests"] += 1
                if node is not None:
                    context.stats["matching_intervals"] += 1
            node = context.interval_cache[interval_key]
            if node is not None:
                if len(result["nodes"]) == policy.max_gap_nodes:
                    result["reason"] = "gap_node_budget"
                    return result
                result["nodes"].append(node)
    nodes = result["nodes"]
    result.update(analyze_gap_graph(context, nodes, start, end, lower, upper))
    context.stats[result["reason"]] += 1
    if result["complete"]:
        context.gap_cache[key] = deepcopy(result)
    return result


def analyze_internal_geometry(context, parent, structural):
    """Pure geometry replay. Its result never authorizes classification."""
    source, policy = context.source, context.policy
    context.parent_work = Counter()
    geometry = _source_valid(source, parent)
    record = {"policy": POLICY, "search_policy": asdict(policy), "source_sha256": context.source_sha256,
              "coordinate_system": "zero_based_half_open", "parent_geometry": geometry,
              "structural_rule": dict(vars(structural)), "status": "abstained",
              "eligibility_decision": "not_authorized_by_geometry", "gaps": []}
    if len(geometry[0]) < 5:
        record["reason"] = "insufficient_original_repeats"
        return None, record
    counts = Counter(parent.list_repeats)
    maximum = max(counts.values())
    motifs = sorted(motif for motif, count in counts.items() if count == maximum)
    if len(motifs) != 1 or maximum * 5 < len(parent.list_repeats) * 4:
        record["reason"] = "no_unique_dominant_source_motif"
        return None, record
    motif = motifs[0]
    ordinary = [len(spacer) for spacer in parent.list_spacers
                if structural.min_avg_spacer_length <= len(spacer) <= structural.max_avg_spacer_length]
    if len(ordinary) < 3:
        record["reason"] = "insufficient_ordinary_spacers"
        return None, record
    center = median(ordinary)
    mad = median(abs(value - center) for value in ordinary)
    radius = max(12, 3 * mad)
    lower = max(structural.min_avg_spacer_length, math.floor(center - radius))
    upper = min(structural.max_avg_spacer_length, math.ceil(center + radius))
    record.update(motif=motif, motif_count=maximum, ordinary_spacer_lengths=ordinary,
                  spacer_median=center, spacer_mad=mad, local_spacer_bounds=(lower, upper))
    anomalous = [(index, end, start) for index, (end, start) in enumerate(zip(geometry[1], geometry[0][1:]))
                 if start - end > structural.max_avg_spacer_length]
    if not anomalous:
        record["reason"] = "no_anomalous_gap"
        return None, record
    record["searched_gap_bases"] = sum(end - start for _, start, end in anomalous)
    if any(end - start > policy.max_gap_bases for _, start, end in anomalous):
        record["reason"] = "gap_base_budget"
        return None, record
    if record["searched_gap_bases"] > policy.max_parent_gap_bases:
        record["reason"] = "parent_gap_base_budget"
        return None, record
    # The full declared domain is known without running the matcher; reject
    # a parent whose interval ceiling cannot cover its complete search.
    domains = [sum(max(0, end - lower - length - (start + lower) + 1)
                   for length in range(max(1, len(motif) - 3), len(motif) + 4))
               for _, start, end in anomalous]
    if sum(domains) > policy.max_parent_interval_tests:
        record["reason"] = "parent_interval_budget"
        return None, record
    inserted = []
    for (index, start, end), domain in zip(anomalous, domains):
        gap = _gap_search(context, motif, start, end, lower, upper)
        gap["original_gap_index"] = index
        record["gaps"].append(gap)
        if gap["reason"] == "unique_minimum_tiling":
            inserted.extend(gap["selected_intervals"])
    record["interval_tests"] = sum(gap["interval_tests"] for gap in record["gaps"])
    record["search_complete"] = all(gap["complete"] for gap in record["gaps"])
    failed = [gap["reason"] for gap in record["gaps"] if gap["reason"] != "unique_minimum_tiling"]
    if failed:
        record["reason"] = failed[0]
        return None, record
    intervals = sorted(list(zip(*geometry)) + inserted)
    child = build(source, None, *zip(*intervals))
    record["inserted_intervals"] = tuple(inserted)
    record["gap_closed_geometry"] = source_geometry(child)
    # Invoke stage 1 only on a newly built, source-attested intermediate.
    # Original eligibility and geometry remain untouched.
    checks = structural.evidence(child)
    if not checks.passes:
        repairs, boundary_records, _ = propose_terminal_repairs(source, [child], structural,
                                                               source_sha256=context.source_sha256)
        if len(repairs) == 1:
            child = repairs[0]
            record["post_gap_boundary_attestation"] = boundary_records[0]
            checks = structural.evidence(child)
    record["resulting_checks"] = checks.as_dict()
    if not checks.passes:
        record["reason"] = "final_structural_filter"
        return None, record
    record.update(status="geometry_supported", reason="unique_complete_internal_path",
                  child_geometry=source_geometry(child))
    return child, record


def propose_internal_repairs(context, candidates, structural):
    by_geometry = {}
    for candidate in candidates:
        by_geometry.setdefault(source_geometry(candidate), []).append(candidate)
    children, records = [], []
    attempted = 0
    for geometry, members in sorted(by_geometry.items(), key=lambda item: (
            -(item[0][1][-1] - item[0][0][0]), -len(item[0][0]), item[0])):
        if len(geometry[0]) < 5 or not any(len(x) > structural.max_avg_spacer_length for x in members[0].list_spacers):
            continue
        context.stats["candidate_geometries"] += 1
        try:
            authorizations = [parent_authorization(context.source, member, structural, context.source_sha256)
                              for member in members]
        except (ValueError, KeyError, TypeError, AttributeError) as error:
            records.append({"parent_geometry": geometry, "status": "abstained",
                            "reason": "parent_authorization_rejected", "detail": str(error)})
            context.stats["parent_authorization_rejected"] += 1
            continue
        motif_counts = Counter(members[0].list_repeats)
        maximum = max(motif_counts.values())
        qualifying = (sum(count == maximum for count in motif_counts.values()) == 1
                      and maximum * 5 >= len(members[0].list_repeats) * 4
                      and sum(structural.min_avg_spacer_length <= len(spacer) <= structural.max_avg_spacer_length
                              for spacer in members[0].list_spacers) >= 3)
        if not qualifying:
            _, record = analyze_internal_geometry(context, members[0], structural)
            records.append(record)
            continue
        if attempted >= context.policy.max_parents_per_component:
            records.append({"parent_geometry": geometry, "status": "abstained", "reason": "component_parent_budget"})
            context.stats["component_parent_budget"] += 1
            continue
        attempted += 1
        context.stats["parents_attempted"] += 1
        child, record = analyze_internal_geometry(context, members[0], structural)
        record["parent_authorizations"] = sorted({_json(item): item for item in authorizations}.values(), key=_json)
        records.append(record)
        context.stats[record["reason"]] += 1
        if child is not None:
            record["eligibility_decision"] = "pending_single_core_transaction"
            child.internal_repair_record = record
            child.source_path_classification_eligible = False
            child.source_path_classification_reason = "pending_internal_completion"
            children.append(child)
    return children, records


def replay_internal_attestation(source, record, structural, context=None):
    context = context or InternalSearchContext(source)
    if (record.get("policy") != POLICY or record.get("source_sha256") != context.source_sha256
            or record.get("search_policy") != asdict(context.policy)
            or record.get("structural_rule") != vars(structural)
            or record.get("eligibility_decision") != "pending_single_core_transaction"):
        raise ValueError("internal completion attestation policy/source differs")
    authorizations = record.get("parent_authorizations")
    if not isinstance(authorizations, list) or not authorizations:
        raise ValueError("missing internal parent authorizations")
    parents = [_restore_authorization(source, record["parent_geometry"], auth, structural, context.source_sha256)
               for auth in authorizations]
    child, replay = analyze_internal_geometry(context, parents[0], structural)
    if child is None:
        raise ValueError("internal geometry is not replayable")
    expected = dict(record)
    expected.pop("parent_authorizations", None)
    expected["eligibility_decision"] = "not_authorized_by_geometry"
    if _json(replay) != _json(expected):
        raise ValueError("internal completion attestation replay differs")
    return child
