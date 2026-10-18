"""Bounded source-sealed CRISPR candidates from existing Vmatch edges."""

from dataclasses import dataclass
import hashlib
import json

from components.components_detection import _within_one_edit
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
)
from components.source_intact_existence_rescue import (
    ELIGIBLE,
    SourceIntactEvaluationContext,
    evaluate_source_intact_existence_in_context,
)


VMATCH_CHAIN_SCHEMA = "crispridentify.vmatch_occurrence_chain@1.1.0"
MINIMUM_CHAIN_OCCURRENCES = 5
MINIMUM_SOURCE_SPACER = 18
MAXIMUM_SOURCE_SPACER = 80
# Vmatch may report the same physical repeat through two seed alignments with
# slightly different boundaries.  Treat them as aliases only when both source
# boundaries agree to three bases and at least 80% of the shorter interval is
# the identical physical source slice.  Complete-link grouping below prevents
# a chain of small shifts from merging distinct repeat loci transitively.
MAXIMUM_ALIAS_BOUNDARY_DELTA = 3
MINIMUM_ALIAS_OVERLAP_FRACTION = 0.80
DEFAULT_EDGE_CAP = 4096
DEFAULT_OCCURRENCE_CAP = 4096
DEFAULT_COMPONENT_CAP = 128
DEFAULT_DRAFT_CAP = 256


def _canonical_json(value):
    return json.dumps(
        value,
        ensure_ascii=True,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _digest(value):
    return hashlib.sha256(_canonical_json(value).encode("ascii")).hexdigest()


def _sequence_sha256(value):
    return hashlib.sha256(value.upper().encode("ascii")).hexdigest()


def _candidate_geometry(candidate):
    return {
        "repeat_starts": list(candidate.list_repeat_starts),
        "repeats": list(candidate.list_repeats),
        "spacers": list(candidate.list_spacers),
        "consensus": candidate.consensus,
    }


@dataclass(frozen=True)
class VmatchChainEdge:
    edge_id: str
    first_start: int
    first_end: int
    first_sha256: str
    second_start: int
    second_end: int
    second_sha256: str
    seed: str

    def as_dict(self):
        return {
            "edge_id": self.edge_id,
            "first_start": self.first_start,
            "first_end": self.first_end,
            "first_sha256": self.first_sha256,
            "second_start": self.second_start,
            "second_end": self.second_end,
            "second_sha256": self.second_sha256,
            "seed": self.seed,
        }


@dataclass(frozen=True)
class VmatchChainAttestation:
    registry_id: str
    source_sha256: str
    component_id: str
    chain_id: str
    cluster_index: int
    occurrence_ids: tuple
    occurrence_groups: tuple
    repeat_starts: tuple
    repeat_ends: tuple
    edge_records: tuple
    geometry_sha256: str
    existence_decision: dict

    def as_dict(self):
        return {
            "schema": VMATCH_CHAIN_SCHEMA,
            "registry_id": self.registry_id,
            "source_sha256": self.source_sha256,
            "component_id": self.component_id,
            "chain_id": self.chain_id,
            "cluster_index": self.cluster_index,
            "occurrence_ids": list(self.occurrence_ids),
            "occurrence_groups": [
                {
                    "occurrence_id": group[0],
                    "representative": list(group[1]),
                    "members": [list(member) for member in group[2]],
                }
                for group in self.occurrence_groups
            ],
            "repeat_starts": list(self.repeat_starts),
            "repeat_ends": list(self.repeat_ends),
            "edge_records": [
                edge.as_dict() for edge in self.edge_records
            ],
            "geometry_sha256": self.geometry_sha256,
            "existence_decision": self.existence_decision,
        }


class VmatchChainRegistry:
    """Private object seals plus immutable source/edge attestations."""

    def __init__(self, source_sequence, registered, telemetry):
        self._source_sha256 = _sequence_sha256(source_sequence)
        self._registered = tuple(registered)
        self._by_registry_id = {
            attestation.registry_id: (attestation, candidate)
            for attestation, candidate in self._registered
        }
        if len(self._by_registry_id) != len(self._registered):
            raise ValueError("Vmatch chain registry IDs are not unique")
        self.telemetry = dict(telemetry)

    @property
    def source_sha256(self):
        return self._source_sha256

    def entries(self):
        return tuple(self._registered)

    def verify(self, source_sequence, attestation, candidate):
        if _sequence_sha256(source_sequence) != self._source_sha256:
            return False
        registered = self._by_registry_id.get(
            getattr(attestation, "registry_id", None)
        )
        if (
            registered is None
            or registered[0] is not attestation
            or registered[1] is not candidate
        ):
            return False
        if _digest(_candidate_geometry(candidate)) != attestation.geometry_sha256:
            return False
        if (
            tuple(candidate.list_repeat_starts) != attestation.repeat_starts
            or tuple(
                start + len(repeat) for start, repeat in zip(
                    candidate.list_repeat_starts, candidate.list_repeats
                )
            ) != attestation.repeat_ends
        ):
            return False
        source = source_sequence.upper()
        if (
            tuple(group[0] for group in attestation.occurrence_groups)
            != attestation.occurrence_ids
            or tuple(group[1][1] for group in attestation.occurrence_groups)
            != attestation.repeat_starts
            or tuple(group[1][2] for group in attestation.occurrence_groups)
            != attestation.repeat_ends
        ):
            return False
        raw_to_group = {}
        for group_id, representative, members in attestation.occurrence_groups:
            if not members or representative not in members:
                return False
            for member in members:
                if len(member) != 4:
                    return False
                raw_id, start, end, sequence_hash = member
                if (
                    start < 0
                    or end <= start
                    or end > len(source)
                    or _sequence_sha256(source[start:end]) != sequence_hash
                    or raw_id != _digest({
                        "start": start,
                        "end": end,
                        "source_slice_sha256": sequence_hash,
                    })
                    or raw_id in raw_to_group
                ):
                    return False
                raw_to_group[raw_id] = group_id
            if any(
                not _same_source_locus(first, second)
                for member_index, first in enumerate(members)
                for second in members[member_index + 1:]
            ):
                return False
            expected_group_id = _digest({
                "representative": list(representative),
                "member_occurrence_ids": sorted(
                    member[0] for member in members
                ),
            })
            if group_id != expected_group_id:
                return False

        connected_group_pairs = set()
        for edge in attestation.edge_records:
            first = source[edge.first_start:edge.first_end]
            second = source[edge.second_start:edge.second_end]
            if (
                _sequence_sha256(first) != edge.first_sha256
                or _sequence_sha256(second) != edge.second_sha256
                or not _within_one_edit(first, second)
            ):
                return False
            first_raw = _occurrence(edge, 0)[0]
            second_raw = _occurrence(edge, 1)[0]
            first_group = raw_to_group.get(first_raw)
            second_group = raw_to_group.get(second_raw)
            if (
                first_group is None
                or second_group is None
                or first_group == second_group
            ):
                return False
            connected_group_pairs.add(frozenset((first_group, second_group)))
        if any(
            frozenset((left, right)) not in connected_group_pairs
            for left, right in zip(
                attestation.occurrence_ids,
                attestation.occurrence_ids[1:],
            )
        ):
            return False
        return True


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def add(self, value):
        self.parent.setdefault(value, value)

    def find(self, value):
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            following = self.parent[value]
            self.parent[value] = root
            value = following
        return root

    def union(self, first, second):
        self.add(first)
        self.add(second)
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root == second_root:
            return
        if first_root > second_root:
            first_root, second_root = second_root, first_root
        self.parent[second_root] = first_root


def _edge_from_repeat(source, repeat):
    first_start = repeat.begin_first
    second_start = repeat.begin_second
    first_end = first_start + repeat.length_first
    second_end = second_start + repeat.length_second
    if (
        first_start < 0
        or second_start < 0
        or first_end > len(source)
        or second_end > len(source)
    ):
        raise ValueError("Vmatch edge lies outside the source")
    first = source[first_start:first_end].upper()
    second = source[second_start:second_end].upper()
    if not first or not second or not _within_one_edit(first, second):
        raise ValueError("Vmatch edge does not replay within one edit")
    edge_payload = {
        "first_start": first_start,
        "first_end": first_end,
        "first_sha256": _sequence_sha256(first),
        "second_start": second_start,
        "second_end": second_end,
        "second_sha256": _sequence_sha256(second),
        "seed": repeat.sequence.upper(),
    }
    return VmatchChainEdge(
        edge_id=_digest(edge_payload),
        **edge_payload
    )


def _occurrence(edge, side):
    if side == 0:
        start, end, sequence_hash = (
            edge.first_start, edge.first_end, edge.first_sha256
        )
    else:
        start, end, sequence_hash = (
            edge.second_start, edge.second_end, edge.second_sha256
        )
    occurrence_id = _digest({
        "start": start,
        "end": end,
        "source_slice_sha256": sequence_hash,
    })
    return occurrence_id, start, end, sequence_hash


def _same_source_locus(first, second):
    """Prove that two raw Vmatch endpoints are boundary aliases."""

    first_start, first_end = first[1], first[2]
    second_start, second_end = second[1], second[2]
    if (
        abs(first_start - second_start) > MAXIMUM_ALIAS_BOUNDARY_DELTA
        or abs(first_end - second_end) > MAXIMUM_ALIAS_BOUNDARY_DELTA
    ):
        return False
    overlap = min(first_end, second_end) - max(first_start, second_start)
    shorter = min(first_end - first_start, second_end - second_start)
    return (
        overlap > 0
        and overlap / float(shorter) >= MINIMUM_ALIAS_OVERLAP_FRACTION
    )


def _coalesce_source_occurrences(occurrences):
    """Return deterministic complete-link physical-repeat alias groups."""

    groups = []
    for occurrence in sorted(occurrences.values(), key=lambda item: (
            item[1], item[2], item[0])):
        compatible = [
            group_index
            for group_index, group in enumerate(groups)
            if all(_same_source_locus(occurrence, member) for member in group)
        ]
        if len(compatible) > 1:
            raise ValueError(
                "a Vmatch occurrence aliases multiple physical-repeat groups"
            )
        if compatible:
            groups[compatible[0]].append(occurrence)
        else:
            groups.append([occurrence])

    records = []
    raw_to_group = {}
    for group in groups:
        members = tuple(sorted(group, key=lambda item: (
            item[1], item[2], item[0]
        )))
        representative = min(members, key=lambda item: (
            -(item[2] - item[1]), item[1], item[2], item[0]
        ))
        group_id = _digest({
            "representative": list(representative),
            "member_occurrence_ids": sorted(member[0] for member in members),
        })
        record = (
            group_id,
            representative[1],
            representative[2],
            representative[3],
            members,
            representative,
        )
        records.append(record)
        for member in members:
            raw_to_group[member[0]] = group_id
    return tuple(records), raw_to_group


def _maximal_source_runs(occurrences, connected_pairs):
    if not occurrences:
        return ()
    runs = []
    current = [occurrences[0]]
    for occurrence in occurrences[1:]:
        previous = current[-1]
        spacer = occurrence[1] - previous[2]
        if (
            occurrence[1] < previous[2]
            or spacer < MINIMUM_SOURCE_SPACER
            or spacer > MAXIMUM_SOURCE_SPACER
            or frozenset((previous[0], occurrence[0])) not in connected_pairs
        ):
            runs.append(tuple(current))
            current = [occurrence]
        else:
            current.append(occurrence)
    runs.append(tuple(current))
    return tuple(runs)


def build_vmatch_chain_registry(
    source_sequence,
    clusters,
    *,
    edge_cap=DEFAULT_EDGE_CAP,
    occurrence_cap=DEFAULT_OCCURRENCE_CAP,
    component_cap=DEFAULT_COMPONENT_CAP,
    draft_cap=DEFAULT_DRAFT_CAP,
    alias_joined_only=True,
):
    """Build only maximal, strict source-intact chains from existing edges."""

    caps = {
        "edge_cap": edge_cap,
        "occurrence_cap": occurrence_cap,
        "component_cap": component_cap,
        "draft_cap": draft_cap,
    }
    if any(
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 1
        for value in caps.values()
    ):
        raise ValueError("Vmatch chain caps must be positive integers")
    if not isinstance(source_sequence, str) or not source_sequence:
        raise ValueError("Vmatch chain source must be a nonempty string")
    if not isinstance(alias_joined_only, bool):
        raise TypeError("alias-joined-only policy must be boolean")

    source = source_sequence.upper()
    source_sha256 = _sequence_sha256(source)
    source_context = SourceIntactEvaluationContext(source)
    registered = []
    telemetry = {
        "schema": VMATCH_CHAIN_SCHEMA,
        "clusters_reviewed": 0,
        "edges_replayed": 0,
        "components_reviewed": 0,
        "maximal_runs_reviewed": 0,
        "short_runs_rejected": 0,
        "conflicting_occurrence_components_rejected": 0,
        "ambiguous_alias_clusters_rejected": 0,
        "raw_occurrences_coalesced": 0,
        "alias_joined_runs_rejected": 0,
        "source_disconnected_runs_split": 0,
        "source_intact_rejected": 0,
        "drafts_registered": 0,
        "cap_abstentions": 0,
        **caps,
    }

    for cluster_index, cluster in enumerate(clusters):
        telemetry["clusters_reviewed"] += 1
        raw_repeats = tuple(getattr(cluster, "list_cluster_repeats", ()))
        if len(raw_repeats) > edge_cap:
            telemetry["cap_abstentions"] += 1
            continue
        try:
            edges = tuple(sorted(
                {_edge_from_repeat(source, repeat) for repeat in raw_repeats},
                key=lambda edge: edge.edge_id,
            ))
        except (AttributeError, TypeError, ValueError):
            telemetry["cap_abstentions"] += 1
            continue
        telemetry["edges_replayed"] += len(edges)

        occurrences = {}
        edge_endpoints = {}
        for edge in edges:
            first = _occurrence(edge, 0)
            second = _occurrence(edge, 1)
            occurrences[first[0]] = first
            occurrences[second[0]] = second
            edge_endpoints[edge.edge_id] = (first[0], second[0])
        if len(occurrences) > occurrence_cap:
            telemetry["cap_abstentions"] += 1
            continue

        try:
            occurrence_groups, raw_to_group = (
                _coalesce_source_occurrences(occurrences)
            )
        except ValueError:
            telemetry["ambiguous_alias_clusters_rejected"] += 1
            continue
        telemetry["raw_occurrences_coalesced"] += (
            len(occurrences) - len(occurrence_groups)
        )
        groups_by_id = {group[0]: group for group in occurrence_groups}
        union = _UnionFind()
        edges_by_group = {}
        connected_pairs = set()
        edge_groups = {}
        for edge in edges:
            first_raw, second_raw = edge_endpoints[edge.edge_id]
            first_group = raw_to_group[first_raw]
            second_group = raw_to_group[second_raw]
            if first_group == second_group:
                continue
            union.union(first_group, second_group)
            connected_pairs.add(frozenset((first_group, second_group)))
            edges_by_group.setdefault(first_group, set()).add(edge)
            edges_by_group.setdefault(second_group, set()).add(edge)
            edge_groups[edge.edge_id] = (first_group, second_group)

        components = {}
        for occurrence_id in sorted(union.parent):
            components.setdefault(
                union.find(occurrence_id), []
            ).append(groups_by_id[occurrence_id])
        if len(components) > component_cap:
            telemetry["cap_abstentions"] += 1
            continue

        for component_occurrences in components.values():
            telemetry["components_reviewed"] += 1
            by_start = {}
            conflict = False
            for occurrence in component_occurrences:
                previous = by_start.get(occurrence[1])
                if previous is not None and previous[2:] != occurrence[2:]:
                    conflict = True
                    break
                by_start[occurrence[1]] = occurrence
            if conflict:
                telemetry[
                    "conflicting_occurrence_components_rejected"
                ] += 1
                continue
            ordered = tuple(sorted(by_start.values(), key=lambda item: (
                item[1], item[2], item[0]
            )))
            component_edge_records = tuple(sorted({
                edge
                for occurrence in ordered
                for edge in edges_by_group.get(occurrence[0], ())
            }, key=lambda edge: edge.edge_id))
            component_id = _digest({
                "source_sha256": source_sha256,
                "occurrence_ids": [item[0] for item in ordered],
                "edge_ids": [
                    edge.edge_id for edge in component_edge_records
                ],
            })

            source_runs = _maximal_source_runs(ordered, connected_pairs)
            if len(source_runs) > 1:
                telemetry["source_disconnected_runs_split"] += (
                    len(source_runs) - 1
                )
            for run in source_runs:
                telemetry["maximal_runs_reviewed"] += 1
                if len(run) < MINIMUM_CHAIN_OCCURRENCES:
                    telemetry["short_runs_rejected"] += 1
                    continue
                if (
                    alias_joined_only
                    and not any(len(occurrence[4]) > 1 for occurrence in run)
                ):
                    telemetry["alias_joined_runs_rejected"] += 1
                    continue
                if len(registered) >= draft_cap:
                    telemetry["cap_abstentions"] += 1
                    continue
                starts = tuple(item[1] for item in run)
                ends = tuple(item[2] for item in run)
                try:
                    candidate = build_source_backed_candidate_from_boundaries(
                        source, None, starts, ends
                    )
                except (TypeError, ValueError):
                    telemetry["source_intact_rejected"] += 1
                    continue
                existence = evaluate_source_intact_existence_in_context(
                    source_context, candidate
                )
                if existence.decision != ELIGIBLE:
                    telemetry["source_intact_rejected"] += 1
                    continue

                occurrence_ids = tuple(item[0] for item in run)
                run_ids = set(occurrence_ids)
                run_edge_records = tuple(
                    edge for edge in component_edge_records
                    if set(edge_groups.get(edge.edge_id, ())).issubset(run_ids)
                )
                if any(
                    frozenset((left, right)) not in {
                        frozenset(edge_groups[edge.edge_id])
                        for edge in run_edge_records
                    }
                    for left, right in zip(
                        occurrence_ids, occurrence_ids[1:]
                    )
                ):
                    telemetry["source_intact_rejected"] += 1
                    continue
                chain_id = _digest({
                    "component_id": component_id,
                    "occurrence_ids": list(occurrence_ids),
                })
                geometry_sha256 = _digest(_candidate_geometry(candidate))
                registry_id = _digest({
                    "schema": VMATCH_CHAIN_SCHEMA,
                    "source_sha256": source_sha256,
                    "chain_id": chain_id,
                    "geometry_sha256": geometry_sha256,
                })
                attestation = VmatchChainAttestation(
                    registry_id=registry_id,
                    source_sha256=source_sha256,
                    component_id=component_id,
                    chain_id=chain_id,
                    cluster_index=cluster_index,
                    occurrence_ids=occurrence_ids,
                    occurrence_groups=tuple(
                        (
                            item[0],
                            tuple(item[5]),
                            tuple(item[4]),
                        )
                        for item in run
                    ),
                    repeat_starts=starts,
                    repeat_ends=ends,
                    edge_records=run_edge_records,
                    geometry_sha256=geometry_sha256,
                    existence_decision=existence.as_dict(),
                )
                registered.append((attestation, candidate))

    registered.sort(key=lambda item: (
        item[0].repeat_starts[0],
        item[0].repeat_ends[-1],
        item[0].registry_id,
    ))
    telemetry["drafts_registered"] = len(registered)
    return VmatchChainRegistry(source, registered, telemetry)
