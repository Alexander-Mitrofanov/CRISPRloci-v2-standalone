"""Bounded, source-replayed alternatives to fuzzy-match boundaries.

The legacy fuzzy detector deliberately follows a leftmost, non-overlapping
``finditer`` path.  That contract can select a deletion-shortened interval
even when a nearby same-length interval is a better biological occurrence.
An exact-only overlay cannot recover the latter when it contains a genuine
substitution.

This module searches a much narrower question: for each already detected
deletion-shortened occurrence, which intervals within an explicit boundary
radius restore its length toward (but never beyond) the candidate repeat
length while satisfying the *unchanged* fuzzy edit limits when their
boundaries are forced?  Occurrences that are not deletion-shortened remain
fixed at their baseline source interval.  The module then enumerates every
source-ordered, non-overlapping assignment containing one interval per
baseline occurrence.

Completeness is local to those declared interval domains.  Seed, interval,
node, and path ceilings are explicit.  If any ceiling is reached, the plan
abstains and emits no alternatives; a sampled prefix can never be mistaken
for a complete candidate universe.  Every complete lattice is bound to the original source and matcher version.
"""

import hashlib
import json
from dataclasses import asdict, dataclass, field, fields
from functools import lru_cache

import regex

from components.fuzzy_alignment import normalize_fuzzy_changes
from components.fuzzy_alignment import snapshot_owned_fuzzy_match_state
from components.fuzzy_alignment_merge import build_pairwise_fuzzy_alignment


FUZZY_BOUNDARY_LATTICE_SCHEMA = (
    "crispridentify.fuzzy-boundary-alternative-lattice.v2"
)
MATCHER_CONTRACT = (
    "regex-enhanced-fullmatch-deletion-restoration-boundary.v2"
)


def _sha256_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


@lru_cache(maxsize=None)
def _dataclass_field_names(record_type):
    """Return the serialized dataclass fields once per record type."""

    return tuple(item.name for item in fields(record_type))


def _dataclass_json_default(value):
    """Expose provenance dataclasses to the native JSON encoder."""

    value_type = type(value)
    if getattr(value_type, "__dataclass_fields__", None) is None:
        raise TypeError(
            "Object of type {} is not JSON serializable".format(
                value_type.__name__
            )
        )
    return {
        name: getattr(value, name)
        for name in _dataclass_field_names(value_type)
    }


def _canonical_sha256(value):
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        default=_dataclass_json_default,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


_PLAIN_ATOMIC_TYPES = frozenset((type(None), bool, int, float, str))


def _plain_value(value):
    """Recursively convert frozen provenance records to canonical JSON data."""

    value_type = type(value)
    if value_type in _PLAIN_ATOMIC_TYPES:
        return value
    if value_type is dict:
        return {
            key: _plain_value(item) for key, item in value.items()
        }
    if value_type is tuple or value_type is list:
        return [_plain_value(item) for item in value]
    record_type = value if isinstance(value, type) else value_type
    if getattr(record_type, "__dataclass_fields__", None) is not None:
        if record_type is value:
            raise TypeError(
                "asdict() should be called on dataclass instances"
            )
        return {
            name: _plain_value(getattr(value, name))
            for name in _dataclass_field_names(record_type)
        }
    if isinstance(value, dict):
        return {
            key: _plain_value(item) for key, item in value.items()
        }
    if isinstance(value, (tuple, list)):
        return [_plain_value(item) for item in value]
    return value


def _positive_integer(value, label):
    if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or value <= 0):
        raise ValueError("{} must be a positive integer".format(label))


def _non_negative_integer(value, label):
    if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or value < 0):
        raise ValueError(
            "{} must be a non-negative integer".format(label)
        )


@dataclass(frozen=True)
class FuzzyEditLimits:
    """The exact edit policy used to test every forced interval."""

    max_substitutions: int = 3
    max_insertions: int = 3
    max_deletions: int = 3
    max_total_errors: int = 6

    def __post_init__(self):
        for label, value in (
                ("substitution limit", self.max_substitutions),
                ("insertion limit", self.max_insertions),
                ("deletion limit", self.max_deletions),
                ("total error limit", self.max_total_errors)):
            _non_negative_integer(value, label)

    def accepts(self, fuzzy_counts):
        if len(fuzzy_counts) != 3:
            return False
        substitutions, insertions, deletions = fuzzy_counts
        if any(
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 0
                for value in fuzzy_counts):
            return False
        return (
            substitutions <= self.max_substitutions
            and insertions <= self.max_insertions
            and deletions <= self.max_deletions
            and sum(fuzzy_counts) <= self.max_total_errors
        )

    def pattern_suffix(self):
        return (
            "{{s<={s},i<={i},d<={d},i+d+s<={total}}}".format(
                s=self.max_substitutions,
                i=self.max_insertions,
                d=self.max_deletions,
                total=self.max_total_errors,
            )
        )


@dataclass(frozen=True)
class FuzzyBoundarySearchPolicy:
    """Finite local domain and materialization ceilings.

    The one-base default is intentionally conservative.  A caller may widen
    it only up to the unchanged total-error ceiling, while the interval cap
    remains an independent bound on actual matcher calls.
    """

    boundary_radius: int = 1
    max_seed_occurrences: int = 512
    max_interval_evaluations: int = 4096
    max_materialized_nodes: int = 512
    max_materialized_paths: int = 512

    def __post_init__(self):
        _non_negative_integer(self.boundary_radius, "boundary radius")
        for label, value in (
                ("seed ceiling", self.max_seed_occurrences),
                ("interval-evaluation ceiling",
                 self.max_interval_evaluations),
                ("node ceiling", self.max_materialized_nodes),
                ("path ceiling", self.max_materialized_paths)):
            _positive_integer(value, label)


@dataclass(frozen=True, repr=False)
class BoundFuzzyBoundarySource:
    """One immutable source string whose digest is computed only once."""

    sequence: str = field(repr=False)
    sha256: str = field(init=False)
    length: int = field(init=False)

    def __post_init__(self):
        if not isinstance(self.sequence, str):
            raise ValueError("fuzzy boundary source must be a string")
        object.__setattr__(self, "sha256", _sha256_text(self.sequence))
        object.__setattr__(self, "length", len(self.sequence))

    def __repr__(self):
        return (
            "BoundFuzzyBoundarySource(length={}, sha256={!r})".format(
                self.length, self.sha256
            )
        )


def bind_fuzzy_boundary_source(source_dna):
    """Hash a source once for reuse by multiple local boundary searches."""

    if isinstance(source_dna, BoundFuzzyBoundarySource):
        return source_dna
    return BoundFuzzyBoundarySource(source_dna)


@dataclass(frozen=True)
class FuzzyBoundarySeed:
    seed_id: str
    seed_index: int
    start: int
    end: int
    sequence: str
    slice_sha256: str
    fuzzy_counts: tuple


@dataclass(frozen=True)
class FuzzyBoundaryDomain:
    seed_id: str
    seed_index: int
    domain_kind: str
    start_min: int
    start_max: int
    end_min: int
    end_max: int
    minimum_length: int
    maximum_length: int
    interval_count: int
    intervals_sha256: str


@dataclass(frozen=True)
class FuzzyBoundaryNode:
    node_id: str
    start: int
    end: int
    sequence: str
    slice_sha256: str
    fuzzy_counts: tuple
    normalized_changes: tuple
    alignment_sha256: str
    seed_indexes: tuple


@dataclass(frozen=True)
class FuzzyBoundaryPath:
    path_id: str
    kind: str
    node_ids: tuple
    geometry_sha256: str


@dataclass(frozen=True)
class FuzzyBoundaryLattice:
    schema: str
    lattice_id: str
    status: str
    status_reasons: tuple
    alternatives_complete: bool
    interval_evaluation_complete: bool
    node_universe_complete: bool
    path_universe_complete: bool
    source_length: int
    source_sha256: str
    window_start: int
    window_end: int
    window_sha256: str
    repeat_candidate: str
    repeat_sha256: str
    matcher_contract: str
    matcher_version: str
    edit_limits: FuzzyEditLimits
    policy: FuzzyBoundarySearchPolicy
    seeds: tuple
    domains: tuple
    planned_interval_count: int
    evaluated_interval_count: int
    accepted_interval_count: int
    evaluation_sha256: str
    nodes: tuple
    compatibility_edges: tuple
    paths: tuple
    baseline_path_id: str

    @property
    def abstained(self):
        return not self.alternatives_complete

    def as_dict(self):
        """Return a JSON-safe canonical provenance object."""

        value = asdict(self)
        value["abstained"] = self.abstained
        return value


class FuzzyBoundaryMatch:
    """A ``regex.Match``-compatible view in the original search window."""

    __slots__ = (
        "_text",
        "_start",
        "_end",
        "_group",
        "_deletion_position_shift",
        "fuzzy_counts",
        "fuzzy_changes",
    )

    def __init__(
            self, text, start, end, fuzzy_counts, fuzzy_changes,
            deletion_position_shift=None):
        self._text = text
        self._start = start
        self._end = end
        self._group = text[start:end]
        self.fuzzy_counts = tuple(fuzzy_counts)
        self.fuzzy_changes = tuple(
            list(positions) for positions in fuzzy_changes
        )
        self._deletion_position_shift = deletion_position_shift

    def _fuzzy_alignment_cache_snapshot(self):
        return snapshot_owned_fuzzy_match_state(
            self.fuzzy_counts, self.fuzzy_changes
        )

    @property
    def deletion_position_shift(self):
        return self._deletion_position_shift

    def start(self):
        return self._start

    def end(self):
        return self._end

    def span(self):
        return self._start, self._end

    def group(self):
        return self._group


@dataclass(frozen=True)
class FuzzyBoundaryMatchPath:
    path_id: str
    node_ids: tuple
    matches: tuple


@dataclass(frozen=True)
class FuzzyBoundaryAlternativePlan:
    lattice: FuzzyBoundaryLattice
    alternatives: tuple


def _seed_id(source_sha256, repeat_sha256, index, start, end, counts):
    return "fuzzy-boundary-seed:" + _canonical_sha256({
        "source_sha256": source_sha256,
        "repeat_sha256": repeat_sha256,
        "seed_index": index,
        "start": start,
        "end": end,
        "fuzzy_counts": counts,
    })


def _source_seeds(
        repeat, text, baseline_matches, sequence_start, source_binding,
        edit_limits):
    seeds = []
    previous_end = sequence_start
    for index, match in enumerate(tuple(baseline_matches)):
        start = match.start()
        end = match.end()
        if (
                not isinstance(start, int)
                or isinstance(start, bool)
                or not isinstance(end, int)
                or isinstance(end, bool)
                or start < 0
                or end <= start
                or end > len(text)
                or start < previous_end - sequence_start):
            raise ValueError(
                "baseline fuzzy occurrences must be ordered, positive, "
                "and non-overlapping"
            )
        sequence = text[start:end]
        if match.group() != sequence:
            raise ValueError(
                "baseline fuzzy occurrence does not replay its search window"
            )
        try:
            counts = tuple(match.fuzzy_counts)
        except (AttributeError, TypeError):
            raise ValueError(
                "baseline fuzzy occurrence has no valid edit counts"
            )
        if not edit_limits.accepts(counts):
            raise ValueError(
                "baseline fuzzy occurrence exceeds the declared edit limits"
            )
        expected_length = (
            len(repeat) + counts[1] - counts[2]
        )
        if len(sequence) != expected_length:
            raise ValueError(
                "baseline fuzzy occurrence length contradicts its edit counts"
            )
        absolute_start = sequence_start + start
        absolute_end = sequence_start + end
        if (
                source_binding.sequence[absolute_start:absolute_end]
                != sequence):
            raise ValueError(
                "baseline fuzzy occurrence does not replay original source"
            )
        seed = FuzzyBoundarySeed(
            _seed_id(
                source_binding.sha256,
                _sha256_text(repeat),
                index,
                absolute_start,
                absolute_end,
                counts,
            ),
            index,
            absolute_start,
            absolute_end,
            sequence,
            _sha256_text(sequence),
            counts,
        )
        seeds.append(seed)
        previous_end = absolute_end
    if not seeds:
        raise ValueError(
            "fuzzy boundary search requires a non-empty baseline path"
        )
    return tuple(seeds)


def _validate_replay_seeds(lattice, source_binding):
    """Re-slice every imported seed instead of trusting recorded indexes."""

    if _sha256_text(lattice.repeat_candidate) != lattice.repeat_sha256:
        raise ValueError("fuzzy boundary lattice repeat identity is stale")
    seeds = lattice.seeds
    if not isinstance(seeds, tuple) or not seeds:
        raise ValueError("fuzzy boundary lattice seeds are malformed")
    previous_end = lattice.window_start
    observed_ids = set()
    for expected_index, seed in enumerate(seeds):
        if not isinstance(seed, FuzzyBoundarySeed):
            raise ValueError("fuzzy boundary lattice seed has the wrong type")
        if seed.seed_index != expected_index:
            raise ValueError(
                "fuzzy boundary lattice seed indexes are not contiguous"
            )
        if (
                not isinstance(seed.start, int)
                or isinstance(seed.start, bool)
                or not isinstance(seed.end, int)
                or isinstance(seed.end, bool)
                or seed.start < previous_end
                or seed.end <= seed.start
                or seed.start < lattice.window_start
                or seed.end > lattice.window_end):
            raise ValueError(
                "fuzzy boundary lattice seed interval is malformed"
            )
        source_slice = source_binding.sequence[seed.start:seed.end]
        if source_slice != seed.sequence:
            raise ValueError(
                "fuzzy boundary lattice seed does not replay original source"
            )
        if _sha256_text(source_slice) != seed.slice_sha256:
            raise ValueError(
                "fuzzy boundary lattice seed slice identity is stale"
            )
        counts = tuple(seed.fuzzy_counts)
        if not lattice.edit_limits.accepts(counts):
            raise ValueError(
                "fuzzy boundary lattice seed exceeds its edit limits"
            )
        expected_length = (
            len(lattice.repeat_candidate) + counts[1] - counts[2]
        )
        if len(source_slice) != expected_length:
            raise ValueError(
                "fuzzy boundary lattice seed length contradicts edit counts"
            )
        expected_id = _seed_id(
            source_binding.sha256,
            lattice.repeat_sha256,
            seed.seed_index,
            seed.start,
            seed.end,
            counts,
        )
        if seed.seed_id != expected_id or seed.seed_id in observed_ids:
            raise ValueError(
                "fuzzy boundary lattice seed identity is inconsistent"
            )
        observed_ids.add(seed.seed_id)
        previous_end = seed.end


def _domain_intervals(
        seed, window_start, window_end, repeat_length, limits, radius):
    del limits
    baseline_length = seed.end - seed.start
    deletion_shortened = (
        baseline_length < repeat_length
        and seed.fuzzy_counts[2] > 0
    )
    if deletion_shortened:
        domain_kind = "deletion_shortened_restoration"
        start_min = max(window_start, seed.start - radius)
        start_max = min(window_end - 1, seed.start + radius)
        end_min = max(window_start + 1, seed.end - radius)
        end_max = min(window_end, seed.end + radius)
        minimum_length = baseline_length
        maximum_length = repeat_length
        intervals = tuple(
            (start, end)
            for start in range(start_min, start_max + 1)
            for end in range(end_min, end_max + 1)
            if minimum_length <= end - start <= maximum_length
        )
    else:
        domain_kind = "fixed_baseline"
        start_min = seed.start
        start_max = seed.start
        end_min = seed.end
        end_max = seed.end
        minimum_length = baseline_length
        maximum_length = baseline_length
        intervals = ((seed.start, seed.end),)
    domain = FuzzyBoundaryDomain(
        seed.seed_id,
        seed.seed_index,
        domain_kind,
        start_min,
        start_max,
        end_min,
        end_max,
        minimum_length,
        maximum_length,
        len(intervals),
        _canonical_sha256(intervals),
    )
    return domain, intervals


def _translated_match(local_match, text, relative_start):
    fuzzy_changes = tuple(
        [relative_start + position for position in positions]
        for positions in local_match.fuzzy_changes
    )
    return FuzzyBoundaryMatch(
        text,
        relative_start,
        relative_start + len(local_match.group()),
        local_match.fuzzy_counts,
        fuzzy_changes,
        getattr(local_match, "deletion_position_shift", None),
    )


def _node_from_match(
        repeat, source_sha256, repeat_sha256, absolute_start,
        absolute_end, sequence, seed_indexes, local_match):
    normalized_changes = tuple(
        tuple(positions)
        for positions in normalize_fuzzy_changes(local_match)
    )
    gapped_pattern, gapped_match = build_pairwise_fuzzy_alignment(
        repeat, sequence, normalized_changes
    )
    fuzzy_counts = tuple(local_match.fuzzy_counts)
    alignment_sha256 = _canonical_sha256({
        "gapped_pattern": gapped_pattern,
        "gapped_match": gapped_match,
        "fuzzy_counts": fuzzy_counts,
        "normalized_changes": normalized_changes,
    })
    node_payload = {
        "source_sha256": source_sha256,
        "repeat_sha256": repeat_sha256,
        "start": absolute_start,
        "end": absolute_end,
        "slice_sha256": _sha256_text(sequence),
        "fuzzy_counts": fuzzy_counts,
        "normalized_changes": normalized_changes,
        "alignment_sha256": alignment_sha256,
        "seed_indexes": seed_indexes,
    }
    return FuzzyBoundaryNode(
        "fuzzy-boundary-node:" + _canonical_sha256(node_payload),
        absolute_start,
        absolute_end,
        sequence,
        node_payload["slice_sha256"],
        fuzzy_counts,
        normalized_changes,
        alignment_sha256,
        seed_indexes,
    )


def _path_count(choices, ceiling):
    count = 0

    def visit(seed_index, previous_end):
        nonlocal count
        if count > ceiling:
            return
        if seed_index == len(choices):
            count += 1
            return
        for node in choices[seed_index]:
            if node.start >= previous_end:
                visit(seed_index + 1, node.end)
                if count > ceiling:
                    return

    visit(0, 0)
    return count


def _enumerate_node_paths(choices):
    paths = []
    current = []

    def visit(seed_index, previous_end):
        if seed_index == len(choices):
            paths.append(tuple(current))
            return
        for node in choices[seed_index]:
            if node.start >= previous_end:
                current.append(node)
                visit(seed_index + 1, node.end)
                current.pop()

    visit(0, 0)
    return tuple(paths)


def _lattice_id_payload(values):
    payload = dict(values)
    payload.pop("lattice_id", None)
    return payload


def _make_lattice(**values):
    values["lattice_id"] = ""
    lattice_id = "fuzzy-boundary-lattice:" + _canonical_sha256(
        _lattice_id_payload(values)
    )
    values["lattice_id"] = lattice_id
    return FuzzyBoundaryLattice(**values)


def _base_lattice_values(
        source_binding, window_start, window_end, repeat, limits, policy,
        seeds, domains):
    return {
        "schema": FUZZY_BOUNDARY_LATTICE_SCHEMA,
        "source_length": source_binding.length,
        "source_sha256": source_binding.sha256,
        "window_start": window_start,
        "window_end": window_end,
        "window_sha256": _sha256_text(
            source_binding.sequence[window_start:window_end]
        ),
        "repeat_candidate": repeat,
        "repeat_sha256": _sha256_text(repeat),
        "matcher_contract": MATCHER_CONTRACT,
        "matcher_version": regex.__version__,
        "edit_limits": limits,
        "policy": policy,
        "seeds": seeds,
        "domains": domains,
    }


def _abstained_plan(
        base_values, status, reason, planned_count,
        evaluated_count=0, accepted_count=0, evaluation_sha256="",
        interval_complete=False, node_complete=False, nodes=()):
    lattice = _make_lattice(
        **base_values,
        status=status,
        status_reasons=(reason,),
        alternatives_complete=False,
        interval_evaluation_complete=interval_complete,
        node_universe_complete=node_complete,
        path_universe_complete=False,
        planned_interval_count=planned_count,
        evaluated_interval_count=evaluated_count,
        accepted_interval_count=accepted_count,
        evaluation_sha256=evaluation_sha256,
        nodes=tuple(nodes),
        compatibility_edges=(),
        paths=(),
        baseline_path_id=None,
    )
    return FuzzyBoundaryAlternativePlan(lattice, ())


def _build_from_seeds(
        repeat, source_binding, window_start, window_end, seeds,
        edit_limits, policy):
    if policy.boundary_radius > edit_limits.max_total_errors:
        raise ValueError(
            "boundary radius cannot exceed the unchanged total-error limit"
        )

    domains = []
    interval_seed_indexes = {}
    if len(seeds) <= policy.max_seed_occurrences:
        for seed in seeds:
            domain, intervals = _domain_intervals(
                seed,
                window_start,
                window_end,
                len(repeat),
                edit_limits,
                policy.boundary_radius,
            )
            domains.append(domain)
            for interval in intervals:
                interval_seed_indexes.setdefault(interval, []).append(
                    seed.seed_index
                )

    base_values = _base_lattice_values(
        source_binding,
        window_start,
        window_end,
        repeat,
        edit_limits,
        policy,
        seeds,
        tuple(domains),
    )
    if len(seeds) > policy.max_seed_occurrences:
        return _abstained_plan(
            base_values,
            "abstained_seed_cap",
            "seed_occurrence_ceiling_reached",
            0,
        )

    intervals = tuple(sorted(interval_seed_indexes))
    planned_count = len(intervals)
    if planned_count > policy.max_interval_evaluations:
        return _abstained_plan(
            base_values,
            "abstained_interval_cap",
            "interval_evaluation_ceiling_reached",
            planned_count,
        )

    source = source_binding.sequence
    text = source[window_start:window_end]
    repeat_sha256 = _sha256_text(repeat)
    expression = (
        "(?e)(?:{}){}".format(
            regex.escape(repeat), edit_limits.pattern_suffix()
        )
    )
    compiled = regex.compile(expression)
    nodes = []
    matches_by_node_id = {}
    evaluations = []
    for start, end in intervals:
        sequence = source[start:end]
        slice_sha256 = _sha256_text(sequence)
        local_match = compiled.fullmatch(sequence)
        if local_match is None:
            evaluations.append({
                "start": start,
                "end": end,
                "slice_sha256": slice_sha256,
                "status": "rejected",
            })
            continue
        counts = tuple(local_match.fuzzy_counts)
        if not edit_limits.accepts(counts):
            raise RuntimeError(
                "forced-boundary matcher exceeded its declared edit limits"
            )
        seed_indexes = tuple(interval_seed_indexes[(start, end)])
        node = _node_from_match(
            repeat,
            source_binding.sha256,
            repeat_sha256,
            start,
            end,
            sequence,
            seed_indexes,
            local_match,
        )
        nodes.append(node)
        relative_start = start - window_start
        matches_by_node_id[node.node_id] = _translated_match(
            local_match, text, relative_start
        )
        evaluations.append({
            "start": start,
            "end": end,
            "slice_sha256": slice_sha256,
            "status": "accepted",
            "node_id": node.node_id,
        })

    nodes = tuple(sorted(nodes, key=lambda node: (
        node.start,
        node.end,
        node.fuzzy_counts,
        node.normalized_changes,
        node.node_id,
    )))
    evaluation_sha256 = _canonical_sha256(evaluations)
    if len(nodes) > policy.max_materialized_nodes:
        return _abstained_plan(
            base_values,
            "abstained_node_cap",
            "node_materialization_ceiling_reached",
            planned_count,
            evaluated_count=planned_count,
            accepted_count=len(nodes),
            evaluation_sha256=evaluation_sha256,
            interval_complete=True,
        )

    nodes_by_geometry = {
        (node.start, node.end): node for node in nodes
    }
    baseline_nodes = tuple(
        nodes_by_geometry.get((seed.start, seed.end)) for seed in seeds
    )
    if any(node is None for node in baseline_nodes):
        return _abstained_plan(
            base_values,
            "abstained_baseline_replay",
            "baseline_geometry_not_reproduced_by_forced_matcher",
            planned_count,
            evaluated_count=planned_count,
            accepted_count=len(nodes),
            evaluation_sha256=evaluation_sha256,
            interval_complete=True,
            node_complete=True,
            nodes=nodes,
        )

    choices = tuple(tuple(
        node for node in nodes if seed.seed_index in node.seed_indexes
    ) for seed in seeds)
    path_count = _path_count(
        choices, policy.max_materialized_paths
    )
    if path_count > policy.max_materialized_paths:
        return _abstained_plan(
            base_values,
            "abstained_path_cap",
            "path_materialization_ceiling_reached",
            planned_count,
            evaluated_count=planned_count,
            accepted_count=len(nodes),
            evaluation_sha256=evaluation_sha256,
            interval_complete=True,
            node_complete=True,
            nodes=nodes,
        )

    node_paths = _enumerate_node_paths(choices)
    baseline_node_ids = tuple(node.node_id for node in baseline_nodes)
    paths = []
    for node_path in node_paths:
        node_ids = tuple(node.node_id for node in node_path)
        geometry = tuple(
            (node.start, node.end) for node in node_path
        )
        path_id = "fuzzy-boundary-path:" + _canonical_sha256({
            "source_sha256": source_binding.sha256,
            "repeat_sha256": repeat_sha256,
            "node_ids": node_ids,
        })
        paths.append(FuzzyBoundaryPath(
            path_id,
            (
                "baseline_geometry"
                if node_ids == baseline_node_ids
                else "forced_boundary_alternative"
            ),
            node_ids,
            _canonical_sha256(geometry),
        ))
    paths = tuple(sorted(paths, key=lambda path: (
        path.kind != "baseline_geometry", path.node_ids, path.path_id
    )))
    baseline_paths = tuple(
        path for path in paths if path.kind == "baseline_geometry"
    )
    if len(baseline_paths) != 1:
        raise RuntimeError(
            "complete fuzzy boundary lattice lost its unique baseline path"
        )
    baseline_path_id = baseline_paths[0].path_id
    compatibility_edges = tuple(sorted({
        (left, right)
        for path in paths
        for left, right in zip(path.node_ids, path.node_ids[1:])
    }))
    status = "complete" if paths else "complete_no_compatible_path"
    lattice = _make_lattice(
        **base_values,
        status=status,
        status_reasons=(),
        alternatives_complete=True,
        interval_evaluation_complete=True,
        node_universe_complete=True,
        path_universe_complete=True,
        planned_interval_count=planned_count,
        evaluated_interval_count=planned_count,
        accepted_interval_count=len(nodes),
        evaluation_sha256=evaluation_sha256,
        nodes=nodes,
        compatibility_edges=compatibility_edges,
        paths=paths,
        baseline_path_id=baseline_path_id,
    )
    alternatives = tuple(
        FuzzyBoundaryMatchPath(
            path.path_id,
            path.node_ids,
            tuple(
                matches_by_node_id[node_id] for node_id in path.node_ids
            ),
        )
        for path in paths
        if path.path_id != baseline_path_id
    )
    return FuzzyBoundaryAlternativePlan(lattice, alternatives)


def build_fuzzy_boundary_alternative_plan(
        repeat, text, baseline_matches, sequence_start=0, *,
        source_dna, edit_limits=None, policy=None):
    """Enumerate complete local boundary alternatives or explicitly abstain.

    ``baseline_matches`` use coordinates relative to ``text``.  All lattice
    intervals are zero-based, half-open coordinates in ``source_dna``.
    Returned match objects remain relative to ``text`` and can therefore be
    passed directly to the existing ``FuzzySearch(precomputed_matches=...)``
    adapter.
    """

    if not isinstance(repeat, str) or not repeat:
        raise ValueError("repeat candidate must be a non-empty string")
    if not isinstance(text, str):
        raise ValueError("fuzzy boundary search window must be a string")
    if (
            not isinstance(sequence_start, int)
            or isinstance(sequence_start, bool)
            or sequence_start < 0):
        raise ValueError("sequence start must be a non-negative integer")
    source_binding = bind_fuzzy_boundary_source(source_dna)
    window_end = sequence_start + len(text)
    if (
            window_end > source_binding.length
            or source_binding.sequence[sequence_start:window_end] != text):
        raise ValueError(
            "fuzzy boundary search window does not replay original source"
        )
    limits = edit_limits or FuzzyEditLimits()
    search_policy = policy or FuzzyBoundarySearchPolicy()
    if not isinstance(limits, FuzzyEditLimits):
        raise TypeError("edit_limits must be a FuzzyEditLimits record")
    if not isinstance(search_policy, FuzzyBoundarySearchPolicy):
        raise TypeError(
            "policy must be a FuzzyBoundarySearchPolicy record"
        )
    seeds = _source_seeds(
        repeat,
        text,
        baseline_matches,
        sequence_start,
        source_binding,
        limits,
    )
    return _build_from_seeds(
        repeat,
        source_binding,
        sequence_start,
        window_end,
        seeds,
        limits,
        search_policy,
    )




__all__ = [
    "BoundFuzzyBoundarySource",
    "FUZZY_BOUNDARY_LATTICE_SCHEMA",
    "FuzzyBoundaryAlternativePlan",
    "FuzzyBoundaryDomain",
    "FuzzyBoundaryLattice",
    "FuzzyBoundaryMatch",
    "FuzzyBoundaryMatchPath",
    "FuzzyBoundaryNode",
    "FuzzyBoundaryPath",
    "FuzzyBoundarySearchPolicy",
    "FuzzyBoundarySeed",
    "FuzzyEditLimits",
    "MATCHER_CONTRACT",
    "bind_fuzzy_boundary_source",
    "build_fuzzy_boundary_alternative_plan",
]
