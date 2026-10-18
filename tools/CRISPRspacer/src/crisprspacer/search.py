"""Exhaustive full-length Hamming search with bounded exact-seed candidates.

For k substitutions, k+1 disjoint nonempty query blocks guarantee at least
one unchanged block. Every nominated start is verified against the full query.
Short queries (length <= k) enumerate all eligible starts explicitly.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import time

from .inputs import CATEGORIES, validate_inventory, validate_records
from .io import InputError, LimitError, canonical, digest, integer
from . import __version__

_ACGT = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Policy:
    substitutions: int = 0
    min_length: int = 1
    max_length: int = 1000
    categories: tuple[str, ...] = ("Bona-fide", "Possible")
    scope: str = "all"
    flank: int = 10
    array_padding: int = 0
    max_hits: int = 100_000
    max_candidates: int = 2_000_000
    max_pairs: int = 100_000
    max_evidence_bytes: int = 64 * 1024 * 1024
    max_work: int = 500_000_000
    timeout: int = 300

    def validate(self):
        for key, low, high in (("substitutions", 0, 2), ("min_length", 1, 10_000),
            ("max_length", self.min_length, 10_000), ("flank", 0, 1000), ("array_padding", 0, 10_000),
            ("max_hits", 1, 1_000_000), ("max_candidates", 1, 20_000_000),
            ("max_pairs", 1, 1_000_000), ("max_evidence_bytes", 1, 256 * 1024 * 1024),
            ("max_work", 1, 100_000_000_000), ("timeout", 1, 86_400)):
            integer(getattr(self, key), key, low, high)
        if self.scope not in {"all", "same-record", "same-assembly"}:
            raise InputError("Invalid target scope")
        if not isinstance(self.categories, (tuple, list)) or not self.categories or any(x not in CATEGORIES for x in self.categories):
            raise InputError("Invalid array category selection")
        if len(set(self.categories)) != len(self.categories):
            raise InputError("Duplicate category policy")

    def to_dict(self):
        result = asdict(self)
        result["categories"] = list(self.categories)
        return result


class Budget:
    def __init__(self, policy: Policy):
        self.policy = policy
        self.deadline = time.monotonic() + policy.timeout
        self.candidates = 0
        self.work = 0
        self.pairs = 0
        self.evidence_bytes = 0

    def check(self):
        if time.monotonic() > self.deadline:
            raise LimitError("Search deadline exceeded; no completed result was published")

    def charge(self, work: int = 0, candidates: int = 0):
        self.work += work
        self.candidates += candidates
        self.check()
        if self.work > self.policy.max_work:
            raise LimitError("Search work budget exceeded; no completed result was published")
        if self.candidates > self.policy.max_candidates:
            raise LimitError("Search candidate budget exceeded; no completed result was published")

    def pair(self, query_count: int):
        self.pairs += 1
        self.charge(work=max(1, query_count))
        if self.pairs > self.policy.max_pairs:
            raise LimitError("Query–target pair budget exceeded; no completed result was published")

    def retain(self, row: dict):
        self.evidence_bytes += len(canonical(row)) + 64
        self.check()
        if self.evidence_bytes > self.policy.max_evidence_bytes:
            raise LimitError("Evidence byte budget exceeded; no completed result was published")


def windows(target: dict, length: int) -> tuple[str, int]:
    sequence = target["sequence"]
    if length > len(sequence):
        return sequence, 0
    if target["topology"] == "circular":
        return sequence + sequence[:length-1], len(sequence)
    return sequence, len(sequence) - length + 1


def eligible_window_count(sequence: str, length: int, starts: int, budget: Budget) -> int:
    """Count ACGT-only windows without constructing a window per start."""
    run = total = 0
    budget.charge(len(sequence))
    for pos, letter in enumerate(sequence):
        run = run+1 if letter in "ACGT" else 0
        if run >= length and pos-length+1 < starts:
            total += 1
        if pos % 65536 == 0:
            budget.check()
    return total


def find_matches(query: str, target: dict, substitutions: int, budget: Budget) -> list[tuple[int, str, int]]:
    """Returns zero-based starts, relative strand (+/-/both), and Hamming distance."""
    sequence, starts = windows(target, len(query))
    if not starts:
        return []
    reverse = reverse_complement(query)
    choices = [(query, "both")] if query == reverse else [(query, "+"), (reverse, "-")]
    hits = []
    exact = substitutions == 0 and bool(query) and set(query) <= _ACGT
    for pattern, strand in choices:
        if exact:
            # A full-query exact seed already proves the complete ACGT window.
            # Stream overlapping starts; no candidate set or window copy is needed.
            budget.charge(len(sequence))
            start = sequence.find(pattern)
            while start != -1:
                budget.charge(len(pattern), candidates=1)
                hits.append((start, strand, 0))
                if len(hits) > budget.policy.max_hits:
                    raise LimitError("Hit budget exceeded; no completed result was published")
                start = sequence.find(pattern, start+1)
            continue
        nominated = set()
        if len(pattern) <= substitutions:
            budget.charge(candidates=starts)
            nominated.update(range(starts))
        else:
            block_count = substitutions + 1
            for block in range(block_count):
                left = block * len(pattern) // block_count
                right = (block+1) * len(pattern) // block_count
                seed = pattern[left:right]
                budget.charge(len(sequence))
                pos = sequence.find(seed)
                while pos != -1:
                    budget.charge(candidates=1)
                    start = pos-left
                    if 0 <= start < starts:
                        nominated.add(start)
                    pos = sequence.find(seed, pos+1)
        for start in sorted(nominated):
            budget.charge(len(pattern))
            candidate = sequence[start:start+len(pattern)]
            if any(base not in "ACGT" for base in candidate):
                continue
            distance = sum(a != b for a, b in zip(pattern, candidate))
            if distance <= substitutions:
                hits.append((start, strand, distance))
                if len(hits) > budget.policy.max_hits:
                    raise LimitError("Hit budget exceeded; no completed result was published")
    return sorted(hits)


def segments(start: int, length: int, size: int) -> list[dict]:
    end = start + length
    if end <= size:
        return [{"start": start+1, "end": end}]
    return [{"start": start+1, "end": size}, {"start": 1, "end": end-size}]


def forward_slice(target: dict, start: int, length: int) -> str:
    seq = target["sequence"]
    if target["topology"] == "circular":
        start %= len(seq)
        first = min(length, len(seq)-start)
        return seq[start:start+first] + seq[:length-first]
    return seq[max(0, start):max(0, min(len(seq), start+length))]


def query_status(query: dict, policy: Policy) -> str:
    if query["kind"] == "deletion":
        return "deletion"
    if query["category"] is not None and query["category"] not in policy.categories:
        return "category_excluded"
    if not policy.min_length <= len(query["sequence"]) <= policy.max_length:
        return "unsupported_length"
    if not set(query["sequence"]) <= _ACGT:
        return "ambiguous_query"
    return "eligible"


def overlaps(parts: list[dict], array: dict, padding: int, size: int, circular: bool) -> bool:
    low, high = array["start"]-padding, array["end"]+padding
    if circular and high-low+1 >= size:
        return True
    ranges = [(max(1, low), min(size, high))]
    if circular and low < 1:
        ranges.append((size+low, size))
    if circular and high > size:
        ranges.append((1, high-size))
    return any(part["start"] <= end and part["end"] >= start for part in parts for start, end in ranges)


def analyze(inventory: dict, targets: list[dict], policy: Policy, source_targets: bool = False) -> dict:
    policy.validate()
    validate_inventory(inventory)
    target_map = validate_records(targets)
    if not target_map:
        raise InputError("At least one target record is required")
    if source_targets and targets != inventory["sources"]:
        raise InputError("Source-target mode must use the exact source record inventory")
    if not source_targets and policy.scope != "all":
        raise InputError("Restricted scopes require self analysis with source-bound queries")
    sources = {r["id"]: r for r in inventory["sources"]}
    annotated_records = set(inventory["annotated_records"])
    ordered_targets = sorted(targets, key=lambda row: row["id"])
    arrays_by_record = {}
    for array in inventory["arrays"]:
        arrays_by_record.setdefault(array["source_id"], []).append(array)
    budget = Budget(policy)
    query_rows = []
    eligible = []
    for query in inventory["queries"]:
        budget.charge(1)
        status = query_status(query, policy)
        query_row = {**query, "status": status}
        budget.retain(query_row)
        query_rows.append(query_row)
        if status == "eligible":
            if policy.scope != "all" and query["source_id"] is None:
                raise InputError("Self search requires source-bound query occurrences")
            if policy.scope == "same-assembly" and sources[query["source_id"]]["assembly_id"] is None:
                raise InputError("Same-assembly scope requires an explicit assembly_id for every eligible query source")
            eligible.append(query)
    # Deduplicate computation, retaining every query occurrence in reported rows.
    groups = {}
    for query in eligible:
        groups.setdefault(query["sequence"], []).append(query)
    rows, coverage = [], []
    coverage_cache = {}
    for sequence, queries in sorted(groups.items()):
        sequence_sha256 = queries[0]["sha256"]
        for target in ordered_targets:
            budget.pair(len(queries))
            selected = queries if policy.scope == "all" else [query for query in queries if
                        (policy.scope == "same-record" and query["source_id"] == target["id"]) or
                        (policy.scope == "same-assembly" and
                         sources[query["source_id"]]["assembly_id"] == target["assembly_id"])]
            if not selected:
                continue
            key = (target["id"], len(sequence))
            if key not in coverage_cache:
                text, starts = windows(target, len(sequence))
                assessed = eligible_window_count(text, len(sequence), starts, budget) if starts else 0
                coverage_cache[key] = starts, assessed
            starts, assessed = coverage_cache[key]
            coverage_row = {"query_sha256": sequence_sha256, "target_id": target["id"],
                             "possible_windows": starts, "assessed_windows": assessed,
                             "ambiguous_windows": starts-assessed, "query_occurrence_count": len(selected)}
            budget.retain(coverage_row)
            coverage.append(coverage_row)
            hits = find_matches(sequence, target, policy.substitutions, budget)
            if len(rows) + len(hits)*len(selected) > policy.max_hits:
                raise LimitError("Occurrence-expanded hit budget exceeded; no completed result was published")
            for start, strand, distance in hits:
                budget.charge(2 * len(arrays_by_record.get(target["id"], [])))
                parts = segments(start, len(sequence), target["length"])
                known = source_targets and target["id"] in annotated_records
                actual = [a for a in arrays_by_record.get(target["id"], []) if known and
                          overlaps(parts, a, 0, target["length"], target["topology"] == "circular")]
                padded = [a for a in arrays_by_record.get(target["id"], []) if known and not actual and
                          overlaps(parts, a, policy.array_padding, target["length"], target["topology"] == "circular")]
                actual_ids = {array["id"] for array in actual}
                overlapping_ids = [array["id"] for array in actual or padded]
                forward = forward_slice(target, start, len(sequence))
                query_frame = reverse_complement(forward) if strand == "-" else forward
                flank_size = min(policy.flank, target["length"]-len(sequence)) if target["topology"] == "circular" else policy.flank
                left_flank = forward_slice(target, start-flank_size, flank_size)
                right_flank = forward_slice(target, (start+len(sequence)) % target["length"] if target["topology"] == "circular" else start+len(sequence), flank_size)
                for query in selected:
                    budget.check()
                    same_record = query["source_id"] == target["id"] if source_targets else None
                    source = sources.get(query["source_id"])
                    assembly = source["assembly_id"] if source else None
                    same_assembly = assembly == target["assembly_id"] if assembly is not None and target["assembly_id"] is not None else None
                    overlap_kind = "unknown"
                    if known:
                        overlap_kind = "source_array" if query["array_id"] in actual_ids else (
                            "other_array" if actual else "array_padding" if padded else "outside_reported_arrays")
                    excluded = bool(actual or padded)
                    row = {"query_id": query["id"], "target_id": target["id"],
                        "query_sha256": query["sha256"], "target_sha256": target["sha256"],
                        "target_start": start+1, "target_end": (start+len(sequence)-1) % target["length"]+1,
                        "target_segments": parts, "wraps_origin": len(parts) == 2,
                        "relative_strand": strand, "query_start": 1, "query_end": len(sequence),
                        "query_length": len(sequence), "query_coverage": 1.0,
                        "identities": len(sequence)-distance, "substitutions": distance, "gaps": 0,
                        "target_sequence_forward": forward,
                        "target_sequence_query_frame": query_frame,
                        "left_flank_forward": left_flank or None, "right_flank_forward": right_flank or None,
                        "same_source_record": same_record, "same_declared_assembly": same_assembly,
                        "array_overlap": overlap_kind, "overlapping_array_ids": overlapping_ids,
                        "excluded": excluded, "pam_status": "not_evaluated", "evalue": None, "confidence": None}
                    budget.retain(row)
                    rows.append(row)
    budget.check()
    rows.sort(key=lambda row: (row["query_id"], row["target_id"], row["target_start"], row["relative_strand"]))
    budget.check()
    for index, row in enumerate(rows, 1):
        row["id"] = f"match-{index:08d}"
    candidate_count = sum(not row["excluded"] for row in rows)
    outcome = "no_eligible_queries" if not eligible else "matches" if candidate_count else "no_candidate_matches"
    counts = {"query_occurrences": len(query_rows), "eligible_queries": len(eligible),
              "unique_eligible_sequences": len(groups), "targets": len(targets), "alignment_occurrences": len(rows),
              "candidate_matches": candidate_count, "excluded_array_matches": len(rows)-candidate_count,
              "assessed_windows": sum(row["assessed_windows"] for row in coverage),
              "ambiguous_windows": sum(row["ambiguous_windows"] for row in coverage)}
    warnings = ["Sequence matches do not establish functional targeting, active infection, immunity or host range.",
                "PAMs are not evaluated. Missing metadata remains unknown.",
                "Array exclusions cover reported arrays only; undetected or omitted arrays may remain."]
    if not counts["assessed_windows"] and eligible:
        outcome = "no_assessable_windows"
    budget.check()
    return {"schema": {"name": "crisprspacer.analysis", "version": "1.0.0"}, "tool_version": __version__,
            "execution_status": "completed", "retrieval": "exhaustive_full_length_hamming",
            "retrieval_scope": "eligible queries and unambiguous windows of the declared targets, at most one circular lap",
            "outcome": outcome, "search_complete": True, "results_truncated": False,
            "coordinate_system": "1-based inclusive forward target; ordered segments for origin-crossing matches",
            "policy": policy.to_dict(), "source_targets": source_targets,
            "inventory_sha256": digest(canonical(inventory)), "targets_sha256": digest(canonical(targets)),
            "counts": counts, "queries": query_rows, "coverage": coverage, "matches": rows, "warnings": warnings}
