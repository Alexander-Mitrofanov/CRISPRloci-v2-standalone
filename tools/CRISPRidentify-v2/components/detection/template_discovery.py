"""Bounded whole-source literal scans projected into source-local tasks."""

import hashlib
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass

import regex

from components.detection import literal_fuzzy
from components.detection.fuzzy_alignment import (
    _deletion_position_shift, MODERN_REGEX_VERSION,
    normalize_fuzzy_changes,
    snapshot_owned_fuzzy_match_state,
)
from components.detection.fuzzy_alignment_merge import build_pairwise_fuzzy_alignment
from components.detection.repeat_templates import (
    DEFAULT_POLICY, POLICY, TemplateDiscoveryError, TemplatePolicy, check_limit, digest,
    panel_manifest, parse_repeat_template_panel, policy_for_engine, reverse_complement,
    validate_template_source,
)
from components.detection.search_contract import FUZZY_ERROR, MAX_LOCAL_CHAIN_SPACER_LENGTH


@dataclass(frozen=True)
class TemplateHit:
    source_sha256: str
    pattern_sha256: str
    canonical: str
    pattern: str
    start: int
    end: int
    sequence: str
    fuzzy_counts: tuple
    fuzzy_changes: tuple
    deletion_position_shift: int


class WindowedLiteralMatch:
    """Owned match adapter. Coordinates, including raw deletions, translate equally."""

    def __init__(self, hit, local_source, offset):
        self._source = local_source
        self._start, self._end = hit.start - offset, hit.end - offset
        if not 0 <= self._start < self._end <= len(local_source):
            raise TemplateDiscoveryError("window clips a template match")
        if local_source[self._start:self._end] != hit.sequence:
            raise TemplateDiscoveryError("window match sequence differs")
        self.fuzzy_counts = hit.fuzzy_counts
        self.fuzzy_changes = tuple([p - offset for p in group] for group in hit.fuzzy_changes)
        self.deletion_position_shift = hit.deletion_position_shift

    def start(self):
        return self._start

    def end(self):
        return self._end

    def span(self):
        return self._start, self._end

    def group(self):
        return self._source[self._start:self._end]

    def _fuzzy_alignment_cache_snapshot(self):
        return snapshot_owned_fuzzy_match_state(self.fuzzy_counts, self.fuzzy_changes)


def validate_hit(hit, source, source_sha256):
    if (hit.source_sha256 != source_sha256
            or hit.pattern_sha256 != hashlib.sha256(hit.pattern.encode("ascii")).hexdigest()
            or hit.canonical != min(hit.pattern, reverse_complement(hit.pattern))
            or not 15 <= len(hit.pattern) <= 63 or set(hit.pattern) - set("ACGT")):
        raise TemplateDiscoveryError("template hit identity differs")
    if (type(hit.start) is not int or type(hit.end) is not int
            or not 0 <= hit.start < hit.end <= len(source)
            or source[hit.start:hit.end] != hit.sequence):
        raise TemplateDiscoveryError("template hit span/source differs")
    if (len(hit.fuzzy_counts) != 3 or len(hit.fuzzy_changes) != 3
            or any(type(c) is not int or not 0 <= c <= 3 for c in hit.fuzzy_counts)
            or sum(hit.fuzzy_counts) > 6
            or any(len(group) != c for c, group in zip(hit.fuzzy_counts, hit.fuzzy_changes))
            or any(type(p) is not int for group in hit.fuzzy_changes for p in group)
            or type(hit.deletion_position_shift) is not int
            or hit.deletion_position_shift not in (0, 1)):
        raise TemplateDiscoveryError("template hit edit contract differs")
    try:
        match = WindowedLiteralMatch(hit, source, 0)
        build_pairwise_fuzzy_alignment(hit.pattern, hit.sequence, normalize_fuzzy_changes(match))
    except ValueError as error:
        raise TemplateDiscoveryError("template hit edit replay failed: " + str(error)) from error


def scan_pattern(source, pattern, engine, cap):
    try:
        if engine == "native":
            return literal_fuzzy.finditer(pattern, source, max_matches=cap)
        if engine != "full" or regex.__version__ != MODERN_REGEX_VERSION:
            raise TemplateDiscoveryError("template full mode requires pinned regex " + MODERN_REGEX_VERSION)
        matches = []
        for match in regex.finditer("(?e)(?:" + pattern + ")" + FUZZY_ERROR, source):
            if len(matches) == cap:
                raise TemplateDiscoveryError("template selected-hit limit exceeded")
            matches.append(match)
        return matches
    except Exception as error:
        raise TemplateDiscoveryError("template {} scan failed: {}".format(engine, error)) from error


@dataclass(frozen=True)
class TemplateTask:
    id: str
    pattern: str
    start: int
    end: int
    hits: tuple

    def matches(self, source):
        local = source[self.start:self.end]
        return tuple(WindowedLiteralMatch(hit, local, self.start) for hit in self.hits)


def plan_windows(source, hits, flank, policy=DEFAULT_POLICY):
    if type(flank) is not int or flank < 0:
        raise TemplateDiscoveryError("template flank must be a non-negative integer")
    if not hits:
        return ()
    chains = [[hits[0]]]
    for previous, hit in zip(hits, hits[1:]):
        if hit.pattern != previous.pattern or hit.start < previous.end:
            raise TemplateDiscoveryError("template hits must be ordered, nonoverlapping and same-oriented")
        if hit.start - previous.end > MAX_LOCAL_CHAIN_SPACER_LENGTH:
            chains.append([])
        chains[-1].append(hit)
    tasks = {}
    starts, ends = [h.start for h in hits], [h.end for h in hits]
    for chain in chains:
        start, end = max(0, chain[0].start - flank), min(len(source), chain[-1].end + flank)
        check_limit("window_span", end - start, policy)
        contained = tuple(hits[bisect_left(starts, start):bisect_right(ends, end)])
        check_limit("hits_per_task", len(contained), policy)
        task_id = digest([hits[0].pattern, start, end, [asdict(h) for h in contained]])
        tasks[task_id] = TemplateTask(task_id, hits[0].pattern, start, end, contained)
        check_limit("tasks", len(tasks), policy)
    return tuple(tasks.values())


def discover_template_tasks(source, panel, policy=None, engine="native", flank=100):
    """Yield bounded tasks one pattern at a time; each scan owns its immutable hits."""
    policy = policy or policy_for_engine(engine)
    validate_template_source(source, panel, policy)
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    total_hits, task_count = 0, 0
    for pattern in panel.patterns:
        try:
            matches = scan_pattern(source, pattern, engine, min(policy.hits_per_pattern, policy.total_hits - total_hits))
            hits = tuple(TemplateHit(
                source_hash, hashlib.sha256(pattern.encode("ascii")).hexdigest(),
                min(pattern, reverse_complement(pattern)), pattern, m.start(), m.end(), m.group(),
                tuple(m.fuzzy_counts), tuple(tuple(g) for g in m.fuzzy_changes),
                _deletion_position_shift(m)) for m in matches)
            del matches
            for hit in hits:
                validate_hit(hit, source, source_hash)
            total_hits += len(hits)
            check_limit("total_hits", total_hits, policy)
            tasks = plan_windows(source, hits, flank, policy)
            task_count += len(tasks)
            check_limit("tasks", task_count, policy)
            yield pattern, hits, tasks
        except Exception as error:
            raise TemplateDiscoveryError(
                "{} source={} panel={} pattern={}: {}".format(
                    POLICY, source_hash, panel.raw_sha256, pattern, error)) from error


def discovery_ledger(source, panel, engine, flank, policy=None):
    policy = policy or policy_for_engine(engine)
    return {"policy": POLICY, "engine": engine, "panel": panel_manifest(panel),
            "source_sha256": hashlib.sha256(source.encode("ascii")).hexdigest(),
            "source_length": len(source), "limits": asdict(policy), "flank": flank,
            "patterns": [], "tasks": [], "births": []}


def validate_discovery_ledger(ledger, source, configuration):
    """Replay direct discovery evidence before canonical report publication."""
    try:
        if configuration is None:
            if ledger is not None:
                raise ValueError("discovery ledger lacks configured panel")
            return
        if not isinstance(ledger, dict):
            raise ValueError("configured panel lacks discovery ledger")
        for key in ("policy", "engine", "panel", "limits"):
            if ledger[key] != configuration[key]:
                raise ValueError("discovery configuration differs: " + key)
        if ledger["policy"] != POLICY or ledger["engine"] not in ("native", "full"):
            raise ValueError("unsupported discovery policy/engine")
        manifest = ledger["panel"]
        if (any(not isinstance(manifest[k], str) or len(manifest[k]) != 64
                or set(manifest[k]) - set("0123456789abcdef")
                for k in ("raw_sha256", "normalized_sha256"))
                or type(manifest["raw_bytes"]) is not int or manifest["raw_bytes"] < 1):
            raise ValueError("panel byte identity is malformed")
        rebuilt = panel_manifest(parse_repeat_template_panel("".join(
            ">{}\n{}\n".format(a["id"], a["sequence"]) for a in manifest["aliases"]).encode("ascii")))
        for key in ("aliases", "canonical_sequences", "patterns", "normalization", "normalized_sha256"):
            if rebuilt[key] != manifest[key]:
                raise ValueError("panel normalization differs")
        policy = TemplatePolicy(**ledger["limits"])
        check_limit("raw_bytes", manifest["raw_bytes"], policy)
        for name in ("source_length", "hit_count", "task_count", "candidate_count", "flank"):
            if type(ledger[name]) is not int or ledger[name] < 0:
                raise ValueError("discovery count is malformed: " + name)
        source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
        if ledger["source_sha256"] != source_hash or ledger["source_length"] != len(source):
            raise ValueError("discovery source differs")
        check_limit("source_length", len(source), policy)
        check_limit("base_pattern_positions", len(source) * len(manifest["patterns"]), policy)
        if [p["pattern"] for p in ledger["patterns"]] != manifest["patterns"]:
            raise ValueError("incomplete discovery pattern inventory")
        expected_tasks, total_hits = [], 0
        for record in ledger["patterns"]:
            hits = tuple(TemplateHit(**dict(h, fuzzy_counts=tuple(h["fuzzy_counts"]),
                                           fuzzy_changes=tuple(tuple(g) for g in h["fuzzy_changes"])))
                         for h in record["hits"])
            check_limit("hits_per_pattern", len(hits), policy)
            total_hits += len(hits)
            for hit in hits:
                if hit.pattern != record["pattern"]:
                    raise ValueError("hit pattern differs")
                validate_hit(hit, source, source_hash)
            for task in plan_windows(source, hits, ledger["flank"], policy):
                expected_tasks.append({"id": task.id, "pattern": task.pattern,
                                       "start": task.start, "end": task.end, "hit_count": len(task.hits)})
        check_limit("total_hits", total_hits, policy)
        check_limit("tasks", len(expected_tasks), policy)
        if ledger["tasks"] != expected_tasks or ledger["hit_count"] != total_hits:
            raise ValueError("discovery task/hit inventory differs")
        if ledger["task_count"] != len(expected_tasks) or ledger["candidate_count"] != len(ledger["births"]):
            raise ValueError("discovery counts differ")
        by_id = {t["id"]: t for t in expected_tasks}
        for birth in ledger["births"]:
            task = by_id[birth["task_id"]]
            starts, repeats, gapped, spacers = birth["representation"]
            if (birth["id"] != digest(birth["representation"]) or birth["pattern"] != task["pattern"]
                    or not starts or len(starts) != len(repeats) or len(gapped) != len(repeats)
                    or len(spacers) != len(repeats) - 1 or len(set(map(len, gapped))) != 1):
                raise ValueError("discovery birth representation differs")
            for i, (start, repeat, aligned) in enumerate(zip(starts, repeats, gapped)):
                end = start + len(repeat)
                if (type(start) is not int or not repeat or not task["start"] <= start < end <= task["end"]
                        or source[start:end] != repeat or aligned.replace("-", "").replace(" ", "") != repeat):
                    raise ValueError("discovery birth does not reconstruct source")
                if i + 1 < len(starts) and (end > starts[i + 1] or source[end:starts[i + 1]] != spacers[i]):
                    raise ValueError("discovery birth spacer differs")
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as error:
        raise TemplateDiscoveryError("invalid repeat discovery ledger: " + str(error)) from error
