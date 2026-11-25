import copy
import dataclasses
import hashlib
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from components.detection import literal_fuzzy
from components.detection.fuzzy_alignment import normalize_fuzzy_changes
from components.detection.repeat_templates import (
    DEFAULT_POLICY, TemplateDiscoveryError, parse_repeat_template_panel, reverse_complement,
)
from components.detection.search import ClusterSequence, FuzzySearch
from components.detection.search_contract import FUZZY_ERROR
from components.detection.stage import Detection, _initialize_fuzzy_source
from components.detection.template_discovery import (
    TemplateHit, WindowedLiteralMatch, discover_template_tasks, plan_windows, scan_pattern, validate_hit,
)
from components.evidence.discovery_origins import inherit_discovery_origins, mark_discovery_birth
from components.pipeline import Pipeline
from components.refinement.template_integration import reconcile_template_candidates

MOTIF = "GTTTTTATCGTACCTATGAGGAATTGAAAC"
ROOT = Path(__file__).resolve().parents[1]


def panel(pattern=MOTIF):
    return parse_repeat_template_panel((">query\n" + pattern + "\n").encode())


def hit(source, start, sequence=MOTIF, changes=((), (), ()), counts=(0, 0, 0)):
    return TemplateHit(hashlib.sha256(source.encode()).hexdigest(), hashlib.sha256(MOTIF.encode()).hexdigest(),
                       min(MOTIF, reverse_complement(MOTIF)), MOTIF, start, start + len(sequence), sequence,
                       counts, changes, 1)


def raw_candidate(starts, template=False, repeats=None, role=True):
    repeats = repeats or [MOTIF] * len(starts)
    c = SimpleNamespace(list_absolute_start=starts, list_repeats=repeats,
                        list_gaped_repeats=repeats, list_spacers=["A" * (b-a-len(r))
                           for a, b, r in zip(starts, starts[1:], repeats)], number_errors=0,
                        source_path_classification_eligible=role)
    if template:
        mark_discovery_birth(c, "a" * 64, "b" * 64, MOTIF, "task")
    return c


class TemplateDiscoveryTests(unittest.TestCase):
    def test_native_regex_and_projection_edit_signatures(self):
        copies = (MOTIF[:8] + "A" + MOTIF[9:], MOTIF[:10] + "A" + MOTIF[10:],
                  MOTIF[:12] + MOTIF[13:], MOTIF[1:], MOTIF[:-1])
        for observed in copies:
            source = "N" * 80 + observed + "N" * 80
            native, regex = (scan_pattern(source, MOTIF, engine, 100) for engine in ("native", "full"))
            self.assertEqual([(m.span(), m.group(), m.fuzzy_counts, normalize_fuzzy_changes(m)) for m in native],
                             [(m.span(), m.group(), m.fuzzy_counts, normalize_fuzzy_changes(m)) for m in regex])
            self.assertTrue(native)
            for m in native:
                h = hit(source, m.start(), m.group(), tuple(tuple(g) for g in m.fuzzy_changes), m.fuzzy_counts)
                validate_hit(h, source, h.source_sha256)
                projected = WindowedLiteralMatch(h, source[50:130], 50)
                self.assertEqual(m.span(), (projected.start() + 50, projected.end() + 50))
                self.assertEqual(m.group(), projected.group())
                self.assertEqual(normalize_fuzzy_changes(m), normalize_fuzzy_changes(projected))
                self.assertEqual(m._fuzzy_alignment_cache_snapshot()[0], projected._fuzzy_alignment_cache_snapshot()[0])

    def test_leading_deletion_at_local_minus_one_is_preserved(self):
        source = "N" * 40 + MOTIF[1:]
        h = hit(source, 40, MOTIF[1:], ((), (), (39,)), (0, 0, 1))
        validate_hit(h, source, h.source_sha256)
        local = WindowedLiteralMatch(h, source[40:], 40)
        self.assertEqual([-1], local.fuzzy_changes[2])
        self.assertEqual([0], normalize_fuzzy_changes(local)[2])

    def test_tampered_hits_fail_source_and_alignment_replay(self):
        source = "N" * 50 + MOTIF
        original = hit(source, 50)
        for replacement in (dict(start=True), dict(end=79), dict(sequence="A" * 30),
                            dict(source_sha256="a" * 64), dict(pattern_sha256="a" * 64),
                            dict(fuzzy_counts=(1, 0, 0), fuzzy_changes=((50,), (), ())),
                            dict(fuzzy_counts=(True, 0, 0))):
            with self.subTest(replacement=replacement), self.assertRaises(TemplateDiscoveryError):
                validate_hit(dataclasses.replace(original, **replacement), source, original.source_sha256)

    def test_native_overflow_never_retries_or_grows_output(self):
        capacities = []
        def overflow(pattern, n, source, length, buffer, capacity):
            capacities.append((capacity, len(buffer)))
            return -1000
        with mock.patch.object(literal_fuzzy, "_FINDITER", side_effect=overflow):
            with self.assertRaises(literal_fuzzy.NativeMatchLimitExceeded):
                literal_fuzzy.finditer(MOTIF, "A" * 1000, max_matches=3)
        self.assertEqual([(4, 4)], capacities)
        source = MOTIF + "N" * 50 + MOTIF
        for engine in ("native", "full"):
            for cap in (1, 2, 3):
                if cap == 1:
                    with self.assertRaises(TemplateDiscoveryError):
                        scan_pattern(source, MOTIF, engine, cap)
                else:
                    self.assertEqual(2, len(scan_pattern(source, MOTIF, engine, cap)))

    def test_regex_stops_at_sentinel_and_engine_failure_is_visible(self):
        seen = []
        def iterator(*args):
            for i in range(100):
                seen.append(i)
                yield object()
        with mock.patch("components.detection.template_discovery.regex.finditer", side_effect=iterator):
            with self.assertRaises(TemplateDiscoveryError):
                scan_pattern("A" * 200, MOTIF, "full", 2)
        self.assertEqual([0, 1, 2], seen)
        with mock.patch.object(literal_fuzzy, "_FINDITER", side_effect=RuntimeError("broken")):
            with self.assertRaisesRegex(TemplateDiscoveryError, "broken"):
                list(discover_template_tasks("A" * 200, panel()))

    def test_work_budgets_checked_before_scanning(self):
        source = "N" * 100
        for field, count in (("source_length", 100), ("base_pattern_positions", 200), ("oriented_patterns", 2)):
            for cap in (count - 1, count, count + 1):
                policy = dataclasses.replace(DEFAULT_POLICY, **{field: cap})
                with mock.patch("components.detection.template_discovery.scan_pattern", return_value=[]) as scan:
                    if cap < count:
                        with self.assertRaisesRegex(TemplateDiscoveryError, field):
                            list(discover_template_tasks(source, panel(), policy))
                        scan.assert_not_called()
                    else:
                        list(discover_template_tasks(source, panel(), policy))
                        self.assertEqual(2, scan.call_count)

    def test_hit_task_and_span_caps(self):
        source = MOTIF + "N" * 20 + MOTIF
        for field, count in (("hits_per_pattern", 2), ("total_hits", 2), ("hits_per_task", 2), ("window_span", 80)):
            for cap in (count - 1, count, count + 1):
                policy = dataclasses.replace(DEFAULT_POLICY, **{field: cap})
                if cap < count:
                    with self.assertRaises(TemplateDiscoveryError):
                        list(discover_template_tasks(source, panel(), policy, flank=0))
                else:
                    list(discover_template_tasks(source, panel(), policy, flank=0))
        source = MOTIF + "N" * 200 + MOTIF
        for cap in (1, 2, 3):
            policy = dataclasses.replace(DEFAULT_POLICY, tasks=cap)
            if cap == 1:
                with self.assertRaisesRegex(TemplateDiscoveryError, "tasks"):
                    list(discover_template_tasks(source, panel(), policy, flank=0))
            else:
                self.assertEqual(2, sum(len(t) for _, _, t in discover_template_tasks(source, panel(), policy, flank=0)))

    def test_flanks_include_complete_neighbors_but_do_not_clip(self):
        source = "N" * 10 + MOTIF + "N" * 160 + MOTIF + "N" * 10
        hits = (hit(source, 10), hit(source, 200))
        tasks = plan_windows(source, hits, 175)
        self.assertEqual([1, 1], [len(t.hits) for t in tasks])
        tasks = plan_windows(source, hits, 190)
        self.assertTrue(all(len(t.hits) == 2 for t in tasks))

    def test_orientation_and_distant_loci_never_mix(self):
        source = MOTIF + "N" * 20 + reverse_complement(MOTIF) + "N" * 200 + MOTIF
        tasks = list(discover_template_tasks(source, panel(), flank=0))
        self.assertEqual(3, sum(len(t) for _, _, t in tasks))
        self.assertTrue(all(len(task.hits) == 1 for _, _, group in tasks for task in group))

    def test_precomputed_window_does_not_repeat_baseline_scan(self):
        source = "N" * 40 + MOTIF + "N" * 25 + MOTIF + "N" * 40
        for pattern, _, tasks in discover_template_tasks(source, panel()):
            for task in tasks:
                with mock.patch.object(FuzzySearch, "compile_pattern", side_effect=AssertionError("rescan")):
                    search = FuzzySearch(source[task.start:task.end], task.start, pattern, FUZZY_ERROR,
                                         source_dna=source, precomputed_matches=task.matches(source),
                                         _build_boundary_paths=False)
                    self.assertEqual([40, 95], search.list_absolute_start)

    def test_duplicate_only_preserves_self_keys_and_preferred_representative(self):
        first, duplicate = raw_candidate([80, 142, 204]), raw_candidate([80, 142, 204], True)
        a, b = ClusterSequence("A", 0, 300, ("self",)), ClusterSequence("A", 5, 305, (MOTIF,))
        result = reconcile_template_candidates({a: [first], b: [duplicate]})
        self.assertEqual([a], list(result))
        self.assertIs(first, result[a][0])
        self.assertTrue(first.discovery_origins)
        self.assertFalse(hasattr(first, "discovery_route"))

    def test_novel_overlap_competes_while_distant_and_touching_loci_stay_apart(self):
        a, b, c = (ClusterSequence("A", i, i+1000, (str(i),)) for i in range(3))
        baseline = raw_candidate([80, 142, 204])
        novel = raw_candidate([90, 152, 214], True)
        distant = raw_candidate([1000, 1062, 1124], True)
        result = reconcile_template_candidates({a: [baseline], b: [novel], c: [distant]})
        self.assertEqual([2, 1], [len(rows) for rows in result.values()])
        touching = raw_candidate([234, 296, 358], True)
        result = reconcile_template_candidates({a: [baseline], b: [touching]})
        self.assertEqual(2, len(result))

    def test_eligibility_roles_do_not_deduplicate(self):
        a, b = (ClusterSequence("A", i, i+1000, (str(i),)) for i in range(2))
        result = reconcile_template_candidates({a: [raw_candidate([80, 142, 204], role=False)],
                                                b: [raw_candidate([80, 142, 204], True)]})
        self.assertEqual(2, sum(map(len, result.values())))

    def test_derived_ancestry_keeps_birth_and_rejects_stale_geometry(self):
        parent = raw_candidate([80, 142, 204], True)
        derived = raw_candidate([81, 143, 205])
        inherit_discovery_origins(derived, (parent,))
        self.assertEqual(parent.discovery_origins[1], derived.discovery_origins[1])
        self.assertNotEqual(parent.discovery_origins[0], derived.discovery_origins[0])
        parent.list_absolute_start[0] += 1
        with self.assertRaises(ValueError):
            inherit_discovery_origins(derived, (parent,))

    def test_descriptor_collision_routing_and_delta_restore(self):
        source = "N" * 40 + MOTIF + "N" * 25 + MOTIF + "N" * 40
        d = Detection.__new__(Detection)
        d.dna, d.fuzzy_search_mode, d.fuzzy_search_flank = source, "native", 100
        d.flag_fuzzy_boundary_paths, d.repeat_template_panel = False, panel()
        key = ClusterSequence(source, 0, len(source), (MOTIF,))
        d.clusters, d.cluster_sequences, d.discovery_sources = [], [key], [[]]
        existing = raw_candidate([40, 95])
        d.dict_fuzzy_crisprs = {key: [existing]}
        _initialize_fuzzy_source(source, False)
        d._run_template_detection()
        self.assertEqual(1, len(d.cluster_sequences))
        self.assertIs(existing, d.dict_fuzzy_crisprs[key][0])
        self.assertEqual({MOTIF}, d.locus_seed_patterns(0, ((40, 125),)))
        self.assertEqual((0,), d.cluster_indices_for_intervals(((40, 125),))[0])
        pipeline = Pipeline.__new__(Pipeline)
        pipeline._raw_detection = d
        selected_candidate = SimpleNamespace(list_repeats=[MOTIF, MOTIF], consensus=MOTIF)
        for scope in ("lineage", "locus", "family"):
            selected = pipeline._fast_rescue_seed_scope(
                [{"source_interval": (40, 125), "trigger_candidate": selected_candidate}], (0,), scope)
            self.assertEqual({0: (MOTIF,)}, selected)
        before = copy.deepcopy(d.discovery_sources)
        d.fast_run_rescue_trim_profile = "full"
        d.compact_fuzzy_topology_telemetry = {}
        with mock.patch.object(d, "_run_array_detection"):
            d.run_start_end_delta((0,))
        self.assertEqual(before, d.discovery_sources)
