import json
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event
import unittest
from unittest.mock import patch

from components import components_evaluated_arrays_enhancement as enhancement
from components.fast_run_rescue import candidate_score_fingerprint
from tests.test_repeat_boundary_repair import synthetic


def source_candidate(repeats, spacer="GATC" * 8, prefix="ACGT"):
    source = prefix
    starts, ends = [], []
    for repeat in repeats:
        starts.append(len(source))
        source += repeat
        ends.append(len(source))
        source += spacer
    return source, enhancement.build_source_backed_candidate_from_boundaries(
        source, None, starts, ends,
    )


def state_bytes(candidate):
    return json.dumps(vars(candidate), sort_keys=True, separators=(",", ":")).encode()


class SourceAlignmentCacheTests(unittest.TestCase):
    def setUp(self):
        enhancement.clear_source_alignment_caches()

    def tearDown(self):
        enhancement.clear_source_alignment_caches()

    def test_cold_warm_and_evicted_state_and_score_inputs_are_identical(self):
        for kwargs in ({"right": False}, {}, {"terminal_gap": True},
                       {"internal_indel": True}, {"left": True}):
            enhancement.clear_source_alignment_caches()
            source, cold = synthetic(**kwargs)
            starts = cold.list_repeat_starts
            ends = [start + len(repeat) for start, repeat in zip(starts, cold.list_repeats)]
            expected = state_bytes(cold)
            expected_fingerprint = candidate_score_fingerprint(cold)
            warm = enhancement.build_source_backed_candidate_from_boundaries(source, None, starts, ends)
            self.assertEqual(expected, state_bytes(warm))
            self.assertEqual(expected_fingerprint, candidate_score_fingerprint(warm))
            self.assertIsNot(cold, warm)
            self.assertIsNot(cold.list_repeats_gaped, warm.list_repeats_gaped)
            # Changing one output must not poison any cached pure alignment.
            warm.list_repeats_gaped[0] = "N" * len(warm.list_repeats_gaped[0])
            again = enhancement.build_source_backed_candidate_from_boundaries(source, None, starts, ends)
            self.assertEqual(expected, state_bytes(again))
            self.assertGreater(enhancement.source_alignment_cache_info()["alignment"]["alignment_hits"], 0)
            # Eviction and reconstruction must preserve exactly the same state.
            with patch.object(enhancement, "_SOURCE_ALIGNMENT_CACHE_MAXSIZE", 1):
                enhancement._source_repeat_alignment(("A" * 20, "C" * 20))
                self.assertNotIn(tuple(cold.list_repeats), enhancement._SOURCE_ALIGNMENT_CACHE)
                evicted = enhancement.build_source_backed_candidate_from_boundaries(source, None, starts, ends)
            self.assertEqual(expected, state_bytes(evicted))
            self.assertEqual(expected_fingerprint, candidate_score_fingerprint(evicted))

    def test_complete_order_and_multiplicity_are_in_the_alignment_key(self):
        first = "GTTACCGGATCGTACCGGTTA"
        second = first[:8] + "C" + first[9:]
        third = first[:6] + first[7:]
        cases = ((first, first, second, third), (second, first, first, third),
                 (first, second, second, third), (first, second, third))
        for repeats in cases:
            _, candidate = source_candidate(repeats)
            self.assertEqual(list(repeats), candidate.list_repeats)
            expected = list(enhancement._source_repeat_alignment_uncached(repeats))
            self.assertEqual(expected, candidate.list_repeats_gaped)
        self.assertEqual(len(cases), enhancement.source_alignment_cache_info()["alignment"]["currsize"])
        self.assertTrue(all(tuple(repeats) in enhancement._SOURCE_ALIGNMENT_CACHE for repeats in cases))

    def test_cached_alignment_never_reuses_source_coordinates_or_spacers(self):
        repeats = ("GTTACCGGATCGTACCGGTTA", "GTTACCGGATGGTACCGGTTA") * 3
        _, first = source_candidate(repeats)
        source, second = source_candidate(repeats, spacer="CTGA" * 9, prefix="TTGGCCAATT")
        self.assertEqual(first.list_repeats_gaped, second.list_repeats_gaped)
        self.assertNotEqual(first.list_repeat_starts, second.list_repeat_starts)
        self.assertNotEqual(first.list_spacers, second.list_spacers)
        enhancement.validate_candidate_source_reconstruction(
            source, second.list_repeats, second.list_repeat_starts, second.list_spacers,
        )

    def test_edit_script_cache_values_are_immutable_and_results_defensive(self):
        for first, second in (("ACGTTAC", "ACGTAC"), ("ACGTAC", "ACGTTAC"),
                              ("ACGTAC", "AGGTAC"), ("ACGTAC", "ACGTAC")):
            expected = enhancement._minimal_edit_changes_uncached(first, second)
            result = enhancement._minimal_edit_changes(first, second)
            self.assertEqual(expected, result)
            result[0].append(999)
            self.assertEqual(expected, enhancement._minimal_edit_changes(first, second))
            self.assertEqual(expected, enhancement._minimal_edit_changes_uncached(first, second))
        cached = enhancement._cached_minimal_edit_changes("ACGTTAC", "ACGTAC")
        self.assertIsInstance(cached, tuple)
        self.assertTrue(all(isinstance(positions, tuple) for positions in cached))

    def test_pair_caches_bound_entries_and_eviction_changes_only_cost(self):
        for index in range(enhancement._SOURCE_DISTANCE_CACHE_MAXSIZE + 1):
            text = str(index)
            self.assertEqual(len(text), enhancement._source_repeat_edit_distance("", text))
        info = enhancement.source_alignment_cache_info()["distance"]
        self.assertEqual(enhancement._SOURCE_DISTANCE_CACHE_MAXSIZE, info["currsize"])
        self.assertEqual(enhancement._SOURCE_DISTANCE_CACHE_MAXSIZE, info["maxsize"])
        for index in range(enhancement._SOURCE_EDIT_CACHE_MAXSIZE + 1):
            text = str(index)
            self.assertEqual(enhancement._minimal_edit_changes_uncached("", text),
                             enhancement._minimal_edit_changes("", text))
        info = enhancement.source_alignment_cache_info()["edit_scripts"]
        self.assertEqual(enhancement._SOURCE_EDIT_CACHE_MAXSIZE, info["currsize"])
        self.assertEqual(enhancement._SOURCE_EDIT_CACHE_MAXSIZE, info["maxsize"])
        # The original evicted pair still has exactly the same result.
        self.assertEqual(1, enhancement._source_repeat_edit_distance("", "0"))
        self.assertEqual([[], [0], []], enhancement._minimal_edit_changes("", "0"))

    def test_alignment_bounds_and_oversize_bypass(self):
        for index in range(enhancement._SOURCE_ALIGNMENT_CACHE_MAXSIZE + 1):
            # Literal sequences are accepted by the pure alignment routine;
            # no biological candidate is admitted or scored by this test.
            enhancement._source_repeat_alignment(("A" + str(index),) * 2)
        info = enhancement.source_alignment_cache_info()["alignment"]
        self.assertEqual(enhancement._SOURCE_ALIGNMENT_CACHE_MAXSIZE, info["currsize"])
        self.assertEqual(1, info["alignment_evictions"])
        cases = (("_SOURCE_ALIGNMENT_CACHE_MAX_BASES", "alignment_oversize_key_bypasses"),
                 ("_SOURCE_ALIGNMENT_CACHE_MAX_REPEATS", "alignment_oversize_key_bypasses"),
                 ("_SOURCE_ALIGNMENT_CACHE_MAX_VALUE_BASES", "alignment_oversize_value_bypasses"))
        repeats = ("ACGT", "AGGT")
        expected = enhancement._source_repeat_alignment_uncached(repeats)
        for constant, counter in cases:
            enhancement.clear_source_alignment_caches()
            with patch.object(enhancement, constant, 1):
                observed = enhancement._source_repeat_alignment(repeats)
            self.assertEqual(expected, observed)
            info = enhancement.source_alignment_cache_info()["alignment"]
            self.assertEqual(0, info["currsize"])
            self.assertEqual(1, info[counter])
        enhancement.clear_source_alignment_caches()
        with patch.object(enhancement, "_SOURCE_PAIR_CACHE_MAX_LENGTH", 1):
            self.assertEqual(1, enhancement._source_repeat_edit_distance("ACGT", "AGGT"))
            self.assertEqual([[1], [], []], enhancement._minimal_edit_changes("ACGT", "AGGT"))
        self.assertEqual(0, enhancement.source_alignment_cache_info()["distance"]["currsize"])
        self.assertEqual(0, enhancement.source_alignment_cache_info()["edit_scripts"]["currsize"])

    def test_payload_exact_boundary_lru_eviction_and_reset(self):
        first, second, third = (("ACGT", "AGGT"), ("TCGT", "TGGT"), ("CCGT", "CGGT"))
        alignment = enhancement._source_repeat_alignment_uncached(first)
        size = enhancement._source_alignment_payload_bytes(first, alignment)
        with patch.object(enhancement, "_SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES", size - 1):
            self.assertEqual(alignment, enhancement._source_repeat_alignment(first))
            info = enhancement.source_alignment_cache_info()["alignment"]
            self.assertEqual(0, info["currsize"])
            self.assertEqual(0, info["retained_payload_bytes"])
            self.assertEqual(1, info["alignment_payload_oversize_bypasses"])
        enhancement.clear_source_alignment_caches()
        with patch.object(enhancement, "_SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES", size):
            enhancement._source_repeat_alignment(first)
            info = enhancement.source_alignment_cache_info()["alignment"]
            self.assertEqual(size, info["retained_payload_bytes"])
            self.assertEqual(size, info["peak_retained_payload_bytes"])
            enhancement._source_repeat_alignment(first)
            self.assertEqual(size, enhancement.source_alignment_cache_info()["alignment"]["retained_payload_bytes"])
        enhancement.clear_source_alignment_caches()
        with patch.object(enhancement, "_SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES", 2 * size):
            for repeats in (first, second, first, third):
                enhancement._source_repeat_alignment(repeats)
            self.assertEqual([first, third], list(enhancement._SOURCE_ALIGNMENT_CACHE))
            info = enhancement.source_alignment_cache_info()["alignment"]
            self.assertEqual(2 * size, info["retained_payload_bytes"])
            self.assertEqual(2 * size, info["peak_retained_payload_bytes"])
            self.assertEqual(1, info["alignment_byte_limit_evictions"])
            self.assertEqual(0, info.get("alignment_entry_limit_evictions", 0))
            enhancement._source_repeat_alignment(second)
            self.assertEqual([third, second], list(enhancement._SOURCE_ALIGNMENT_CACHE))
            self.assertEqual(2 * size, enhancement.source_alignment_cache_info()["alignment"]["retained_payload_bytes"])
        enhancement.clear_source_alignment_caches()
        info = enhancement.source_alignment_cache_info()["alignment"]
        self.assertEqual(0, info["currsize"])
        self.assertEqual(0, info["retained_payload_bytes"])
        self.assertEqual(0, info["peak_retained_payload_bytes"])
        self.assertEqual(0, info.get("alignment_evictions", 0))

    def test_duplicate_concurrent_insertions_count_payload_once(self):
        repeats = ("ACGT", "AGGT")
        original = enhancement._source_repeat_alignment_uncached
        expected = original(repeats)
        size = enhancement._source_alignment_payload_bytes(repeats, expected)
        barrier = Barrier(2)

        def synchronize(key):
            barrier.wait(timeout=5)
            return original(key)

        with patch.object(enhancement, "_source_repeat_alignment_uncached", side_effect=synchronize):
            with ThreadPoolExecutor(max_workers=2) as executor:
                results = list(executor.map(enhancement._source_repeat_alignment, (repeats, repeats)))
        self.assertEqual([expected, expected], results)
        info = enhancement.source_alignment_cache_info()["alignment"]
        self.assertEqual(1, info["currsize"])
        self.assertEqual(size, info["retained_payload_bytes"])
        self.assertEqual(size, info["peak_retained_payload_bytes"])
        self.assertEqual(1, info["alignment_duplicate_insertions"])

    def test_reset_during_computation_does_not_restore_stale_accounting(self):
        entered, release = Event(), Event()
        original = enhancement._source_repeat_alignment_uncached

        def synchronize(key):
            entered.set()
            if not release.wait(timeout=5):
                raise RuntimeError("test alignment was not released")
            return original(key)

        with patch.object(enhancement, "_source_repeat_alignment_uncached", side_effect=synchronize):
            with ThreadPoolExecutor(max_workers=1) as executor:
                future = executor.submit(enhancement._source_repeat_alignment, ("ACGT", "AGGT"))
                self.assertTrue(entered.wait(timeout=5))
                enhancement.clear_source_alignment_caches()
                release.set()
                future.result(timeout=5)
        info = enhancement.source_alignment_cache_info()["alignment"]
        self.assertEqual(0, info["currsize"])
        self.assertEqual(0, info["retained_payload_bytes"])
        self.assertEqual(0, info["peak_retained_payload_bytes"])
        self.assertEqual(1, info["alignment_reset_race_bypasses"])

    def test_entry_and_payload_limits_are_independent(self):
        with patch.object(enhancement, "_SOURCE_ALIGNMENT_CACHE_MAXSIZE", 1):
            enhancement._source_repeat_alignment(("ACGT", "AGGT"))
            enhancement._source_repeat_alignment(("TCGT", "TGGT"))
            info = enhancement.source_alignment_cache_info()["alignment"]
            self.assertEqual(1, info["currsize"])
            self.assertEqual(1, info["alignment_entry_limit_evictions"])
            self.assertEqual(0, info.get("alignment_byte_limit_evictions", 0))
            self.assertEqual(sum(entry[1] for entry in enhancement._SOURCE_ALIGNMENT_CACHE.values()),
                             info["retained_payload_bytes"])


if __name__ == "__main__":
    unittest.main()
