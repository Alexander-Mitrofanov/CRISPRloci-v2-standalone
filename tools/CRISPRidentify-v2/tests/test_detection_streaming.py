"""Scheduling bounds must preserve every ordered candidate and native match."""
import os
import pickle
from types import SimpleNamespace
import unittest
from unittest import mock

from components import literal_fuzzy
from components.components_detection import ClusterSequence, FuzzySearch, RejectedFuzzyCandidateTopology
from components.module_detection import Detection, _initialize_fuzzy_source, _initialize_fuzzy_worker, _native_worker_pool, _ordered_bounded_results


def worker_runtime_probe(unused):
    from components import module_detection
    return {'cpus': sorted(os.sched_getaffinity(0)),
            'omp': os.environ['OMP_NUM_THREADS'],
            'blas': os.environ['OPENBLAS_NUM_THREADS'],
            'hashseed': os.environ.get('PYTHONHASHSEED'),
            'boundary_paths': module_detection._fuzzy_build_boundary_paths(),
            'source': module_detection._fuzzy_source_dna()}


class RecordingPool:
    def __init__(self):
        self.submitted = []
        self.collected = []
        self.maximum_pending = 0

    def apply_async(self, worker, arguments):
        task, = arguments
        self.submitted.append(task)
        self.maximum_pending = max(self.maximum_pending, len(self.submitted) - len(self.collected))

        def get():
            self.collected.append(task)
            return worker(task)
        return SimpleNamespace(get=get)


class DetectionStreamingTests(unittest.TestCase):
    def test_slow_first_task_cannot_enqueue_unbounded_completed_results(self):
        pool = RecordingPool()
        seen = []

        def worker(task):
            # A pool may complete every submitted later task while waiting
            # for this first result; none beyond the window may be submitted.
            if task == 0:
                self.assertEqual([0, 1], pool.submitted)
            return task * 7

        results = _ordered_bounded_results(pool, worker, range(17), window=2)
        for result in results:
            seen.append(result)
            self.assertLessEqual(len(pool.submitted) - len(pool.collected), 1)
        self.assertEqual([i * 7 for i in range(17)], seen)
        self.assertEqual(list(range(17)), pool.collected)
        self.assertEqual(2, pool.maximum_pending)
        self.assertEqual([], list(_ordered_bounded_results(pool, worker, [], window=2)))
        with self.assertRaises(ValueError):
            list(_ordered_bounded_results(pool, worker, [1], window=0))

    def test_worker_failure_is_propagated_without_scheduling_remaining_tasks(self):
        pool = RecordingPool()
        def fail(task):
            if task == 1:
                raise RuntimeError('worker failed')
            return task
        stream = _ordered_bounded_results(pool, fail, range(20), window=2)
        self.assertEqual(0, next(stream))
        with self.assertRaisesRegex(RuntimeError, 'worker failed'):
            next(stream)
        self.assertEqual([0, 1, 2], pool.submitted)

    def test_batch_collection_preserves_cluster_pattern_and_topology_order(self):
        clusters = [ClusterSequence('A' * 100, 0, 100, tuple(range(259))),
                    ClusterSequence('C' * 100, 100, 200, ()),
                    ClusterSequence('G' * 100, 200, 300, ('last',))]
        marker = RejectedFuzzyCandidateTopology(((2, 5),), (), 3, 1)

        def worker(arguments):
            return [[pattern, marker] for pattern in arguments[0]]

        outputs = []
        for parallel in (False, True):
            detector = Detection.__new__(Detection)
            detector.cluster_sequences = clusters
            detector.dict_fuzzy_crisprs = {}
            detector.dna = 'A' * 300
            detector.flag_parallel = parallel
            detector.flag_cpu = '2'
            pool = RecordingPool()
            manager = mock.MagicMock()
            manager.__enter__.return_value = pool
            with mock.patch('components.module_detection._native_worker_pool', return_value=manager), \
                    mock.patch.object(Detection, '_parallel_run_fuzzy_native_batch', side_effect=worker):
                detector._run_native_array_detection('unchanged')
            outputs.append(detector.dict_fuzzy_crisprs)
            self.assertEqual({'raw_rejected_candidates': 780,
                              'reduced_rejected_candidates': 260,
                              'interval_runs': 260}, detector.compact_fuzzy_topology_telemetry)
            self.assertEqual([], detector.dict_fuzzy_crisprs[clusters[1]])
            self.assertEqual([value for pattern in range(259) for value in (pattern, marker)],
                             detector.dict_fuzzy_crisprs[clusters[0]])
        self.assertEqual(outputs[0], outputs[1])

    @unittest.skipUnless(literal_fuzzy.AVAILABLE, 'native matcher not built')
    def test_real_spawn_matches_serial_and_cleans_workers_after_failure(self):
        repeat = 'ACGTACGTTGCAACGTTACGTTAGCT'
        source = 'T' * 40 + (repeat + 'C' * 30) * 4 + 'T' * 40
        tasks = [((repeat, 'ACGT'), source, 0, '{i<=3,d<=3,s<=3,i+d+s<=6}'),
                 ((repeat[1:],), source, 0, '{i<=3,d<=3,s<=3,i+d+s<=6}')]
        for boundary_paths in (False, True):
            with self.subTest(boundary_paths=boundary_paths), mock.patch.dict(os.environ, PYTHONHASHSEED='0'):
                _initialize_fuzzy_source(source, boundary_paths)
                expected = [Detection._parallel_run_fuzzy_native_batch(task) for task in tasks]
                with _native_worker_pool(2, initializer=_initialize_fuzzy_worker,
                                         initargs=(source, boundary_paths)) as pool:
                    workers = list(pool._pool)
                    actual = list(_ordered_bounded_results(
                        pool, Detection._parallel_run_fuzzy_native_batch, tasks, window=2))
                    runtime = pool.apply(worker_runtime_probe, (None,))
                    self.assertEqual({'cpus': sorted(os.sched_getaffinity(0)), 'omp': '1', 'blas': '1',
                                      'hashseed': '0', 'boundary_paths': boundary_paths, 'source': source}, runtime)
                    self.assertEqual(len(expected), len(actual))
                    for before_batch, after_batch in zip(expected, actual):
                        self.assertEqual(len(before_batch), len(after_batch))
                        for before, after in zip(before_batch, after_batch):
                            self.assertEqual([c.__getstate__() for c in before],
                                             [c.__getstate__() for c in after])
                    with self.assertRaisesRegex(ValueError, 'invalid arity'):
                        list(_ordered_bounded_results(pool, Detection._parallel_run_fuzzy_native_batch,
                                                      [(1, 2, 3, 4, 5)], window=2))
                self.assertTrue(all(not worker.is_alive() for worker in workers))

    @unittest.skipUnless(literal_fuzzy.AVAILABLE, 'native matcher not built')
    def test_native_iterator_matches_regex_and_owns_its_scratch_copies(self):
        patterns = ('ACGTACGTTGCAACGTTACG', 'TCGATCGATGCTAGCTACGA', 'G' * 20)
        source = 'TTTT' + (patterns[0] + 'C' * 30) * 12 + patterns[1] + 'T' * 30
        expected = [literal_fuzzy.signature(list(FuzzySearch.compile_pattern(p, '{i<=3,d<=3,s<=3,e<=6}').finditer(source)))
                    for p in patterns]
        iterator = literal_fuzzy.iter_match_batches(patterns, source)
        first = next(iterator)
        saved = literal_fuzzy.signature(first)
        tail = list(iterator)
        self.assertEqual(saved, literal_fuzzy.signature(first))
        observed = [saved] + [literal_fuzzy.signature(matches) for matches in tail]
        self.assertEqual(expected, observed)
        self.assertEqual(observed, [literal_fuzzy.signature(m) for m in literal_fuzzy.finditer_many(patterns, source)])
        self.assertEqual([], literal_fuzzy.finditer_many([], source))
        with self.assertRaises(literal_fuzzy.NativeMatcherUnavailable):
            list(literal_fuzzy.iter_match_batches([patterns[0], 'unsupported!'], source))

    @unittest.skipUnless(literal_fuzzy.AVAILABLE, 'native matcher not built')
    def test_mixed_native_fallback_worker_matches_materialized_baseline(self):
        repeat = 'ACGTACGTTGCAACGTTACGTTAGCT'
        source = 'T' * 50 + (repeat + 'C' * 30) * 4 + 'T' * 50
        # The short pattern exercises the regex fallback interleaved with native.
        patterns = (repeat, 'ACGT', repeat[1:], 'G' * 25)
        arguments = (patterns, source, 0, '{i<=3,d<=3,s<=3,e<=6}')
        _initialize_fuzzy_source(source, True)
        original = literal_fuzzy.iter_match_batches
        with mock.patch.object(literal_fuzzy, 'iter_match_batches',
                               side_effect=lambda pats, text: list(original(pats, text))):
            expected = Detection._parallel_run_fuzzy_native_batch(arguments)
        actual = Detection._parallel_run_fuzzy_native_batch(arguments)
        self.assertEqual(pickle.dumps(expected), pickle.dumps(actual))
        self.assertEqual(len(patterns), len(actual))

    @unittest.skipUnless(literal_fuzzy.AVAILABLE, 'native matcher not built')
    def test_flat_worker_releases_sibling_cache_without_losing_paths_or_provenance(self):
        repeat = 'A' * 21
        region = 'C' + repeat + 'T' * 20 + repeat
        source = 'N' * 137 + region + 'N' * 20
        error = '{i<=3,d<=3,s<=3,i+d+s<=6}'
        arguments = ((repeat,), region, 137, error)
        _initialize_fuzzy_source(source, True)
        with mock.patch.object(FuzzySearch, 'detach_materialized_path_alternatives',
                               side_effect=lambda: None):
            before = Detection._parallel_run_fuzzy_native_batch(arguments)[0]
        self.assertTrue(any(c.occurrence_path_alternatives or c.boundary_path_alternatives
                            for c in before))
        after = Detection._parallel_run_fuzzy_native_batch(arguments)[0]
        self.assertEqual(len(before), len(after))
        for original, consumed in zip(before, after):
            expected = original.__getstate__()
            observed = consumed.__getstate__()
            for name in ('occurrence_path_alternatives', 'boundary_path_alternatives'):
                expected.pop(name)
                self.assertFalse(observed.pop(name))
            self.assertEqual(pickle.dumps(expected), pickle.dumps(observed))
            consumed._validate_source_path_provenance()
            if consumed.occurrence_path_lattice is not None:
                consumed._validate_occurrence_path_link()
            if consumed.boundary_path_lattice is not None:
                consumed._validate_boundary_path_link()
        direct = FuzzySearch(region, 137, repeat, error, source_dna=source)
        Detection._array_candidates_from_search(direct)
        self.assertGreater(len(direct.complete_occurrence_path_family()), 1)


if __name__ == '__main__':
    unittest.main()
