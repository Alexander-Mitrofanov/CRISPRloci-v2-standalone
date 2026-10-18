"""Idle lattice storage must preserve scientific state and fail-closed replay."""
import dataclasses
import gc
import hashlib
import pickle
import unittest
import weakref
import zlib
from unittest import mock

from components.components_detection import ClusterSequence, FuzzySearch
from components.components_detection_refinement import SameStartEndFilter
from components.fuzzy_boundary_lattice import _canonical_sha256
from components.fuzzy_lattice_storage import CompressedLattice, LatticeDecodeCache
from components.module_detection import Detection, _initialize_fuzzy_source
from components.module_detection_refinement import DetectionRefinement

ERROR = '{i<=3,d<=3,s<=3,i+d+s<=6}'
REPEAT = 'ACGTACGTTGCAACGTTACGTTAGCT'
DELETED = REPEAT[:12] + REPEAT[13:]
SOURCE = 'T' * 20 + (DELETED + 'C' * 30) * 3 + 'T' * 20
PARAMETERS = {
    'param_min_repeats': 3, 'param_min_avg_repeat_length': 21,
    'param_max_avg_repeat_length': 55, 'param_min_avg_spacer_length': 18,
    'param_max_avg_spacer_length': 78, 'param_max_identical_spacers': 4,
    'param_max_identical_cluster_spacers': 3,
}


def family():
    _initialize_fuzzy_source(SOURCE, True)
    return Detection._parallel_run_fuzzy_native_batch(((REPEAT,), SOURCE, 0, ERROR))[0]


def occurrence_lattice():
    repeat = 'A' * 21
    source = 'C' + repeat + 'T' * 20 + repeat
    return FuzzySearch(source, 0, repeat, ERROR, source_dna=source).occurrence_path_lattice


def compact(candidates, capacity=4):
    cache, memo = LatticeDecodeCache(capacity), {}
    for candidate in candidates:
        candidate._compact_lattice_storage(cache, memo)
    return cache


class LatticeStorageTests(unittest.TestCase):
    def test_occurrence_alternatives_and_abstained_boundary_records_remain_complete(self):
        repeat = 'A' * 21
        source = 'C' + repeat + 'T' * 20 + repeat
        search = FuzzySearch(source, 0, repeat, ERROR, source_dna=source)
        candidates = list(search.complete_occurrence_path_family())
        self.assertGreater(len(candidates), 1)
        for candidate in candidates:
            candidate.detach_materialized_path_alternatives()
        expected = [candidate.__getstate__() for candidate in candidates]
        compact(candidates, capacity=1)
        for before, candidate in zip(expected, candidates):
            candidate._validate_source_path_provenance()
            self.assertEqual(before, candidate.__getstate__())
            self.assertIsNone(candidate._occurrence_path_validation_seal)
            self.assertIsNone(candidate._source_provenance_validation_seal)
        policy = dataclasses.replace(FuzzySearch._BOUNDARY_SEARCH_POLICY, max_materialized_paths=1)
        with mock.patch.object(FuzzySearch, '_BOUNDARY_SEARCH_POLICY', policy):
            abstained = family()
        self.assertEqual(1, len(abstained))
        self.assertEqual('abstained_path_cap', abstained[0].boundary_path_lattice.status)
        before = abstained[0].__getstate__()
        compact(abstained)
        abstained[0]._validate_source_path_provenance()
        self.assertEqual(before, pickle.loads(pickle.dumps(abstained[0])).__getstate__())

    def test_typed_boundary_and_occurrence_roundtrip_and_lru_eviction(self):
        original = family()[0]
        cache = LatticeDecodeCache(capacity=1)
        lattices = (original.boundary_path_lattice, occurrence_lattice())
        self.assertTrue(all(lattice is not None for lattice in lattices))
        references = [cache.compact(lattice) for lattice in lattices]
        first = references[0].resolve()
        for lattice, reference in zip(lattices, references):
            recovered = reference.resolve()
            self.assertIs(type(lattice), type(recovered))
            self.assertEqual(lattice, recovered)
            self.assertEqual(_canonical_sha256(lattice), _canonical_sha256(recovered))
            self.assertEqual(1, len(cache._decoded))
        again = references[0].resolve()
        self.assertEqual(first, again)
        self.assertIsNot(first, again)
        self.assertLess(len(references[0].payload), len(pickle.dumps(lattices[0])))

    def test_records_and_cache_die_with_last_reference(self):
        candidates = family()
        lattice = weakref.ref(candidates[0].boundary_path_lattice)
        cache = compact(candidates)
        gc.collect()
        self.assertIsNone(lattice(), 'batch memo or validation seal retained original lattice')
        reference = candidates[0].__dict__['boundary_path_lattice']
        decoded = weakref.ref(reference.resolve())
        token = weakref.ref(reference)
        self.assertTrue(cache._decoded)
        del candidates, reference
        gc.collect()
        self.assertIsNone(token())
        self.assertIsNone(decoded())
        self.assertEqual(0, len(cache._decoded))

    def test_corruption_and_invalid_records_fail_closed_even_after_cache_hit(self):
        original = family()[0].boundary_path_lattice
        cache = LatticeDecodeCache()
        reference = cache.compact(original)
        reference.resolve()
        object.__setattr__(reference, 'payload', reference.payload + b'x')
        with self.assertRaisesRegex(ValueError, 'digest'):
            reference.resolve()
        for serialized in (pickle.dumps(original) + b'extra', pickle.dumps({'unexpected': []}),
                           pickle.dumps(Exception('unsupported global'))):
            payload = zlib.compress(serialized)
            reference = CompressedLattice(payload, hashlib.sha256(payload).hexdigest(), type(original), cache)
            with self.assertRaises((ValueError, pickle.UnpicklingError)):
                reference.resolve()
        for payload in (zlib.compress(pickle.dumps(original)) + b'trailing',
                        zlib.compress(pickle.dumps(original))[:-2]):
            reference = CompressedLattice(payload, hashlib.sha256(payload).hexdigest(), type(original), cache)
            with self.assertRaisesRegex(ValueError, 'framing'):
                reference.resolve()
        with self.assertRaises(TypeError):
            cache.compact(dataclasses.replace(original, nodes=list(original.nodes)))
        with self.assertRaises(TypeError):
            cache.compact({'not': 'a lattice'})
        for capacity in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                LatticeDecodeCache(capacity)

    def test_every_candidate_serializes_legacy_logical_state_and_replays_after_eviction(self):
        candidates = family()
        expected = [candidate.__getstate__() for candidate in candidates]
        cache = compact(candidates, capacity=1)
        self.assertGreater(len(candidates), 2)
        shared = candidates[0].__dict__['boundary_path_lattice']
        self.assertTrue(all(candidate.__dict__['boundary_path_lattice'] is shared for candidate in candidates))
        for before, candidate in zip(expected, candidates):
            candidate._validate_source_path_provenance()
            self.assertEqual(before, candidate.__getstate__())
            for name in ('_source_provenance_validation_seal', '_boundary_path_validation_seal',
                         '_occurrence_path_validation_seal'):
                self.assertIsNone(getattr(candidate, name))
            loaded = pickle.loads(pickle.dumps(candidate))
            self.assertEqual(before, loaded.__getstate__())
            self.assertNotIsInstance(loaded.__dict__['boundary_path_lattice'], CompressedLattice)
        self.assertLessEqual(len(cache._decoded), 1)

    def test_geometry_counts_node_identity_and_record_replacement_still_rejected(self):
        mutations = [lambda c: c.list_absolute_start.__setitem__(0, c.list_absolute_start[0] + 1),
                     lambda c: c.list_repeats.__setitem__(0, 'T' + c.list_repeats[0][1:]),
                     lambda c: c.list_fuzzy_counts.__setitem__(0, (3, 3, 3)),
                     lambda c: setattr(c, 'boundary_path_node_ids', ('missing',)),
                     lambda c: setattr(c, 'boundary_path_lattice',
                                       dataclasses.replace(c.boundary_path_lattice, paths=()))]
        for mutation in mutations:
            candidate = family()[0]
            cache = compact([candidate], capacity=1)
            candidate._validate_source_path_provenance()
            # Evict the boundary record before corrupting authoritative fields.
            other = cache.compact(occurrence_lattice())
            other.resolve()
            mutation(candidate)
            with self.assertRaises(ValueError):
                candidate._validate_source_path_provenance()
            self.assertLessEqual(len(cache._decoded), 1)

    def test_local_chain_shares_only_immutable_storage_and_keeps_independent_geometry(self):
        parent = family()[0]
        expected = parent._slice_local_chain(0, 2).__getstate__()
        compact([parent])
        child = parent._slice_local_chain(0, 2)
        self.assertEqual(expected, child.__getstate__())
        for name in ('boundary_path_lattice', 'occurrence_path_lattice'):
            self.assertIs(parent.__dict__[name], child.__dict__[name])
        start = parent.list_absolute_start[0]
        child.list_absolute_start[0] += 1
        self.assertEqual(start, parent.list_absolute_start[0])
        with self.assertRaises(ValueError):
            child._validate_source_path_provenance()
        parent._validate_source_path_provenance()

    def test_same_start_end_and_refinement_keep_full_scientific_state(self):
        plain, stored = family(), family()
        compact(stored)
        for candidates in (plain, stored):
            # Exercise provenance merging, not just distinct candidate keys.
            candidates.append(pickle.loads(pickle.dumps(candidates[0])))
        key = ClusterSequence(SOURCE, 0, len(SOURCE), (REPEAT,))
        before = SameStartEndFilter({key: plain}).output()
        after = SameStartEndFilter({key: stored}).output()
        self.assertEqual([c.__getstate__() for c in before[key]], [c.__getstate__() for c in after[key]])
        results = [DetectionRefinement(value, PARAMETERS, False, source_dna=SOURCE) for value in (before, after)]
        def scientific_state(result):
            return [(tuple(key), [candidate.__dict__ for candidate in candidates])
                    for key, candidates in result.output().items()]
        self.assertEqual(scientific_state(results[0]), scientific_state(results[1]))
        self.assertEqual(results[0].boundary_repair_records, results[1].boundary_repair_records)
        self.assertEqual(results[0].internal_repair_records, results[1].internal_repair_records)

    def test_failed_compaction_leaves_both_fields_and_seals_unchanged(self):
        candidate = family()[0]
        candidate.occurrence_path_lattice = {'invalid': 'mutable record'}
        before = candidate.__dict__.copy()
        with self.assertRaises(TypeError):
            candidate._compact_lattice_storage(LatticeDecodeCache(), {})
        for key, value in before.items():
            self.assertIs(value, candidate.__dict__[key])

    def test_native_collector_matches_disabled_compression_serial_and_spawn(self):
        results = []
        for parallel in (False, True):
            for enabled in (False, True):
                detector = Detection.__new__(Detection)
                detector.cluster_sequences = [ClusterSequence(SOURCE, 0, len(SOURCE), (REPEAT, REPEAT[:-1]))]
                detector.dict_fuzzy_crisprs = {}
                detector.dna = SOURCE
                detector.flag_parallel = parallel
                detector.flag_cpu = '2'
                _initialize_fuzzy_source(SOURCE, True)
                method = FuzzySearch._compact_lattice_storage if enabled else lambda *_: None
                with mock.patch.object(FuzzySearch, '_compact_lattice_storage', method):
                    detector._run_native_array_detection(ERROR)
                candidates = next(iter(detector.dict_fuzzy_crisprs.values()))
                self.assertGreater(len(candidates), 2)
                self.assertEqual(enabled, any(c._has_compressed_lattices() for c in candidates))
                results.append([candidate.__getstate__() for candidate in candidates])
        for result in results[1:]:
            self.assertEqual(results[0], result)


if __name__ == '__main__':
    unittest.main()
