"""Owned rescue storage may shrink only after complete evidence transfer."""
from contextlib import ExitStack
import gc
import pickle
import unittest
import weakref
from unittest import mock

from components.components_detection import ClusterSequence, FuzzySearch, RejectedFuzzyCandidateTopology
from components.module_detection_refinement import DetectionRefinement
from components import module_detection_refinement as module
from tests.test_fuzzy_lattice_storage import family, compact, PARAMETERS, REPEAT, SOURCE


def clone(value):
    return pickle.loads(pickle.dumps(value))


def scientific_state(result):
    return {
        'output': [(tuple(key), [candidate.__dict__ for candidate in candidates])
                   for key, candidates in result.output().items()],
        'fallback': [(key.start, key.end) for key in result.fallback_keys],
        'empty_witness': [(key.start, key.end) for key in result.compact_witness_empty_keys],
        'boundary_records': result.boundary_repair_records,
        'boundary_stats': result.boundary_repair_stats,
        'internal_records': result.internal_repair_records,
        'parents': [candidate.__dict__ for candidate in result.internal_parent_candidates],
        'pending': [candidate.__dict__ for candidate in result.pending_internal_repairs],
    }


class RefinementConsumptionTests(unittest.TestCase):
    def compare(self, mapping, parameters=PARAMETERS, source=SOURCE):
        observations, results = [], []
        original_proposals = module.propose_internal_repairs
        stages = ('_filter_out_same_start_end_cases', '_split_disjoint_source_loci',
                  '_filter_out_non_crispr_cases', '_collect_pending_internal_repairs',
                  '_reformat_ac_crispr_candidates')
        originals = {name:getattr(DetectionRefinement,name) for name in stages}
        for consume in (False, True):
            current = clone(mapping)
            compact([candidate for rows in current.values() for candidate in rows
                     if type(candidate) is FuzzySearch])
            before_counts = [len(values) for values in current.values()]
            calls, order = [], []
            def proposals(context, candidates, predicate):
                calls.append([candidate.__dict__.copy() for candidate in candidates])
                return original_proposals(context, candidates, predicate)
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(module, 'propose_internal_repairs', side_effect=proposals))
                for name in stages:
                    def wrapped(instance, name=name):
                        order.append(name)
                        return originals[name](instance)
                    stack.enter_context(mock.patch.object(DetectionRefinement,name,wrapped))
                result = DetectionRefinement(current, parameters, False, source_dna=source,
                                             consume_input=consume)
            self.assertEqual(list(stages), order)
            self.assertEqual({} if consume else before_counts,
                             current if consume else [len(values) for values in current.values()])
            observations.append(calls)
            results.append(scientific_state(result))
        self.assertEqual(observations[0], observations[1])
        self.assertEqual(results[0], results[1])
        return results[1]

    def test_complete_refinement_parity_with_all_global_parents_and_duplicate_origins(self):
        first = ClusterSequence(SOURCE, 0, len(SOURCE), (REPEAT,))
        second = ClusterSequence(SOURCE, 1, len(SOURCE), (REPEAT,))
        rows = family()
        compact(rows)
        mapping = {first: rows, second: [rows[0], rows[-1]]}
        result = self.compare(mapping)
        self.assertGreater(len(result['parents']), 100)
        self.assertTrue(result['output'])

    def test_actual_terminal_and_pending_internal_repairs_keep_complete_records(self):
        from tests.test_repeat_boundary_repair import synthetic as terminal_fixture, fuzzy, Locus
        from tests.test_internal_repeat_completion import synthetic as internal_fixture
        source, parent = terminal_fixture()
        terminal = self.compare({Locus(0, len(source)): [fuzzy(parent)]}, source=source)
        self.assertGreater(terminal['boundary_stats']['children_generated'], 0)
        source, parent, _ = internal_fixture(missing=(3, 4))
        internal = self.compare({Locus(0, len(source)): [fuzzy(parent)]}, source=source)
        self.assertEqual(1, len(internal['pending']))

    def test_no_pass_fallback_and_topology_connectivity_remain_equal(self):
        key = ClusterSequence(SOURCE, 0, len(SOURCE), (REPEAT,))
        rows = family()
        rows.append(RejectedFuzzyCandidateTopology(((0, len(SOURCE)),), (), 1, 1))
        result = self.compare({key: rows}, dict(PARAMETERS, param_min_repeats=100))
        self.assertTrue(result['fallback'])
        self.assertEqual(1, sum(len(rows) for _, rows in result['output']))

    def test_delta_raw_rows_die_incrementally_while_original_detector_rows_stay_live(self):
        key = ClusterSequence(SOURCE, 0, len(SOURCE), (REPEAT,))
        rows = family()
        compact(rows)
        original_map = {key: [rows[0]]}
        merged_map = {key: list(rows)}
        refs = [weakref.ref(candidate) for candidate in rows]
        original_geometry = clone(rows[0].__getstate__())
        del rows
        seen = []
        normalize = DetectionRefinement._normalize_fuzzy_candidate
        proposals = module.propose_internal_repairs
        def checked_normalize(candidate):
            if seen:
                gc.collect()
                # The first row belongs to the original detector and must
                # survive; later delta-only rows must die before this call.
                self.assertTrue(all(ref() is None for ref in seen[1:]))
            seen.append(weakref.ref(candidate))
            return normalize(candidate)
        def checked_proposals(context, candidates, predicate):
            self.assertEqual(len(refs), len(candidates))
            self.assertTrue(all(ref() is None for ref in refs[1:]))
            self.assertIs(refs[0](), original_map[key][0])
            return proposals(context,candidates,predicate)
        with mock.patch.object(DetectionRefinement, '_normalize_fuzzy_candidate', new=staticmethod(checked_normalize)), \
                mock.patch.object(module, 'propose_internal_repairs', side_effect=checked_proposals):
            result = DetectionRefinement(merged_map, PARAMETERS, False, source_dna=SOURCE, consume_input=True)
        self.assertEqual({}, merged_map)
        self.assertEqual({}, result.dict_fuzzy_crisprs_refined_st_end)
        self.assertEqual(original_geometry, original_map[key][0].__getstate__())
        self.assertEqual(len(refs), len(result.internal_parent_candidates))

    def test_default_preserves_raw_inventory_and_boolean_flag_is_explicit(self):
        key = ClusterSequence(SOURCE, 0, len(SOURCE), (REPEAT,))
        rows = family()
        original = {key: rows}
        result = DetectionRefinement(original, PARAMETERS, False, source_dna=SOURCE)
        self.assertIs(rows, original[key])
        self.assertIs(original, result.dict_fuzzy_crisprs)
        self.assertTrue(all(type(candidate) is FuzzySearch for candidate in rows))
        for flag in ('True', 1, None):
            with self.assertRaises(TypeError):
                DetectionRefinement(original, PARAMETERS, False, source_dna=SOURCE, consume_input=flag)


if __name__ == '__main__':
    unittest.main()
