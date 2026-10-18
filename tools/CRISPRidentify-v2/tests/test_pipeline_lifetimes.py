"""Obsolete owners must die without erasing shared scientific evidence."""
from contextlib import ExitStack
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest import mock
import weakref

from components.pipeline import Pipeline


class Payload:
    pass


class PipelineLifetimeTests(unittest.TestCase):
    def test_rescue_releases_refinement_before_scoring_and_keeps_transferred_evidence(self):
        pipeline = Pipeline.__new__(Pipeline)
        pipeline.flags = {'flag_fast_run': True, 'flag_fast_run_rescue': True}
        pipeline.parameters = {}
        pipeline.flag_dev_mode = False
        pipeline.genome_sequence = 'A' * 100
        candidate, pending, internal, contextual = Payload(), Payload(), Payload(), Payload()
        rows = {(10, 30): [[.9, candidate, []]]}
        pipeline._initial_scored_rows = rows
        pipeline.categories = [rows, {}, {}, {}, {}]
        raw = {(10, 30): [Payload()]}
        pipeline._raw_detection = SimpleNamespace(
            cluster_indices_for_intervals=lambda intervals: ([0], {}),
            run_start_end_delta=lambda *args, **kwargs: (raw, {}),
        )
        refs = {}
        test = self

        class Refinement:
            def __init__(self, **kwargs):
                test.assertIs(kwargs['consume_input'], True)
                kwargs['dict_fuzzy_crisprs'].clear()
                self.workspace = Payload()
                refs['workspace'] = weakref.ref(self.workspace)
                refs['owner'] = weakref.ref(self)
                self.pending_internal_repairs = [pending]
                self.pending_contextual_spacers = [contextual]
                self.internal_repair_records = [internal]
                self.boundary_repair_records = [{'proof': 'preserved'}]
                self.boundary_repair_stats = {'proposed': 1}

            def output(self):
                return {(10, 30): [candidate]}

        def score(candidates, previous, scorer):
            self.assertIsNone(refs['workspace']())
            self.assertIsNone(refs['owner']())
            self.assertIs(candidates[(10, 30)][0], candidate)
            self.assertIs(pipeline.pending_internal_repairs[0], pending)
            self.assertIs(pipeline.pending_contextual_spacers[0], contextual)
            self.assertIs(pipeline.internal_repair_records[0], internal)
            self.assertEqual(pipeline.repeat_boundary_repair_records,
                             [{'proof': 'preserved', 'stage': 'fast_rescue'}])
            return rows, {'recomputed': 0, 'reused': 1}

        with ExitStack() as stack:
            def patch(name, **kwargs):
                return stack.enter_context(mock.patch('components.pipeline.' + name, **kwargs))
            patch('DetectionRefinement', new=Refinement)
            patch('build_advanced_fuzzy_search_filter')
            patch('select_unresolved_components', return_value=([
                {'source_interval': (10, 30), 'trigger_score': .9, 'features': {}}
            ], {}))
            patch('collapse_source_backed_phase_equivalent_candidates',
                  side_effect=lambda source, candidates, **kwargs: (candidates, (), {}))
            patch('candidate_score_fingerprint', return_value='same-candidate')
            patch('ArrayEvaluation.categorize_scored_rows', return_value=pipeline.categories)
            patch('CandidateScoreCache')
            for name in ('_record_source_phase_collapse', '_validate_refined_candidates'):
                stack.enter_context(mock.patch.object(pipeline, name))
            stack.enter_context(mock.patch.object(pipeline, '_accepted_structural_signatures',
                                                  return_value=Counter()))
            stack.enter_context(mock.patch.object(pipeline, '_accepted_intervals_by_identity',
                                                  return_value={}))
            stack.enter_context(mock.patch.object(pipeline, '_score_rows_with_structural_cache',
                                                  side_effect=score))
            pipeline._run_fast_routed_start_end_rescue()
        self.assertEqual({}, pipeline.dict_fuzzy_crisprs)
        self.assertIs(pipeline.dict_crispr_candidates[(10, 30)][0], candidate)

    def test_context_capture_is_separate_and_only_enabled_by_conservative_policy(self):
        for policy in ("legacy", "conservative"):
            pipeline = Pipeline.__new__(Pipeline)
            pipeline.parameters = {"param_reconstruction_policy": policy}
            pipeline.genome_sequence = "ACGT" * 25
            pipeline.dict_fuzzy_crisprs = {}
            pipeline.flag_dev_mode = False
            original, captured = Payload(), Payload()
            refinement = SimpleNamespace(
                output=lambda: {(0, 100): [original]},
                pending_contextual_spacers=[captured] if policy == "conservative" else [],
                boundary_repair_records=[], boundary_repair_stats={},
            )
            with self.subTest(policy=policy), mock.patch(
                    "components.pipeline.DetectionRefinement", return_value=refinement) as detect:
                pipeline._run_detection_refinement()
            self.assertEqual(policy == "conservative", detect.call_args.kwargs[
                "capture_contextual_spacers"])
            self.assertEqual({(0, 100): [original]}, pipeline.dict_crispr_candidates)
            self.assertEqual([captured] if policy == "conservative" else [],
                             pipeline.pending_contextual_spacers)

    def test_context_transaction_respects_completed_acceptance_and_releases_pending_owner(self):
        pipeline = Pipeline.__new__(Pipeline)
        pipeline.parameters = {"param_reconstruction_policy": "conservative"}
        pipeline.genome_sequence = "ACGT"
        pipeline.categories = [{}, {}, {}, {}, {}]
        pipeline._source_cds = None
        pending = Payload()
        pending_ref = weakref.ref(pending)
        pipeline.pending_contextual_spacers = [pending]
        del pending
        order = []
        result = [{}, {}, {(0, 4): [[.8, Payload(), []]]}, {}, {}]

        def ordinary(name):
            def apply(source, categories, parameters, scorer, source_cds):
                order.append(name)
                self.assertIsNotNone(pending_ref())
                return categories, {"stage": name}
            return apply

        def context(source, categories, captured, parameters, scorer, source_cds):
            order.append("context")
            self.assertEqual(["BF", "Possible", "context"], order)
            self.assertIs(captured[0], pending_ref())
            self.assertEqual({"stage": "Possible"}, pipeline.possible_reconstruction_evidence)
            return result, {"stage": "context"}

        with mock.patch("components.pipeline.reconstruct_accepted_arrays", side_effect=ordinary("BF")), \
                mock.patch("components.pipeline.reconstruct_possible_arrays", side_effect=ordinary("Possible")), \
                mock.patch("components.pipeline.rescue_contextual_spacer_arrays", side_effect=context) as rescue:
            pipeline._reconstruct_accepted_arrays()
            # Mock call history deliberately owns the argument until released.
            rescue.reset_mock()
        self.assertIsNone(pending_ref())
        self.assertEqual([], pipeline.pending_contextual_spacers)
        self.assertIs(result, pipeline.categories)
        self.assertEqual({"stage": "context"}, pipeline.contextual_spacer_evidence)

    def test_constructor_releases_workspace_only_after_its_last_consumer(self):
        refs = {}
        retained = Payload()
        categories = [{(10, 30): [[.9, retained, []]]}, {}, {}, {}, {}]
        test = self

        def payload(name):
            value = Payload()
            refs[name] = weakref.ref(value)
            return value

        def detect(pipeline):
            raw_map = {(10, 30): [payload('raw_witness')]}
            pipeline._raw_detection = SimpleNamespace(dict_fuzzy_crisprs=raw_map)
            pipeline.dict_fuzzy_crisprs = raw_map

        def refine(pipeline):
            pipeline.dict_crispr_candidates = {(10, 30): [payload('refinement_workspace')]}

        def evaluate(pipeline):
            pipeline.categories = categories
            pipeline._initial_scored_rows = {(10, 30): [payload('initial_workspace')]}
            pipeline._candidate_score_cache = payload('score_cache')

        def rescue(pipeline):
            test.assertIsNotNone(refs['refinement_workspace']())
            test.assertIsNotNone(refs['initial_workspace']())
            test.assertIsNotNone(refs['raw_witness']())

        def enhance(pipeline):
            test.assertIsNone(refs['refinement_workspace']())
            test.assertIsNone(refs['initial_workspace']())
            # Rebinding the pipeline map must not clear this shared raw map.
            test.assertIs(pipeline._raw_detection.dict_fuzzy_crisprs[(10, 30)][0], refs['raw_witness']())
            test.assertIsNotNone(refs['score_cache']())

        def reevaluate(pipeline, **kwargs):
            test.assertIsNotNone(refs['score_cache']())
            pipeline._post_enhancement_scored_rows = {(10, 30): [payload('late_workspace')]}
            pipeline._final_scored_rows = pipeline._post_enhancement_scored_rows

        def two_repeat(pipeline):
            test.assertIsNone(refs['score_cache']())
            test.assertIsNotNone(refs['raw_witness']())

        def family(pipeline):
            test.assertIsNotNone(refs['raw_witness']())
            test.assertIs(pipeline._fast_accepted_family_anchor_categories[0][(10, 30)][0][1], retained)

        def phase(pipeline):
            test.assertIsNone(refs['raw_witness']())
            test.assertIsNone(pipeline._raw_detection)
            test.assertIsNone(pipeline._fast_accepted_family_anchor_categories)
            pipeline._candidate_score_cache = payload('late_cache')

        def augment(pipeline):
            test.assertIsNotNone(refs['late_workspace']())
            return False

        def reconstruct(pipeline):
            test.assertIsNotNone(refs['late_workspace']())
            pipeline.reconstruction_evidence = {'preserved': True}

        def quarantine(pipeline):
            test.assertIsNone(refs['late_workspace']())
            test.assertIsNone(refs['late_cache']())
            test.assertIs(pipeline.categories[0][(10, 30)][0][1], retained)
            test.assertEqual({'preserved': True}, pipeline.reconstruction_evidence)

        effects = {
            '_run_detection': detect, '_run_detection_refinement': refine,
            '_run_evaluation': evaluate, '_run_fast_routed_start_end_rescue': rescue,
            '_results_enhancement': enhance, '_reevaluate_enhanced_candidates': reevaluate,
            '_run_fast_two_repeat_terminal_repair': two_repeat,
            '_run_fast_accepted_family_rescue': family, '_run_fast_bounded_phase_repair': phase,
            '_run_configured_source_backed_augmentation': augment,
            '_reconstruct_accepted_arrays': reconstruct,
            '_quarantine_likely_tandems_transactionally': quarantine,
        }
        unused = (
            '_get_header', '_validate_refined_candidates',
            '_collapse_source_backed_phase_equivalent_candidates',
            '_validate_candidate_categories', '_normalize_post_enhancement_source_phase',
            '_complete_attested_terminal_repairs', '_complete_pending_internal_repeats',
            '_annotate_model_score_classification_bases',
            '_demote_accepted_family_tandems_transactionally',
            '_complete_accepted_core_boundaries_transactionally',
            '_recover_source_intact_candidates_transactionally',
            '_repair_source_periodic_lattices_transactionally',
            '_recover_source_regular_subpaths_transactionally',
            '_evaluate_circular_origin_joins', '_run_non_crispr_computation',
            '_validate_and_write_output',
        )
        with TemporaryDirectory() as directory, ExitStack() as stack:
            for name in unused:
                stack.enter_context(mock.patch.object(Pipeline, name, autospec=True))
            for name, effect in effects.items():
                stack.enter_context(mock.patch.object(Pipeline, name, autospec=True, side_effect=effect))
            pipeline = Pipeline(directory, '', '', str(Path(directory)/'source.fa'),
                                {'param_accepted_only': False},
                                {'flag_fast_accepted_family_rescue': True}, False,
                                retain_candidate_evidence=True)
        self.assertIs(pipeline._candidate_evidence_categories[0][(10, 30)][0][1], retained)
        self.assertIs(pipeline.categories[0][(10, 30)][0][1], retained)


if __name__ == '__main__':
    unittest.main()
