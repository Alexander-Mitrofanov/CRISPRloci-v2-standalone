"""One chosen ordinary observation moves; unused alignment evidence survives."""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from components.candidate_evidence import stage_candidate_evidence
from components.candidate_origins import candidate_origin_evidence, inherit_candidate_origins
from components.compact_model import predict_compact_scores
from components.components_detection_refinement import CrisprCandidate
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.output_transaction import OutputTransaction
from components.reconstruction_transaction import reconstruct_accepted_arrays
from components.report_model import ReportModelBuilder
from tests.test_candidate_evidence import PARAMETERS
from tests.test_reconstruction_boundary_family import fixture as family_fixture
from tests.test_reconstruction_transaction import fixture as ordinary_fixture


FIXTURE = Path(__file__).with_name('fixtures') / 'reconstruction_alias_cp108133.json'
FIXTURE_SHA = '10ed84b5a4f934a474498b7da1ce5f13ecd5530dfcc6d4a219a76254977ba7c5'


def cp108_fixture():
    """Exact observed local strings and vectors; no new feature extraction."""
    raw = FIXTURE.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == FIXTURE_SHA
    data = json.loads(raw)
    source = data['source_sequence']
    assert hashlib.sha256(source.encode()).hexdigest() == data['source_sha256']
    rows = []
    for saved in data['candidates']:
        iv = saved['repeat_intervals']
        candidate = CrisprCandidate(list_repeat_starts=[a for a, _ in iv],
            list_repeats=[source[a:b] for a, b in iv], list_repeats_gaped=saved['gapped_repeats'],
            list_spacers=[source[x[1]:y[0]] for x, y in zip(iv, iv[1:])])
        # Final evidence has compact derivation counts, not the recursive raw
        # history. These test records exercise attribute lifetime only.
        candidate.source_path_provenance = tuple(
            SimpleNamespace(derivation=kind, fixture_record={'ordinal': i, 'name': saved['name']})
            for kind, count in (saved['source_path_provenance_summary'] or {}).get('derivations', {}).items()
            for i in range(count))
        features = np.array([[saved['features'][n] for n in COMPACT_FEATURE_SCHEMA]])
        assert float(predict_compact_scores(features)[0]) == saved['score']
        rows.append([saved['score'], candidate, [features]])
    inherit_candidate_origins(rows[2][1], [rows[1][1]])
    features_by_candidate = {id(r[1]): r[2][0] for r in rows}

    def scorer(requests):
        return {key: [[float(predict_compact_scores(features_by_candidate[id(c)])[0]),
                       c, [features_by_candidate[id(c)].copy()]] for c in values]
                for key, values in requests.items()}

    key = (0, len(source))
    return source, [{key: [rows[0]]}, {key: rows[1:]}, {}, {}, {}], scorer


class OrdinaryAliasRetentionTests(unittest.TestCase):
    def assert_row_unchanged(self, row, before):
        self.assertEqual(row[0], before[0])
        self.assertEqual(row[1].__dict__, before[1].__dict__)
        self.assertEqual(len(row), len(before))
        self.assertEqual(len(row[2]), len(before[2]))
        for a, b in zip(row[2], before[2]):
            np.testing.assert_array_equal(a, b)
        self.assertEqual(row[3:], before[3:])

    def test_real_cp108_alignments_origins_and_provenance_survive_permutations(self):
        expected = None
        for reverse in (False, True):
            source, categories, scorer = cp108_fixture()
            key = next(iter(categories[0]))
            if reverse:
                categories[1][key].reverse()
            originals = categories[1][key]
            core = categories[0][key][0]
            snapshots = [copy.deepcopy(r) for r in [core]+originals]
            output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
            selected = output[0][key][0]
            decision = audit['decisions'][0]
            chosen_index = decision['selected']['row_index']
            unused = originals[1-chosen_index]
            self.assertEqual(audit['actions'], {'selected': 1})
            self.assertEqual(decision['evidence']['same_rank_peer_geometries'], 1)
            self.assertEqual(decision['evidence']['selected_alignment_tie_count'], 2)
            self.assertEqual(len(decision['selected_aliases']), 2)
            self.assertEqual(len(output[0][key]), 1)
            self.assertEqual(len(output[1][key]), 2)
            self.assertIs(output[1][key][0], unused)
            self.assertIs(output[1][key][1], core)
            self.assertTrue(all(r is not originals[chosen_index] for r in output[1][key]))
            self.assertEqual(selected[0], .7420351112150484)
            self.assertEqual(len(selected[1].list_repeats), 5)
            signature = (selected[0], selected[1].list_gaped_repeats,
                         selected[3]['classification_basis']['reconstruction']['selection']['relationship'])
            if expected is None:
                expected = copy.deepcopy(signature)
            else:
                self.assertEqual(signature, expected)
            for row, before in zip([core]+originals, snapshots):
                self.assert_row_unchanged(row, before)
            origins = [candidate_origin_evidence(r[1]) for r in originals]
            self.assertNotEqual(origins[0]['current_representation_sha256'], origins[1]['current_representation_sha256'])
            self.assertEqual(origins[0]['ordinary_birth_sha256'], origins[1]['ordinary_birth_sha256'])
            with tempfile.TemporaryDirectory() as temporary:
                fasta = Path(temporary)/'cp108.fa'
                fasta.write_text('>cp108\n'+source+'\n')
                parameters = dict(PARAMETERS, param_reconstruction_policy='conservative')
                report = ReportModelBuilder(fasta, output, {}, [COMPACT_FEATURE_SCHEMA], parameters=parameters).build()
                path = Path(temporary)/'candidate_evidence.jsonl.gz'
                with OutputTransaction() as publication:
                    stage_candidate_evidence(publication, path, source_header='cp108', source_sequence=source,
                        categories=output, parameters=parameters, feature_names=[COMPACT_FEATURE_SCHEMA],
                        canonical_report=report)
                    publication.publish()
                records = [json.loads(line) for line in gzip.decompress(path.read_bytes()).splitlines()]
            candidates = [r for r in records if r['record_type']=='candidate']
            self.assertEqual(len(candidates), 3)
            self.assertEqual(sum(r['accepted'] for r in candidates), 1)
            aliases = [r for r in candidates if r['repeat_count']==5]
            self.assertEqual(len(aliases), 2)
            self.assertEqual(len({r['source_geometry_id'] for r in aliases}), 1)
            self.assertEqual(len({r['alignment_id'] for r in aliases}), 2)
            self.assertEqual(sorted(r['source_path_provenance']['record_count'] for r in aliases), [0, 1])
            retained = next(r for r in aliases if not r['accepted'])
            self.assertEqual(retained['category'], 'Alternative')
            self.assertEqual(retained['gapped_repeats'], unused[1].list_gaped_repeats)
            self.assertEqual(retained['candidate_origins'], candidate_origin_evidence(unused[1]))

    def test_identical_row_object_in_two_positions_moves_only_one_observation(self):
        source, categories, scorer = ordinary_fixture()
        key = next(iter(categories[1]))
        donor = categories[1][key][0]
        categories[1][key].append(donor)
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit['stats']['selected'], 1)
        self.assertEqual(len(audit['decisions'][0]['selected_aliases']), 2)
        self.assertEqual(len(output[0][key]), 1)
        self.assertEqual(len(output[1][key]), 2)
        self.assertIs(output[1][key][0], donor)
        self.assertIs(output[1][key][1], categories[0][key][0])
        self.assertEqual(len(categories[1][key]), 2)

    def test_distinct_identical_observations_and_lower_score_alias_stay_unchanged(self):
        source, categories, original_scorer = ordinary_fixture()
        key = next(iter(categories[1]))
        donor = categories[1][key][0]
        duplicate, lower = copy.deepcopy(donor), copy.deepcopy(donor)
        lower[2][0][0, 2] = .7
        lower[0] = float(predict_compact_scores(lower[2][0])[0])
        self.assertLess(lower[0], donor[0])
        categories[1][key].extend([duplicate, lower])
        snapshots = [copy.deepcopy(r) for r in categories[1][key]]

        def scorer(requests):
            result = original_scorer(requests)
            for rows in result.values():
                for row in rows:
                    if row[1] is lower[1]:
                        row[0], row[2] = lower[0], copy.deepcopy(lower[2])
            return result

        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        decision = audit['decisions'][0]
        self.assertEqual(audit['stats']['selected'], 1)
        self.assertEqual(decision['reconstruction_score'], donor[0])
        self.assertEqual(decision['evidence']['selected_alignment_tie_count'], 2)
        chosen = decision['selected']['row_index']
        remaining = [r for i, r in enumerate(categories[1][key]) if i != chosen]
        self.assertEqual(len(output[1][key]), 3)
        self.assertTrue(all(a is b for a,b in zip(output[1][key], remaining+[categories[0][key][0]])))
        for row, before in zip(categories[1][key], snapshots):
            self.assert_row_unchanged(row, before)

    def test_family_still_retains_all_original_observations(self):
        source, categories, scorer = family_fixture()
        key = next(iter(categories[1]))
        categories[1][key].append(copy.deepcopy(categories[1][key][0]))
        originals = list(categories[1][key])
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit['stats']['boundary_uncertain_selected'], 1)
        self.assertEqual(len(output[0][key]), 1)
        self.assertEqual(len(output[1][key]), len(originals)+1)
        self.assertTrue(all(a is b for a,b in zip(output[1][key], originals+[categories[0][key][0]])))
        proof = output[0][key][0][3]['classification_basis']['boundary_family']
        self.assertEqual(len(proof['members']), 2)
        self.assertEqual(len(proof['additional_peer_representations']), 1)


if __name__ == '__main__':
    unittest.main()
