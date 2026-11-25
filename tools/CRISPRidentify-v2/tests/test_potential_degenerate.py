import random
import tempfile
import unittest
from pathlib import Path

from components.refinement.enhancement_search import build_source_backed_candidate_from_boundaries
from components.refinement.potential_degenerate import annotate_members, complete, geometry, propose
from components.reporting.html_report import _component_table
from components.reporting.report_model import ReportModelBuilder


class PotentialDegenerateTests(unittest.TestCase):
    motif = "GCTTCAGCTGTTAGGATCGTCAACGCTCAG"
    # Seven separated substitutions beyond the default six-edit budget.
    mutant = "".join({"A":"C", "C":"G", "G":"T", "T":"A"}[b] if i in (3,6,10,13,17,21,25)
                     else b for i, b in enumerate("GCTTCAGCTGTTAGGATCGTCAACGCTCAG"))

    def fixture(self, missing, strand=False, sequence=None):
        copies = [self.motif] * 4
        copies[missing] = sequence if sequence is not None else self.mutant
        rng = random.Random(42)
        dna = ""
        starts = []
        for repeat in copies:
            if dna:
                dna += "".join(rng.choices("ACGT", k=31))
            starts.append(len(dna))
            dna += repeat
        intervals = [(s, s+len(r)) for s, r in zip(starts, copies)]
        if strand:
            dna = dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            intervals = [(len(dna)-b, len(dna)-a) for a, b in intervals]
        target = intervals.pop(missing)
        core = build_source_backed_candidate_from_boundaries(dna, None, *zip(*sorted(intervals)))
        return dna, core, target

    def test_ends_internal_and_reverse_with_indels(self):
        for reverse in (False, True):
            for missing in (0, 1, 3):
                for seq in (self.mutant, self.mutant[:15] + "A" + self.mutant[15:]):
                    with self.subTest(reverse=reverse, missing=missing, seq=seq):
                        dna, core, loc = self.fixture(missing, reverse, seq)
                        found = propose(dna, core, (0, len(dna)), 6)
                        self.assertEqual({loc}, set(found))
                        ev = found[loc]
                        self.assertEqual("potential_degenerate_repeat", ev["status"])
                        self.assertEqual(ev["edit_distance"], sum(ev[k] for k in
                                         ("substitutions", "insertions", "deletions")))
                        self.assertLessEqual(ev["edit_distance"], ev["relaxed_edit_limit"])

    def test_negative_controls(self):
        for sequence in (self.motif, "A"*30, self.mutant[:4]+"N"+self.mutant[5:]):
            dna, core, _ = self.fixture(3, sequence=sequence)
            self.assertFalse(propose(dna, core, (0, len(dna)), 6))
        dna, core, loc = self.fixture(3)
        self.assertFalse(propose(dna, core, (0, loc[0]), 6))

    def test_scoring_rollback_and_public_membership(self):
        dna, core, loc = self.fixture(0)
        key = (0, len(dna))
        row = [.9, core, []]
        categories = [{key:[row]}, {}, {}, {}, {}]
        def score(requests):
            return {key: [[.8, x, []] for x in requests[key]]}
        scored = {key:[row]}
        unchanged, _, accepted = complete(dna, categories, scored, 6, score, lambda _: categories)
        self.assertFalse(accepted)
        self.assertIs(unchanged[0][key][0][1], core)
        changed, _, committed = complete(dna, categories, scored, 6, score,
                                         lambda rows: [rows, {}, {}, {}, {}])
        annotate_members(changed, committed)
        winner = changed[0][key][0][1]
        self.assertEqual(set(geometry(core)) | {loc}, set(geometry(winner)))
        with tempfile.TemporaryDirectory() as wd:
            fasta = Path(wd)/"test.fa"
            fasta.write_text(">record\n"+dna+"\n")
            model = ReportModelBuilder(fasta, changed, {}, []).build()
        array, = model["arrays"]
        self.assertEqual(4, array["repeat_count"])
        self.assertEqual(3, array["spacer_count"])
        self.assertEqual(dna, array["source_sequence"])
        first = array["repeats"][0]
        self.assertEqual(self.mutant, first["sequence"])
        self.assertEqual("potential_degenerate_repeat", first["status"])
        self.assertEqual(dict(type="repeat", **first), array["architecture"][0])
        self.assertIn("Potential degenerate repeat", _component_table(array))
