"""Core-only phase evidence, legacy replay and coupled-member transactions.

The real CP109078 fixture supplies source/geometry arithmetic, not new native
features or predicted acceptance. Transaction vectors below are deliberately
synthetic, scored by the bundled compact model without external predictors.
"""
import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest import mock

import numpy as np

from components.compact_model import predict_compact_scores
from components.components_detection_refinement import CrisprCandidate
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.reconstruction_boundary_family import METHOD as FAMILY_METHOD, digest
from components.reconstruction_reference import (
    ANCHOR_REFERENCE, MAPPED_CORE_REFERENCE, assess_added_repeat_reference,
    recorded_reference_rule,
)
from components.reconstruction_selection import (
    candidate_geometry, monotone_repeat_mapping, relationship, select_reconstructions,
)
from components.reconstruction_transaction import (
    reconstruct_accepted_arrays, validate_boundary_family_basis, validate_reconstruction_basis,
)
from components.report_model import ReportModelBuilder
from tests.test_candidate_evidence import PARAMETERS
from tests.test_reconstruction_boundary_family import qualified, fixture as legacy_family_fixture
from tests.test_reconstruction_transaction import fixture as legacy_fixture


CORE = "ACCACCCCCACGGGCGTGGGGAGCAC"
PHASE = "CGG" + CORE
WEAK_OLD = "CGGACCACCCCCACGCGTGAGGAGCAC"
FIXTURE = Path(__file__).with_name("fixtures") / "reconstruction_phase_cp109078.json"
FIXTURE_SHA = "a5a51cee95230ff4daa857f40e6aeb3c214856100f771cc93ba62ef10f1c1838"


def from_source(source, intervals):
    repeats = [source[a:b] for a, b in intervals]
    width = max(map(len, repeats))
    return CrisprCandidate(list_repeats=repeats,
        list_repeats_gaped=[r + "-" * (width-len(r)) for r in repeats],
        list_repeat_starts=[a for a, b in intervals],
        list_spacers=[source[a[1]:b[0]] for a, b in zip(intervals, intervals[1:])])


def phase_fixture(*, family=False, aliases=False):
    """Real-source transaction fixture with independently model-scored rows."""
    # The long and trimmed fourth copies have legacy identities 26/30 and
    # 22/28 respectively. Their shared start permits a mixed-reference family.
    fourth = "CGGACCACCCCCACGGGCGTGGGGGAGCAC" if family else WEAK_OLD
    repeats = [PHASE] * 3 + [fourth] + [PHASE] * 2
    source = "ACGT" * 8
    intervals = []
    rng = random.Random(5077)
    for index, repeat in enumerate(repeats):
        start = len(source)
        source += repeat
        intervals.append((start, len(source)))
        source += "".join(rng.choice("ACGT") for _ in range(31 + index))
    core = from_source(source, [(a+3,b) for a,b in intervals[:3]])
    donors = [from_source(source, intervals)]
    if family:
        trimmed = list(intervals)
        trimmed[3] = (trimmed[3][0], trimmed[3][1]-2)
        donors.append(from_source(source, trimmed))
        if aliases:
            donors.append(copy.deepcopy(donors[-1]))

    def score(candidate):
        n = len(candidate.list_repeats)
        # Synthetic feature rows intentionally make the two full geometries
        # equally ranked. This is binding/replay coverage, not extraction.
        features = np.array([[40, n, 1 if n == 3 else .97, .45, 34.5,
                              .254305, 0, 0, -8.2, 15.9]])
        return [float(predict_compact_scores(features)[0]), candidate, [features]]

    def score_batch(requests):
        return {key: [score(candidate) for candidate in values] for key, values in requests.items()}

    key = (0, len(source))
    return source, [{key: [score(core)]}, {key: [score(d) for d in donors]}, {}, {}, {}], score_batch


class SourcePhaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = FIXTURE.read_bytes()
        cls.saved = json.loads(cls.raw)
        cls.candidates = {}
        for item in cls.saved["candidates"]:
            cls.candidates[item["label"]] = CrisprCandidate(
                list_repeats=item["repeat_sequences"],
                list_repeats_gaped=item["gapped_repeats"],
                list_spacers=item["spacer_sequences"],
                list_repeat_starts=[a for a, b in item["repeat_intervals"]])

    def test_immutable_fixture_and_all_source_slices(self):
        self.assertEqual(hashlib.sha256(self.raw).hexdigest(), FIXTURE_SHA)
        source = self.saved["source"]
        text = source["slice_sequence"]
        offset, end = source["slice_interval"]
        self.assertEqual(len(text), end-offset)
        self.assertEqual(hashlib.sha256(text.encode()).hexdigest(), source["slice_sha256"])
        for item in self.saved["candidates"]:
            geometry = item["repeat_intervals"]
            self.assertEqual([text[a-offset:b-offset] for a,b in geometry], item["repeat_sequences"])
            self.assertEqual([text[a[1]-offset:b[0]-offset] for a,b in zip(geometry,geometry[1:])],
                             item["spacer_sequences"])

    def test_real_59_legacy_pass_65_66_fallback_and_67_weak_copy(self):
        core = self.candidates["original_core"]
        published = self.candidates["published59"]
        self.assertEqual(relationship(core,published),
                         relationship(core,published,reference_rule=ANCHOR_REFERENCE))
        self.assertNotIn("added_repeat_reference", relationship(core,published))
        for name, expected_min in (("retained65", 1-3/29), ("retained66", 1-4/30)):
            donor = self.candidates[name]
            self.assertFalse(relationship(core,donor,reference_rule=ANCHOR_REFERENCE)["qualified"])
            result = relationship(core,donor)
            self.assertTrue(result["qualified"])
            ref = result["added_repeat_reference"]
            self.assertEqual(ref["mapped_core_indexes"], [0,1,2])
            self.assertEqual(ref["boundary_offsets"], [-3,0])
            self.assertEqual(ref["original_consensus"], CORE)
            self.assertEqual(ref["mapped_core_sequence"], PHASE)
            self.assertEqual(result["minimum_added_repeat_identity"], expected_min)
            last = next(x for x in ref["added_repeat_checks"]
                        if candidate_geometry(donor)[x["index"]] == (5670899,5670926))
            self.assertEqual((last["original_edit_distance"],last["original_denominator"]), (6,27))
            self.assertEqual((last["mapped_edit_distance"],last["mapped_denominator"]), (3,29))
        donor = self.candidates["retained67"]
        result = self.assess(core, donor)
        self.assertFalse(result["qualified"])
        self.assertEqual(result["reference_status"], "mapped_reference_identity_below_limit")
        rejected = [x for x in result["reference_evidence"]["added_repeat_checks"] if x["mapped_identity"] < .80]
        self.assertEqual(len(rejected), 1)
        self.assertEqual(candidate_geometry(donor)[rejected[0]["index"]], (5670545,5670571))
        self.assertEqual((rejected[0]["mapped_edit_distance"], rejected[0]["mapped_denominator"]), (7,29))

    @staticmethod
    def assess(core, donor):
        left,right = candidate_geometry(core),candidate_geometry(donor)
        mapping = monotone_repeat_mapping(left,right)
        added = [i for i in range(len(right)) if i not in mapping.mapping]
        return assess_added_repeat_reference(core.list_repeats,core.consensus,left,
            donor.list_repeats,right,mapping.mapping,added,.80)

    def test_longer_count_does_not_hide_prior_occurrence_loss(self):
        old = candidate_geometry(self.candidates["published59"])
        self.assertEqual(old[56], (5670714,5670744))
        self.assertEqual(monotone_repeat_mapping(old,candidate_geometry(self.candidates["retained65"])).status,
                         "no_mapping")
        sixty_six = candidate_geometry(self.candidates["retained66"])
        mapped = monotone_repeat_mapping(old,sixty_six)
        self.assertEqual(mapped.status,"unique")
        changed = [(i,a,sixty_six[j]) for i,(a,j) in enumerate(zip(old,mapped.mapping)) if a != sixty_six[j]]
        self.assertEqual(len(changed),1)
        self.assertEqual(changed[0][0],14)
        self.assertEqual(changed[0][1][0],changed[0][2][0])
        self.assertEqual(changed[0][1][1]-changed[0][2][1],2)
        sixty_seven = candidate_geometry(self.candidates["retained67"])
        mapped67 = monotone_repeat_mapping(old,sixty_seven)
        self.assertEqual(mapped67.status,"unique")
        self.assertTrue(all(a == sixty_seven[j] for a,j in zip(old,mapped67.mapping)))


class ReferenceContractTests(unittest.TestCase):
    def setup_values(self):
        left = [(3,29),(63,89),(123,149)]
        right = [(0,29),(60,89),(120,149),(180,207)]
        return [CORE]*3, CORE, left, [PHASE]*3+[WEAK_OLD], right, [0,1,2], [3], .80

    def test_authority_is_mapped_original_core_never_added_majority(self):
        values = list(self.setup_values())
        expected = assess_added_repeat_reference(*values)
        self.assertTrue(expected["qualified"])
        for mutation,status in (("original", "original_core_not_unanimous"),
                                ("mapped", "mapped_core_sequence_not_unanimous"),
                                ("offset", "mapped_core_offsets_not_uniform")):
            changed = copy.deepcopy(values)
            if mutation == "original":
                changed[0][1] = "T" + CORE[1:]
            elif mutation == "mapped":
                changed[3][1] = "T" + PHASE[1:]
            else:
                changed[4][1] = (61,89)
            result = assess_added_repeat_reference(*changed)
            self.assertFalse(result["qualified"])
            self.assertEqual(result["reference_status"],status)
        values[3].extend([WEAK_OLD]*7)
        values[4].extend([(240+60*i,267+60*i) for i in range(7)])
        values[6].extend(range(4,11))
        result = assess_added_repeat_reference(*values)
        self.assertEqual(result["reference_evidence"]["mapped_core_sequence"],PHASE)
        self.assertEqual(result["reference_evidence"]["boundary_offsets"],[-3,0])

    def test_all_copies_use_one_reference_not_per_copy_max(self):
        values = list(self.setup_values())
        # Four substitutions to the old core pass .80 against that core but
        # cost seven edits against the29-base phase; the other copy needs phase.
        old_only = "ACCACGCCCAGGGGCCTGGGCAGCAC"
        values[3].append(old_only)
        values[4].append((240,266))
        values[6].append(4)
        result = assess_added_repeat_reference(*values)
        checks = result["reference_evidence"]["added_repeat_checks"]
        self.assertLess(checks[0]["original_identity"],.80)
        self.assertGreaterEqual(checks[0]["mapped_identity"],.80)
        self.assertGreaterEqual(checks[1]["original_identity"],.80)
        self.assertLess(checks[1]["mapped_identity"],.80)
        self.assertFalse(result["qualified"])

    def test_uniform_source_extensions_trims_and_shifts(self):
        outer = "TTAACGTTACGATCAGTCCGATGCAATGCGCT"
        complement = str.maketrans("ACGT","TGCA")
        for name,left,right,expected in (("extension",(3,29),(0,29),[-3,0]),
                                         ("trim",(0,29),(3,29),[3,0]),
                                         ("shift",(3,29),(0,26),[-3,-3])):
            reference = outer[slice(*right)]
            added = list(reference)
            for index in (5,11,17,23):
                added[index] = added[index].translate(complement)
            added = "".join(added)
            source = (outer + "GTAC"*9)*3 + added
            step = len(outer)+36
            anchor_intervals = [(i*step+left[0],i*step+left[1]) for i in range(3)]
            donor_intervals = [(i*step+right[0],i*step+right[1]) for i in range(3)]
            donor_intervals.append((3*step,3*step+len(added)))
            anchor,donor = from_source(source,anchor_intervals),from_source(source,donor_intervals)
            with self.subTest(name=name):
                self.assertFalse(relationship(anchor,donor,reference_rule=ANCHOR_REFERENCE)["qualified"])
                relation = relationship(anchor,donor)
                self.assertTrue(relation["qualified"])
                ref = relation["added_repeat_reference"]
                self.assertEqual(ref["boundary_offsets"],expected)
                self.assertEqual(ref["mapped_core_sequence"],reference)
                self.assertEqual(ref["added_repeat_checks"][0]["mapped_edit_distance"],4)

    def test_exact_point_eight_identity_boundary(self):
        reference = "ACGTTACGATCAGTCCGATG"
        for added,edits,expected in (("ACCTTAGGATGAGTGCGATG",4,True),
                                     ("ACCTTAGGATGAGTGCGAAG",5,False)):
            source = (("TTA"+reference+"GTAC"*9)*3)+added
            step = 59
            left = [(i*step,i*step+23) for i in range(3)]
            right = [(i*step+3,i*step+23) for i in range(3)]+[(3*step,3*step+20)]
            anchor,donor = from_source(source,left),from_source(source,right)
            result = SourcePhaseTests.assess(anchor,donor)
            self.assertEqual(result["qualified"],expected)
            check = result["reference_evidence"]["added_repeat_checks"][0]
            self.assertEqual((check["mapped_edit_distance"],check["mapped_denominator"]),(edits,20))
            self.assertEqual(check["mapped_identity"],.80 if expected else .75)

    def test_noncanonical_and_strict_mapping_limits_remain(self):
        source,categories,_ = phase_fixture()
        core = next(iter(categories[0].values()))[0][1]
        donor = next(iter(categories[1].values()))[0][1]
        donor.consensus = "T"*29
        self.assertEqual(relationship(core,donor)["added_repeat_reference"]["mapped_core_sequence"],PHASE)
        intervals = list(candidate_geometry(core))
        intervals[0] = (intervals[0][0]+1,intervals[0][1])
        shifted = from_source(source,intervals)
        self.assertEqual(relationship(shifted,donor)["reason"],"no_mapping")
        values = list(self.setup_values())
        values[3][-1] = "N"+WEAK_OLD[1:]
        self.assertEqual(assess_added_repeat_reference(*values)["reason"],"added_repeat_not_canonical")

    def test_legacy_absence_and_unknown_explicit_rule(self):
        self.assertEqual(recorded_reference_rule({}),ANCHOR_REFERENCE)
        good = assess_added_repeat_reference(*self.setup_values())["reference_evidence"]
        self.assertEqual(recorded_reference_rule({"added_repeat_reference":good}),MAPPED_CORE_REFERENCE)
        for bad in (None,{},dict(good,rule="future-unapproved"),dict(good,extra=True)):
            with self.subTest(bad=bad),self.assertRaises(ValueError):
                recorded_reference_rule({"added_repeat_reference":bad})
        with self.assertRaises(ValueError):
            assess_added_repeat_reference(*self.setup_values(),reference_rule="unknown")


class PhaseTransactionTests(unittest.TestCase):
    def run_fixture(self, **kwargs):
        source,categories,scorer = phase_fixture(**kwargs)
        output,audit = reconstruct_accepted_arrays(source,categories,PARAMETERS,scorer,lambda:())
        row = next(iter(output[0].values()))[0]
        self.assertEqual(audit["stats"]["selected"],1,audit)
        return source,categories,output,audit,row

    def test_new_fallback_replay_and_accepted_only_report(self):
        source,_,output,_,row = self.run_fixture()
        basis = row[3]["classification_basis"]
        relation = basis["reconstruction"]["selection"]["relationship"]
        self.assertEqual(relation["added_repeat_reference"]["rule"],MAPPED_CORE_REFERENCE)
        self.assertEqual(basis,validate_reconstruction_basis(source,row[1],basis,row[0],0,row,output))
        projected = [output[0],{},output[2],{},{}]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"phase.fa"
            path.write_text(">phase\n"+source+"\n")
            report = ReportModelBuilder(path,projected,{},[COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS,param_accepted_only=True),provenance_categories=output).build()
        self.assertEqual(report["arrays"][0]["classification_basis"],basis)

    def test_legacy_proof_stays_old_shaped_and_replays(self):
        for fixture in (legacy_fixture,legacy_family_fixture):
            source,categories,scorer = fixture()
            output,audit = reconstruct_accepted_arrays(source,categories,PARAMETERS,scorer,lambda:())
            with mock.patch("components.reconstruction_transaction.MAPPED_CORE_REFERENCE",ANCHOR_REFERENCE):
                legacy,legacy_audit = reconstruct_accepted_arrays(source,categories,PARAMETERS,scorer,lambda:())
            row = next(iter(output[0].values()))[0]
            basis = row[3]["classification_basis"]
            self.assertEqual(basis,next(iter(legacy[0].values()))[0][3]["classification_basis"])
            self.assertEqual(audit,legacy_audit)
            if fixture is legacy_fixture:
                self.assertNotIn("added_repeat_reference",basis["reconstruction"]["selection"]["relationship"])
                validator = validate_reconstruction_basis
            else:
                proof = basis["boundary_family"]
                self.assertTrue(all("added_repeat_reference" not in item["relationship"]
                                    for item in proof["members"]+proof["additional_peer_representations"]))
                validator = validate_boundary_family_basis
            self.assertEqual(basis,validator(source,row[1],basis,row[0],0,row,output))

    def test_missing_forged_and_stale_reference_reject_even_with_new_digest(self):
        source,_,output,_,row = self.run_fixture()
        for mutation in ("remove","rule","sequence","original","offset","map","geometry","distance","denominator","identity"):
            basis = copy.deepcopy(row[3]["classification_basis"])
            proof = basis["reconstruction"]
            relation = proof["selection"]["relationship"]
            ref = relation["added_repeat_reference"]
            if mutation == "remove":
                del relation["added_repeat_reference"]
            elif mutation == "rule":
                ref["rule"] = "unapproved"
            elif mutation == "sequence":
                ref["mapped_core_sequence"] = "T"+PHASE[1:]
            elif mutation == "original":
                ref["original_consensus"] = "T"+CORE[1:]
            elif mutation == "offset":
                ref["boundary_offsets"][0] += 1
            elif mutation == "map":
                ref["mapped_core_indexes"][1] = 2
            elif mutation == "geometry":
                ref["mapped_core_intervals"][0][0] += 1
            else:
                key = {"distance":"mapped_edit_distance","denominator":"mapped_denominator",
                       "identity":"mapped_identity"}[mutation]
                ref["added_repeat_checks"][0][key] += 1
            basis["evidence_sha256"] = digest(proof)
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                validate_reconstruction_basis(source,row[1],basis,row[0],0,row,output)
        with self.assertRaises(ValueError):
            validate_reconstruction_basis("T"+source[1:],row[1],row[3]["classification_basis"],row[0],0,row,output)

    def test_fallback_relationship_rejects_boolean_and_float_coordinate_aliases(self):
        source,_,output,_,row = self.run_fixture()
        for field in ("nested_mapping","nested_offsets","top_mapping","top_added"):
            for convert in (bool,float):
                basis = copy.deepcopy(row[3]["classification_basis"])
                proof = basis["reconstruction"]
                relation = proof["selection"]["relationship"]
                if field == "nested_mapping":
                    target,index = relation["added_repeat_reference"]["mapped_core_indexes"],1
                elif field == "nested_offsets":
                    target,index = relation["added_repeat_reference"]["boundary_offsets"],1
                elif field == "top_mapping":
                    target,index = relation["correspondence"]["anchor_to_candidate_indexes"],1
                else:
                    target,index = relation["added_repeat_indexes"],0
                    # Index3 is not bool-equal; use integral float here to
                    # independently catch top-level numeric type coercion.
                    if convert is bool:
                        continue
                target[index] = convert(target[index])
                basis["evidence_sha256"] = digest(proof)
                with self.subTest(field=field,type=convert.__name__),self.assertRaises(ValueError):
                    validate_reconstruction_basis(source,row[1],basis,row[0],0,row,output)

    def test_fresh_score_and_negative_failures_still_block_fallback(self):
        _source,categories,_ = phase_fixture()
        donor = next(iter(categories[1].values()))[0]
        for negative,score in ((False,donor[0]),(None,donor[0]),(True,.59)):
            def validate(row):
                if row is donor:
                    return replace(qualified([score,row[1],row[2]]),negatives_passed=negative)
                return qualified(row)
            decision = select_reconstructions(categories,validate).decisions[0]
            self.assertNotEqual(decision.action,"selected")
        # A valid original-core correspondence cannot erase a disjoint accepted
        # core that the proposed extension would encompass.
        changed = copy.deepcopy(categories)
        accepted = from_source(_source,candidate_geometry(donor[1])[3:])
        changed[2] = {"other":[[.65,accepted,[]]]}
        decision = select_reconstructions(changed,qualified).decisions[0]
        self.assertNotEqual(decision.action,"selected")

    def test_fallback_cannot_erase_an_original_s1_occurrence_witness(self):
        from components.candidate_origins import mark_repair_origin
        from components.initial_competition import initial_competition_groups
        source,categories,_ = phase_fixture()
        key = next(iter(categories[0]))
        core = categories[0][key][0]
        geometry = candidate_geometry(core[1])
        # Each original witness intersects a26-base core by13 bases. The
        # left-extended29-base donor intersects it by the same13 (< one half).
        witness = [.76,from_source(source,[(a+13,b+13) for a,b in geometry]),[]]
        mark_repair_origin(core[1],{
            "policy":"terminal-only-trim-v1","status":"generated",
            "coordinate_system":"zero_based_half_open",
            "source_sha256":hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry":[[a for a,b in geometry],[b for a,b in geometry]],
        })
        groups,conflicts = initial_competition_groups([witness,core])
        self.assertFalse(conflicts)
        bound = next(r for r in groups[0].rows if r[1] is core[1])
        categories[0][key] = [bound]
        categories[1][key].append(witness)
        decision = select_reconstructions(categories,qualified).decisions[0]
        self.assertEqual(decision.reason,"no_supported_extension")
        self.assertEqual(decision.evidence["geometry_rejections"]["protected_ordinary_occurrence_conflict"],1)

    def test_mixed_family_members_and_alias_replay_their_own_reference(self):
        source,_,output,audit,row = self.run_fixture(family=True,aliases=True)
        basis = row[3]["classification_basis"]
        self.assertEqual(basis["method"],FAMILY_METHOD)
        self.assertEqual(audit["stats"]["boundary_uncertain_selected"],1)
        proof = basis["boundary_family"]
        self.assertEqual(len(proof["members"]),2)
        self.assertEqual(len(proof["additional_peer_representations"]),1)
        records = proof["members"]+proof["additional_peer_representations"]
        rules = [recorded_reference_rule(item["relationship"]) for item in records]
        self.assertEqual(rules.count(ANCHOR_REFERENCE),1)
        self.assertEqual(rules.count(MAPPED_CORE_REFERENCE),2)
        restored = json.loads(json.dumps(basis))
        self.assertEqual(restored,validate_boundary_family_basis(source,row[1],restored,row[0],0,row,output))
        for location in ("members","additional_peer_representations"):
            altered = copy.deepcopy(output)
            current = next(iter(altered[0].values()))[0]
            corrupted = current[3]["classification_basis"]
            bad_proof = corrupted["boundary_family"]
            item = next(m for m in bad_proof[location] if "added_repeat_reference" in m["relationship"])
            del item["relationship"]["added_repeat_reference"]
            corrupted["evidence_sha256"] = digest(bad_proof)
            with self.subTest(location=location),self.assertRaises(ValueError):
                validate_boundary_family_basis(source,current[1],corrupted,current[0],0,current,altered)

    def test_family_member_and_alias_typed_fallback_relationship_replay(self):
        source,_,output,_,_ = self.run_fixture(family=True,aliases=True)
        for location in ("members","additional_peer_representations"):
            for field in ("nested_mapping","top_mapping"):
                altered = copy.deepcopy(output)
                row = next(iter(altered[0].values()))[0]
                basis = row[3]["classification_basis"]
                proof = basis["boundary_family"]
                relation = next(m["relationship"] for m in proof[location]
                                if "added_repeat_reference" in m["relationship"])
                if field == "nested_mapping":
                    relation["added_repeat_reference"]["mapped_core_indexes"][1] = True
                else:
                    relation["correspondence"]["anchor_to_candidate_indexes"][1] = 1.0
                basis["evidence_sha256"] = digest(proof)
                with self.subTest(location=location,field=field),self.assertRaises(ValueError):
                    validate_boundary_family_basis(source,row[1],basis,row[0],0,row,altered)


if __name__ == "__main__":
    unittest.main()
