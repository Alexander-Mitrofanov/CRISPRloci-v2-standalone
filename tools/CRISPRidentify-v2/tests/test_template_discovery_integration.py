import copy
import hashlib
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from components.detection.repeat_templates import (
    TemplateDiscoveryError, load_repeat_template_panel, panel_configuration, reverse_complement,
)
from components.detection.template_discovery import validate_discovery_ledger
from components.sequence.validation import read_single_fasta
from tests.test_batch_runner import batch

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data/template_discovery"
MANIFEST = json.loads((DATA / "manifest.json").read_text())


def scientific(report):
    return [{key: a.get(key) for key in ("category", "certainty_score", "feature_sets", "repeats", "spacers",
                                   "source_interval", "classification_basis", "initial_competition")}
            for a in report["arrays"]]


class TemplateIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.folder = Path(cls.temporary.name)
        cls.serial = 0
        cls.reports = {}
        for source, panel in (("three_polymorphic_columns_30nt", "three_copy"), ("four_copy_known_like", "four_copy")):
            cls.reports[source] = cls.run_cli(DATA / (source + ".fa"), DATA / (panel + "_panel.fa"))

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    @classmethod
    def run_cli(cls, source, panel=None, extra=(), expect=0):
        cls.serial += 1
        output = cls.folder / ("out" + str(cls.serial))
        command = [str(ROOT / "crispridentify"), "--file", str(source), "--result_folder", str(output),
                   "--parallel", "False", "--cpu", "1", "--strand", "True" if "full" in extra and "--preset" in extra else "False", "--cas", "False"]
        if panel:
            command += ["--repeat-template-fasta", str(panel)]
        result = subprocess.run(command + list(extra), cwd=ROOT, capture_output=True, text=True)
        if expect == 0 and result.returncode != 0:
            raise AssertionError(result.stdout + result.stderr)
        if expect != 0:
            if result.returncode == 0 or list(output.rglob("report.json")):
                raise AssertionError("failed template input published success")
            return result
        paths = list(output.rglob("report.json"))
        return json.loads(paths[0].read_text())

    def test_frozen_fixture_hashes_and_exact_positive_countercontrol(self):
        for name, record in MANIFEST["fixtures"].items():
            self.assertEqual(record["sha256"], hashlib.sha256((DATA / name).read_bytes()).hexdigest())
        for key, source in (("template-injection-result", "three_polymorphic_columns_30nt"),
                            ("four-copy-template-injection-result", "four_copy_known_like")):
            expected = MANIFEST["injection_probes"][key]["arrays"][0]
            report = self.reports[source]
            self.assertEqual(1, len(report["arrays"]))
            actual = report["arrays"][0]
            for field in ("repeats", "spacers", "category", "feature_sets"):
                self.assertEqual(expected[field], actual[field], field)
            self.assertAlmostEqual(expected["certainty_score"], actual["certainty_score"], places=14)
            self.assertEqual(int(source.startswith("four")), report["summary"]["accepted_array_count"])
            self.assertEqual(32, actual["spacers"][0]["source_interval"]["length"])

    def test_balanced_sensitive_full_and_regex_controls(self):
        source = DATA / "four_copy_known_like.fa"
        baseline = self.run_cli(source)
        self.assertEqual(0, baseline["summary"]["array_count"])
        sensitive = self.run_cli(source, extra=("--fast_run_seed_profile", "sensitive"))
        self.assertEqual(1, sensitive["summary"]["accepted_array_count"])
        self.assertAlmostEqual(.6204423621476783, sensitive["arrays"][0]["certainty_score"], places=14)
        self.assertEqual(21, sensitive["arrays"][0]["repeats"][0]["source_interval"]["length"])
        for extra in (("--preset", "full"), ("--fuzzy_search_mode", "full")):
            report = self.run_cli(source, DATA / "four_copy_panel.fa", extra)
            self.assertEqual(1, report["summary"]["accepted_array_count"])
            self.assertEqual([81, 143, 205, 267], [r["source_interval"]["start"] for r in report["arrays"][0]["repeats"]])
        mixed = self.run_cli(source, DATA / "four_copy_panel.fa", ("--fast_run_seed_profile", "sensitive"))
        self.assertEqual(1, mixed["summary"]["accepted_array_count"])

    def test_existing_count_and_low_score_contracts(self):
        report = self.run_cli(DATA / "eleven_base_core_34nt.fa", extra=("--fast_run_seed_profile", "sensitive"))
        self.assertEqual("Low score", report["arrays"][0]["category"])
        self.assertAlmostEqual(.37789706497487846, report["arrays"][0]["certainty_score"], places=14)
        baseline = self.run_cli(DATA / "two_known_repeats.fa")
        allowed = self.run_cli(DATA / "two_known_repeats.fa", extra=("--min_repeats", "2"))
        self.assertEqual("Low score", baseline["arrays"][0]["category"])
        self.assertEqual("Possible", allowed["arrays"][0]["category"])
        self.assertEqual(baseline["arrays"][0]["repeats"], allowed["arrays"][0]["repeats"])
        self.assertAlmostEqual(.6571611015115605, allowed["arrays"][0]["certainty_score"], places=14)
        invalid = self.run_cli(DATA / "four_copy_known_like.fa", DATA / "four_copy_panel.fa", ("--min_len_rep", "40"))
        self.assertEqual(0, invalid["summary"]["accepted_array_count"])

    def test_reverse_source_and_reverse_only_duplicate_panels(self):
        _, source = read_single_fasta(DATA / "four_copy_known_like.fa")
        original = load_repeat_template_panel(DATA / "four_copy_panel.fa").aliases[0].sequence
        reverse_panel = self.folder / "reverse duplicate panel.fa"
        reverse_panel.write_text(">reverse\n" + reverse_complement(original) + "\n>duplicate\n" + original.lower() + "\n")
        forward = self.run_cli(DATA / "four_copy_known_like.fa", reverse_panel)
        self.assertEqual(scientific(self.reports["four_copy_known_like"]), scientific(forward))
        reverse_source = self.folder / "reverse.fa"
        reverse_source.write_text(">reverse\n" + reverse_complement(source) + "\n")
        reverse_panel.write_text(">reverse\n" + reverse_complement(original) + "\n")
        report = self.run_cli(reverse_source, reverse_panel)
        self.assertEqual(1, report["summary"]["accepted_array_count"])
        observed = [(r["source_interval"]["start"] - 1, r["source_interval"]["end"]) for r in report["arrays"][0]["repeats"]]
        self.assertEqual(sorted((len(source)-b, len(source)-a) for a, b in ((80,110),(142,172),(204,234),(266,296))), observed)

    def test_no_panel_does_not_call_template_lane(self):
        # Execute the real CLI in-process so accidental loader/scanner calls fail.
        script = """import runpy,sys
from unittest.mock import patch
with patch('components.detection.repeat_templates.load_repeat_template_panel', side_effect=AssertionError('loader')), patch('components.detection.stage.Detection._run_template_detection', side_effect=AssertionError('scanner')), patch('components.refinement.template_integration.reconcile_template_candidates', side_effect=AssertionError('grouping')):
    sys.argv = sys.argv[1:]
    runpy.run_path(sys.argv[0], run_name='__main__')
"""
        env = dict(os.environ, PYTHONHASHSEED="0",
                   PATH=str(ROOT / ".envs/crispridentify-core/bin") + os.pathsep + os.environ["PATH"],
                   GTHDATADIR=str(ROOT / ".envs/crispridentify-core/share/genometools"))
        result = subprocess.run([sys.executable, "-c", script, str(ROOT / "CRISPRidentify.py"), "--file",
            str(DATA / "four_copy_known_like.fa"), "--result_folder", str(self.folder / "disabled"),
            "--parallel", "False", "--strand", "False", "--cpu", "1"], env=env, cwd=ROOT,
            capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)

    def test_invalid_panel_or_workflow_never_publishes(self):
        invalid = self.folder / "empty.fa"
        invalid.write_text("")
        self.run_cli(DATA / "four_copy_known_like.fa", invalid, expect=1)
        for extra in (("--preset", "region-screen"), ("--accepted_only", "True"), ("--fuzzy_search_mode", "family")):
            self.run_cli(DATA / "four_copy_known_like.fa", DATA / "four_copy_panel.fa", extra, expect=1)
        self.run_cli(DATA / "four_copy_known_like.fa", DATA / "four_copy_panel.fa",
                     ("--_repeat_template_sha256", "0" * 64), expect=1)

    def test_negative_controls_have_no_new_accepted_events(self):
        motif = load_repeat_template_panel(DATA / "four_copy_panel.fa").aliases[0].sequence
        rng = random.Random(1007)
        controls = {
            "isolated": "N" * 80 + motif + "N" * 80,
            "inverted": "N" * 80 + motif + "N" * 32 + reverse_complement(motif) + "N" * 80,
            "low_complexity": "A" * 2000,
            "random": "".join(rng.choices("ACGT", k=2000)),
            "identical_spacers": "N" * 80 + (motif + "ACGT" * 8) * 8 + motif + "N" * 80,
        }
        for name, source in controls.items():
            path = self.folder / (name + ".fa")
            path.write_text(">" + name + "\n" + source + "\n")
            report = self.run_cli(path, DATA / "four_copy_panel.fa")
            self.assertEqual(0, report["summary"]["accepted_array_count"], name)
            audit = report["run"]["parameters"]["param_repeat_template_discovery"]
            self.assertEqual(len(audit["births"]), audit["candidate_count"])
        shuffled = list(motif)
        rng.shuffle(shuffled)
        query = self.folder / "shuffled.fa"
        query.write_text(">shuffled\n" + "".join(shuffled) + "\n")
        report = self.run_cli(DATA / "four_copy_known_like.fa", query)
        self.assertEqual(0, report["summary"]["accepted_array_count"])

    def test_launcher_panel_path_is_relative_to_callers_directory(self):
        caller = self.folder / "caller with spaces"
        caller.mkdir()
        (caller / "panel.fa").write_bytes((DATA / "four_copy_panel.fa").read_bytes())
        for option in ("--repeat-template-fasta", "--repeat_template_fasta", "--repeat-t", "--repeat_t"):
            for args in ((option, "panel.fa"), (option + "=panel.fa",)):
                result = subprocess.run([str(ROOT / "crispridentify"), "--print-effective-config", *args],
                                        cwd=caller, capture_output=True, text=True)
                self.assertEqual(0, result.returncode, result.stderr)
                config = json.loads(result.stdout)["effective_config"]["repeat_template_discovery"]
                self.assertEqual(load_repeat_template_panel(caller / "panel.fa").raw_sha256,
                                 config["panel"]["raw_sha256"])

    def test_mutating_original_after_parent_snapshot_does_not_change_folder_children(self):
        query = self.folder / "mutable panel.fa"
        query.write_bytes((DATA / "four_copy_panel.fa").read_bytes())
        expected = load_repeat_template_panel(query).raw_sha256
        source = self.folder / "snapshot records.fa"
        text = (DATA / "four_copy_known_like.fa").read_text()
        source.write_text(text + text.replace(">four_copy_known_like", ">second"))
        output = self.folder / "snapshot children"
        script = """import runpy,sys
from pathlib import Path
from unittest.mock import patch
import components.detection.repeat_templates as panels
original = panels.snapshot_repeat_template_panel
query = Path(sys.argv[2])
def snapshot(panel, destination):
    result = original(panel, destination)
    query.write_text('>changed\\n' + 'A'*30 + '\\n')
    return result
with patch.object(panels, 'snapshot_repeat_template_panel', snapshot):
    sys.argv = [sys.argv[1]] + sys.argv[3:]
    runpy.run_path(sys.argv[0], run_name='__main__')
"""
        env = dict(os.environ, PYTHONHASHSEED="0",
                   PATH=str(ROOT / ".envs/crispridentify-core/bin") + os.pathsep + os.environ["PATH"])
        result = subprocess.run([sys.executable, "-c", script, str(ROOT / "CRISPRidentify.py"), str(query),
            "--file", str(source), "--result_folder", str(output), "--repeat-template-fasta", str(query),
            "--folder_jobs", "2", "--cpu", "2", "--strand", "False", "--parallel", "False"],
            env=env, cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)
        reports = [json.loads(p.read_text()) for p in output.rglob("report.json")]
        self.assertEqual(2, len(reports))
        for report in reports:
            self.assertEqual(expected, report["run"]["configuration"]["effective_config"]["repeat_template_discovery"]["panel"]["raw_sha256"])
            self.assertEqual(1, report["summary"]["accepted_array_count"])
        self.assertNotEqual(expected, load_repeat_template_panel(query).raw_sha256)

    def test_search_cap_and_corrupt_ledger_fail_before_canonical_publication(self):
        for failure in ("cap", "ledger"):
            output = self.folder / ("failed-" + failure)
            script = """import runpy,sys,dataclasses
from unittest.mock import patch
import components.detection.repeat_templates as panels
failure = sys.argv[2]
if failure == 'cap':
    original = panels.policy_for_engine
    def policy(engine):
        return dataclasses.replace(original(engine), hits_per_pattern=3)
    context = patch.object(panels, 'policy_for_engine', policy)
else:
    from components.detection.stage import Detection
    original = Detection._run_template_detection
    def corrupt(self):
        original(self)
        self.template_discovery_audit['source_sha256'] = '0'*64
    context = patch.object(Detection, '_run_template_detection', corrupt)
with context:
    sys.argv = [sys.argv[1]] + sys.argv[3:]
    runpy.run_path(sys.argv[0], run_name='__main__')
"""
            env = dict(os.environ, PYTHONHASHSEED="0",
                       PATH=str(ROOT / ".envs/crispridentify-core/bin") + os.pathsep + os.environ["PATH"])
            result = subprocess.run([sys.executable, "-c", script, str(ROOT / "CRISPRidentify.py"), failure,
                "--file", str(DATA / "four_copy_known_like.fa"), "--result_folder", str(output),
                "--repeat-template-fasta", str(DATA / "four_copy_panel.fa"), "--parallel", "False",
                "--strand", "False", "--cpu", "1"], env=env, cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(0, result.returncode)
            self.assertFalse(list(output.rglob("report.json")))
            self.assertIn("limit" if failure == "cap" else "discovery source differs", result.stderr)

    def test_ledger_replays_and_rejects_tampered_geometry_policy_or_missing_evidence(self):
        report = self.reports["four_copy_known_like"]
        audit = report["run"]["parameters"]["param_repeat_template_discovery"]
        config = report["run"]["configuration"]["effective_config"]["repeat_template_discovery"]
        _, source = read_single_fasta(DATA / "four_copy_known_like.fa")
        validate_discovery_ledger(audit, source, config)
        for key, value in (("hit_count", 0), ("source_sha256", "0" * 64), ("tasks", []),
                           ("patterns", []), ("policy", "unknown")):
            changed = copy.deepcopy(audit)
            changed[key] = value
            with self.assertRaises(TemplateDiscoveryError):
                validate_discovery_ledger(changed, source, config)
        changed = copy.deepcopy(audit)
        changed["births"][0]["representation"][0][0] += 1
        with self.assertRaises(TemplateDiscoveryError):
            validate_discovery_ledger(changed, source, config)
        with self.assertRaises(TemplateDiscoveryError):
            validate_discovery_ledger(None, source, config)

    def test_all_rescue_seed_scopes_and_folder_record_isolation(self):
        for scope in ("lineage", "locus", "family", "cluster"):
            report = self.run_cli(DATA / "three_polymorphic_columns_30nt.fa", DATA / "three_copy_panel.fa",
                                   ("--fast_run_rescue_seed_scope", scope))
            self.assertEqual(scientific(self.reports["three_polymorphic_columns_30nt"]), scientific(report))
        multifasta = self.folder / "two records.fa"
        multifasta.write_text((DATA / "four_copy_known_like.fa").read_text() + ">empty\n" + "N" * 100 + "\n")
        hashes = []
        for lanes in (1, 2):
            output = self.folder / ("lanes" + str(lanes))
            command = [str(ROOT / "crispridentify"), "--file", str(multifasta), "--result_folder", str(output),
                       "--repeat_template_fasta=" + str(DATA / "four_copy_panel.fa"),
                       "--folder_jobs", str(lanes), "--cpu", "2", "--strand", "False", "--parallel", "False"]
            result = subprocess.run(command, capture_output=True, text=True, cwd=ROOT)
            self.assertEqual(0, result.returncode, result.stderr)
            reports = [json.loads(p.read_text()) for p in output.rglob("report.json")]
            self.assertEqual(2, len(reports))
            self.assertEqual([0, 1], sorted(r["summary"]["accepted_array_count"] for r in reports))
            hashes.extend(r["run"]["parameters"]["param_repeat_template_discovery"]["panel"]["raw_sha256"] for r in reports)
        self.assertEqual(1, len(set(hashes)))

    def test_batch_managed_panel_snapshot_resume_and_wrong_child_identity(self):
        for option in ("--repeat-template-fasta", "--repeat_template_fasta", "--repeat-t", "--repeat_t", "--_repeat_template"):
            for args in ((option, "query.fa"), (option + "=query.fa",)):
                with self.assertRaises(batch.InputError):
                    batch._validate_extra_arguments(args)
        report = copy.deepcopy(self.reports["four_copy_known_like"])
        config = panel_configuration(load_repeat_template_panel(DATA / "four_copy_panel.fa"), "native")
        report_path = self.folder / "child.json"
        report["run"]["configuration"]["effective_config"]["repeat_template_discovery"]["panel"]["raw_sha256"] = "0" * 64
        report_path.write_text(json.dumps(report))
        with self.assertRaisesRegex(batch.CompletionError, "panel identity"):
            batch.validate_report(report_path, None, None, False, expected_repeat_template_config=config)
        inputs = self.folder / "batch inputs"
        inputs.mkdir()
        (inputs / "sample.fa").write_bytes((DATA / "four_copy_known_like.fa").read_bytes())
        panel_path = self.folder / "batch panel.fa"
        panel_path.write_bytes((DATA / "four_copy_panel.fa").read_bytes())
        manifest = self.folder / "batch.tsv"
        command = [sys.executable, str(ROOT / "scripts/run_crispridentify_batch.py")]
        subprocess.run(command + ["manifest", "--input-root", str(inputs), "--output", str(manifest)],
                       check=True, capture_output=True)
        run = self.folder / "batch run"
        argv = command + ["run", "--manifest", str(manifest), "--run-dir", str(run), "--preset", "balanced",
                          "--cpus-per-worker", "1", "--repeat-template-fasta", str(panel_path), "--min-free-gib", "0"]
        for _ in range(2):
            result = subprocess.run(argv, capture_output=True, text=True, cwd=ROOT)
            self.assertEqual(0, result.returncode, result.stderr + result.stdout)
        self.assertEqual(1, len(list(run.rglob("completion.json"))))
        self.assertEqual(1, len(list(run.rglob("attempt.json"))))
        panel_path.write_text(">changed\n" + "A" * 30 + "\n")
        result = subprocess.run(argv, capture_output=True, text=True, cwd=ROOT)
        self.assertNotEqual(0, result.returncode)
        self.assertIn("different configuration", result.stderr)
        self.assertEqual(1, len(list(run.rglob("attempt.json"))))
