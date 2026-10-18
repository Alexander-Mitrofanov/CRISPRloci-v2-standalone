import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_crispridentify_batch.py"
SPEC = importlib.util.spec_from_file_location("crispridentify_batch_under_test", SCRIPT)
batch = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = batch
SPEC.loader.exec_module(batch)


FAKE_RUNNER = r"""#!/usr/bin/env python3
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

if "--help" in sys.argv:
    print("fake detector --preset {balanced,region-screen,full}")
    raise SystemExit(0)

def value(name):
    return sys.argv[sys.argv.index(name) + 1]

fasta = Path(value("--file"))
result = Path(value("--result_folder"))
counter_root = Path(os.environ["FAKE_COUNTER_ROOT"])
counter_root.mkdir(parents=True, exist_ok=True)
counter = counter_root / (fasta.name + ".count")
count = int(counter.read_text() if counter.exists() else "0") + 1
counter.write_text(str(count))

text = fasta.read_text(encoding="ascii")
header = text.splitlines()[0][1:].split()[0]
if header == "RETRY_ONCE" and count == 1:
    raise SystemExit(7)
if header == "FAIL_ALWAYS":
    raise SystemExit(9)
if header == "SLEEP":
    ready = counter_root / "sleep-child.ready"
    child = subprocess.Popen([
        sys.executable,
        "-c",
        "import pathlib,signal,sys,time; "
        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
        "pathlib.Path(sys.argv[1]).write_text('ready'); time.sleep(30)",
        str(ready),
    ])
    deadline = time.monotonic() + 1
    while not ready.exists() and time.monotonic() < deadline:
        time.sleep(0.005)
    (counter_root / "sleep-child.pid").write_text(str(child.pid))
    time.sleep(30)

records = []
current = None
parts = []
for raw in text.splitlines():
    if raw.startswith(">"):
        if current is not None:
            records.append((current, "".join(parts).upper()))
        current = raw[1:].split()[0]
        parts = []
    elif raw.strip():
        parts.append("".join(raw.split()))
if current is not None:
    records.append((current, "".join(parts).upper()))

for source_id, sequence in records:
    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", source_id).strip("_")
    target = result / safe
    target.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
        "source": {
            "id": source_id,
            "length": len(sequence),
            "sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        },
        "run": {
            "parameters": {
                "param_runtime_preset": "region-screen",
                "param_accepted_only": True,
            },
            "configuration": {"preset": "region-screen"},
        },
        "summary": {
            "array_count": 0,
            "accepted_array_count": 0,
            "category_counts": {
                "Bona-fide": 0,
                "Possible": 0,
                "Low score": 0,
            },
        },
        "validation": {
            "status": "valid",
            "method": "independent_original_sequence_slicing",
            "validated_array_count": 0,
            "invalid_array_count": 0,
        },
        "arrays": [],
        "extensions": {
            "likely_tandem": {
                "schema": {
                    "name": "CRISPRidentify-likely-tandem-quarantine",
                    "version": "1.1.0",
                },
                "source_sha256": hashlib.sha256(
                    sequence.encode("ascii")
                ).hexdigest(),
                "disposition": "likely_tandem",
                "tag": "tandem",
                "reported_category": "Low score",
                "status": "applied",
                "external_tool_support_used_for_decision": False,
                "rule": {
                    "name": "complete-link-low-spacer-diversity",
                    "version": "1.0.0",
                    "linkage": "complete_link_all_spacer_pairs",
                    "identity": (
                        "1 - global_levenshtein/max(sequence_lengths)"
                    ),
                    "policy": {
                        "minimum_spacer_count": 2,
                        "maximum_spacer_count": 8,
                        "minimum_pairwise_edit_identity": 0.95,
                        "canonical_alphabet": "ACGT",
                    },
                },
                "transaction_stats": {
                    "total_rows_before": 0,
                    "total_rows_after": 0,
                    "accepted_rows_reviewed": 0,
                    "accepted_rows_before": 0,
                    "accepted_rows_after": 0,
                    "likely_tandem_count": 0,
                    "retained_count": 0,
                    "bona_fide_quarantined": 0,
                    "possible_quarantined": 0,
                    "diagnostic_rows_before": 0,
                    "diagnostic_rows_after": 0,
                },
                "array_count": 0,
                "arrays": [],
            }
        },
    }
    (target / "report.json").write_text(
        json.dumps(report, sort_keys=True) + "\n", encoding="utf-8"
    )
"""


def write_fasta(path, source_id, sequence="ACGTACGTACGT"):
    path.write_text(">{} synthetic\n{}\n".format(source_id, sequence), encoding="ascii")


def _current_tandem_extension(source_sequence, stats, records=()):
    return {
        "schema": {
            "name": "CRISPRidentify-likely-tandem-quarantine",
            "version": "1.1.0",
        },
        "source_sha256": hashlib.sha256(
            source_sequence.encode("ascii")
        ).hexdigest(),
        "disposition": "likely_tandem",
        "tag": "tandem",
        "reported_category": "Low score",
        "status": "applied",
        "external_tool_support_used_for_decision": False,
        "rule": {
            "name": "complete-link-low-spacer-diversity",
            "version": "1.0.0",
            "linkage": "complete_link_all_spacer_pairs",
            "identity": "1 - global_levenshtein/max(sequence_lengths)",
            "policy": {
                "minimum_spacer_count": 2,
                "maximum_spacer_count": 8,
                "minimum_pairwise_edit_identity": 0.95,
                "canonical_alphabet": "ACGT",
            },
        },
        "transaction_stats": stats,
        "array_count": len(records),
        "arrays": list(records),
    }


def _retained_stats(accepted_count):
    return {
        "total_rows_before": accepted_count,
        "total_rows_after": accepted_count,
        "accepted_rows_reviewed": accepted_count,
        "accepted_rows_before": accepted_count,
        "accepted_rows_after": accepted_count,
        "likely_tandem_count": 0,
        "retained_count": accepted_count,
        "bona_fide_quarantined": 0,
        "possible_quarantined": 0,
        "diagnostic_rows_before": 0,
        "diagnostic_rows_after": 0,
    }


def _canonical_array(first_spacer, second_spacer, category="Possible"):
    repeat = "ACGTAC"
    source_sequence = (
        repeat + first_spacer + repeat + second_spacer + repeat
    )
    first_repeat_end = len(repeat)
    first_spacer_end = first_repeat_end + len(first_spacer)
    second_repeat_end = first_spacer_end + len(repeat)
    second_spacer_end = second_repeat_end + len(second_spacer)
    repeats = [
        {
            "index": 1,
            "source_interval": {
                "start": 1,
                "end": first_repeat_end,
                "length": len(repeat),
            },
            "sequence": repeat,
        },
        {
            "index": 2,
            "source_interval": {
                "start": first_spacer_end + 1,
                "end": second_repeat_end,
                "length": len(repeat),
            },
            "sequence": repeat,
        },
        {
            "index": 3,
            "source_interval": {
                "start": second_spacer_end + 1,
                "end": len(source_sequence),
                "length": len(repeat),
            },
            "sequence": repeat,
        },
    ]
    spacers = [
        {
            "index": 1,
            "source_interval": {
                "start": first_repeat_end + 1,
                "end": first_spacer_end,
                "length": len(first_spacer),
            },
            "sequence": first_spacer,
        },
        {
            "index": 2,
            "source_interval": {
                "start": second_repeat_end + 1,
                "end": second_spacer_end,
                "length": len(second_spacer),
            },
            "sequence": second_spacer,
        },
    ]
    architecture = []
    for position, repeat_record in enumerate(repeats):
        architecture.append(dict(repeat_record, type="repeat"))
        if position < len(spacers):
            architecture.append(dict(spacers[position], type="spacer"))
    return source_sequence, {
        "id": "candidate",
        "category": category,
        "certainty_score": 0.55,
        "source_interval": {
            "start": 1,
            "end": len(source_sequence),
            "length": len(source_sequence),
        },
        "source_sequence": source_sequence,
        "repeat_count": 3,
        "spacer_count": 2,
        "non_deleted_spacer_count": 2,
        "spacer_deletion_count": 0,
        "repeats": repeats,
        "spacers": spacers,
        "spacer_deletions": [],
        "architecture": architecture,
        "validation": {
            "status": "valid",
            "source_reconstructed": True,
            "method": "independent_original_sequence_slicing",
        },
    }


def _current_report(source_sequence, arrays, accepted_only=True):
    counts = {
        category: sum(array["category"] == category for array in arrays)
        for category in ("Bona-fide", "Possible", "Low score")
    }
    report = {
        "schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
        "source": {
            "id": "A",
            "length": len(source_sequence),
            "sha256": hashlib.sha256(
                source_sequence.encode("ascii")
            ).hexdigest(),
        },
        "run": {
            "parameters": {"param_accepted_only": accepted_only},
            "configuration": {"preset": "region-screen"},
        },
        "summary": {
            "array_count": len(arrays),
            "accepted_array_count": counts["Bona-fide"] + counts["Possible"],
            "category_counts": counts,
        },
        "validation": {
            "status": "valid",
            "method": "independent_original_sequence_slicing",
            "validated_array_count": len(arrays),
            "invalid_array_count": 0,
        },
        "arrays": arrays,
        "extensions": {
            "likely_tandem": _current_tandem_extension(
                source_sequence,
                _retained_stats(counts["Bona-fide"] + counts["Possible"]),
            )
        },
    }
    return report


def _identical_tandem_evidence(spacer):
    return {
        "decision": "likely_tandem",
        "category": "likely_tandem",
        "reason_code": "all_spacer_pairs_at_or_above_identity_threshold",
        "rule_name": "complete-link-low-spacer-diversity",
        "rule_version": "1.0.0",
        "coordinate_system": "source_forward_strand",
        "policy": {
            "minimum_spacer_count": 2,
            "maximum_spacer_count": 8,
            "minimum_pairwise_edit_identity": 0.95,
            "canonical_alphabet": "ACGT",
        },
        "spacer_count": 2,
        "nonempty_spacer_count": 2,
        "canonical_spacer_count": 2,
        "distinct_spacer_count": 1,
        "exact_duplicate_excess_count": 1,
        "shortest_spacer_length": len(spacer),
        "longest_spacer_length": len(spacer),
        "pairwise_comparison_count": 1,
        "minimum_pairwise_edit_identity": 1.0,
        "maximum_pairwise_edit_identity": 1.0,
        "maximum_pairwise_edit_distance": 0,
        "minimum_identity_pair_1based": [1, 2],
        "minimum_identity_pair_edit_distance": 0,
        "minimum_identity_pair_maximum_length": len(spacer),
        "complete_link_threshold_satisfied": True,
    }


def _accepted_only_tandem_report():
    repeat = "ACGTAC"
    spacer = "A" * 20
    source_sequence = repeat + spacer + repeat + spacer + repeat
    starts = (0, 26, 52)

    def ledger_element(index, start, sequence):
        end = start + len(sequence)
        return {
            "index": index,
            "sequence": sequence,
            "sequence_sha256": hashlib.sha256(
                sequence.encode("ascii")
            ).hexdigest(),
            "source_interval": {"start": start, "end": end},
        }

    repeats = [
        ledger_element(index, start, repeat)
        for index, start in enumerate(starts, 1)
    ]
    spacers = [
        ledger_element(1, 6, spacer),
        ledger_element(2, 32, spacer),
    ]
    record = {
        "id": "LT1",
        "action": "quarantined",
        "disposition": "likely_tandem",
        "tag": "tandem",
        "source_category_index": 2,
        "source_category_name": "Possible",
        "target_category_index": 4,
        "target_category_name": "Low score",
        "row_ordinal": 0,
        "score": 0.55,
        "locus_key": [0, len(source_sequence)],
        "repeat_count": 3,
        "source_interval": {
            "start": 1,
            "end": len(source_sequence),
            "length": len(source_sequence),
        },
        "spacers_source_forward": [spacer, spacer],
        "evidence": _identical_tandem_evidence(spacer),
        "source_architecture": {
            "coordinate_system": "zero_based_half_open",
            "reference": "input_forward_strand",
            "repeats": repeats,
            "spacers": spacers,
            "reconstructed_source_span_sha256": hashlib.sha256(
                source_sequence.encode("ascii")
            ).hexdigest(),
        },
        "reported_array_id": None,
        "core_projection": "omitted_by_accepted_only",
    }
    stats = {
        "total_rows_before": 1,
        "total_rows_after": 1,
        "accepted_rows_reviewed": 1,
        "accepted_rows_before": 1,
        "accepted_rows_after": 0,
        "likely_tandem_count": 1,
        "retained_count": 0,
        "bona_fide_quarantined": 0,
        "possible_quarantined": 1,
        "diagnostic_rows_before": 0,
        "diagnostic_rows_after": 1,
    }
    report = _current_report(source_sequence, [])
    report["extensions"]["likely_tandem"] = _current_tandem_extension(
        source_sequence, stats, (record,)
    )
    return source_sequence, report


def _full_tandem_report():
    source_sequence, accepted_only_report = _accepted_only_tandem_report()
    record = copy.deepcopy(
        accepted_only_report["extensions"]["likely_tandem"]["arrays"][0]
    )
    rebuilt_sequence, array = _canonical_array(
        "A" * 20, "A" * 20, category="Low score"
    )
    if rebuilt_sequence != source_sequence:
        raise AssertionError("full tandem test fixture source differs")
    record["score"] = 1.0
    record["reported_array_id"] = array["id"]
    record["core_projection"] = "reported_as_low_score"
    array["certainty_score"] = 1.0
    array["tags"] = ["tandem"]
    array["legacy_interval_key"] = copy.deepcopy(record["locus_key"])
    array["classification_basis"] = {
        "method": "low_spacer_diversity_tandem_demotion",
        "disposition": "likely_tandem",
        "tag": "tandem",
        "likely_tandem": copy.deepcopy(record["evidence"]),
        "transition": {
            "rule_name": "likely-tandem-category-quarantine",
            "rule_version": "1.1.0",
            "source_category_index": 2,
            "source_category_name": "Possible",
            "target_category_index": 4,
            "target_category_name": "Low score",
            "row_ordinal": 0,
            "locus_key": copy.deepcopy(record["locus_key"]),
            "reported_disposition": "likely_tandem",
            "reported_tag": "tandem",
        },
    }
    report = _current_report(source_sequence, [array], accepted_only=False)
    report["run"]["configuration"]["preset"] = "balanced"
    report["extensions"]["likely_tandem"] = _current_tandem_extension(
        source_sequence,
        copy.deepcopy(
            accepted_only_report["extensions"]["likely_tandem"][
                "transaction_stats"
            ]
        ),
        (record,),
    )
    return source_sequence, report


class BatchTestCase(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.inputs = self.root / "inputs"
        self.inputs.mkdir()
        self.runner = self.root / "fake detector"
        self.runner.write_text(FAKE_RUNNER, encoding="utf-8")
        self.runner.chmod(0o755)
        self.counter_root = self.root / "counters"

    def make_manifest(self):
        manifest = self.root / "manifest.tsv"
        batch.create_manifest(self.inputs, manifest, (".fna",))
        return manifest

    def run_batch(self, manifest, run_dir=None, **overrides):
        run_dir = run_dir or (self.root / "run")
        arguments = [
            "run",
            "--manifest",
            str(manifest),
            "--run-dir",
            str(run_dir),
            "--runner",
            str(self.runner),
            "--workers",
            str(overrides.pop("workers", 2)),
            "--cpus-per-worker",
            "1",
            "--timeout-seconds",
            str(overrides.pop("timeout", 2)),
            "--termination-grace-seconds",
            "0.1",
            "--max-attempts",
            str(overrides.pop("max_attempts", 2)),
            "--scratch-root",
            str(self.root / "scratch"),
            "--min-free-gib",
            "0",
        ]
        for key, value in overrides.items():
            arguments.extend(("--" + key.replace("_", "-"), str(value)))
        environment = {"FAKE_COUNTER_ROOT": str(self.counter_root)}
        with (
            mock.patch.dict(os.environ, environment, clear=False),
            mock.patch.object(batch, "runtime_tree_sha256", return_value="a" * 64),
            mock.patch.object(batch, "_runtime_commit", return_value="b" * 40),
        ):
            return batch.main(arguments)


class ManifestTests(BatchTestCase):
    def test_manifest_is_sorted_hashed_and_repeatable(self):
        nested = self.inputs / "z"
        nested.mkdir()
        write_fasta(nested / "b.fna", "B")
        write_fasta(self.inputs / "a.fna", "A", "TTTT")
        ignored = self.inputs / "notes.txt"
        ignored.write_text("not FASTA", encoding="utf-8")
        manifest = self.make_manifest()
        first = manifest.read_bytes()
        self.assertEqual(2, batch.create_manifest(self.inputs, manifest, ("fna",)))
        self.assertEqual(first, manifest.read_bytes())
        with manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(["a", "z/b"], [row["sample_id"] for row in rows])
        for row in rows:
            source = (manifest.parent / row["fasta"]).resolve()
            self.assertEqual(source.stat().st_size, int(row["bytes"]))
            self.assertEqual(
                hashlib.sha256(source.read_bytes()).hexdigest(), row["sha256"]
            )

    def test_sharding_is_stable_and_disjoint(self):
        samples = ["sample-{:04d}".format(index) for index in range(1000)]
        assignments = {sample: batch.shard_for_sample(sample, 4) for sample in samples}
        self.assertEqual(
            assignments,
            {sample: batch.shard_for_sample(sample, 4) for sample in samples},
        )
        selected = [
            {sample for sample in samples if assignments[sample] == shard}
            for shard in range(4)
        ]
        self.assertEqual(set(samples), set().union(*selected))
        self.assertTrue(
            all(
                not (left & right)
                for index, left in enumerate(selected)
                for right in selected[index + 1 :]
            )
        )

    def test_manifest_duplicate_selected_id_fails(self):
        fasta = self.inputs / "one.fna"
        write_fasta(fasta, "ONE")
        digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
        manifest = self.root / "duplicate.tsv"
        manifest.write_text(
            "\t".join(batch.MANIFEST_FIELDS)
            + "\n"
            + "one\t{}\t{}\t{}\n".format(fasta, fasta.stat().st_size, digest) * 2,
            encoding="utf-8",
        )
        with self.assertRaisesRegex(batch.InputError, "duplicate"):
            batch.index_manifest(manifest, self.root / "state.sqlite", 1, 0)


class PresetTests(BatchTestCase):
    def test_region_screen_prefers_public_detector_preset(self):
        parser = batch.build_parser(ROOT)
        args = parser.parse_args(
            [
                "run",
                "--manifest",
                "unused",
                "--run-dir",
                "unused",
                "--runner",
                str(self.runner),
            ]
        )
        with mock.patch.object(batch, "runtime_tree_sha256", return_value="a" * 64):
            config = batch.build_tool_configuration(args, ROOT)
        self.assertEqual("detector-preset", config.preset_mode)
        self.assertEqual(("--preset", "region-screen"), config.analysis_arguments)

    def test_region_screen_fallback_is_explicit_and_lean(self):
        parser = batch.build_parser(ROOT)
        args = parser.parse_args(
            [
                "run",
                "--manifest",
                "unused",
                "--run-dir",
                "unused",
                "--runner",
                str(self.runner),
            ]
        )
        with (
            mock.patch.object(batch, "_probe_preset_support", return_value=False),
            mock.patch.object(batch, "runtime_tree_sha256", return_value="a" * 64),
        ):
            config = batch.build_tool_configuration(args, ROOT)
        arguments = list(config.analysis_arguments)
        self.assertEqual("compatibility-expansion", config.preset_mode)
        for option in (
            "--output_profile",
            "--accepted_only",
            "--fuzzy_boundary_paths",
            "--fast_run_rescue",
            "--degenerated",
            "--strand",
            "--cas",
        ):
            self.assertIn(option, arguments)
        self.assertEqual("False", arguments[arguments.index("--degenerated") + 1])


class RestartAndFailureTests(BatchTestCase):
    def test_successful_reports_are_validated_and_safely_skipped(self):
        write_fasta(self.inputs / "a.fna", "A")
        write_fasta(self.inputs / "b.fna", "B")
        manifest = self.make_manifest()
        self.assertEqual(0, self.run_batch(manifest))
        summary = json.loads((self.root / "run" / "summary.json").read_text())
        self.assertEqual(2, summary["counts"]["succeeded"])
        self.assertEqual(0, summary["counts"]["failed"])
        self.assertGreater(summary["elapsed_seconds"], 0)
        self.assertGreater(summary["completed_genomes_per_hour"], 0)
        self.assertEqual(
            "fresh_succeeded_only_excludes_resume_skips",
            summary["throughput_basis"],
        )
        self.assertEqual(1.0, summary["completion_fraction"])
        self.assertEqual(1, int((self.counter_root / "a.fna.count").read_text()))
        self.assertEqual(0, self.run_batch(manifest))
        resumed = json.loads((self.root / "run" / "summary.json").read_text())
        self.assertEqual(2, resumed["counts"]["skipped"])
        self.assertEqual(0, resumed["counts"]["succeeded"])
        self.assertEqual(0, resumed["completed_genomes_per_hour"])
        self.assertEqual(1, int((self.counter_root / "a.fna.count").read_text()))

    def test_non_region_validation_allows_internal_candidate_categories(self):
        report = self.root / "report.json"
        sequence = "ACGT"
        payload = {
            "schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
            "source": {
                "id": "A",
                "length": len(sequence),
                "sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            },
            "summary": {"array_count": 1},
            "arrays": [{"id": "candidate", "category": "Alternative"}],
        }
        report.write_text(json.dumps(payload), encoding="utf-8")
        batch.validate_report(report, None, None, region_screen=False)
        with self.assertRaisesRegex(batch.CompletionError, "non-accepted"):
            batch.validate_report(report, None, None, region_screen=True)

    def test_region_screen_rejects_an_accepted_likely_tandem_survivor(self):
        sequence, array = _canonical_array("A" * 20, "A" * 20)
        payload = _current_report(sequence, [array])
        payload["extensions"] = {}
        report = self.root / "report.json"
        report.write_text(json.dumps(payload), encoding="utf-8")

        with self.assertRaisesRegex(
            batch.CompletionError, "still qualifies as likely tandem"
        ):
            batch.validate_report(
                report,
                {"A": sequence},
                None,
                region_screen=True,
            )

    def test_completion_replays_array_architecture_against_tables(self):
        sequence, array = _canonical_array("A" * 20, "C" * 20)
        array["spacers"][0]["sequence"] = "T" * 20
        payload = _current_report(sequence, [array])
        report = self.root / "report.json"
        report.write_text(json.dumps(payload), encoding="utf-8")

        with self.assertRaisesRegex(
            batch.CompletionError, "architecture/table differs from source"
        ):
            batch.validate_report(
                report,
                {"A": sequence},
                None,
                region_screen=True,
            )

    def test_region_screen_requires_current_likely_tandem_ledger(self):
        sequence = "ACGTACGT"
        payload = _current_report(sequence, [])
        payload["extensions"].pop("likely_tandem")
        report = self.root / "report.json"
        report.write_text(json.dumps(payload), encoding="utf-8")

        with self.assertRaisesRegex(
            batch.CompletionError, "lacks the likely-tandem audit ledger"
        ):
            batch.validate_report(
                report,
                {"A": sequence},
                None,
                region_screen=True,
            )

    def test_accepted_only_tandem_ledger_is_source_replayed(self):
        sequence, payload = _accepted_only_tandem_report()
        report = self.root / "report.json"
        report.write_text(json.dumps(payload), encoding="utf-8")
        batch.validate_report(
            report,
            {"A": sequence},
            None,
            region_screen=True,
        )

        forged_evidence = copy.deepcopy(payload)
        forged_evidence["extensions"]["likely_tandem"]["arrays"][0][
            "evidence"
        ]["minimum_pairwise_edit_identity"] = 0.0
        report.write_text(json.dumps(forged_evidence), encoding="utf-8")
        with self.assertRaisesRegex(
            batch.CompletionError, "evidence does not replay"
        ):
            batch.validate_report(
                report,
                {"A": sequence},
                None,
                region_screen=True,
            )

        for path, value in (
            (("evidence", "minimum_pairwise_edit_identity"), True),
            (("evidence", "maximum_pairwise_edit_distance"), False),
            (("source_architecture", "repeats", 0, "index"), True),
        ):
            with self.subTest(strict_json_path=path):
                forged_type = copy.deepcopy(payload)
                target = forged_type["extensions"]["likely_tandem"][
                    "arrays"
                ][0]
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                report.write_text(json.dumps(forged_type), encoding="utf-8")
                with self.assertRaises(batch.CompletionError):
                    batch.validate_report(
                        report,
                        {"A": sequence},
                        None,
                        region_screen=True,
                    )

        forged_source = copy.deepcopy(payload)
        ledger_spacer = forged_source["extensions"]["likely_tandem"][
            "arrays"
        ][0]["source_architecture"]["spacers"][1]
        ledger_spacer["sequence"] = "C" * 20
        ledger_spacer["sequence_sha256"] = hashlib.sha256(
            ledger_spacer["sequence"].encode("ascii")
        ).hexdigest()
        report.write_text(json.dumps(forged_source), encoding="utf-8")
        with self.assertRaisesRegex(
            batch.CompletionError, "ledger spacer differs from source"
        ):
            batch.validate_report(
                report,
                {"A": sequence},
                None,
                region_screen=True,
            )

    def test_full_tandem_ledger_requires_exact_provenance_and_json_types(self):
        sequence, payload = _full_tandem_report()
        report = self.root / "report.json"
        report.write_text(json.dumps(payload), encoding="utf-8")
        batch.validate_report(
            report,
            {"A": sequence},
            None,
            region_screen=False,
        )

        provenance_mutations = (
            (("method",), "forged-method"),
            (("transition", "rule_name"), "forged-rule"),
            (("transition", "rule_version"), "0.0.0"),
            (("transition", "source_category_name"), "Bona-fide"),
            (("transition", "target_category_index"), 3),
            (("transition", "reported_tag"), "forged-tag"),
            (("transition", "row_ordinal"), False),
        )
        for path, value in provenance_mutations:
            with self.subTest(provenance_path=path):
                forged = copy.deepcopy(payload)
                target = forged["arrays"][0]["classification_basis"]
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                report.write_text(json.dumps(forged), encoding="utf-8")
                with self.assertRaisesRegex(
                    batch.CompletionError, "ledger does not match"
                ):
                    batch.validate_report(
                        report,
                        {"A": sequence},
                        None,
                        region_screen=False,
                    )

        boolean_score = copy.deepcopy(payload)
        boolean_score["arrays"][0]["certainty_score"] = True
        report.write_text(json.dumps(boolean_score), encoding="utf-8")
        with self.assertRaisesRegex(batch.CompletionError, "ledger does not match"):
            batch.validate_report(
                report,
                {"A": sequence},
                None,
                region_screen=False,
            )

    def test_tampered_report_invalidates_completion_and_creates_new_attempt(self):
        write_fasta(self.inputs / "a.fna", "A")
        manifest = self.make_manifest()
        self.assertEqual(0, self.run_batch(manifest))
        report = next((self.root / "run" / "samples").rglob("report.json"))
        report.write_text("{}\n", encoding="utf-8")
        self.assertEqual(0, self.run_batch(manifest))
        self.assertEqual(2, int((self.counter_root / "a.fna.count").read_text()))
        item = next((self.root / "run" / "samples").iterdir())
        self.assertEqual(1, len(list((item / "superseded").glob("completion-*.json"))))
        self.assertEqual(2, len(list((item / "attempts").iterdir())))

    def test_retry_succeeds_on_second_isolated_attempt(self):
        write_fasta(self.inputs / "retry.fna", "RETRY_ONCE")
        manifest = self.make_manifest()
        self.assertEqual(0, self.run_batch(manifest, max_attempts=2))
        summary = json.loads((self.root / "run" / "summary.json").read_text())
        self.assertEqual(1, summary["counts"]["succeeded"])
        item = next((self.root / "run" / "samples").iterdir())
        attempts = sorted((item / "attempts").iterdir())
        self.assertEqual(2, len(attempts))
        statuses = [
            json.loads((path / "attempt.json").read_text())["status"]
            for path in attempts
        ]
        self.assertEqual(["exit_nonzero", "success"], statuses)

    def test_timeout_is_process_group_bounded_and_returns_nonzero(self):
        write_fasta(self.inputs / "sleep.fna", "SLEEP")
        manifest = self.make_manifest()
        started = time.monotonic()
        self.assertEqual(1, self.run_batch(manifest, timeout=0.15, max_attempts=1))
        self.assertLess(time.monotonic() - started, 3.0)
        summary = json.loads((self.root / "run" / "summary.json").read_text())
        self.assertEqual("failed", summary["status"])
        self.assertEqual(1, summary["counts"]["failed"])
        failures = (self.root / "run" / "failures.jsonl").read_text().splitlines()
        self.assertEqual(1, len(failures))
        self.assertEqual("timeout", json.loads(failures[0])["error_code"])
        child_pid = int((self.counter_root / "sleep-child.pid").read_text())
        process_status = Path("/proc") / str(child_pid) / "stat"
        deadline = time.monotonic() + 1
        while process_status.exists() and time.monotonic() < deadline:
            fields = process_status.read_text().split()
            if len(fields) > 2 and fields[2] == "Z":
                break
            time.sleep(0.01)
        if process_status.exists():
            self.assertEqual("Z", process_status.read_text().split()[2])

    def test_permanent_failure_has_deterministic_ledgers_and_nonzero_exit(self):
        write_fasta(self.inputs / "fail.fna", "FAIL_ALWAYS")
        manifest = self.make_manifest()
        self.assertEqual(1, self.run_batch(manifest, max_attempts=2))
        status = (self.root / "run" / "status.tsv").read_text().splitlines()
        self.assertEqual(2, len(status))
        failure = json.loads(
            (self.root / "run" / "failures.jsonl").read_text().splitlines()[0]
        )
        self.assertEqual("fail", failure["sample_id"])
        self.assertEqual(2, failure["attempts"])


class BoundedSchedulingTests(unittest.TestCase):
    def test_scheduler_never_consumes_more_than_fixed_window(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            state = root / "state.sqlite"
            connection = sqlite3.connect(str(state))
            connection.execute(
                """CREATE TABLE jobs (
                    ordinal INTEGER PRIMARY KEY, sample_id TEXT, fasta TEXT,
                    source_bytes INTEGER, source_sha256 TEXT, status TEXT,
                    attempts INTEGER, elapsed_seconds REAL, error_code TEXT,
                    error_message TEXT, completion TEXT)"""
            )
            for ordinal in range(1, 51):
                connection.execute(
                    "INSERT INTO jobs VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                    (
                        ordinal,
                        str(ordinal),
                        "unused",
                        1,
                        "a" * 64,
                        "pending",
                        0,
                        0.0,
                        "",
                        "",
                        "",
                    ),
                )
            connection.commit()
            connection.close()

            lock = threading.Lock()
            consumed = 0
            completed = 0
            maximum_outstanding = 0

            def jobs():
                nonlocal consumed, maximum_outstanding
                for ordinal in range(1, 51):
                    with lock:
                        consumed += 1
                        maximum_outstanding = max(
                            maximum_outstanding, consumed - completed
                        )
                    yield batch.Job(ordinal, str(ordinal), Path("unused"), 1, "a" * 64)

            def worker(job, _configuration):
                nonlocal completed
                time.sleep(0.002)
                with lock:
                    completed += 1
                return batch.JobResult(job.ordinal, job.sample_id, "succeeded", 1, 0.0)

            with mock.patch.object(batch, "execute_job", side_effect=worker):
                batch.run_bounded(jobs(), object(), 2, state, max_in_flight=4)
            self.assertLessEqual(maximum_outstanding, 4)
            self.assertEqual(50, completed)


if __name__ == "__main__":
    unittest.main()
