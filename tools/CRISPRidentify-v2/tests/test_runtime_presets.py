import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from components.runtime_presets import BALANCED_DEFAULTS
from components.runtime_presets import FULL_CONFIG
from components.runtime_presets import REGION_SCREEN_CONFIG
from components.runtime_presets import PresetConflictError
from components.runtime_presets import explicit_option_destinations
from components.runtime_presets import preset_provenance
from components.runtime_presets import resolve_preset


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "NC_017040.1.fasta"


class RuntimePresetUnitTests(unittest.TestCase):
    def test_balanced_preserves_backward_compatible_direct_flags(self):
        values = dict(BALANCED_DEFAULTS)
        values.update({
            "cas": True,
            "fast_run": False,
            "fuzzy_boundary_paths": False,
        })

        self.assertEqual(
            values,
            resolve_preset(
                "balanced",
                values,
                explicit_destinations={
                    "cas", "fast_run", "fuzzy_boundary_paths"
                },
            ),
        )

    def test_region_screen_resolves_the_validated_atomic_contract(self):
        resolved = resolve_preset("region-screen", dict(BALANCED_DEFAULTS))

        for destination, expected in REGION_SCREEN_CONFIG.items():
            self.assertEqual(expected, resolved[destination], destination)
        self.assertEqual("machine", resolved["output_profile"])
        self.assertTrue(resolved["accepted_only"])
        self.assertFalse(resolved["fuzzy_boundary_paths"])
        self.assertFalse(resolved["degenerated"])
        self.assertFalse(resolved["fast_run_rescue"])
        self.assertEqual("legacy", resolved["reconstruction_policy"])

    def test_architecture_rescue_keeps_legacy_unless_incompatibility_is_explicit(self):
        values = dict(BALANCED_DEFAULTS, architecture_rescue=True)
        resolved = resolve_preset("balanced", values,
                                  explicit_destinations={"architecture_rescue"})
        self.assertEqual("legacy", resolved["reconstruction_policy"])
        self.assertEqual("conservative", values["reconstruction_policy"])
        with self.assertRaisesRegex(PresetConflictError, "conservative reconstruction requires"):
            resolve_preset("balanced", values, explicit_destinations={
                "architecture_rescue", "reconstruction_policy"})

    def test_region_screen_rejects_explicit_scientific_conflicts(self):
        values = dict(BALANCED_DEFAULTS)
        values["cas"] = True

        with self.assertRaisesRegex(
                PresetConflictError, r"--cas=true .*requires false"):
            resolve_preset(
                "region-screen", values, explicit_destinations={"cas"}
            )

        with self.assertRaisesRegex(
                PresetConflictError, r"--fuzzy_boundary_paths=true"):
            resolve_preset(
                "region-screen",
                dict(BALANCED_DEFAULTS),
                explicit_destinations={"fuzzy_boundary_paths"},
            )

    def test_region_screen_accepts_matching_explicit_confirmation(self):
        values = dict(BALANCED_DEFAULTS)
        values["fuzzy_boundary_paths"] = False
        resolved = resolve_preset(
            "region-screen",
            values,
            explicit_destinations={"fuzzy_boundary_paths"},
        )
        self.assertFalse(resolved["fuzzy_boundary_paths"])

    def test_full_keeps_cas_and_is_annotations_independently_opt_in(self):
        self.assertNotIn("cas", FULL_CONFIG)
        self.assertNotIn("is_element", FULL_CONFIG)
        values = dict(BALANCED_DEFAULTS)
        values.update({"cas": True, "is_element": True})

        resolved = resolve_preset(
            "full",
            values,
            explicit_destinations={"cas", "is_element"},
        )

        self.assertFalse(resolved["fast_run"])
        self.assertFalse(resolved["fast_run_rescue"])
        self.assertFalse(resolved["fast_run_phase_repair"])
        self.assertFalse(resolved["split_competitive_components"])
        self.assertTrue(resolved["cas"])
        self.assertTrue(resolved["is_element"])

    def test_full_rejects_an_explicit_fast_search_conflict(self):
        with self.assertRaisesRegex(PresetConflictError, r"--fast_run=true"):
            resolve_preset(
                "full",
                dict(BALANCED_DEFAULTS),
                explicit_destinations={"fast_run"},
            )

    def test_explicit_option_detection_supports_both_cli_forms(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--fast_run")
        parser.add_argument("--cas")

        observed = explicit_option_destinations(
            parser, ("--fast_run=False", "--cas", "True")
        )

        self.assertEqual(frozenset(("fast_run", "cas")), observed)

    def test_provenance_is_deterministic_and_records_selection(self):
        config = dict(REGION_SCREEN_CONFIG)
        first = preset_provenance(
            "region-screen",
            config,
            explicit_destinations={"preset", "cas"},
            preset_explicit=True,
        )
        second = preset_provenance(
            "region-screen",
            dict(reversed(tuple(config.items()))),
            explicit_destinations={"cas", "preset"},
            preset_explicit=True,
        )

        self.assertEqual(first, second)
        self.assertEqual("explicit", first["selection"])
        self.assertEqual(["cas"], first["explicit_controlled_options"])
        self.assertEqual(64, len(first["effective_config_sha256"]))


class RuntimePresetCliTests(unittest.TestCase):
    @staticmethod
    def _run_cli(*arguments):
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = "0"
        environment["CRISPRIDENTIFY_ENV"] = sys.prefix
        return subprocess.run(
            [str(ROOT / "crispridentify")] + list(arguments),
            cwd=str(ROOT),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=False,
        )

    def test_print_effective_region_screen_config_needs_no_input(self):
        completed = self._run_cli(
            "--preset", "region-screen", "--print-effective-config"
        )

        self.assertEqual(0, completed.returncode, completed.stderr)
        document = json.loads(completed.stdout)
        self.assertEqual("region-screen", document["preset"])
        effective = document["effective_config"]
        for destination, expected in REGION_SCREEN_CONFIG.items():
            self.assertEqual(expected, effective[destination], destination)

    def test_default_reconstruction_matches_explicit_conservative(self):
        default = self._run_cli("--print-effective-config")
        enabled = self._run_cli(
            "--reconstruction-policy=conservative", "--print-effective-config"
        )
        self.assertEqual(0, default.returncode, default.stderr)
        self.assertEqual(0, enabled.returncode, enabled.stderr)
        default_config = json.loads(default.stdout)
        enabled_config = json.loads(enabled.stdout)
        self.assertEqual(default_config["effective_config"], enabled_config["effective_config"])
        self.assertEqual("conservative", default_config["effective_config"]["reconstruction_policy"])
        self.assertEqual(default_config["preset_provenance"]["effective_config_sha256"],
                         enabled_config["preset_provenance"]["effective_config_sha256"])
        self.assertIn("reconstruction_policy",
                      enabled_config["preset_provenance"]["explicit_controlled_options"])

    def test_legacy_reconstruction_remains_an_explicit_opt_out(self):
        legacy = self._run_cli(
            "--reconstruction_policy", "legacy", "--print-effective-config"
        )
        self.assertEqual(0, legacy.returncode, legacy.stderr)
        self.assertEqual("legacy", json.loads(legacy.stdout)["effective_config"]["reconstruction_policy"])

    def test_other_workflows_keep_legacy_and_reject_explicit_conservative(self):
        for extra in (("--preset", "full"), ("--preset", "region-screen"),
                      ("--architecture_rescue=True",)):
            with self.subTest(extra=extra):
                implicit = self._run_cli(*extra, "--print-effective-config")
                self.assertEqual(0, implicit.returncode, implicit.stderr)
                self.assertEqual("legacy", json.loads(implicit.stdout)["effective_config"]["reconstruction_policy"])
                explicit_legacy = self._run_cli(
                    *extra, "--reconstruction-policy=legacy", "--print-effective-config")
                self.assertEqual(0, explicit_legacy.returncode, explicit_legacy.stderr)
                rejected = self._run_cli(
                    "--reconstruction_policy", "conservative",
                    *extra, "--print-effective-config"
                )
                self.assertEqual(2, rejected.returncode)
                self.assertIn("reconstruction", rejected.stderr)

    def test_cli_reports_atomic_preset_conflict(self):
        completed = self._run_cli(
            "--preset", "region-screen",
            "--fuzzy_boundary_paths", "True",
            "--print-effective-config",
        )

        self.assertEqual(2, completed.returncode)
        self.assertIn(
            "explicit option conflicts with --preset region-screen",
            completed.stderr,
        )
        self.assertIn("--fuzzy_boundary_paths=true", completed.stderr)

    def test_default_runs_reconstruction_and_legacy_removes_its_sidecar(self):
        with tempfile.TemporaryDirectory() as directory:
            common = ("--file", str(EXAMPLE), "--result_folder", directory,
                      "--cpu", "1", "--parallel", "False", "--quiet")
            default = self._run_cli(*common)
            self.assertEqual(0, default.returncode, default.stderr)
            reports = list(Path(directory).rglob("report.json"))
            self.assertEqual(1, len(reports))
            report = json.loads(reports[0].read_text())
            self.assertEqual("conservative", report["run"]["parameters"]["param_reconstruction_policy"])
            evidence = reports[0].with_name("reconstruction_evidence.json")
            audit = json.loads(evidence.read_text())
            self.assertFalse(audit["experimental"])
            self.assertFalse(audit["biological_accuracy_established"])

            legacy = self._run_cli(*common, "--reconstruction_policy", "legacy")
            self.assertEqual(0, legacy.returncode, legacy.stderr)
            self.assertFalse(evidence.exists())
            legacy_report = json.loads(reports[0].read_text())
            self.assertEqual("legacy", legacy_report["run"]["parameters"]["param_reconstruction_policy"])

    def test_quiet_and_verbose_are_mutually_exclusive(self):
        completed = self._run_cli(
            "--quiet", "--verbose", "--print-effective-config"
        )

        self.assertEqual(2, completed.returncode)
        self.assertIn("not allowed with argument --quiet", completed.stderr)

    def test_default_is_concise_and_verbose_restores_pipeline_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            common = (
                "--file", str(EXAMPLE),
                "--preset", "region-screen",
                "--cpu", "1",
            )
            concise = self._run_cli(
                *common, "--result_folder", str(Path(temporary) / "concise")
            )
            verbose = self._run_cli(
                *common, "--result_folder", str(Path(temporary) / "verbose"),
                "--verbose",
            )

        self.assertEqual(0, concise.returncode, concise.stderr)
        self.assertIn("Processing 1/1", concise.stdout)
        self.assertIn("Finished in", concise.stdout)
        self.assertNotIn("1. Run initial array detection", concise.stdout)
        self.assertEqual(0, verbose.returncode, verbose.stderr)
        self.assertIn("1. Run initial array detection", verbose.stdout)

    def test_quiet_region_screen_is_silent_and_records_configuration(self):
        with tempfile.TemporaryDirectory() as temporary:
            completed = self._run_cli(
                "--file", str(EXAMPLE),
                "--result_folder", temporary,
                "--preset", "region-screen",
                "--quiet",
                "--cpu", "1",
            )

            self.assertEqual(0, completed.returncode, completed.stderr)
            self.assertEqual("", completed.stdout)
            reports = list(Path(temporary).rglob("report.json"))
            self.assertEqual(1, len(reports))
            report = json.loads(reports[0].read_text(encoding="utf-8"))
            configuration = report["run"]["configuration"]
            self.assertEqual("region-screen", configuration["preset"])
            effective = configuration["effective_config"]
            self.assertFalse(effective["fuzzy_boundary_paths"])
            self.assertFalse(effective["degenerated"])
            self.assertFalse(effective["strand"])
            self.assertFalse(effective["cas"])
            self.assertEqual(
                False,
                report["run"]["parameters"]["param_fuzzy_boundary_paths"],
            )


if __name__ == "__main__":
    unittest.main()
