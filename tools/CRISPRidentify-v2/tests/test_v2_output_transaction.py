from pathlib import Path
import tempfile
import unittest
from unittest import mock

from tests.test_v2_casandra_reporting import _casandra_result
from tests.test_v2_casandra_reporting import cas_output


publication_module = cas_output._output_publication


def _hidden_transaction_files(directory):
    return sorted(
        path.name for path in Path(directory).iterdir()
        if path.name.startswith((
            ".crispridentify-stage-",
            ".crispridentify-backup-",
        ))
    )


class CasPairRollbackTests(unittest.TestCase):
    def test_second_csv_replace_failure_restores_the_old_pair(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            systems = root / "Cas_systems.csv"
            genes = root / "Cas_genes.csv"
            systems.write_text("old systems\n", encoding="utf-8")
            genes.write_text("old genes\n", encoding="utf-8")
            real_replace = publication_module.os.replace

            def fail_second_csv(source, target):
                source = Path(source)
                target = Path(target)
                if (
                        source.name.startswith(".crispridentify-stage-")
                        and target.name == "Cas_genes.csv"):
                    raise OSError("injected second CSV replacement failure")
                return real_replace(source, target)

            with mock.patch.object(
                    publication_module.os,
                    "replace",
                    side_effect=fail_second_csv):
                with self.assertRaisesRegex(OSError, "injected"):
                    cas_output.write_cas_outputs(
                        root,
                        _casandra_result(),
                        source_id="record",
                        source_length=1_000,
                    )

            self.assertEqual(
                "old systems\n", systems.read_text(encoding="utf-8")
            )
            self.assertEqual(
                "old genes\n", genes.read_text(encoding="utf-8")
            )
            self.assertEqual([], _hidden_transaction_files(root))

    def test_duplicate_stage_target_does_not_leak_second_stage_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with cas_output.OutputTransaction() as publication:
                publication.stage_text(root / "artifact.txt", "first")
                self.assertEqual(1, len(_hidden_transaction_files(root)))
                with self.assertRaisesRegex(ValueError, "more than once"):
                    publication.stage_text(root / "artifact.txt", "second")
                self.assertEqual(1, len(_hidden_transaction_files(root)))
            self.assertEqual([], _hidden_transaction_files(root))


class CompletionMarkerTests(unittest.TestCase):
    def _stage_complete_run(self, publication, root):
        cas_output.stage_cas_outputs(
            publication,
            root,
            _casandra_result(),
            source_id="record",
            source_length=1_000,
        )
        publication.stage_text(root / "report.html", "new html\n")
        publication.stage_text(root / "report.json", '{"complete":true}\n')

    def test_failure_rolls_back_cas_and_html_without_completion_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            systems = root / "Cas_systems.csv"
            genes = root / "Cas_genes.csv"
            html = root / "report.html"
            completion = root / "report.json"
            systems.write_text("old systems\n", encoding="utf-8")
            genes.write_text("old genes\n", encoding="utf-8")
            html.write_text("old html\n", encoding="utf-8")
            completion.write_text("stale success\n", encoding="utf-8")

            publication_module.remove_completion_marker(completion)
            self.assertFalse(completion.exists())
            real_replace = publication_module.os.replace

            def fail_second_csv(source, target):
                source = Path(source)
                target = Path(target)
                if (
                        source.name.startswith(".crispridentify-stage-")
                        and target.name == "Cas_genes.csv"):
                    raise OSError("injected shared publication failure")
                return real_replace(source, target)

            with cas_output.OutputTransaction() as publication:
                self._stage_complete_run(publication, root)
                with mock.patch.object(
                        publication_module.os,
                        "replace",
                        side_effect=fail_second_csv):
                    with self.assertRaisesRegex(OSError, "injected"):
                        publication.publish(completion_target=completion)

            self.assertFalse(completion.exists())
            self.assertEqual(
                "old systems\n", systems.read_text(encoding="utf-8")
            )
            self.assertEqual(
                "old genes\n", genes.read_text(encoding="utf-8")
            )
            self.assertEqual("old html\n", html.read_text(encoding="utf-8"))
            self.assertEqual([], _hidden_transaction_files(root))

    def test_successful_publication_replaces_report_json_strictly_last(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            completion = root / "report.json"
            stage_targets = []
            real_replace = publication_module.os.replace

            def record_replace(source, target):
                source = Path(source)
                target = Path(target)
                if source.name.startswith(".crispridentify-stage-"):
                    stage_targets.append(target.name)
                return real_replace(source, target)

            with cas_output.OutputTransaction() as publication:
                self._stage_complete_run(publication, root)
                with mock.patch.object(
                        publication_module.os,
                        "replace",
                        side_effect=record_replace):
                    publication.publish(completion_target=completion)

            self.assertEqual("report.json", stage_targets[-1])
            self.assertEqual(
                ["Cas_systems.csv", "Cas_genes.csv", "report.html"],
                stage_targets[:-1],
            )
            self.assertEqual(
                '{"complete":true}\n',
                completion.read_text(encoding="utf-8"),
            )
            self.assertEqual([], _hidden_transaction_files(root))


if __name__ == "__main__":
    unittest.main()

