import dataclasses
import tempfile
import unittest
from pathlib import Path

from components.cli.parser import build_parser
from components.detection.repeat_templates import (
    DEFAULT_POLICY, TemplateDiscoveryError, load_repeat_template_panel,
    panel_manifest, parse_repeat_template_panel, reverse_complement,
    snapshot_repeat_template_panel, validate_template_workflow,
)
from components.runtime.folder_scheduler import rewrite_lane_argv

MOTIF = "GTTTTTATCGTACCTATGAGGAATTGAAAC"


class TemplateLibraryTests(unittest.TestCase):
    def test_normalization_preserves_aliases_but_scans_unique_orientations(self):
        first = parse_repeat_template_panel((">z\n" + MOTIF[:10].lower() + "\n" + MOTIF[10:] +
                                            "\n>a\n" + reverse_complement(MOTIF) + "\n").encode())
        second = parse_repeat_template_panel((">a\n" + reverse_complement(MOTIF).lower() +
                                             "\n>z\n" + MOTIF + "\n").encode())
        self.assertEqual(first.normalized_sha256, second.normalized_sha256)
        self.assertNotEqual(first.raw_sha256, second.raw_sha256)
        self.assertEqual(2, len(first.patterns))
        self.assertEqual(1, len(first.canonical_sequences))
        self.assertEqual(["a", "z"], [a.id for a in first.aliases])
        palindrome = parse_repeat_template_panel(b">pal\nACGTACGTACGTACGT\n")
        self.assertEqual(1, len(palindrome.patterns))

    def test_rejects_invalid_panels(self):
        for raw in (b"", b">\nACGT\n", b">id\n", b"ACGT\n", b">id\n" + b"A" * 14,
                    b">id\n" + b"A" * 64, b">id\n" + b"A" * 15 + b" N", b">id\n" + b"A" * 15 + b"-",
                    b">id\n" + b"A" * 15 + b"\xff", b">id\n" + b"A" * 15 + b"\n>id\n" + b"C" * 15,
                    b">id\nAAAAAAAA AAAAAAAA\n", b">id\nAAAAAAAAAAAAAAAN\n"):
            with self.subTest(raw=raw), self.assertRaises(TemplateDiscoveryError):
                parse_repeat_template_panel(raw)

    def test_panel_limits_at_boundaries(self):
        raw = (">one\n" + MOTIF + "\n").encode()
        for count in (len(raw) - 1, len(raw), len(raw) + 1):
            policy = dataclasses.replace(DEFAULT_POLICY, raw_bytes=count)
            if count < len(raw):
                with self.assertRaisesRegex(TemplateDiscoveryError, "raw_bytes"):
                    parse_repeat_template_panel(raw, policy)
            else:
                parse_repeat_template_panel(raw, policy)
        for field, sample, count in (
                ("records", ">a\n" + MOTIF + "\n>b\n" + MOTIF + "\n", 2),
                ("canonical_motifs", ">a\n" + MOTIF + "\n>b\n" + "A" * 30 + "\n", 2),
                ("oriented_patterns", ">a\n" + MOTIF + "\n", 2)):
            for cap in (count - 1, count, count + 1):
                policy = dataclasses.replace(DEFAULT_POLICY, **{field: cap})
                if cap < count:
                    with self.assertRaisesRegex(TemplateDiscoveryError, field):
                        parse_repeat_template_panel(sample.encode(), policy)
                else:
                    parse_repeat_template_panel(sample.encode(), policy)

    def test_snapshot_owns_bytes_and_refuses_changed_existing_file(self):
        with tempfile.TemporaryDirectory() as folder:
            original, snapshot = Path(folder) / "source panel.fa", Path(folder) / "frozen panel.fa"
            original.write_text(">one\n" + MOTIF + "\n")
            panel = load_repeat_template_panel(original)
            snapshot_repeat_template_panel(panel, snapshot)
            original.write_text(">changed\n" + "A" * 30 + "\n")
            self.assertEqual(panel_manifest(panel), panel_manifest(load_repeat_template_panel(snapshot)))
            snapshot.write_text("tampered")
            with self.assertRaises(TemplateDiscoveryError):
                snapshot_repeat_template_panel(panel, snapshot)

    def test_cli_and_folder_aliases_and_abbreviations(self):
        for alias in ("--repeat-template-fasta", "--repeat_template_fasta", "--repeat-t", "--repeat_t"):
            for args in ([alias, "original panel.fa"], [alias + "=original panel.fa"]):
                parsed = build_parser().parse_args(args)
                self.assertEqual("original panel.fa", parsed.repeat_template_fasta)
                rewritten = rewrite_lane_argv(args, "/inputs", "/outputs", 1,
                    repeat_template_snapshot="/tmp/frozen panel.fa", repeat_template_sha256="a" * 64)
                self.assertNotIn("original panel.fa", rewritten)
                self.assertIn("/tmp/frozen panel.fa", rewritten)
                self.assertEqual(1, rewritten.count("--repeat-template-fasta"))

    def test_unsupported_workflows_fail(self):
        for preset, accepted, engine in (("region-screen", False, "native"),
                                         ("balanced", True, "native"), ("full", False, "family")):
            with self.assertRaises(TemplateDiscoveryError):
                validate_template_workflow(preset, accepted, engine)
