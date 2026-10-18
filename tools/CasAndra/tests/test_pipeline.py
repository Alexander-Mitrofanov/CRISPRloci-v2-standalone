from types import SimpleNamespace

from casandra.cli import build_parser, main


def test_cli_exposes_cas_only_genome_protein_and_cassette_modes():
    parser = build_parser()
    genome = parser.parse_args(["predict-genome", "--model", "m", "--genome", "g.fa", "--output", "out"])
    proteins = parser.parse_args(["predict-proteins", "--model", "m", "--input", "p.fa", "--output", "out.jsonl"])
    annotation = parser.parse_args(["annotate-proteins", "--model", "m", "--input", "p.fa", "--output", "out"])
    cassette = parser.parse_args(["classify-cassette", "--model", "m", "--input", "p.fa", "--output", "out"])
    assert genome.command == "predict-genome"
    assert proteins.command == "predict-proteins"
    assert annotation.command == "annotate-proteins"
    assert cassette.command == "classify-cassette"


def test_cli_dispatches_atomic_protein_annotation(monkeypatch, capsys):
    calls = []

    class FakeAnalyzer:
        def annotate_proteins(self, input_path, output_path, *, threads):
            calls.append((input_path, output_path, threads))
            return SimpleNamespace(summary={"analysis": "annotate_cas_genes"})

    monkeypatch.setattr(
        "casandra.cli.Analyzer.from_bundle", lambda _model: FakeAnalyzer()
    )

    assert (
        main(
            [
                "annotate-proteins",
                "--model",
                "bundle",
                "--input",
                "proteins.faa",
                "--output",
                "annotation",
                "--threads",
                "3",
            ]
        )
        == 0
    )
    assert [(str(source), str(output), threads) for source, output, threads in calls] == [
        ("proteins.faa", "annotation", 3)
    ]
    assert '"analysis": "annotate_cas_genes"' in capsys.readouterr().out
