from dataclasses import replace
from pathlib import Path

import pytest

from crisprtracrrna.config import Settings
from crisprtracrrna.errors import ContractError, TracrError
from crisprtracrrna.evidence import short_rescue_admission
from crisprtracrrna.models import AntiRepeat, Array, Candidate, Hit, Interval, Sequence
from crisprtracrrna.search import interactions, parse_alignment_table, parse_intarna_batch
from crisprtracrrna.sequence import read_fasta

HEADER = "id1;id2;start1;end1;start2;end2;E;hybridDP\n"
ROW = "t;q;1;4;1;4;-4;((((&))))\n"


def candidate(identity="c", repeat="T" * 20, strand="+"):
    array = Array("array", Interval(250, 300), repeat, "+", "Possible", 0.8)
    anti = AntiRepeat(array, repeat, "+", Interval(100, 120), strand, 0.9, 0.5, True)
    return Candidate(identity, anti, Interval(95, 180), "A" * 85, Interval(5, 25))


def parse(tmp_path, rows):
    path = tmp_path / "out.csv"
    path.write_text(HEADER + rows)
    return parse_intarna_batch(path, {"t": 100, "absent": 100}, "q", 20)


def test_all_rows_validated_and_order_independent(tmp_path):
    other = "t;q;2;5;2;5;-4;((((&))))\n"
    assert parse(tmp_path, ROW + other) == parse(tmp_path, other + ROW)
    assert parse(tmp_path, ROW)["t"]["target_paired_positions"] == [1, 2, 3, 4]
    assert parse(tmp_path, ROW)["absent"] is None
    assert parse(tmp_path, "") == {"t": None, "absent": None}


@pytest.mark.parametrize(
    "bad",
    [
        "x;q;1;4;1;4;-3;((((&))))\n",
        "t;x;1;4;1;4;-3;((((&))))\n",
        "t;q;1;400;1;4;-3;((((&))))\n",
        "t;q;4;1;1;4;-3;((((&))))\n",
        "t;q;0;3;1;4;-3;((((&))))\n",
        "t;q;1;4;1;4;nan;((((&))))\n",
        "t;q;1;4;1;4;-3;(((.&))))\n",
        "t;q;1;4;1;4;-3;....&....\n",
        "t;q;1;4;1;4;-3;(((x&))))\n",
        "t;q;1;4;1;4;-3;(((&)))\n",
        "t;q;1;4;1;4;-3;((((&&))))\n",
        "t;q;1;4;1;4;-3;\n",
        ROW,
        "t;q;1;4;1;4;-3;((((&))))\n",
    ],
)
def test_bad_nonwinning_rows_fail(tmp_path, bad):
    with pytest.raises(ContractError):
        parse(tmp_path, ROW + bad)


def test_duplicate_columns_and_truncated_rows_fail(tmp_path):
    path = tmp_path / "out.csv"
    for value in (HEADER.strip() + ";E\n", HEADER + "t;q;1\n"):
        path.write_text(value)
        with pytest.raises(ContractError):
            parse_intarna_batch(path, {"t": 100}, "q", 20)


class PairRunner:
    def __init__(self):
        self.pairs = []
        self.calls = []

    def run(self, args):
        self.calls.append(args)
        targets = list(read_fasta(Path(args[args.index("-t") + 1])))
        queries = list(read_fasta(Path(args[args.index("-q") + 1])))
        assert len(queries) == 1
        rows = []
        for target in reversed(targets):
            self.pairs.append((target.bases, queries[0].bases))
            rows.append(f"{target.id};{queries[0].id};1;4;1;4;-4;((((&))))\n")
        Path(args[args.index("--out") + 1]).write_text(HEADER + "".join(rows))


def test_batch_routes_exact_pairs_chunks_and_empty(tmp_path):
    inputs = [replace(candidate(f"c{i}"), sequence="A" * (85 + i)) for i in range(5)]
    inputs += [candidate("reverse", "A" * 20), candidate("duplicate_sequence")]
    runner = PairRunner()
    result = interactions(inputs[::-1], tmp_path, Settings(intarna_batch_size=2), runner)
    assert set(result) == {c.id for c in inputs}
    assert sorted(runner.pairs) == sorted((c.sequence, c.anti.repeat) for c in inputs)
    assert len(runner.calls) == 4
    assert all(len(list(read_fasta(Path(c[c.index("-t") + 1])))) <= 2 for c in runner.calls)
    assert interactions([], tmp_path, Settings(), runner) == {}
    with pytest.raises(ContractError):
        interactions([inputs[0], inputs[0]], tmp_path, Settings(), runner)
    assert len(runner.calls) == 4


def test_anchored_calls_use_single_target_and_own_seed_range(tmp_path):
    runner = PairRunner()
    inputs = [candidate("a"), replace(candidate("b"), anti_local=Interval(10, 30))]
    interactions(inputs, tmp_path, Settings(), runner, anchored=True)
    assert [c[c.index("--seedTRange") + 1] for c in runner.calls] == ["6-25", "11-30"]
    assert all(len(list(read_fasta(Path(c[c.index("-t") + 1])))) == 1 for c in runner.calls)


@pytest.mark.parametrize("value", [True, False, 0, -1, 1.5, "128"])
def test_batch_size_is_strict_positive_integer(value):
    with pytest.raises(TracrError):
        Settings(intarna_batch_size=value)


def test_experimental_scope_and_boolean():
    for options in ({"model_type": "V"}, {"mode": "model"}, {"short_anti_repeat_rescue": 1}):
        with pytest.raises(TracrError):
            Settings(**{"short_anti_repeat_rescue": True, **options})


def test_short_proposals_are_native_short_complementary_and_array_external(tmp_path):
    c = candidate()
    source = Sequence("s", "s", "A" * 500)
    repeat = "T" * 40
    queries = {"q": (c.anti.array, repeat, "+")}
    path = tmp_path / "align.tsv"
    path.write_text(
        "q s 80 12 0 0 1 12 30 41 1 5\n"  # accepted equality
        "q s 100 11 0 0 1 11 50 60 1 5\n"  # too short
        "q s 79 12 0 0 1 12 70 81 1 5\n"  # low identity
        "q s 100 24 0 0 1 24 100 123 1 5\n"  # standard boundary
        "q s 100 12 0 0 1 12 249 260 1 5\n"  # array overlap
    )
    hits = parse_alignment_table(path, queries, source, [c.anti.array], Settings(), short=True)
    assert len(hits) == 1 and hits[0].complementary and hits[0].strand == "-"
    assert hits[0].coverage == 0.3


def duplex(positions=None):
    return {
        "target_interval": Interval(5, 25).public(),
        "energy_kcal_mol": -4,
        "target_paired_positions": list(range(6, 26)) if positions is None else positions,
        "repeat_paired_positions": list(range(1, 21))
        if positions is None
        else list(range(1, len(positions) + 1)),
    }


@pytest.mark.parametrize("strand", ["+", "-"])
def test_short_gate_pairs_anchored_support_thresholds_and_strands(strand):
    c = candidate(strand=strand)
    model = Hit("c", Interval(25, 60), "+", 5, 0.01)
    term = Hit("c", Interval(25, 50), "+", 5, 0.001)
    accepted = short_rescue_admission(c, [model], [], duplex(), Settings())
    assert accepted["accepted"] and accepted["anchored_paired_count"] == 20
    assert accepted["antiparallel_pairs"][0] == {"target": 6, "repeat": 20, "bases": "AU"}
    assert short_rescue_admission(c, [], [term], duplex(), Settings())["accepted"]
    assert not short_rescue_admission(c, [replace(model, strand="-")], [], duplex(), Settings())[
        "accepted"
    ]
    assert not short_rescue_admission(
        c, [replace(model, evalue=0.01001)], [], duplex(), Settings()
    )["accepted"]
    assert not short_rescue_admission(
        c, [], [replace(term, interval=Interval(24, 50))], duplex(), Settings()
    )["accepted"]
    assert not short_rescue_admission(
        c, [replace(model, interval=Interval(15, 60))], [], duplex(), Settings()
    )["accepted"]
    assert not short_rescue_admission(c, [model], [], None, Settings())["accepted"]
    bad = replace(c, sequence="C" * 85)
    assert (
        "noncanonical_or_ambiguous_pair"
        in short_rescue_admission(bad, [model], [], duplex(), Settings())["rejection_reasons"]
    )
    # An overlapping bounding span cannot stand in for actual paired anchor bases.
    unanchored = duplex(list(range(26, 46)))
    assert not short_rescue_admission(c, [model], [], unanchored, Settings())["accepted"]
    assert short_rescue_admission(c, [model], [], duplex(list(range(6, 16))), Settings())[
        "accepted"
    ]
    assert not short_rescue_admission(c, [model], [], duplex(list(range(6, 15))), Settings())[
        "accepted"
    ]


def test_python_module_provenance_distinguishes_interpreter_checkout(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from crisprtracrrna.pipeline import command_provenance

    interpreter = tmp_path / "interpreter" / "python"
    interpreter.parent.mkdir()
    interpreter.write_text("python")
    package = tmp_path / "actual" / "casandra" / "__init__.py"
    package.parent.mkdir(parents=True)
    package.write_text("raise RuntimeError('must never import')")

    def fake_run(args, **kwargs):
        return SimpleNamespace(
            returncode=0,
            stdout="module-revision\n" if "actual" in args[2] else "interpreter-revision\n",
        )

    monkeypatch.setattr("crisprtracrrna.pipeline.subprocess.run", fake_run)
    result = command_provenance(
        (str(interpreter), "-m", "casandra"), {"PYTHONPATH": str(tmp_path / "actual")}
    )
    assert result["git_commit"] == "module-revision"
    assert result["interpreter_checkout_git_commit"] == "interpreter-revision"
    assert result["module_source"] == str(package)
    unresolved = command_provenance(
        (str(interpreter), "-m", "missing"), {"PYTHONPATH": str(tmp_path / "actual")}
    )
    assert "git_commit" not in unresolved and "unresolved" in unresolved["module_resolution"]


def test_harmless_empty_extra_native_column(tmp_path):
    path = tmp_path / "extra.csv"
    path.write_text(HEADER.strip() + ";future\n" + ROW.strip() + ";\n")
    assert parse_intarna_batch(path, {"t": 100}, "q", 20)["t"]["energy_kcal_mol"] == -4


def test_batch_contract_failure_publishes_failure_not_complete(tmp_path, monkeypatch):
    import crisprtracrrna.pipeline as pipeline

    source = tmp_path / "input.fasta"
    source.write_text(">s\n" + "A" * 500 + "\n")
    c = candidate()
    monkeypatch.setattr(pipeline, "doctor", lambda settings: {"ready": True})
    monkeypatch.setattr(pipeline, "run_identify", lambda *args: ([c.anti.array], {}))
    monkeypatch.setattr(pipeline, "run_casandra", lambda *args: {"cassettes": []})
    monkeypatch.setattr(pipeline, "search_anti_repeats", lambda *args: [c.anti])
    monkeypatch.setattr(pipeline, "scan_models", lambda *args: [])
    monkeypatch.setattr(pipeline, "terminators", lambda *args: {})

    def fail(*args):
        raise ContractError("malformed native batch")

    monkeypatch.setattr(pipeline, "interactions", fail)
    with pytest.raises(ContractError):
        pipeline.run(source, tmp_path / "output", Settings())
    assert (tmp_path / "output/failure.json").exists()
    assert not (tmp_path / "output/report.json").exists()


def test_legacy_single_interactions_preserve_every_calls_evidence(tmp_path):
    from crisprtracrrna.search import interaction

    runner = PairRunner()
    for identity in ("first", "second", "first", "../unsafe"):
        interaction(candidate(identity), tmp_path, Settings(), runner)
    assert len(list(tmp_path.glob("intarna-single-*/*.csv"))) == 4
    assert len(list(tmp_path.glob("intarna-single-*/*.ids.json"))) == 4
    assert len({str(Path(c[c.index("--out") + 1]).parent) for c in runner.calls}) == 4
