import hashlib
import json
from pathlib import Path

import pytest
from cctk.comparison import compare, run_compare, plot_comparison, table_cell
from cctk.errors import CCTKError
from cctk.io import safe_open, read_arrays, read_fasta, load_json


def document():
    return json.loads((Path(__file__).parents[1] / "examples/arrays.json").read_text())


def test_occurrences_gaps_threshold_and_isolates():
    result = compare(document())
    assert result["edges"][0]["nshared"] == 1
    assert result["edges"][0]["jaccard"] == pytest.approx(1 / 3)
    assert len(result["arrays"][0]["spacers"]) == 4
    assert (
        result["arrays"][0]["spacers"][1]["symbol"] == result["arrays"][0]["spacers"][2]["symbol"]
    )
    assert result["arrays"][0]["spacers"][-1]["symbol"] is None
    assert result["arrays"][0]["strand"] is None
    assert result["isolates"] == ["a2"]
    assert compare(document(), min_shared=2)["isolates"] == ["a0", "a1", "a2"]


def test_reverse_complement_explicit_not_strand_inference():
    doc = document()
    doc["arrays"] = doc["arrays"][:2]
    doc["arrays"][0]["spacers"] = doc["arrays"][0]["spacers"][:1]
    doc["arrays"][1]["spacers"] = doc["arrays"][1]["spacers"][:1]
    doc["arrays"][1]["spacers"][0]["sequence"] = "TGCAACGT"
    doc["arrays"][1]["strand"] = "-"
    assert compare(doc)["edges"] == []
    assert compare(doc, identity="reverse-complement")["edges"][0]["jaccard"] == 1
    assert compare(doc, identity="reverse-complement")["arrays"][1]["strand"] == "-"


@pytest.mark.parametrize(
    "change",
    [
        lambda d: d["arrays"][0]["spacers"][0].update(sequence="NNNN"),
        lambda d: d["arrays"][0]["spacers"][0].update(sequence=None),
        lambda d: d["arrays"][0]["spacers"][1].update(id="occ0"),
        lambda d: d["arrays"][0]["spacers"][1].update(ordinal=4),
        lambda d: d["arrays"][0].update(strand="unknown"),
        lambda d: d["arrays"][0].update(untrusted="field"),
        lambda d: d["arrays"][1].update(source_id="source0"),
    ],
)
def test_malformed_array_documents(change):
    doc = document()
    change(doc)
    with pytest.raises(CCTKError):
        compare(doc)


def test_singleton_and_all_deleted_no_artificial_edge():
    doc = document()
    doc["arrays"] = doc["arrays"][:1]
    assert compare(doc)["isolates"] == ["a0"]
    doc["arrays"][0]["spacers"] = [
        {"id": "gap", "ordinal": 1, "sequence": None, "status": "deletion"}
    ]
    assert compare(doc)["edges"] == []


def test_transaction_and_digest(tmp_path):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(document()))
    result = run_compare(
        path, tmp_path / "out", identity="exact", min_shared=1, order="supplied", plot=True
    )
    assert result["input_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert (tmp_path / "out/comparison.png").is_file()
    with pytest.raises(CCTKError, match="new"):
        run_compare(path, tmp_path / "out", identity="exact", min_shared=1, order="supplied")
    assert not list(tmp_path.glob(".cctk-*"))


def test_deleted_marker_not_dna_rectangle(tmp_path, monkeypatch):
    from cctk import plotting
    import matplotlib.axes

    rectangles = []
    labels = []
    original = plotting.plot_spacer
    original_text = matplotlib.axes.Axes.text

    def spacer(*args, **kwargs):
        rectangles.append(args)
        return original(*args, **kwargs)

    def text(self, *args, **kwargs):
        labels.append(args[2])
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(plotting, "plot_spacer", spacer)
    monkeypatch.setattr(matplotlib.axes.Axes, "text", text)
    result = compare(document())
    plot_comparison(result, tmp_path / "plot.png")
    assert len(rectangles) == 6
    assert "×" in labels


def test_bounded_symlink_and_duplicate_inputs(tmp_path):
    path = tmp_path / "arrays"
    path.write_text("A 1\nA 2\n")
    with pytest.raises(CCTKError):
        read_arrays(path)
    alias = tmp_path / "alias"
    alias.symlink_to(path)
    with pytest.raises(CCTKError):
        safe_open(alias)
    path.write_text('{"arrays":[],"arrays":[]}')
    with pytest.raises(CCTKError):
        load_json(path)
    path.write_text(">same\nACGT\n>same\nACGT\n")
    with pytest.raises(CCTKError):
        read_fasta(path)


def test_append_budget(tmp_path, monkeypatch):
    import cctk.io as io

    monkeypatch.setattr(io, "MAX_BYTES", 8)
    path = tmp_path / "append"
    path.write_text("123456")
    with io.safe_open(path, "a") as handle:
        with pytest.raises(CCTKError):
            handle.write("789")
    assert path.read_text() == "123456"


def test_formula_text_escaped_only_in_table():
    assert table_cell("=HYPERLINK()").startswith("'")
    doc = document()
    doc["arrays"][0]["id"] = "=unsafe"
    assert compare(doc)["arrays"][0]["id"] == "=unsafe"


def test_order_deterministic_and_source_scoped_ids():
    doc = document()
    doc["arrays"].reverse()
    assert [a["id"] for a in compare(doc)["arrays"]] == ["a2", "a1", "a0"]
    assert [a["id"] for a in compare(doc, order="id")["arrays"]] == ["a0", "a1", "a2"]
    # The repeated original array_1 and occurrence IDs are scoped, not conflated.
    assert len(compare(doc)["arrays"]) == 3


@pytest.mark.parametrize("all_deleted", [False, True])
def test_gap_and_label_visible_on_canvas(tmp_path, monkeypatch, all_deleted):
    import matplotlib.figure

    doc = document()
    doc["arrays"] = doc["arrays"][:1]
    doc["arrays"][0]["spacers"] = [
        {"id": "gap", "ordinal": 1, "sequence": None, "status": "deletion"}
    ]
    if not all_deleted:
        doc["arrays"][0]["spacers"].append(
            {"id": "dna", "ordinal": 2, "sequence": "ACGT", "status": "observed"}
        )
    original = matplotlib.figure.Figure.savefig
    tested = []

    def save(self, *args, **kwargs):
        self.canvas.draw()
        renderer = self.canvas.get_renderer()
        for axes in self.axes:
            for label in axes.texts:
                if label.get_text() in {"×", "a0"}:
                    bbox = label.get_window_extent(renderer)
                    assert self.bbox.contains(bbox.x0, bbox.y0)
                    assert self.bbox.contains(bbox.x1, bbox.y1)
                    tested.append(label.get_text())
        return original(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", save)
    result = compare(doc)
    assert result["arrays"][0]["evidence_status"] == (
        "no_observed_spacers" if all_deleted else "observed"
    )
    plot_comparison(result, tmp_path / "visible.png")
    assert set(tested) == {"×", "a0"}


def test_duplicate_connection_work_is_bounded(tmp_path, monkeypatch):
    from cctk import plotting

    def forbidden(*args, **kwargs):
        raise AssertionError("unbounded connection work was entered")

    monkeypatch.setattr(plotting, "draw_diffplot_lines", forbidden)
    plotting.plot_diffplot(
        {"a": ["s"] * 101, "b": ["s"] * 101},
        ["a", "b"],
        ["s"],
        {"s": ("#123456", "#123456")},
        outfile=str(tmp_path / "bounded.png"),
        plot_width=8,
        plot_height=3,
        text_size=8,
        dpi=80,
    )


@pytest.mark.parametrize("strand", [[], {}, 0, False])
def test_malformed_strand_is_domain_error(strand):
    doc = document()
    doc["arrays"][0]["strand"] = strand
    with pytest.raises(CCTKError):
        compare(doc)


@pytest.mark.parametrize("options", [{"identity": []}, {"order": {}}, {"min_shared": True}])
def test_malformed_options_are_domain_errors(options):
    with pytest.raises(CCTKError):
        compare(document(), **options)


def test_invalid_unicode_identifier_is_rejected_before_publication():
    doc = document()
    doc["arrays"][0]["id"] = "\ud800"
    with pytest.raises(CCTKError, match="identifier"):
        compare(doc)


@pytest.mark.parametrize("payload", [b'{"schema":', b'\xff', b'[' * 2000])
def test_malformed_json_is_domain_error_without_partial_output(tmp_path, payload):
    source = tmp_path / "input.json"
    source.write_bytes(payload)
    with pytest.raises(CCTKError, match="invalid JSON"):
        run_compare(source, tmp_path / "out", identity="exact", min_shared=1, order="supplied")
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob(".cctk-*"))


@pytest.mark.parametrize("competitor", ["empty", "file", "symlink"])
def test_atomic_publication_does_not_replace_competing_output(tmp_path, monkeypatch, competitor):
    import cctk.comparison as comparison

    source = tmp_path / "input.json"
    source.write_text(json.dumps(document()))
    target = tmp_path / "out"
    original = comparison.publish_directory
    other = tmp_path / "other"
    other.mkdir()

    def compete(staging, destination):
        if competitor == "empty":
            destination.mkdir()
        elif competitor == "file":
            destination.write_text("competitor")
        else:
            destination.symlink_to(other, target_is_directory=True)
        original(staging, destination)

    monkeypatch.setattr(comparison, "publish_directory", compete)
    with pytest.raises(CCTKError, match="appeared"):
        run_compare(source, target, identity="exact", min_shared=1, order="supplied")
    if competitor == "empty":
        assert target.is_dir() and not list(target.iterdir())
    elif competitor == "file":
        assert target.read_text() == "competitor"
    else:
        assert target.is_symlink() and target.resolve() == other
    assert not list(tmp_path.glob(".cctk-*"))


def test_output_budget_failure_leaves_no_partial_comparison(tmp_path, monkeypatch):
    import cctk.io as io

    source = tmp_path / "input.json"
    source.write_text(json.dumps(document()))
    monkeypatch.setattr(io, "MAX_BYTES", source.stat().st_size + 1)
    with pytest.raises(CCTKError, match="bounded"):
        run_compare(source, tmp_path / "out", identity="exact", min_shared=1, order="supplied")
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob(".cctk-*"))


@pytest.mark.parametrize("outlined", [False, True])
def test_repeated_spacer_connectors_draw_each_occurrence_pair_once(outlined):
    from cctk.plotting import draw_diffplot_lines

    class Axes:
        def __init__(self):
            self.lines = []

        def plot(self, xs, ys, **kwargs):
            self.lines.append((tuple(xs), tuple(ys), kwargs["zorder"]))

    axes = Axes()
    draw_diffplot_lines(
        axes, 0, ["a", "b"], {"a": ["s"] * 2, "b": ["s"] * 3},
        {"s": ("#123456", "#000000")}, 2, 0.35, 0.5, 0.3, 0.1, 1, outlined,
    )
    assert len(axes.lines) == 6 * (2 if outlined else 1)
    assert len(set(axes.lines)) == len(axes.lines)


def test_connector_budget_counts_pairs_without_duplicate_overdraw(tmp_path, monkeypatch):
    from cctk import plotting

    rows = []
    monkeypatch.setattr(plotting, "draw_diffplot_lines", lambda *a, **k: rows.append(a[1]))
    plotting.plot_diffplot(
        {"a": ["s"] * 30, "b": ["s"] * 30}, ["a", "b"], ["s"],
        {"s": ("#123456", "#123456")}, outfile=str(tmp_path / "pairs.png"),
        plot_width=8, plot_height=3, text_size=8, dpi=80,
    )
    assert rows == [0]  # 900 pairs fit the limit; 27,000 overdrawn lines did not.


def test_observed_labels_are_literal_and_do_not_use_legacy_label_geometry(tmp_path, monkeypatch):
    from cctk import plotting

    def forbidden(*args, **kwargs):
        raise AssertionError("legacy label geometry is unused in the observed viewport")

    monkeypatch.setattr(plotting, "calc_vh_ratio_and_label_pad_diffplot", forbidden)
    doc = document()
    doc["arrays"][0]["id"] = r"$\unknowncommand$"
    plot_comparison(compare(doc), tmp_path / "literal.png")
    assert (tmp_path / "literal.png").is_file()


def test_transitive_components_preserve_scoped_occurrences_and_absence():
    slots = [["A", "A", None, "G"], ["A", "C"], ["C"], [None]]
    doc = {"schema": {"name": "cctk.arrays", "version": "1.0.0"}, "arrays": [
        {
            "id": ident, "source_id": "source", "array_id": ident, "strand": None,
            "spacers": [
                {"id": f"occ{n}", "ordinal": n, "sequence": sequence,
                 "status": "deletion" if sequence is None else "observed"}
                for n, sequence in enumerate(sequences, 1)
            ],
        }
        for ident, sequences in zip(["a", "b", "c", "d"], slots)
    ]}
    result = compare(doc)
    assert [(e["a"], e["b"], e["nshared"]) for e in result["edges"]] == [
        ("a", "b", 1), ("b", "c", 1)
    ]
    assert [e["jaccard"] for e in result["edges"]] == [1 / 3, 1 / 2]
    assert result["components"] == [["a", "b", "c"], ["d"]]
    assert result["isolates"] == ["d"]
    assert result["arrays"][0]["observed_spacer_count"] == 3
    assert result["arrays"][-1]["evidence_status"] == "no_observed_spacers"
    assert compare(doc, min_shared=2)["components"] == [["a"], ["b"], ["c"], ["d"]]


def test_multiline_fasta_preserves_record_order_case_and_ambiguity(tmp_path):
    source = tmp_path / "multiline.fa"
    source.write_text(">b description\n  aCgT  \nRYSW\n\nKMBDHVN\n>a\nAC\nGT\n")
    records = read_fasta(source)
    assert list(records) == ["b description", "a"]
    assert records == {"b description": "aCgTRYSWKMBDHVN", "a": "ACGT"}


@pytest.mark.parametrize("identity", ["exact", "reverse-complement"])
def test_repeated_dna_equivalence_retains_every_occurrence(identity):
    doc = document()
    doc["arrays"] = doc["arrays"][:1]
    sequences = ["ACGTTGCA", "TGCAACGT", "ACGTTGCA", None]
    doc["arrays"][0]["spacers"] = [
        {"id": f"occ{n}", "ordinal": n, "sequence": dna,
         "status": "observed" if dna is not None else "deletion"}
        for n, dna in enumerate(sequences, 1)
    ]
    result = compare(doc, identity=identity)
    slots = result["arrays"][0]["spacers"]
    assert [(s["id"], s["ordinal"], s["sequence"]) for s in slots] == [
        (f"occ{n}", n, dna) for n, dna in enumerate(sequences, 1)
    ]
    expected = "s_" + hashlib.sha256(b"ACGTTGCA").hexdigest()
    assert slots[0]["symbol"] == slots[2]["symbol"] == expected
    assert (slots[1]["symbol"] == expected) == (identity == "reverse-complement")
    assert slots[-1]["symbol"] is None
    assert result["arrays"][0]["observed_spacer_count"] == 3
