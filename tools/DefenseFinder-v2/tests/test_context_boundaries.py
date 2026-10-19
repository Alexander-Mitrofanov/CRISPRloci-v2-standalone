"""Independent native-row controls for source identity and parsing limits."""
import csv
import hashlib
import logging
from types import SimpleNamespace

import pytest

from defense_context import pipeline


@pytest.fixture
def work(tmp_path):
    for family in pipeline.GROUPS:
        folder = tmp_path / family
        (folder / "hmmer_results").mkdir(parents=True)
        (folder / "best_solution.tsv").write_text("# No Systems found\n")
    return tmp_path


def gene(identity="gene", position=1, source="source", replicon="replicon"):
    return dict(id=identity, position=position, source_id=source, replicon=replicon,
                protein="M" * 100, start=100 * position, end=100 * position + 99)


def extract(work, rows, family="DefenseFinder_1", name="profile"):
    path = work / family / "hmmer_results" / f"{name}.res_hmm_extract"
    with path.open("w", newline="") as stream:
        csv.writer(stream, delimiter="\t").writerows(rows)
    return path


def hit(score="45", start="1"):
    return ["gene", "replicon", "1", "100", "profile", "1e-12", score, "0.9", "0.8", start, "80"]


def system_row(**changes):
    row = dict(hit_id="gene", replicon="replicon", hit_pos="1", hit_seq_len="100",
               sys_id="system_a", model_fqn="defense-finder-models/DefenseFinder_1/RM/RM_II",
               sys_wholeness="1.000", sys_score="2.000", gene_name="profile",
               hit_status="mandatory", hit_score="45", hit_i_eval="1e-12", hit_profile_cov="0.9")
    return row | changes


def systems(work, rows):
    with (work / "DefenseFinder_1/best_solution.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.parametrize("text", [">\nACGT\n", ">   \nACGT\n", ">good\nACGU\n"])
def test_invalid_headers_and_non_dna_fail_with_value_error(tmp_path, text):
    path = tmp_path / "input.fna"
    path.write_text(text)
    with pytest.raises(ValueError):
        pipeline.read_fasta(path)


def test_many_short_lines_join_once_and_keep_source_order(tmp_path):
    path = tmp_path / "input.fna"
    path.write_text(">second description\n" + "a\n" * 100_000 + ">first\nryswkmbdhvn\n")
    assert list(pipeline.read_fasta(path).items()) == [("second", "A" * 100_000), ("first", "RYSWKMBDHVN")]


@pytest.mark.parametrize("text", [">a\n" + "A" * 2_000_001, "".join(f">r{i}\nA\n" for i in range(21))])
def test_source_size_and_count_limits(tmp_path, text):
    path = tmp_path / "input.fna"
    path.write_text(text)
    with pytest.raises(ValueError):
        pipeline.read_fasta(path)


@pytest.mark.parametrize("column,value", [(0, "unknown"), (1, "wrong"), (2, "2"), (3, "99"),
    (4, ""), (5, "nan"), (5, "-1"), (6, "inf"), (7, "nan"), (8, "-inf"), (9, "0"), (10, "101")])
def test_malformed_or_unbound_extracts_rejected(work, column, value):
    row = hit()
    row[column] = value
    extract(work, [row])
    with pytest.raises(ValueError):
        pipeline.collect(work, [gene()])


@pytest.mark.parametrize("kind", ["symlink", "oversized"])
def test_extracts_have_same_file_bounds_as_system_tables(work, kind):
    path = extract(work, [hit()])
    if kind == "symlink":
        original = path.with_suffix(".original")
        path.rename(original)
        path.symlink_to(original)
    else:
        with path.open("r+b") as stream:
            stream.truncate(pipeline.LIMIT + 1)
    with pytest.raises(ValueError):
        pipeline.collect(work, [gene()])


def test_profile_winner_and_tie_preserve_original_traversal(work):
    extract(work, [hit(score="45", start="1"), hit(score="46", start="2")], name="a")
    extract(work, [hit(score="46", start="3")], name="b")
    extract(work, [hit(score="46", start="4")], family="DefenseFinder_2")
    _, _, hits = pipeline.collect(work, [gene()])
    assert [(row["score"], row["protein_start"]) for row in hits] == [(46, 2)]


def test_shared_protein_systems_and_duplicate_profile_assignments_survive(work):
    systems(work, [system_row(), system_row(gene_name="other_profile"), system_row(sys_id="system_b")])
    result, components, _ = pipeline.collect(work, [gene()])
    assert [(row["id"], row["gene_ids"], row["start"], row["end"]) for row in result] == [
        ("system_a", ["gene"], 100, 199), ("system_b", ["gene"], 100, 199)]
    assert len(components) == 3


@pytest.mark.parametrize("changes", [dict(model_fqn="defense-finder-models/DefenseFinder_1/RM/other"),
    dict(sys_wholeness="0.500"), dict(sys_score="3.000"),
    dict(hit_id="other", replicon="other_replicon", hit_pos="2")])
def test_conflicting_system_metadata_rejected(work, changes):
    systems(work, [system_row(), system_row(**changes)])
    with pytest.raises(ValueError):
        pipeline.collect(work, [gene(), gene("other", source="other_source", replicon="other_replicon")])


@pytest.fixture
def stub_engine(monkeypatch):
    monkeypatch.setattr(pipeline, "version", lambda name: "2.1.4")
    monkeypatch.setattr(pipeline.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(stdout="HMMER 3.4 "))
    monkeypatch.setattr(pipeline, "verify_models", lambda path: dict(pipeline.MODEL))


def test_input_hash_uses_the_same_bytes_as_called_genes(tmp_path, monkeypatch, stub_engine):
    path = tmp_path / "source.fna"
    original = b">source\nACGTACGTACGT\n"
    path.write_bytes(original)

    def changed_after_read(records):
        assert records == {"source": "ACGTACGTACGT"}
        path.write_text(">changed\nTGCATGCA\n")
        return []

    monkeypatch.setattr(pipeline, "call_genes", changed_after_read)
    report = pipeline.run(path, tmp_path / "result", tmp_path, tmp_path / "hmmsearch", 1)
    assert report["input"]["sha256"] == hashlib.sha256(original).hexdigest()
    assert report["sources"] == [dict(id="source", length=12, sha256=hashlib.sha256(b"ACGTACGTACGT").hexdigest())]


@pytest.mark.parametrize("failure", [RuntimeError, KeyboardInterrupt])
def test_native_handlers_close_on_failure_or_interruption(tmp_path, monkeypatch, stub_engine, failure):
    path = tmp_path / "source.fna"
    path.write_text(">source\nACGT\n")
    monkeypatch.setattr(pipeline, "call_genes", lambda records: [gene()])
    logger = logging.getLogger("macsypy")
    previous = logger.handlers[:]
    logger.handlers = []
    handler = logging.StreamHandler()
    closed = []
    monkeypatch.setattr(handler, "close", lambda: closed.append(True))

    def fail(args):
        logger.addHandler(handler)
        raise failure("native failure")

    monkeypatch.setattr(pipeline.macsyfinder, "main", fail)
    try:
        with pytest.raises(failure):
            pipeline.run(path, tmp_path / "result", tmp_path, tmp_path / "hmmsearch", 1)
        assert closed == [True]
        assert not logger.handlers
        assert not (tmp_path / "result/report.json").exists()
    finally:
        logger.handlers = previous


def test_native_global_rank_preserves_exported_contig_local_position(work):
    first, second = gene(), gene("second", source="second_source", replicon="second_replicon")
    systems(work, [system_row(), system_row(hit_id="second", replicon="second_replicon", hit_pos="2", sys_id="system_b")])
    row = hit()
    row[:4] = ["second", "second_replicon", "2", "100"]
    extract(work, [row])
    result, components, hits = pipeline.collect(work, [first, second])
    assert second["position"] == 1
    assert [(system["source_id"], system["gene_ids"]) for system in result] == [
        ("second_source", ["second"]), ("source", ["gene"])]
    assert len(components) == 2 and hits[0]["gene_id"] == "second"


@pytest.mark.parametrize("workers", [0, 9, True, 1.0])
def test_invalid_workers_rejected_before_native_execution(tmp_path, monkeypatch, workers):
    def unexpected(*args, **kwargs):
        pytest.fail("Invalid workers must fail before native execution")
    monkeypatch.setattr(pipeline.subprocess, "run", unexpected)
    with pytest.raises(ValueError, match="Worker count"):
        pipeline.run(tmp_path, tmp_path, tmp_path, tmp_path, workers)


@pytest.mark.parametrize("location", ["environment", "virtualenv", "home", "cwd", "model", "symlink"])
def test_undeclared_configuration_rejected_before_analysis(tmp_path, monkeypatch, location):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("VIRTUAL_ENV", raising=False)
    monkeypatch.setenv("MACSY_CONF", str(tmp_path / "ambient.conf"))
    paths = {"environment": tmp_path / "ambient.conf", "virtualenv": tmp_path / "venv/etc/macsyfinder/macsyfinder.conf", "home": tmp_path / "home/.macsyfinder/macsyfinder.conf", "cwd": tmp_path / "macsyfinder.conf", "model": tmp_path / "models/defense-finder-models/RM/model_conf.xml", "symlink": tmp_path / "ambient.conf"}
    if location == "virtualenv":
        monkeypatch.setenv("VIRTUAL_ENV", str(tmp_path / "venv"))
    path = paths[location]
    path.parent.mkdir(parents=True, exist_ok=True)
    if location == "symlink":
        path.symlink_to(tmp_path / "missing")
    else:
        path.write_text("[score_opt]\nmandatory_weight = 3\n")
    monkeypatch.setattr(pipeline, "call_genes", lambda _: pytest.fail("configuration must fail before analysis"))
    with pytest.raises(ValueError, match="Undeclared MacSyFinder"):
        pipeline.run(tmp_path / "absent.fna", tmp_path / "result", tmp_path / "models", tmp_path / "unused", 1)


def test_config_resolution_matches_virtualenv_precedence(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("VIRTUAL_ENV", str(tmp_path / "venv"))
    poison = tmp_path / "unused.conf"
    poison.write_text("[score_opt]\nmandatory_weight = 3\n")
    monkeypatch.setenv("MACSY_CONF", str(poison))
    pipeline.reject_ambient_configuration(tmp_path / "models")


@pytest.mark.parametrize('column,value', [(7, '-0.1'), (7, '1.01'), (8, '-0.1'), (8, '1.01')])
def test_profile_coverages_are_fractions(work, column, value):
    row = hit()
    row[column] = value
    extract(work, [row])
    with pytest.raises(ValueError, match='coverage'):
        pipeline.collect(work, [gene()])


@pytest.mark.parametrize('changes', [dict(sys_wholeness='-0.1'),
    dict(hit_i_eval='-1'), dict(hit_profile_cov='-0.1'), dict(hit_profile_cov='1.01')])
def test_system_evidence_rejects_impossible_fractions_and_evalues(work, changes):
    systems(work, [system_row(**changes)])
    with pytest.raises(ValueError):
        pipeline.collect(work, [gene()])
