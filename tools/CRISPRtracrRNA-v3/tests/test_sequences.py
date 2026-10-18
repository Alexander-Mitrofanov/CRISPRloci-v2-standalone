import gzip

import pytest

from crisprtracrrna.config import Settings, load_settings
from crisprtracrrna.errors import TracrError
from crisprtracrrna.models import AntiRepeat, Array, Interval, Sequence
from crisprtracrrna.search import make_candidates
from crisprtracrrna.sequence import input_files, read_fasta, reverse_complement, subsequence


def test_iupac_complement_and_lowercase():
    assert reverse_complement("acgtryswkmbdhvn") == "NBDHVKMWSRYACGT"
    assert reverse_complement(reverse_complement("ACGTRYSWKMBDHVN")) == "ACGTRYSWKMBDHVN"


def test_fasta_multirecord_gzip_and_full_identifiers(tmp_path):
    path = tmp_path / "a.b.fasta.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">seq.1 description, preserved\nacgu\nNN\n>seq.2\nGCTA\n")
    records = list(read_fasta(path))
    assert [(r.id, r.bases) for r in records] == [("seq.1", "ACGTNN"), ("seq.2", "GCTA")]
    assert records[0].description == "seq.1 description, preserved"
    (tmp_path / "ignore.txt").write_text("text")
    assert input_files(tmp_path) == [path]


@pytest.mark.parametrize("text", ["", "ACTG\n", ">x\n", ">x\nACZ\n", ">x\nAC\n>x\nGT\n", ">\nAC\n"])
def test_fasta_rejects_ambiguous_or_malformed_records(tmp_path, text):
    path = tmp_path / "test.fasta"
    path.write_text(text)
    with pytest.raises(TracrError):
        list(read_fasta(path))


def test_closed_coordinates_include_final_base_and_containment():
    interval = Interval.from_public(2, 4)
    assert subsequence("ACGTA", interval) == "CGT"
    assert interval.length == 3
    assert Interval(0, 100).overlap(Interval(30, 40)) == 10
    assert Interval(0, 5).distance(Interval(5, 10)) == 0
    assert Interval(0, 5).distance(Interval(8, 10)) == 3
    with pytest.raises(ValueError):
        Interval.from_public(0, 5)
    with pytest.raises(TracrError):
        subsequence("AA", Interval(0, 3))


@pytest.mark.parametrize(
    "strand,start,end", [("+", 0, 3), ("-", 0, 3), ("+", 17, 20), ("-", 17, 20)]
)
def test_candidate_clipping_and_roundtrip(strand, start, end):
    source = Sequence("s", "s", "ACGTTGCAACGTTGCAACGT")
    array = Array("a", Interval(5, 10), "ACGT", "+", "Possible", 0.8)
    anti = AntiRepeat(array, "ACGT", "+", Interval(start, end), strand, 0.9, 0.8, True)
    candidate = make_candidates(source, [anti], Settings())[0]
    assert candidate.to_genome(candidate.anti_local) == anti.interval
    assert candidate.sequence[candidate.anti_local.start : candidate.anti_local.end] == subsequence(
        source.bases, anti.interval, strand
    )
    assert len(candidate.sequence) == candidate.interval.length
    assert candidate.interval.end <= len(source.bases)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"threads": 0},
        {"coverage": 1.1},
        {"similarity": float("nan")},
        {"timeout": float("inf")},
        {"tail_length": -1},
        {"identify_command": "echo"},
        {"model_type": "X"},
        {"type_v_motifs": "False"},
        {"weights": {"madeup": 3}},
    ],
)
def test_invalid_settings(kwargs):
    with pytest.raises(TracrError):
        Settings(**kwargs)


def test_config_paths_and_typo_rejection(tmp_path):
    path = tmp_path / "settings.toml"
    path.write_text('models = ["data/custom.cm"]\nidentify_command = ["tools/identify"]\n')
    settings = load_settings(path, {})
    assert settings.models == (str(tmp_path / "data/custom.cm"),)
    assert settings.identify_command == (str(tmp_path / "tools/identify"),)
    path.write_text("threds = 2\n")
    with pytest.raises(TracrError, match="Unknown"):
        load_settings(path, {})


def test_command_path_keeps_virtual_environment_symlink(tmp_path):
    import sys

    command = tmp_path / "env/bin/python"
    command.parent.mkdir(parents=True)
    command.symlink_to(sys.executable)
    path = tmp_path / "settings.toml"
    path.write_text('casandra_command = ["env/bin/python", "-m", "casandra"]\n')
    settings = load_settings(path, {})
    assert settings.casandra_command[0] == str(command)
    assert settings.casandra_command[0] != str(command.resolve())
