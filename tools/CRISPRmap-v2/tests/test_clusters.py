import json
import sqlite3

import pytest

from crisprmap import clusters
from crisprmap.catalog_index import CandidateIndex, build_index, file_sha256
from crisprmap.cli import modern_command, parser
from crisprmap.clusters import assemble_clusters, build_forest, orientation_distance
from crisprmap.io import InputError, LimitError
from test_catalog_index import make_catalog


def core_fixture(tmp_path, monkeypatch):
    # A--B--C path: endpoints are outside the edge threshold. X is only Possible.
    seqs = ['AAAAAAAAAACCCCCCCCCC', 'GAAAAAAAAACCCCCCCCCC', 'GGAAAAAAAACCCCCCCCCC',
            'GGGAAAAAAACCCCCCCCCC', 'TTGCGTACGATGACGTGTAG']
    source = make_catalog(tmp_path, monkeypatch, seqs)
    with sqlite3.connect(source / 'catalog.sqlite3') as db:
        db.execute('UPDATE nodes SET bona_fide_consensus_count=1')
        db.execute('UPDATE nodes SET bona_fide_consensus_count=0 WHERE sequence=?', (seqs[3],))
    manifest = json.loads((source / 'manifest.json').read_text())
    manifest['database']['sha256'] = file_sha256(source / 'catalog.sqlite3')
    (source / 'manifest.json').write_text(json.dumps(manifest))
    return source


def memberships(path):
    with sqlite3.connect(path / 'clusters.sqlite3') as db:
        return db.execute('SELECT * FROM members ORDER BY threshold,node_id').fetchall()


def test_complete_components_equal_split_lane_forests_and_do_not_bridge_possible(tmp_path, monkeypatch):
    source = core_fixture(tmp_path, monkeypatch)
    index = tmp_path / 'index'
    build_index(source, index, min_free_bytes=0)
    one = tmp_path / 'one'
    build_forest(source, index, one, thresholds=['0.95'])
    whole = tmp_path / 'whole'
    manifest = assemble_clusters(source, [one], whole)
    parts = [tmp_path / f'part{i}' for i in range(3)]
    for i, part in enumerate(parts):
        build_forest(source, index, part, lane=i, lanes=3, thresholds=['0.95'])
    split = tmp_path / 'split'
    assemble_clusters(source, parts, split)
    assert memberships(whole) == memberships(split)
    assert manifest['core_nodes'] == 4
    assert manifest['thresholds']['0.95']['largest_cluster'] == 3
    assert manifest['thresholds']['0.95']['members_outside_representative_threshold'] == 1


def test_missing_lanes_and_changed_edges_fail(tmp_path, monkeypatch):
    source = core_fixture(tmp_path, monkeypatch)
    index = tmp_path / 'index'
    build_index(source, index, min_free_bytes=0)
    a = tmp_path / 'a'
    build_forest(source, index, a, lane=0, lanes=2)
    with pytest.raises(InputError, match='missing or duplicated'):
        assemble_clusters(source, [a], tmp_path / 'bad')
    b = tmp_path / 'b'
    build_forest(source, index, b)
    with sqlite3.connect(b / 'forest.sqlite3') as db:
        db.execute('UPDATE edges SET distance=0')
    manifest = json.loads((b / 'manifest.json').read_text())
    manifest['database']['sha256'] = file_sha256(b / 'forest.sqlite3')
    (b / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(InputError, match='Invalid forest edge'):
        assemble_clusters(source, [b], tmp_path / 'corrupt')


def test_memory_budget_and_reverse_complement_distance(tmp_path, monkeypatch):
    source = core_fixture(tmp_path, monkeypatch)
    index = tmp_path / 'index'
    build_index(source, index, min_free_bytes=0)
    with pytest.raises(LimitError, match='memory budget'):
        build_forest(source, index, tmp_path / 'bounded', max_core_nodes=1)
    assert orientation_distance('AACCT', 'AGGTT') == 0


def test_resume_refuses_unowned_or_symlink_database(tmp_path, monkeypatch):
    source = core_fixture(tmp_path, monkeypatch)
    index = tmp_path / 'index'
    build_index(source, index, min_free_bytes=0)
    root = tmp_path / 'unowned'
    root.mkdir()
    unrelated = root / 'forest.sqlite3'
    with sqlite3.connect(unrelated) as db:
        db.execute('CREATE TABLE precious(value TEXT)')
        db.execute("INSERT INTO precious VALUES ('preserve')")
    before = unrelated.read_bytes()
    with pytest.raises(InputError, match='does not belong'):
        build_forest(source, index, root)
    assert unrelated.read_bytes() == before
    alias = tmp_path / 'alias'
    alias.mkdir()
    (alias / 'forest.sqlite3').symlink_to(unrelated)
    with pytest.raises(InputError, match='symlink'):
        build_forest(source, index, alias)
    assert unrelated.read_bytes() == before


@pytest.mark.parametrize('checkpoint_queries', [1, 2])
@pytest.mark.parametrize('interruption', [LimitError, KeyboardInterrupt])
def test_forest_failure_rolls_back_batch_and_resume_reconstructs_union_state(tmp_path, monkeypatch,
                                                                          checkpoint_queries, interruption):
    source = core_fixture(tmp_path, monkeypatch)
    index = tmp_path / 'index'
    build_index(source, index, min_free_bytes=0)
    original = CandidateIndex.iter_candidates
    calls = 0

    def interrupted(self, sequence, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 4:
            raise interruption('injected query interruption')
        yield from original(self, sequence, **kwargs)

    monkeypatch.setattr(CandidateIndex, 'iter_candidates', interrupted)
    resumed = tmp_path / 'resumed'
    with pytest.raises(interruption, match='injected query'):
        build_forest(source, index, resumed, thresholds=['0.95'], checkpoint_queries=checkpoint_queries)
    retained = 3 if checkpoint_queries == 1 else 2
    with sqlite3.connect(resumed / 'forest.sqlite3') as db:
        assert db.execute("SELECT value FROM metadata WHERE key='next_position'").fetchone()[0] == str(retained)
        # Query 3 connects the chain's last member. Its edge must roll back
        # with a 2-query checkpoint and persist with legacy per-query commits.
        assert db.execute('SELECT count(*) FROM edges').fetchone()[0] == (2 if retained == 3 else 1)
    assert not (resumed / 'manifest.json').exists()
    monkeypatch.setattr(CandidateIndex, 'iter_candidates', original)
    result = build_forest(source, index, resumed, thresholds=['0.95'], checkpoint_queries=16)
    clean = tmp_path / 'clean'
    clean_result = build_forest(source, index, clean, thresholds=['0.95'], checkpoint_queries=1)
    with sqlite3.connect(resumed / 'forest.sqlite3') as db:
        assert db.execute("SELECT value FROM metadata WHERE key='next_position'").fetchone()[0] == '4'
        resumed_edges = db.execute('SELECT * FROM edges ORDER BY threshold,a,b').fetchall()
    with sqlite3.connect(clean / 'forest.sqlite3') as db:
        assert db.execute('SELECT * FROM edges ORDER BY threshold,a,b').fetchall() == resumed_edges
    assert result['binding'] == clean_result['binding']
    assert result['invocation_policy']['checkpoint_queries'] == 16
    assemble_clusters(source, [resumed], tmp_path / 'resumed_clusters')
    assemble_clusters(source, [clean], tmp_path / 'clean_clusters')
    assert memberships(tmp_path / 'resumed_clusters') == memberships(tmp_path / 'clean_clusters')


@pytest.mark.parametrize('value', [0, 1025, True, 1.5])
def test_forest_rejects_invalid_checkpoint_bound_before_access(tmp_path, value):
    with pytest.raises(InputError, match='checkpoint_queries'):
        build_forest(tmp_path / 'absent', tmp_path / 'absent_index', tmp_path / 'output', checkpoint_queries=value)
    assert not (tmp_path / 'output').exists()


def test_forest_cli_forwards_explicit_and_default_checkpoint_size(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(clusters, 'build_forest', lambda *args, **kwargs: seen.append(kwargs))
    command = ['catalog-cluster-forest', '--catalog', str(tmp_path / 'source'), '--index', str(tmp_path / 'index'),
               '--output', str(tmp_path / 'output')]
    modern_command(parser().parse_args(command))
    modern_command(parser().parse_args(command + ['--checkpoint-queries', '32']))
    assert [kwargs['checkpoint_queries'] for kwargs in seen] == [16, 32]
