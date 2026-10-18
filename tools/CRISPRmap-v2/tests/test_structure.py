import importlib.util
import json
import sqlite3

import pytest

from crisprmap import structure
from crisprmap.catalog_index import file_sha256
from crisprmap.cli import modern_command, parser
from crisprmap.io import InputError, LimitError
from crisprmap.structure import fold_catalog, fold_sequence
from test_catalog_index import make_catalog


@pytest.mark.skipif(importlib.util.find_spec('RNA') is None, reason='optional pinned ViennaRNA absent')
def test_two_structure_hypotheses_preserve_orientation_and_thermodynamic_evidence():
    result = fold_sequence('GCGCGCAAAAGCGCGC')
    assert [r['orientation_hypothesis'] for r in result] == [1, -1]
    for row in result:
        assert len(row['dot_bracket']) == 16
        assert 0 <= row['mfe_structure_probability'] <= 1
        assert 0 <= row['normalized_ensemble_defect_to_mfe'] <= 1
    with pytest.raises(InputError, match='version differs'):
        fold_sequence('ACGT', expected_version='0.0.0')


def mock_predictions(sequence, **_):
    return [{'orientation_hypothesis': orientation, 'dot_bracket': '.' * len(sequence),
             'mfe_kcal_mol': 0., 'ensemble_free_energy_kcal_mol': 0.,
             'mfe_structure_probability': 1., 'ensemble_mean_base_pair_distance': 0.,
             'normalized_ensemble_defect_to_mfe': 0.} for orientation in (1, -1)]


def prediction_rows(directory):
    with sqlite3.connect(directory / 'structure.sqlite3') as db:
        return db.execute('SELECT * FROM predictions ORDER BY node_id,orientation').fetchall()


@pytest.mark.parametrize('checkpoint_nodes', [1, 2])
@pytest.mark.parametrize('interruption', [LimitError, KeyboardInterrupt])
def test_fold_failure_rolls_back_batch_and_resume_keeps_prior_nodes(tmp_path, monkeypatch,
                                                                  checkpoint_nodes, interruption):
    source = make_catalog(tmp_path, monkeypatch, ['A' * 19 + base for base in 'ACGT'] + ['C' * 20])
    with sqlite3.connect(source / 'catalog.sqlite3') as db:
        db.execute('UPDATE nodes SET bona_fide_consensus_count=1')
    manifest = json.loads((source / 'manifest.json').read_text())
    manifest['database']['sha256'] = file_sha256(source / 'catalog.sqlite3')
    (source / 'manifest.json').write_text(json.dumps(manifest))
    calls = []

    def interrupted(sequence, **kwargs):
        calls.append(sequence)
        if len(calls) == 4:
            raise interruption('injected folding interruption')
        return mock_predictions(sequence, **kwargs)

    monkeypatch.setattr(structure, 'fold_sequence', interrupted)
    resumed = tmp_path / 'resumed'
    with pytest.raises(interruption, match='injected folding'):
        fold_catalog(source, resumed, checkpoint_nodes=checkpoint_nodes, min_free_bytes=0)
    retained = 3 if checkpoint_nodes == 1 else 2
    rows = prediction_rows(resumed)
    assert len(rows) == retained * 2
    assert {row[0] for row in rows} == set(range(1, retained + 1))
    assert not (resumed / 'manifest.json').exists()
    calls.clear()

    def resumed_fold(sequence, **kwargs):
        calls.append(sequence)
        return mock_predictions(sequence, **kwargs)

    monkeypatch.setattr(structure, 'fold_sequence', resumed_fold)
    result = fold_catalog(source, resumed, checkpoint_nodes=100, min_free_bytes=0)
    assert len(calls) == 5 - retained
    assert len(prediction_rows(resumed)) == 10  # Final short batch is committed.
    clean = tmp_path / 'clean'
    clean_result = fold_catalog(source, clean, checkpoint_nodes=2, min_free_bytes=0)
    assert prediction_rows(resumed) == prediction_rows(clean)
    assert result['policy'] == clean_result['policy']
    assert result['coverage'] == clean_result['coverage']
    assert result['invocation_policy']['checkpoint_nodes'] == 100
    with sqlite3.connect(resumed / 'structure.sqlite3') as db:
        resumed_binding = db.execute("SELECT value FROM metadata WHERE key='binding'").fetchone()[0]
    with sqlite3.connect(clean / 'structure.sqlite3') as db:
        assert db.execute("SELECT value FROM metadata WHERE key='binding'").fetchone()[0] == resumed_binding


@pytest.mark.parametrize('value', [0, 10001, True, 1.5])
def test_fold_rejects_invalid_checkpoint_bound_before_access(tmp_path, value):
    with pytest.raises(InputError, match='checkpoint_nodes'):
        fold_catalog(tmp_path / 'absent', tmp_path / 'output', checkpoint_nodes=value)
    assert not (tmp_path / 'output').exists()


def test_fold_cli_forwards_explicit_and_default_checkpoint_size(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(structure, 'fold_catalog', lambda *args, **kwargs: seen.append(kwargs))
    command = ['catalog-fold', '--catalog', str(tmp_path / 'source'), '--output', str(tmp_path / 'output')]
    modern_command(parser().parse_args(command))
    modern_command(parser().parse_args(command + ['--checkpoint-nodes', '250']))
    assert [kwargs['checkpoint_nodes'] for kwargs in seen] == [100, 250]
