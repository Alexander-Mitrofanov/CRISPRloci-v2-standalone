"""Separate secondary-structure predictions for both orientation hypotheses.

RNA folding is supporting evidence. It neither assigns transcription strand nor
validates a shared structural motif. No sequence-group label becomes a motif.
"""
from __future__ import annotations

import importlib.metadata
import math
from pathlib import Path
import sqlite3
import shutil
import time

from .catalog_index import file_sha256
from .core import reverse_complement
from .io import InputError, LimitError, atomic_json, integer

PINNED_VIENNA = '2.7.0'


def fold_sequence(sequence, *, expected_version=PINNED_VIENNA):
    try:
        import RNA
        version = importlib.metadata.version('ViennaRNA')
    except ImportError as exc:
        raise InputError('Structure evidence requires the optional ViennaRNA dependency') from exc
    if version != expected_version or getattr(RNA, '__version__', None) != expected_version:
        raise InputError('ViennaRNA version differs from the explicitly pinned version')
    if not isinstance(sequence, str) or not 1 <= len(sequence) <= 200 or set(sequence) - set('ACGT'):
        raise InputError('Folding requires 1–200 unambiguous DNA bases')
    RNA.params_load_RNA_Turner2004()
    # Constructor keywords rebuild pair tables and ignore mutable process-wide
    # defaults. Copy the reviewed Turner2004 tables before creating compounds.
    model = RNA.md(
        temperature=37.0, dangles=2, noLP=0, salt=1.021,
        noGU=0, noGUclosure=0, gquad=0, circ=0, min_loop_size=3,
        compute_bpp=1, energy_set=0, logML=0, special_hp=1,
        betaScale=1.0, pf_smooth=1, sfact=1.07, backtrack=1,
        backtrack_type="F", max_bp_span=-1, window_size=-1,
    )
    parameters = RNA.param(model)
    partition_parameters = RNA.exp_param(model)
    predictions = []
    for orientation, dna in ((1, sequence), (-1, reverse_complement(sequence))):
        compound = RNA.fold_compound(dna.replace('T', 'U'), model)
        compound.params_subst(parameters)
        compound.exp_params_subst(partition_parameters)
        structure, mfe = compound.mfe()
        compound.exp_params_rescale(mfe)
        _, ensemble_free_energy = compound.pf()
        probability = compound.pr_structure(structure)
        diversity = compound.mean_bp_distance()
        defect = compound.ensemble_defect(structure)
        values = (mfe, ensemble_free_energy, probability, diversity, defect)
        if not all(math.isfinite(v) for v in values) or not 0 <= probability <= 1.000001 or not 0 <= defect <= 1.000001:
            raise InputError('RNA folding returned invalid thermodynamic values')
        predictions.append({'orientation_hypothesis': orientation, 'dot_bracket': structure,
            'mfe_kcal_mol': mfe, 'ensemble_free_energy_kcal_mol': ensemble_free_energy,
            'mfe_structure_probability': probability, 'ensemble_mean_base_pair_distance': diversity,
            'normalized_ensemble_defect_to_mfe': defect})
    return predictions


def fold_catalog(catalog_dir, destination, *, lane=0, lanes=1, max_nodes=10_000_000,
                 expected_version=PINNED_VIENNA, deadline_seconds=172800, min_free_bytes=2 * 1024**3,
                 checkpoint_nodes=100):
    """Fold core nodes in bounded transactions; resume committed node batches."""
    from .catalog import load_catalog
    integer(lanes, 'lanes', 1, 1024)
    integer(lane, 'lane', 0, lanes - 1)
    integer(max_nodes, 'max_nodes', 1, 100_000_000)
    integer(min_free_bytes, 'min_free_bytes', 0, 2**63 - 1)
    integer(checkpoint_nodes, 'checkpoint_nodes', 1, 10000)
    if type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds) or deadline_seconds <= 0:
        raise InputError('Deadline must be positive and finite')
    root = Path(destination).absolute()
    if any(p.is_symlink() for p in (root, *root.parents)):
        raise InputError('Output paths must not traverse symlinks')
    root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(root).free < min_free_bytes:
        raise LimitError('Structure free-space reserve would be exhausted')
    if (root / 'manifest.json').exists():
        raise InputError('Structure partition already complete')
    catalog, source = load_catalog(catalog_dir, verify_hash=True)
    catalog_hash = file_sha256(Path(catalog_dir) / 'manifest.json')
    dbpath = root / 'structure.sqlite3'
    binding = f'{catalog_hash}:{lane}:{lanes}:{expected_version}:37:dangles2:core_all'
    if dbpath.is_symlink():
        catalog.close()
        raise InputError('Structure checkpoint must not be a symlink')
    if dbpath.exists():
        probe = sqlite3.connect(dbpath.as_uri() + '?mode=ro', uri=True)
        try:
            old = probe.execute("SELECT value FROM metadata WHERE key='binding'").fetchone()
            if old is None or old[0] != binding:
                raise InputError('Checkpoint does not belong to this structure build')
        except sqlite3.Error as exc:
            catalog.close()
            raise InputError('Checkpoint does not belong to a structure build') from exc
        except BaseException:
            catalog.close()
            raise
        finally:
            probe.close()
    db = sqlite3.connect(dbpath)
    started = time.monotonic()
    try:
        db.executescript('''CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT);
          CREATE TABLE IF NOT EXISTS predictions(node_id INTEGER,orientation INTEGER,dot_bracket TEXT,
          mfe REAL,ensemble_free_energy REAL,mfe_probability REAL,ensemble_diversity REAL,ensemble_defect REAL,
          PRIMARY KEY(node_id,orientation)) WITHOUT ROWID;''')
        binding = f'{catalog_hash}:{lane}:{lanes}:{expected_version}:37:dangles2:core_all'
        old = db.execute("SELECT value FROM metadata WHERE key='binding'").fetchone()
        if old and old[0] != binding:
            raise InputError('Cannot resume changed structure source or policy')
        db.execute("INSERT OR IGNORE INTO metadata VALUES ('binding',?)", (binding,))
        db.commit()
        count = 0
        total = 0
        pending_nodes = 0
        for position, (node, sequence) in enumerate(catalog.execute(
                'SELECT node_id,sequence FROM nodes WHERE eligible=1 AND bona_fide_consensus_count>0 ORDER BY node_id')):
            total += 1
            if position % lanes != lane:
                continue
            count += 1
            if count % 1000 == 0 and shutil.disk_usage(root).free < min_free_bytes:
                raise LimitError('Structure free-space reserve would be exhausted')
            if count > max_nodes or time.monotonic() - started > deadline_seconds:
                raise LimitError('Structure budget reached; committed nodes can resume')
            if db.execute('SELECT count(*) FROM predictions WHERE node_id=?', (node,)).fetchone()[0] == 2:
                continue
            results = fold_sequence(sequence, expected_version=expected_version)
            if time.monotonic() - started > deadline_seconds:
                raise LimitError('Structure deadline reached; committed nodes can resume')
            for r in results:
                db.execute('INSERT INTO predictions VALUES (?,?,?,?,?,?,?,?)', (node, r['orientation_hypothesis'],
                    r['dot_bracket'], r['mfe_kcal_mol'], r['ensemble_free_energy_kcal_mol'],
                    r['mfe_structure_probability'], r['ensemble_mean_base_pair_distance'], r['normalized_ensemble_defect_to_mfe']))
            pending_nodes += 1
            if pending_nodes >= checkpoint_nodes:
                db.commit()
                pending_nodes = 0
        # The final partial batch must be durable before its completion marker.
        db.commit()
    except BaseException:
        db.rollback()
        raise
    finally:
        db.close()
        catalog.close()
    result = {'schema': {'name': 'crisprmap.structure-evidence', 'version': '1.0.0'},
        'status': 'complete', 'catalog': {'manifest_sha256': catalog_hash, 'database_sha256': source['database']['sha256']},
        'policy': {'engine': 'ViennaRNA', 'version': expected_version, 'temperature_celsius': 37,
            'dangles': 2, 'parameters': 'ViennaRNA_default_Turner2004', 'orientations': [1, -1],
            'scope': 'all_eligible_bona_fide_consensus_nodes', 'lane': lane, 'lanes': lanes},
        'coverage': {'core_nodes': total, 'nodes_in_lane': count, 'orientation_predictions': count * 2},
        'invocation_policy': {'checkpoint_nodes': checkpoint_nodes, 'batch_unit': 'newly_folded_nodes',
            'failure_recovery': 'rollback_uncommitted_batch'},
        'interpretation': 'Predicted secondary structures; no strand selection, validated motif, or covariation claim',
        'database': {'file': dbpath.name, 'sha256': file_sha256(dbpath), 'bytes': dbpath.stat().st_size},
        'elapsed_seconds_this_invocation': time.monotonic() - started}
    if shutil.disk_usage(root).free < min_free_bytes:
        raise LimitError('Structure free-space reserve would be exhausted')
    if time.monotonic() - started > deadline_seconds:
        raise LimitError('Structure deadline reached before publication; committed nodes can resume')
    atomic_json(root / 'manifest.json', result)
    return result
