"""Reproducible operational sequence groups, never biological family labels.

Connected components use a lossless candidate enumerator and verified global
unit edit distance. Each worker retains a spanning forest, not a full graph.
Possible-only and observed-variant-only nodes cannot connect core components.
"""
from __future__ import annotations

from array import array
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time

from .catalog_index import CandidateIndex, distance_threshold, file_sha256, similarity_fraction
from .core import levenshtein, reverse_complement
from .io import InputError, LimitError, atomic_json, integer, read_json

THRESHOLDS = ("0.85", "0.90", "0.95")


def edit_distance(a, b, cutoff=None):
    try:
        from rapidfuzz.distance.Levenshtein import distance
    except ImportError:
        return levenshtein(a, b)
    return distance(a, b, score_cutoff=cutoff)


def orientation_distance(a, b, cutoff=None):
    return min(edit_distance(a, b, cutoff), edit_distance(reverse_complement(a), b, cutoff))


class UnionFind:
    def __init__(self, ids, position=None):
        self.position = position if position is not None else {node: i for i, node in enumerate(ids)}
        self.parent = array('I', range(len(ids)))
        self.size = array('I', [1]) * len(ids)

    def find(self, node):
        i = self.position[node]
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]
        return True


def _destination(path):
    path = Path(path).absolute()
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise InputError('Output path traverses a symlink')
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_forest(catalog_dir, index_dir, destination, *, lane=0, lanes=1,
                 thresholds=THRESHOLDS, max_candidates=2_000_000,
                 max_postings=50_000_000, deadline_seconds=172800, max_core_nodes=2_000_000,
                 checkpoint_queries=16):
    """Build/resume one deterministic lane; all lanes are required for assembly.

    A pair belongs to the lane of its smaller node ID. Skipping verification
    inside an existing component preserves exact components at that threshold.
    A clean interruption rolls back the uncommitted query batch; committed
    forests reconstruct union state on resume. Hard termination can leave a hot
    SQLite journal requiring operator recovery before the ownership probe can
    read the checkpoint. Checkpoint frequency can change on resume.
    """
    integer(lanes, 'lanes', 1, 1024)
    integer(lane, 'lane', 0, lanes - 1)
    integer(max_core_nodes, 'max_core_nodes', 1, 10_000_000)
    integer(checkpoint_queries, 'checkpoint_queries', 1, 1024)
    if type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds) or deadline_seconds <= 0:
        raise InputError('Deadline must be positive and finite')
    if not isinstance(thresholds, (list, tuple)) or not 1 <= len(thresholds) <= 8:
        raise InputError('Choose one to eight sequence similarity thresholds')
    fractions = sorted({similarity_fraction(t) for t in thresholds})
    if not fractions:
        raise InputError('At least one threshold is required')
    labels = [format(float(t), '.4f').rstrip('0').rstrip('.') for t in fractions]
    root = _destination(destination)
    if (root / 'manifest.json').exists():
        raise InputError('Forest already complete')
    started = time.monotonic()
    with CandidateIndex(catalog_dir, index_dir) as index:
        binding = {'catalog_manifest_sha256': file_sha256(Path(catalog_dir) / 'manifest.json'),
                   'catalog_database_sha256': index.catalog_manifest['database']['sha256'],
                   'index_manifest_sha256': file_sha256(Path(index_dir) / 'manifest.json'),
                   'thresholds': labels, 'lane': lane, 'lanes': lanes,
                   'core': 'eligible_and_bona_fide_consensus_count_positive',
                   'method': 'global_unit_edit_both_orientations_connected_components'}
        count = index.catalog.execute('SELECT count(*) FROM nodes WHERE eligible=1 AND bona_fide_consensus_count>0').fetchone()[0]
        if count > max_core_nodes:
            raise LimitError('Core node count exceeds the configured memory budget')
        rows = index.catalog.execute('SELECT node_id FROM nodes WHERE eligible=1 AND '
                                     'bona_fide_consensus_count>0 ORDER BY node_id').fetchall()
        ids = [r[0] for r in rows]
        positions = {node: i for i, node in enumerate(ids)}
        forests = {label: UnionFind(ids, positions) for label in labels}
        del rows
        dbpath = root / 'forest.sqlite3'
        encoded = json.dumps(binding, sort_keys=True)
        if dbpath.is_symlink():
            raise InputError('Forest checkpoint must not be a symlink')
        if dbpath.exists():
            probe = sqlite3.connect(dbpath.as_uri() + '?mode=ro', uri=True)
            try:
                old = probe.execute("SELECT value FROM metadata WHERE key='binding'").fetchone()
                if old is None or old[0] != encoded:
                    raise InputError('Checkpoint does not belong to this forest build')
            except sqlite3.Error as exc:
                raise InputError('Checkpoint does not belong to a forest build') from exc
            finally:
                probe.close()
        db = sqlite3.connect(dbpath)
        try:
            db.executescript('''CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT);
                CREATE TABLE IF NOT EXISTS edges(threshold TEXT,a INTEGER,b INTEGER,distance INTEGER,
                  PRIMARY KEY(threshold,a,b)) WITHOUT ROWID;''')
            encoded = json.dumps(binding, sort_keys=True)
            previous = db.execute("SELECT value FROM metadata WHERE key='binding'").fetchone()
            if previous and previous[0] != encoded:
                raise InputError('Cannot resume forest with changed source or policy')
            db.execute("INSERT OR IGNORE INTO metadata VALUES ('binding',?)", (encoded,))
            db.execute("INSERT OR IGNORE INTO metadata VALUES ('next_position','0')")
            db.commit()
            for label, a, b in db.execute('SELECT threshold,a,b FROM edges'):
                forests[label].join(a, b)
            position = int(db.execute("SELECT value FROM metadata WHERE key='next_position'").fetchone()[0])
            if not 0 <= position <= len(ids):
                raise InputError('Invalid forest resume position')
            pending_queries = 0
            for pos in range(position, len(ids)):
                a = ids[pos]
                if pos % lanes != lane:
                    continue
                if time.monotonic() - started > deadline_seconds:
                    raise LimitError('Forest build deadline reached; committed progress can resume')
                sequence = index.catalog.execute('SELECT sequence FROM nodes WHERE node_id=?', (a,)).fetchone()[0]
                for other in index.iter_candidates(sequence, max_distance=None,
                        min_similarity=float(fractions[0]), core_only=True, min_node_id=a,
                        max_candidates=max_candidates, max_postings=max_postings,
                        deadline_seconds=min(300, max(.01, deadline_seconds - (time.monotonic() - started)))):
                    if time.monotonic() - started > deadline_seconds:
                        raise LimitError('Forest deadline reached; committed progress can resume')
                    b, target = other['node_id'], other['sequence']
                    active = [(label, f) for label, f in zip(labels, fractions)
                              if forests[label].find(a) != forests[label].find(b)]
                    if not active:
                        continue
                    cutoff = distance_threshold(len(sequence), len(target), fractions[0])
                    d = orientation_distance(sequence, target, cutoff)
                    for label, fraction in active:
                        if d <= distance_threshold(len(sequence), len(target), fraction):
                            forests[label].join(a, b)
                            db.execute('INSERT INTO edges VALUES (?,?,?,?)', (label, a, b, d))
                db.execute("UPDATE metadata SET value=? WHERE key='next_position'", (str(pos + 1),))
                pending_queries += 1
                if pending_queries >= checkpoint_queries:
                    db.commit()
                    pending_queries = 0
            # Commit final partial batches, including the terminal resume marker.
            db.execute("UPDATE metadata SET value=? WHERE key='next_position'", (str(len(ids)),))
            db.commit()
            counts = dict(db.execute('SELECT threshold,count(*) FROM edges GROUP BY threshold'))
        except BaseException:
            db.rollback()
            raise
        finally:
            db.close()
        result = {'schema': {'name': 'crisprmap.cluster-forest', 'version': '1.0.0'},
                  'status': 'complete', 'binding': binding, 'core_nodes': len(ids),
                  'forest_edge_counts': counts, 'database': {'file': 'forest.sqlite3',
                  'sha256': file_sha256(dbpath), 'bytes': dbpath.stat().st_size},
                  'invocation_policy': {'checkpoint_queries': checkpoint_queries,
                      'batch_unit': 'completed_queries_in_lane', 'failure_recovery': 'rollback_uncommitted_batch'},
                  'elapsed_seconds_this_invocation': time.monotonic() - started}
        if time.monotonic() - started > deadline_seconds:
            raise LimitError('Cluster deadline reached before publication')
        atomic_json(root / 'manifest.json', result)
        return result


def assemble_clusters(catalog_dir, forest_dirs, destination, *, max_core_nodes=2_000_000, deadline_seconds=172800):
    """Validate all lanes and derive deterministic membership and chaining metrics."""
    from .catalog import load_catalog
    integer(max_core_nodes, 'max_core_nodes', 1, 10_000_000)
    if type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds) or deadline_seconds <= 0:
        raise InputError('Deadline must be positive and finite')
    started = time.monotonic()
    dirs = [Path(p).absolute() for p in forest_dirs]
    manifests = [read_json(p / 'manifest.json') for p in dirs]
    if not manifests:
        raise InputError('No forest partitions')
    for manifest in manifests:
        if (not isinstance(manifest, dict) or manifest.get('schema') != {'name': 'crisprmap.cluster-forest', 'version': '1.0.0'}
                or manifest.get('status') != 'complete' or not isinstance(manifest.get('binding'), dict)):
            raise InputError('Invalid forest completion manifest')
        needed = {'catalog_manifest_sha256','catalog_database_sha256','index_manifest_sha256','thresholds','lane','lanes','core','method'}
        if set(manifest['binding']) != needed:
            raise InputError('Invalid forest policy binding')
        integer(manifest['binding']['lanes'], 'forest lanes', 1, 1024)
        integer(manifest['binding']['lane'], 'forest lane', 0, manifest['binding']['lanes'] - 1)
        values = manifest['binding']['thresholds']
        if not isinstance(values, list) or not 1 <= len(values) <= 8:
            raise InputError('Invalid forest thresholds')
        expected = [format(float(f), '.4f').rstrip('0').rstrip('.') for f in sorted({similarity_fraction(v) for v in values})]
        if expected != values or manifest['binding']['core'] != 'eligible_and_bona_fide_consensus_count_positive' or manifest['binding']['method'] != 'global_unit_edit_both_orientations_connected_components':
            raise InputError('Unsupported forest policy')
    binding = manifests[0]['binding']
    common = {k: v for k, v in binding.items() if k != 'lane'}
    if len(manifests) != binding['lanes'] or {m['binding']['lane'] for m in manifests} != set(range(binding['lanes'])):
        raise InputError('Cluster forest lanes are missing or duplicated')
    if any(m['status'] != 'complete' or {k: v for k, v in m['binding'].items() if k != 'lane'} != common for m in manifests):
        raise InputError('Forest partitions have incompatible sources or policies')
    if file_sha256(Path(catalog_dir) / 'manifest.json') != binding['catalog_manifest_sha256']:
        raise InputError('Cluster source differs from forest source')
    catalog, catalog_manifest = load_catalog(catalog_dir, verify_hash=True)
    root = _destination(destination)
    if any(root.iterdir()):
        catalog.close()
        raise InputError('Cluster output directory must be empty')
    try:
        count = catalog.execute('SELECT count(*) FROM nodes WHERE eligible=1 AND bona_fide_consensus_count>0').fetchone()[0]
        if count > max_core_nodes:
            raise LimitError('Core node count exceeds the configured memory budget')
        nodes = {r[0]: (r[1], r[2]) for r in catalog.execute(
            'SELECT node_id,sequence,sequence_sha256 FROM nodes WHERE eligible=1 AND bona_fide_consensus_count>0')}
        labels = binding['thresholds']
        positions = {node: i for i, node in enumerate(sorted(nodes))}
        forests = {label: UnionFind(sorted(nodes), positions) for label in labels}
        for path, manifest in zip(dirs, manifests):
            if file_sha256(path / 'forest.sqlite3') != manifest['database']['sha256']:
                raise InputError('Forest checksum mismatch')
            source = sqlite3.connect((path / 'forest.sqlite3').as_uri() + '?mode=ro&immutable=1', uri=True)
            try:
                stored = dict(source.execute('SELECT key,value FROM metadata'))
                if (stored.get('binding') != json.dumps(manifest['binding'], sort_keys=True)
                        or stored.get('next_position') != str(len(nodes)) or manifest.get('core_nodes') != len(nodes)):
                    raise InputError('Incomplete or inconsistent forest checkpoint')
                for label, a, b, d in source.execute('SELECT threshold,a,b,distance FROM edges'):
                    if time.monotonic() - started > deadline_seconds:
                        raise LimitError('Cluster assembly deadline exceeded')
                    if label not in forests or a not in nodes or b not in nodes or not a < b or type(d) is not int or d < 0:
                        raise InputError('Invalid forest edge')
                    if d > distance_threshold(len(nodes[a][0]), len(nodes[b][0]), label) or orientation_distance(nodes[a][0], nodes[b][0], d) != d:
                        raise InputError('Invalid forest edge')
                    forests[label].join(a, b)
            finally:
                source.close()
        dbpath = root / 'clusters.sqlite3'
        db = sqlite3.connect(dbpath)
        try:
            db.executescript('''CREATE TABLE members(node_id INTEGER,threshold TEXT,cluster_id TEXT,
                PRIMARY KEY(threshold,node_id)) WITHOUT ROWID;
                CREATE TABLE clusters(cluster_id TEXT PRIMARY KEY,threshold TEXT,representative_node_id INTEGER,
                member_count INTEGER,member_sha256 TEXT,max_representative_distance INTEGER,
                outside_representative_threshold INTEGER,representative_distance_histogram TEXT);
                CREATE INDEX members_cluster ON members(cluster_id);''')
            summaries = {}
            for label in labels:
                groups = defaultdict(list)
                for node in sorted(nodes):
                    groups[forests[label].find(node)].append(node)
                histogram = Counter()
                outside_total = 0
                for members in groups.values():
                    members.sort(key=lambda n: nodes[n][0])
                    rep = members[0]
                    identity = hashlib.sha256()
                    identity.update(('crisprmap-sequence-cluster-v1:' + label + '\n').encode())
                    for node in members:
                        identity.update(bytes.fromhex(nodes[node][1]))
                    member_sha = identity.hexdigest()
                    cluster_id = 'MAP2-S' + str(int(Fraction(label) * 10000)) + '-' + member_sha[:24]
                    distances = Counter()
                    outside = 0
                    for node in members:
                        if time.monotonic() - started > deadline_seconds:
                            raise LimitError('Cluster summary deadline exceeded')
                        d = orientation_distance(nodes[rep][0], nodes[node][0])
                        distances[d] += 1
                        outside += d > distance_threshold(len(nodes[rep][0]), len(nodes[node][0]), label)
                    db.execute('INSERT INTO clusters VALUES (?,?,?,?,?,?,?,?)',
                        (cluster_id, label, rep, len(members), member_sha, max(distances), outside, json.dumps(dict(distances), sort_keys=True)))
                    db.executemany('INSERT INTO members VALUES (?,?,?)', [(n, label, cluster_id) for n in members])
                    outside_total += outside
                    histogram[len(members)] += 1
                db.commit()
                summaries[label] = {'clusters': len(groups), 'singleton_clusters': histogram[1],
                    'largest_cluster': max(histogram, default=0), 'members_outside_representative_threshold': outside_total,
                    'size_histogram': dict(sorted(histogram.items()))}
        finally:
            db.close()
        result = {'schema': {'name': 'crisprmap.sequence-clusters', 'version': '1.0.0'},
            'status': 'complete', 'catalog': {'manifest_sha256': binding['catalog_manifest_sha256'],
            'database_sha256': catalog_manifest['database']['sha256']}, 'policy': common,
            'interpretation': 'Operational sequence connected components; not validated biological families or structural motifs',
            'representative_policy': 'lexicographically_smallest_canonical_sequence',
            'edge_storage': 'verified_spanning_forests_not_full_graph', 'core_nodes': len(nodes),
            'thresholds': summaries, 'forest_manifest_sha256': [file_sha256(p / 'manifest.json') for p in dirs],
            'database': {'file': 'clusters.sqlite3', 'sha256': file_sha256(dbpath), 'bytes': dbpath.stat().st_size}}
        if time.monotonic() - started > deadline_seconds:
            raise LimitError('Cluster deadline reached before publication')
        atomic_json(root / 'manifest.json', result)
        return result
    finally:
        catalog.close()
