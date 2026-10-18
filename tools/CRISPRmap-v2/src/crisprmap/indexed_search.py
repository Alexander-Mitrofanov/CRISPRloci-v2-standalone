"""Stream modern-catalog search results without expanding occurrence inventories."""

from __future__ import annotations

import math
import os
from pathlib import Path
import shutil
import sqlite3
import time

from . import __version__
from .catalog_index import (COUNT_FIELDS, NODE_COLUMNS, CandidateIndex, distance_threshold,
                            file_sha256, similarity_fraction)
from .core import canonical_sequence, global_alignment, levenshtein, reverse_complement
from .io import InputError, LimitError, atomic_json, canonical, digest, integer, normalize_records, read_json

SCHEMA = {"name": "crisprmap.catalog-search", "version": "1.0.0"}


def _distance_backend():
    """Select equivalent unit-cost engines; record the engine actually loaded."""
    try:
        import rapidfuzz
        from rapidfuzz.distance.Levenshtein import distance as native_distance
    except ImportError:
        def python_distance(a, b, cutoff=None):
            return levenshtein(a, b)

        return python_distance, {"name": "crisprmap.myers", "version": __version__,
                          "insertion_cost": 1, "deletion_cost": 1, "substitution_cost": 1}

    def accelerated(a, b, cutoff=None):
        return native_distance(a, b, weights=(1, 1, 1), score_cutoff=cutoff)

    return accelerated, {"name": "rapidfuzz.Levenshtein", "version": rapidfuzz.__version__,
                         "insertion_cost": 1, "deletion_cost": 1, "substitution_cost": 1}


def search_catalog(records, catalog_dir, index_dir, destination, *, max_distance=3, top_k=20,
                   deadline_seconds=3600, max_candidates=2_000_000, max_postings=10_000_000,
                   max_total_candidates=100_000_000, max_total_postings=1_000_000_000,
                   max_queries=1_000_000, max_output_bytes=8 * 1024**3, verify_hash=True,
                   input_sha256=None, clusters_dir=None, max_cluster_groups=100_000,
                   cluster_similarity=0.9) -> dict:
    """Publish JSONL, all-best-node SQLite evidence and a completion manifest.

    Every closest node is retained in matches; JSONL only displays top_k nodes.
    Occurrence counts sum every tied node before this display cap. No historical
    families, motifs, Cas types or organism labels are assigned by this method.
    """
    integer(max_distance, "max_distance", 0, 200)
    integer(top_k, "top_k", 1, 1000)
    integer(max_queries, "max_queries", 1, 100_000_000)
    integer(max_output_bytes, "max_output_bytes", 1024 * 1024, 256 * 1024**3)
    integer(max_total_candidates, "max_total_candidates", 1, 10_000_000_000)
    integer(max_total_postings, "max_total_postings", 1, 100_000_000_000)
    integer(max_candidates, "max_candidates", 1, 100_000_000)
    integer(max_postings, "max_postings", 1, 1_000_000_000)
    integer(max_cluster_groups, "max_cluster_groups", 1, 1_000_000)
    cluster_fraction = similarity_fraction(cluster_similarity)
    cluster_label = format(float(cluster_fraction), ".4f").rstrip("0").rstrip(".")
    if type(verify_hash) is not bool:
        raise InputError("verify_hash must be boolean")
    if type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds) or deadline_seconds <= 0:
        raise InputError("Deadline must be positive and finite")
    if input_sha256 is not None and (not isinstance(input_sha256, str) or len(input_sha256) != 64
                                     or set(input_sha256) - set("0123456789abcdef")):
        raise InputError("input_sha256 must be null or lowercase SHA256")
    destination = Path(destination).absolute()
    if any(p.is_symlink() for p in (destination, *destination.parents)):
        raise InputError("Output paths must not traverse symlinks")
    started = time.monotonic()
    distance_fn, distance_backend = _distance_backend()
    database = None
    owned = False
    try:
        with CandidateIndex(catalog_dir, index_dir, verify_hash=verify_hash) as reader:
            if max_distance > reader.max_distance:
                raise InputError("Search distance exceeds index fixed-distance coverage")
            cluster_binding = None
            if clusters_dir is not None:
                clusters_dir = Path(clusters_dir).absolute()
                cluster_manifest = read_json(clusters_dir / "manifest.json")
                if (not isinstance(cluster_manifest, dict)
                        or cluster_manifest.get("schema") != {"name": "crisprmap.sequence-clusters", "version": "1.0.0"}
                        or cluster_manifest.get("status") != "complete"):
                    raise InputError("Unsupported or incomplete sequence-cluster release")
                if cluster_label not in cluster_manifest.get("thresholds", {}):
                    raise InputError("Requested cluster similarity is absent from the cluster release")
                if cluster_fraction < reader.min_similarity:
                    raise InputError("Requested cluster similarity is below index coverage")
                if (cluster_manifest.get("catalog", {}).get("manifest_sha256") != file_sha256(Path(catalog_dir) / "manifest.json")
                        or cluster_manifest.get("catalog", {}).get("database_sha256") != reader.catalog_manifest["database"]["sha256"]
                        or cluster_manifest.get("policy", {}).get("index_manifest_sha256") != file_sha256(Path(index_dir) / "manifest.json")):
                    raise InputError("Sequence clusters belong to a different catalog or index")
                cluster_path = clusters_dir / "clusters.sqlite3"
                descriptor = cluster_manifest.get("database", {})
                if (descriptor.get("file") != cluster_path.name or not cluster_path.is_file()
                        or cluster_path.stat().st_size != descriptor.get("bytes")
                        or file_sha256(cluster_path) != descriptor.get("sha256")):
                    raise InputError("Sequence-cluster database checksum or size mismatch")
                reader.connection.execute("ATTACH DATABASE ? AS clusters", (cluster_path.as_uri() + "?mode=ro&immutable=1",))
                cluster_binding = {"manifest_sha256": file_sha256(clusters_dir / "manifest.json"),
                                   "database_sha256": descriptor["sha256"]}
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.mkdir()
            owned = True
            database_path = destination / "results.sqlite3"
            database = sqlite3.connect(database_path)
            database.execute("PRAGMA journal_mode=DELETE")
            database.execute("PRAGMA synchronous=FULL")
            database.execute(f"PRAGMA max_page_count={max_output_bytes // 4096}")
            database.executescript("""
                CREATE TABLE queries(query_id TEXT PRIMARY KEY,sequence TEXT NOT NULL,
                  sequence_sha256 TEXT NOT NULL,status TEXT NOT NULL,best_distance INTEGER,
                  best_node_count INTEGER NOT NULL,summary_json TEXT NOT NULL);
                CREATE TABLE matches(query_id TEXT NOT NULL,node_id INTEGER NOT NULL,
                  distance INTEGER NOT NULL,orientation TEXT NOT NULL,
                  PRIMARY KEY(query_id,node_id)) WITHOUT ROWID;
                CREATE INDEX matches_node ON matches(node_id);
                CREATE TABLE cluster_matches(query_id TEXT NOT NULL,threshold TEXT NOT NULL,
                  cluster_id TEXT NOT NULL,best_node_count INTEGER NOT NULL,
                  PRIMARY KEY(query_id,threshold,cluster_id)) WITHOUT ROWID;
                CREATE TABLE cluster_candidates(query_id TEXT NOT NULL,threshold TEXT NOT NULL,
                  cluster_id TEXT NOT NULL,best_distance INTEGER NOT NULL,
                  qualifying_core_node_count INTEGER NOT NULL,
                  PRIMARY KEY(query_id,threshold,cluster_id)) WITHOUT ROWID;
            """)
            reader.connection.execute("CREATE TEMP TABLE best_nodes(node_id INTEGER PRIMARY KEY,distance INTEGER,orientation TEXT)")
            reader.connection.execute("CREATE TEMP TABLE qualifying_core_nodes(node_id INTEGER PRIMARY KEY,distance INTEGER)")
            counts = {"queries": 0, "exact": 0, "near": 0, "no_hit": 0, "unsupported": 0,
                      "best_node_matches": 0, "candidate_nodes": 0, "candidate_postings": 0,
                      "distance_comparisons": 0, "cluster_mapping_queries": 0,
                      "cluster_mapping_qualifying_core_nodes": 0}
            jsonl_bytes = 0
            jsonl_path = destination / "results.jsonl"

            def budget():
                remaining = deadline_seconds - (time.monotonic() - started)
                if remaining <= 0:
                    raise LimitError("Catalog search exceeded deadline")
                if jsonl_bytes + database_path.stat().st_size > max_output_bytes:
                    raise LimitError("Catalog search exceeded output byte budget")
                return remaining

            def account_candidates():
                counts["candidate_nodes"] += reader.last_stats["candidates"]
                counts["candidate_postings"] += reader.last_stats["postings"]
                if counts["candidate_nodes"] > max_total_candidates:
                    raise LimitError("Catalog search exceeded aggregate candidate budget")
                if counts["candidate_postings"] > max_total_postings:
                    raise LimitError("Catalog search exceeded aggregate posting budget")

            with jsonl_path.open("wb") as stream:
                for raw_record in records:
                    record = normalize_records([raw_record])[0]
                    counts["queries"] += 1
                    if counts["queries"] > max_queries:
                        raise LimitError("Catalog search exceeded query budget")
                    budget()
                    sequence = record["sequence"]
                    reader.connection.execute("DELETE FROM temp.best_nodes")
                    best = max_distance + 1
                    unsupported = bool(set(sequence) - set("ACGT"))
                    if not unsupported:
                        exact = reader.connection.execute(
                            "SELECT " + NODE_COLUMNS + " FROM catalog.nodes WHERE eligible=1 AND sequence=?",
                            (canonical_sequence(sequence),)).fetchone()
                        if exact is not None:
                            target = exact["sequence"]
                            orientation = "both" if target == reverse_complement(target) else (
                                "forward" if sequence == target else "reverse")
                            reader.connection.execute("INSERT INTO temp.best_nodes VALUES (?,0,?)",
                                                      (exact["node_id"], orientation))
                            best = 0
                        else:
                            for node in reader.iter_candidates(
                                sequence, max_distance=max_distance, max_candidates=max_candidates,
                                max_postings=max_postings, deadline_seconds=budget()
                            ):
                                forward = distance_fn(sequence, node["sequence"], max_distance)
                                target_rc = reverse_complement(node["sequence"])
                                reverse = forward if target_rc == node["sequence"] else distance_fn(sequence, target_rc, max_distance)
                                counts["distance_comparisons"] += 1 if target_rc == node["sequence"] else 2
                                current = min(forward, reverse)
                                if current > max_distance or current > best:
                                    continue
                                if current < best:
                                    reader.connection.execute("DELETE FROM temp.best_nodes")
                                    best = current
                                orientation = "both" if forward == reverse else ("forward" if forward < reverse else "reverse")
                                reader.connection.execute("INSERT INTO temp.best_nodes VALUES (?,?,?)",
                                                          (node["node_id"], best, orientation))
                            account_candidates()
                    aggregate = reader.connection.execute(
                        "SELECT count(*)," + ",".join(f"coalesce(sum(n.{field}),0)" for field in COUNT_FIELDS)
                        + " FROM temp.best_nodes b JOIN catalog.nodes n ON n.node_id=b.node_id").fetchone()
                    best_count = aggregate[0]
                    status = "unsupported" if unsupported else ("no_hit" if not best_count else ("exact" if best == 0 else "near"))
                    counts[status] += 1
                    counts["best_node_matches"] += best_count
                    mapping_status = "not_requested"
                    cluster_candidates = []
                    reader.connection.execute("DELETE FROM temp.qualifying_core_nodes")
                    if cluster_binding is not None:
                        if unsupported:
                            mapping_status = "unsupported_sequence"
                        elif status == "exact" and exact["bona_fide_consensus_count"] > 0:
                            # Every threshold-qualified neighbor of a core member
                            # belongs to that same connected component. Direct
                            # membership already establishes the complete set of
                            # qualifying cluster IDs without another graph query.
                            members = reader.connection.execute(
                                "SELECT cluster_id FROM clusters.members WHERE threshold=? AND node_id=?",
                                (cluster_label, exact["node_id"])).fetchall()
                            if len(members) != 1:
                                raise InputError("Core node has missing or conflicting cluster membership")
                            mapping_status = "not_needed_exact_core_member"
                        else:
                            counts["cluster_mapping_queries"] += 1
                            for node in reader.iter_candidates(
                                sequence, max_distance=None, min_similarity=cluster_similarity,
                                core_only=True, max_candidates=max_candidates, max_postings=max_postings,
                                deadline_seconds=budget()
                            ):
                                threshold = distance_threshold(len(sequence), node["length"], cluster_fraction)
                                forward = distance_fn(sequence, node["sequence"], threshold)
                                target_rc = reverse_complement(node["sequence"])
                                reverse = forward if target_rc == node["sequence"] else distance_fn(sequence, target_rc, threshold)
                                counts["distance_comparisons"] += 1 if target_rc == node["sequence"] else 2
                                distance = min(forward, reverse)
                                if distance <= threshold:
                                    reader.connection.execute("INSERT INTO temp.qualifying_core_nodes VALUES (?,?)",
                                                              (node["node_id"], distance))
                            account_candidates()
                            qualifying = reader.connection.execute("SELECT count(*) FROM temp.qualifying_core_nodes").fetchone()[0]
                            mapped = reader.connection.execute(
                                "SELECT count(*) FROM temp.qualifying_core_nodes q JOIN clusters.members m "
                                "ON m.node_id=q.node_id WHERE m.threshold=?", (cluster_label,)).fetchone()[0]
                            if mapped != qualifying:
                                raise InputError("Qualifying core nodes have missing or conflicting cluster memberships")
                            counts["cluster_mapping_qualifying_core_nodes"] += qualifying
                            for group in reader.connection.execute(
                                "SELECT m.cluster_id,min(q.distance),count(*) FROM temp.qualifying_core_nodes q "
                                "JOIN clusters.members m ON m.node_id=q.node_id WHERE m.threshold=? "
                                "GROUP BY m.cluster_id ORDER BY m.cluster_id", (cluster_label,)
                            ):
                                if len(cluster_candidates) >= max_cluster_groups:
                                    raise LimitError("Catalog search exceeded cluster-group budget")
                                cluster_candidates.append({"threshold": cluster_label, "cluster_id": group[0],
                                                           "best_distance": group[1], "qualifying_core_node_count": group[2]})
                                database.execute("INSERT INTO cluster_candidates VALUES (?,?,?,?,?)",
                                                 (record["id"], cluster_label, *group))
                                if len(cluster_candidates) % 1000 == 0:
                                    budget()
                            mapping_status = "completed_with_candidates" if cluster_candidates else "completed_no_candidates"
                    cluster_groups = []
                    if cluster_binding is not None:
                        for group in reader.connection.execute(
                            "SELECT m.threshold,m.cluster_id,count(*) FROM temp.best_nodes b "
                            "JOIN clusters.members m ON m.node_id=b.node_id "
                            "GROUP BY m.threshold,m.cluster_id ORDER BY m.threshold,m.cluster_id"
                        ):
                            if len(cluster_groups) >= max_cluster_groups:
                                raise LimitError("Catalog search exceeded cluster-group budget")
                            cluster_groups.append({"threshold": group[0], "cluster_id": group[1],
                                                   "best_node_count": group[2]})
                            database.execute("INSERT INTO cluster_matches VALUES (?,?,?,?)", (record["id"], *group))
                            if len(cluster_groups) % 1000 == 0:
                                budget()
                    hits = []
                    for node in reader.connection.execute(
                        "SELECT n.*,b.distance,b.orientation FROM temp.best_nodes b JOIN catalog.nodes n "
                        "ON n.node_id=b.node_id ORDER BY n.node_id LIMIT ?", (top_k,)):
                        aligned_orientation = "reverse" if node["orientation"] == "reverse" else "forward"
                        target = reverse_complement(node["sequence"]) if aligned_orientation == "reverse" else node["sequence"]
                        hits.append({"node_id": node["node_id"], "sequence_sha256": node["sequence_sha256"],
                                     "sequence": node["sequence"], "distance": node["distance"],
                                     "orientation": node["orientation"], "aligned_orientation": aligned_orientation,
                                     "similarity": 1 - node["distance"] / max(len(sequence), node["length"]),
                                     "occurrence_counts": {field: node[field] for field in COUNT_FIELDS},
                                     **global_alignment(sequence, target)})
                    result = {**record, "sequence_sha256": digest(sequence), "status": status,
                              "reason": "ambiguous_bases" if unsupported else (
                                  "no_reference_within_edit_budget" if status == "no_hit" else None),
                              "best_distance": best if best_count else None, "best_node_count": best_count,
                              "best_occurrence_counts": dict(zip(COUNT_FIELDS, aggregate[1:])),
                              "reported_node_count": len(hits), "hits_truncated": len(hits) < best_count,
                              "biological_family": None, "structural_motif": None,
                              "operational_clusters": cluster_groups,
                              "operational_cluster_candidates": cluster_candidates,
                              "cluster_mapping_status": mapping_status,
                              "cluster_mapping_threshold": cluster_label if cluster_binding is not None else None,
                              "cluster_relation": ("exact_core_membership" if status == "exact" else
                                                   "nearest_reference_cluster_evidence") if cluster_groups else None,
                              "hits": hits}
                    encoded = canonical(result)
                    jsonl_bytes += len(encoded)
                    if jsonl_bytes > max_output_bytes:
                        raise LimitError("Catalog search exceeded JSONL byte budget")
                    database.execute("INSERT INTO queries VALUES (?,?,?,?,?,?,?)",
                                     (record["id"], sequence, digest(sequence), status,
                                      result["best_distance"], best_count, encoded.decode().rstrip()))
                    for position, row in enumerate(reader.connection.execute(
                        "SELECT node_id,distance,orientation FROM temp.best_nodes ORDER BY node_id"
                    ), 1):
                        database.execute("INSERT INTO matches VALUES (?,?,?,?)", (record["id"], *row))
                        if position % 1000 == 0:
                            budget()
                    stream.write(encoded)
                    if counts["queries"] % 1000 == 0:
                        database.commit()
                    budget()
                if counts["queries"] == 0:
                    raise InputError("No query records supplied")
                stream.flush()
                os.fsync(stream.fileno())
            database.commit()
            budget()
            database.close()
            database = None
            manifest = {"schema": dict(SCHEMA), "status": "complete", "execution_completed": True,
                        "distance_backend": distance_backend,
                        "catalog_manifest_sha256": file_sha256(Path(catalog_dir) / "manifest.json"),
                        "index_manifest_sha256": file_sha256(Path(index_dir) / "manifest.json"),
                        "clusters": cluster_binding, "reference_hashes_verified": verify_hash,
                        "input_sha256": input_sha256, "counts": counts,
                        "policy": {"max_distance": max_distance, "top_k": top_k,
                                   "deadline_seconds": deadline_seconds, "max_candidates": max_candidates,
                                   "max_postings": max_postings, "max_total_candidates": max_total_candidates,
                                   "max_total_postings": max_total_postings, "max_queries": max_queries,
                                   "max_output_bytes": max_output_bytes,
                                   "max_cluster_groups": max_cluster_groups,
                                   "cluster_mapping_threshold": cluster_label,
                                   "cluster_mapping": "all_threshold_qualified_core_clusters_no_membership_changes",
                                   "normalization": "uppercase_U_to_T",
                                   "retrieval": "all_nearest_nodes_within_edit_budget_both_orientations",
                                   "annotation": "source_occurrence_counts_no_family_or_motif_assignment"},
                        "files": {name: {"sha256": file_sha256(destination / name),
                                         "bytes": (destination / name).stat().st_size}
                                  for name in ("results.jsonl", "results.sqlite3")}}
            budget()
            atomic_json(destination / "manifest.json", manifest)
            return manifest
    except sqlite3.IntegrityError as exc:
        raise InputError("Duplicate query identifier or inconsistent search evidence") from exc
    except (sqlite3.Error, OSError) as exc:
        raise InputError(f"Catalog search failed: {type(exc).__name__}") from exc
    finally:
        if database is not None:
            database.close()
        if owned and not (destination / "manifest.json").exists():
            shutil.rmtree(destination)
