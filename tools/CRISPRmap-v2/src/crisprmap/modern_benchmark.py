"""Genome-heldout repeat coverage with full exact-overlap and scoped homology audits.

Reference membership is rebuilt from source occurrences, never from global node
counts. This is a retrieval evaluation, not a biological family classifier test.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
import os
from importlib.metadata import version
from pathlib import Path
import shutil
import sqlite3
import tempfile
import time

from .catalog_index import CandidateIndex, distance_threshold, file_sha256
from .core import canonical_sequence, levenshtein, reverse_complement
from .io import InputError, LimitError, atomic_json, integer
from .reference import load_reference
from .source_import import check_free, validate_free_reserve

SCHEMA = {"name": "crisprmap.modern-benchmark", "version": "1.0.0"}
SIMILARITY = Fraction(9, 10)
MAX_ANCHORS = 100_000
ORACLE_QUERIES = 12
ORACLE_REFERENCES = 512


class _RunBudget:
    def __init__(self, output, deadline_seconds, min_free_bytes, max_candidates, max_postings):
        if type(deadline_seconds) not in (int, float) or not math.isfinite(deadline_seconds) or deadline_seconds <= 0:
            raise InputError("Benchmark deadline must be finite and positive")
        validate_free_reserve(min_free_bytes)
        integer(max_candidates, "max_total_candidates", 1, 10**10)
        integer(max_postings, "max_total_postings", 1, 10**11)
        self.until = time.monotonic() + deadline_seconds
        self.minimum = min_free_bytes
        self.directories = (output, Path(os.environ.get("SQLITE_TMPDIR") or tempfile.gettempdir()))
        self.max_candidates, self.max_postings = max_candidates, max_postings
        self.candidates = self.postings = 0
        self.last_disk_check = 0
        self.failure = None

    def check(self, force_disk=False):
        if self.failure is not None:
            raise self.failure
        now = time.monotonic()
        if now >= self.until:
            raise LimitError("Benchmark exceeded its total deadline")
        if force_disk or now - self.last_disk_check >= 1:
            for directory in self.directories:
                if directory.exists():
                    check_free(directory, self.minimum)
            self.last_disk_check = now

    def progress(self):
        try:
            self.check()
            return 0
        except (InputError, OSError) as error:
            self.failure = error
            return 1

    def install(self, connection):
        connection.set_progress_handler(self.progress, 10000)

    def candidate_options(self):
        self.check()
        candidates = self.max_candidates - self.candidates
        postings = self.max_postings - self.postings
        if min(candidates, postings) <= 0:
            raise LimitError("Benchmark exhausted its total candidate/posting budget")
        return {"max_candidates": min(2_000_000, candidates), "max_postings": min(10_000_000, postings),
                "deadline_seconds": min(300, self.until - time.monotonic())}

    def record(self, stats):
        self.candidates += stats["candidates"]
        self.postings += stats["postings"]
        self.check()


def accession_fold(stem, seed="20261001"):
    """Keep all versions of an accession together; fold zero is held out."""
    return int.from_bytes(hashlib.sha256((seed + "\0assembly\0" + stem).encode()).digest()[:8], "big") % 5


def _sample_key(sha, seed):
    return hashlib.sha256((seed + "\0sequence\0" + sha).encode()).hexdigest()


def _scalar(connection, sql, parameters=()):
    return connection.execute(sql, parameters).fetchone()[0]


def _fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def _distance_backend():
    try:
        from rapidfuzz.distance.Levenshtein import distance
    except ImportError:
        return levenshtein, "crisprmap.core.levenshtein"
    return distance, "rapidfuzz==" + version("rapidfuzz")


def _both_distance(query, target, distance):
    forward = distance(query, target)
    reverse = reverse_complement(target)
    return min(forward, distance(query, reverse)) if target != reverse else forward


def _oracle_distance(query, target):
    """Independent two-row dynamic programming, without the production kernel."""
    previous = list(range(len(target) + 1))
    for i, a in enumerate(query, 1):
        current = [i]
        for j, b in enumerate(target, 1):
            current.append(min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (a != b)))
        previous = current
    return previous[-1]


def _cohorts(connection, seed, maximum):
    """Materialize full source-to-fold membership and leakage sets on disk."""
    connection.create_function("benchmark_fold", 1, lambda value: accession_fold(value, seed), deterministic=True)
    connection.create_function("benchmark_sample", 1, lambda value: _sample_key(value, seed), deterministic=True)
    connection.executescript("""
        CREATE TEMP TABLE bm_assemblies(assembly_id INTEGER PRIMARY KEY,heldout INTEGER NOT NULL,qualified INTEGER NOT NULL);
        INSERT INTO bm_assemblies SELECT assembly_id,benchmark_fold(accession_stem)=0,
          CASE WHEN status='complete' AND validation_status='valid' THEN 1 ELSE 0 END FROM catalog.assemblies;
        CREATE INDEX bm_assemblies_qualified ON bm_assemblies(qualified,assembly_id);
        CREATE TEMP TABLE bm_shared_contigs(sequence_sha256 TEXT PRIMARY KEY) WITHOUT ROWID;
        INSERT INTO bm_shared_contigs
          SELECT c.sequence_sha256 FROM catalog.contigs c JOIN bm_assemblies f USING(assembly_id)
          WHERE f.qualified=1 AND c.sequence_sha256 IS NOT NULL GROUP BY c.sequence_sha256 HAVING min(f.heldout)<>max(f.heldout);
        CREATE TEMP TABLE bm_contigs(contig_id INTEGER PRIMARY KEY,heldout INTEGER NOT NULL);
        INSERT INTO bm_contigs SELECT c.contig_id,f.heldout FROM catalog.contigs c
          JOIN bm_assemblies f USING(assembly_id) WHERE f.qualified=1 AND c.sequence_sha256 IS NOT NULL
          AND NOT EXISTS(SELECT 1 FROM bm_shared_contigs x WHERE x.sequence_sha256=c.sequence_sha256);
        CREATE INDEX bm_contigs_fold ON bm_contigs(heldout,contig_id);
        CREATE TEMP TABLE bm_train(node_id INTEGER PRIMARY KEY);
        INSERT OR IGNORE INTO bm_train SELECT a.consensus_node_id FROM bm_contigs c
          JOIN catalog.arrays a USING(contig_id) JOIN catalog.nodes n ON n.node_id=a.consensus_node_id
          WHERE c.heldout=0 AND a.category='Bona-fide' AND n.eligible=1;
        CREATE TEMP TABLE bm_heldout_nodes(node_id INTEGER PRIMARY KEY);
        INSERT OR IGNORE INTO bm_heldout_nodes SELECT a.consensus_node_id FROM bm_contigs c
          JOIN catalog.arrays a USING(contig_id) WHERE c.heldout=1 AND a.consensus_node_id IS NOT NULL;
        INSERT OR IGNORE INTO bm_heldout_nodes SELECT r.node_id FROM bm_contigs c
          JOIN catalog.arrays a USING(contig_id) JOIN catalog.repeat_occurrences r USING(array_id)
          WHERE c.heldout=1 AND r.node_id IS NOT NULL;
        CREATE TEMP TABLE bm_queries(node_id INTEGER PRIMARY KEY,occurrences INTEGER NOT NULL);
        INSERT INTO bm_queries SELECT a.consensus_node_id,count(*) FROM bm_contigs c
          JOIN catalog.arrays a USING(contig_id) JOIN catalog.nodes n ON n.node_id=a.consensus_node_id
          WHERE c.heldout=1 AND a.category='Bona-fide' AND n.eligible=1 GROUP BY a.consensus_node_id;
        CREATE TEMP TABLE bm_sample(node_id INTEGER PRIMARY KEY,rank_key TEXT NOT NULL);
    """)
    connection.execute("INSERT INTO bm_sample SELECT q.node_id,benchmark_sample(n.sequence_sha256) "
                       "FROM bm_queries q JOIN catalog.nodes n USING(node_id) ORDER BY 2,n.sequence_sha256 LIMIT ?", (maximum,))
    connection.executescript("""
        CREATE TEMP TABLE bm_anchor_evidence(node_id INTEGER PRIMARY KEY);
        INSERT INTO bm_anchor_evidence SELECT node_id FROM bm_sample;
        INSERT OR IGNORE INTO bm_anchor_evidence SELECT r.node_id FROM bm_sample s
          JOIN catalog.arrays a ON a.consensus_node_id=s.node_id JOIN bm_contigs c USING(contig_id)
          JOIN catalog.repeat_occurrences r USING(array_id)
          WHERE c.heldout=1 AND r.node_id IS NOT NULL;
        CREATE TEMP TABLE bm_anchors(node_id INTEGER PRIMARY KEY);
        INSERT INTO bm_anchors SELECT e.node_id FROM bm_anchor_evidence e JOIN catalog.nodes n USING(node_id) WHERE n.eligible=1;
        CREATE TEMP TABLE bm_purge(node_id INTEGER PRIMARY KEY);
        CREATE TEMP TABLE bm_matches(query_node_id INTEGER NOT NULL,node_id INTEGER NOT NULL,
          distance INTEGER NOT NULL,PRIMARY KEY(query_node_id,node_id)) WITHOUT ROWID;
        CREATE TEMP TABLE bm_oracle(node_id INTEGER PRIMARY KEY);
    """)
    counts = {
        "all_assemblies": _scalar(connection, "SELECT count(*) FROM bm_assemblies"),
        "heldout_assemblies_before_contig_exclusions": _scalar(connection, "SELECT count(*) FROM bm_assemblies WHERE heldout=1"),
        "training_assemblies_before_contig_exclusions": _scalar(connection, "SELECT count(*) FROM bm_assemblies WHERE heldout=0"),
        "failed_or_unavailable_assemblies": _scalar(connection, "SELECT count(*) FROM catalog.assemblies WHERE status<>'complete'"),
        "excluded_nonvalidated_assemblies": _scalar(connection, "SELECT count(*) FROM bm_assemblies WHERE qualified=0"),
        "excluded_nonvalidated_contigs": _scalar(connection, "SELECT count(*) FROM bm_assemblies f JOIN catalog.contigs c USING(assembly_id) WHERE f.qualified=0"),
        "excluded_nonvalidated_arrays": _scalar(connection, "SELECT count(*) FROM bm_assemblies f JOIN catalog.contigs c USING(assembly_id) JOIN catalog.arrays a USING(contig_id) WHERE f.qualified=0"),
        "excluded_nonvalidated_repeat_occurrences": _scalar(connection, "SELECT count(*) FROM bm_assemblies f JOIN catalog.contigs c USING(assembly_id) JOIN catalog.arrays a USING(contig_id) JOIN catalog.repeat_occurrences r USING(array_id) WHERE f.qualified=0"),
        "qualified_training_assemblies_before_contig_exclusions": _scalar(connection, "SELECT count(*) FROM bm_assemblies WHERE qualified=1 AND heldout=0"),
        "qualified_heldout_assemblies_before_contig_exclusions": _scalar(connection, "SELECT count(*) FROM bm_assemblies WHERE qualified=1 AND heldout=1"),
        "cross_fold_contig_hashes": _scalar(connection, "SELECT count(*) FROM bm_shared_contigs"),
        "excluded_cross_fold_contigs": _scalar(connection, "SELECT count(*) FROM catalog.contigs c JOIN bm_shared_contigs s USING(sequence_sha256) JOIN bm_assemblies f USING(assembly_id) WHERE f.qualified=1"),
        "excluded_missing_hash_contigs": _scalar(connection, "SELECT count(*) FROM catalog.contigs c JOIN bm_assemblies f USING(assembly_id) WHERE f.qualified=1 AND c.sequence_sha256 IS NULL"),
        "retained_training_contigs": _scalar(connection, "SELECT count(*) FROM bm_contigs WHERE heldout=0"),
        "retained_heldout_contigs": _scalar(connection, "SELECT count(*) FROM bm_contigs WHERE heldout=1"),
        "heldout_arrays_all_tiers": _scalar(connection, "SELECT count(*) FROM bm_contigs c JOIN catalog.arrays a USING(contig_id) WHERE c.heldout=1"),
        "heldout_observed_repeats_all_tiers": _scalar(connection, "SELECT count(*) FROM bm_contigs c JOIN catalog.arrays a USING(contig_id) JOIN catalog.repeat_occurrences r USING(array_id) WHERE c.heldout=1"),
        "training_core_nodes": _scalar(connection, "SELECT count(*) FROM bm_train"),
        "heldout_all_evidence_nodes": _scalar(connection, "SELECT count(*) FROM bm_heldout_nodes"),
        "full_cohort_exact_rc_overlapping_training_nodes": _scalar(connection, "SELECT count(*) FROM bm_train JOIN bm_heldout_nodes USING(node_id)"),
        "heldout_core_query_nodes": _scalar(connection, "SELECT count(*) FROM bm_queries"),
        "heldout_core_query_occurrences": _scalar(connection, "SELECT coalesce(sum(occurrences),0) FROM bm_queries"),
        "sampled_query_nodes": _scalar(connection, "SELECT count(*) FROM bm_sample"),
        "sample_anchor_evidence_nodes": _scalar(connection, "SELECT count(*) FROM bm_anchor_evidence"),
        "heldout_core_arrays_before_sequence_eligibility": _scalar(connection, "SELECT count(*) FROM bm_contigs c JOIN catalog.arrays a USING(contig_id) WHERE c.heldout=1 AND a.category='Bona-fide'"),
        "heldout_missing_consensus_sequences": _scalar(connection, "SELECT count(*) FROM bm_contigs c JOIN catalog.arrays a USING(contig_id) WHERE c.heldout=1 AND a.consensus_node_id IS NULL"),
        "heldout_missing_observed_repeat_sequences": _scalar(connection, "SELECT count(*) FROM bm_contigs c JOIN catalog.arrays a USING(contig_id) JOIN catalog.repeat_occurrences r USING(array_id) WHERE c.heldout=1 AND r.node_id IS NULL"),
        "sample_scoped_homology_anchors": _scalar(connection, "SELECT count(*) FROM bm_anchors"),
        "assemblies_with_tax_id": _scalar(connection, "SELECT count(*) FROM catalog.assemblies WHERE tax_id IS NOT NULL"),
        "cross_fold_tax_ids": _scalar(connection, "SELECT count(*) FROM (SELECT a.tax_id FROM catalog.assemblies a JOIN bm_assemblies f USING(assembly_id) WHERE a.tax_id IS NOT NULL GROUP BY a.tax_id HAVING min(f.heldout)<>max(f.heldout))"),
    }
    counts["training_nodes_after_full_exact_rc_purge"] = counts["training_core_nodes"] - counts["full_cohort_exact_rc_overlapping_training_nodes"]
    counts["sample_anchor_nodes_excluded_as_ineligible"] = counts["sample_anchor_evidence_nodes"] - counts["sample_scoped_homology_anchors"]
    counts["heldout_core_arrays_excluded_as_missing_or_ineligible"] = counts["heldout_core_arrays_before_sequence_eligibility"] - counts["heldout_core_query_occurrences"]
    counts["query_nodes_not_sampled"] = counts["heldout_core_query_nodes"] - counts["sampled_query_nodes"]
    return counts


def _verified_matches(reader, sequence, distance, budget):
    try:
        for node in reader.iter_candidates(sequence, max_distance=None, min_similarity=0.9, core_only=True,
                                           **budget.candidate_options()):
            budget.check()
            observed = _both_distance(sequence, node["sequence"], distance)
            if observed <= distance_threshold(len(sequence), node["length"], SIMILARITY):
                yield node["node_id"], observed
        budget.record(reader.last_stats)
    finally:
        # CandidateIndex installs a shorter per-query handler while generating
        # candidates, then clears it. Restore the total-run SQL watchdog.
        budget.install(reader.connection)


def _metrics(connection, sample, name, exclusion=""):
    rows = []
    for node_id, sequence, sha, occurrences in sample:
        best, hits, exact = connection.execute(
            "SELECT min(m.distance),count(*),coalesce(sum(m.distance=0),0) FROM bm_matches m "
            "WHERE m.query_node_id=? " + exclusion, (node_id,)).fetchone()
        rows.append({"query_sequence_sha256": sha, "heldout_array_occurrences": occurrences,
                     "covered": bool(hits), "exact_rc_hit": bool(exact),
                     "best_distance_within_90_percent": best, "qualifying_reference_nodes": hits})
    covered = sum(row["covered"] for row in rows)
    exact = sum(row["exact_rc_hit"] for row in rows)
    return {"name": name, "evaluated_unique_sequences": len(rows), "covered_sequences": covered,
            "coverage": _fraction(covered, len(rows)), "exact_rc_sequences": exact,
            "near_only_sequences": covered - exact, "no_hit_sequences": len(rows) - covered,
            "rows": rows}


def _legacy_coverage(sample, path, distance, budget):
    if path is None:
        return {"status": "not_run", "reason": "No historical reference supplied"}
    reference, file_hash = load_reference(Path(path))
    sequences = sorted({canonical_sequence(row["sequence"]) for row in reference["entries"]
                        if set(row["sequence"]) <= set("ACGT")})
    rows = []
    for _, query, sha, _ in sample:
        best = None
        hits = exact = 0
        for target in sequences:
            budget.check()
            limit = distance_threshold(len(query), len(target), SIMILARITY)
            if abs(len(query) - len(target)) > limit:
                continue
            observed = _both_distance(query, target, distance)
            if observed <= limit:
                hits += 1
                exact += observed == 0
                best = observed if best is None else min(best, observed)
        rows.append({"query_sequence_sha256": sha, "covered": bool(hits), "exact_rc_hit": bool(exact),
                     "best_distance_within_90_percent": best, "qualifying_reference_nodes": hits})
    covered = sum(row["covered"] for row in rows)
    return {"status": "complete", "namespace": reference["provenance"]["namespace"],
            "reference_file_sha256": file_hash, "reference_content_sha256": reference["content_sha256"],
            "eligible_canonical_sequences": len(sequences), "evaluated_unique_sequences": len(rows),
            "covered_sequences": covered, "coverage": _fraction(covered, len(rows)),
            "exact_rc_sequences": sum(row["exact_rc_hit"] for row in rows), "rows": rows,
            "interpretation": "Frozen historical-reference coverage; its training overlap with the heldout genomes is unknown. Labels are not scored."}


def _oracle(reader, sample, budget):
    """A bounded reference-subset oracle; no claim of full-corpus DP validation."""
    connection = reader.connection
    connection.execute("INSERT OR IGNORE INTO bm_oracle SELECT node_id FROM bm_train "
                       "JOIN catalog.nodes USING(node_id) ORDER BY benchmark_sample(sequence_sha256),sequence_sha256 LIMIT ?", (ORACLE_REFERENCES,))
    # Include known positives as well as a deterministic background, so a rare-hit
    # corpus cannot yield a vacuous all-negative engineering check.
    for node_id, *_ in sample[:ORACLE_QUERIES]:
        connection.execute("INSERT OR IGNORE INTO bm_oracle SELECT node_id FROM bm_matches "
                           "WHERE query_node_id=? ORDER BY node_id LIMIT 32", (node_id,))
    references = list(connection.execute("SELECT n.node_id,n.sequence FROM bm_oracle o JOIN catalog.nodes n USING(node_id) ORDER BY n.node_id"))
    mismatches = []
    pairs = positives = indexed_positives = 0
    for node_id, sequence, sha, _ in sample[:ORACLE_QUERIES]:
        for target_id, target in references:
            budget.check()
            pairs += 1
            target_rc = target.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            observed = min(_oracle_distance(sequence, target), _oracle_distance(sequence, target_rc))
            expected = 10 * observed <= max(len(sequence), len(target))
            actual_row = connection.execute("SELECT distance FROM bm_matches WHERE query_node_id=? AND node_id=?", (node_id, target_id)).fetchone()
            actual = actual_row is not None
            positives += expected
            indexed_positives += actual
            if actual != expected or (actual and actual_row[0] != observed):
                mismatches.append({"query_sequence_sha256": sha,
                                   "target_sequence_sha256": hashlib.sha256(target.encode()).hexdigest(),
                                   "expected_qualifying": expected, "indexed_qualifying": actual,
                                   "oracle_distance": observed, "indexed_distance": actual_row[0] if actual else None})
    return {"status": ("passed" if pairs else "not_run") if not mismatches else "failed", "scope": "bounded training-reference subset, not full-corpus exhaustive alignment",
            "reason": (None if pairs else ("no_eligible_heldout_queries" if not sample else "no_eligible_training_references")),
            "query_count": min(len(sample), ORACLE_QUERIES), "reference_count": len(references),
            "compared_pairs": pairs, "oracle_qualifying_pairs": positives,
            "indexed_qualifying_pairs": indexed_positives, "mismatches": mismatches,
            "sampling": "512 smallest seeded sequence hashes plus at most 32 indexed qualifying references per query"}


def _report(metrics):
    audit = metrics["overlap_audit"]
    lines = ["# Modern repeat reference benchmark", "", "This measures retrieval coverage, not biological family accuracy.", "",
             f"Unique heldout Bona-fide consensus queries: {audit['sampled_query_nodes']} sampled from {audit['heldout_core_query_nodes']}.",
             f"Training core nodes: {audit['training_core_nodes']}; {audit['full_cohort_exact_rc_overlapping_training_nodes']} overlap the full heldout consensus/variant evidence set.",
             f"Excluded contigs: {audit['excluded_nonvalidated_contigs']} from assemblies without complete/valid outcomes, {audit['excluded_cross_fold_contigs']} with cross-fold hashes, and {audit['excluded_missing_hash_contigs']} without hashes.", "",
             "| Evaluation | Queries | Covered | Exact/RC | Coverage |", "|---|---:|---:|---:|---:|"]
    modes = list(metrics["evaluations"].values())
    if metrics["legacy"]["status"] == "complete":
        modes.append({"name": "Frozen historical reference (overlap unknown)", **metrics["legacy"]})
    for row in modes:
        percent = "NA" if row["coverage"] is None else f"{100 * row['coverage']:.2f}%"
        lines.append(f"| {row['name']} | {row['evaluated_unique_sequences']} | {row['covered_sequences']} | {row['exact_rc_sequences']} | {percent} |")
    lines += ["", "Genome holdout keeps accession versions together and admits only complete, source-validated assemblies. References contain only eligible Bona-fide training consensuses; query evidence comes from retained heldout contigs.", "",
              "Exact/RC purge uses every retained heldout consensus and observed repeat from both source tiers, including unsampled arrays. Homology90 stress additionally purges every training node meeting the same 90% similarity cutoff against sampled queries and their associated observed repeat variants. Zero coverage is expected by construction; this is a sample-scoped policy stress test, not a full-cohort homology-independent split.", "",
              "Taxonomy is audited but not used for splitting: this is not a species-heldout evaluation. Upstream detector training membership is unknown. Historical-reference overlap is unknown, so coverage differences do not establish predictive superiority.", "",
              f"Independent dynamic-programming oracle: {metrics['oracle']['status']}, {metrics['oracle']['compared_pairs']} pairs in a bounded reference subset. This does not replace full index-inventory validation.", ""]
    return "\n".join(lines)


def run_benchmark(catalog_dir, index_dir, output, legacy_reference=None, max_queries=1000, seed="20261001",
                  *, deadline_seconds=172800, min_free_bytes=2 * 1024**3,
                  max_total_candidates=100_000_000, max_total_postings=1_000_000_000):
    """Create a new benchmark directory, publishing its completion marker last."""
    integer(max_queries, "max_queries", 1, 10000)
    if not isinstance(seed, str) or not 1 <= len(seed) <= 100 or any(ord(c) < 32 or ord(c) > 126 for c in seed):
        raise InputError("Seed must be bounded printable text")
    output = Path(output).absolute()
    if any(p.is_symlink() for p in (output, *output.parents)):
        raise InputError("Benchmark output must not traverse symlinks")
    budget = _RunBudget(output, deadline_seconds, min_free_bytes, max_total_candidates, max_total_postings)
    started = time.monotonic()
    distance, backend = _distance_backend()
    owned = False
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.mkdir()
        owned = True
        budget.check(force_disk=True)
        with CandidateIndex(catalog_dir, index_dir) as reader:
            connection = reader.connection
            budget.install(connection)
            budget.check()
            print(json.dumps({"event": "benchmark_full_cohort_audit_started"}), flush=True)
            audit = _cohorts(connection, seed, max_queries)
            print(json.dumps({"event": "benchmark_full_cohort_audit_complete", "counts": audit}, sort_keys=True), flush=True)
            if audit["sample_scoped_homology_anchors"] > MAX_ANCHORS:
                raise LimitError("Sample-associated homology anchors exceed 100000; reduce max_queries")
            sample = [tuple(row) for row in connection.execute(
                "SELECT s.node_id,n.sequence,n.sequence_sha256,q.occurrences FROM bm_sample s "
                "JOIN catalog.nodes n USING(node_id) JOIN bm_queries q USING(node_id) ORDER BY s.rank_key,n.sequence_sha256")]
            retrieval_started = time.monotonic()
            for completed, (node_id, sequence, _, _) in enumerate(sample, 1):
                for target_id, observed in _verified_matches(reader, sequence, distance, budget):
                    connection.execute("INSERT OR IGNORE INTO bm_matches SELECT ?,node_id,? FROM bm_train WHERE node_id=?",
                                       (node_id, observed, target_id))
                if completed % 25 == 0 or completed == len(sample):
                    print(json.dumps({"event": "benchmark_retrieval_progress", "completed": completed, "total": len(sample)}), flush=True)
            retrieval_seconds = time.monotonic() - retrieval_started
            stress_started = time.monotonic()
            anchors = connection.execute("SELECT a.node_id,n.sequence FROM bm_anchors a JOIN catalog.nodes n USING(node_id) ORDER BY a.node_id")
            for completed, (node_id, sequence) in enumerate(anchors, 1):
                if _scalar(connection, "SELECT count(*) FROM bm_sample WHERE node_id=?", (node_id,)):
                    connection.execute("INSERT OR IGNORE INTO bm_purge SELECT node_id FROM bm_matches WHERE query_node_id=?", (node_id,))
                else:
                    for target_id, _ in _verified_matches(reader, sequence, distance, budget):
                        connection.execute("INSERT OR IGNORE INTO bm_purge SELECT node_id FROM bm_train WHERE node_id=?", (target_id,))
                if completed % 25 == 0 or completed == audit["sample_scoped_homology_anchors"]:
                    print(json.dumps({"event": "benchmark_homology_stress_progress", "completed": completed, "total": audit["sample_scoped_homology_anchors"]}), flush=True)
            stress_seconds = time.monotonic() - stress_started
            audit["sample_homology_purged_training_nodes"] = _scalar(connection, "SELECT count(*) FROM bm_purge")
            audit["training_nodes_after_exact_and_sample_homology_purge"] = _scalar(connection,
                "SELECT count(*) FROM bm_train t WHERE NOT EXISTS(SELECT 1 FROM bm_heldout_nodes h WHERE h.node_id=t.node_id) "
                "AND NOT EXISTS(SELECT 1 FROM bm_purge p WHERE p.node_id=t.node_id)")
            exact_exclusion = "AND NOT EXISTS(SELECT 1 FROM bm_heldout_nodes h WHERE h.node_id=m.node_id)"
            evaluations = {
                "genome_heldout": _metrics(connection, sample, "Genome-heldout coverage"),
                "full_exact_rc_purged": _metrics(connection, sample, "Full-cohort exact/RC-purged coverage", exact_exclusion),
                "sample_homology90_stress": _metrics(connection, sample, "Sample-scoped homology90 purge stress", exact_exclusion +
                    " AND NOT EXISTS(SELECT 1 FROM bm_purge p WHERE p.node_id=m.node_id)"),
            }
            oracle_started = time.monotonic()
            oracle = _oracle(reader, sample, budget)
            oracle_seconds = time.monotonic() - oracle_started
            if oracle["mismatches"]:
                raise InputError("Independent benchmark oracle disagrees with indexed retrieval")
            if evaluations["sample_homology90_stress"]["covered_sequences"]:
                raise InputError("Homology stress retained a forbidden qualifying neighbor")
            legacy_started = time.monotonic()
            legacy = _legacy_coverage(sample, legacy_reference, distance, budget)
            legacy_seconds = time.monotonic() - legacy_started
            metrics = {"schema": dict(SCHEMA), "status": "complete", "seed": seed,
                       "catalog_database_sha256": reader.catalog_manifest["database"]["sha256"],
                       "index_database_sha256": reader.manifest["database"]["sha256"],
                       "analysis_profile_sha256": reader.catalog_manifest["profile_sha256"],
                       "policy": {"similarity": "1-global_edit_distance/max(query_length,reference_length)",
                                  "minimum_similarity_fraction": [9, 10], "orientations": "both",
                                  "training_reference": "eligible Bona-fide consensus nodes with retained training-contig occurrences",
                                  "assembly_eligibility": "status=complete AND validation_status=valid",
                                  "fold": "sha256(seed+NUL+assembly+NUL+accession_stem),first8bytes_big_endian_mod5;zero heldout",
                                  "sampling_unit": "unique canonical heldout Bona-fide consensus sequence",
                                  "sampling": "smallest SHA256(seed+NUL+sequence+NUL+sequence_sha256), tie sequence_sha256",
                                  "max_queries": max_queries, "taxonomy_used_for_split": False,
                                  "full_exact_rc_audit": True, "full_cohort_homology_audit": False,
                                  "homology_stress_scope": "sampled query consensuses and all eligible observed variants of retained heldout arrays bearing those consensuses",
                                  "upstream_detector_training_overlap": "unknown", "biological_family_accuracy": "unavailable"},
                       "overlap_audit": audit, "evaluations": evaluations,
                       "legacy": legacy, "oracle": oracle,
                       "budgets": {"deadline_seconds": deadline_seconds, "min_free_bytes": min_free_bytes,
                                   "max_total_candidates": max_total_candidates, "max_total_postings": max_total_postings,
                                   "retrieved_candidates": budget.candidates, "visited_postings": budget.postings},
                       "runtime": {"distance_backend": backend, "retrieval_seconds": retrieval_seconds,
                                   "sample_homology_stress_seconds": stress_seconds, "bounded_oracle_seconds": oracle_seconds,
                                   "legacy_reference_coverage_seconds": legacy_seconds,
                                   "comparison_scope": "Different workloads; no original-tool runtime or speedup is estimated",
                                   "total_seconds_before_publication": time.monotonic() - started}}
            budget.check(force_disk=True)
            atomic_json(output / "metrics.json", metrics)
            (output / "REPORT.md").write_text(_report(metrics), encoding="utf-8")
            manifest = {"schema": dict(SCHEMA), "status": "complete",
                        "catalog_database_sha256": metrics["catalog_database_sha256"],
                        "index_database_sha256": metrics["index_database_sha256"],
                        "files": {name: {"sha256": file_sha256(output / name), "bytes": (output / name).stat().st_size}
                                  for name in ("metrics.json", "REPORT.md")},
                        "evaluated_unique_sequences": len(sample), "oracle_status": oracle["status"]}
            budget.check(force_disk=True)
            atomic_json(output / "manifest.json", manifest)
            return manifest
    except sqlite3.OperationalError as error:
        if budget.failure is not None:
            raise budget.failure from error
        raise InputError("Benchmark SQL operation failed") from error
    finally:
        if owned and not (output / "manifest.json").exists():
            shutil.rmtree(output)
