"""Small synthetic leakage fixtures; production runs belong on source compute."""
import gzip
import hashlib
import json
import sqlite3

import pytest

from crisprmap.catalog import build_partition
from crisprmap.catalog_index import CandidateIndex, build_index, file_sha256
from crisprmap.core import reverse_complement
from crisprmap.io import InputError, LimitError, atomic_json
from crisprmap.modern_benchmark import (_cohorts, _oracle_distance, accession_fold, run_benchmark)
from crisprmap.reference import build_reference

Q1 = "AAAAACCCCCGGGGGTTTTT"
VARIANT = "AAAAACCCCCGGGGGTTTTA"
NEAR = "AAAAACCCCCGGGGGTTTAA"
Q2 = "ACGTACGTACGTACGTACGT"
Q3 = "AGAGAGAGAGAGAGAGAGAG"
EXCLUDED = "AGCTAGCTTTAGCTAGCTTT"
PROFILE = "1" * 64


def fixture_catalog(tmp_path):
    stems = {True: [], False: []}
    for i in range(1, 1000):
        stem = f"GCF_{i:09d}"
        stems[accession_fold(stem) == 0].append(stem)
        if min(map(len, stems.values())) >= 4:
            break
    sources = []

    def source(stem, arrays, contig_hash, label, accession_version=1):
        accession = f"{stem}.{accession_version}"
        calls = [{"call_id": str(i), "interval": {"start": 1 + 100 * i, "end": 70 + 100 * i},
                  "strand": None, "repeat_consensus": sequence, "repeat_count": len(variants), "spacer_count": 0,
                  "confidence": {"category": category},
                  "repeats": [{"index": j, "sequence": variant,
                               "interval": {"start": 1 + 100 * i + 25 * j, "end": len(variant) + 100 * i + 25 * j}}
                              for j, variant in enumerate(variants, 1)]}
                 for i, (sequence, category, variants) in enumerate(arrays)]
        payload = {"assembly": {"accession": accession, "ncbi_metadata": {"organism": {"taxId": 123}}},
                   "status": "complete", "validation": {"status": "valid"},
                   "provenance": {"analysis_profile_sha256": PROFILE}, "summary": {"array_count": len(calls)},
                   "sequences": [{"id": "contig", "length": 10000, "sha256": contig_hash,
                                  "tools": [{"name": "CRISPRidentify", "version": "2.3.2", "arrays": calls}]}]}
        path = tmp_path / f"{label}.json.gz"
        path.write_bytes(gzip.compress(json.dumps(payload).encode(), mtime=0))
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        sources.append({"kind": "raw", "path": str(path), "sha256": sha,
                        "entries": [{"accession": accession, "status": "complete", "sha256": sha}]})

    source(stems[False][0], [(s, "Bona-fide", [s]) for s in [VARIANT, NEAR, Q2, "C" * 20]] +
           [(Q3, "Possible", [Q3])], "a" * 64, "training")
    source(stems[True][0], [(Q1, "Bona-fide", [Q1]), (Q2, "Bona-fide", [reverse_complement(Q2)]),
                          (Q3, "Bona-fide", [Q3]), ("AT" * 10, "Possible", [reverse_complement(VARIANT)])],
           "b" * 64, "heldout")
    source(stems[True][0], [(Q1, "Bona-fide", [Q1])], "c" * 64, "heldout-version", 2)
    source(stems[False][1], [(EXCLUDED, "Bona-fide", [EXCLUDED])], "d" * 64, "shared-training")
    source(stems[True][1], [(EXCLUDED, "Bona-fide", [EXCLUDED])], "d" * 64, "shared-heldout")
    source(stems[True][2], [("GATTACA" * 3, "Bona-fide", ["GATTACA" * 3])], None, "missing-hash")
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"profile_sha256": PROFILE, "sources": sources}))
    catalog = tmp_path / "catalog"
    build_partition(plan, catalog, PROFILE, min_free_bytes=0)
    index = tmp_path / "index"
    build_index(catalog, index, min_free_bytes=0)
    return catalog, index


def test_full_variant_purge_contig_exclusions_and_core_only_membership(tmp_path):
    catalog, index = fixture_catalog(tmp_path)
    output = tmp_path / "benchmark"
    result = run_benchmark(catalog, index, output, min_free_bytes=0)
    metrics = json.loads((output / "metrics.json").read_text())
    audit = metrics["overlap_audit"]
    assert result["evaluated_unique_sequences"] == 3
    assert audit["excluded_cross_fold_contigs"] == 2
    assert audit["excluded_missing_hash_contigs"] == 1
    assert audit["training_core_nodes"] == 4
    assert audit["heldout_core_query_occurrences"] == 4
    assert audit["full_cohort_exact_rc_overlapping_training_nodes"] == 2
    assert audit["training_nodes_after_full_exact_rc_purge"] == 2
    assert audit["cross_fold_tax_ids"] == 1
    assert metrics["policy"]["taxonomy_used_for_split"] is False
    modes = metrics["evaluations"]
    assert modes["genome_heldout"]["covered_sequences"] == 2
    assert modes["genome_heldout"]["exact_rc_sequences"] == 1
    assert modes["full_exact_rc_purged"]["covered_sequences"] == 1
    assert modes["sample_homology90_stress"]["covered_sequences"] == 0
    assert metrics["oracle"]["status"] == "passed"
    assert metrics["oracle"]["oracle_qualifying_pairs"] > 0
    assert str(tmp_path) not in (output / "metrics.json").read_text()
    assert "family accuracy" in (output / "REPORT.md").read_text()


def test_unsampled_heldout_variants_still_purge_training_nodes(tmp_path):
    catalog, index = fixture_catalog(tmp_path)
    with CandidateIndex(catalog, index) as reader:
        counts = _cohorts(reader.connection, "20261001", 1)
        assert counts["sampled_query_nodes"] == 1
        assert counts["full_cohort_exact_rc_overlapping_training_nodes"] == 2
        assert counts["query_nodes_not_sampled"] == 2


def test_optional_historical_coverage_has_separate_overlap_caveat(tmp_path):
    catalog, index = fixture_catalog(tmp_path)
    legacy = build_reference([{"id": "old", "sequence": Q3, "description": "synthetic historical record"}], {},
                             {"namespace": "synthetic-v1", "source_url": "https://example.invalid/fixture", "source_sha256": "2" * 64})
    path = tmp_path / "legacy.json"
    atomic_json(path, legacy)
    output = tmp_path / "benchmark"
    run_benchmark(catalog, index, output, legacy_reference=path, min_free_bytes=0)
    metrics = json.loads((output / "metrics.json").read_text())
    assert metrics["legacy"]["covered_sequences"] == 1
    assert metrics["legacy"]["exact_rc_sequences"] == 1
    assert "unknown" in metrics["legacy"]["interpretation"]
    before = (output / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(catalog, index, output, min_free_bytes=0)
    assert (output / "manifest.json").read_bytes() == before


def test_sample_is_deterministic_and_not_all_occurrences(tmp_path):
    catalog, index = fixture_catalog(tmp_path)
    results = []
    for label in ("one", "two"):
        output = tmp_path / label
        run_benchmark(catalog, index, output, max_queries=2, min_free_bytes=0)
        metrics = json.loads((output / "metrics.json").read_text())
        results.append(metrics["evaluations"]["genome_heldout"]["rows"])
        assert metrics["overlap_audit"]["query_nodes_not_sampled"] == 1
    assert results[0] == results[1]


def test_oracle_distance_and_invalid_limits(tmp_path):
    assert _oracle_distance("", "ACGT") == 4
    assert _oracle_distance("AACG", "ACG") == 1
    assert _oracle_distance("ACGT", "AGGT") == 1
    with pytest.raises(InputError):
        run_benchmark(tmp_path, tmp_path, tmp_path / "output", max_queries=0)
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("deadline", [0, -1, True, float("nan"), float("inf")])
def test_invalid_deadline_rejected_before_output(tmp_path, deadline):
    with pytest.raises(InputError, match="deadline"):
        run_benchmark(tmp_path, tmp_path, tmp_path / "output", deadline_seconds=deadline)
    assert not (tmp_path / "output").exists()


def test_total_candidate_limit_cleans_unpublished_result(tmp_path):
    catalog, index = fixture_catalog(tmp_path)
    output = tmp_path / "benchmark"
    with pytest.raises(LimitError, match="budget"):
        run_benchmark(catalog, index, output, max_total_candidates=1, min_free_bytes=0)
    assert not output.exists()


@pytest.mark.parametrize("status,validation", [("failed", "valid"), ("complete", "invalid"), ("complete", None)])
def test_failed_or_unvalidated_assembly_features_do_not_enter_cohorts(tmp_path, status, validation):
    catalog, _ = fixture_catalog(tmp_path)
    database = catalog / "catalog.sqlite3"
    with sqlite3.connect(database) as connection:
        training_id = next(row[0] for row in connection.execute("SELECT assembly_id,accession_stem FROM assemblies")
                           if accession_fold(row[1]) != 0 and connection.execute(
                               "SELECT count(*) FROM contigs c JOIN arrays a USING(contig_id) WHERE c.assembly_id=?", (row[0],)).fetchone()[0] == 5)
        connection.execute("UPDATE assemblies SET status=?,validation_status=? WHERE assembly_id=?", (status, validation, training_id))
    manifest = json.loads((catalog / "manifest.json").read_text())
    manifest["database"].update({"sha256": file_sha256(database), "bytes": database.stat().st_size})
    (catalog / "manifest.json").write_text(json.dumps(manifest))
    index = tmp_path / "filtered-index"
    build_index(catalog, index, min_free_bytes=0)
    with CandidateIndex(catalog, index) as reader:
        audit = _cohorts(reader.connection, "20261001", 1000)
        assert audit["excluded_nonvalidated_assemblies"] == 1
        assert audit["excluded_nonvalidated_contigs"] == 1
        assert audit["excluded_nonvalidated_arrays"] == 5
        assert audit["excluded_nonvalidated_repeat_occurrences"] == 5
        assert audit["training_core_nodes"] == 0
        assert audit["heldout_core_query_nodes"] == 3
