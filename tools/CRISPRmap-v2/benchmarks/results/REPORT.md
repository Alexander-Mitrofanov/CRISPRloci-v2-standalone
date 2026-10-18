# Benchmark against recovered original CRISPRmap evidence

Generated: 2026-09-30T23:49:15.701894+00:00

This is an original-output comparison, not a rerun of the original executable. The original service was offline and no local source release was recovered. Original runtime and speedup are unavailable.

Reference: `crisprmap-2013-v1`; 3,527 records. SHA256: `af386aed0c55f3c2976de02aec68eb6723478ef93ae4c8c2d8f89c42f8690108`.

## Original published 16-query example

| Label | Original known labels | v2 assigned | Assigned agreement |
| --- | ---: | ---: | ---: |
| Sequence family | 11 | 10 | 10/10 |
| Structural motif | 12 | 10 | 10/10 |

v2 statuses: `{"exact": 14, "near": 2}`. v2 elapsed: 0.179 seconds.

Known-label coverage is reported because abstention is not a correct positive assignment. Null-to-null matches are retained in JSON but are not counted as positive-label agreement. Original superclass inference is unsupported and excluded. Nonexact structural motifs are intentionally not assigned.

The two positive-label differences are explicit: Ecoli_K12 matches reference occurrences with competing families F2/F4 and motifs M1/M2, so v2 abstains; the original selected F4/M2. Stepidermidis has original motif M28 but is nonexact, and v2 does not infer its structural motif. This tool does not reproduce the original structure-aware assignment capability.

## Reference recovery

All 7,054 reference and reverse-complement queries: 7,054 exact matches. Annotation evidence was orientation-invariant for 3,527/3,527 pairs. This is an engineering self-recovery check, not independent validation.

## Sequence-group holdout

Evaluated 2031 canonical sequence/reverse-complement groups. Split: Five folds: SHA256('crisprmap-v2-heldout-2026-v1:' + canonical_sequence) modulo 5; every group evaluated once. Training record counts per fold: [2908, 2857, 2744, 2862, 2737]. Each query group and every occurrence of its reverse complement were removed from its training fold. Closely related sequences may occur on both sides; species, homology and family separation were not imposed.

Counts: `{"conflicting_or_partial_historical_label_queries": 56, "missing_historical_label_queries": 951, "queries_with_candidates": 810, "sole_candidate_agrees": 619, "sole_candidate_queries": 676, "true_label_in_candidates": 750, "unambiguous_label_queries": 1024}`.

Candidate coverage for unambiguous historical families: 0.791. True historical family present among candidates: 0.732. Agreement among single-candidate queries: 0.916.

Macro family metrics (equal family weight): `{"agreement_among_sole_candidates": 0.8628791180850678, "candidate_coverage": 0.727142736969777, "families_contributing_sole_candidate_agreement": 38, "families_with_unambiguous_query_labels": 40, "true_label_candidate_recall": 0.6550640631923665}`.

These are historical candidate-agreement metrics, not calibrated biological accuracy. The single-candidate metric requires complete unanimous candidate evidence; incomplete or conflicting reference evidence is excluded. Missing or contradictory query labels are excluded from those denominators and counted separately. The maximum edit distance was fixed before evaluation; no test-set tuning was performed. Per-family counts and per-query hashes are in the JSON.

## Family-absent stress test

Removed each of 40 historical families in turn, excluding complete duplicate/reverse-complement groups. Across 1024 unambiguous query groups, 218 retrieved evidence for another family (0.213). Query family assignments: 0; zero assigned labels follows from the nonexact-abstention design and is not evidence of biological specificity.

Cross-family candidates are not validated family assignments and may reflect genuine similarity across historical partitions. This stress test quantifies candidate retrieval when the source family is absent; it does not calibrate biological novelty detection.

## Synthetic controls and timing

Control outcomes: `{"homopolymer": {"no_hit": 4, "query_assignment_count": 0, "query_count": 4, "with_family_candidate_evidence": 0}, "random": {"no_hit": 100, "query_assignment_count": 0, "query_count": 100, "with_family_candidate_evidence": 0}, "shuffled": {"no_hit": 100, "query_assignment_count": 0, "query_count": 100, "with_family_candidate_evidence": 0}}`.

Synthetic sequences are engineering controls, not experimentally established non-CRISPR negatives. Random controls match length, and shuffles preserve mononucleotide composition; no dinucleotide preservation is claimed.

| Workload | v2 wall seconds |
| --- | ---: |
| original_comparison | 0.179 |
| reference_self_recovery | 3.074 |
| heldout | 66.336 |
| family_absent | 48.337 |
| synthetic_controls | 8.742 |

Timing uses `time.perf_counter` around the in-process search including reference validation; it excludes download, file loading and interpreter startup. One measured run per workload; no confidence interval or original-tool speedup is claimed. Hardware can differ between runs.

## Reproduce

```bash
python scripts/fetch_legacy.py
python benchmarks/run_benchmark.py --output-dir benchmarks/results
```

The raw third-party assets stay in ignored data/. Public result files contain derived scientific comparison counts and query identifiers/hashes, not downloaded reference sequences. See docs/LEGACY_PROVENANCE.md and docs/SCIENCE.md for source and interpretation details.
