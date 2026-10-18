# Natural-locus failure analysis

Stage: `final_fixed/complete`. All 10 misses ran successfully. None is solely a ranking failure: no returned candidate meets the fixed localization criterion. Every locus has an accepted array; all truth intervals lie outside array masks.

Diagnostic searches keep both repeat orientations. Isolated-RNA searches are truth-assisted explanations, not additional benchmark detections. Native genome BLAST uses the same masked input, 70% identity and 60% query-span coverage gates.

| Locus | Repeat evidence | Native FASTA alignments at truth | Isolated FASTA span; identity | CM at truth | Interpretation |
|---|---|---:|---|---|---|
| Esp2 | Exact published repeat (either strand) | 0 (passing: 0) | 10/36; 100.0% | No localizing hit reported | short repeat similarity below coverage gate |
| Bbo | Inferred repeat differs from published repeat | 0 (passing: 0) | 19/38; 84.2% | No localizing hit reported | short repeat similarity below coverage gate |
| Cca2 | Exact published repeat (either strand) | 2 (passing: 0) | 16/36; 93.8% | No localizing hit reported | short repeat similarity below coverage gate |
| Cca1 | Inferred repeat differs from published repeat | 2 (passing: 0) | 15/47; 100.0% | E=0.0045, fails 1e-5 gate | short repeat similarity below coverage gate |
| Cga | Exact published repeat (either strand) | 2 (passing: 0) | 26/47; 96.2% | No localizing hit reported | short repeat similarity below coverage gate |
| Cme1 | Exact published repeat (either strand) | 2 (passing: 0) | 16/46; 100.0% | E=0.012, fails 1e-5 gate | short repeat similarity below coverage gate |
| Ffr | Inferred repeat differs from published repeat | 2 (passing: 0) | 23/46; 100.0% | No localizing hit reported | short repeat similarity below coverage gate |
| Ghy2 | Exact published repeat (either strand) | 2 (passing: 0) | 20/46; 94.7% | No localizing hit reported | short repeat similarity below coverage gate |
| Orh | Exact published repeat (either strand) | 2 (passing: 0) | 23/47; 95.7% | E=0.0038, fails 1e-5 gate | short repeat similarity below coverage gate |
| Psp | Inferred repeat differs from published repeat | 2 (passing: 0) | 24/49; 100.0% | E=0.00077, fails 1e-5 gate | short repeat similarity below coverage gate |

The isolated FASTA column reports the highest-scoring identity-passing alignment (or highest-scoring available alignment), not the maximum coverage across every alignment. Full per-alignment evidence, native BLAST comparisons, exact coordinates, configuration and source hashes are in `failure_diagnosis.json`.

Eight loci have native FASTA alignments at the known RNA, rejected for insufficient repeat coverage. Esp2 and Bbo have no native FASTA genomic alignment at the RNA; BLAST and isolated-RNA diagnostics expose only subthreshold matches. Bbo additionally has a 26/38 nt isolated match at 69.23% identity, failing the 70% identity gate. The accepted consensus differs from the published repeat in Bbo, Cca1, Ffr and Psp, but substituting the published repeat (both orientations) rescues none at default gates. Array inference differences therefore do not explain these failures by themselves. No correct-strand CM matches are confused with wrong-strand hits: Cca1, Cme1, Orh and Psp have correct-strand hits that fail the rescue threshold; the other six have no localizing hit in the raw default cmscan output.

## Partial boundaries from model rescue

| Locus | Model-supported interval | Covered reference | Rank | E-value |
|---|---|---|---:|---:|
| Lan | 1033–1090 (-) | 58/99 (58.6%) | 1 | 6.1e-26 |
| Rsp | 14845–14904 (-) | 60/97 (61.9%) | 3 | 1.5e-06 |
| Cme3 | 5035–5092 (+) | 58/93 (62.4%) | 1 | 5.5e-10 |
| Ghe | 8132–8192 (+) | 61/96 (63.5%) | 1 | 4.6e-10 |

These four loci pass localization because a conserved subregion is detected, but do not pass the separate 90% reference-coverage metric. Their `boundary_kind=model_supported_region` is deliberate: extending to the complete transcript without anti-repeat or boundary evidence would invent precision.

## Model-threshold diagnostic

The following is a rescore of raw genomic CM hits, not a rerun or a recommended cutoff. The union with anti-repeat predictions is potential localization before rescue suppression/deduplication, not measured final pipeline recall. Counts before array exclusion/overlap deduplication expose the increasing evidence burden; unmatched hits are not established false positives.

| Maximum E-value | Localized references by CM / 78 | Anti-repeat or raw CM union / 78 | Raw CM hits / 78 loci |
|---:|---:|---:|---:|
| 1e-05 | 34 | 68 | 69 |
| 0.0001 | 35 | 68 | 78 |
| 0.001 | 38 | 69 | 89 |
| 0.01 | 42 | 71 | 110 |
| 0.1 | 45 | 72 | 153 |
| 1 | 45 | 72 | 206 |
| 10 | 45 | 72 | 307 |

Canonical Type II/unknown sequence controls combine the retained batches by sequence hash, retaining the first batch for duplicates. Membership is exactly II, II-* or unknown; Type III is excluded. They provide a separate descriptive check. E-values depend on search-space size, so these values cannot calibrate the per-genome rescue threshold. Shuffles retain mononucleotide composition and are not validated biological negatives.

| Maximum E-value | Correct-strand recognized / unique inputs | Shuffles with any reported hit / same denominator |
|---:|---:|---:|
| 1e-05 | 7063/42320 | 1/42320 |
| 0.0001 | 9680/42320 | 8/42320 |
| 0.001 | 12393/42320 | 32/42320 |
| 0.01 | 14409/42320 | 132/42320 |
| 0.1 | 15048/42320 | 382/42320 |
| 1 | 15116/42320 | 527/42320 |
| 10 | 15124/42320 | 532/42320 |

## Follow-up decisions

Do not relabel these misses as ranking failures or silently lower gates. If an alternate backend reports an otherwise eligible genomic alignment, evaluate it on the complete reference set and suitable negatives before changing the default. Shorter anti-repeat coverage and weaker CM support call for an explicitly calibrated sensitivity profile and independently held-out families. Missing native transcript boundaries require separate evidence. This positive-only, training-overlapping collection cannot establish specificity or select a new cutoff.

Reproduce with `PYTHONPATH=src python scripts/diagnose_literature_failures.py`; add `--rerun-alignments` to rerun the cached native diagnostics.
