# Modern repeat reference benchmark

This measures retrieval coverage, not biological family accuracy.

Unique heldout Bona-fide consensus queries: 1000 sampled from 58477.
Training core nodes: 174285; 31878 overlap the full heldout consensus/variant evidence set.
Excluded contigs: 0 from assemblies without complete/valid outcomes, 1037029 with cross-fold hashes, and 0 without hashes.

| Evaluation | Queries | Covered | Exact/RC | Coverage |
|---|---:|---:|---:|---:|
| Genome-heldout coverage | 1000 | 761 | 469 | 76.10% |
| Full-cohort exact/RC-purged coverage | 1000 | 661 | 0 | 66.10% |
| Sample-scoped homology90 purge stress | 1000 | 0 | 0 | 0.00% |
| Frozen historical reference (overlap unknown) | 1000 | 185 | 12 | 18.50% |

Genome holdout keeps accession versions together and admits only complete, source-validated assemblies. References contain only eligible Bona-fide training consensuses; query evidence comes from retained heldout contigs.

Exact/RC purge uses every retained heldout consensus and observed repeat from both source tiers, including unsampled arrays. Homology90 stress additionally purges every training node meeting the same 90% similarity cutoff against sampled queries and their associated observed repeat variants. Zero coverage is expected by construction; this is a sample-scoped policy stress test, not a full-cohort homology-independent split.

Taxonomy is audited but not used for splitting: this is not a species-heldout evaluation. Upstream detector training membership is unknown. Historical-reference overlap is unknown, so coverage differences do not establish predictive superiority.

Independent dynamic-programming oracle: passed, 7068 pairs in a bounded reference subset. This does not replace full index-inventory validation.
