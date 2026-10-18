# Modern reference build and benchmark — 1 October 2026

The full frozen ClusterCRISPRCasDB repeat catalogue was built and checked on eight campaign VMs. All scientific extraction, indexing, clustering, folding, tests and benchmark computation ran on VMs. The local workstation handled source edits, transport, repository operations and review of aggregate reports.

This is a computational reference release, not an experimentally validated family classification. The NCBI GenBank bacterial/archaeal source catalogue was frozen on 2 September 2026. CRISPRidentify 2.3.2 predictions retain their original evidence tier. Source validation checks software output and coordinates; it is not biological confirmation.

## Reconciled catalogue

| Inventory | Verified count |
|---|---:|
| Assembly outcomes | 3,331,847 |
| Complete assemblies | 3,331,684 |
| Failed assemblies | 151 |
| Unavailable assemblies | 12 |
| Contigs supporting arrays | 7,559,090 |
| Arrays | 9,583,665 |
| Bona-fide arrays | 3,370,380 |
| Possible arrays | 6,213,285 |
| Observed repeat occurrences | 55,023,435 |

| Sequence groups, after exact/RC deduplication | Count |
|---|---:|
| All canonical groups | 2,184,322 |
| Groups with a consensus occurrence | 1,340,023 |
| Observed-variant-only groups | 844,299 |
| Eligible unambiguous search groups | 2,175,735 |
| Eligible Bona-fide consensus core groups | 206,946 |

The group categories overlap where stated; consensus and observed-variant-only groups partition the total. The core is a subset of eligible search groups. Consensus and individual observed repeats have separate occurrence counters. Possible-only and variant-only evidence cannot connect core clusters.

Every source count reconciles. Missing consensus sequences: 0; missing observed repeat sequences: 0; arrays attached to noncomplete or unvalidated assemblies: 0. The 43,738 independently pinned source receipts include expansion of four raw source groups. The complete SQLite catalogue occupies 8,251,097,088 bytes.

## Search index and operational groups

The lossless segment index contains 47,885,443 postings for 2,175,735 eligible nodes. Deep validation checked every expected posting and core flag against the catalogue. Fixed-distance coverage is 0–6 edits; relative profiles support the declared 85%, 90% and 95% policies.

Groups are connected components under full-length unit edit distance, considering both orientations. They use a new MAP2 vocabulary. They do not inherit historical F-family or M-motif labels. Chaining can place distant members in one component; the table reports that effect.

| Minimum similarity | Components | Singletons | Largest component | Members outside representative threshold |
|---|---:|---:|---:|---:|
| 85% | 68,750 | 60,805 | 55,338 | 122,373 |
| 90% | 88,549 | 77,189 | 21,681 | 96,838 |
| 95% | 126,433 | 111,922 | 15,932 | 61,530 |

The representative is the lexicographically smallest canonical member. All query lanes were assembled and retained forest edge distances recomputed during assembly. Forests preserve connected components, but do not provide full-graph degrees or diameter.

## Separate structure evidence

ViennaRNA 2.7.0 computed 413,892 predictions: both orientations of all 206,946 core nodes. Every node/orientation, partition hash, dot-bracket length and balance, and finite thermodynamic value was checked. The fixed model uses Turner 2004 parameters, 37 °C and dangles=2. MFE and ensemble evidence remain separate from sequence groups; no transcription strand or structural motif is assigned.

## Modern retrieval benchmark

The deterministic accession-stem split samples 1,000 unique heldout Bona-fide consensus sequences from 58,477 eligible query groups. A query is covered when at least one retained training core consensus has similarity `1 - global_edit_distance / max(query_length, reference_length) >= 0.90`, considering both orientations. The threshold was fixed before evaluation.

Training contains 174,285 eligible Bona-fide consensus groups. This benchmark evaluates that core-consensus subset; it does not measure recovery across all 2,175,735 eligible search groups, including Possible and observed-only groups. Sampling weights unique sequence groups equally, not their genomic abundance.

| Reference or overlap policy | Queries | Covered | Exact/RC | Coverage |
|---|---:|---:|---:|---:|
| Full-cohort exact/RC-purged coverage | 1,000 | 661 | 0 | 66.10% |
| Genome-heldout coverage | 1,000 | 761 | 469 | 76.10% |
| Sample-scoped homology90 purge stress | 1,000 | 0 | 0 | 0.00% |
| Frozen historical 2013 reference | 1,000 | 185 | 12 | 18.50% |

Identical contig checksums occurring across folds exclude 1,037,029 contigs; 0 lack a usable checksum. The exact/RC purge uses every retained heldout consensus and observed variant from both evidence tiers, including unsampled arrays, and removes 31,878 overlapping training groups.

The sample-scoped homology stress uses sampled queries plus all eligible observed variants associated with their retained heldout arrays: 1,289 anchors; 3 anchor-evidence nodes are ineligible. Purging every training neighbor at the same 90% retrieval threshold forces zero coverage by construction. That row checks enforcement of the policy and does not establish family generalization.

The historical-reference coverage row uses the new mapper against the frozen 2013 reference; it is not an original-executable run.

Accession stems keep versions together. Different GCA/GCF accessions, related strains, taxa and homologous sequences can occur across folds after identical-contig screening. Taxonomy is audited but not used to define a species-disjoint split. Upstream detector training overlap and historical-reference overlap are unknown. Coverage differences therefore do not establish improved biological prediction. Modern family and motif accuracy require independent labels.

An independent dynamic-programming oracle checked 7,068 query/reference pairs in a bounded training subset, including 81 qualifying pairs: zero mismatches. This complements the full index-inventory check; it is not an exhaustive all-pairs alignment of the whole catalogue.

## Original CRISPRmap comparison

The fresh VM run uses authentic original 2013 example outputs and the pinned historical reference. All 7,054 reference/RC queries recover exact matches. For known original example labels, v2 assigns 10/11 sequence-family labels and 10/12 motif labels; all assigned labels agree. It abstains for conflicting Ecoli_K12 reference evidence and does not infer the nonexact Stepidermidis motif. The original executable was not recovered, so its runtime and any original-tool speedup are unavailable. See the separate historical report for heldout annotation agreement and controls.

## Execution and reproducibility

All 168 automated tests and lint passed in the isolated VM environment. CI also passed on Python 3.10, 3.12 and 3.13 and with pinned scientific extras. The final full-catalogue CLI smoke check covers core membership, reverse-complement invariance, noncore evidence, a one-edit control and ambiguous-input abstention.

| Measured stage | Wall seconds |
|---|---:|
| Merge plus catalogue integrity/publication | 483.698 |
| Index build plus deep inventory validation | 2525.703 |
| Final grouping resume and assembly; prior invocations excluded | 639.062 |
| Historical benchmark process | 129.069 |

Modern benchmark retrieval: 25.821 seconds; full benchmark before publication: 154.782 seconds. These are single-run operational measurements on a 32-vCPU analysis VM. Stages can share resources; timings have no confidence intervals and do not compare equivalent original-tool runtimes. Extraction used four workers on each of eight VMs; core graph construction used 16 query lanes and folding used four lanes.

Forest construction and the modern benchmark use commit `87fb812`; folding completed with `87041e9`; indexing used `191d18a` and extraction used `28d6fd8`. The unchanged historical benchmark ran from an earlier source snapshot, whose Python source hashes are recorded separately. Catalogue merging used the tested intermediate source archive recorded in the validation manifest; its merge implementation matches the final implementation. An independent whole-catalogue check confirms the later assembly-eligibility rule holds for every imported array. Code archive hashes, catalogue/index hashes, group manifests, structure coverage and raw aggregate benchmark metrics accompany this report.

Analysis was cleanly interrupted to apply three measured optimizations: batched disk syncs, candidate-first database lookups and a tighter necessary bound on segment offsets. Paired full-catalogue queries preserved all verified neighbor IDs and distances; repeating the complete modern benchmark preserved every scientific field. Owned checkpoints were verified before resuming; the benchmark restarted from its full-cohort audit. The coordinator timing above covers the final grouping resume and assembly only; it excludes earlier attempts and reuses already-completed folding and benchmark results. The changes preserve scientific policies, verified qualifying neighbor sets and sequence-group definitions; the raw candidate superset can shrink. A hard termination can require operator recovery of SQLite journals before resuming.

The full scientific SQLite files remain on the campaign VMs. This repository contains code, reproducibility commands, manifests and aggregate reports. It does not redistribute raw historical reference assets or infrastructure access details. Follow `docs/BUILD_MODERN.md` using the frozen source campaign, and keep source attribution and upstream data terms.
