# Scientific tool refinement — 13 October 2026

Each of the 15 included scientific tools received a separate code review and
refinement. CasAndra and CRISPRidentify-v2 were excluded; their source pins remain
unchanged. The work removes repeated computation, simplifies maintained code and
adds concise navigation/invariants for future agents. Scientific model bytes,
thresholds, schemas, coordinates, source identities, ordering and uncertainty
remain compatible; demonstrated engineering defects have targeted regressions.
Tool package versions are unchanged. Exact deployment pins belong in
[tools.lock.json](../tools/tools.lock.json).

Backend/frontend release sources are 2.12.2; Integration remains 0.14.0.

| Release step | Status |
| --- | --- |
| Per-tool review and local tests | Complete for all 15 tools |
| Source publication | All 15 tool commits published privately; public frontend commit `6c93629a28cf` pushed |
| Fresh pinned scientific environments and native differential controls | Passed, including rebuilt SpacePHARER, exact ISEScan/RNAalifold runtimes and full-pVOG Phigaro |
| Platform and frontend regression | 1,484 platform tests passed, zero skips; 396 frontend tests and production build passed |
| Packaged HTTP release gate | 37 native cases across 20 modes and 5 predispatch rejections passed |
| Production backend | 2.12.2 healthy and committed; 9 isolated and 9 live scientific cases passed |
| Operator source verification | All 17 tool checkout HEADs checked; CasAndra and CRISPRidentify pins unchanged |
| Existing data retention | All 40 prior jobs and all 40 prior job directories retained; no deletions |
| Pages deployment and browser verification | Exact tested commit deployed; all 47 frontend file hashes, real RNA-family evidence, authenticated downloads, handoff, recovery and desktop/mobile checks passed; no JavaScript/console/CSP errors |

The table combines completed module suites and fresh candidate-native controls.
These engineering comparisons preserve baseline scientific outputs and do not
establish new biological accuracy or calibration.

| Tool | Reviewed commit | Refinement and significant fixes | Validation and remaining limits |
| --- | --- | --- | --- |
| CRISPRtracrRNA-v3 | `9811f27d2025` | Stream native evidence parsers and bound Type V motif pairing; preserve scientific errors when command journaling also fails. | 393 passed (386 portable, 7 native); exact native parser/motif evidence. |
| SpacerPlacer | `0e4c22494829` | Reuse order-graph indegrees, iterate reachability and make summaries reusable; correct independent-gain diagnostic accounting. | 154 passed; 48 scientific columns and 28 exports match frozen controls. |
| evOr | `5120465837e3` | Vectorize orientation selection while retaining fit/tree provenance; repair diagnostics and trusted-pickle/export option propagation. | 131 passed; fitted numeric, tree and result controls match. |
| CRISPRleader-v2 | `b53101638356` | Reuse PWM/window work and serialize scientific snapshots in one traversal. | 320 passed; 120 randomized cases match. Leader prediction remains unavailable. |
| CRISPRrepeat | `602155c32968` | Index oriented references once; own ViennaRNA model and thermodynamic snapshots so other native calls cannot alter recorded folds. | 86 passed in source and dedicated runtime; complete native control report unchanged. |
| CRISPRspacer | `15298183a59b` | Stream exact matches and reuse occurrence metadata; replay verified HSP bytes and strictly bind evidence types. | 243 passed with real BLAST; 17 reports and 17 bundles byte-identical. |
| CRISPRmap-v2 | `4b35029a3600` | Reuse exact/tied alignments and linear FASTA accounting; import biological ordinals correctly and reject malformed source data atomically. | 230 passed, 1 optional fold skip covered by 10 separate pinned-native tests; 3,527 reference queries match. |
| RepeatTyper | `d55ebf535f21` | Precompute canonical 4-mer columns, batch feature counts and deduplicate sequence work; enforce input bounds before materialization. | 68 passed in source and wheel; full 291×37 score matrix and 138-column features retain exact parity. |
| SpacePHARER | `4f264252846a` | Index native HSP rows in original order; clean probe process groups and preserve worker-thread signal compatibility. | Fresh rebuilt engine: 97 passed, 1 optional original-upstream oracle skip; 2 CTest controls passed. Frozen/pre-refinement parity and sanitizer controls passed. |
| CCTK | `c2f619fb3dde` | Cache overlaps/memberships, join FASTA and slice deletions; repair singleton, large and zero-sharing ordering and impossible deletion states. | 160 passed in source and wheel; frozen exports/plots plus 650 ordering and 500 seeded RNG controls match. BLAST/MinCED databases unprovisioned locally. |
| clinker-v2 | `87dd54db5787` | Reuse ordered protein alignments and symmetric synteny work; repair stale group overlap, mutating deserialization and ranged GenBank extraction. | 710 passed; four frozen reports/CLI exports and five real examples retain parity. Large unbounded/browser paths not benchmarked. |
| ISEScan-v2 | `74d7969e3eb7` | Share HMM completion checks, sweep intervals and index ranked boundaries; reject invalid BLAST/TIR evidence and close native export handles. | 76 passed; 750 selection controls and eight native runs match. Fresh pinned Python 3.11/BLAST 2.16 release controls passed. |
| DefenseFinder-v2 | `ff7c62f12051` | Stream bounded tables and reuse memberships; bind parsed input bytes, finite HMM evidence, source proteins and global FASTA ranks. | 40 context tests and 44 retained upstream tests passed; native positive/negative/multisource parity. Two unchanged classic-CLI missing-model failures remain outside the bounded distribution. |
| RNAalifold-v2 | `d5bc9ffdc117` | Own comparative ViennaRNA parameters with family-size scaling; bound native alignment/evidence checks and reuse column/pair work. | 54 passed; MFE, partition energy and every probability match for seven families. Fresh pinned ViennaRNA 2.7.0/MAFFT 7.526 controls passed. |
| Phigaro-v2 | `b4f32844b831` | Group genes per source, join FASTA and vectorize the existing triangular reduction; reject orphan protein evidence and correct Python wheel tagging. | 71 passed in source and wheel, plus 36 integration evidence tests; recorded full-model and fresh Prodigal parity. Fresh complete pinned pVOG native comparisons passed. |

An independent source review found no actionable regression in SpacerPlacer,
evOr, CRISPRspacer, CRISPRmap-v2, RepeatTyper, CCTK, clinker-v2, ISEScan-v2 and
DefenseFinder-v2. Additional controls passed: 3,000 spacer DAGs, 1,000 orientation
selections, exhaustive 4-mers plus 1,500 random feature vectors, 2,400 interval
sets, 5,000 exact searches including work accounting, and 1,600 synteny matrices.
This review complements the module suites and does not replace release validation.

Local measurements demonstrate reductions in specific repeated work. They are
synthetic component measurements; no end-to-end production speedup is claimed.
CRISPRmap's pure-Python alignment fallback was tested without RapidFuzz; its full
modern catalogue campaign and original executable timing were not rerun. Explicit
modern catalogue rebuilds may change physical node identifiers/digests with the
corrected ordinal traversal; source occurrences and immutable existing releases
are preserved. Existing scientific scope and documented uncertainty remain in force, including unavailable
leader prediction, uncalibrated model evidence and heuristic viral associations.

Fresh scientific controls used the rebuilt candidate engine and exact pinned
native environments. They include SpacePHARER's full candidate suite and both
CTest controls, ISEScan with Python 3.11/BLAST 2.16, RNAalifold with MAFFT 7.526,
and Phigaro full-pVOG baseline/candidate comparisons. Source, installed-package,
native executable and model receipts were verified before the packaged release
gate and production activation. Operational receipts remain in private storage;
public publication exports only the frontend.
