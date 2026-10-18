# Scientific tool review — 11 October 2026

Ten separate tool agents reviewed the independent checkouts. CRISPRidentify and CasAndra were excluded from source review and retain their existing revisions. The initial review produced local commits without production operations. A subsequent explicit production request authorized publishing all ten commits to their private origins and staging deployment; execution evidence remains private.

## Reviewed revisions

| Tool | Version | Commit | Main fixes | Source tests |
| --- | --- | --- | --- | --- |
| CRISPRtracrRNA-v3 | 3.0.0 | `ca225a2a3be3` | Duplicate JSON, invalid numeric support, repeated cassette identity, Infernal strand/coordinate conflicts | 378 passed |
| SpacerPlacer | 1.0.1.post3 | `87993c9f928f` | Group and alignment identity, supplied trees, fitted-tree exports, finite CLI options | 138 passed |
| evOr | 1.0.3 | `369241d16a72` | Group/cache identity, supplied/selected tree integrity, finite options, event node identities | 124 passed |
| CRISPRleader-v2 | 2.0.0a1 | `5051d07ab827` | Context settings and source provenance, bounded JSON reads, atomic annotation publication | 316 passed |
| CRISPRrepeat | 1.0.1 | `eb26716becfc` | Fixed folding parameters, native numeric/callback checks, substrate budgets, bounded reference matching | 67 passed |
| CRISPRspacer | 0.3.1 | `0e4951664862` | Typed replay, corrupt gzip, source/topology binding, SQLite integrity, cancellation/deadlines | 225 passed |
| CRISPRmap-v2 | 2.1.0a2 | `b1990e3fa08b` | Strict annotation/JSON imports and portable read-only catalogue partition merging | 199 passed |
| RepeatTyper | 2.0.1 | `29727e2456d3` | Atomic publication, anchored model provenance, checked resource snapshots and prediction matrices | 59 passed |
| SpacePHARER | 2.1.1 | `f11ee0d85c70` | Duplicate/gapped evidence, calibration populations, input snapshots, bounded catalogue/FDR lookup and packaging | 51 passed |
| CCTK | 2.0.2 | `ae7e97b21042` | Repeated-spacer plots, safe publication, input/traceback checks and ambiguity-code reverse complements | 138 passed |

CRISPRmap skipped one optional ViennaRNA test. SpacePHARER skipped one separately built original-upstream oracle; three exact prepatch native calibration controls and frozen biological evidence checks passed. CCTK passed its 138-test suite against both source and the installed wheel. RepeatTyper preserved exact upstream CPU scores on frozen, random and reverse-complement controls. tracrRNA passed seven real native integrations; spacer searches passed 109 focused controls with pinned BLAST 2.17 and installed-wheel round trips.

## Server integration

Backend/frontend versions are 2.11.1; Integration is 0.13.1. The lock and frontend tool documentation revisions now match all twelve clean checkouts. Existing JSON/report schemas, scientific policies and defaults remain compatible. Validators explicitly support the reviewed patch versions alongside supported older outputs. SpacePHARER also accepts existing 2.1.0 reference receipts; the database format did not change.

Integration rejects duplicate/nonfinite annotation JSON and bounds object/JSONL/checksum reads. Its repeat gate now independently verifies context coverage, orientation, offsets, consensus folds, finite energies and pair/accessibility arithmetic, plus unavailable comparison evidence. Both repeat validators reject booleans masquerading as coordinates/counts. Both association validators reject double-gap columns. Startup independently accepts only the explicitly supported spacer and association versions when references are configured. Backend validation continues to import neither Integration validators nor scientific packages.

Provisioning verifies completed toolsets rather than returning unverified cached installations, and validates source/launcher/environment receipts before resuming. The scientific release gate includes reconstruction and orientation jobs with branch extension, checks API source revisions against the scientific manifest, and binds download provenance to the same exact identities.

## Verification status

`scripts/check.sh` passed the complete source gate; after the final startup compatibility correction, `scripts/check.sh backend` passed **1,130 Python tests with no skipped tests**, using the new immutable tool executables. Frontend checks passed **298 tests**, lint, formatting, source boundaries, production build and public export scan.

The fresh twelve-tool toolset is `1d22e92ee2c1d41437f58830e255bd5897ac8f83c8d850bb4888fb47444489c0`. The final backend/Integration release is `cffc2de6e3c2682c0990308eae7bb11cc4df67303456181c14897f0b404e1dee`. Source, installed package/native/model, launcher and dependency receipts verified. Its packaged HTTP gate passed **26 real scientific jobs across all 15 modes**, including positive/negative controls, both RNA model families, both branch-extension cases, reference-dependent workflows, exact source identities and authenticated complete archives. This gated release is selected at `.runtime/current`; `.runtime/previous` retains the prior release. No service was started in production.

Headless Chromium against that exact packaged release passed real locus submission, validated evidence table counts, repository provenance, authenticated ZIP download, reload recovery, mobile layout and no JavaScript exceptions. Evidence is local-only under `output/tool-audit-20261011/browser-real-release/`. It displayed five Cas rows and four tracrRNA rows.

The packaged gate exposed and repaired two startup checks that still rejected the reviewed spacer/association versions. Eight configuration regressions cover supported old/new patches and unsupported versions. The browser's old 97-candidate assertion also described an obsolete Identify reconstruction: the existing e8a3d7c pin masks the complete 5481–7366 array instead of the old 6142–6442 fragment. Original and reviewed tracrRNA parsers produce the same four candidates from the new retained evidence; the UI check now binds to validated API counts.

One repeated backend run hit the four-gigabyte worker admission reserve during runtime construction. Stored heartbeats confirmed `storage_low`; the guards were unchanged. Review-created test environments/builds were removed, and the final gate passed using the immutable executables with adequate reserve.

Earlier local checks passed 235 Integration/release/option tests, 101 focused independent backend tests including real BLAST (one opt-in lifecycle skipped at that stage), 21 real RepeatTyper server tests, 20 CCTK server tests and 36 SpacePHARER server tests. The association and comparison checks include actual HTTP queue → isolated worker → authenticated export flows. The final repository gate enables those real executables and the previously opt-in spacer lifecycle.

Headless Chromium sequence-exchange acceptance passed upload/rejection, authenticated compact/sequence downloads, direct continuation, recovery preservation and mobile layout against an intercepted fixture API. The standalone audit found no executable/configuration references or runtime symlinks into the two retired source trees.

## Availability and scope

Exact source commits exist in the local independent repositories and offline lock verification passes. The subsequent production request published all ten reviewed commits to their private GitHub origins, making the exact pins available to fresh restores. The parent source workspace has no accessible Git metadata; install its .githooks in a real checkout before publication. Privacy checks passed in the workspace and every affected tool commit.

The web workflows retain their existing scientific scope: leader extraction supplies adjacent windows, modern CRISPRmap catalogue commands remain operator-only standalone tools, viral retrieval remains heuristic, and model scores retain their documented uncertainty. This engineering review does not establish new biological calibration or large-reference performance claims. Legacy trusted-pickle/simulation paths were not comprehensively reaudited.

## Subsequent production rollout

An explicit production request deployed backend/frontend 2.11.1 and Integration 0.13.1 on October 12. All ten reviewed commits are published in their private origins. The host passed 26 native release cases, four isolated service checks and four live checks; all 17 prior jobs, configuration and references were preserved. The public frontend publication contains only the three reviewed tool-reference/version files. Live browser asset identity, native typing→folding, authenticated handoff, recovery and mobile checks passed. Operational identifiers, receipts and rollback records remain exclusively in ignored private operations storage.
