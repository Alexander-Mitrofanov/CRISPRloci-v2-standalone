# Scientific tool modernization and integration

Requested 10 October 2026: modernize and integrate RepeatTyper, SpacePHARER and
CCTK sequentially, each with a dedicated agent and a separate repository under
Alexander-Mitrofanov. Create each repository before refactoring; publish the
validated standalone refactor before integrating its pinned commit into the
server. Repositories default to private, matching existing scientific tools.

The scientific algorithms, original attribution and applicable licenses must be
preserved. Modernization covers modular first-party code, typed interfaces,
validated I/O, consistent failures, reproducible dependencies/model identities,
offline execution, tests, CI, packaging and documentation. Third-party code must
remain explicitly attributed and versioned. Scientific behavior changes require
documented evidence rather than being hidden inside a refactor.

| Stage | Repository | State |
| --- | --- | --- |
| 1. RepeatTyper | https://github.com/Alexander-Mitrofanov/RepeatTyper | Standalone v2.0.0 published; server integration accepted |
| 2. SpacePHARER | https://github.com/Alexander-Mitrofanov/SpacePHARER | Standalone v2.0.2 published; server integration accepted |
| 3. CCTK | https://github.com/Alexander-Mitrofanov/CCTK | Standalone v2.0.1 published; server integration accepted |

Each stage requires standalone scientific parity checks and a pushed GitHub
revision, then isolated worker integration, independent backend validation,
source-bound machine-readable results, frontend workflow, immutable runtime
provisioning and release fixtures. Run repository checks and local browser flows.
Record measured outcomes and remaining limitations below as stages finish.

Production deployment and activation of production services are outside this
request. Candidate builds and disposable local verification remain in scope.
Never publish backend configuration, user jobs, output or private operations
records to the new tool repositories. Install Git publication guards and run
privacy scanning before each push.

Initial verification: the sandboxed baseline stopped progressing in the known
FastAPI threaded-test environment and was interrupted. The Python check passed
in the approved local execution environment: 707 passed, 1 skipped. Frontend:
242 tests passed, source checks and production build passed. Ruff and private-data
checks also passed. These are the pre-integration baseline, not verification of
the new modes; no test assertions were weakened.

RepeatTyper standalone checkpoint: commit
`7a153a08ae31682c99b5a255d01e5408ae267d66`, tag `v2.0.0`. Agent reported
22 installed-wheel tests passing, including exact label and native-score parity
with the original classifier on frozen upstream examples, 100 seeded random
valid repeats and reverse complements. Ruff, mypy, source/wheel builds, CLI
publication and privacy checks passed before server integration began.
Parent review independently ran the installed CLI on the five-record example:
three predicted, one uncertain and one unsupported result; duplicate occurrence
IDs and explicit nulls were preserved. The upstream commit matches the GitHub
`v1.9.0` tag. GitHub CI passed on Python 3.11, 3.12 and 3.13 for the published
commit: https://github.com/Alexander-Mitrofanov/RepeatTyper/actions/runs/38056465118.

RepeatTyper integration accepted: `repeat_type` mode; backend/frontend 2.8.0,
Integration 0.10.0, scientific manifest 7.0.0. Full `scripts/check.sh` passed:
735 Python tests, one existing skip; 248 frontend tests, lint, source boundaries,
production build and privacy check. A subsequent focused runtime/manifest check
passed 34 tests. Real HTTP submission, isolated native classification, independent
validation and authenticated downloads passed, including a detection repeat-port
handoff and retained occurrence origins. Both server wheels were built, inspected
and executed in a native smoke test. Parent reviewed desktop/mobile browser
screenshots and receipt. See [repeat typing documentation](../docs/repeat-typing.md).
The ten-tool immutable runtime was not rebuilt or activated; provisioning and
candidate gate support were implemented and tested separately.

SpacePHARER source baseline: upstream commit
`d1e94c04025be6e04fd8204f4f8bb1ea782e993c`. Original native engine built locally
for comparison before standalone modernization. Audit identified an unsafe
FDR slope-vector access for sparse controls and an empty-control FDR-zero claim;
the agent is treating these as explicit correctness changes, separately from
parity-preserving refactoring. The order-sensitive combined-score fallback is
being retained and tested as upstream behavior. Validation is still in progress.

SpacePHARER standalone published before server edits: commit
`6cae31af73b0aee0542d968a0490a03b17935e17`, tag `v2.0.0`. Agent reports
21 Python tests and native Release CTest passed, including byte-exact calibrated
FDR parity on populated target/control score distributions and all six biological
fixture alignments. Ruff, mypy, sdist/wheel builds, privacy and clean-PATH
installed-wheel execution passed. Parent independently ran the complete native
CLI pipeline: combined scores 12.24 and 74.56, six source-bound alignments,
explicitly unavailable calibration and null FDRs, with input groups retained.
Calibration is read from a versioned native receipt. Boundary rows above the
requested FDR remain candidates. Reference preparation currently supports targeted
panels of at most 100 single-contig genomes and 50 million total bases.

Publication gate caught a source packaging defect in the initial SpacePHARER tag:
GitHub clean-checkout CMake failed because required vendored build files were
omitted by ignore rules. The agent identified 135 upstream regular source/data
files affected by imported ignore rules. These must be included and archive
completeness tested before the server pin is accepted; v2.0.0 is not the accepted
integration revision. No scientific server edits had begun when this was caught.

Corrected standalone checkpoint: commit
`82ac8222aea4f81b64757e5bcb6ef6267494a53c`, tag `v2.0.1`. A clean Git archive
configured and built successfully; native Release CTest and 22 Python tests
passed from the archive, including an upstream-file completeness regression.
An installed-wheel execution with a clean PATH also passed. This supersedes
v2.0.0 for integration; its original tag was retained. The server integration
began only after this corrected standalone checkpoint.
The corrected tag also passed GitHub's clean-checkout CI, including the native
build and scientific parity checks:
https://github.com/Alexander-Mitrofanov/SpacePHARER/actions/runs/38058589673.

Integration review prompted standalone v2.0.2, commit
`3121a8eb3fd5b79996c7175552a19aa465486458`: reference loading now checks the
complete database inventory, aggregate size limits and symlinked index paths.
The native scientific calculations are unchanged. This new tag is the server's
integration pin; earlier tags remain intact.
Standalone v2.0.2 passed 25 tests and active native Release CTest. GitHub CI also
passed for this exact pin:
https://github.com/Alexander-Mitrofanov/SpacePHARER/actions/runs/38059479794.

Parent packaged-server acceptance passed for backend 2.9.0 and Integration
0.11.0. Both wheels were built and extracted into a separate directory; imports
were verified to come from the wheels. The real v2.0.2 engine processed 71
occurrences in three explicit groups: two uncalibrated candidates, six support
alignments and one unsupported occurrence. Reordered inputs and duplicate
occurrence identities survived the full compact export, while preview limits
remained separate. The independently loaded backend rejected a deliberately
altered alignment. Local receipt: `/tmp/spacepharer-wheel-smoke.json`.
After review corrected the ZIP exclusion policy and reference-manifest readiness
check, both final wheels were rebuilt and the smoke repeated successfully. All
116 packaged Python/JSON files match the current source, and the packaged
Integration CLI reports 0.11.0.

SpacePHARER integration accepted: `spacer_association` mode, backend/frontend
2.9.0, Integration 0.11.0, scientific manifest 8.0.0 and eleven pinned tools.
The full Python check passed 775 tests with one existing skip; Ruff passed.
After correcting a missing mode-schema entry, the focused contract/runtime check
passed 32 tests. Final frontend check passed all 254 tests, formatting, lint,
source boundaries, production build and public scan. The real HTTP worker tests
passed positive, negative, unsupported and detection-port inputs, per-source
grouping, duplicate projection and protected exports. Parent reviewed final
desktop/mobile browser evidence and refreshed grouping controls; browser checks
passed without page errors, CSP violations or horizontal overflow. Final browser
receipt: `/tmp/spacepharer-browser-final/acceptance.json`. The immutable eleven-tool
runtime was not rebuilt or activated; see [integration documentation](../docs/spacer-association.md).

CCTK began only after SpacePHARER acceptance. The sub-agent service refused a new
thread at its thread limit, so an existing council agent was assigned exclusively
to CCTK implementation. The private repository was created before source
refactoring. Original algorithm/oracle review is in progress; scientific server
integration must wait for a tested, published standalone checkpoint.

CCTK standalone v2.0.0 was published privately at
`54a7d6638dfb10d5ac490e3df6ab68df10682d9f` after a clean Git archive, sdist,
wheel and 96 installed-package tests. The original first-party modules were
modernized in place within a new installable package, retaining the original
alignment, parsimony and simulation algorithms. Original GPL attribution and
full upstream source identity are recorded. A reproducible oracle generator
verifies upstream commit `fe94fde10a2855cbcc3bf118a253ec097d3a5e46` before
regeneration. Parent independently regenerated and compared eight fixture files,
exact CRISPRdiff image pixels and the expected tree/constrained-tree outputs.
Observed-array plots separately preserve deletion slots as visible gaps and
retain arrays without observed spacers as explicit isolates.

Parent process review found that the optional BLAST thread-pool path used
`preexec_fn` for output limits, which is unsafe in threaded callers. A corrected
standalone patch release is required before integration; v2.0.0 is not the final
server pin. The correction must retain bounded output and descendant cleanup.
No CCTK scientific server edits had begun at this checkpoint.

CCTK corrected standalone checkpoint: tag `v2.0.1`, commit
`6fb81fec60f4951bf23259199ea3d9bc8f8a3d22`, published before server integration.
A fresh interpreter exec helper safely applies process output limits when called
from threads, and invalid external UTF-8 produces an explicit tool error.
Clean Git archive → sdist → installed wheel verification passed 99 tests,
Ruff and type checking of the five public boundary modules. Tests cover concurrent
launches, deadlines, SIGTERM, descendants and output limits. Native and observed
plots cap repeated connector rendering at 10,000 with an explicit omission note;
complete JSON and TSV comparisons remain available. Privacy checks passed.
Server integration began only after this corrected publication checkpoint.
GitHub clean-checkout CI passed for this exact CCTK pin:
https://github.com/Alexander-Mitrofanov/CCTK/actions/runs/38062075753.

CCTK runtime wiring implemented: twelve-tool lock, separate locked Python
environment, `cctk` runner and scientific manifest 9.0.0. Focused runtime/release
verification passed 64 tests plus Ruff. Tests export the exact published Git
revision, include the installed exec helper, preserve old optional toolsets and
reject modified installed code or incomplete release identities. The candidate
HTTP gate now includes an existing detection array export and a new synthetic
sharing fixture. Parent ran the actual installed CCTK on that fixture: exact
matching produces one edge, three components and two isolates; explicit reverse
complement matching produces two edges, two components and one isolate. Duplicate
DNA occurrences, scoped array IDs, deletion-only arrays and unknown/reverse strand
remain intact. These native checks do not constitute a consolidated runtime gate;
the complete API/worker/frontend acceptance remains in progress.

Parent CCTK packaged-server check passed: backend 2.10.0 and Integration 0.12.0
wheels were built and loaded from isolated extracted packages. All 122 packaged
Python/JSON files match source; the installed Integration CLI reports 0.12.0.
Real CCTK 2.0.1 execution passed exact, explicit reverse-complement and
minimum-shared-two controls. A one-array preview still yielded complete
four-array JSON and TSV exports, including duplicate occurrences, gaps and
unchanged reusable input arrays. Independent validation rejected altered symbols
for every case. Local receipt: `/tmp/cctk-server-wheel-smoke.json`.
Admission checks also measure the actual in-flight input tree and confirm that
array maps fit the inode reservation and existing operator record limits.

CCTK integration accepted: `array_compare`, backend/frontend 2.10.0,
Integration 0.12.0 and scientific manifest 9.0.0. The full backend check passed
817 tests with one existing skip, plus Ruff. A final annotation-bound validation
fix passed 22 focused CCTK/admission tests, including actual HTTP detection
handoff, isolated worker execution and authenticated exports. The backend wheel
was rebuilt after that fix; the complete packaged native smoke passed again.
All 122 packaged Python/JSON files still match source.

Final frontend checks passed 262 tests, formatting, lint, source boundaries,
production build and public scan. Parent reviewed desktop and 390-pixel mobile
screenshots after the mobile table adjustment; all pairwise metrics are visible
without horizontal scrolling. Browser acceptance used real native public
evidence with an intercepted API, separately from the real HTTP worker tests.
It passed graph interaction, recovery, header-authenticated download, explicit
truncation and missing-evidence checks with no page errors, CSP violations or
page overflow. Receipt: `/tmp/cctk-browser/acceptance.json`. Privacy scanning
reported zero findings. All three scientific tool checkouts are clean and match
their published integration pins.

The requested repository-first, standalone-first implementation sequence is
complete for all three tools. Server changes are local; scientific repositories
are published privately. No production deployment occurred. The consolidated
twelve-tool immutable runtime still requires an operator build and candidate
gate on a host with adequate storage. SpacePHARER additionally requires its
operator-provisioned immutable reference panel before the mode is advertised.

Local verification constraint: approximately 8 GB free while the existing
immutable runtimes occupy 14 GB. Use isolated tool environments for real local
worker checks; do not duplicate all existing environments or remove user data.
Record separately whether a consolidated immutable release was built.
