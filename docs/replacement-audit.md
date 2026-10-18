# Replacement audit — 28 September 2026

Baseline: the local **WebServer-Identify-SpacerPlacer-evOr** project, not the
historical 2021 CRISPRloci server reviewed in `improvenments/01-parity-audit.md`.

CRISPRloci v2 is an additive superset of this predecessor's public application
contract and scientific workflows. This establishes application compatibility
and broader functionality, not greater scientific accuracy or faster execution.
Full production retirement and a claim that private data was never previously
pushed are not established by this local audit.

| Capability | Predecessor | Current v2 and evidence |
| --- | --- | --- |
| Detection, reconstruction, orientation | Three modes | All retained; backend, Integration and scientific-result UI regressions |
| Public API | Eight routes | All retained unchanged; seventeen declared routes; frozen predecessor-contract test plus live OpenAPI contract tests |
| Additional scientific workflows | Not in its public modes | Cas, tracrRNA, complete loci, adjacent leader context, supplied repeat folding, contextual repeat folding, target-to-spacer search, spacer-to-virus search: eleven modes total |
| Reference exploration | No reference catalogue routes | Accession, arrays, repeats, spacers, Cas systems/genes, linked records and exports; catalogue API/UI tests |
| Job privacy and lifecycle | Capabilities, queue, cancellation, recovery, retention, downloads | Retained backend hardening tests; browser submission, authenticated download and reload recovery passed |
| Scientific uncertainty | Explicit missing values and orientation uncertainty | Retained; flank context is not leader prediction, folds are model evidence, sequence matches do not establish infection |
| Scientific execution | Separate workers and environments | Retained; immutable snapshots, source/tool inventories, exact-release gates and activation locks |
| Existing data | Historical queue and artifacts | Atomic migrations through schema 6; capability, artifact and salt preservation tests; incompatible rollback requires matching data snapshot |
| Production resource and backup controls | Bounded units and authenticated SQLite backups | Resource policies retained; restored predecessor policy/backup tests; worker temporary mounts require sticky writable mode |
| Publication privacy | Frontend-only exporter and ignore rules | Added private operations journal, broader ignores, staged-blob/outgoing-history hooks, exact private identifier scanning and guarded frontend export |

All predecessor backend, Integration and frontend test files were found in v2.
Of 132 named backend test functions, 131 names remain; the old additive-digest
rollback test became `test_additive_digest_migration_preserves_legacy_insert_compatibility`.
The newer schema migration tests explicitly preserve existing jobs, capability
digests and downloads. All 35 named Integration test functions remain. Test
counts alone are not proof of parity: the frozen public contract and behavioral
checks provide the executable baseline without reading the old directory.

## Changes made by this audit

- Preserved the complete deployment record and SSH configuration privately;
  wrote a dated connection procedure/journal, record index and integrity hashes.
  Private directory/file modes are 0700/0600. Existing output/evidence folders
  are ignored and their top-level directories now have mode 0700.
- Replaced publishable deployment details and the frontend's hard-coded
  infrastructure hostname with generic operator instructions.
- Added pre-commit and pre-push guards. Tests exercise real Git commits and
  local pushes, force-added ignored files, staged/worktree differences,
  removed-in-later-commit secrets, non-commit tags, copied VM identifiers,
  symlinks and actual frontend export. Logs redact matched private identifiers.
- Restored private-monorepo CI, frozen predecessor API compatibility, service
  resource-policy checks and the unchanged backup helper's security tests.
  CI runs application/fixture checks; real scientific release approval remains
  the separate operator gate with pinned tool environments.
- Updated the README to eight tools, eleven modes and seventeen gate cases;
  corrected schema migration/startup/rollback instructions; restored an
  operations runbook and removed broken references to predecessor-only reports.

## Verification and boundaries

The full gate passed **577 Python tests** (one opt-in installed-tool test skipped)
and **192 frontend tests**. All **12 final privacy tests** passed, including two
installer cases added after the full gate. Results are recorded in
[replacement-verification.json](replacement-verification.json).
The complete `scripts/check.sh` gate checks Python, frontend, lint, formatting,
public asset scanning and the production build. The additional browser runs use
a compiled build: a disposable local API/fixture worker for recovery and
intercepted synthetic evidence for both spacer-search directions. They passed
mobile layout and JavaScript/CSP checks. They are not new scientific validation.

The standalone audit inspected 221 runtime code/configuration files and 22,553
runtime/tool symlinks without finding a dependency on the predecessor directory.
The tests and setup do not require that old source folder.

The retained deployment evidence reports seventeen passing real-science cases
for the eight-tool release. That evidence was inspected, not rerun; the selected
local immutable release was not rebuilt or activated. A future release still
needs its exact real-science gate. Earlier dated evidence in `verification.md`
describes earlier releases, not the latest workspace.

The predecessor's installed edge bootstrap/transaction/redaction controls are
not all implemented by v2's local release manager. Existing host controls must
remain until a separately authorized deployment validates their replacement.
No live service or SSH connection was used in this audit. Do not infer remote
health, successful production rollback or permission to remove installed guards.

Git metadata is inaccessible in this workspace. Hook behavior was verified in
disposable repositories, but installation in the real checkout and its prior
history audit remain required before publication. `.gitignore` cannot erase old
commits, and hooks can be bypassed. See [private operations](private-operations.md)
for installation, signature registration and the limits of these safeguards.
