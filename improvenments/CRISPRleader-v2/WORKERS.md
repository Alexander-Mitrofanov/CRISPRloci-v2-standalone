# Independent workers and coordination

Status: owner confirmation received; three implementation workers were launched
from a shared contract commit in separate worktrees. The table below preserves
the approved responsibilities. Engineering components were delivered; biological
reference construction, inference and evaluation remain incomplete after safety
restrictions stopped that work. See the current implementation status in the
standalone repository. There are four concurrent slots: coordinator plus three workers.

## Council responsibilities

| Council member | Bounded task |
| --- | --- |
| `leader_science_council` | Original and recent primary literature; algorithm scope; biological interpretation. |
| `leader_reference_council` | Source/model inventory; licenses/provenance; independent labels and benchmark design. |
| `leader_architecture_council` | Standalone interfaces; coordinate correctness; worker ownership; runtime and release gates. |
| Coordinator | Verified GitHub identity, source boundary cross-checks, proposed choices and reviewable plan. |

## Contract commit before parallel coding

The coordinator first defines and reviews schemas and interfaces for
`SourceRecord`, `ArrayAnnotation`, `FlankCandidate`, `RawEvidence`,
`LeaderPrediction` and `ReferenceManifest`, plus fixture IDs and coordinate
semantics. Freeze signatures for input loading, candidate extraction, scoring,
decision, independent validation and bundle writing. These are proposed names,
not implemented APIs.

Workers start from that same commit in isolated branches/worktrees, with a
separate temporary directory and no writes into another tool checkout. Share
contract changes through a coordinator-reviewed commit, not concurrent edits.
Communicate assumptions early; do not invent incompatible local schema variants.

## Implementation ownership

| Worker | Exclusive modules and outputs | Work possible independently | Handoff gate |
| --- | --- | --- | --- |
| W1 — Genomic inputs and mapping | `sequence_io.py`, `coordinates.py`, `adapters/`, corresponding unit/property tests | FASTA/native annotation parser, canonical Identify adapter, segments/strand/window extraction using frozen fixtures | Source binding, coordinate conversions, reverse complements, circular and clipped records all verified. |
| W2 — Scientific inference and references | `inference/`, `references/`, reference manifests/build recipes, scientific unit tests | Candidate-family retrieval, raw profile/motif evidence, provenance verification and abstention using fixture candidates/mock runner | Usable reference provenance, declared score semantics, tied/unsupported cases, decision-policy and calibration records. |
| W3 — Runtime and independent evaluation | `runtime/`, `benchmarks/`, process tests and independent biological/CLI fixtures | Safe process runner, benchmark accession ledger, split/leakage audit and oracle tests from frozen interfaces | Failure/timeout/cancellation tests, independent labels/controls, untouched final evaluation and reproducible report. |
| Coordinator | Shared schemas/domain types, CLI composition, result validation/serialization/export, package metadata, dependency locks, CI, docs | Compose narrow interfaces; review/reference rights tracking; packaging/export tests | Reviewed merges, complete end-to-end tool, reference/model identities and justified release claims. |

Use mocks only to develop interfaces; real scientific smoke tests are required
before integration is considered done. W3 must not derive expected biological
answers solely from W2's implementation output. W2 may inspect development
partitions but does not tune against final evaluation labels.

Each handoff includes a commit, affected contract/version, tests and results,
evidence provenance, known limitations and any unresolved question. No worker
updates dependencies, reference snapshots or another owner's modules without
coordinator review. Keep complete stage order visible in composition.

## Repository quality gates

- Typed immutable domain records, `src/` packaging, minimal declared runtime
  dependencies, a supported Python matrix and deterministic effective config.
- Lint/type checks and tests for real failure modes, not tests that restate the
  implementation. Fixtures use synthetic or explicitly permitted biological data.
- Source reconstruction independent of predictor code, malformed/stale input
  rejection, stable IDs, truthful null/zero semantics and schema round trips.
- Shell-free external programs, bounded resources/output, cancellation and
  timeouts; no network/model updates during prediction.
- Fresh wheel installation and offline CLI use from a directory outside any
  CRISPRloci/predecessor checkout; no absolute-path or editable-install reliance.
- Scientific result consistency across supported thread counts; ignore expected
  run-time/provenance timestamps when comparing deterministic scientific content.
- Explicit model/data licenses, hashes, build recipe and reproducible benchmark
  manifests. No fake production model made from synthetic fixtures.

## Boundary with existing tools and server

Identify supplies optional canonical array evidence; CasAndra supplies optional
context. Neither must be reimplemented, modified or imported by the new tool.
evOr evidence, when available later, remains a separate group-level hypothesis.
CRISPRtracrRNA and SpacerPlacer keep their current responsibilities/environments.

The new repository is independently versioned. Creating a local checkout under
`tools/` does not add it to the active scientific release. Later integration
requires an explicit contract/release plan, independent backend validation,
public export rules and relevant browser tests. Workers cannot operate
production services or silently activate a new server runtime.

## Implementation handoff record

- W1: commit `65253eb`, 106 input/mapping/adapter tests; later added 44 export/validation regression tests.
- W2: commit `288748f`, 44 reference-container and synthetic numerical tests. No biological predictor or default bundle.
- W3: commit `9a9aaf3`, 19 runtime tests; benchmark protocol `f3c3daa` and pre-output correction `183dbf9`. No biological evaluation result.
- Coordinator reviewed and merged these into the private repository, added explicit extraction-only behavior, CLI, packaging and checks. The final aggregate check record is maintained in that repository.

Worker counts above are component results, not independent biological validation.
No worker's stopped biological research was rerouted to bypass the restriction.
