# Pipeline and reporting council consensus plan

Council: `pipeline_architect` (orchestration and module architecture) and
`science_contracts` (independent data-integrity and failure-contract review).
Both members independently reviewed `pipeline.py` and its tests, exchanged
findings, reproduced the defects below, and explicitly agreed to this bounded
plan on 2026-09-27. A separate implementation agent will execute the plan.

## Objective and observed defects

Separate result publication from biological orchestration, make retained records
consistent with the final report, and preserve the original error when publication
fails. Preserve valid scientific predictions, their ordering and IDs, the CLI,
report schemas, export columns and bytes, and the completion-marker convention.

The council independently reproduced these problems using temporary outputs and
mocked scientific computation:

- A failed `configuration.json` write occurs outside the failure boundary and
  leaves no `failure.json`, despite an output directory having been created.
- If scientific computation raises `ToolError` and writing `failure.json` then
  raises `OSError`, callers receive the latter error instead of the scientific
  error. The diagnostic failure replaces the primary exception.
- `work/<key>/record.json` is saved before `write_exports` adds `export_id` to
  each prediction. Consequently, the retained record differs from its final
  `report.json` counterpart. Accepted short-rescue proposals can share these
  prediction objects and need consistent IDs as well.
- A record is appended before its JSON is successfully published. When record
  publication fails, `completed_sequences` includes a record that was not retained.
- `atomic_json` uses a fixed sibling `.tmp` name and does not clean it up if
  writing or replacing fails. This couples unrelated attempts to the same
  temporary file and can leave incomplete artifacts behind.

## Agreed implementation

### 1. Give publication a focused module

Move `atomic_json` and `write_exports` into a small `reporting.py` module. Keep
their existing imports from `crisprtracrrna.pipeline` working through explicit
re-exports. In particular, `pipeline.run` must continue resolving its
`atomic_json` and `write_exports` names through the pipeline module so existing
monkeypatch seams and failure-injection tests remain effective.

Introduce one private deterministic helper for assigning export IDs from each
record key and prediction ordinal. Call it before saving each record in `run`.
Also use it inside `write_exports` to preserve that function's existing behavior
for standalone callers, including its mutation of prediction dictionaries. Keep
the exact `<record key>.p<ordinal padded to six digits>` identity and the original
prediction ordering. Do not duplicate the ID construction rule or introduce a
new public report type/framework.

Move a record into the completed records list only after `record.json` has been
successfully published. The final report and retained record must agree, including
export IDs on accepted short-rescue proposals that share prediction objects.

### 2. Preserve publication and error contracts

After successfully creating a fresh output directory, protect work-directory
creation, settings serialization and configuration publication with the same
failure handling as analysis and final exports. Keep preflight and refusal to
overwrite existing output unchanged; failures must never modify an existing run.

On any existing supported failure or interruption, attempt to publish
`failure.json` with the existing keys. If that diagnostic operation fails, log
the secondary failure and re-raise the original exception unchanged. Do not let
failure-marker creation or temporary cleanup replace the primary exception.
`completed_sequences` must count successfully persisted per-record results.

Retain the existing order: record evidence and exports first, manifest next,
`report.json` last. Preserve the manifest's documented exclusion of itself and
the final report. A failed export, manifest or report publication must not leave
a successful top-level `report.json`.

### 3. Make atomic JSON attempts independent

Use a unique temporary sibling in the destination directory, followed by atomic
replacement. Serialize with the current UTF-8 encoding, indentation, trailing
newline and `allow_nan=False` behavior. Preserve ordinary file-creation
permissions where feasible; an exclusively opened unique filename obeys the
existing umask, while standard temporary-file helpers may change permissions.

Close and remove abandoned temporary files after serialization, write or replace
failures while preserving the original exception. A previously published
destination must remain unchanged if the new attempt does not reach replacement.
This is an atomic single-file publication helper; do not add fsync guarantees,
multi-file transactions or retry machinery.

## Ownership and affected files

The separate implementation agent owns:

- `new_version/src/crisprtracrrna/pipeline.py`
- New `new_version/src/crisprtracrrna/reporting.py`
- A dedicated focused test file such as `new_version/tests/test_publication.py`
- `new_version/docs/code_quality/03_pipeline_implementation.md`

Runtime/configuration and scientific parsing/search/evidence modules belong to
other councils and must not be changed by this implementation. Existing report
tests and native integrations are regression targets. Do not edit `old_version`,
external/native dependencies, previous benchmark evidence, local configuration or
model assets.

## Meaningful validation

1. A small mocked successful run covering multiple records, empty predictions,
   reused original sequence IDs and accepted proposal aliases must demonstrate
   that each persisted record equals its corresponding final report record.
   Verify prediction export IDs across JSON, CSV, FASTA and GFF3 and verify every
   manifest digest. Preserve CSV quoting and GFF/FASTA escaping checks.
2. Compare valid export bytes before and after extraction, including model-mode
   defaults and complete-mode evidence fields. The intended artifact difference
   is that saved per-record JSON now includes the already-public export IDs;
   associated manifest checksums therefore legitimately change.
3. Inject failures at work/configuration setup, a later record save, exports,
   manifest publication and final-report publication. Each failure must remain
   observable, must not publish a completion report, and should publish a failure
   marker when storage permits. A later record-save failure must report only the
   number of records successfully saved before it.
4. Inject a second failure while publishing `failure.json` and confirm callers
   receive the original exception, including interruption semantics, with a
   logged diagnostic about the secondary failure.
5. Exercise atomic JSON serialization/write/replace failures against an existing
   destination. Assert the prior content survives and abandoned temporary files
   are removed. Verify strict nonfinite-number rejection and cleanup-error
   handling do not mask the original failure. Verify a pre-existing sibling
   named `<destination>.tmp` is not overwritten by a new publication attempt.
6. Run the focused tests, existing process/report and council provenance tests,
   the complete dependency-free test suite, Ruff, formatting and mypy checks.
   The lead will also compare fresh native integration results with the frozen
   pre-change baseline; no benchmark evidence may be rewritten.

## Resolved tradeoffs and non-goals

- Preserve `write_exports` mutation for compatibility rather than silently
  changing its API. Explicit pre-publication assignment makes the lifecycle clear.
- The saved per-record artifact gains matching `export_id` values; final report
  schema and values remain stable. This intentional consistency fix is not a
  scientific behavior change.
- Keep `command_provenance` and its static module-resolution policy unchanged.
  Its existing tests cover interpreter versus module checkout distinctions, and
  the council found no reproduced defect requiring a provenance refactor here.
- Leave `analyze_record` biological stage order, rescue admission, scoring,
  sorting and native tool calls unchanged. Extracting those stages now would
  increase regression risk without addressing the demonstrated publication gaps.
- No new dependencies, durable filesystem transaction promises, parallel
  execution, process policy changes or broad report-schema typing migration.

No unresolved dissent remains. The architecture expert proposed separating
reporting and fixing error masking and snapshot inconsistency; the integrity
expert independently reproduced those issues, added completed-record accounting,
and explicitly accepted the implementation and validation plan above.
