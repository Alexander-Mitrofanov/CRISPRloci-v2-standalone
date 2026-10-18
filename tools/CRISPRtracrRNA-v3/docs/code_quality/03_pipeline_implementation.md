# Pipeline and reporting implementation

Implemented by the separate `pipeline_implementer` executor on 2026-09-27,
following the consensus in `03_pipeline_plan.md`.

## Changes

- Extracted JSON and CSV/FASTA/GFF3 publication into `reporting.py`. Explicit
  re-exports retain `pipeline.atomic_json` and `pipeline.write_exports`; `run`
  continues calling those module names so existing failure-injection seams work.
- Centralized deterministic export-ID assignment in one private helper. `run`
  assigns IDs before writing each retained record; standalone `write_exports`
  preserves its existing mutation behavior. Accepted rescue proposal aliases
  therefore have the same IDs in retained and final JSON.
- Count records as complete only after their individual JSON publication succeeds.
- Moved work-directory creation, settings serialization and configuration writing
  inside the failure boundary. Preflight and refusal to overwrite remain outside.
- Failure-marker errors are logged without replacing the original exception,
  including `KeyboardInterrupt` and `SystemExit`.
- Atomic JSON now exclusively creates a UUID-named sibling using normal umask
  permissions, closes it before replacement, and cleans abandoned attempts. Cleanup
  errors are logged without masking the primary error. An old fixed `.tmp` sibling
  is never touched. Strict nonfinite-number rejection and serialization bytes are
  unchanged.
- Added 23 focused tests in `tests/test_publication.py`, covering consistency,
  publication ordering, manifest digests, export bytes and escaping, setup and
  publication fault injection, retained-record accounting, secondary diagnostic
  failures, atomic serialization/write/replace/cleanup failures, unique siblings,
  ordinary creation permissions, and compatibility re-exports.

`analyze_record`, `doctor`, and `command_provenance` have identical ASTs to their
frozen pre-change implementations. Biological stage order and scoring were not
modified. No legacy, upstream, native, model, configuration, or historical
benchmark files were edited by this executor.

## Validation

Executed from `new_version`:

| Check | Result |
| --- | --- |
| `python3 -m pytest -q tests/test_publication.py tests/test_process_and_reports.py tests/test_council_controls.py tests/test_council_round2.py` | 70 passed |
| `python3 -m pytest -q -m 'not integration'` | 349 passed, 7 deselected at executor check time |
| `python3 -m ruff check src tests/test_publication.py` | Passed |
| `python3 -m ruff format --check src tests/test_publication.py` | 15 files already formatted |
| `python3 -m mypy` | Passed; 14 source files |

Additional direct comparisons loaded the original functions from
`/tmp/tracrrna-quality-baseline/src/crisprtracrrna/pipeline.py`:

1. Compared all three export files byte for byte, plus the exporter-mutated
   prediction dictionaries, using mixed model-mode defaults and complete evidence,
   empty output, and reused sequence IDs. All comparisons matched. The mixed
   evidence bytes are independently retained as literals in the regression test.
2. Ran the frozen pipeline and changed pipeline on the same mocked scientific
   result with timestamps and durations pinned. `report.json`,
   `configuration.json`, `summary.csv`, `predictions.fasta`, and
   `predictions.gff3` were byte-identical. The only retained-record difference was
   the intended addition of export IDs to predictions and aliased accepted
   proposals. Their manifest digests consequently change as intended.
3. Compared parsed function ASTs for `analyze_record`, `doctor`, and
   `command_provenance`; all were unchanged.

The lead developer owns the final combined checks, native integration comparison,
and independent validation council. Whole-suite counts above are a snapshot while
other implementation agents were adding their tests.

## Deviations and limitations

No implementation scope deviations. A repository-wide `ruff check .` also scanned
protected native and downloaded dependency code and produced unrelated findings;
no automatic fixes were applied there. Scoped source checks passed, and the lead
will use the project's documented first-party check scope.

Atomic publication does not add fsync or a transaction across the export files.
Storage failures can prevent failure-marker or temporary-file cleanup; these
secondary failures are logged while preserving the original error. The final
report remains the completion marker published after exports and manifest.
