# Independent tooling and packaging validation

Reviewer: `validation_tooling`, an independent final-council member who did not
plan or implement these changes. Review date: 2026-09-27.

## Scope and evidence

Compared `scripts/dependencies.py`, the new `scripts/check_wheel.py`, focused
maintenance/resource tests, CI, and maintenance/development documentation with
`/tmp/tracrrna-quality-baseline`. Packaging configuration is unchanged. The
validation used only temporary local repositories and local built distributions;
no real upstream checkout or dependency lock was changed.

Independent checks:

- `python3 -m pytest -q tests/test_dependency_updates.py tests/test_wheel_check.py`:
  **25 passed**. The tests cover detached and branch rollback, staging and
  publication failures, a checkout that changes HEAD before failing, explicit
  restoration failure, annotated tag resolution and subsequent maintenance,
  invalid refs, mismatched checkout/remote, failed sync and retry, existing
  destination preservation, and resource corruption/missing files.
- Ruff passed for both maintenance scripts and both focused test files.
- Direct inspection of the final wheel archive verified all **27 resource
  hashes**, the typing marker, the retained `py3-none-linux_x86_64` tag, and absence
  of external checkout content. The source archive contains the maintenance
  scripts and focused tests byte-for-byte as in the working source.
- The isolated-wheel smoke check was independently repeated while `PYTHONPATH`
  pointed to a shadow `crisprtracrrna` package whose import immediately raises.
  It passed all 27 resources and both entry points. The installed probes therefore
  ignored this inherited source-shadow location. The command log is
  `/tmp/tracrrna-wheel-shadow/check.txt`.
- The lead's integrated unit, native, typing, formatting, build, wheel and source
  distribution results were inspected as supporting evidence, rather than
  repeatedly rerunning unchanged integrated checks.

The wheel check proves installation, entry points and packaged resource integrity;
it does not execute scientific tools or certify compatibility with arbitrary
upstream revisions. The lack of concurrent-update/crash recovery is explicitly
documented and agrees with the council's bounded plan.

## Reproduced finding requiring correction

The initial reviewed implementation unconditionally called
`temporary.unlink(missing_ok=True)` from `update()`'s `finally` block. A cleanup
failure replaces the exception already being propagated, including the recovery
message when checkout restoration also fails.

Independent reproducer: `/tmp/tracrrna-tooling-cleanup-probe.py`, with output in
`/tmp/tracrrna-tooling-cleanup-probe.txt`. The probe creates a local upstream with
two commits, syncs the first into an isolated project, and updates toward the
second. It injects an `OSError` at lock replacement and a `PermissionError` when
removing the owned staged lock. A second run also injects failure of the checkout
restoring the first commit. No network or real dependency is involved.

Observed before correction:

| Injected failures | Propagated exception message | Checkout and lock |
| --- | --- | --- |
| Publication and cleanup | `CLEANUP: cannot remove staged lock` | HEAD restored; disk lock and caller mapping unchanged |
| Publication, restoration and cleanup | `CLEANUP: cannot remove staged lock` | HEAD at candidate; disk lock and caller mapping unchanged |

The full Python traceback retains the prior error in exception context. However,
the propagated exception's type/message loses the primary error and, in the
second case, the actionable restoration diagnostic. This violates the agreed
original-error/recovery contract. The appropriate bounded correction is to retain
an already propagating exception and add secondary cleanup details including the
staged path, with regressions for both combinations. Cleanup must not imply that
the checkout was restored when restoration failed.

## Council discussion and resolution

`validation_runtime` independently inspected the affected lines and probe output
and agrees that this is a concrete blocker; runtime/reporting review has no
additional blocker. `validation_science` accepts the scientific changes within
the tested software-equivalence scope and agrees final acceptance must await this
repair. All three reviewers agree on the bounded remedy above.

The initial disposition was conditional approval pending this correction. The
implementation agent subsequently preserved the pending exception and added a
cleanup-failure note containing the staged path. Successful publication clears
the consumed staging path, avoiding redundant cleanup.

The independent reviewer reran the original local-Git probe against the correction:

- Publication plus cleanup failure now propagates the original publication
  `OSError`; its note reports the cleanup failure and owned staged path. HEAD is
  restored, and disk lock/caller mapping remain unchanged.
- Publication plus restoration plus cleanup failure now propagates the combined
  `RuntimeError`, retaining the publication error as its cause and the original
  revision/manual recovery message. A note reports the cleanup failure and path.
  HEAD remains at the candidate as honestly reported; disk lock/caller mapping
  remain unchanged.
- The expanded focused suite passes **30 tests**; Ruff check and format check pass
  for both maintenance scripts and both focused test files. Added regression
  coverage exercises attached/detached states, successful/failed restoration,
  exception identity/cause, residue notes and successful publication without a
  redundant unlink.

The recheck output is `/tmp/tracrrna-tooling-cleanup-fixed.txt`. The earlier wheel
proof remains applicable because this correction changes maintenance scripts and
tests, not installed package code or packaging configuration. The lead refreshes
the source distribution and integrated evidence after these edits.

Final tooling disposition: **ACCEPTED**, with no remaining actionable blocker in
the reviewed scope. `validation_runtime` explicitly concurs in acceptance of the complete change
set after receiving the independent recheck. `validation_science` was notified
that its agreed tooling condition is satisfied; scientific evidence remains
applicable because the correction only affects maintenance tooling. Acceptance
remains within the documented software-validation scope.
