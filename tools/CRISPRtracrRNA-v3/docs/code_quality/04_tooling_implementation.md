# Developer tooling implementation

Executed the agreed [tooling council plan](04_tooling_plan.md) as a separate
implementation agent on 2026-09-27.

## Changes

- `scripts/dependencies.py` resolves annotated tags to commits, rejects empty
  refs, options and branch-writing refspecs before Git access, and stages a copied
  lock in a unique sibling temporary file before changing HEAD. It publishes the
  lock before updating the caller's mapping. A caught failure restores either the
  original local branch or the original detached HEAD without a destructive
  reset. Failed restoration reports both errors and the original recovery state.
- A missing upstream is cloned and checked out inside an owned temporary
  directory, then moved to its final path only when ready. Failed clones and
  failed pinned checkouts clean up their staging directory. Existing destinations
  and untracked build products in already pinned checkouts retain prior behavior.
- `scripts/check_wheel.py` creates a fresh temporary virtual environment and
  installs one wheel with `--no-index --no-deps`. Its installed-package probe and
  module entry point run with isolated Python from outside the source checkout.
  The installed console launcher runs directly with inherited Python path
  overrides removed and safe-path mode enabled. The check verifies environment
  import location, package/distribution version agreement, both entry points'
  help and version output, `py.typed`, the manifest, all 27 declared resource
  checksums, all five model/index sets and ERPIN executable permission on Linux.
- The existing CI Python matrix now invokes the isolated wheel check after its
  build. README development instructions document the same command; dependency
  documentation states the recovery guarantees and their limits.

## Focused validation

Commands run from `new_version/`:

```bash
python3 -m pytest -q tests/test_dependency_updates.py tests/test_wheel_check.py
python3 -m ruff check scripts/dependencies.py scripts/check_wheel.py tests/test_dependency_updates.py tests/test_wheel_check.py
python3 -m ruff format --check scripts/dependencies.py scripts/check_wheel.py tests/test_dependency_updates.py tests/test_wheel_check.py
```

Initial results: **25 tests passed**; Ruff check passed; all four Python files pass the
format check. Failure-injection tests cover partial staging writes, failed final
replacement, checkout failure after HEAD moves, and failed restoration. Both
attached-branch and detached-HEAD recovery are exercised. Local repository tests
also cover annotated tags and subsequent maintenance, dirty checkout refusal,
remote/HEAD mismatch, invalid refs before Git access, failed clone/checkout with
successful retry, existing-directory preservation, and untracked build products.
Resource tests detect missing assets/typing information, incomplete manifests,
corrupted models and lost executable permissions.

All dependency tests use temporary local Git repositories. No real upstream was
fetched, checked out, cleaned, reset, patched or rebuilt. Dependency pins,
scientific resources, packaging configuration and native sources are unchanged.

Integrated source/native/build validation and the actual built-wheel smoke result
are recorded by the lead in the overall validation record after all independent
implementations are ready; this implementer did not repeat those runs.

## Independent review correction

The final independent council reproduced a failure in the initial cleanup path:
an `OSError` while removing the staged lock replaced the original publication
error or the combined publication/restoration diagnostic. The bounded correction
preserves the pending exception and adds a note containing the failed cleanup
path and error. Successful publication clears the consumed staging path, so no
redundant cleanup can interrupt publication or prevent updating the caller's
mapping.

New regression cases combine publication and cleanup failures with both successful
and failed restoration, each from attached and detached HEAD states. They assert
the primary exception identity or combined diagnostic/cause, cleanup note and
residue, unchanged disk lock and caller mapping, and the resulting branch/HEAD.
A successful-publication case verifies that blocked unlink cannot leave the
caller mapping behind the saved lock and checked-out revision.

After the correction, the same focused command passes **30 tests**. Ruff check
and formatting pass for the two changed Python files. The independent tooling
reviewer reruns the original reproduction separately; the lead reruns the
integrated unit, lint and distribution checks. Production package code is
unchanged by this correction.

## Limits and deviations

No scope deviations. Recovery covers ordinary caught exceptions; it does not
provide crash atomicity or concurrent dependency maintenance. The packaging check
does not execute scientific tools or certify upstream compatibility. No packaging
configuration edits were needed during focused implementation.
