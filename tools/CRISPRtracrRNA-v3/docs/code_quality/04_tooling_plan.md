# Developer tooling council consensus plan

Council: `tooling_architect` (maintainability and packaging) and
`tooling_reliability` (dependency reproducibility and failure handling).
Both experts inspected the development tooling independently, exchanged concrete
findings, and agreed to the bounded plan below on 2026-09-27. A separate
implementation agent executes this plan.

## Objective and observed defects

Make explicit upstream maintenance recoverable after ordinary failures and make
the installed-wheel check exercise the actual distribution. Keep scientific
behavior, models, native binaries, dependency pins and benchmark evidence intact.

The architecture expert reproduced two defects using temporary local Git
repositories, without contacting or changing the real upstream checkouts:

- An injected lock-file replacement failure leaves the checkout at the new
  revision while the on-disk lock retains the old revision. The passed lock
  mapping is already changed and the fixed-name temporary file remains behind.
  The reliability expert independently identified the same failure ordering.
- Updating to an annotated release tag records the tag-object ID from
  `FETCH_HEAD`, although `git checkout` resolves it to a commit. A successful
  command therefore produces a checkout/lock mismatch that prevents the next
  ordinary sync or update.

Both experts also found that `sync` clones directly into its final destination;
failure to check out an unavailable locked commit can leave a poisoned checkout
there. Existing CI installs a wheel in the development environment and checks
only `--version`; it does not prove isolation from the source checkout or verify
that scientific resources survive packaging. Baseline sdist and wheel builds,
sdist unit tests, and an ordinary wheel installation already pass. The new check
adds proof of isolation and installed assets; it does not claim those builds
were broken.

## Agreed implementation

### 1. Preserve dependency-update consistency

Keep the current CLI and its explicit candidate-update semantics. Continue to
require a clean checkout, matching remote and matching pinned HEAD. Reject empty
refs, option-like refs and colon-containing fetch refspecs before fetching: the
documented argument accepts a tag, branch or commit, and must not move another
local branch as a side effect. Resolve the fetched result to a verified commit
with `FETCH_HEAD^{commit}` before recording or checking it out.

Prepare a copy of the proposed lock, including a copied dependency entry. Write
its complete JSON to a uniquely named temporary file beside the final lock
before changing the checkout. Preserve the original HEAD and whether it was
attached to a local branch. Check out the new commit, atomically replace the
lock file, then update the caller's mapping only after successful publication.

For caught failures before publication, remove only the temporary file owned by
this invocation and restore the original branch/HEAD if checkout was attempted.
Never use a destructive reset. Keep the original lock bytes and caller mapping
unchanged on failure. If restoration itself fails, report the original failure
and the restoration failure together with enough state for manual recovery;
never imply the update completed or that restoration succeeded.

This is recovery for caught failures, not crash atomicity across Git and a JSON
file. Do not introduce a transaction framework, persistent journal or concurrent
update protocol. Document this limitation concisely.

### 2. Publish synchronized checkouts only when ready

For a missing destination, clone into an owned temporary directory under the
same parent, check out the requested pinned commit there, and publish the
complete checkout at the final destination only after success. Clean up that
owned temporary directory on clone or checkout failure so a corrected subsequent
sync can succeed normally. Preserve any pre-existing destination.

Keep current existing-checkout protections and the deliberate behavior that
untracked upstream build products do not invalidate a successful no-op sync.
Do not clean, reset, patch or rebuild real upstream repositories. Broad lock
schema redesign and Git configuration hardening are outside this change.

### 3. Check the installed wheel independently

Add a small reusable distribution smoke check and wire CI to run it against the
built wheel in a fresh temporary virtual environment. Install the wheel without
dependencies or index access, execute from outside the source tree with isolated
Python (`-I`), and verify:

- The imported package comes from that environment, with package and distribution
  versions agreeing.
- Both the installed console command and `python -m crisprtracrrna` accept
  `--version` and `--help`.
- The installed manifest is present and every declared resource exists with its
  expected SHA-256 digest; `py.typed` and model/index paths are present.
- The bundled ERPIN executable retains executable permission on the supported
  Linux platform.

Use the simplest reusable script/CI split that keeps these checks legible. Keep
the existing Linux x86-64 distribution tag and supported Python matrix. Do not
run scientific tools or download native dependencies in this packaging check.
Document the exact local command in the development instructions. Full source
tests remain separate because some intentionally inspect maintenance scripts and
benchmark fixtures that are not wheel contents.

## Ownership

The implementation agent may edit:

- `new_version/scripts/dependencies.py` and a focused new packaging-check script.
- `new_version/tests/test_dependency_updates.py` and focused tooling tests if
  needed.
- `new_version/.github/workflows/ci.yml`.
- The development section of `new_version/README.md` and relevant maintenance
  guarantees in `new_version/docs/DEPENDENCIES.md`.
- A new `new_version/docs/code_quality/04_tooling_implementation.md` record.

Production package modules belong to the other councils. Packaging configuration
was reviewed; no planned edits to `pyproject.toml`, `setup.cfg` or `MANIFEST.in` are
needed unless the isolated check exposes an actual omission. Coordinate any such
addition with the lead before expanding file ownership. Leave existing benchmark
scripts, artifacts, upstream pins, local configuration and native sources alone.

## Acceptance checks

1. Use temporary local Git repositories for all dependency tests. Preserve the
   existing successful update and dirty-source refusal cases.
2. Update through an annotated tag and assert that recorded `commit`, checked-out
   HEAD and the underlying commit ID agree. Subsequent sync/update must work.
3. Inject both staging-write and final-replacement failures. Assert unchanged
   disk lock bytes and caller mapping, the original branch or detached HEAD,
   clean tracked worktree, and no owned temporary-file residue. Cover restoration
   failure diagnostics separately without mutating a real checkout.
4. Verify unsupported refs fail before a Git fetch. Preserve remote-mismatch and
   locked-HEAD mismatch refusals with focused coverage.
5. Fail sync using a nonexistent pinned commit; assert the final destination is
   absent and a corrected retry succeeds. Verify pre-existing directories remain
   intact, and an existing pinned checkout with untracked build products is a
   no-op that preserves those products.
6. The implementer runs focused tooling tests and relevant Ruff checks and
   formatting. The lead runs the integrated dependency-free suite, strict mypy,
   native regressions, sdist/wheel builds, and the new isolated wheel check.
   Record the actual wheel and check results without repeating successful
   integrated runs unless a later change warrants it.

## Consensus and limits

Both experts prioritize demonstrated consistency defects over broad script
refactoring. The reliability expert added restoration of attached branches,
explicit dual-failure diagnostics and rejection of branch-writing refspecs; the
architecture expert reproduced the annotated-tag object/commit mismatch. Both
accepted isolated wheel verification and deferred speculative lock-schema work,
benchmark changes, crash recovery and concurrent maintenance support. No
unresolved dissent remains.
