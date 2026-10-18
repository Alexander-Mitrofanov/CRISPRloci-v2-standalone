# Runtime council consensus plan

Council: `runtime_architect` (architecture and API contracts) and
`runtime_reliability` (failure handling and regression coverage).
Both members reviewed the configuration, CLI, process, sequence and domain modules,
exchanged findings and explicitly agreed to this bounded plan on 2026-09-27.

## Objective and evidence

Improve validation and process reproducibility without changing valid scientific
calculations, defaults, command-line interfaces or report schemas. The council
reproduced these defects in the current implementation:

- `Settings(similarity=True)` and `Settings(timeout=True)` are accepted;
  `Settings(weights=None)` raises an incidental `TypeError`.
- TOML `erpin_model = 4` and `cmscan = 4` fail inside path manipulation;
  `models = [""]` silently resolves to the configuration directory.
- `Interval(0.0, 2.0)` and `Interval(False, True)` are accepted. An unknown
  `subsequence` strand silently selects the reverse complement.
- `Runner.run(["./tool"], cwd=directory)` looks up the executable in the caller's
  working directory. When both directories contain that name, it can hash one
  executable and launch another. An empty argv raises `IndexError`.
- Reporting a failed process reads its entire stderr file before slicing the last
  3,000 characters, despite scientific tools potentially producing large logs.

## Agreed implementation

### 1. Explicit configuration boundaries

Use small private validation helpers where they remove repetition in `Settings`;
avoid a new schema framework or dependency. Numeric settings and weights must
reject booleans, nonnumeric values and nonfinite numbers consistently. Keep
integer settings strictly integral. Validate weights as the expected mapping before
examining its keys and values, and preserve the existing required keys.

Ensure configured path values have their intended nonempty string/list shapes
before path operations; never turn an empty model or model-file path into a
directory. Invalid settings should raise `TracrError` with the setting name, rather
than leaking incidental `TypeError` or `AttributeError` exceptions. Keep the
existing distinction between data-path resolution and executable paths: executable
and interpreter symlinks must remain intact. Keep command arrays as argv arrays,
including flags and `--key=value` arguments, without interpreting them as shell
text. Preserve valid override precedence and config-relative paths.

No production CLI refactor is required. Cover its public behavior with focused
tests, including legacy direct-flag invocation, aliases and explicit booleans.

### 2. Validate coordinate and strand primitives

Require actual integer coordinates, excluding `bool`, at both `Interval`
construction and `Interval.from_public`; validate public inputs before arithmetic
can coerce them. Retain `ValueError` as the interval error contract and the existing
zero-based, half-open versus one-based, inclusive conversions. Reject strands
outside `+` and `-` in `subsequence` with `TracrError`, keeping both valid results
unchanged. Do not add broad validation to every domain dataclass in this change.

### 3. Resolve and report the executable that is actually launched

Validate `Runner.run` argv as a nonempty list/tuple of nonempty strings and report
malformed commands with `ToolError`. Compute the effective child working directory
once. Resolve explicit relative executable paths and relative/empty `PATH` entries
against that directory; retain normal absolute paths and search order. Pass the
same absolute, symlink-preserving executable path to hashing, the record and
`Popen`. In particular, do not call `resolve()` on a venv interpreter path.

Read failure diagnostics from a bounded tail of the stderr file using seek/read
and replacement decoding. Retain the existing 3,000-character diagnostic limit,
full stderr evidence on disk, return-code behavior and `ToolError` context. A byte
window large enough for the final 3,000 UTF-8 characters is acceptable; include
decoding-boundary coverage.

## Ownership and affected files

The separate implementation agent owns these production files:

- `new_version/src/crisprtracrrna/config.py`
- `new_version/src/crisprtracrrna/models.py`
- `new_version/src/crisprtracrrna/sequence.py`
- `new_version/src/crisprtracrrna/process.py`

Prefer dedicated focused test files (for example `test_runtime_contracts.py`) to
avoid conflicts with other implementation agents. Existing sequence/process tests
are regression targets; small additions there are also permitted if coordinated.
The CLI and errors modules were reviewed and need no planned production changes.

## Meaningful validation

1. Parameterize invalid numeric booleans, nonfinite values, weights shapes and
   malformed TOML paths; assert `TracrError`. Check valid overridden values,
   relative config paths, executable symlinks and environment-value redaction.
2. Assert interval constructor/public-constructor rejection of floats and booleans,
   invalid-strand rejection, and existing plus/minus clipping and coordinate
   round trips. Preserve exact values for valid coordinate operations.
3. Run harmless temporary executables with `./tool` and relative/empty `PATH`
   entries under a requested child cwd. Include different same-named tools in
   caller and child directories; verify output, recorded absolute executable path
   and SHA-256 all identify the child tool. Cover an executable symlink without
   replacing its path with the resolved target. Test malformed argv inputs.
4. A failing subprocess with a large multibyte stderr log must retain the final
   diagnostic marker and full evidence file while exposing at most the existing
   3,000-character diagnostic tail. Verify bounded reading directly if needed with
   a read wrapper that rejects unbounded reads.
5. Preserve existing literal shell-metacharacter, nonzero-exit and descendant
   timeout tests. CLI tests should patch `run`/`doctor` and assert settings reaching
   them, avoiding any external scientific runtime.
6. Run focused tests, then the complete dependency-free suite and repository lint
   checks. Record external-runtime skips honestly; no scientific benchmark rerun
   or dependency installation is required for these input-boundary changes.

## Non-goals and resolved tradeoffs

- No change to scientific algorithms, thresholds, scoring, model data, CLI flags,
  report/command-record keys, process cancellation policy or log filenames.
- No edits to `old_version`, external/native dependencies, existing benchmarks or
  local runtime configuration.
- Deep immutability of `Settings` dictionaries was considered but deferred: it
  would change public container types and serialization for little immediate
  benefit to this fix.
- Atomic command-manifest writes and cancellation refactors were considered but
  deferred to keep failure semantics stable. Existing timeout coverage remains.
- Rejecting previously tolerated malformed settings and coordinates is intentional;
  valid input behavior must remain unchanged. Both experts accepted this boundary.
- The relative-executable fix takes priority over stylistic refactoring because it
  repairs a demonstrated mismatch between provenance and execution.

No unresolved dissent remains. The architecture expert proposed the bounded
validation plan; the reliability expert added and reproduced the cwd/provenance
defect, and both explicitly accepted its inclusion and the deferrals above.
