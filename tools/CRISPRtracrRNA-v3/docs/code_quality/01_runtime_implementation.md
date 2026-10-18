# Runtime implementation report

Implemented the agreed runtime council plan on 2026-09-27 in a separate
implementation agent, without changing CLI flags, scientific calculations,
configuration defaults, dependency profiles or command-record keys.

## Changes

- `config.py`: small validation helpers consistently reject boolean, nonnumeric
  and nonfinite numeric settings and weights. Weights are checked as a mapping
  before inspecting keys or values. Configured model, model-file, executable and
  command values are validated before path operations, so malformed TOML values
  produce a setting-specific `TracrError`. Valid override precedence,
  config-relative paths, executable symlinks and literal flag arguments remain
  unchanged. Integer settings retain their strict integer checks.
- `models.py`: `Interval` and `Interval.from_public` reject noninteger coordinates,
  including booleans, with `ValueError` before comparisons or arithmetic.
- `sequence.py`: `subsequence` rejects unknown strands with `TracrError` before
  extracting sequence. Valid plus/minus behavior remains unchanged.
- `process.py`: `Runner.run` validates argv and computes one absolute child working
  directory. Explicit relative executable paths and relative/empty `PATH`
  entries are interpreted in that directory. The same absolute executable path
  is hashed, recorded and launched, retaining symlink identity. Failed commands
  read at most 12,000 bytes from the stderr tail, decode with replacement and
  expose at most 3,000 diagnostic characters while keeping the complete log.
  Empty argument strings remain valid; only the executable must be nonempty.
- `tests/test_runtime_contracts.py`: 131 focused parameterized regression cases
  cover validation failures, valid config paths and overrides, environment-value
  redaction, exact coordinate results, malformed argv, executable lookup/search
  order/digests, interpreter symlinks, bounded multibyte stderr reading and CLI
  aliases/booleans with patched scientific entry points.

## Validation

Run from `new_version` using the existing development interpreter, with no
dependency installation:

| Command | Result |
| --- | --- |
| `python3 -m pytest -q tests/test_runtime_contracts.py tests/test_sequences.py tests/test_process_and_reports.py` | 161 passed |
| `python3 -m pytest -q -m 'not integration'` | 243 passed, 7 deselected |
| `python3 -m ruff check src tests scripts CRISPRtracrRNA.py` | Passed |
| `python3 -m ruff format --check src tests scripts CRISPRtracrRNA.py` | Passed; 39 files formatted |
| `python3 -m mypy` | Passed; 12 source files |

Existing literal-shell-argument, failure-code, descendant-timeout and both-strand
candidate-clipping regressions passed. The new process tests use harmless
temporary executables with different caller/child implementations and assert
both observed output and the recorded SHA-256. A bounded-reader wrapper rejects
unbounded reads and verifies the exact final 3,000 characters when the byte window
starts inside a multibyte UTF-8 character.

The seven external-runtime integration cases were intentionally deselected for
this component's validation; the coordinating lead owns final integrated native
validation. Counts describe the shared workspace at the time of these checks.

## Deviations and scope

The lead's implementation review corrected one overly restrictive plan detail:
generic subprocess argument strings may legitimately be empty, for example
`[python, "-c", ""]` or an option's empty value. Runner therefore requires a
nonempty executable and string argument elements, preserving valid empty
arguments. Two regression cases verify empty Python code and empty payload
arguments. Settings command validation retains its existing nonempty-item
contract. This bounded correction preserves the plan's valid-command objective.

Existing test modules and CLI production code were left unchanged; all new
runtime coverage is in the dedicated test file. Protected
legacy, external/native, benchmark, local configuration and model asset files
were not edited. Deferred immutability, manifest atomicity and cancellation
refactors remain deferred as agreed by the council.
