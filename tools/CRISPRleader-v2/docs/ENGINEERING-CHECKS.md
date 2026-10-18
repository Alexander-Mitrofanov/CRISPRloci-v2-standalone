# Engineering verification — 27 September 2026

Verified locally on Linux x86_64, Python 3.12.7. These are software checks, not
measurements of biological accuracy or scientific release readiness.

| Check | Result |
| --- | --- |
| `python -m pytest -q` | 293 passed |
| `python -m ruff check src tests scripts` | Passed |
| `python -m mypy src` | Passed, 16 source files |
| `python -m build --no-isolation` | sdist and wheel built |
| `python scripts/check_wheel.py` | Passed, offline installation into a fresh environment outside the checkout |
| Scientific prediction / independent biological gates | Unavailable / not run |

Tests cover input binding, strict formats, circular arcs and strands, observed
junction projection, numerical PWM nulls, full multiple-testing scope, reference
file integrity, process failure/timeout/cancellation, atomic publication races,
export reconstruction, malformed records and explicit CLI failure behavior.
The installed-wheel check verifies extraction, result validation, file receipts,
packaged data/schema resources, descriptive orientation reporting and prediction's
explicit unavailable error. The added evaluator has 59 synthetic tests, including
hand-computed denominators, paired resampling, missing-run eligibility and
imbalanced/zero-call group accounting.

`requirements-dev.lock` records the engineering dependency versions. No biological
packages or models are required by the runtime package. GitHub Actions is
configured for Python 3.11, 3.12 and 3.13; its remote results are separate from the
local results recorded here.

The parent CRISPRloci application and its active tool lock were not changed.
This snapshot was not integrated or deployed.

The original 228-test engineering snapshot preceded the scientific-value follow-up.
The current 293-test count includes the evaluator and CLI/publication regressions.
Independent methodology review checked label semantics, paired denominators and
missing-run gates. No biological prediction improvement was measured.


## Refinement verification — 13 October 2026

Verified on the existing Linux/Python 3.12.7 development runtime. The baseline
passed 316 tests; the refined version passes 320. Lint, static types (16 source
files), whitespace checks, sdist and wheel builds pass. An extracted-wheel CLI
smoke test ran outside the checkout under isolated Python with all source-tree
import paths removed. It verified extraction, independent validation, output
receipts, packaged resources, descriptive evaluation and prediction's explicit
unavailable error without creating another environment.

The scanner now prepares matrices/backgrounds once, retains the bounded null-table
cache, encodes candidate bases once, and uses prefix counts for ambiguous windows.
It slices reported occurrences without recreating full junction strings. Record
serialization copies containers once while preserving null values, scores and
source-forward segments. The code navigation entry point is [CODE_MAP.md](CODE_MAP.md).

Exact comparisons with the previous implementation passed 120 randomized synthetic
cases comprising 2678 eligible tests and 298 reported hits, including ambiguity,
both strands, all domains and restricted junction offsets. Full hit records,
statistical scope receipts and JSON serialization matched. The synthetic extraction
report matched after removing its creation timestamp. New exhaustive-enumeration
and record-snapshot checks supplement existing circular-coordinate/export controls.
These engineering results do not evaluate the biological acceptance gates.
