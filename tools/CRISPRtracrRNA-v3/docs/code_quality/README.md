# Code-quality council review

**Completed and accepted by all three independent final reviewers.** See the
[final outcome and validation results](08_final_outcome.md).

This engineering review covers the maintained v3 tool in `new_version/`.
It preserves scientific thresholds, ranking, default behavior, public command
interfaces, and report schemas. Historical code, native and upstream dependency
sources, existing benchmark evidence, models, and local runtime configurations
are outside the editing scope.

## Method and ownership

Each area has two expert reviewers who inspect independently, exchange findings,
resolve tradeoffs, and write a consensus plan. A separate implementation agent
then executes that plan. After integration, a new council reviews the combined
changes independently of the authors and planning councils.

| Area | Consensus plan | Separate implementation |
| --- | --- | --- |
| Runtime and inputs | [Architecture + reliability](01_runtime_plan.md) | [Runtime changes and checks](01_runtime_implementation.md) |
| Scientific contracts | [Architecture + contract validation](02_science_plan.md) | [Scientific changes and equivalence](02_science_implementation.md) |
| Pipeline and reporting | [Architecture + data integrity](03_pipeline_plan.md) | [Publication changes and checks](03_pipeline_implementation.md) |
| Developer tooling | [Architecture + reproducibility](04_tooling_plan.md) | [Maintenance and package checks](04_tooling_implementation.md) |

The work is sequenced to avoid overlapping file ownership. Every agreed plan
records its acceptance checks. The final evaluation records any unresolved limits.

## Starting baseline

- Unit suite: **112 passed, 7 integration tests deselected**.
- Ruff checks and formatting: passed (38 Python files).
- Strict mypy: passed (12 package modules).
- Source, scripts, tests, and build configuration were copied and hashed before
  edits in `/tmp/tracrrna-quality-baseline/`. This workspace does not expose a
  usable Git repository, so before/after review uses that snapshot.
- Native integration baseline: **7 passed**, using `local.reproducible.toml`.
- Source distribution and wheel builds passed; all 112 unit tests also passed
  from the extracted source distribution. An isolated wheel install reported
  version 3.0.0 outside the checkout.
- [Baseline evidence](baseline/summary.json) includes logs and source hashes.
  [The pre-change source archive](baseline/source.tar.gz) preserves the reviewed
  implementation without local configurations or dependency assets.

This review does not establish new biological accuracy or independent scientific
generalization; it evaluates software correctness and preserved behavior.

## Independent final council

Three new reviewers inspected the integrated work independently of all planning
and implementation agents:

- [Runtime and publication review](05_validation_runtime.md).
- [Scientific equivalence review](06_validation_science.md).
- [Dependency tooling and packaging review](07_validation_tooling.md).

The council found a secondary-cleanup error that could obscure dependency-update
recovery instructions. A separate executor corrected it and the tooling reviewer
independently verified the repair; all three reviewers then accepted the result.
The original reproduction is preserved in
[`validation/tooling-cleanup-before.txt`](validation/tooling-cleanup-before.txt).

Integrated check logs, native-result comparison, preservation checks and the
source diff are retained in [`validation/`](validation/).
