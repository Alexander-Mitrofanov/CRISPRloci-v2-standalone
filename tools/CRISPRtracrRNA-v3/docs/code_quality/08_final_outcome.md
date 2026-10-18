# Code-quality council: final outcome

**Accepted.** Four expert councils agreed on bounded plans, separate implementation
agents executed each plan, and a new three-member validation council reviewed the
integrated changes. All three final reviewers explicitly accepted the completed
change set after one independently discovered error-handling defect was corrected
and rechecked. No actionable blocker remains in the reviewed scope.

## Completed changes

| Area | Result |
| --- | --- |
| Runtime and inputs | Consistent settings, coordinate and strand validation; child-directory-aware executable lookup with matching provenance; bounded stderr diagnostics. Valid empty subprocess arguments and interpreter symlinks remain supported. |
| Scientific contracts | Shared strict JSON/JSONL decoding, contextual upstream errors, early validation of consumed identifiers/counts, and one immutable source for existing experimental rescue thresholds. |
| Pipeline and reporting | A focused reporting module; consistent prediction IDs in retained records and final exports; accurate completed-record counts; independent atomic writes; original failures preserved when diagnostic publication also fails. |
| Development tooling | Annotated tags resolve to commits; failed dependency updates restore checkout state where possible; failed syncs leave no incomplete destination; isolated wheel checks verify imports, entry points and packaged resources. |

All [plans and implementation records](README.md) identify participants, ownership,
decisions, tests and scope limits. Scientific thresholds, ranking, defaults, native
command behavior, public CLI flags and final report schemas are preserved. Invalid
settings and malformed upstream data now fail earlier and more consistently.

## Independent council findings

- [Runtime/publication reviewer](05_validation_runtime.md): accepted after
  independent executable-path, hash, timeout, concurrent-publication and
  multi-failure checks.
- [Scientific reviewer](06_validation_science.md): accepted after independently
  reproducing the original differential comparison and adding 697 comparison
  cases for search boundaries, admission, cassette selection and numeric values.
- [Tooling/packaging reviewer](07_validation_tooling.md): found that temporary-file
  cleanup could mask a dependency-update error and its recovery instructions.
  The executor corrected this and added five regression cases. The reviewer
  reran the original reproduction and accepted the repair. Runtime and science
  reviewers explicitly concurred in final acceptance.

The initial finding and corrected behavior are retained in
[before](validation/tooling-cleanup-before.txt) and
[after](validation/tooling-cleanup-after.txt) probe outputs. This validation was
not limited to reviewing passing test counts.

## Final verification

| Check | Result |
| --- | --- |
| Dependency-free unit suite | **360 passed**; baseline was 112 |
| Native-tool integrations | **7 passed**, using the existing reproducible runtime profile |
| Native baseline comparison | **104 prediction records across five reports unchanged**; FASTA/GFF3 bytes unchanged; CSV unchanged after normalizing temporary input-directory prefixes |
| Retained output integrity | Per-record JSON matches final records; every retained manifest digest verifies |
| Scientific differential probe | 1,156 experimental admission decisions plus standard ranking, adapter, search and model-rescue cases produce byte-identical baseline/current output |
| Ruff and formatting | Passed for all 45 maintained Python files |
| Strict mypy | Passed for all 14 package modules |
| Source distribution | Rebuilt after the final correction; **360 unit tests passed** from the extracted archive |
| Wheel | Isolated installation and both entry points pass; **27 resources** verify; all 14 installed module files match final source bytes |

The last correction changed only maintenance code, its tests and documentation.
Runtime module hashes remained identical to the native-tested source, so native
checks were not needlessly repeated. The independent wheel check also passed
with a deliberately broken shadow package on the inherited `PYTHONPATH`.

Detailed [validation summary](validation/summary.json),
[native comparison](validation/native-comparison.json),
[source hashes](validation/source-hashes.json), and
[before/after source diff](validation/changes.patch) accompany the check logs.
The [pre-change source archive](baseline/source.tar.gz) preserves the review
baseline; no usable Git repository was exposed in this workspace.

## Preservation and limits

The final [preservation audit](validation/protection.json) found no additions,
missing files, or size/mtime changes among **113,233 pre-existing protected
files** in historical code, upstream/native dependencies and prior benchmarks.
SHA-256 checks also verified **39 model, fixture and configuration files**.
Five newly generated Ruff cache files in an archived directory were removed;
their paths are recorded. Real dependency pins and checkouts were not updated.

This establishes code correctness and preserved behavior within the tested local
scope, not new biological accuracy or generalization. The full historical
literature benchmark was not rerun. Existing scientific dependency/platform
limits remain. Dependency-update recovery handles ordinary caught failures;
crash-atomic Git/filesystem transactions and concurrent maintenance remain outside
the agreed scope. The preserved CI Python matrix was not executed remotely.
