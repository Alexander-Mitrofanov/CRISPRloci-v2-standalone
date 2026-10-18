# Independent final council: scientific contracts and regression review

Reviewer: `validation_science`, a new council member who did not author the plans
or implementation. Reviewed on 2026-09-27.

Decision: **ACCEPTED; final council concurrence confirmed**, with no remaining
blocker in the reviewed scope. This is evidence of software equivalence for valid inputs
and improved rejection of malformed upstream records. It is not a validation of
biological sensitivity, specificity, or the calibration of experimental rules.

## Scope and code review

Compared `adapters.py`, `search.py`, and `evidence.py` against the saved original
source in `/tmp/tracrrna-quality-baseline/src` and inspected the new
`short_rescue_policy.py`, the science consensus plan, implementation report,
focused tests, and differential probe. The original source is also archived at
[baseline/source.tar.gz](baseline/source.tar.gz).

- Every experimental threshold literal and public threshold key, key order, and
  integer/float representation matches the baseline. Query span and identity
  minima remain inclusive; standard-coverage equality still excludes a short
  proposal; energy and CM/duplex overlap remain strict; model and terminator
  E-value maxima remain inclusive. Anchor membership, antiparallel pairing,
  support-position rules, and rejection-reason order are unchanged.
- Ranking expressions, tie ordering, model-rescue calculations, default settings,
  and native scientific command construction are untouched by these diffs. Runtime
  executable-resolution changes were separately reviewed by the runtime member.
- The shared JSON decoder rejects named nonfinite numbers and exponent overflow
  wherever they occur, including additive fields. Valid finite additive fields
  remain preserved. JSONL is checksum-verified before decoding, retains row order,
  skips blank physical lines, and gives file/line context for record failures.
- Consumed identifiers and counts now fail early with `ContractError` when invalid.
  Strict coordinate failures from the runtime change remain inside the adapter
  contract. Invalid repeat normalization also becomes a path-scoped contract error.
  Existing schema-version checks and ignoring unaccepted array categories remain
  unchanged. The scope intentionally does not become full upstream schema validation.

## Independently executed checks

1. Ran `python3 -m pytest -q tests/test_science_contracts.py
   tests/test_upstream_contracts.py tests/test_scientific_parsers.py
   tests/test_council_round2.py`: **118 passed**. This includes malformed numeric
   fields, consumed IDs/counts, checksum precedence, invalid coordinates/repeats,
   valid schema extensions, short-policy boundaries, and existing native-parser
   contracts.
2. Reran [science_differential_probe.py](science_differential_probe.py) against
   isolated baseline and current source imports. The output was byte-identical:
   **1,511,671 bytes**, SHA-256
   `d4ab4b54b629e53cfdeea4ca1aceb4dc65df584304eea8574d9c2729d1654d0d`.
   This independently reproduces the implementation evidence: 116 alignment hits,
   16 ordinary evidence/ranking records, 1,156 experimental decisions, valid
   CRISPRidentify/CasAndra adapters, and model-rescue cases.
3. Wrote and ran a separate reviewer probe with literal expected decisions,
   [validation_science_probe.py](validation_science_probe.py). It covers **672
   alignment cases** including both query/genome roles, all query/subject direction
   combinations, values immediately around 70%/80% identity, and span boundaries;
   **18 admission cases** with odd/even anchor lengths, half-anchor equality, and
   signed-zero/negative energy; **2 cassette-ID tie cases**; and **5 finite JSON
   extreme values**, including subnormal values and negative zero. Both baseline
   and current runs pass all assertions and produce byte-identical output:
   **121,027 bytes**, SHA-256
   `360c878d46f436feaa8f593c990a8e0cefbfb2f172e0588d19bf4211a3176165`.
   The reviewer probe passes Ruff lint and formatting.
4. Inspected the lead's completed native comparison and its comparison script.
   The five saved reports contain **104 predictions**; predictions and raw model
   evidence compare exactly, FASTA/GFF3 bytes match, CSV matches after temporary
   fixture-root normalization, and final manifest hashes and per-record JSON
   consistency checks pass. The native log records **7 integration tests passed**.
   This review did not independently rerun those external tools; the separate
   nonnative probes above were independently executed.

To reproduce the reviewer probe from `new_version`:

```sh
PYTHONPATH=/tmp/tracrrna-quality-baseline/src python3 docs/code_quality/validation_science_probe.py > /tmp/tracrrna-independent-science-before.json
PYTHONPATH="$PWD/src" python3 docs/code_quality/validation_science_probe.py > /tmp/tracrrna-independent-science-after.json
cmp /tmp/tracrrna-independent-science-before.json /tmp/tracrrna-independent-science-after.json
sha256sum /tmp/tracrrna-independent-science-before.json /tmp/tracrrna-independent-science-after.json
```

## Council discussion and limits

The runtime reviewer exchanged findings and independently accepted executable
resolution, failure handling, and atomic publication, with additional probes for
PATH/cwd/symlinks, process-group timeout, concurrent writers, and publication
failure. I explicitly accepted those findings and their stated limits: per-file
atomic replacement is not a multi-file transaction or an fsync durability promise.
No scientific concern arises from that limitation in this change set.

The tooling reviewer identified a reproduced blocker outside the scientific
modules: dependency-update cleanup masked the primary publication or rollback
failure. All three reviewers agreed that the primary error and recovery
diagnostics must survive cleanup failure. The bounded correction was implemented
and the tooling reviewer independently reran the original local-Git reproducer:
primary exceptions, their causes, rollback state and recovery diagnostics now
survive, with cleanup details and the residue path retained in exception notes.
The expanded 30-test focused suite and Ruff/format checks passed. I read the
completed [tooling re-review](07_validation_tooling.md) and explicitly concur with
its acceptance. **The agreed condition is satisfied; I accept the final council
conclusion with no remaining blocker in the documented review scope.** Scientific
modules and package configuration did not change during that correction, so the
independent scientific results above remain applicable without repeated testing.

The protected-tree audit has no modified or missing protected files and no asset
or configuration hash mismatch. After I flagged five newly generated Ruff cache
files under an archived benchmark directory, the lead removed exactly those
new files and recorded the cleanup in the audit. The updated audit has zero
remaining protected-tree additions; scientific assets remain hash-identical.

Finite fixture/native regressions cannot prove equivalence over all possible
inputs or establish biological accuracy. The experimental gate remains an
experimental hypothesis with the same thresholds as before; no improved predictive
performance or calibrated probability is claimed. Malformed data rejected earlier
is an intentional contract change, so equivalence applies to valid inputs.
