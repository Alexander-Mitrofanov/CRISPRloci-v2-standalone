# Scientific adapters and evidence implementation

Implemented the consensus in [02_science_plan.md](02_science_plan.md).

## Changes

- `adapters.py` uses one JSON-object decoder for structured reports and CasAndra
  JSONL records. It rejects non-object records, named nonfinite constants and
  finite-literal overflow, including values nested in additive fields. Finite
  additive fields remain intact.
- CasAndra JSONL is read in UTF-8 one physical line at a time after checksum
  verification. Record validation failures identify the file and one-based line.
  Input FASTA hashing occurs once per parse invocation. Cassette IDs and accepted
  CRISPRidentify array IDs must be nonempty strings; consumed CasAndra counts must
  be nonnegative integers, excluding booleans.
- Invalid upstream consensus repeats now produce a path-scoped `ContractError`.
  Strict interval errors introduced by the runtime implementation stay inside the
  same adapter boundary. Unaccepted CRISPRidentify categories and unused optional
  upstream fields retain their previous handling.
- A frozen `ShortRescuePolicy` in `short_rescue_policy.py` supplies the proposal
  filter, experimental admission comparisons and public threshold dictionary.
  Every threshold, comparison boundary, report-key order, numeric representation
  and rejection-reason order is unchanged. Unrelated ranking constants remain
  untouched.
- Added 65 focused cases in `tests/test_science_contracts.py` for the external
  contracts, preserved extensions, coordinate error translation, and literal
  experimental gate output.

## Validation

Run from `new_version`, using `python3` from the available conda runtime:

- Focused upstream/scientific/round-two suite plus new tests: **118 passed**.
- Full suite with `-m 'not integration'`: **308 passed, 7 deselected** at the time
  of this implementation. This includes the completed runtime additions; the
  pipeline implementation was still in progress.
- Strict mypy: **no issues in 13 source files** at the time of the check.
- Ruff lint over `src tests`: **passed**.
- Ruff lint and formatting over all six science-owned Python files, including the
  differential probe: **passed**. The full formatting check identified the new
  pipeline-owned `reporting.py` while its executor was still editing it; this
  unrelated file was left to that owner for final formatting.
- Native integration reruns are coordinated by the lead developer.

## Valid-output differential evidence

[science_differential_probe.py](science_differential_probe.py) executes against
isolated source imports without requiring native tools. The same probe was run
against the original source snapshot and the changed package. Its JSON output was
**byte-identical**, preserving dictionary order and integer/float representation,
across:

- Three CRISPRidentify fixtures covering both strands and unknown orientation,
  normalized repeats, additive fields and ignored categories.
- A valid multirow CasAndra report with UTF-8 metadata, finite additive fields,
  blank lines, and three nearest-cassette selection results.
- Four standard/experimental alignment matrices, using default and customized
  standard thresholds, producing 116 total hits across both strands.
- Sixteen standard candidate evidence/ranking records, including absent evidence
  and forward/reverse native hits on both genomic strands.
- 1,156 experimental admission decisions spanning both genomic strands, disabled
  and enabled experimental settings, exact paired-base, energy, model E-value,
  duplex-overlap, terminator E-value and position boundaries, and absent duplexes.
- Four model-rescue result lists, including empty and nonempty outputs on both
  genomic strands.

The common output is 1,511,671 bytes with SHA-256
`d4ab4b54b629e53cfdeea4ca1aceb4dc65df584304eea8574d9c2729d1654d0d`.
Counts and digest are saved in
[science_differential_summary.json](science_differential_summary.json).
The complete local outputs are `/tmp/tracrrna-quality-science-before.json` and
`/tmp/tracrrna-quality-science-after.json`.

Reproduce from `new_version`:

```sh
PYTHONPATH=/tmp/tracrrna-quality-baseline/src python3 docs/code_quality/science_differential_probe.py > /tmp/tracrrna-quality-science-before.json
PYTHONPATH="$PWD/src" python3 docs/code_quality/science_differential_probe.py > /tmp/tracrrna-quality-science-after.json
cmp /tmp/tracrrna-quality-science-before.json /tmp/tracrrna-quality-science-after.json
sha256sum /tmp/tracrrna-quality-science-before.json /tmp/tracrrna-quality-science-after.json
```

## Deviations and limits

No implementation-scope deviation. JSONL field-contract failures receive the same
file/line context as decoding failures, improving diagnostics without extending
validation to unused schema fields. Rejection of malformed nonfinite values,
consumed identifiers/counts and coordinate types is intentional; valid scientific
outputs match the baseline probe. The differential fixture probe complements,
but does not replace, the lead developer's native regression rerun.
