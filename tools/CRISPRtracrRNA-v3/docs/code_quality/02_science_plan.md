# Scientific adapters and evidence council consensus plan

Council: `science_architect` (architecture and maintainability) and
`science_contracts` (external data contracts and scientific regression protection).
Both experts independently reviewed `adapters.py`, `search.py`, `evidence.py` and
their tests, exchanged findings, and explicitly agreed to this plan on 2026-09-27.

## Objective and reproduced findings

Improve external-data error handling and reduce drift between the existing
experimental gate and its reported thresholds. Preserve all valid biological
predictions, command arguments, thresholds, defaults, ranking, ordering and public
report schemas. This is an engineering change; no scientific rules are being
recalibrated.

The council reproduced these defects:

- `load_json` rejects literal `NaN` and `Infinity`, but accepts `1e999` as an infinite
  float. CasAndra JSONL bypasses even the existing constant check.
- A checksum-valid CasAndra cassette without `cassette_id` is accepted, then
  `nearest_cas` raises an incidental `KeyError` when the cassette is used.
- Invalid or empty CRISPRidentify consensus repeats raise `TracrError` from
  normalization rather than the adapter's path-scoped `ContractError`.

Separately, the frozen short-rescue gate repeats numerical policy across the
search filter, admission decision and public `thresholds` dictionary. These values
currently agree, but maintaining separate copies permits future drift.

## Agreed implementation

### 1. A consistent, bounded external JSON boundary

Factor a small strict JSON-object decoder inside `adapters.py`, shared by the
existing `load_json` entry point and the CasAndra JSONL reader. Reject non-object
records, nonfinite named constants, and floating-point overflow such as `1e999`.
Keep regular finite numbers and all legitimate additive fields intact; do not
whitelist report keys, reject unknown fields, or introduce a schema dependency.

Read JSONL with explicit UTF-8 and line iteration rather than materializing and
splitting the whole file. Preserve row order and blank-line handling. Keep checksum
verification before parsing. Malformed JSONL errors should identify the source
file and one-based physical line number, while preserving `ContractError` as the
public adapter failure. Checksum failures must remain identifiable as such.

Compute the input FASTA SHA-256 once per `parse_casandra` call rather than once for
each inspected input-manifest entry. Preserve the existing requirement that some
input entry matches both digest and `genome_fasta` kind.

Validate only identifiers and counts the tool consumes: accepted CRISPRidentify
array IDs and CasAndra cassette IDs must be nonempty strings; the existing
`genes`/`cas_cassettes` counts must be nonnegative integers excluding booleans.
Maintain existing source, coordinate, prediction-boolean, classification-object
and count-equality checks. Do not add requirements for unused protein IDs or
optional upstream fields. Continue ignoring unaccepted CRISPRidentify categories
as before.

Translate normalization failures for upstream repeat sequences into
`ContractError` with the CRISPRidentify report path. The runtime council is making
interval coordinates strictly integral while retaining `ValueError`; ensure
upstream coordinate failures still reach callers through the adapter's
`ContractError` boundary. Avoid a broad catch of every exception.

### 2. Single-source the frozen short-rescue policy

Introduce one small immutable internal policy value in a dependency-free module
(for example `short_rescue_policy.py`). A frozen dataclass with literal field
defaults and a method producing the public dictionary is suitable. Both
`search.parse_alignment_table(short=True)` and
`evidence.short_rescue_admission` must use this same policy.

The exact public fields and values remain:

| Field | Value |
| --- | --- |
| `min_query_span` | `12` |
| `min_identity` | `0.8` |
| `min_anchored_paired_bases` | `10` |
| `min_anchored_paired_fraction` | `0.5` |
| `energy_strictly_below` | `0` |
| `cm_evalue_max` | `0.01` |
| `cm_duplex_overlap_strictly_below` | `10` |
| `terminator_evalue_max` | `0.001` |

Preserve every comparison's inclusivity, rejection-reason ordering, paired-base
calculation, report key order and integer/float representation. Standard coverage
and identity remain governed by the existing `Settings`. Do not repurpose these
experimental-policy fields for unrelated ranking bonuses merely because some
numbers happen to match. Do not expose a new user option or modify defaults.

## Ownership

The separate implementation agent owns `adapters.py`, `search.py`, `evidence.py`,
the small new policy module, and focused new tests such as
`tests/test_science_contracts.py`. Existing scientific tests may be extended if
needed, but separate files minimize overlap with other agents.

Do not edit runtime-owned `config.py`, `cli.py`, `process.py`, `sequence.py` or
`models.py`, or pipeline/reporting modules assigned to another council. Do not
change `old_version`, external/native dependencies, benchmarks, model assets or
local configuration files.

## Meaningful validation

1. Exercise checksum-valid CasAndra reports with `NaN`, `Infinity`, negative
   infinity and exponent overflow, including a nonwinning/additive record field.
   These malformed numerical values must fail at the input boundary. Confirm
   valid finite additive fields remain present and unchanged in returned data.
2. Cover non-object JSONL records, malformed JSON on a later physical line, blank
   lines and UTF-8 data; assert filename/line context for malformed records.
3. Cover missing, empty and non-string consumed IDs, boolean/fractional/negative
   counts, malformed repeats and invalid coordinate types. Assert `ContractError`
   rather than downstream incidental errors. Retain valid additive schema updates,
   checksum-tampering detection, row ordering and existing schema-version rules.
4. Preserve existing short-proposal and short-admission boundary tests. Add only
   missing meaningful boundaries, such as nonnegative energy at exactly zero, and
   assert the reported thresholds against the literal pre-change dictionary.
   Assertions should test externally visible decisions, not merely compare the
   new policy object to itself.
5. Run focused upstream/scientific/round-two tests, then the complete nonintegration
   suite, Ruff lint/format checks and strict mypy. The environment's `python3`
   (`/home/alex/miniconda3/bin/python3`) provides these tools; the project `.venv`
   does not contain pytest. Run from `new_version`.
6. Demonstrate valid-output equivalence against the saved source baseline at
   `/tmp/tracrrna-quality-baseline`, including default and experimental admission
   outputs, thresholds, ranking values, both strands, and valid adapter fixtures.
   Root will coordinate the complete native regression rerun; the baseline passed
   all seven available integration tests and stored native artifacts there.

## Discussion and resolved tradeoffs

The architecture expert proposed shared strict JSON handling and frozen policy;
the contracts expert independently found the JSONL hole, confirmed missing-ID and
normalization failures, and added early validation of consumed IDs/counts and
literal threshold regression checks. Both explicitly accepted this scope.

We deferred splitting the scientific module, streaming all native text-table
parsers, generic ranking helpers, and full upstream-schema validation. These would
increase churn without addressing the demonstrated failures. Rejecting malformed
nonfinite JSON and invalid consumed IDs/counts is intentional; valid native
outputs must remain unchanged. No unresolved dissent remains.
