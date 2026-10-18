# Descriptive transcription-orientation evaluation

`crisprleader.evaluation.evaluate(dataset, predictions, bootstrap_replicates=2000,
seed=0)` reports supplied annotation tables. It does not run a predictor,
collect biological labels, validate a scientific reference bundle, or establish
leader boundaries. This is secondary reporting infrastructure: conserved-core
localization and its scientific validation require a separate study.

All tests for this module use synthetic fixtures. Every output explicitly carries
`data_kind`, `synthetic`, and `scientific_superiority_established: false`. Setting
`data_kind` to `biological` is a caller assertion, not evidence of biological
validity. A positive confidence interval does not change the superiority flag.

## Dataset document

The exact top-level fields are:

```json
{
  "schema_version": "1.0.0",
  "dataset_id": "example-synthetic",
  "data_kind": "synthetic",
  "orientation_semantics": "transcription_orientation",
  "cases": [
    {
      "id": "case-1",
      "group_id": "independence-group-1",
      "stratum": "example-stratum",
      "truth": "+",
      "label_kind": "experimental",
      "label_source": "Synthetic fixture; no biological evidence"
    }
  ]
}
```

`data_kind` is `synthetic` or `biological`. `truth` is `+`, `-`, or null.
`label_kind` is `experimental`, `reviewed_comparative`, `computational`, or
`unknown`. Every field is required. Unknown fields fail closed at every level.

Only cases with a known truth and an `experimental` or `reviewed_comparative`
label are eligible for correctness scoring. Eligible cases require a nonempty
`label_source`. Other cases remain in `case_audit`, every tool's `case_results`,
and separate `unscored` summaries; a matching computational label never becomes
correctness evidence. A missing truth remains unscored even if its label kind is
experimental. The software does not independently verify the source, label, or
independence of a group.

The fixed semantics are **transcription orientation**. Acquisition orientation,
leader side, and transcription orientation are different annotation targets.
Tables using another target must not be relabeled to fit this schema.

## Prediction document

```json
{
  "schema_version": "1.0.0",
  "dataset_id": "example-synthetic",
  "tools": [
    {
      "id": "example-tool",
      "version": "fixture-v1",
      "reference_id": null,
      "predictions": [
        {
          "case_id": "case-1",
          "status": "abstained",
          "orientation": null,
          "reason": "Synthetic abstention example"
        }
      ]
    }
  ]
}
```

`status` is `called`, `abstained`, `unsupported`, or `failed`. `called` requires
`+` or `-`; every other status requires a null orientation. `reason` is required
but may be an empty string. `reference_id` is a nonempty string or null. At least
one tool is required. Missing rows become the distinct report outcome `missing`;
`missing` cannot be supplied as a prediction status. Unknown case IDs, repeated
case IDs, repeated tool IDs, and repeated predictions for a tool/case fail closed.
All dataset cases, including unscored cases, are valid prediction targets.

## Common denominators and summaries

Each tool's `eligible.overall`, `eligible.by_stratum`, and
`eligible.by_label_kind` contain the same fields. For the relevant eligible case
set, let `N` be its size and `C`/`I` be correct/incorrect called predictions:

| Field | Definition |
| --- | --- |
| `denominator` | `N`, fixed by the dataset, independent of tool output |
| `counts` | Separate correct, incorrect, abstained, unsupported, failed, missing counts; sum equals `N` |
| `coverage` | `(C + I) / N` |
| `conditional_accuracy` | `C / (C + I)` |
| `conditional_error` | `I / (C + I)` |
| `correct_yield` | `C / N` |
| `incorrect_yield` | `I / N` |
| `completion` | `(C + I + abstained + unsupported) / N` |

A zero denominator produces null, never zero or NaN. Conditional accuracy/error
are null when there are no calls. `unsupported` denotes a completed response
outside a tool's supported scope; it increases completion, not call coverage.
Failed and missing predictions remain in the denominator. They are never silently
removed to make an accuracy or completion result look better. An entirely
missing or failed eligible run has no conditional accuracy and is not an
available comparator. Its bookkeeping correct yield is zero on the fixed
denominator, but that is not evidence of zero biological accuracy.

Breakdowns include every declared stratum and label kind, with zero eligible
counts where appropriate. They do not constitute separate hypothesis tests.
Unscored summaries report `called_unscored`, abstained, unsupported, failed, and
missing counts, with `accuracy: null`. They never compare predictions to their
weak or missing labels.

## Equal-group comparison and paired bootstrap

The primary **descriptive comparison in this report** is the difference in
equal-group mean correct-call yield. It is not the primary scientific endpoint
of the proposed leader-localization study.

For each eligible group `g` and tool `t`, compute
`y(t,g) = correct(t,g) / eligible_cases(g)`. A tool's
`primary_group_correct_yield` is the unweighted mean of these group yields.
Each comparable `pairwise_comparisons` entry reports the first tool's mean minus
the second tool's mean, ordered by tool ID. This differs from case-weighted correct yield
when groups have different sizes. Missing, failed, unsupported, and abstained
predictions contribute zero correct calls while keeping their cases in each
group denominator. All tools are compared on exactly the same eligible cases
and groups. Whole-group membership may span strata.

The companion `group_weighted` summary keeps the comparison's weighting visible
for every metric. `group_count` and `eligible_cases` give its scope. For each
eligible group, `group_yields` records counts of correct calls, incorrect calls,
all calls, and completed outcomes, plus their fractions of eligible group cases.
`group_weighted.coverage`, `correct_yield`, `incorrect_yield`, and `completion`
are the equal-group means of those fractions, including groups with zero calls.
Its conditional accuracy is `mean(correct fraction) / mean(call fraction)`;
conditional error is `mean(incorrect fraction) / mean(call fraction)`. These
ratios are **not** the average of per-group accuracies/errors among called groups
and need not equal the pooled conditional metrics in `eligible.overall`.
If mean call coverage is zero, both conditional metrics are null. With no
eligible groups, every group-weighted rate is null.

For example, groups containing one correctly called case, four cases with one
incorrect call and three non-calls, and two cases with no calls have equal-group
coverage `(1 + 1/4 + 0) / 3 = 5/12`, correct yield `1/3`, incorrect yield `1/12`,
and conditional accuracy `(1/3)/(5/12) = 4/5`. Their pooled conditional accuracy
is `1/2`. The zero-call group remains in the equal-group denominators.

A tool with no completed eligible outcome (`called`, `abstained`, or
`unsupported`) has `comparison_eligible: false` and a null
`primary_group_correct_yield`. Its group-level counts and `group_weighted`
bookkeeping remain visible, including zero correct yield; these do not override
the comparator eligibility gate or establish zero biological accuracy. Any pair
involving it has a null difference and interval, `comparison_status:
not_comparable`, and `comparison_unavailable_reason:
no_completed_eligible_outcomes_for_tool`. Partially missing/failed runs remain
descriptively comparable on the full common denominator, and
`all_eligible_predictions_completed` exposes their incompleteness. No eligible
cases instead produce the reason `no_eligible_cases`. These declarations check
submitted rows only: they do not verify which external comparator tools should
have been submitted or whether the comparator registry is complete.

For `G >= 2` eligible groups, each bootstrap replicate samples `G` whole group IDs
with replacement. Every case in a sampled group is retained, and repeated groups
retain their sampling multiplicity. The same sampled groups are used for every
tool and every pair. Group IDs and tool IDs are sorted before using a local
`random.Random(seed)` generator. Percentile endpoints use linear interpolation
at `p * (replicates - 1)` for `p = 0.025, 0.975`. They are unadjusted descriptive
95% intervals. There are no hypothesis-test p-values or superiority decisions.

For otherwise comparable tools, fewer than two eligible groups produce a null
interval and an explicit reason; no eligible groups also produce a null point
difference. A single tool produces
no pairwise comparisons. `bootstrap_replicates=0` explicitly disables intervals.
Small positive replicate counts are permitted for diagnostics, but fewer than
1000 replicates produce a Monte Carlo uncertainty warning. Such intervals,
especially those from very few groups or replicates, are not sufficient for
confirmatory inference. Degenerate data can legitimately yield a zero-width
percentile interval; that does not establish certainty or scientific superiority.

The report distinguishes the requested replicate count (`bootstrap_replicates`)
from `bootstrap_replicates_computed`, which is zero when no interval can be
computed. Case IDs, tools, strata, and groups are deterministically sorted. Reordering
input lists does not change the report, semantic input hashes, or intervals for
the same seed. This reproducibility contract applies to the same implementation
and Python runtime; the report does not claim identical random-number behavior
across all possible future implementations.

## Bounds, receipts, and interpretive limits

The module enforces the following limits without silently selecting a subset:

- 10,000 cases, 16 tools, and 100,000 tool-by-case cells.
- IDs, versions, strata, and reference IDs: 256 characters each; label sources
  and reasons: 4096 characters each.
- 16 MiB of UTF-8 text across parsed input fields; invalid Unicode is rejected.
- Integer bootstrap counts from 0 to 10,000; integer seeds in `[0, 2**64)`.
  Boolean values are not accepted as integers.
- At most 20,000,000 declared resampling work units, calculated as
  `replicates * (groups * (tools + 1) + tool_pairs)` when intervals can be
  computed for at least one comparable pair. This is a conservative declared
  work budget covering all submitted tools and pairs. Exceeding this cap raises
  an actionable error. The caller may
  explicitly request fewer replicates or zero; no group is dropped automatically.

The API receives in-memory JSON-style dictionaries; callers remain responsible
for bounded file reading and duplicate-key detection before parsing. All fields
have explicit scalar or collection types, and nonfinite numerical substitutes
are rejected. Input objects are not mutated. No external dependency, network
access, inference model, or global random state is used.

The report contains SHA256 receipts of canonical semantic JSON with sorted case
and tool IDs, UTF-8 text, sorted object keys, and compact separators. These are
not raw file-byte hashes, and they do not authenticate biological labels or
upstream tool executions. The CLI additionally records `input_document_receipts`
with `canonical_json_sha256` for each parsed input document: sorted object keys,
compact separators, and ASCII-escaped JSON preserve the supplied list order.
Those receipts intentionally differ from the evaluator's semantic receipts,
which normalize case/tool list order and use UTF-8 text. Neither receipt type
claims a raw file-byte hash or requires rereading a changing input file. CLI
metadata also names and versions the evaluator and identifies its analysis kind.

Any scientific comparison still needs curated evidence-backed labels, an audit
of independence/grouping and training/calibration/test leakage, preregistered
estimands and margins with power or precision targets, multiplicity decisions,
and actual complete comparator runs with frozen configurations and references.
The report lists these prerequisites. Merely filling `label_source`, assigning
group IDs, declaring `biological`, or observing a favorable unadjusted interval
does not satisfy them.
