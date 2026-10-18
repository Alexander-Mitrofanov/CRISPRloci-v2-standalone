# Reference container schema 1.0.0

This document defines a software container and numerical scanner. It does not
establish the biological validity or permission status of any supplied data.
There is intentionally no synthetic default reference pack. Missing production
references cause an actionable `ReferenceError`, not a biological no-call.

`references.load_bundle(path=None)` takes a directory or its `manifest.json`.
Omitting the path addresses installed `crisprleader/data/default/`. Prediction
never fetches data. The directory must contain exactly the regular data files
listed in `files`, plus `manifest.json`; symlinks, unsafe paths, extra files,
missing files, executable formats, duplicate JSON keys, and nonfinite JSON
numbers fail closed. Version 1.0.0 accepts at most 8 MiB of manifest JSON,
64 MiB per member, 256 MiB total, 4,096 files, 4,096 families, and 256 motifs.
Every member's SHA-256 is verified and its verified bytes retained in memory.
Checksums provide content identity and corruption detection, not authentication.

The manifest has exactly these top-level keys:

| Field | Contract |
| --- | --- |
| `schema_version` | String `1.0.0`. |
| `bundle_id`, `version` | Portable identifiers matching `[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}`. |
| `families` | List of family objects below, with unique IDs. |
| `motifs` | List of generic PWM objects below, with unique IDs. |
| `provenance` | Object with `build` object, `sources` list, `rights` list, and `limitations` string list. Additional fields allowed here. |
| `files` | Object mapping relative POSIX paths to lowercase SHA-256 digests. The manifest does not list itself. |
| `policy` | Explicit policy metadata object; the container does not invent or calibrate a scientific decision rule. |

Family objects require `id`, nonempty unique uppercase IUPAC DNA string lists
`repeat_sequences` and `core_sequences`, `boundary_kind`, nonempty unique
`training_locus_ids` and `citations` string lists, and a `metadata` object.
`boundary_kind` is one of `comparative_conserved_segment`,
`experimentally_interrogated_segment`, or `unresolved`. These are declared
semantics, not an independent verification of the curator's claim.
Optional fields are `description`, `calibration` object, and paired
`model_file`/`model_format`. Model files must be checksum-listed data;
formats are `hmmer3-dna` or `json`. Loading never executes or interprets them.
`metadata` can carry per-reference IDs, accession versions, coordinates,
sequence hashes, boundary meaning, reviewed evidence, calibration provenance,
declared applicability and limitations. Neither model membership nor a motif
match establishes experimental function.

Motif objects require:

| Field | Contract |
| --- | --- |
| `id` | Portable identifier. |
| `domain` | `leader`, `junction`, or `repeat`. This only sets the eligible sequence coordinates. |
| `pwm` | 1–128 position rows, each an object with exactly `A`, `C`, `G`, `T`, positive probabilities summing to one. Zero counts require explicit, provenance-recorded smoothing before import. |
| `background` | Exactly `A`, `C`, `G`, `T` probabilities, positive and summing to one. Fixed IID null on the candidate's forward-oriented sequence. |
| `strands` | Unique nonempty list containing `+`, `-`, or both, relative to the candidate coordinate frame. |
| `citations` | Nonempty unique string list. |
| `metadata` | Object for source, applicability, restrictions, and other annotations. |

Optional `score_scale` is an integer from 1 through 1,000 (default 100).
Optional `report_pvalue` is in [0,1] (default 0.001); this is a reporting
threshold after the full test universe is scored. Optional
`junction_end_range: [minimum, maximum]` restricts the interval's end offset
from the flank/repeat junction to that inclusive range. It is a predeclared
coordinate restriction; no post hoc positional selection occurs.

The returned plain dictionary adds two reserved keys that are not allowed in
on-disk manifests: `_receipt`, containing manifest SHA-256, file checksums,
schema/bundle/version identities and verified bytes; and `_members`, mapping
member names to their verified immutable bytes. Use `reference_receipt(bundle)`
for an independent JSON-safe receipt copy. Exclude reserved keys when exporting
the original manifest. Do not serialize `_members` as scientific output.

## Generic motif scanner

`scan_motifs(candidates, motifs, *, max_tests=1000000) -> list[MotifHit]` scans
all explicitly supplied candidates and models in one multiple-testing family.
`candidates` contains the frozen `FlankCandidate` records; `motifs` is the
validated list of PWM objects. Missing actual observed terminal-repeat sequence
fails; consensus repeats are never substituted. Hit coordinates index
`candidate.sequence + candidate.observed_repeat`. Candidate orientation and
motif strand are separate.

Leader-domain windows are wholly inside the flank. Repeat-domain windows are
wholly inside the observed terminal repeat. Junction windows must contain at
least one base on each side. The optional end-offset restriction applies in
that same candidate coordinate frame. A reverse-strand match uses the
reverse-complemented PWM with the same forward-sequence background, so an
asymmetric background receives a distinct correct null distribution.

Each position's log2 likelihood ratio is rounded to the nearest integer
`score_scale` unit (ties to even). Dynamic programming enumerates the entire
distribution of this declared integer score under the fixed IID background.
`score` is the quantized score in bits; the unquantized log2 score is retained
separately in `details`. `pvalue` is the upper-tail probability **of that
quantized statistic**, exact up to floating-point arithmetic. It is not a
probability of function, activity, a leader, or correct orientation. Matrices
requiring more than 100,000 null states fail rather than silently approximate.
Combinations with a minimum possible null-string probability below `1e-300`
fail rather than claim underflowed zero p-values.

Benjamini–Hochberg adjusted `qvalue` values include all models, candidates,
eligible positions, and selected strands. Ambiguous IUPAC windows remain in
the family with p=1 but are not emitted as occurrences. The scanner does not
assume that overlapping windows meet any particular independence condition;
the values are explicitly BH adjustments under the declared background, not
an empirically demonstrated false-discovery guarantee for genomic sequence.
Reporting thresholds apply after adjustment. A test-bound overflow fails
before scoring, so output caps cannot shrink the statistical denominator.

Every reported hit includes its matched sequence, background, raw and
quantized scores, rounding rule, model citation, scope digest, number of
tests, ambiguous-test count, and model/candidate content identities. Empty
hit lists establish no biological absence. This scanner never makes a
leader decision. Exhaustive short-string synthetic tests validate numerical
and engineering behavior only; biological benchmarks are a separate task.
