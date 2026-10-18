# Implementation status — 27 September 2026

Owner approval was received. The private repository is
<https://github.com/Alexander-Mitrofanov/CRISPRleader-v2>.
This is an incomplete development snapshot, not a validated leader predictor.

## Implemented and tested

- Source-bound FASTA/native input and canonical Identify adapter.
- Both context windows, circular arcs, clipping and strand-aware projection.
- Independent source/array/flank/result validation and explicit format conversion.
- Atomic no-clobber result export with output file hashes.
- Offline data-only reference containers and generic PWM numerical scanning.
- Shell-free bounded subprocess groups, timeouts and cancellation.
- CLI, packaging, synthetic examples, software tests and CI configuration.
- Descriptive evaluation of supplied transcription-orientation labels/calls,
  with complete outcome accounting and paired equal-group bootstrap intervals.

## Unfinished requirements

- The scientific family retrieval, scoring and abstention policy.
- A useful, provenance-reviewed default biological reference bundle.
- Reproducible reference construction, calibration and comparison of methods.
- Independently curated biological labels and execution of all frozen gates.
- Claims supported by measured biological performance and coverage.

Automated safety restrictions stopped biological research in two worker turns.
The coordinator stopped that work and did not route it through another worker or
acquisition method. Software infrastructure was completed where unaffected.
The restriction was not a failed biological benchmark or a withdrawal of the
owner's approval. No benchmark result or passing scientific gate is claimed.

`predict` therefore fails explicitly. The successful `extract` command states
that no scientific inference was performed; it never labels an extracted window
as a leader. Synthetic fixtures and passing software tests cannot fulfill the
approved scientific completion criteria.

## Readiness decisions

| Gate | Status |
| --- | --- |
| Owner confirmation | Received |
| Private repository | Created |
| Software components | Implemented; see [engineering checks](ENGINEERING-CHECKS.md) |
| Default biological references | Missing |
| Leader prediction | Unavailable |
| Biological evaluation | Not run |
| Scientific release | Not ready |
| Server integration / deployment | Separate scope; not performed |

Retain these criteria until independently satisfied. Do not replace missing
biological references with a synthetic default or a permissive fallback.

## Scientific-value follow-up

The comparison literature audit is in `benchmarks/COMPARATORS.md`.
`benchmarks/SUPERIORITY_PROTOCOL.md` strengthens the original prerelease gates
without relaxing them. The proposed primary scientific task remains conserved-core
localization; the implemented orientation reporter is secondary infrastructure.
No comparator has been run and no biological accuracy improvement is demonstrated.
The study design remains incomplete, with practical margins, power/sample size,
label review and final comparator/reference identities unresolved. Descriptive
bootstrap intervals and software tests cannot satisfy those requirements.
