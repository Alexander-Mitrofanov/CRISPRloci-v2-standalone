# CRISPRtracrRNA v3 maintenance

This checkout owns the typed, standard-library tracrRNA orchestrator. Scientific
programs run as external processes; upstream applications live in their own
repositories and environments.

## Code map

- `cli.py`, `config.py`: public commands, effective settings and preflight values.
- `sequence.py`, `models.py`: normalized DNA, immutable domain values and intervals.
- `adapters.py`: independently validate CRISPRidentify and CasAndra evidence.
- `search.py`: native commands, strict streamed parsers and bounded IntaRNA batches.
- `evidence.py`, `short_rescue_policy.py`: evidence selection, ranking and frozen gates.
- `pipeline.py`: per-record orchestration, provenance and completion publication.
- `process.py`: subprocess isolation, cancellation, logs and executable digests.
- `reporting.py`: atomic JSON and consistent FASTA/CSV/GFF3 export identities.
- `tests/`: portable contract regressions; opt-in real native-tool integrations.

## Scientific invariants

- Internal intervals are zero-based half-open; public intervals are one-based
  end-inclusive on the input forward strand. Strand remains a separate value.
- Candidate sequences and local evidence are in transcript orientation. Reverse
  model/terminator hits remain raw evidence and cannot support that hypothesis.
- Missing evidence stays null; absent hits are distinct from failed computation.
  Ranking scores are uncalibrated. Search windows are not native transcript bounds.
- Preserve thresholds, weights, model assets, default branches, deterministic
  tie/order rules, command arguments and schemas during engineering refinements.
  Experimental short rescue stays off by default and retains its frozen gates.
- Validate every native/upstream row, including nonwinning evidence. Never hide
  failures with retries or treat malformed tool output as an empty prediction.
- Hash retained evidence, publish `report.json` last, and refuse existing outputs.
  Diagnostics must not obscure the original scientific failure or cancellation.
- Do not merge scientific environments, download dependencies during a job, or
  add runtime configurations/private operations evidence to version control.

## Checks

See README Development for portable pytest, Ruff, strict mypy and wheel checks.
Native tests require `CRISPRTRACRRNA_TEST_CONFIG` pointing to explicit installed
commands: `pytest -q -m integration`. Reuse existing isolated environments.
For parser/evidence optimization, assert exact ordered scientific output against
reference fixtures; retain uncertainty, coordinate reconstruction and zero-hit cases.
