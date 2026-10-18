# Pre-result protocol clarifications

2026-09-27, before new native panel results.

- Detection council accepts the frozen contract. Freeze the full cross-locus repeat
  donor pool in the control manifest. Shuffle seeds are specimen-stable; changing
  the donor universe is a different derangement panel and must not be silent.
- Engineering council accepts the frozen contract. Plan section 2's default batch
  limit of 128 means the configured positive `intarna_batch_size`; users may set
  other positive sizes. The frozen comparison still uses 1 versus 128.
- Installed IntaRNA 3.2.0 supports `seedTRange` only with one target. The executor
  therefore uses individual anchored experimental calls, as plan section 3 permits.
  Standard unanchored candidates remain batched by exact oriented repeat.
- Existing export IDs encode rank. Default IDs and order remain exactly equivalent.
  Experimental insertions may shift rank-derived export IDs; stable candidate IDs
  and existing biological fields remain preserved, and export agreement is audited.

These clarify implementation/provenance, not thresholds, data, or acceptance gates.

Validation council also finds no blocking protocol error. Its pre-result
clarification is adopted: the genomic-control review examines EVERY newly admitted
short-rescue candidate supported by CM and/or terminator, including
`sequence_duplex_cm`, `sequence_duplex_terminator`, and
`sequence_duplex_cm_terminator`, rather than choosing a 'strongest' tier after
looking at results. Any such new corroborated control admission triggers explicit
council review and keeps default promotion blocked. This does not imply it is a
biological false positive, and no default promotion is planned.

Coordinator found an existing provenance defect before new panel results:
`command_provenance()` assigns a Python interpreter's enclosing Git checkout
revision to `python -m casandra`, although the imported CasAndra source is selected
through the configured PYTHONPATH. Historical reports show interpreter checkout
`784cb6d...`; the protected loaded CasAndra checkout is `ba046f90242ac91dc075183ea51dcbedfac28484`.
The executor is asked to correct or explicitly relabel this metadata and retain
actual module provenance under the child environment, with tests. This is a
metadata correctness correction within the provenance gate; biological settings,
models, candidate rules and evaluation thresholds remain frozen.
