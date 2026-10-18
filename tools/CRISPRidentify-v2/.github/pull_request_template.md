## Purpose

<!-- What problem does this solve? Link the issue. -->

## Behavior and scientific contract

<!--
Describe user-visible behavior, affected categories/coordinates/schema fields,
and whether the target is region recall, representation, or exact borders.
Write "not applicable" when this is not a scientific-result change.
-->

## Validation

<!-- List exact commands and summarized results. -->

- [ ] Complete unit suite passed.
- [ ] Relevant bundled-genome or integration smoke test passed.
- [ ] git diff --check passed.
- [ ] New or changed behavior has regression coverage.
- [ ] Documentation and schema material match the implementation.

## Reproducibility and performance

<!--
For scientific/performance changes, identify inputs and hashes, environment
lock, hardware, timing method, baseline, and output-equivalence contract.
Do not mix measurements with projections.
-->

## Safety, provenance, and licensing

- [ ] No credentials, private sequence data, environments, results, or caches are included.
- [ ] Imported code, data, and model artifacts have documented origin, license, and checksums.
- [ ] Vendored CasAndra manifests/notices are unchanged or were deliberately updated and validated.
- [ ] Security-sensitive details are being handled through the private process in SECURITY.md.

## Reviewer notes

<!-- Migration concerns, known limitations, or focused areas for review. -->
