# Contributing to CRISPRidentify v2

Thank you for improving CRISPRidentify. Contributions are expected to preserve
scientific traceability, deterministic reports, and fail-closed behavior at
external-tool and model boundaries.

This repository is currently private. These instructions apply to authorized
collaborators and do not imply public contribution access.

## Before starting

Choose the right route:

- Read [SUPPORT.md](SUPPORT.md) for installation and usage questions.
- Follow [SECURITY.md](SECURITY.md) for vulnerabilities. Never put exploit
  details, credentials, or sensitive sequence data in a public issue.
- Open a scientific-result issue for an unexpected array or Cas call.
- Open a performance issue for a reproducible slowdown.
- Discuss broad algorithm, schema, model, or dependency changes before
  implementing them.

Small fixes can go directly to a pull request when their scope and expected
behavior are unambiguous.

## Development setup

Use the locked Linux x86-64 environment:

```bash
./install.sh
.envs/crispridentify-core/bin/python scripts/check_environment.py --profile core
```

Create a focused branch from the current target branch. Do not mix generated
benchmark data, unrelated formatting, and functional changes in one pull
request.

The vendored runtime and large model files are intentionally ordinary Git
files. Do not replace them with Git LFS pointers.

## Tests

Run the complete tracked suite after installation:

```bash
.envs/crispridentify-core/bin/python -m unittest discover -v
git diff --check
```

Also run the bundled-genome command when changing the CLI, installation,
native code, reporting, CasAndra integration, or output transactions:

```bash
./crispridentify \
  --file examples/NC_017040.1.fasta \
  --result_folder runs/contributor-smoke \
  --parallel False \
  --cas True \
  --html_report True
```

Use a fresh result directory for the smoke test. A change is not ready merely
because a narrow unit test passes. Add regression coverage for the failure mode
and run the complete suite.

Documentation-only changes should still pass `git diff --check`, link review,
and any documentation contract tests.

## Scientific-result changes

A proposed change to candidate generation, region competition, array
boundaries, category assignment, Cas classification, or confidence semantics
must include:

1. the biological question and decision contract;
2. source accessions or redistributable minimal fixtures plus sequence hashes;
3. exact commands, software commit, environment lock hash, and hardware;
4. before/after results using a declared coordinate and overlap convention;
5. separate counts for regions, candidate representations, and exact borders;
6. evidence for both recovered calls and newly introduced calls;
7. tests that fail closed on malformed or incomplete evidence;
8. a clear statement of whether evaluation records overlap training or
   threshold-development data.

Do not describe source-reconstruction validation as biological validation.
Do not describe `certainty_score` as a calibrated probability. Report the
five canonical categories—Bona-fide, Possible, Alternative, Possible
Discarded, and Low score—separately where they are emitted, and state which
categories define the comparison target.

Any new trained artifact needs a model card, training/evaluation provenance,
license, immutable checksum, and a leakage audit. Deployment refits must not be
presented as independent evaluation models.

## Performance changes

Performance claims must identify:

- input manifest and total sequence size;
- cold/warm cache policy and number of repetitions;
- wall time, CPU time, and peak resident memory;
- CPU model, allocated cores, RAM, storage, and operating environment;
- the exact command and concurrency layout;
- output-equivalence or region-recall contract;
- failures, timeouts, and excluded inputs.

Do not multiply speedup ratios that use different baselines. Keep prototypes
and projected throughput clearly separate from measured end-to-end results.

## Code and documentation style

- Follow the existing Python style and keep changes narrowly scoped.
- Prefer explicit validation and typed failure over silent truncation or
  fallback.
- Keep user-facing coordinates one-based and end-inclusive unless a field
  explicitly documents another convention.
- Preserve deterministic ordering and transactional publication of
  `report.json`.
- Wrap user-facing prose for readable diffs.
- Do not mechanically reformat `vendor/casandra/`.

The repository may add automated formatters over time. The checked-in
configuration and CI are authoritative when they differ from this prose.

## Vendored CasAndra changes

`vendor/casandra/` is a modified downstream snapshot covered by
`SOURCE_MANIFEST.sha256`. Its supplied source had no immutable upstream VCS
revision or published artifact locator.

Coordinate any edit under that directory with a maintainer. A valid change must
preserve the Apache-2.0 materials, update the downstream modification record,
update all affected manifests and checksums, and pass the package-integrity
tests. Never invent or infer an upstream commit.

## Pull requests

A reviewable pull request:

- has one clear purpose and links the relevant issue;
- explains user-visible and scientific behavior changes;
- lists exact validation commands and results;
- updates README, schema, and changelog material when contracts change;
- contains no environments, result directories, caches, logs, or private data;
- preserves third-party notices and model/data attribution.

Reviewers may request an independent reproduction for biological or
performance claims. Maintainers decide release inclusion after CI and review;
opening a pull request does not guarantee acceptance.

## Conduct and licensing

Participation is governed by [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
Contributions must be legally redistributable and compatible with the licenses
recorded in [LICENSE](LICENSE) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
State the origin and license of any imported code, model, or dataset.
