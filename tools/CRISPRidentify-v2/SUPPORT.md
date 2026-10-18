# Support

This page explains how to get useful help without exposing sensitive sequence
data or losing the information needed to reproduce a result.

## Supported environment

The supported installation target is Linux x86-64 using the release's
locked Conda environment. Other operating systems, architectures, dependency
solves, and modified model bundles may be useful experiments, but maintainers
may ask you to reproduce a problem in the locked environment.

## Check locally first

From the extracted release directory, collect:

```bash
./crispridentify --version
.envs/crispridentify-core/bin/python scripts/check_environment.py --profile core
./crispridentify --help
sha256sum locks/conda-linux-64-core.lock
```

Also retain the SHA-256 supplied with the release archive.

Then retry the smallest affected FASTA with a new result directory. Do not
overwrite or delete the original evidence while diagnosing a transactional
output failure.

Common first checks:

- Confirm that the input is readable FASTA and that folder inputs have unique
  result keys.
- Confirm that at least 5 GiB was available during installation.
- Confirm that `CRISPRIDENTIFY_ENV`, if set, points to the intended absolute
  environment prefix.
- For CasAndra, include `--cas True` and check the configured
  `--cas_timeout_seconds`.
- For IS-element annotation, install and validate the `is-elements` profile.
- For a machine-profile run, remove HTML, FASTA, pickle, and legacy JSON export
  options; those outputs are intentionally incompatible with machine mode.

## Choose the right issue

| Need | Route |
| --- | --- |
| Installation, crash, or incorrect file handling | Bug report |
| Unexpected or missing biological call | Scientific-result report |
| Slow run, memory pressure, or poor scaling | Performance report |
| New behavior or interface | Feature request |
| Security vulnerability | Private route in [SECURITY.md](SECURITY.md) |

The issue forms request the minimum useful evidence. A question about biological
interpretation is not automatically a software defect; provide comparison
coordinates and confidence categories so it can be investigated.

## Information to include

For most reports, include:

- software version, release-archive SHA-256, and lock-file SHA-256;
- exact command with secrets and private paths removed;
- operating system, CPU allocation, RAM, and storage type;
- input accession, sequence SHA-256, record count, and total bases;
- complete error text or the smallest relevant log tail;
- whether the problem reproduces with the bundled example;
- whether the output directory was new or reused.

For result disagreements, additionally state the coordinate system, overlap
rule, target categories, external-tool version, and whether different calls
may be representations of the same biological region.

For performance questions, follow the measurement checklist in
[CONTRIBUTING.md](CONTRIBUTING.md) and read
[docs/PERFORMANCE.md](docs/PERFORMANCE.md).

## Data privacy

Do not attach confidential, clinical, controlled-access, or personally
identifiable sequence data to an issue. Prefer a public accession. If a minimal
synthetic reproducer is possible, use it and state how it was constructed.
Checksums identify data without publishing it, but a checksum alone does not
allow maintainers to reproduce a private input.

## Scope of support

Maintainers can investigate installation, runtime, report integrity, and
reproducible scientific discrepancies. They cannot provide clinical advice,
guarantee that a predicted array is active in vivo, or treat the model score as
a calibrated biological probability.

See the [README](README.md), [cluster guide](docs/CLUSTER.md),
[performance notes](docs/PERFORMANCE.md), and
[output contract](docs/OUTPUTS.md) before opening an issue.
