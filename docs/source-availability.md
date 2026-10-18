# Complete source availability and reproduction

The publication repository is
[Alexander-Mitrofanov/CRISPRloci-v2](https://github.com/Alexander-Mitrofanov/CRISPRloci-v2).
It contains the application and all 17 scientific source modules. An ordinary
clone or GitHub source download includes the module files directly: there are
no submodules, Git LFS pointers or private-origin downloads required to inspect
or restore the scientific source.

## Source inventory

- `frontend/`: Vue user interface and frontend tests.
- `backend/`: FastAPI API, durable SQLite queue, isolated worker, independent
  scientific output validation and protected result exports.
- `Integration/`: scientific command orchestration, adapters and validation.
- `tools/`: all 17 scientific modules, including their native source, tests,
  bundled models, licenses and upstream provenance.
- `scripts/`: setup, dependency provisioning, immutable runtime releases,
  scientific acceptance gates, database preparation and source/privacy checks.
- `tests/`, `docs/`, `deploy/`: repository regression tests, contracts, sample
  input/output and generic deployment templates.

The [root README](../README.md) maps the 20 analysis modes to their scientific
modules. The [tool inventory](../tools/README.md) lists every exact source pin.
The [license inventory](../THIRD_PARTY.md) explains the separately licensed
components and third-party dependencies.

## Verify a downloaded source snapshot

```bash
python3 scripts/sync_tools.py --check
python3 scripts/check_private_data.py workspace
```

`tools/tools.lock.json` preserves the original upstream repository identities,
versions and commits. In the complete publication repository,
`tools/bundled-sources.json` additionally records their original Git tree IDs and the exported source tree
IDs. Git applies each module's own `.gitattributes` line-ending rules when
exporting sources; these rules and their source files remain included.
Verification recomputes each complete exported tree, including executable permissions,
and rejects missing, changed or additional source files. It works offline without
nested `.git` directories. Python/test caches do not change source identity.
Runtime provisioning verifies these same trees before copying sources into an
isolated immutable toolset. The development layout without a bundled manifest
continues to use separate pinned Git checkouts.

## Install and run

Use Linux x86-64, Python 3.12, Node 24, Git, Conda/Mamba and a C++ toolchain.
Follow [Run locally](../README.md#run-locally). Network access is needed for
public package/native dependencies and explicitly acquired reference datasets;
access to the private module origins is unnecessary for the bundled snapshot.
SpacerPlacer and evOr must retain separate environments because both own
`sp_model`. API processes never import or execute scientific tools.

After setup, run `bash scripts/check.sh` for application regression tests and
source/build checks. Scientific runtime provisioning and the real-science
release gate are separate, as documented in the README. The tests and documented
[precomputed examples](precomputed-examples.md) provide reproducible inputs and
outputs; the [verification guide](verification.md) states their coverage.

This is a source distribution. It excludes credentials, operator configuration,
connection records, live jobs, runtime environments, local deployment evidence
and production databases. Generic configuration examples remain in `deploy/`.
Large external reference collections are acquired explicitly using the documented
import/preparation scripts, with their own source and checksum identities.
No production deployment is performed by creating or cloning this repository.

## Manuscript archive and availability statement

The repository provides one complete source link for editors. A GitHub URL alone
does not provide a permanent DOI. NAR's current
[data deposition instructions](https://academic.oup.com/nar/pages/data_deposition_and_standardization?login=false)
require a persistent archive and describe Figshare import for code hosted on GitHub.

For submission, freeze the manuscript's exact source version as a GitHub release,
import that release into Figshare, verify the full source archive, publish the
archive and put the assigned DOI plus the GitHub URL in the manuscript's Data
Availability statement. Complete author/title/version/license metadata for that
specific manuscript version. No archive DOI or accepted manuscript version is
claimed by this source repository.

Availability wording once a DOI has been issued:

> The complete CRISPRloci v2 source code, including its web interface, API,
> isolated worker, integration adapters and 17 pinned scientific modules, is
> available at https://github.com/Alexander-Mitrofanov/CRISPRloci-v2. The exact
> version used for this study is archived in Figshare [insert assigned DOI].
> Individual scientific tools retain their original licenses and attribution.

Replace the DOI placeholder only with an actual published archive identifier.
