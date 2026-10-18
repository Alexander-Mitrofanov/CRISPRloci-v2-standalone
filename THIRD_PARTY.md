# Licensing and third-party source inventory

The root [MIT license](LICENSE) covers original CRISPRloci application code,
including the frontend, backend, worker, Integration adapters and repository
scripts. It does not relicense scientific tools, bundled third-party code,
models, fixtures or reference data. Their original copyright, license, NOTICE,
UPSTREAM and provenance files remain in place and govern those components.

This is an aggregation of separately licensed tools, invoked through isolated
processes. In particular, GPL components retain their GPL terms. Each module's
nested notices also cover its bundled dependencies; the table below is an entry
point, not a replacement for those notices.

| Module | Version | Source commit | License entry point |
| --- | --- | --- | --- |
| [CRISPRidentify-v2](tools/CRISPRidentify-v2/README.md) | 2.3.2 | `e8a3d7c038ea` | [MIT; bundled CasAndra notices retained](tools/CRISPRidentify-v2/LICENSE) |
| [CasAndra](tools/CasAndra/README.md) | 0.3.0.dev2 | `ba046f90242a` | [Apache-2.0](tools/CasAndra/LICENSE) |
| [CRISPRtracrRNA-v3](tools/CRISPRtracrRNA-v3/README.md) | 3.0.0 | `9307e2667733` | [MIT; native FASTA notices retained](tools/CRISPRtracrRNA-v3/LICENSE) |
| [SpacerPlacer](tools/SpacerPlacer/README.md) | 1.0.1.post3 | `2a91510a8ec9` | [GPL-3.0](tools/SpacerPlacer/LICENSE) |
| [evOr](tools/evOr/README.md) | 1.0.4 | `505dafe5bc11` | [GPL-3.0](tools/evOr/LICENSE) |
| [CRISPRleader-v2](tools/CRISPRleader-v2/README.md) | 2.0.0a1 | `b53101638356` | [MIT](tools/CRISPRleader-v2/LICENSE) |
| [CRISPRrepeat](tools/CRISPRrepeat/README.md) | 1.0.1 | `602155c32968` | [MIT](tools/CRISPRrepeat/LICENSE) |
| [CRISPRspacer](tools/CRISPRspacer/README.md) | 0.3.1 | `ddea3f5da9d5` | [MIT](tools/CRISPRspacer/LICENSE) |
| [CRISPRmap-v2](tools/CRISPRmap-v2/README.md) | 2.1.0a2 | `4b35029a3600` | [MIT](tools/CRISPRmap-v2/LICENSE) |
| [RepeatTyper](tools/RepeatTyper/README.md) | 2.0.1 | `d55ebf535f21` | [MIT](tools/RepeatTyper/LICENSE) |
| [SpacePHARER](tools/SpacePHARER/README.md) | 2.1.1 | `8dabb0d0f4d7` | [GPL-3.0-or-later; vendored dependency notices retained](tools/SpacePHARER/LICENCE.md) |
| [CCTK](tools/CCTK/README.md) | 2.0.2 | `c2f619fb3dde` | [GPL-3.0](tools/CCTK/LICENSE) |
| [clinker-v2](tools/clinker-v2/README.md) | 2.0.1 | `1526377ef54b` | [MIT; vendored dependency notices retained](tools/clinker-v2/LICENSE) |
| [ISEScan-v2](tools/ISEScan-v2/README.md) | 2.0.0 | `30877d2ca43f` | [Apache-2.0; supplied model notices retained](tools/ISEScan-v2/LICENSE) |
| [DefenseFinder-v2](tools/DefenseFinder-v2/README.md) | 2.0.1 | `f579941ca35f` | [GPL-3.0; dependency/model notices retained](tools/DefenseFinder-v2/LICENSE) |
| [RNAalifold-v2](tools/RNAalifold-v2/README.md) | 1.0.1 | `0498c888203b` | [MIT wrapper; separate ViennaRNA and MAFFT terms](tools/RNAalifold-v2/LICENSE) |
| [Phigaro-v2](tools/Phigaro-v2/README.md) | 2.4.1 | `0fae62b856a6` | [MIT; external executable notices retained](tools/Phigaro-v2/LICENSE.md) |

Additional notices include `frontend/public/licenses/igv.txt`,
`Integration/crispr_integration/licenses/`, and fixture-specific license files.
Scientific citation and upstream attribution files are retained within each
module. Cite the original methods as well as the CRISPRloci platform where used.

External dependencies installed by Conda/Mamba or pip retain their own terms.
ViennaRNA has a custom license; the RNAalifold wrapper's MIT license does not
relicense ViennaRNA. Reference datasets and databases retain their own access,
attribution and redistribution terms. No independent ISfinder database, complete
CRISPRcasDb catalogue, spacerome or production reference collection is included.
See [context integration licensing](docs/context-integrations.md).
