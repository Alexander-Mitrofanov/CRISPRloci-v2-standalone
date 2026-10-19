# Scientific tools

The complete [publication repository](https://github.com/Alexander-Mitrofanov/CRISPRloci-v2-standalone)
bundles all 17 module source trees directly. Their original GitHub repositories
are independent development origins; [tools.lock.json](tools.lock.json) records
exact commits. In a publication snapshot, `bundled-sources.json` records the
original and exported Git tree IDs and allows offline verification without private access.
Every source tree retains its original license and upstream provenance; see
[the license inventory](../THIRD_PARTY.md).

Standalone bug fixes retain the original upstream commit and `upstream_tree` as
provenance. Patched entries also pin `bundled_tree` in `tools.lock.json`, matching
the exported `tree` in `bundled-sources.json`. This pin participates in the
runtime identity, so patched sources build a new isolated toolset instead of
reusing a runtime built from older bundled code.

The [19 October standalone audit](../docs/standalone-tool-audit-2026-10-19.md)
lists per-tool fixes, commits and verification limits.

| Tool | Version | Pinned commit | Command |
| --- | --- | --- | --- |
| [CRISPRidentify-v2](CRISPRidentify-v2/README.md) | 2.3.2 | `e8a3d7c038ea` | `./crispridentify` |
| [CasAndra](CasAndra/README.md) | 0.3.0.dev2 | `ba046f90242a` | `casandra` |
| [CRISPRtracrRNA-v3](CRISPRtracrRNA-v3/README.md) | 3.0.0 | `9307e2667733` | `crisprtracrrna` |
| [SpacerPlacer](SpacerPlacer/README.md) | 1.0.1.post3 | `2a91510a8ec9` | `spacerplacer` |
| [evOr](evOr/README.md) | 1.0.4 | `505dafe5bc11` | `crispr-evor` |
| [CRISPRleader-v2](CRISPRleader-v2/README.md) | 2.0.0a1 | `b53101638356` | `crisprleader` |
| [CRISPRrepeat](CRISPRrepeat/README.md) | 1.0.1 | `602155c32968` | `crisprrepeat` |
| [CRISPRspacer](CRISPRspacer/README.md) | 0.3.1 | `ddea3f5da9d5` | `crisprspacer` |
| [CRISPRmap-v2](CRISPRmap-v2/README.md) | 2.1.0a2 | `4b35029a3600` | `crisprmap` |
| [RepeatTyper](RepeatTyper/README.md) | 2.0.1 | `d55ebf535f21` | `repeattyper` |
| [SpacePHARER](SpacePHARER/README.md) | 2.1.1 | `8dabb0d0f4d7` | `spacepharer` |
| [CCTK](CCTK/README.md) | 2.0.2 | `c2f619fb3dde` | `cctk` |
| [clinker-v2](clinker-v2/README.md) | 2.0.1 | `1526377ef54b` | `clinker-v2` |
| [ISEScan-v2](ISEScan-v2/README.md) | 2.0.0 | `30877d2ca43f` | `isescan-v2` |
| [DefenseFinder-v2](DefenseFinder-v2/README.md) | 2.0.1 | `f579941ca35f` | `defense-context` |
| [RNAalifold-v2](RNAalifold-v2/README.md) | 1.0.1 | `0498c888203b` | `rnaalifold-v2` |
| [Phigaro-v2](Phigaro-v2/README.md) | 2.4.1 | `0fae62b856a6` | `phigaro-v2` |

The five context additions have separate license/model decisions and bounded
workflow contracts in [context integrations](../docs/context-integrations.md).

## Restore or verify

From the CRISPRloci-v2 project root, with Python 3 and Git available:

```bash
# Only for the separate development layout without bundled sources:
# gh auth login
# gh auth setup-git

# Clone missing tools or check out the exact commits in the lock file.
python3 scripts/sync_tools.py

# Offline verification: origins, commits, and clean worktrees.
python3 scripts/sync_tools.py --check

# Optional: operate on one tool only.
python3 scripts/sync_tools.py --check --tool SpacerPlacer
```

In a publication snapshot, both commands verify the bundled source trees
without contacting any origin. In the separate development layout, the restore
script refuses to overwrite local changes or reuse a checkout with a different origin. It never pulls a moving branch or automatically advances the
lock file. Git authentication is required only for cloning/fetching.

The development workspace uses independent checkouts rather than submodules;
its `.gitignore` excludes those nested repositories. The publication repository
tracks exported source trees instead, so ordinary GitHub source archives contain
every module. Runtime source export supports both layouts and preserves the
original lock identities.

## Tool responsibilities

- **CRISPRidentify v2** detects CRISPR arrays and produces canonical `report.json`
  output. It also bundles an optional CasAndra runtime.
- **CasAndra** detects Cas proteins and classifies nearby Cas cassettes using the
  bundled cap15 model. Its standalone checkout lets the server run Cas analysis
  independently of CRISPR array detection.
- **CRISPRtracrRNA v3** predicts tracrRNA candidates. Model mode uses Infernal;
  complete mode also consumes independently installed CRISPRidentify and CasAndra.
- **SpacerPlacer** reconstructs ancestral CRISPR spacer arrays and estimates
  evolutionary parameters for related array groups.
- **evOr** compares forward and reverse evolutionary reconstructions to report
  `Forward`, `Reverse`, or `ND` acquisition orientation.

The identify pin includes the 6 October 2026 reconstruction and memory update;
its upstream `VERSION` remains 2.3.2. The exact commit distinguishes it from the
published 2.3.2 archive. The matching tracrRNA dependency lock retains that identify pin and the
CasAndra pin; its October 13 refinements preserve their compatibility.

CRISPRtracrRNA's [upstream lock](CRISPRtracrRNA-v3/upstream-lock.json) references
exactly the CRISPRidentify and CasAndra commits pinned here.

Use separate processes and explicit executable paths when combining tools.
**SpacerPlacer and evOr both install `sp_model`; never install them in the same
environment or import them together.** Avoid adding all tool directories to one
`PYTHONPATH`. The worker gives each job its own output directory and preserves scientific
warnings, missing values, and uncertainty from tool outputs.

CRISPRleader v2 extracts both adjacent contexts for accepted arrays. Its current
release does **not** predict leaders. The server uses `adapt identify`, `extract`
and `validate` in an isolated environment. See [server integration](../docs/leader-integration.md).

## Runtime installation

Run `python3.12 scripts/provision_tools.py` from the project root. It exports
exact Git commits, creates separate environments, verifies command entry points,
and records source, installed-file and dependency receipts. See the
[root README](../README.md) for application builds, real scientific gates,
activation and updates. IntaRNA uses its own native environment because its
ViennaRNA dependency is incompatible with the Python 3.12 RNA orchestrator.

## Extracted source provenance

SpacerPlacer and evOr were imported from the fixed working trees in
`/home/alex/CursorProjects/web_server_identify_placer_evor`, including uncommitted
fixes and new regression tests. Each new repository has a clean initial commit,
retains the upstream GPL license, and records the upstream commit, original file
hashes, and package verification results:

- [SpacerPlacer provenance](SpacerPlacer/PROVENANCE.md) and
  [validation](SpacerPlacer/provenance/validation.json).
- [evOr provenance](evOr/PROVENANCE.md) and
  [validation](evOr/provenance/validation.json).

Ignored caches, bytecode, build output, and installed package metadata were
excluded. Original scientific runtime files, data, examples, and tests were
preserved byte-for-byte. The original source directories were not modified.

CRISPRspacer checks uploaded target DNA against an operator-installed immutable
spacer snapshot. Full reference data stay on the VM; only a three-record sample
is bundled. See [protospacer integration](../docs/protospacer-integration.md).

CRISPRmap-v2 is the independent 2026 prokaryotic repeat reference mapper. It
provides offline exact historical annotation recovery and near-reference family
candidates with uncertainty. Its [council record](CRISPRmap-v2/docs/COUNCIL.md)
and [original-output benchmark](CRISPRmap-v2/benchmarks/results/REPORT.md)
describe the scope and measured limitations. Its historical `classify` contract is exposed as queued `repeat_map`; modern catalogue/index/group-search commands remain standalone. Operator provisioning gives the web mapper a separate standard-library-only Python environment. The original 2013 reference is acquired explicitly with its
SHA-pinned importer and is not bundled or automatically fetched during jobs.

CRISPRmap has a queued historical-reference web workflow; its newer standalone modern catalogue contract is not yet exposed. See the [integration audit](../docs/integration-audit-2026-10-03.md).

RepeatTyper is integrated as the queued `repeat_type` mode, consuming repeat DNA
FASTA or the existing repeat occurrence handoff. See [repeat typing](../docs/repeat-typing.md) for immutable model identity,
uncertainty semantics, the isolated runtime and validated exports.

SpacePHARER consumes explicit spacer groups in queued `spacer_association` jobs.
Its isolated native engine uses operator-provisioned target/control panels and
reports unavailable empirical calibration explicitly. See
[spacer association](../docs/spacer-association.md) for panel limits and preparation.

CCTK consumes detection `arrays.json` in queued `array_compare` jobs. It retains
ordered occurrences and deletion slots while computing observed sharing networks;
its standalone evolutionary commands remain separate from this web workflow.
See [array comparison](../docs/array-comparison.md) for contracts and web limits.

All three repositories were created, modernized and published before their server
integration pins were added. The seventeen-tool runtime recipes cover 20 modes. The full release gate defines
37 native cases and 5 predispatch input rejections. Current refinement and
verification status are recorded in [the October 13 refinement](../docs/tool-refinement-2026-10-13.md);
the earlier twelve-tool release remains documented in [the October 11 review](../docs/tool-review-2026-10-11.md); the earlier modernization snapshot remains in [the implementation record](../improvenments/08-tool-modernization-progress.md).
