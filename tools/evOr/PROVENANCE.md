# Source provenance

Imported on 2026-09-27 from `web_server_identify_placer_evor/CRISPR-evOr`.
The import includes the actual fixed working tree, including uncommitted changes
and new tests; cloning the upstream commit alone would lose those fixes.

- Upstream: https://github.com/fbaumdicker/SpacerPlacer
- Upstream base commit: `f2424aef2f79d6ee5b0f7c222abf33aa2d464aea`
- Imported files: 146
- Source inventory: [provenance/import.json](provenance/import.json)
- Original source SHA-256 hashes: [provenance/source-files.sha256](provenance/source-files.sha256)

Runtime Python modules, scientific data, examples, test fixtures, and tests are
copied byte-for-byte. Executable permissions on native helpers are preserved.
The original GPL license and authorship are retained. This is a new standalone
repository with a clean initial commit; the upstream history is referenced above.

`README.md` from the source is retained as `README_SPACERPLACER.md`; `README_CRISPR_EVOR.md` is unchanged. The package readme now points to the standalone `README.md`.
Only standalone documentation, packaging's readme reference (evOr), and ignore
rules were adjusted. Original Git metadata, build artifacts, installed package
metadata, caches, and bytecode are excluded. Upstream bundled scientific assets
are retained, including the likelihood files required by ODE mode.

The hash inventory describes the original source files, before these documented
README, packaging, and ignore-rule adjustments. It is a provenance record, not a
checksum list for the final repository root.

Version 1.0.3 adds a source audit after the import: exact alignment-cache/input identity checks, preserved dotted group
names, finite CLI parameter validation, supplied-tree validation, collision-free internal-node naming, actual fitted
tree provenance, and stricter selected-tree and rho-fallback output checks. The original inventory and validation JSON
remain the historical 1.0.2 import records. The current regression suite contains 124 passing tests; Ruff passes.
The obsolete, unreferenced `additional_data/230329_death_lh_up_to_68_lambdifyed_old.pickle` is removed; supported ODE
entrypoints continue to use both likelihood resources under `sp_model/bdm_likelihood_computation/sympy_bdm_lh_fct/`.
Its original hash remains in the historical import inventory.
