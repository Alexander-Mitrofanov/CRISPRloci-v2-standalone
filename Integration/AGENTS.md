# Integration agent guide

Integration adapts isolated scientific tools into versioned machine-readable
contracts. It is part of the scientific release identity.

Use `crispr_integration/README.md` as the module and hotspot map before opening
the large runner or parser.

- Keep runner commands operator-owned and shell-free.
- Validate each upstream tool output before emitting adapter or pipeline output.
- Preserve schema tolerance only where explicitly covered by fixtures.
- A contract change must propagate through backend validation, public summary,
  frontend models/views, fixtures, documentation, version, wheel, and release
  manifests.
- Do not place both SpacerPlacer and CRISPR-evOr `sp_model` trees on one import
  path. Tests run in isolated working directories/processes.
- Source movement or wheel changes are release-scoped; do not regenerate
  distributions during an unrelated refactor.
- `runner.py` deliberately retains the visible stage order and independent
  output gates; `parser.py` keeps canonical/legacy source-binding state in one
  place. Split either only as a versioned scientific release with regenerated
  identities and the full contract chain.

Focused check: `python -m pytest -q -p no:cacheprovider Integration/tests`.
