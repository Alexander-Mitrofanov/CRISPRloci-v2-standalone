# evOr maintenance

Read [README.md](README.md) and [docs/CODE_MAP.md](docs/CODE_MAP.md) before editing.

- This independent repository installs `sp_model`. Keep its environment and imports separate from SpacerPlacer.
- Preserve array/group order, spacer identities, exact likelihood arithmetic, confidence boundaries and explicit uncertainty.
- Estimate forward and reverse trees separately unless the caller supplies an independent tree. Export the actual fitted trees.
- Keep effective alpha/rho correction provenance and per-group fallback warnings together with the selected fit.
- Pickle input remains disabled unless a trusted standalone caller explicitly passes `--allow_unsafe_pickle`.
- Preserve public imports and legacy output filenames. Historical provenance inventories describe the original import.
- Validate changes with `python -m pytest -q` and `ruff check .`; scientific changes also need bundled-example or frozen output comparisons.
- Keep generated outputs, native assets, model resources and production deployment outside routine refactors.
