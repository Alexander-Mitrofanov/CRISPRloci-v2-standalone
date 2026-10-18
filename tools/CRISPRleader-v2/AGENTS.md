# CRisprleader v2

The owner approved the council plan on 27 September 2026. Work only in this
independent repository/its assigned worktree. Never edit or activate another
scientific tool or the CRISPRloci server as part of standalone development.

- Python package `src/crisprleader`; source-forward zero-based half-open native coordinates.
- Reference-supported conserved cores, extracted windows, motifs and full
  biological leaders are distinct. Never manufacture confidence or absence.
- Inference is offline; references are immutable data, never executable pickle.
- Source hashes, coordinates, sequences and exports must independently agree.
- Read `docs/CONTRACT.md` before implementation. Propose contract changes to the
  coordinator; do not silently change shared domain types.
- Workers use isolated branches/worktrees and exclusive ownership. W1 owns
  coordinates/sequence_io/adapters; W2 owns inference/references; W3 owns
  runtime/benchmarks; coordinator owns shared models, CLI, validation, export,
  packaging, schemas and documentation.
- Include meaningful tests with each scientific/runtime change. Biological
  benchmarks and software tests answer different questions. No training/test leakage.
- A default biological reference bundle and measured independent evaluation
  are required; synthetic models alone do not complete the tool.
- Do not copy third-party code, models or tables with undocumented reuse terms.
  Independently authored methods and provenance-recorded source sequences are
  the approved route; preserve source-specific terms and experimental scope.
- Run `python -m pytest`, `ruff check`, type checks and the documented release
  validation applicable to changes. Keep errors actionable and logs bounded.
- No production services, public release, package registry publication or server
  runtime changes without an explicit separate request.
