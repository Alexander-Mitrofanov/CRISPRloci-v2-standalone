# Public-summary module map

`crispr_web.summary` is the stable façade. Implementation is split by the
kind of contract an agent changes:

- `contracts.py` — size limits, allowlists, and schema constants
- `sanitize.py` — recursive public-output privacy and size enforcement
- `bias.py` — bias-correction policy and warning normalization
- `bounded_io.py` — safe paths, bounded reads, JSON, and checkpoints
- `detection.py` — CRISPRidentify report summary
- `adapter.py` — group membership and eligibility summary
- `pipeline.py` — sanitized Integration stage summary
- `result_tables.py` — orientation and SpacerPlacer result tables
- `reconstruction_trees.py` — tree/source parsing
- `reconstruction_events.py` — canonical spacer and event identities
- `reconstruction_discovery.py` — bounded reconstruction discovery/publication
- `builder.py` — final schema composition

`common.py`, `reconstruction.py`, and `stages.py` are narrow compatibility
import façades. New implementation code belongs in the responsibility module,
not a façade.

Run:

```bash
cd backend
.venv/bin/pytest -q tests/test_output_hardening.py \
  tests/test_resource_exhaustion_regressions.py
```

