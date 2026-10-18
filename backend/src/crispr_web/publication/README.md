# Artifact publication module map

`crispr_web.artifacts` is the stable compatibility façade.

- `contracts.py` owns allowlists, denied path parts, and size limits.
- `safety.py` owns private-path detection, manifest sanitization, safe file
  resolution, hashing, registration, and fsync.
- `discovery.py` owns bounded output traversal and public text/JSON screening.
- `archive.py` owns the exact streaming ZIP byte cap.
- `collector.py` composes provenance, safe downloads, and the result bundle.
- `array_data.py` reduces full array reports to scientific data and verifies
  observed sequences against original source intervals.
- `sequence_ports.py` publishes bounded, complete module handoffs using the
  versioned models in `crispr_web.exchange`.
- `module_result.py` publishes the common versioned `result.json` envelope for
  every successful mode. `compact/` projects full validated evidence into
  module-specific data, deduplicating folds, references and source metadata.
  It never uses capped browser previews or substitutes empty findings when
  an export is unavailable.

Keep these checks independent from summary sanitization. Changes require output
hardening, resource exhaustion, worker, and privacy tests.
