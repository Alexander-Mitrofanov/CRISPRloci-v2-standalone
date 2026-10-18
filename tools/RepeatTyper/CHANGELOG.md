# Changes

## 2.0.1

- Anchor bundled provenance before validating model and label bytes, load only
  checksum-verified bytes, and report missing/incompatible resources as model errors.
- Validate the multiclass objective and finite normalized score matrices before
  producing repeat evidence; preserved model predictions and policy are unchanged.
- Reject malformed API records, digests and arbitrarily large thresholds with
  typed input errors.
- Publish with an atomic operation that refuses concurrently created outputs,
  including empty directories; clean up failed temporary publications.

## 2.0.0

- Extract the standalone CRISPRCasTyper 1.9.0 repeat classifier into a typed src-layout package.
- Preserve exact canonical 4-mer/GC/length features, label order, model and iteration range.
- Replace implicit error exits with explicit nonzero status; publish atomic JSON/TSV outputs.
- Retain occurrence identity, deduplicate inference, expose unsupported bases and uncertainty.
- Pin and validate source/model identities; add parity, hostile-input, API and CLI tests and CI.
