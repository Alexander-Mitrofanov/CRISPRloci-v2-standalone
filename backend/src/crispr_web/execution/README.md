# Worker execution module map

`crispr_web.worker` remains the stable CLI/test façade. It owns the worker
loop and process-containment sequence because failure to kill an attributed
process must stop future claims in the same orchestration context.

Focused helpers live here:

- `contracts.py` — mode identifiers and successful pipeline statuses
- `files.py` — bounded output reads and log capture
- `storage.py` — storage accounting and fail-closed attempt quarantine
- `input_integrity.py` — canonical FASTA binding verification
- `provenance.py` — sanitized public provenance
- `validation_common.py` — schema, number, and result-table primitives
- `validation_newick.py` — bounded generated-tree grammar/topology checks
- `validation_adapter.py` — adapter membership and orientation records
- `validation_identify.py` — CRISPRidentify report/input binding
- `validation_manifest.py` — adapter counts and exact group membership
- `validation_reconstruction.py` — SpacerPlacer tables and generated trees
- `validation_orientation.py` — evOr metadata, likelihoods, and tree policy
- `validation_association.py` — SpacePHARER groups, calibration and source-bound alignments
- `validation_association_reference.py` — independent immutable catalogue/index verification and bounded target reads
- `validation_pipeline.py` — small end-to-end status/composition gate
- `validation.py` — compatibility imports used by `worker.py` and tests

Backend validation intentionally does not import Integration validators.
Changes here require worker, hardening, Integration-contract, output, and
resource-exhaustion tests.
