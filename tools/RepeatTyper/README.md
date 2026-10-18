# RepeatTyper 2

Offline repeat-sequence CRISPR-Cas subtype inference, extracted and refactored
from [CRISPRCasTyper 1.9.0](https://github.com/Russel88/CRISPRCasTyper/tree/v1.9.0).
This independently maintained MIT package preserves the original classifier and
37-label vocabulary. It does not perform array detection, retrain the model,
or replace protein-based system classification.

```sh
python -m venv .venv
.venv/bin/pip install -r requirements.lock .
.venv/bin/repeattyper classify --input examples/repeats.fasta --output results
```

Supply FASTA with unique ASCII identifiers and 1–200 IUPAC DNA bases per record.
Repeated sequences retain their distinct occurrence IDs and are classified once.
Ambiguous IUPAC sequences remain explicit `unsupported` rows with null predictions
and scores. Empty FASTA produces a completed empty result. Limits: 10,000 records,
2 MiB input, 32 MiB JSON. Output directory must not exist; complete JSON and TSV
are published together by atomic directory rename. Existing output is never overwritten.
No network, external commands or model downloads occur during analysis.
Atomic publication without replacement is supported on Linux, macOS and Windows;
other platforms fail explicitly if that operation is unavailable.

`report.json` uses `repeattyper.classification` schema 1.0.0. `predictions.tsv`
contains one row per occurrence. `winner` and `native_score` retain the original
closed-set result even below `--threshold` (default 0.75); `prediction` is null
below threshold. Scores are not calibrated confidence that a biological locus
belongs to a subtype. Unknown/out-of-vocabulary families can receive high scores.
Canonical reverse-complement 4-mers provide no orientation evidence.

Python API: `from repeattyper import Record, classify`; call
`classify([Record("copy1", "ACGTACGT")], input_sha256=<FASTA digest>)`.
The return value matches the JSON schema and contains immutable model provenance.
CLI exits: 0 completed, 2 invalid input/options, 3 incompatible/changed model,
4 filesystem/publication failure. Reports never fabricate missing scientific values.

The refactor separates bounded parsing, typed occurrence records, feature construction,
model validation, inference, atomic publication and CLI. Only numpy and CPU XGBoost
are runtime dependencies. Provenance, modern JSON model bytes and labels are SHA-256
verified before each run, and inference loads the same verified bytes. See
`src/repeattyper/data/provenance.json` for exact upstream commit,
wheel and model checksums. `requirements.lock` records the tested CPU runtime.

## Source navigation

The stable flow is `cli.py` → `input.py` → `inference.py` → `output.py` under
`src/repeattyper/`. `features.py` fixes the 138 feature columns and precomputes
all 256 canonical 4-mer mappings; `model.py` verifies package resources and loads
those same bytes for each invocation. Inference computes features and sequence
digests once per distinct sequence, then restores occurrence IDs and order.
FASTA record and sequence limits are enforced while parsing as well as in the API.

`tests/upstream/xgb.py` is the frozen scientific oracle. Keep its implementation,
model resources, vocabulary, float32 features, zero sentinel row, tie selection
and exclusive iteration range unchanged. `test_classifier.py` compares the full
37-column score matrix, winners and scores on frozen/random/reverse-complement
controls; `test_refinement.py` independently checks every 4-mer and boundary policy.
`schemas/classification-1.0.0.json` defines the public result contract.

## Scientific preservation and validation

The upstream exclusive iteration range `(0, 41)` is intentionally preserved, even
though `best_ntree_limit` is 42; changing this would change predictions. Tests execute
the original upstream `XGB.predict_repeats` implementation against the bundled model,
comparing every score and label on frozen upstream repeat examples and their reverse
complements. Those examples are parity fixtures drawn from the upstream data, **not**
independent biological validation or a performance benchmark. Generalization and
subtype-specific accuracy remain subject to the original model's limitations.

```sh
pip install -e '.[test]'
pytest -q
ruff check src tests --exclude tests/upstream
mypy
python -m build
```

## Attribution and citation

Original CRISPRCasTyper software and model: Jakob Russel and contributors.
The original MIT license is preserved in `LICENSE`; the parity oracle is the
unmodified upstream `cctyper/xgb.py`, also MIT licensed.
Please cite Russel et al., *CRISPRCasTyper: Automated Identification, Annotation,
and Classification of CRISPR-Cas Loci*, The CRISPR Journal (2020),
[doi:10.1089/crispr.2020.0059](https://doi.org/10.1089/crispr.2020.0059).
This fork is not an endorsement by upstream authors.
