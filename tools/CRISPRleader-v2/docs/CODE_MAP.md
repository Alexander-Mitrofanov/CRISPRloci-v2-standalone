# Code map

Runtime package: `src/crisprleader/`. Read [CONTRACT.md](CONTRACT.md) and
[STATUS.md](STATUS.md) before changing scientific behavior. Prediction remains
unavailable; context extraction and descriptive evaluation are operational.

| Responsibility | Source | Main checks |
| --- | --- | --- |
| Frozen records and JSON snapshots | `models.py` | `test_models.py`, `test_schemas.py` |
| Bounded FASTA/native binding | `sequence_io.py` | `test_sequence_io.py` |
| Canonical Identify conversion | `adapters/identify.py` | `test_identify.py` |
| Circular arcs, flanks and projection | `coordinates.py` | `test_coordinates.py` |
| Immutable offline references | `references/__init__.py` | `test_references.py` |
| Exact discrete PWM tails and complete BH scope | `inference/motifs.py` | `test_motifs.py` |
| Extraction and independent result validation | `pipeline.py`, `validation.py` | `test_outputs.py` |
| Atomic exports and offline CLI | `export.py`, `cli.py` | `test_outputs.py`, `test_cli.py` |
| Bounded subprocess groups | `runtime/__init__.py` | `test_runtime.py` |
| Supplied-label descriptive evaluation | `evaluation.py` | `test_evaluation.py`, `test_benchmark_cli.py` |

Extraction flows from `load_inputs` through `extract_candidates`, `json_record`
and independent `validate_result`, then `write_bundle`. The independent validator
reconstructs source slices without importing the coordinate implementation.
Native segments remain source-forward; strand is applied after concatenation.
Candidate-local motif coordinates index the oriented flank plus observed repeat.

PWM tables depend on the matrix, background and score scale, never the candidate.
Each scan prepares matrices and backgrounds once per model, retrieves null tables
through the existing bounded cache, encodes each context once, and counts
ambiguous windows with a prefix sum. All eligible windows remain in the
BH family, including ambiguous windows assigned p=1. Reporting follows adjustment.
Score summation order, quantization, null distributions and scope receipts are
preserved; statistical occurrences do not assert biological function.

`json_record` copies nested record/list/dict containers in one traversal. Preserve
null values, full score precision and isolation between an export and its input.
Do not collapse clipped candidates or fill missing scientific values.

Run `scripts/check.sh` for tests, lint, types, packaging and installed-wheel CLI
checks. Numerical and export fixtures are synthetic software controls; the
biological gates in `benchmarks/GATES.md` remain separate and unevaluated.
