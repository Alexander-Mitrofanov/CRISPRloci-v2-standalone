# Code map

Paths below are relative to the repository root. The `tool/` checkout contains
the application; surrounding benchmark and research work is separate.

## Execution and data flow

`crispridentify` → `scripts/run_crispridentify.sh` → `CRISPRidentify.py` →
`components/pipeline.py:Pipeline`.

The launcher activates the locked environment. The Python entry point sets
hash determinism, resolves CLI inputs and CPU limits, and dispatches files or
folder lanes. `components/cli/parser.py:build_parser` only builds the parser;
`presets.py` resolves scientific defaults. `Pipeline.__init__` coordinates
detection, refinement, scoring, enhancement/rescue, reconstruction, annotation,
and report publication. Individual features are controlled by the resolved
configuration; preserve the stage order when changing orchestration.

`crispridentify-batch` launches `scripts/run_crispridentify_batch.py`. This
stdlib-only supervisor owns cohort manifests, bounded workers, immutable
attempts, completion validation, and resume. Its main navigation points are
`index_manifest`, `build_tool_configuration`, `validate_report`, `execute_job`,
`run_bounded`, and `command_run`.

## Find an implementation and its tests

Test names below are modules under `tests/`; run them with
`.envs/crispridentify-core/bin/python -m unittest tests.<name> -v`.

| Concern | Implementation under `components/` | Tests to start with |
| --- | --- | --- |
| CLI options and presets | `cli/parser.py`, `cli/validation.py`, `cli/presets.py`; runtime plumbing in `CRISPRidentify.py` | `test_runtime_presets`, `test_v2_casandra_timeout_plumbing`, `test_cli_worker_import` |
| CPU lanes and subprocess cleanup | `runtime/folder_scheduler.py`, `runtime/runtime_resources.py`, `runtime/bounded_subprocess.py`, `runtime/external_tools.py` | `test_runtime_resources`, `test_v2_folder_scheduler_casandra`, `test_v2_bounded_subprocess` |
| FASTA and source correctness | `sequence/fasta.py`, `sequence/coordinates.py`, `sequence/validation.py`, `sequence/topology.py` | `test_reconstruction_validation`, `test_contextual_spacer_circular` |
| Discovery orchestration | `detection/stage.py:Detection`, `detection/search.py` | `test_detection_streaming`, `test_v2_fuzzy_boundary_paths_flag` |
| Executed discovery configuration | `detection/search_contract.py` | `test_search_contract` |
| Supplied repeat panel, bounded scanning and source projection | `detection/repeat_templates.py`, `detection/template_discovery.py` | `test_repeat_template_library`, `test_template_discovery`, `test_template_discovery_integration` |
| Cross-route duplicate competition and observational ancestry | `refinement/template_integration.py`, `evidence/discovery_origins.py` | `test_template_discovery`, `test_template_discovery_integration`, `test_candidate_origins` |
| Native matching and occurrence paths | `detection/literal_fuzzy.py`, `detection/fuzzy_*.py`, `detection/vmatch_occurrence_chain.py`; C implementation in `native/` | `test_fuzzy_lattice_storage`, `test_detection_streaming`, `test_v2_distribution` |
| Candidate representation and structural filtering | `refinement/candidates.py:CrisprCandidate`, `refinement/stage.py:DetectionRefinement` | `test_refinement_consumption`, `test_candidate_origins` |
| Enhancement and source alignment | `refinement/enhancement.py:EvaluatedArraysEnhancement`, `refinement/enhancement_search.py` | `test_source_alignment_cache`, `test_source_repeat_edit_distance`, `test_source_weighted_medoid`, `test_source_unit_edit_kernel` |
| Local spacer-period reconstruction and curated boundaries | `refinement/periodic_completion.py`, `refinement/periodic_transaction.py`, `refinement/boundary_panel.py` | `test_periodic_completion` |
| Corroborated fragment recovery and optional source-derived extent proposals | `refinement/extent_trials.py` | `test_extent_trials` |
| Convergent longer candidates and competing extents | `refinement/extent_trials.py` | `test_extent_trials` |
| Comparison of independently supported boundaries, memberships and terminal placements | `refinement/comparative_transaction.py`, `refinement/reviewed_boundary_selection.py`, `refinement/reviewed_membership_selection.py`, `refinement/reviewed_terminal_selection.py` | `test_comparative_transaction`, `test_reviewed_boundary_selection`, `test_reviewed_membership_selection`, `test_reviewed_terminal_selection` |
| Research comparisons of repeat and spacer boundaries | `evidence/review_hypotheses.py`, `reporting/candidate_evidence.py` | `test_review_hypotheses`, `test_candidate_evidence` |
| Divergent repeat membership | `refinement/potential_degenerate.py` | `test_potential_degenerate` |
| Terminal and internal repeat completion | `refinement/repeat_boundary_repair.py`, `refinement/internal_*.py`, `refinement/terminal_repair_completion.py` | `test_repeat_boundary_repair`, `test_internal_repeat_completion`, `test_internal_local_core_bridge`, `test_terminal_repair_completion` |
| Feature/model inference | `scoring/feature_schema.py`, `scoring/compact_feature_extraction.py`, `scoring/compact_model.py`, `scoring/compact_ml_runtime.py` | `test_scale_fast_path`, `test_machine_phase_score_gate` |
| EDeN-compatible features | `scoring/sequence_path_vectorizer.py`, `scoring/legacy_hash.py` | `test_scale_fast_path` |
| Evaluation and score reuse | `scoring/evaluation.py:ArrayEvaluation`, `scoring/cache.py:CandidateScoreCache` | `test_scale_fast_path`, `test_source_phase_score_cache` |
| Candidate lineage, initial competition and Low-score output selection | `evidence/candidate_origins.py`, `evidence/initial_competition.py`, `evidence/occurrence_protection.py` | `test_candidate_origins`, `test_initial_competition`, `test_occurrence_protection`, `test_low_score_selection` |
| Boundary/locus evidence | `evidence/boundary_assessment.py`, `evidence/biological_loci.py` | `test_reconstruction_boundary_family`, `test_extent_family_exports` |
| Bona-fide reconstruction | `reconstruction/reconstruction_selection.py`, `reconstruction/reconstruction_validation.py`, `reconstruction/reconstruction_transaction.py` | `test_reconstruction_selection`, `test_reconstruction_validation`, `test_reconstruction_transaction` |
| Phase references and tied boundary families | `reconstruction/reconstruction_reference.py`, `reconstruction/reconstruction_boundary_family.py`, `reconstruction/family_terminal_extension.py` | `test_reconstruction_phase_reference`, `test_reconstruction_alias_retention`, `test_reconstruction_boundary_family`, `test_family_terminal_extension` |
| Possible reconstruction | `reconstruction/possible_reconstruction.py` | `test_possible_reconstruction`, `test_possible_reconstruction_exports` |
| Identical-spacer context and competing extents | `reconstruction/localized_spacer_context.py`, `reconstruction/contextual_spacer_rescue.py`, `reconstruction/competing_terminal_extents.py` | `test_contextual_spacer_rescue`, `test_contextual_spacer_exports`, `test_competing_terminal_extents`, `test_competing_terminal_extent_exports` |
| Fast rescue | `rescue/fast_run_rescue.py`, `rescue/fast_two_repeat_terminal_repair.py`, `rescue/fast_accepted_family_rescue.py` | `test_fast_two_repeat_terminal_repair`, `test_fast_accepted_family_rescue`, `test_scale_fast_path` |
| Architecture rescue and source recovery | `rescue/architecture_rescue_*.py`, `rescue/source_*.py`, `rescue/shared_repeat_augmentation.py` | `test_source_phase_score_cache`, `test_source_alignment_cache`, `test_scale_fast_path` |
| Accepted-boundary completion and tandem checks | `rescue/accepted_core_boundary_*.py`, `rescue/accepted_family_tandem_*.py`, `rescue/accepted_source_intact_conflict.py`, `rescue/likely_tandem*.py` | `test_likely_tandem`, `test_fast_accepted_family_rescue` |
| Circular edge joining | `rescue/circular_array_join.py`, `rescue/circular_array_join_integration.py` | `test_contextual_spacer_circular`, `test_contextual_spacer_exports` |
| Cas, strand, leader and IS annotation | `annotation/stage.py:NonArrayComputations`, `annotation/casandra_adapter.py`, `annotation/sequence_annotations.py`, `annotation/portable_strand.py` | `test_v2_casandra_adapter`, `test_v2_casandra_classification`, `test_v2_casandra_timeout_plumbing` |
| Canonical report construction and replay | `reporting/stage.py:OutputMaker`, `reporting/report_model.py`, `reporting/report_extensions.py`, `reporting/report_provenance_validation.py` | `test_reconstruction_transaction`, `test_contextual_spacer_exports`, `test_extent_family_exports` |
| HTML, CSV, FASTA, and legacy exports | `reporting/html_report.py`, `reporting/legacy_exports.py`, `reporting/cas_output.py`, `reporting/cas_html.py`, `reporting/*_outputs.py` | `test_report_streaming`, `test_v2_casandra_reporting`, `test_competing_terminal_extent_exports` |
| Evidence retention and atomic publication | `reporting/candidate_evidence.py`, `reporting/output_transaction.py` | `test_candidate_evidence`, `test_v2_output_transaction` |
| Object lifetimes and stage sequencing | `pipeline.py` | `test_pipeline_lifetimes`, `test_refinement_consumption` |
| Cohort manifests and resume | `../scripts/run_crispridentify_batch.py` | `test_batch_runner` |
| Package/import compatibility and asset paths | `_compat.py`, `paths.py` | `test_package_layout`, `test_cli_worker_import`, `test_v2_distribution` |

Tests often cover more than one package because provenance must survive the
whole pipeline. The table is a starting point, not a substitute for the complete
suite before delivery. Fixtures live in `tests/data/` and `tests/fixtures/`.

## Naming and imports

`stage.py` is the public stage coordinator within detection, refinement,
annotation, and reporting. Helpers live beside it. Scientific feature names
remain explicit so, for example, a search for `internal_repeat` finds its
proposal, transaction, and tests directly. A `*_transaction.py` module applies
validated decisions; proposal and validation code remain separately named.

Use canonical imports in new code:

```python
from components.detection.stage import Detection
from components.refinement.candidates import CrisprCandidate
from components.scoring.cache import CandidateScoreCache
```

All 93 historical flat module names remain lazy import aliases in
`components/_compat.py`. The aliases share the canonical module object, so old
pickle class references and patches do not create duplicate class/module state.
The mapping also provides an exact old-to-new filename lookup. The old source
file paths have moved; code loading source files directly must use the new
paths. `components.pipeline` and its `CandidateScoreCache` export remain valid.

Package initializers intentionally avoid eager stage imports. Some scientific
modules have local imports to break dependency cycles; preserve their timing.
CLI workers must remain importable before scientific dependencies are loaded.

## Assets and contracts

- `components/paths.py` is the single anchor for repository assets, regardless
  of the caller's working directory or a module's package depth.
- `components/native/` holds the C kernels; `scripts/build_native_matcher.sh`
  builds the ignored shared library in that directory.
- `trained_models/` holds first-party inference artifacts; `tools/` holds BLAST
  and insertion-sequence reference data. These are data, not application code.
- `vendor/casandra/` is a pinned third-party distribution with its own manifest.
- `schemas/report-1.2.0.schema.json` and [OUTPUTS.md](OUTPUTS.md) define report
  structure and interpretation. [CLUSTER.md](CLUSTER.md) covers cohort execution.
- `.github/workflows/ci.yml` contains release/install checks.

Batch resume binds the runtime tree to the completion record. A source-layout
change therefore needs a fresh run directory, even when scientific results
are unchanged; do not bypass that integrity check.
