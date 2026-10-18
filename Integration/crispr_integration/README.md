# Integration module map

This package is the versioned, independently tested boundary around the pinned
scientific tools.

- `contract.py` owns schema identifiers and contract-level constants.
- `models.py` owns adapter records and JSON serialization shapes.
- `parser.py` owns bounded CRISPRidentify report discovery plus canonical and
  legacy source binding. Start at `load_identified_arrays`; the canonical and
  legacy branches are isolated below it.
- `adapter.py` owns array grouping and the downstream FASTA/manifest boundary.
- `runner.py` owns shell-free stage sequencing. Start at `run_pipeline`; output
  gates are grouped as SpacerPlacer, CRISPR-evOr, bias policy, and Newick/leaf
  identity validators.
- `cli.py` is the narrow command-line façade.
- `protospacer.py` owns the optional CRISPRspacer database-search stage and its
  source/reference/alignment gate. It bypasses array detection and comparative
  grouping; it accepts an operator path and pinned SQLite digest only.

`runner.py` and `parser.py` exceed the normal agent-context target because the
scientific release hashes their source and because validation order is itself a
security boundary. Refactor them only in a release-scoped change: preserve the
public imports, update the wheel and release identities, and run the complete
Integration -> backend validator -> public summary -> frontend contract chain.

- `repeat_map.py` owns the isolated CRISPRmap reference stage and input/label gates. The backend independently recomputes nearest ties.

- `repeat_type.py` owns the isolated RepeatTyper stage and occurrence/model/score checks; the backend independently validates its native report.

- `association.py` runs the isolated SpacePHARER engine against operator-pinned target/control references. It validates explicit groups, input occurrence identity, native calibration, target identity and oriented alignments; backend validation remains independent.
- `association_reference.py` verifies the immutable index/catalogue bytes and loads only matched target sequences into a bounded cache. Public evidence carries reference counts and scope, not the full catalogue.

- `array_compare.py` consumes a pinned observed-array map and real spacer FASTA, runs CCTK, and independently checks identities, metrics and components. Backend validation reimplements this gate without imports.

- `locus_compare.py` exports all source-verified called proteins and runs bounded clinker comparisons.
- `prophage_context.py` joins source-derived exact spacer matches to Phigaro intervals.
- `insertion_context.py` and `defense_context.py` retain native evidence and add a separate `context_annotations.py` Cas stage.
- `rna_family.py` validates user-supplied homologous DNA/RNA families, duplicate crosswalks, alignments and comparative folding evidence.
