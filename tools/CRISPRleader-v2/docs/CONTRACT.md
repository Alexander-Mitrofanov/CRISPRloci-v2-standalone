# Contract freeze v1 — 27 September 2026

Owner confirmed the council plan. This commit is the shared worker baseline.
Implementation is standalone; no server tool lock/runtime changes.

Shared records are in `models.py`. Propose changes to the coordinator before
editing this file. Native JSON uses zero-based half-open source-forward ordered
segments. Strands are separate. Arrays/anchors may cross a circular origin;
all segments form one forward arc. Core intervals in RawEvidence/LeaderDecision
are offsets in the oriented flank string. Motif intervals are offsets in the
oriented flank plus observed terminal-repeat string. Intervals must reconstruct.

## Native annotations v1

Object keys: `schema_version` = `1.0.0`, `coordinate_system` =
`zero_based_half_open`, `sources` and `arrays`.

Source objects: `id`, `length`, `sha256`, `topology` (linear/circular/unknown),
`topology_provenance` (required for circular), optional `metadata`.
Array objects match `ArrayAnnotation` serialized with `json_record`. Every
terminal-repeat anchor's sequence is source-forward. Consensus is source-forward
relative to the array arc, independently of a tool's transcription strand.
Unknown input fields fail closed unless within explicit metadata.

## W1 interfaces

`sequence_io.load_inputs(fasta_path, annotation_path, config=None)` returns
`(dict[str, SourceRecord], list[ArrayAnnotation])`. Strict bounded source binding,
IDs, SHA256, geometry and terminal repeats. Raise InputError for invalid inputs.
`sequence_io.load_fasta(path, max_bases=50000000)` returns source dict.
`adapters.identify.adapt_identify(fasta_path, report_paths)` returns native
annotation dict (accept Path list). Do not import other scientific packages.

`coordinates.extract_candidates(source, array, config=None)` returns both
FlankCandidates; default window500 is extraction policy, not inferred boundary.
`coordinates.reverse_complement(sequence)` supports IUPAC DNA.
`coordinates.extract_segments(source_sequence, segments, strand='+')` reconstructs.
`coordinates.project_interval(candidate, start, end, include_repeat=False)`
returns tuple[Segment,...] in source-forward arc order. For include_repeat,
candidate.junction_sequence is the coordinate frame; otherwise candidate.sequence.
Right candidate has negative strand; its source-forward junction arc is anchor
then flank (reverse-complement yields flank then observed repeat).

W1 owns tests/test_coordinates.py, test_sequence_io.py, test_identify.py only.

## W2 interfaces

`references.load_bundle(path=None)` returns a validated plain mapping. None
loads package default bundle. Required missing/corrupt/incompatible references
raise ReferenceError, never a biological unsupported result. Bundle schema v1
contains `schema_version`, `bundle_id`, `version`, `families`, `motifs`,
`provenance`, `files` checksum mapping, `policy` calibrated controls.
Reference worker chooses nested family/motif fields with documented schema.

`inference.analyze_array(array, candidates, bundle, config=None, runner=None)`
returns LeaderDecision. It preserves all candidate raw evidence, motifs and
no-call reasons, and supplied vs supported orientation. Root maps intervals.
No calls based solely on flank extraction. Scores retain units. Motifs can
span leader-repeat junctions; actual observed repeat is mandatory there.
Inference never downloads. Reference build is a maintainer workflow.

W2 owns src/crisprleader/inference/, references/, data/, tests/test_inference.py,
test_references.py, test_motifs.py, scripts/build_references.py. Coordinate and
CLI implementation remains elsewhere. Required biological default reference
needs source terms, reproducible sequence acquisition and independent labels.

## W3 interfaces

`runtime.run_command(argv, *, cwd, timeout_seconds=120.0, max_output_bytes=1048576,
cancel_event=None, env=None)` returns CompletedProcess-like record with stdout,
stderr, returncode. No shell. New process group; timeout/cancellation/output
overflow kill group and raise ExecutionError. Inputs need not use global cwd.
Output bytes bound includes both pipes. Explain process cleanup race handling.

W3 owns runtime/, benchmarks/, tests/test_runtime.py and independent end-to-end
biological fixtures. W3 designs development/test ledger and predeclares scientific
acceptance gates before final runs, independent of W2 prediction output. Share
training/calibration source IDs with W2; keep final labels out of tuning.

## Coordinator interfaces

CLI predict/adapt/validate/references verify; independent validator; JSON/TSV/
FASTA/GFF3/BED exports; schemas; packaging/CI/docs. Result root v1 has version,
execution_status, coordinate_system, sources (identity without entire genome),
arrays (native evidence), candidates, predictions, provenance/reference receipt,
warnings and effective_config. Selected cores and motifs have explicit source
segments and sequence. GFF3/BED wrap features use linked pieces. Core feature
labels never imply full functional leader or activity.

Each worker commits only owned files and reports test commands, limitations,
new dependencies and contract requests. Root merges reviewed commits. No
external asset with undocumented terms is copied into the default bundle.

## Implemented snapshot versus intended contract

The contract above records the approved target interfaces. As of this snapshot,
`inference.analyze_array` is **not implemented** and no default biological bundle
is installed. `pipeline.predict` raises an explicit execution error. See STATUS.md.

`pipeline.extract(fasta, annotations, config=None)` produces a validated report
with `analysis_mode=context_extraction`, no reference and no predictions. Both
flanks remain present even when one is empty because of clipping. A future
prediction result must declare `analysis_mode=prediction` and retain one decision
per array. Validation checks coordinate consistency; it does not independently
verify external biological claims. CLI `extract` and its exports are operational.
