# CRISPR tool integration

Integration 0.14.0 adds source-bound locus comparison (clinker), source-spacer/prophage context (Phigaro), insertion-sequence context (ISEScan), defense-repertoire context (DefenseFinder) and homologous RNA-family folding (RNAalifold). Each worker invocation is offline and independently validated again by the backend. Cas-context stages preserve complete called-gene evidence; RNA families retain input orientation and explicit user-supplied homology and boundary assumptions.

Integration 0.13.1 additionally rejects ambiguous and nonfinite annotation JSON/JSONL and bounds file reads independently of the native tools. It preserves SpacePHARER's version-2 reference and report contracts, with streamed integrity checks and bounded lazy source-genome validation. See [spacer associations](../docs/spacer-association.md).

Integration 0.6.0 additionally exposes `--tools protospacer`. This invokes the
isolated `crisprspacer` runner command against an operator-installed immutable
spacer database and treats the submitted DNA as target sequences. It performs
exact full-length comparisons on both strands without running array detection
or evolutionary reconstruction. Configure `--spacer-database /absolute/bundle`
and `--spacer-database-sha256 SHA256`; the digest identifies the exact bytes of
`spacers.sqlite3`. Anonymous submissions never choose these arguments.

The native search deadline is capped at 300 seconds and cannot exceed the
configured stage deadline. A completed native result must bind every target,
coordinate, sequence, strand, reference digest, alignment count and retained
input. Input and native output folders remain private. Missing reference data,
digest mismatch, resource limits or invalid output fail the stage; they are not
reported as zero matches. Array overlap remains unknown and PAMs are not
evaluated. The backend separately validates these semantics before publication.

This package is the CRISPRloci v2 scientific orchestration boundary. It converts CRISPRidentify v2 machine output into the shared `ccdb`
input layout and can run the complete workflow:

```text
FASTA → CRISPRidentify → adapter
                              ├─ identify: stop after detection and normalization
                              ├─ spacerplacer: input-orientation reconstruction
                              └─ crispr-evor/both: orientation comparison plus selected SpacerPlacer reconstruction
```

The adapter prefers each sequence's canonical `report.json`. Canonical reports must name the
`CRISPRidentify-report` schema; the explicit major-version adapter accepts major 1 and fails closed on an unknown
major. It also accepts the legacy per-sequence JSON and `Complete_json_summary.json` only when original-input binding
was not requested. It deliberately does not scrape CSV or HTML. That keeps the connection stable when presentation
output changes and makes an incompatible machine-schema change fail clearly instead of silently changing the
scientific input.

## Install and run

The integration package itself has no runtime dependencies. The server release
manager builds its wheel and installs it beside the backend; scientific tools
remain in separate environments. Use the [project release workflow](../README.md).

The worker passes validated per-tool CLI overrides through the bounded private
`CRISPRLOCI_TOOL_OPTIONS` JSON environment variable. Integration consumes and
removes that variable, independently validates its exact flag names, values and
workflow applicability, then builds shell-free argument vectors. This transport
does not add flags to scientific CLIs. Direct Python callers can pass the same
object as `RunnerConfig.tool_options`. See [tool options and modes](../docs/tool-options.md).

Deploying source changes still requires building and gating a new immutable
backend/Integration release; running jobs and existing runtime snapshots are
never modified by an options request.
For standalone adapter development:

```bash
python3.12 -m venv Integration/.venv
Integration/.venv/bin/pip install ./Integration
Integration/.venv/bin/crispr-tools run genomes/ job-output/ \
  --runner-config /absolute/toolset/runners.json --tools loci \
  --tracr-model-type II --cpu 8 --stage-timeout 7200
```

Runner JSON keys are `crispridentify`, `casandra`, `crispr_tracrrna`,
`spacerplacer`, and `crispr_evor`. Values are command argument arrays, never shell
strings. `scripts/provision_tools.py` generates exact isolated executable paths.
The new selections `casandra`, `tracrrna` and `loci` validate CasAndra schema-5
output or CRISPRtracrRNA schema-3 output. Complete loci retain and validate native
array/Cas evidence; model screening uses `II` or `V` (Type V-K). See
[the scientific contract](../docs/scientific-contract.md) for evidence semantics.

Choose the biological analysis with `--tools`: `identify` is the detection-only route suitable for a single genome;
`spacerplacer` reconstructs related arrays in their supplied orientations; `crispr-evor` predicts group-wide
acquisition orientation and includes the orientation-selected SpacerPlacer reconstruction. `both` is an alias for
that scientifically coherent CRISPR-evOr result, not two competing reconstructions. Detection-only jobs finish
successfully even when the adapter finds no comparative group.

The `--cpu` budget is forwarded to the selected evolutionary worker as its independent-group alignment limit.
CRISPR-evOr estimates separate forward and reverse trees when no independent phylogeny is supplied, as required by its
method. The runner never passes SpacerPlacer's input-order-derived estimated tree as though it were independent tree
evidence. A reusable strict-JSON alignment cache remains a supported worker artifact, but the runner does not execute a
complete provisional SpacerPlacer reconstruction merely to produce it.

For FASTA batches, `--identify-folder-jobs N` opts into CRISPRidentify's persistent lanes for independent records.
It does not increase the global `--cpu` budget: the worker clips the effective lane count to the number of input records
and CPUs available within that budget, then partitions those CPUs into disjoint lanes. The Integration default remains
the worker's conservative single lane until deployment benchmarks select an operator value; benchmark 1, 2, and 4 by
passing this option explicitly.

For an already completed CRISPRidentify job:

```bash
crispr-tools run identify-results/ job-output/ \
  --input-kind identify-results \
  --runner-config runners.json
```

To inspect or hand off only the normalized input:

```bash
crispr-tools adapt identify-results/ adapter-output/ --source-fasta original-input.fasta
```

Every fresh FASTA pipeline run requires a one-to-one match between submitted records and canonical reports. The
adapter recomputes each record's ID, length, and SHA-256 over its whitespace-free uppercase sequence, then compares
all three with `source.id`, `source.length`, and `source.sha256`; `source.sha256_scope` must be
`uppercase_sequence`. Stale, duplicated, missing, or cross-job reports fail before any evolutionary worker starts.
For a deliberate legacy-results import where the original FASTA is unavailable, omit `--source-fasta`; the manifest
records `source_binding.status` as `not_requested` rather than implying verification.


The downstream tools' alpha/rho bias corrections remain enabled by default. Some small or sparse groups have no
finite optimum for the optional unobserved-spacer rho correction even though the primary likelihood fit is valid.
Integrated SpacerPlacer and CRISPR-evOr runs therefore default to
`--rho-bias-failure-policy fallback_uncorrected`: correction is attempted independently for each group and, for
CRISPR-evOr, each direction. Only an affected fit uses uncorrected rho estimates, and the requested alpha correction
is retained. This does not change the BDM likelihood used for orientation, and only the affected rho/deletion-rate estimate is uncorrected. Optimizer likelihoods and the IDM/BDM LRT are unchanged by rho correction. SpacerPlacer's CSV, `crispr_evor_results.json`, and `integration_result.json` expose the
requested policy, effective status, stable warning code, and full warning. Select `--rho-bias-failure-policy fail` for
strict all-corrections-or-error behavior. Direct worker CLIs also retain strict mode as their own default. Other
optimizer failures are never eligible for fallback. `--disable-bias-corrections` remains the explicit global opt-out
and is recorded in the worker command.
The default accepted categories are `Bona-fide` and `Possible`. `Low score` is excluded unless explicitly named with
`--categories`. Groups require at least two arrays because the downstream tools perform comparative evolutionary
reconstruction.

For FASTA jobs, the default category policy is pushed upstream as CRISPRidentify `--accepted_only True`; this reduced
the bundled example's canonical JSON from 77,318,382 bytes to 41,589 bytes without changing its two accepted arrays.
Explicitly including `Low score` automatically switches CRISPRidentify back to the full report.

## Output contract

Each job writes:

```text
job-output/
├── crispridentify/             # omitted when existing results are reused
├── adapter/
│   ├── ccdb/
│   └── manifest.json
├── spacerplacer/                # only for --tools spacerplacer
├── crispr_evor/                 # for --tools crispr-evor or both
├── logs/
└── integration_result.json
```

`adapter/manifest.json` records every included/skipped array, source-binding status, the source report, orientation
handling, grouping policy, and output file. Its `input` and per-array `report_path` values are portable: they are relative to the adapter input
root (a directory input is represented by `.`), never host-absolute `/srv` or `/opt` paths. Schema 1.1.0 records this
rule in `path_contract`; the defensive fallback for a report found outside that root is
`external/<safe-basename>`. `integration_result.json` records stage status, commands, logs, and the CRISPR-evOr tree policy. For
`--tools both`, its SpacerPlacer stage has status `fulfilled_by_crispr_evor` and points to
`crispr_evor/0_results.csv`; no provisional `spacerplacer/` result directory is created.

A zero worker exit status is necessary but not sufficient for stage completion. Standalone SpacerPlacer must produce a
parseable `0_results.csv` with the stable reconstruction columns and exactly one uniquely named row for every group in
the adapter manifest. Its rates, tree length, and runtime must be finite and nonnegative; its BDM log likelihood must
be finite; its aligned-spacer, reconstructed-insertion, and reconstructed-deletion counts must be nonnegative
integers; and its preferred deletion model must be exactly `BDM` or `IDM`. CRISPR-evOr applies the same validation to
the orientation-selected reconstruction and must additionally produce parseable `crispr_evor_results.json` and
`additional_data/orientation_trees.json`. Its CSV rows, JSON results, and forward/reverse/selected tree maps must name
exactly the manifest groups. Every generated tree is parsed as restricted Newick with at least two unique leaves,
finite nonnegative branch lengths, and exactly the leaf identities declared by the adapter. The selected tree must be
the reverse tree only for a `Reverse` prediction and the forward tree for `Forward` or `ND`. Its serialized branch sum
must agree with the CSV tree length within the cumulative five-decimal Newick rounding bound.

The CSV orientation, recommendation, threshold, and forward/reverse BDM likelihoods are cross-checked against the JSON
values; the likelihood difference is independently recomputed and the prediction is independently derived from that
difference and the threshold. For `Reverse`, the main CSV's `ln_lh_bdm` is the selected reverse reconstruction and
`reversed_ln_lh_bdm` is the original forward reconstruction, while the difference column remains in the canonical
forward-minus-reverse coordinates recorded by the compact JSON. This deliberate selected-row convention is validated
explicitly. The contract rejects header-only, partial, duplicated, stale, cross-job, or scientifically inconsistent
output even when the worker returned zero. Such failures change the stage to `failed_output_validation` with a
stage-specific error. The backend should use `crispr_evor_results.json` as the canonical compact orientation response
and expose the CSV as the detailed selected-reconstruction artifact.

Standalone SpacerPlacer rows must additionally report the exact requested rho-failure policy and the effective
alpha/rho correction status. A fallback row is accepted only under the explicit fallback policy and only with the
stable code `rho_bias_fit_unavailable_used_uncorrected` plus a non-empty warning; warning fields on a non-fallback row
are rejected. The stage summary exposes `bias_correction_policy`, `bias_fallback_count`, and `bias_warnings`, using
`input_order` as the orientation label because this route does not infer orientation.

CRISPRidentify v2 reports source-forward spacer sequences and stores predicted biological orientation separately.
Source sequences are therefore written unchanged under `pos_strand` or `neg_strand`; the downstream CCDB reader does
the one required reverse-complement/order operation. If a future report explicitly supplies biologically oriented
sequences, they are stored under `pos_strand` to prevent a double reversal. A source-frame array whose orientation was
not determined is excluded with manifest reason `orientation_not_determined`: CRISPR-evOr reverses a group globally
and cannot repair one independently misoriented member. An explicitly biologically oriented sequence remains usable
without a separate strand label because its sequence frame is already resolved.

## Grouping policy

Arrays are grouped only when they have the same direct repeat after reverse-complement canonicalization and belong to
the same transitive spacer-overlap component. Spacer identity for this grouping is also reverse-complement neutral.
Downstream `--cluster_spacers` then merges near-identical spacers with the configured edit distance (default 1). This
mirrors the conservative CRISPRCasFinder conversion policy and avoids comparing unrelated CRISPR families.

## Web-server boundary

- Allocate a new empty output directory for every job; resuming into a non-empty directory is rejected.
- Keep the runner configuration under operator control. Do not accept command arrays from public requests.
- Validate/extract uploads in a separate unprivileged job sandbox and enforce CPU, memory, file, process, and wall-time
  limits outside this Python process as well as the per-stage deadline.
- Treat `completed_no_eligible_groups` as a valid comparative-analysis no-result, not a server failure. Detection-only
  jobs use `completed` even when there is no comparative group.
- Keep `adapter/manifest.json`, `integration_result.json`, CRISPR-evOr's JSON result, and the tree artifacts for
  reproducibility.

## Repeat evidence

Version 0.5.0 adds `--tools repeats --molecule DNA|RNA` for supplied repeats and
`--tools repeat-context` for detected genomic arrays. Both invoke an isolated
`crisprrepeat` command from runner configuration; no ViennaRNA import occurs in
Integration. The report is source-validated before the completed marker is emitted.
See [web/tool contract](../docs/repeat-analysis.md).

RepeatTyper is integrated as the queued `repeat_type` mode, consuming repeat DNA
FASTA or the existing repeat occurrence handoff. See [repeat typing](../docs/repeat-typing.md) for immutable model identity,
uncertainty semantics, the isolated runtime and validated exports.

SpacePHARER is integrated as `spacer_association`: explicitly grouped spacer occurrences, an operator-owned targeted phage panel, native empirical calibration and protected evidence exports. See [spacer association](../docs/spacer-association.md).

CCTK is integrated as queued `array_compare`: detection arrays.json feeds observed
spacer order and sharing comparisons. See [array comparison](../docs/array-comparison.md).
