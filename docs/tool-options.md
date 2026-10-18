# Tool options and modes

This reference describes the scientific command lines pinned in
[`tools/tools.lock.json`](../tools/tools.lock.json). An option's existence in a
standalone CLI does not mean it is available in every web workflow. The advanced
settings menu shows the allowlisted options applicable to the selected workflow;
the server validates them again before an isolated worker constructs an argument
vector. There is no shell command or arbitrary-arguments input.

The tables below document the **full scientific CLI inventory**, including
options that are not editable in the server. The precise editable subset, its
workflow applicability, defaults and server bounds are defined in
[`tool-options-v1.json`](../frontend/src/contracts/tool-options-v1.json).
Additional existing workflow controls select model family, flank length, viral
mismatch allowance, repeat molecule, spacer clustering distance and bias policy.
The remaining CLI entries are reference information, not promises that the web
service accepts them.

Defaults below are the **shipped CLI defaults**. The web service can apply a
documented workflow or operator policy, including stricter resource limits.
Leaving an advanced field unchanged preserves the workflow's normal defaults.
Changing a scientific parameter changes the invoked analysis configuration; it
does not guarantee a different biological result for every input.

## Sources and version policy

The following pinned source files define accepted flag spellings, values and
defaults. Read the specified commit, not an independently updated checkout's HEAD.

| Tool | Locked version / commit | Authoritative CLI source |
| --- | --- | --- |
| CRISPRidentify-v2 | 2.3.2 / `e8a3d7c` | [CRISPRidentify.py](https://github.com/Alexander-Mitrofanov/CRISPRidentify-v2/blob/e8a3d7c038eab27bdb955f97c5cd6bc3666a17dc/CRISPRidentify.py) |
| CasAndra | 0.3.0.dev2 / `ba046f9` | [src/casandra/cli.py](https://github.com/Alexander-Mitrofanov/CasAndra/blob/ba046f90242ac91dc075183ea51dcbedfac28484/src/casandra/cli.py) |
| CRISPRtracrRNA-v3 | 3.0.0 / `1e173bc` | [src/crisprtracrrna/cli.py](https://github.com/Alexander-Mitrofanov/CRISPRtracrRNA-v3/blob/1e173bcd4e7d59b86412dac17a0f6ebcb4353eef/src/crisprtracrrna/cli.py) and [config.py](https://github.com/Alexander-Mitrofanov/CRISPRtracrRNA-v3/blob/1e173bcd4e7d59b86412dac17a0f6ebcb4353eef/src/crisprtracrrna/config.py) |
| SpacerPlacer | 1.0.1.post2 / `2742b6e` | [spacerplacer.py](https://github.com/Alexander-Mitrofanov/SpacerPlacer/blob/2742b6ec5336fcbc40add772e1128441807ae239/spacerplacer.py) |
| evOr | 1.0.2 / `5310ba9` | [crispr_evor.py](https://github.com/Alexander-Mitrofanov/evOr/blob/5310ba9e75836717c58d037606de17f98fc5d7b1/crispr_evor.py) and [spacerplacer.py](https://github.com/Alexander-Mitrofanov/evOr/blob/5310ba9e75836717c58d037606de17f98fc5d7b1/spacerplacer.py) |
| CRISPRleader-v2 | 2.0.0a1 / `da7b4d6` | [src/crisprleader/cli.py](https://github.com/Alexander-Mitrofanov/CRISPRleader-v2/blob/da7b4d659afaa8b0eb8873ebc2ccec00aad1000f/src/crisprleader/cli.py) |
| CRISPRrepeat | 1.0.0 / `4a98cc4` | [src/crisprrepeat/cli.py](https://github.com/Alexander-Mitrofanov/CRISPRrepeat/blob/4a98cc465656c4db3ed0faccaafe476707d20ddf/src/crisprrepeat/cli.py) and [folding.py](https://github.com/Alexander-Mitrofanov/CRISPRrepeat/blob/4a98cc465656c4db3ed0faccaafe476707d20ddf/src/crisprrepeat/folding.py) |
| CRISPRspacer | 0.3.0 / `c331f8b` | [src/crisprspacer/cli.py](https://github.com/Alexander-Mitrofanov/CRISPRspacer/blob/c331f8b3b26bb4cea968d288f024c0e3aadd3a77/src/crisprspacer/cli.py) |
| CRISPRmap-v2 | 2.0.0a1 / `b35c245` | [src/crisprmap/cli.py](https://github.com/Alexander-Mitrofanov/CRISPRmap-v2/blob/b35c245e1103598fb418ebd8126f713126cda08c/src/crisprmap/cli.py) |

The online documentation audit on 2026-10-09 also consulted the accessible primary
[CRISPRidentify documentation](https://github.com/BackofenLab/CRISPRidentify),
[CRISPRtracrRNA documentation](https://github.com/BackofenLab/CRISPRtracrRNA),
[SpacerPlacer documentation](https://github.com/fbaumdicker/SpacerPlacer), and the
[CRISPR-evOr paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12638026/).
The newer fork pages were not retrievable by that browser; their pinned local Git
objects supplied the authoritative definitions. Older upstream documentation is
background, not permission to introduce unshipped flags. For example, an upstream
SpacerPlacer README mentions `--orientation_threshold`, but the locked executable
accepts `--orientation_decision_boundary` instead.

## CRISPRidentify: array detection

The executable takes flags directly, without an analysis subcommand. It is used
by detection and workflows that first require detected arrays. The complete
tracrRNA workflow invokes its own CRISPRidentify stage; that stage accepts the
tracrRNA `--identify-preset` interface, not arbitrary outer detection settings.

| Flag | CLI default / values | Effect |
| --- | --- | --- |
| `--preset` | `balanced`; `balanced`, `region-screen`, `full` | Selects the shipped analysis configuration. Region screening changes accepted representations and boundary behavior; full uses the historical full search. |
| `--min_len_rep`, `--max_len_rep` | `21`, `55` | Lower and upper accepted average repeat lengths. |
| `--min_len_spacer`, `--max_len_spacer` | `18`, `78` | Lower and upper accepted average spacer lengths. |
| `--min_repeats` | `3` | Minimum repeat count. |
| `--strand` | `True` | Requests array orientation prediction. |
| `--degenerated` | `True` | Enables degenerate and terminal-repeat enhancement; disabling also disables architecture rescue. |
| `--enhancement_max_min` | `True` | Enables enhancement with the length filter. |
| `--enhancement_start_end` | `True` | Enables start/end omission enhancement. |
| `--enhancement_intermediate` | `True` | Enables intermediate-substring enhancement. |
| `--max_identical_spacers` | `4` | Maximum identical spacers in an array. |
| `--max_identical_cluster_spacers` | `3` | Maximum consecutive identical spacers. |
| `--margin_degenerated` | `30` | Spacer-margin length considered during degenerate search. |
| `--max_edit_distance_enhanced` | `6` | Edit-distance limit for enhancement candidates. |
| `--fast_run` | `True` | Selects the high-recall fast path; false selects historical full search. |
| `--fast_run_seed_profile` | `balanced`; `balanced`, `legacy`, `sensitive` | Selects the fast-path seed/gap search domain. Applies to fast mode. |
| `--fast_run_rescue` | `True` | Enables bounded classifier-routed start/end rescue. Applies to fast mode. |
| `--fast_run_rescue_seed_scope` | `cluster`; `lineage`, `locus`, `family`, `cluster` | Chooses which observed seed hypotheses bounded rescue can consider. |
| `--fast_run_rescue_trim_profile` | `full`; `full`, `total3`, `one_sided` | Chooses terminal-trimming hypotheses considered by rescue. |
| `--fast_run_phase_repair` | `True` | Enables bounded source-backed phase repair. Applies to fast mode. |
| `--fast_two_repeat_terminal_repair` | `False` | Experimental bounded terminal completion for narrow two-repeat accepted-only candidates. |
| `--fast_accepted_family_rescue` | `False` | Experimental accepted-family rescue; the CLI says independent biological holdout is pending. |
| `--split_competitive_components` | `True` | Selects independent accepted winners for disjoint source components; required by default fast rescue. |
| `--reconstruction_policy` / `--reconstruction-policy` | `conservative` in balanced; `legacy` in other presets | Controls evidence-supported array reconstruction and preservation of ambiguous phases. |
| `--architecture_rescue` | `False` | Enables source-backed assembly and coding-context adjudication; requires degenerated enhancement and legacy reconstruction policy. |
| `--fuzzy_search_mode` | `native`; `full`, `family`, `native` | Chooses fuzzy repeat-search strategy. |
| `--fuzzy_search_flank` | `100` | Bases added around each Vmatch cluster. |
| `--fuzzy_family_padding` | `15` | Bases added around family-anchor windows. Applies to the relevant family search. |
| `--fuzzy_boundary_paths` | `True` | Includes forced-boundary alternatives in fuzzy search. |
| `--accepted_only` | `False` | Retains only Bona-fide and Possible arrays in reports. |
| `--biological_loci` | `False` | Adds heuristic locus presentation groups; it does not supply new confidence values. |
| `--cas`, `--cas_gene_mode` | `False`, `auto`; gene mode `auto`, `single`, `meta` | Optional embedded CasAndra annotation and its gene-calling mode. Web workflows select their own Cas stage. |
| `--is_element` | `False` | Optional IS-element annotation, requiring a separately provisioned environment. |

Boolean flags here take explicit values, for example `--degenerated False`.
They are not presence-only switches. `balanced` permits explicit customization;
`region-screen` and `full` are atomic presets and reject contradictory controlled
flags. Repeat/spacer length and repeat-count bounds are not atomic preset fields.
The advanced web allowlist must reject incompatible combinations before queueing.

`--output_profile`, `--fasta_report`, `--html_report`, `--json_report`,
`--pickle_report` and `--retain_candidate_evidence` govern output rather than
scientific search. The service owns output paths and its validated report format.
`--sequence_topology_manifest` requires independently verified source provenance.
`--cpu`, `--parallel`, `--folder_jobs` and `--cas_timeout_seconds` are scheduling
controls. Hidden `--_folder_*` switches and the suppressed `--cas_completeness`
compatibility option are not public scientific controls. `--print-effective-config`,
`--version`, verbosity flags and help do not run analysis.

## CasAndra: Cas annotation

| CLI mode | Purpose / web applicability |
| --- | --- |
| `predict-genome` | Calls or reuses CDSs and predicts Cas cassettes; used by nucleotide Cas workflows. |
| `predict-proteins` | Classifies supplied amino-acid sequences; separate standalone input contract. |
| `annotate-proteins` | Publishes provenance-bound amino-acid annotations; separate standalone input contract. |
| `classify-cassette` | Classifies one ordered protein cassette; separate standalone input contract. |
| `inspect-model` | Operator inspection of a local model bundle. |

`predict-genome --gene-mode` accepts `auto` (default), `single` and `meta`, changing
how Pyrodigal selects/trains its gene caller. `--translation-table` defaults to
11. Its effect depends on the CDS source and selected caller: single-genome
training uses the requested table, whereas the fast metagenomic path selects
codes per gene. It must not be offered as a universally effective metagenomic
control. GenBank can carry source-specific translation annotations.

`--model` is an operator-owned frozen bundle path. `--genome`, `--genbank`,
`--gff3`, `--input` and `--output` select input contracts or locations, not tunable
scientific thresholds. `--threads` defaults to 6 but the worker owns allocation.
The shipped protein commands have no sensitivity, score-cutoff or subtype-threshold
flag; the menu must not invent one.

## CRISPRtracrRNA: model screening and complete analysis

`run` performs analysis; `doctor` inspects the runtime. The shipped `--mode` values
are `complete` and `model`, with compatibility aliases `complete_run` and
`model_run` (also `--run_type`). The standalone default is `complete`; the web
tracrRNA workflow selects `model`, and the combined loci workflow selects
`complete`. A mode change requires the corresponding validated workflow.

| Flag | CLI/config default / values | Effective analysis |
| --- | --- | --- |
| `--model-type` / `--model_type` | `II`; `II`, `V` | Both modes: selects the packaged covariance-model family. |
| `--merge-gap` | `150` | Model mode: maximum interval gap for merging model hits. |
| `--type-v-motifs` | `False` | Both modes with model type V: adds Type V anti-repeat motif evidence. Accepts true/false. |
| `--similarity` | `0.7`, from 0 to 1 | Complete mode: minimum anti-repeat sequence similarity. |
| `--coverage` | `0.6`, from 0 to 1 | Complete mode: minimum anti-repeat coverage. |
| `--tail-length` | `150`, positive integer | Complete mode: candidate tail extension length. |
| `--upstream-length` | `5`, nonnegative integer | Complete mode: upstream candidate extension length. |
| `--search-backend` | `fasta36`; `fasta36`, `blastn` | Complete mode: selects the anti-repeat search executable. Requires a provisioned backend. |
| `--model-rescue` / `--no-model-rescue` | enabled | Complete mode: adds independent genomic covariance-model rescue. |
| `--model-rescue-evalue` | `1e-5`, positive | Complete mode with rescue enabled: rescue-hit E-value ceiling. |
| `--mask-arrays` / `--no-mask-arrays` | enabled | Complete mode: masks known arrays during anti-repeat search. |
| `--identify-preset` | `balanced`; `balanced`, `region-screen`, `full` | Complete mode: CRISPRidentify preset used by the internal detector. |
| `--gene-mode` | `auto`; `auto`, `single`, `meta` | Complete mode: CasAndra gene-calling mode. |
| `--short-anti-repeat-rescue` / `--no-short-anti-repeat-rescue` | disabled | Experimental rescue; requires complete mode and model type II. |

Model mode returns immediately after covariance-model scanning and hit merging.
It does not execute anti-repeat searching, array detection, Cas calling, tail
extension or complete-mode rescue. Such controls must be hidden or rejected for
model-only submissions. `--merge-gap` is not used by the complete-mode candidate
path. The compatibility aliases `--anti_repeat_similarity_threshold`,
`--anti_repeat_coverage_threshold` and `--perform_type_v_anti_repeat_analysis` are
shipped, but the menu can use the canonical spellings above.

`--model` adds local model paths and `--config` loads a local TOML configuration;
both remain operator-owned. The server owns input/output paths, `--threads`,
`--timeout` and `--intarna-batch-size` (default 128). Ranking weights exist in the
configuration object, but this pinned CLI does not ship individual weight flags.

## SpacerPlacer and evOr: reconstruction and orientation

Both tools accept positional input/output paths and the scientific flags below.
SpacerPlacer reconstructs ancestral arrays. `crispr-evor` uses its own isolated
environment and always enables forward/reverse orientation comparison. The two
packages each own `sp_model` and cannot share one import environment.

| Flag | CLI default / values | Effect |
| --- | --- | --- |
| `--cluster_spacers` | absent / false | Enables spacer clustering before alignment. The web comparative pipeline uses clustering. |
| `--cluster_spacers_max_distance` | `1`, nonnegative | Maximum Levenshtein distance for clustered spacer identities; requires clustering. |
| `--seed` | `2357` | Random seed for MAFFT alignment, not a general biological-confidence parameter. |
| `--insertion_rate` | `0.00001`, positive | Insertion rate used during ancestral reconstruction. |
| `--deletion_rate` | `0.1`, positive | Deletion rate used during ancestral reconstruction. |
| `--extend_branches` | `0`, nonnegative | Adds branch length, allowing event placement on short or zero-length branches. |
| `--tree_construction_method` | `upgma`; `upgma`, `nj` | Tree estimation method when no tree is supplied. |
| `--tree_insertion_rate`, `--tree_deletion_rate`, `--tree_alpha` | automatic, positive when supplied | Overrides tree-estimation likelihood parameters. |
| `--tree_lh_fct` | `simple`; `simple`, `ode_based` | Likelihood function used in tree estimation. |
| `--combine_non_unique_arrays` | absent / false | Combines arrays with identical spacer content before estimation/reconstruction; changes leaf identity and requires matching output validation. |
| `--lh_fct` | `simple`; `simple`, `ode_based` | Likelihood function for deletion-model fitting and comparison. |
| `--no_alpha_bias_correction` | absent / false | Disables alpha bias correction. Corrections are available for the simple likelihood. |
| `--no_rho_bias_correction` | absent / false | Disables rho bias correction. Corrections are available for the simple likelihood. |
| `--rho_bias_failure_policy` | `fail`; `fail`, `fallback_uncorrected` | Fails or explicitly retries affected fits without rho correction, recording the effective policy. Operator-owned in the web service. |
| `--significance_level` | `0.05`, between 0 and 1 | Significance level for the independent/block deletion likelihood-ratio test. |
| `--determine_orientation` | false for SpacerPlacer; forced true by evOr | Runs both orientations. The web workflow owns this selection. |
| `--orientation_decision_boundary` | `5`, nonnegative | Absolute log-likelihood difference required for an orientation call; applies only with comparison enabled. Below the threshold, orientation remains undetermined. |

`--rec_model gtr` and `--tree_distance_function likelihood` each have only one
implemented value; there is no useful mode selector to expose. `--min_evidence_level`
(default 4; values 1–4) applies only to CRISPRCasFinder input, not the web adapter's
CRISPRCasDB-style array groups (`--input_type ccdb`). `--input_type`, `--tree_path` and `--aligned_groups_path` alter input
contracts or reuse operator-owned data. When a validated alignment is reused,
MAFFT's seed does not change that alignment.

`--no_plot_reconstruction`, `--no_plot_order_graph`, `--dpi_rec`, `--figsize_rec`,
`--show_spacer_name_row`, `--do_show` and `--save_reconstructed_events` affect
exports/visualization, not inference. The worker uses its fixed headless output
policy. `--alignment_workers` is resource allocation. `--allow_unsafe_pickle`
and `--allow_existing_output` are deliberately unavailable to anonymous jobs.

## CRISPRleader: adjacent flank extraction

The pinned modes are `extract`, `predict`, `adapt identify`, `validate` and
`references verify`. Only `extract` performs the supported web scientific task:
source-verified extraction of both adjacent windows. `predict` is an explicitly
unavailable interface and must not appear as a working prediction mode.

`extract --flank-length` defaults to **500 bases** and changes how much adjacent
sequence is retained on each side. It does not select a leader side or predict a
biological boundary. The web service limits the value to 1–5000 bases. Input
FASTA, array annotations and output paths are worker-owned. `predict --references`
and `predict --threads` do not activate a missing scientific implementation.

## CRISPRrepeat: isolated repeat and array-context structure

This CLI has one analysis entry point. Without `--reports`, it folds supplied
repeat sequences. With canonical CRISPRidentify `--reports`, it extracts repeats
from original genome FASTA and compares isolated and local array-context folding.

| Flag | CLI default / values | Effect |
| --- | --- | --- |
| `--molecule` | `DNA`; `DNA`, `RNA` | Interprets supplied sequence and determines default orientation. RNA inputs must already be transcribed 5′ to 3′. |
| `--orientation` | `both` for DNA, `forward` for RNA; `both`, `forward`, `reverse` | Selects repeat orientations to evaluate. Source/report validation still applies. |
| `--temperature` | `37`; 0–100 °C | Changes the thermodynamic folding model temperature in both workflows. |
| `--no-lonely-pairs` | absent / false | Enables ViennaRNA's no-lonely-pair policy in both workflows. |
| `--window` | `150`; 10–300 | Local-context folding window only; no effect on isolated supplied-repeat folding. |
| `--span` | `80`; 5 through window | Maximum local base-pair span only; no effect on isolated supplied-repeat folding. |
| `--categories` | `Bona-fide Possible` | Selects reported array categories when genome reports are supplied. |

For short contexts, the effective window/span are clipped to the available
sequence length and reported. The engine uses pinned ViennaRNA 2.7.0 and
RNA Turner2004 parameters; DNA repeats are evaluated as RNA sequence hypotheses.
Structure evidence does not establish transcription orientation or functional
activity. `--reference`, `--reports`, input/output paths and `--timeout` (default
300 seconds, accepted 1–3600) remain part of the managed worker contract.

## CRISPRspacer: exact target search and viral reference search

The web **protospacer** workflow invokes `database-search`: submitted target DNA
is searched against an immutable observed-spacer SQLite database using exact
full-length matching. There is no mismatch, substitution or sensitivity flag for
this command. Its shipped knobs are the output limits `--max-hits` (10,000),
`--max-candidates` (2,000,000), `--timeout` (300), and managed database/target/output
paths and expected database digest. Resource limits are not scientific thresholds.

The web **viral search** workflow invokes `viral-search` on supplied spacers.
`--max-mismatches` accepts **0, 1 or 2**, default **2**, and filters full-query
matches by substitution count. BLAST retrieval remains heuristic, even at zero;
this mode is not the exact observed-spacer database search. `--threads` (2),
`--timeout` (300), `--max-hits` (10,000), `--max-hsps` (200,000), the BLAST executable,
reference path/digest and output location are managed by the worker.

Other shipped standalone analysis modes have different input/output contracts:

| CLI mode / flags | Default / choices | Effect |
| --- | --- | --- |
| `prepare` | Source FASTA plus canonical identify reports | Validates and exports a spacer-occurrence inventory. |
| `self --scope` | `same-record`; `same-record`, `same-assembly`, `all` | Chooses source target scope for self-targeting analysis. |
| `search` | Supplied spacers/inventory plus targets/reference | Searches explicitly supplied target sequences. |
| `self` / `search --substitutions` | `0`; 0, 1, 2 | Allowed full-length substitutions in these commands only. |
| `--min-length`, `--max-length` | `1`, `1000` | Spacer length eligibility for self/search. |
| `--category` | Bona-fide and Possible | Repeatable array-category selection for inventories. |
| `--flank` | `10` | Forward-target context retained around hits; no PAM inference. |
| `--array-padding` | `0` | Additional exclusion bases around known arrays. |

The self/search budgets are `--max-hits` 100,000, `--max-candidates` 2,000,000,
`--max-pairs` 100,000, `--max-evidence-bytes` 67,108,864, `--max-work` 500,000,000 and
`--timeout` 300 seconds. These modes are not substitutions for `database-search`.
They need their own web input and validation contracts before exposure.

`index`, `database-index`, `viral-plan`, `viral-download` and `viral-index` are
operator reference-preparation commands. `validate`, `validate-reference`,
`validate-database`, `validate-database-result`, `validate-viral-reference` and
`validate-viral-result` are validation commands. A submitted job cannot download,
replace or rebuild a reference.

## CRISPRmap: frozen repeat-reference mapping

The locked CLI has `classify`, `build-reference` and `inspect-reference` commands.
The web workflow invokes `classify` against its pinned historical reference.
Newer independent checkout commands such as `catalog-search` are not part of
this locked server release.

| Flag | CLI default / accepted range | Effect |
| --- | --- | --- |
| `--max-distance` | `3`; standalone 0–200 | Global unit-cost edit cutoff; zero gives exact recovery. The server can impose a lower maximum. Near matches remain uncalibrated candidates. |
| `--top-k` | `20`; 1–1000 | Maximum displayed nearest hits. All equal-best ties still contribute to annotation evidence before this cap. |
| `--max-comparisons` | `5,000,000` | Work reservation; exceeding the budget fails explicitly. Operator-owned. |
| `--deadline-seconds` | `300` | Analysis deadline; operator-owned. |

The input/output paths, reference and provenance are managed by the worker.
`build-reference` takes FASTA, metadata, namespace and source URL and is an
operator-only publishing operation. `inspect-reference` verifies the reference.
No shipped flag assigns a Cas subtype, converts similarity into confidence,
infers transcription strand or replaces conflicting/missing annotations.

## Verification requirements

Every exposed setting must be accepted by its locked CLI, retained through job
storage and passed to the corresponding worker argument vector. Tests must cover
the end-to-end submission/claim/command path, invalid values and unsupported
tool/workflow combinations. Check real behavior on bounded synthetic fixtures
where available, such as repeat folding at two temperatures, two flank lengths,
or exact versus one-edit repeat reference retrieval. Synthetic behavior tests
establish wiring and parameter effects, not improved biological accuracy.

Independent backend output validation must use the effective submitted policy;
it must not keep hard-coded defaults that reject valid customized results. A
recorded command alone cannot compensate for a setting ignored by the selected
upstream mode. Browser checks should confirm controls are collapsed initially,
only relevant tool settings are shown, changed values are submitted and reset
returns to defaults.

CCTK (`array_compare`) exposes `--identity` (exact or explicit reverse-complement)
and `--min-shared` (1–500 distinct spacers). The worker fixes supplied array order
and seed 2, owns input/output paths, and keeps the native PNG private.
