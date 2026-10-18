# Detection-science council review, phase 1

Date: 2026-09-27. Scope: a bounded next implementation in `new_version`, followed by review by this same council member. No core code, model, external repository, or `old_version` was changed for this review.

## Decision

Implement an **opt-in Type II short anti-repeat rescue branch**, with an anchored duplex check and separate supporting evidence, while retaining the current standard branch and `model_rescue_evalue=1e-5`. Add explicit evidence classes and a non-destructive locus view. Do not make a lower global coverage or CM cutoff the new default. The experiment should test a mechanism that is currently excluded before scoring; success is not defined as recovering all ten known misses.

The first version should reuse native genome-wide alignment output. A seed absent from that output remains absent: finding a match in a truth-cropped reference RNA does not constitute genomic discovery. A new search backend, unrestricted duplex scanning, retrained models, and transcript-end inference belong in subsequent experiments.

## Evidence inspected

Reviewed `docs/LITERATURE_BENCHMARK.md`, `docs/IMPROVEMENTS.md`, `benchmarks/results/failure_analysis.md`, `failure_diagnosis.json`, `metric.json`, retained alignment/interaction examples, and the current configuration, search, evidence, model, and orchestration code.

- The completed baseline has 79/79 successful runs, 68/78 localization, 53/78 top-1, 63/78 top-5, 64/78 at least 90% reference coverage, and 3,395 raw candidate hypotheses. The known model-source overlap means this is a regression collection.
- All ten misses have accepted arrays; masking does not erase their reference coordinates. Eight have native genomic FASTA36 alignments at truth that fail the 60% query-span gate. These native spans include 15/47 for Cca1, 16/46 for Cme1, 16/36 for Cca2, 20/46 for Ghy2, 23/47 for Orh, 24/49 for Psp, and 26/47 and 26/46 for Cga and Ffr. The shorter isolated-RNA summary is not interchangeable with the native genomic alignment.
- Esp2 and Bbo have no corresponding native FASTA36 alignment. Esp2's 10-base isolated match is insufficient evidence to broaden genome-wide seed discovery indiscriminately. Bbo also has an isolated 26/38 alignment at 69.23% identity; adjusting a threshold to include that single alignment is not justified.
- Replacing inferred repeats with published repeats does not resolve the misses under the present gates. Thus repeat-consensus improvement is not the immediate dominant bottleneck in this set.
- Four misses have localizing genomic CM hits above the rescue threshold: Cca1 0.0045, Cme1 0.012, Orh 0.0038, Psp 0.00077. Six have no localizing hit in the retained default scan. Increasing the raw genomic CM threshold to 0.1 raises the potential union only to 72/78 while increasing raw CM hits from 69 to 153. Those counts are before final suppression and cannot be described as implemented recall or precision.
- In `parse_alignment_table`, the coverage gate prevents short matches from becoming candidates at all. In `interaction`/`parse_intarna`, the returned minimum-energy interaction can occur anywhere in the extension window. `candidate_evidence` then rewards its energy without requiring paired residues to coincide with the proposed anti-repeat. An interaction near the far end of a window must not validate a different seed near its beginning.
- A CM-only rescue reports 58–61 bases at Lan, Rsp, Cme3, and Ghe. This is valid subregion localization, not full-transcript reconstruction.

A prospective proposal-volume calculation used the fixed rule below on all 79 retained native genomic alignment files, without running a downstream rescue or reading truth to select seeds. It yields **66 additional complementary hypotheses**, 64 unique coordinate/strand anchors, over 30 inputs; median 0, p95 2, maximum 4 per input. The 3,390 standard anti-repeat hypotheses are preserved (the baseline's additional five predictions are genomic model rescues). Per-input counts and raw-file hashes are in [01_proposal_volume.json](01_proposal_volume.json). Esp2 and Bbo each have two extra proposals elsewhere; this calculation does **not** recover their known references or change the finding that no native FASTA36 alignment exists at their truth intervals.

## Biological rationale and limits

Repeat/anti-repeat pairing is imperfect and structural; a short local sequence alignment can seed a longer interaction with mismatches and bulges. The Chyou and Brown study reports interrupted helices in almost all its reference hybrids and uses structural features alongside sequence/context evidence. This supports assessing a duplex after a short seed; it does not validate a universal energy or short-seed threshold. [Chyou and Brown, 2019](https://doi.org/10.1080/15476286.2018.1498281)

The original CRISPRtracrRNA method separates repeat interaction from a distal-region CM, and combines several signals. The currently bundled models should therefore supply complementary regional evidence rather than be treated as models of complete native transcripts. [CRISPRtracrRNA, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9486595/)

The existing IntaRNA adapter can be extended without editing IntaRNA. Its documented `seedTRange` limits where a target seed can occur; `bpList` or validated `hybridDP` provides actual paired positions. A seed constraint preserves full candidate context for accessibility while avoiding a strongest interaction wholly elsewhere in the window. Verify supported syntax against the installed version and a real-tool fixture before depending on it. [Official IntaRNA documentation](https://github.com/BackofenLab/IntaRNA#usage-and-parameters)

## Prioritized implementation contract

### P0: preserve the reference behavior and expose the experiment

Keep current default `similarity=0.7`, `coverage=0.6`, array masking, independent genomic CM rescue, extension lengths, and heuristic weights. Add `short_anti_repeat_rescue=false`; enable it explicitly for a Type II experiment. Reject or clearly disable this profile for the bundled Type V model pool. Custom-model use must record that the supporting CM's family validity is user supplied.

Record configuration, model/input/tool hashes, branch origin, proposal and rejection counts, and every rescue admission reason. Ordinary long anti-repeat candidates must remain in the standard branch when the option is enabled. Avoid changing established candidate IDs merely because extra proposals are present.

### P1: short-seed proposals with anchored duplex support

Recommended **first experimental settings**, to be frozen before evaluation:

| Setting | Proposed value | Meaning |
| --- | --- | --- |
| `short_anti_repeat_rescue` | `false` by default | Experimental Type II branch only |
| `short_min_query_span` | 12 nt | Absolute query span, not gapped alignment length |
| `short_min_identity` | 0.80 | Additional proposals must meet this, even if the standard branch permits 0.70 |
| Coverage | below the configured standard gate | No new fractional threshold; the absolute span is the proposal criterion |
| Duplex seed location | within the proposed anti-repeat's local interval | Require an actual interaction at the proposal |
| Duplex admission | at least 10 paired target residues within the anti-repeat and at least 50% of its target span paired there; total energy < 0 | Structural plausibility screen, not probability |
| Additional regional evidence | correctly oriented CM at E <= 0.01 **or** downstream forward ERPIN terminator at E <= 0.001 | Must occur in the same oriented candidate window |
| Genomic CM-only rescue | unchanged E <= 1e-5 | Weak CM evidence cannot independently enter through this experiment |
| Extension | current 5 nt upstream / 150 nt downstream | Still a search window, not a transcript |

These round values are an initial engineering hypothesis, not literature-derived calibration and not an optimum selected from the ten misses. In particular, 12 nt deliberately permits proposals outside the observed 15–26 nt native examples. The independent support gates may reject some or all of the ten misses, and that result must be retained rather than followed by same-panel cutoff tuning.

For the short branch, retain both tested repeat orientations and require the complementary transcript hypothesis for admission. Do not hard-gate on predicted array strand, Cas distance, coding overlap, or a new local search radius: these could introduce new family-specific sensitivity losses without evaluation. Record those context features. Continue excluding accepted array overlaps and preserving original sequence coordinates.

Validate the duplex encoding: paired counts on both molecules agree, local offsets are within bounds, paired positions refer to canonical or allowed wobble pairs, and malformed output fails explicitly. Prefer `bpList` for coordinates if the installed version supports the documented contract. Counting a broad interaction span as paired bases is insufficient. A minimum-energy interaction entirely downstream cannot admit a short seed. If no anchored interaction is reported, reject and record the reason.

The exact `hybridDP` contract is: one `&`; left part contains only `(` and `.`, right only `)` and `.`; left and right string lengths equal `end1-start1+1` and `end2-start2+1`, respectively; equal nonzero paired-residue counts. Both intervals must have positive ordered coordinates within their respective input sequence. Enumerate left paired positions in ascending target order and right paired positions in descending query order to reconstruct the antiparallel base pairs. With T treated as U, only AU, UA, GC, CG, GU and UG may be paired; ambiguous residues cannot support an admitted pair. Count target paired positions intersecting `candidate.anti_local`, using consistent zero-based half-open intervals after conversion. The two support inequalities are `anchored_paired_count >= 10` and `anchored_paired_count / anti_local.length >= 0.5`. Validate **every CSV row**, including rows that lose energy selection, unknown target/query IDs, nonfinite energies, out-of-bounds coordinates, and duplicated or inconsistent output. A valid no-hit CSV remains a distinct empty result. If `bpList` is also emitted, its pair set must agree with the dot-bracket reconstruction.

For batching, the standard branch can group identical query repeats when all IntaRNA options agree. The anchored short branch must additionally group by identical `seedTRange`/other interaction constraints, or retain individual calls. Applying one candidate's anti-repeat range to another same-repeat target would change the scientific search. Use isolated output paths and preserve an unambiguous target-ID map.

The additional CM must support the downstream region, not merely a wrong-strand or unrelated upstream hit. Allow less than 10 nt overlap with the duplex and require the CM to start at or after the duplex start. A terminator must be on the candidate's forward strand, start downstream of the duplex end, and lie within the candidate window. The regional CM E-value must be labeled with its actual candidate-window search context: it is not numerically interchangeable with a genomic scan E-value.

### P2: interpretable evidence classes and consolidation

Use descriptive classes such as `sequence_only`, `sequence_duplex`, `sequence_duplex_cm`, `sequence_duplex_terminator`, `sequence_duplex_cm_terminator`, and `cm_only`. Add `admission_branch=standard|short_rescue|genomic_cm_rescue`. Report a missing interaction or terminator as missing evidence, not as experimental evidence of absence. Keep `calibrated_probability=false`; avoid labels like 95% confidence or validated.

Keep raw hypotheses intact. Add a separate same-contig, same-strand locus view to reduce duplicated orientations and windows. Prefer an evidence anchor (anti-repeat/paired region or genomic CM interval) rather than broad 150-base window overlap. Group exact duplicate anchors first; a broader rule can merge anchors with at least 50% reciprocal overlap, anchored to a representative rather than transitive chaining. Preserve alternative repeat IDs and all member IDs. Never combine opposite strands or represent a large union as a newly inferred transcript. Report raw and consolidated burdens separately. Changing the main ranking or replacing raw output with this locus view is not required in the first implementation.

### P3: defer boundary and family expansion

Maintain `boundary_kind=anti_repeat_search_window` for all seed branches and `model_supported_region` for CM-only rescues. Additional coverage created by a long window is not accurate endpoint inference. Do not extend the four partial CM-only references just to pass the 90% metric. Terminators can support a proposed downstream region; native versus mature endpoints require separate evidence.

Describe the current `V` models as having demonstrated **Type V-K** scope. The 0/14 Cas12f recognition result is a family-scope limit, not a universal absence of tracrRNA or an instruction to relabel models. Future family models need independent training/selection/test partitions and native/engineered sequence distinctions.

## Fixed evaluation protocol

Before any implementation run, freeze the agreed profile, baseline source hashes, metric file, and control-generation script in a council protocol. The parent council has frozen the 79 reports, source, and metrics in `benchmarks/council_round2/baseline/`. Any deviation must be written before inspecting the affected result. Do not search a parameter grid on this regression panel.

1. Reproduce unchanged-default results for all 79 inputs, using the fixed 78 exactly mapped references. Ghc2 remains unassessed for localization. Preserve execution failures in the denominator. Cached upstream arrays/Cas calls are acceptable only with verified input/configuration/tool identity; all newly affected search, duplex, evidence admission, ranking, and export stages must actually run. Replay claims must be labeled separately from fresh end-to-end runs.
2. Run the frozen short-rescue profile on **all** 79 inputs, not merely the ten misses. Report standard and added candidates, localization/top-1/top-5, >=90% reference coverage, strand errors, per-reference coverage, and the exact admission evidence for every newly recovered reference. Include median/p95/max and per-input candidate counts, unique loci, count/Mb, runtime, and rejected proposal counts. Assess whether correct hits are newly proposed or merely reranked.
3. Run two deterministic conditional null panels over every locus with accepted arrays: (a) dinucleotide-shuffle each accepted repeat independently using seed tags 20260927 and 20260928 while leaving source DNA and context fixed; derive each RNG seed from the tag plus stable source/array/repeat identity (for example SHA-256), so adding or reordering inputs never changes previous decoys; (b) a reproducible cross-locus repeat derangement, excluding identical or reverse-complement-identical repeats and recording every donor/recipient pair. Keep length/order/hash manifests and use isolated output paths for every panel/replicate. Search genome-wide using each decoy repeat and run the same downstream branch. Count **newly admitted short-rescue candidates** per input, array, and Mb, before and after consolidation. Degenerate repeats that cannot produce a distinct shuffle must be flagged, not silently included as valid decoys. Cross-locus pairing can retain legitimate family compatibility; these are background controls, not known false positives.
4. Retain raw CM/shuffle sequence controls as a separate model-only diagnostic; do not reuse their E-value counts to calibrate genomic precision. Optionally add frozen dinucleotide-shuffled array-external genomic segments with retained upstream metadata as a further conditional control, explicitly labeled as an artificial sequence-background test. Do not call disappearance of arrays in a wholly shuffled genome evidence of good rescue specificity.
5. Publish results by evidence class and admission branch. Compare the added rescued-reference count to added candidate burden, and show controls even when the experiment fails. Avoid converting unannotated natural candidates, decoy admissions, or candidate burden into biological precision/FDR.
6. Independent generalization remains unmeasured. A subsequent external test must exclude model-training homologs at clade/RNA-family level, freeze parameter selection before unblinding, retain native-versus-engineered labels, and document annotation completeness. Neither the 303,646 prediction records nor the later date of a paper creates an independent biological truth set.

## Acceptance and review gates

**Required for merge as an experimental feature:** unchanged-default run status and original prediction projection reproduce 79/79 success and the 68/53/63/64 reference metrics; every original prediction remains available under the experimental profile; added candidates never have absent/wrong-strand/unanchored evidence misrepresented as support; seed/duplex/model coordinates reconstruct on both strands; array masks remain honored; configuration and exports expose branch and boundary semantics; all existing tests plus meaningful new adapter/strand/admission tests pass. The branch may be merged disabled even if it does not improve sensitivity, provided that negative outcome is explicit and the option remains experimental.

**Required to call the experiment a useful regression sensitivity improvement:** at least one new mapped reference localized, none of the 68 baseline-localized references lost, all 79 runs complete, and no more than a two-reference decrease in top-1 or top-5. As predeclared operational guardrails, total raw candidate burden must stay <= 1.5 times baseline (<= 5,092), and no single input may add over 100 accepted short-rescue hypotheses. In each conditional null panel the median added candidate count must be <=1 per locus and the maximum <=10. These round guardrails limit triage burden; they are not biological specificity thresholds. Failure means report the result and revise on a separately designated development set, not alter gates until the current panel passes.

**Required for a default change or accuracy claim:** stronger independent validation with suitable negatives, family-aware separation, and explicit treatment of candidate burden and native endpoint uncertainty. Passing the present regression/control protocol alone does not authorize a claim of calibrated confidence, precision, complete-transcript recovery, universal Type V coverage, or an automatic default switch.

I will review the executor's implementation against the final council-agreed settings and protocol, emphasizing seed anchoring, coordinate/strand validity, scope labels, and the distinction between a real genomic discovery and a truth-assisted diagnostic.
