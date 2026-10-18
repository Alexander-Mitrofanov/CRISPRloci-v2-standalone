# Comparable tasks and literature audit — 27 September 2026

This is a targeted primary-literature audit, not a proof that every existing tool
has been found. Recheck the inventory before freezing an actual comparison.
No competitor was run in this audit. Article results from different datasets and
metrics are not a leaderboard for v2.

| Method / primary source | Published task | Fair comparison with leader v2 | Important limitation |
| --- | --- | --- | --- |
| [CRISPRleader, 2016](https://doi.org/10.1093/bioinformatics/btw454); [author code](https://github.com/BackofenLab/CRISPRleader) | Conserved leader-core localization with repeat-related reference families; array orientation in the pipeline | Main direct comparator for conserved-core boundaries, with identical input arrays where an interface permits | A core boundary is not a full functional leader boundary. Freeze executable and model identities; legacy reuse terms remain unresolved. |
| [CRISPRstrand, 2014](https://doi.org/10.1093/bioinformatics/btu459) | Repeat-based crRNA-encoding strand prediction | Transcription-orientation track with the same arrays and repeat evidence | Does not establish leader-core endpoints. Its reported ROC AUC is not directly comparable with another study's accuracy. |
| [CRISPRDirection, 2014](https://doi.org/10.1093/bioinformatics/btu114) | Transcription direction using repeat/flanking evidence and weighted voting | Orientation track; preserve its no-prediction outcomes | It is not a conserved-core boundary predictor. Reference labels include expansion from a smaller experimentally oriented repeat set; audit overlap and evidence level. |
| [CRISPRDetect, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4869251/); [CRISPRCasFinder, 2018](https://doi.org/10.1093/nar/gky425) | Array detection/annotation including orientation outputs | Deployment tracks for full pipelines, separately from fixed-array comparisons | Shared CRISPRDirection ancestry does not make these independent algorithms. Array detection failures cannot disappear from pipeline denominators. |
| Potential Orientation in [CRISPRCasFinder](https://doi.org/10.1093/nar/gky425), assessed in the [2019 comparison](https://doi.org/10.3389/fmicb.2019.02054) | A leader-context orientation heuristic | Explicit context baseline, separately configured from the full voting method | Shared component ancestry and incomplete applicability; no independent core-boundary claim. |
| [Cas Orientation, 2019](https://doi.org/10.3389/fmicb.2019.02054) | Direction inferred from associated cas-gene context | Context-enabled orientation baseline on eligible loci | Applicability depends on nearby cas context and subtype; not a universal rule for orphan or bidirectional arrays. |
| [Leader motif annotation, Buyukyoruk et al., 2025](https://doi.org/10.1089/crispr.2024.0093) | Taxon/subtype-associated motif evidence and annotation | Motif/subtype evidence track using equivalent labels and annotation inputs | The study found no universal leader motif. Additional annotations are not, by themselves, measured sensitivity or boundary accuracy. |
| [CRISPR-evOr, 2025](https://doi.org/10.1371/journal.pcbi.1013706) | Acquisition orientation from related spacer-array groups | A separate acquisition-orientation / discordance track | Acquisition and transcription orientation are distinct targets. Group-based evidence needs related arrays; do not rank it as a single-genome core-boundary predictor. |
| [CRISPRloci, 2021](https://doi.org/10.1093/nar/gkab456) | Integrated array, leader and cas annotation | Deployment workflow and historical integration comparison | An integrated wrapper sharing underlying methods is not an independent scientific baseline. |

Identify, CasAndra, CRISPRtracrRNA, evOr and SpacerPlacer already belong to the
server. This project must not reimplement them to inflate its feature count.
Identify's existing orientation can be an explicit baseline, and evOr an optional
external evidence source with its own orientation semantics. CasAndra context is
an input only in a separately declared context-enabled track.

## Required run registry

Before final data are evaluated, record each included implementation's repository,
commit/release, executable/container digest, reference hashes, license/availability,
configuration, orientation convention, supported output meanings and input budget.
Distinguish algorithm families from pipeline implementations. Registry entries in
`comparators.json` are candidate entries, not frozen runnable environments.

Use two separate comparison tracks:

1. **Released-tool deployment:** each tool with its documented reference resources;
   differences in resource coverage are part of measured practical performance.
2. **Controlled-method:** equal input information, development/reference data and
   tuning budgets where compatible retraining is actually supported. Otherwise
   mark this comparison unavailable; do not silently substitute a reimplementation.

Fixed windows, simple sequence matching, motif-only scoring and component removal
are internal baselines/ablations. They cannot replace a head-to-head test against
the original tool or another eligible independently developed method. A tool that
cannot be installed, licensed or supplied compatible inputs remains untested;
it cannot be counted as a defeated competitor.

## Search record

Search date: 2026-09-27. Queries included CRISPR leader prediction tools,
CRISPRleader/CRISPRDirection/CRISPRstrand benchmarks, leader prediction in
2024–2026, Cas Orientation, and the 2025 motif/evOr studies. Author publications
and repositories were used for the task mappings above. General reviews were
used only for discovery, not as evidence of individual method performance.
Some full-text retrievals were blocked by publisher errors or CAPTCHA; searchable
primary text, abstracts and the accessible PLOS article supplied the scope evidence.
No biological reference sequences, model files or assay optimization procedures
were acquired for this audit.
