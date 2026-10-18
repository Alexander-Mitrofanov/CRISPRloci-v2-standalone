# Literature and scientific interpretation register

Targeted primary-literature review, checked 27 September 2026. This is not an
exhaustive systematic review. Publications establish precedent and constraints;
the proposed v2 designs and release gates are council recommendations. No
published accuracy result is being transferred to v2 without a local benchmark.
Author documentation establishes software availability, not biological validity.

## Historical scope and existing methods

**E1 — CRISPRloci baseline.** Alkhnbashi et al. (2021), NAR 49(W1):W125–W130,
[CRISPRloci](https://doi.org/10.1093/nar/gkab456).
The [supplementary ZIP](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8265192/supplementaryFiles)
was recovered and Tables S1/S2 visually inspected. Use the historical column in
the parity audit as the compact capability inventory. The
[official help](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=CRISPRloci)
and [original repository](https://github.com/BackofenLab/CRISPRloci) corroborate
operational details. They are living resources, not frozen publication snapshots.

**E2 — Detector evidence.** Mitrofanov et al. (2021; online 2020), NAR 49:e20,
[CRISPRidentify](https://doi.org/10.1093/nar/gkaa1158).
Candidate generation, feature extraction and learned classification support
retaining categories and evidence, including atypical arrays. This paper is
scientific lineage for the existing detector; it does not validate all changes
in the pinned v2 code. Evaluate its current release on independent fixtures.

## Orientation, leaders and repeats

**E3 — Conserved leader prediction.** Alkhnbashi et al. (2016), Bioinformatics
32:i576–i585, [Characterizing leader sequences of CRISPR loci](https://doi.org/10.1093/bioinformatics/btw454).
CRISPRleader supplies precedent for conserved-core boundary prediction.
Heuristic/conservation exclusions may miss divergent functional leaders.
**Consequence:** distinguish detected core, insufficient flank, unsupported
family and no detection. A fixed upstream window is context, not this inference.

**E4 — Transcribed-strand evidence.** Alkhnbashi et al. (2014), Bioinformatics
30:i489–i496, [CRISPRstrand](https://doi.org/10.1093/bioinformatics/btu459).
Repeat-based strand inference supports a transcription-oriented annotation.
**Consequence:** preserve that evidence separately from evOr's group-level
acquisition-order comparison; agreement of correlated methods is not a
calibrated probability.

**E5 — Repeat families.** Lange et al. (2013), NAR 41:8034–8044,
[CRISPRmap](https://doi.org/10.1093/nar/gkt606).
Comparative repeat sequence/structure families have subtype-dependent
relationships. **Consequence:** retain vocabulary and reference versions;
family assignment, similarity retrieval and Cas-subtype inference are different
outputs. Test unfamiliar families and allow unassigned results.

**E6 — Thermodynamic tools.** Lorenz et al. (2011), Algorithms for Molecular
Biology 6:26, [ViennaRNA Package 2.0](https://doi.org/10.1186/1748-7188-6-26).
MFE and ensemble/local pairing predictions motivate RNAfold/RNAplfold evaluation.
**Consequence:** retain software, temperature, energy parameters, sequence and
context; present predicted support, not experimentally measured structure.

**E7 — Experimental counterexample to a hairpin requirement.** Shao et al.
(2016), Structure 24:547–554,
[A Non-Stem-Loop CRISPR RNA Is Processed by Dual Binding Cas6](https://doi.org/10.1016/j.str.2016.02.009).
Structural/biochemical work supports processing of a non-stem-loop repeat in a
protein-stabilized conformation. **Consequence:** poor isolated-repeat folding
must not exclude an otherwise supported array.

**E8 — Repeat-subtype comparator.** Russel et al. (2020), The CRISPR Journal,
[CRISPRCasTyper](https://doi.org/10.1089/crispr.2020.0059).
Repeat-based subtyping provides a comparator for repeat-only and orphan-array
evaluation. **Consequence:** assess family-held-out performance and abstention;
this is not a decision to add a second full locus pipeline beside CasAndra.

**E8b — Processing and structural context.** Reimann et al. (2017), NAR
45:915–925, [Structural constraints and enzymatic promiscuity in the
Cas6-dependent generation of crRNAs](https://doi.org/10.1093/nar/gkw786).
Experiments in a specific Cas6 context motivate inspecting repeat-instance
differences. **Consequence:** per-repeat context deserves visibility, but a
folding metric does not establish universal processing efficiency.

## Spacer matches and host association

**E9 — Search-engine benchmark.** Neri et al. (2026), Bioinformatics
42(7):btag394, published 29 June,
[Computational tool choice impacts CRISPR spacer-protospacer detection](https://doi.org/10.1093/bioinformatics/btag394).
Synthetic, semi-synthetic and real comparisons show that tool choice, distance
metric and reporting behavior affect recall and cost. Its practical settings
motivate testing Bowtie1 and alternatives, not adopting a universal mismatch
threshold. Large comparisons lack exhaustive ground truth, and the paper's
chance-hit discussion contains inconsistent normalized ratios.
**Consequence:** use a small exhaustive oracle and independent biological
evaluation; import no universal false-positive percentage. Author
[workflows](https://github.com/UriNeri/spacer_matching_bench) and
[archived runs](https://doi.org/10.5281/zenodo.15171878) provide reproducible
benchmark starting points. Match correctness and host attribution need separate
measurements.

**E10 — Aggregating host evidence.** Zhang et al. (2021), Bioinformatics
37:3364–3366, [SpacePHARER](https://doi.org/10.1093/bioinformatics/btab222).
Translated search and null-database evaluation support an alternative host
association approach. **Consequence:** useful comparator; its association score
is not a nucleotide/PAM-supported interference prediction.

**E11 — Host prediction benchmark.** Roux et al. (2023), PLOS Biology
21:e3002083, [iPHoP](https://doi.org/10.1371/journal.pbio.3002083).
Integrated evidence was evaluated for host-genus prediction, including novel
viruses. **Consequence:** evaluate taxonomic rank, coverage and uncertainty;
reported integrated performance does not transfer to a CRISPR-only feature.

**E12 — Self-targeting associations.** Nobrega et al. (2020), NAR
48:12074–12084, [Prophages are associated with extensive CRISPR–Cas
auto-immunity](https://doi.org/10.1093/nar/gkaa1071).
Genome-scale associations motivate prophage context around candidate self-hits.
**Consequence:** a match and overlap do not establish its mechanism in an
individual locus, an active immune response or anti-CRISPR presence.

**E12b — Experimental dependence on system/context.** Rollie et al. (2020),
Nature 578:149–153, [study of CRISPR–Cas immunity and temperate phages](https://doi.org/10.1038/s41586-020-1936-2),
demonstrated costs of imperfect targeting and anti-CRISPR suppression in the
studied type I system. Keith et al. (2025), Cell Host & Microbe 33:330–340.e6,
[Cas9 and temperate-phage immunity](https://doi.org/10.1016/j.chom.2025.02.001),
demonstrated prophage curing under different experimental conditions.
**Consequence:** preserve context and hypotheses; do not assign the same
biological interpretation to every genomic spacer match.

**E13 — Local prophage annotation candidate.** Wishart et al. (2023), NAR
51:W443–W450, [PHASTEST](https://doi.org/10.1093/nar/gkad382).
The authors describe improved prophage annotation and a local implementation;
[official downloads](https://www.phastest.ca/databases) provide backend and
database assets. **Consequence:** evaluate installation, licensing, runtime,
reference vintage and supported host lineages before selection. Pin database
and executable independently. Validate predicted regions and boundaries,
not only downstream interval overlaps; bacterial results do not establish
archaeal-virus coverage.

## Cas interpretation and conditional extensions

**E14 — Boundary/module precedent.** Padilha et al. (2021; online 2020),
Bioinformatics 37:1352–1359,
[Casboundary](https://doi.org/10.1093/bioinformatics/btaa984).
Learned cassette boundaries and module decomposition supply a historical
comparison target. **Consequence:** evaluate what current CasAndra already
provides and add missing interpretation as a separate evidence layer; no
automatic legacy-tool reinstatement.

**E15 — Missing-component hypotheses.** Padilha et al. (2020), GigaScience
9:giaa062, [CRISPRcasIdentifier](https://doi.org/10.1093/gigascience/giaa062).
Regression of potentially missing protein-family scores supported cassette
classification and follow-up hypotheses. **Consequence:** a predicted score
does not detect a gene. Test masked-component cases and retain observed versus
inferred status. Replacing this feature with observed coverage alone is a
deliberate capability change that remains on the parity ledger.

**E16 — Current nomenclature context.** Makarova et al. (2025), Nature
Microbiology 10:3346–3361,
[An updated evolutionary classification of CRISPR–Cas systems including rare
variants](https://doi.org/10.1038/s41564-025-02180-8).
The analysis describes 2 classes, 7 types and 46 subtypes, including rare and
functionally varied systems. **Consequence:** version nomenclature and model
coverage separately; do not relabel unsupported calls or equate partial/mixed
cassettes with newly discovered effectors. A taxonomy update does not retrain
the existing classifier.

**E17 — Flanking context, with mechanism limits.** Biswas et al. (2013), RNA
Biology 10:817–827, [CRISPRTarget](https://doi.org/10.4161/rna.24046),
provides precedent for target searches and inspected flanks. Meeske and
Marraffini (2018), Molecular Cell 71:791–801.e3,
[type VI-A target-flank self-recognition](https://doi.org/10.1016/j.molcel.2018.07.013),
demonstrates a distinct RNA-target context.
**Consequence:** annotate oriented flanks first; infer candidate motifs only
with independent support/background controls. A single DNA-PAM filter is
inappropriate across all systems, and motif enrichment may reflect acquisition
or survival rather than experimentally measured interference requirements.

## Evidence strength and remaining work

Historical functionality is supported directly by A/S and corroborated by H.
Leaders, family assignment, folding and matching have published computational
methods; that is sufficient to prioritize evaluations, not to promise their
performance on every lineage. E7, E8b and E12b provide experimental reasons to
avoid overly broad interpretation, with scope limited to the tested systems.

Before implementation selection, prepare accession-level benchmark manifests,
independent labels, training/reference overlap checks, licenses and resource
estimates. Literature disagreement or uncertainty must produce an explicit
evaluation task rather than a silently selected default. Current pinned
CasAndra/tracrRNA/Identify revisions are software facts from the repository;
this review does not imply each revision has separate peer-reviewed validation.
