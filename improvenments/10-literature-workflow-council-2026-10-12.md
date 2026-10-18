**Why these five additions: literature and biological workflows**

12 October 2026. Research and planning only; no integration, installation or deployment.

The five candidates remain **clinker, Phigaro, ISEScan, DefenseFinder and RNAalifold**, but the justification and delivery priorities change after reviewing application papers. The purpose is to let a researcher move from a detected CRISPR locus to an inspectable explanation of its organization, sequence matches, surrounding mobile elements, other defenses and RNA conservation.

This is a biological-use follow-up to the [engineering assessment](09-tool-expansion-council-2026-10-12.md). It supersedes that report's ease-first priority and repeat-first RNAalifold recommendation. Its dependency, licensing and containment findings remain relevant. The earlier [connected-pipeline proposal](06-connected-scientific-pipelines.md) supplies the product direction.

The coordinator and two specialist agents reviewed primary application studies and challenged the shortlist. The comparative reviewer emphasized complete locus comparisons and independent organism phylogeny. The locus/RNA reviewer emphasized interrupted genes, incomplete defense predictions and homologous tracrRNA families. We separated three questions: whether a biological question matters, whether researchers demonstrably use this exact tool for it, and whether this repository can support the necessary inputs and interpretation.

This is a focused qualitative literature review, not a systematic review, popularity ranking or user survey. The papers document real workflows; expected community benefits below are our product inferences. We have not measured demand, time saved, accuracy or runtime in CRISPRloci. Exact software versions in a paper establish historical use, not validation of a different current release.

| Proposed addition | Biological question after initial detection | Finished result worth integrating | Main beneficiaries |
| --- | --- | --- | --- |
| **clinker** | Are related arrays associated with conserved, rearranged or partial gene neighborhoods? | Linked comparison of arrays, Cas cassettes and surrounding genes | Comparative microbiologists and evolutionary bioinformaticians |
| **Phigaro** | What genomic context surrounds a spacer's sequence match? | Spacer-to-target evidence joined to candidate prophage intervals and target genes | Microbial ecologists and researchers interpreting CRISPR–mobile-element relationships |
| **ISEScan** | Could an insertion element help explain an unusual or interrupted locus? | A locus-integrity review with IS evidence, source DNA and comparative context | Researchers studying natural variation and apparently incomplete systems |
| **DefenseFinder** | Which other known defense systems occur in the assessed sequence? | A defense-repertoire comparison with component evidence and separate CRISPR annotations | Researchers comparing microbial immune repertoires |
| **RNAalifold** | Do homologous CRISPR-associated RNAs retain a common structural core? | A homologous RNA-family alignment and consensus-structure report | CRISPR RNA biologists and RNA comparative-genomics researchers |

**What researchers actually do to complete these analyses**

The selected papers connect detection to several further stages: extract and curate sequences; align homologues; compare local gene order; search spacer matches; annotate the target's genomic context; interpret the findings against organism or Cas-protein phylogeny; and export figures and evidence for further inspection. Experimental or transcriptomic evidence is then used where a claim concerns activity rather than sequence organization. The individual studies below show which parts each tool supplied.

The implication for CRISPRloci is practical: the useful integration unit is a reproducible analysis with stable links between its evidence. Users should not have to reconstruct which protein came from which CDS, which spacer occurrence produced a match, or which cropped window corresponds to the original contig.

**1. clinker: connect spacer differences to differences in locus architecture**

**Published use.** Siozios et al. used clinker to compare CRISPR-Cas gene order in *Arsenophonus*, alongside spacer-sharing networks. Their analysis also used ISEScan, DefenseFinder, prophage annotation and distinct phylogenetic analyses. Figure 6 connects neighborhood organization with spacer repertoire comparisons during a study of the transition to endosymbiosis. This is a close precedent for combining clinker with the server's existing CCTK analyses. The study's spacer similarity calculation is not identical to CCTK's current Jaccard output. [Current Biology, 2024; Figure 6 and relevant STAR Methods](https://www.sciencedirect.com/science/article/pii/S0960982224014337).

An independent application by Naki and Gophna compared CRISPR-bearing provirus neighborhoods in halophilic archaea using clinker, with separate Cas-protein phylogenies and spacer/repeat analyses. This supports use beyond one bacterial study and highlights the interpretation of partial cassettes. It does not make neighborhood similarity proof of horizontal transfer or functional cooperation. [microLife, 2025; synteny methods and Figure 1](https://academic.oup.com/microlife/article/doi/10.1093/femsml/uqaf033/8297144).

**Why add it.** The current comparison methods can show that arrays share or differ in spacers. A biologist still needs to determine whether those arrays sit beside comparable Cas machinery and surrounding genes. Without that context, an apparent difference in spacer history may be interpreted without knowing that one neighborhood is truncated or structurally different.

**Community benefit — proposed.** A joined comparison would reduce repeated locus extraction, identifier reconciliation and figure preparation. A useful result might show conserved gene order beside divergent spacer repertoires, or a partial cassette beside an otherwise similar neighborhood. The user can inspect that observation before advancing an evolutionary explanation.

**Minimum useful integration.** Select bounded genomic neighborhoods; include all called CDSs, arrays and Cas annotations; retain their original coordinates and verified translations; provide protein-similarity evidence and downloadable sequence/feature files. Let users align these views with CCTK results. CasAndra already retains all-called CDS/protein intermediates, but publishing them requires the independent validation and sanitized export described in the engineering report.

**Limit.** Clinker is a comparison and visualization tool. Its ordering is not an organism tree; protein similarity is not automatically orthology. Retain partial genes, unclassified neighboring proteins and crop/contig boundaries.

**2. Phigaro: turn a spacer hit into an interpretable target-context result**

**Published use.** Roberts et al. identified spacers, searched them against their source genomes, inspected match differences and flanking sequences, then used Phigaro to investigate prophage context in *Ligilactobacillus salivarius*. Selected genomes also received PHASTER analysis. Their definition combined region overlap, proximity and phage-related gene annotation; those are distinct evidence types. The authors did not establish that the CRISPR systems actively targeted or regulated the mobile elements. [mSphere, 2024; “Prediction of self-targeting spacers and prophages” and Discussion](https://journals.asm.org/doi/10.1128/msphere.00171-24).

**Why add it.** A spacer match alone leaves a basic ecological question unanswered: does it lie in a candidate integrated viral region, elsewhere in the genome, or in an uncertain boundary region? Providing a match table without this context pushes users into another annotation service and a manual coordinate merge.

**Community benefit — proposed.** A researcher could inspect the source array, exact match, target gene and candidate prophage region together. This would make candidate self-matches and mobile-element associations easier to review and make the distinction between an observed sequence relationship and a functional interpretation visible.

**Minimum useful integration.** Connect source-bound detected spacers to a search of the source genome or explicitly supplied target; validate the resulting intervals; annotate the appropriate target genome with Phigaro; report containment, boundary overlap and distance separately. Flag or exclude matches inside source CRISPR arrays when interpreting same-genome matches. Preserve strand, mismatch and source-identity evidence.

This requires exposing an existing native CRISPRspacer capability through the web workflow. Today's public protospacer mode searches the uploaded target using historical reference spacers. It cannot, by itself, establish that the submitted genome contains the matching spacer. Merely adding a prophage track would leave the central literature-supported workflow incomplete.

**Limit.** Phigaro predicts candidate regions, not inducibility or ongoing CRISPR activity. Proximity is not containment; a study-specific proximity cutoff should not silently become a biological default. The choice is motivated by documented use and practical local integration, not a claim that Phigaro is universally the most accurate prophage predictor.

**3. ISEScan: help distinguish locus presence from evidence of intact organization**

**Published use.** Sheng et al. used ISfinder searches and manual boundary evidence, then ISEScan, to investigate natural insertion sequences within Cas-encoding loci. Their expanded computational survey reported 163 candidate natural events across 57 genera; this is not 163 experimentally validated losses of activity. The paper separately tested the biological mechanism in selected contexts. [Nature Communications, 2023; “Naturally occurring transpositions of ISs into Cas-encoding genes”](https://www.nature.com/articles/s41467-023-39964-7).

**Why add it.** A detected array and Cas-like proteins do not establish an intact system. An insertion element may provide an explanation for an unusual neighborhood or split gene. That matters when users interpret retained arrays, incomplete cassettes or differences between related genomes.

**Community benefit — proposed.** This makes a common interpretation problem inspectable: “the annotation says this locus exists, but why does its organization differ?” Users gain a connected review of the unusual locus and its sequence evidence instead of a separate inventory of mobile elements.

**Minimum useful integration.** Place ISEScan intervals and available boundary evidence alongside Cas, array and RNA annotations; offer a corresponding neighborhood comparison through clinker; preserve source DNA around the suspected interruption. Report observed overlaps, nearby features and partial calls separately. Terminal inverted repeats and target-site duplications are different evidence; do not imply ISEScan automatically supplies both.

The interface must also expose intergenic sequence and split/partial CDS evidence. An interrupted gene can be split or missed during gene calling, so overlap with an already predicted Cas CDS is not a sufficient discovery strategy. Systematic detection of interrupted homologues would require additional DNA-level comparative analysis; an initial ISEScan adapter must not claim to provide that automatically.

**Limit and integration gate.** An IS prediction or overlap is not a verdict of functional inactivation. Retain this candidate despite its lower engineering convenience: the biological use is particularly direct. The FragGeneScan version discrepancy documented in the engineering report still requires a reproducible dependency pilot before implementation is committed.

**4. DefenseFinder: interpret CRISPR as part of a broader immune repertoire**

**Published use.** Leprince et al. combined DefenseFinder and PADLOC on 263 *Streptococcus thermophilus* genomes. Each tool missed systems called by the other, including examples for which both had models. Their joint analysis illustrates both the value of repertoire annotation and the danger of treating one detector as exhaustive. [Nature Communications, 2025; “S. thermophilus defensome expands beyond CRISPR-Cas and RM systems” and “Defense system prediction”](https://www.nature.com/articles/s41467-025-62408-3).

Maestri et al. combined protein searches, neighborhood examination and defense annotations in their MADS study. At discovery, DefenseFinder and PADLOC did not identify the complete system, while PADLOC's isolated component hits remained informative. This is a useful counterexample: unclassified neighboring genes should remain visible. Functional interaction required further evidence. [Cell Host & Microbe, 2024; Results and STAR Methods](https://doi.org/10.1016/j.chom.2024.07.005).

**Why add it.** A researcher comparing two CRISPR loci may also need to know how the other predicted defense systems differ. CRISPR-only annotation cannot answer that question, and CRISPR absence should not be presented as absence of antiviral defenses.

**Community benefit — proposed.** A reproducible repertoire view would support comparative immune-system studies and help users identify which observations require further validation. Stable exports would also make external tool comparison possible without forcing researchers to reconcile protein names and genomic positions again.

**Minimum useful integration.** Use the bounded CPU HMM/system-rule workflow with frozen software and model releases. Scan complete ordered protein evidence, preserve contig boundaries, and map calls back to verified CDSs. Display system calls, component genes and isolated profile hits separately. Retain CasAndra and DefenseFinder labels as distinct evidence rather than silently replacing or merging classifications.

**Limit.** “No system detected in the assessed sequence” is the appropriate negative result. A short locus window is not a complete host immune repertoire. Co-occurrence or proximity does not establish synergistic immunity. The papers support making independent results interoperable; they do not establish that adding a second caller would eliminate missed systems. PADLOC comparison/export is a valuable later extension rather than a sixth required integration in this proposal.

**5. RNAalifold: compare homologous RNA families, especially tracrRNAs**

**Published use.** Gasiunas et al. recruited tracrRNA homologues, grouped sequences, aligned them with MAFFT and used RNAalifold to create group structure models; Infernal supplied subsequent sequence/structure comparisons. This was a comparative RNA analysis within a broader Cas9 study, not simply another fold of one repeat. [Nature Communications, 2020; “Computational analysis of Cas9 tracrRNAs”](https://www.nature.com/articles/s41467-020-19344-1).

Faure et al. also combined genomic context, homology searches and RNA structure methods. Their RNAalifold analysis included aligned nonhybridizing tracrRNA regions, with other methods used for repeat–anti-repeat duplexes. This distinction matters: a single-sequence-family consensus fold does not directly model an arbitrary two-strand RNA interaction. [RNA Biology, online 2018 / 2019 issue; comparative methods](https://doi.org/10.1080/15476286.2018.1493331).

Liu et al. independently used MAFFT, RNAalifold and Infernal for comparative tracrRNA analysis in a type V-K study, extending the exact-tool precedent beyond type II. Candidate activity was assessed separately. This is relevant to the server's existing type V-K tracrRNA branch. [Nature Communications, 2025; computational discovery methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11906591/).

**Why add it.** CRISPRtracrRNA identifies RNA candidates and existing repeat analysis describes individual sequences. Researchers may then ask whether homologous candidates share structural features despite sequence variation, or whether an apparent difference comes from inconsistent boundaries or orientation.

**Community benefit — proposed.** A reproducible RNA-family report would carry candidate provenance through alignment and structure comparison, giving RNA biologists an inspectable starting point for curation and further analysis. It would replace manual sequence extraction and detached folds with linked evidence.

**Minimum useful integration.** Let users assemble or import a justified homologous family, review strand and candidate boundaries, and export the alignment, original sequence crosswalk and consensus pairing evidence. Use the existing ViennaRNA engine through an explicit versioned adapter and a controlled alignment step. Preserve sequence redundancy and weighting decisions in the report. Do not group RNAs solely because the loci share a subtype.

**What changes from the previous proposal.** A repeats-only implementation may be an inexpensive pilot, but it does not fulfill the strongest use case found here. Prioritize a homologous tracrRNA-family workflow, with repeat-family comparison as another eligible branch. Existing candidate sequence boundaries are not necessarily mature transcript boundaries. Family selection and validated RNA export are substantive work even though the folding dependency already exists.

**Limit.** A consensus structure does not prove transcript processing, functional compatibility or statistically significant covariation. An alignment with little independent variation should show that limitation. Optional later covariation testing is useful, but it should not be confused with RNAalifold's own scoring or made a prerequisite merely because another tool is available.

**Council decisions and delivery order**

The recommendation is to retain all five with the narrower commitments above. ISEScan moves up in biological relevance; RNAalifold moves away from being an ease-first addition. Two discussions shaped the result:

- **RNA families versus organism phylogeny.** Independent host phylogeny is important when interpreting evolutionary relationships. However, a similarity screen such as skani is not a replacement for that phylogeny, and a new core-genome inference workflow is a larger undertaking. Retain RNAalifold for its direct RNA-biological use; add an optional imported Newick tree and sample metadata, with explicit leaf-to-source mappings, inference-method provenance and branch-support information, as complementary comparison support. Never label a spacer history or similarity network as the host tree. If only a repeats-only folding wrapper is feasible, prioritize host-context support first.
- **Biological priority versus integration ease.** An ISEScan adapter has a stronger specific integrity-review rationale than its engineering rank alone suggested. Nevertheless, its dependency pilot remains a gate. Likewise, an easy RNAalifold call without family curation would not complete the published workflow.

The recommended first release unit is **locus comparison with clinker**, built on a shared, source-verified CDS/protein export. Develop **spacer-to-target context with Phigaro** alongside the missing public source-spacer search connection. Run the **ISEScan feasibility pilot** early, then add its evidence to the locus view. **DefenseFinder** can reuse the same complete protein/coordinate foundation once assessed-genome scope is clear. Deliver **RNAalifold** when the homologous-family and RNA-boundary contracts are ready. These are dependency-based priorities, not estimated implementation durations.

Three shared improvements are necessary for the community benefit to materialize:

1. **Stable evidence joins and usable exports.** Carry sample, source sequence, locus, CDS, RNA and spacer identities through all results. Export sequences, feature tables, alignments, coordinates, parameters and tool/model versions. Treat this as a scientific feature.
2. **Explicit analysis scope.** The current default total input limit is 2 Mb. Many complete bacterial genomes exceed it. A whole-genome repertoire or prophage workflow needs measured capacity changes, assembly/sample grouping and explicit handling of fragmentary inputs. Do not infer organism-wide absence from an uploaded fragment.
3. **Results organized by biological question.** Present locus organization, target context, integrity evidence, immune repertoire and RNA conservation as connected views. Keep detected features, model predictions and biological interpretations distinguishable.

A future integration should be judged against reference examples with expected coordinate joins and known limitations, not merely whether its executable returns success. Useful acceptance questions are: can a user reconstruct the source of every displayed feature; compare the relevant loci; export the evidence; and identify where the result stops short of a functional conclusion? Scientific sensitivity, host coverage and runtime still require separate evaluation.

**Evidence access and remaining uncertainty**

The review checked relevant primary Methods, Results and figure text rather than relying only on software descriptions. For Siozios et al., the publisher's relevant methods passages and Figure 6 caption were retrieved through expanded indexed extraction; direct full-article opening failed. That is a narrower access claim than reading its uninterrupted full text. Publisher/PMC text or a primary-author repository copy was available for the other cited application evidence. We excluded unverified workflow claims and did not use citation counts as evidence of community preference.

The next evidence gap is user experience and local scientific evaluation: whether the proposed joins save researchers effort, whether the selected versions reproduce appropriate published examples, and how predictions behave on the server's supported organisms and incomplete assemblies. Those questions are not settled by literature examples or council agreement.
