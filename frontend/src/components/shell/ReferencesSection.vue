<script setup>
import { ref } from "vue";
import AppIcon from "../common/AppIcon.vue";

const disclosure = ref(null);
defineExpose({
  open() {
    disclosure.value.open = true;
    disclosure.value.scrollIntoView({ block: "start" });
    disclosure.value.querySelector("summary")?.focus({ preventScroll: true });
  },
});

const citations = [
  {
    tool: "clinker",
    authors: "Gilchrist and Chooi",
    venue: "Bioinformatics · 2021",
    title: "clinker & clustermap.js: automatic generation of gene cluster comparison figures",
    doi: "https://doi.org/10.1093/bioinformatics/btab007",
    source: "https://github.com/gamcil/clinker",
    note: "Protein similarity and gene order support neighborhood comparison; neither establishes orthology or organism phylogeny.",
  },
  {
    tool: "Phigaro",
    authors: "Starikova et al.",
    venue: "Bioinformatics · 2020",
    title: "Phigaro: high-throughput prophage sequence annotation",
    doi: "https://doi.org/10.1093/bioinformatics/btaa250",
    source: "https://github.com/bobeobibo/phigaro",
    note: "Candidate prophage intervals provide context for source-spacer matches. They do not establish viral activity or CRISPR targeting.",
  },
  {
    tool: "ISEScan",
    authors: "Xie and Tang",
    venue: "Bioinformatics · 2017",
    title:
      "ISEScan: automated identification of insertion sequence elements in prokaryotic genomes",
    doi: "https://doi.org/10.1093/bioinformatics/btx433",
    source: "https://github.com/xiezhq/ISEScan",
    note: "Insertion-element and boundary evidence support locus inspection; predicted overlaps do not prove gene inactivation.",
  },
  {
    tool: "DefenseFinder",
    authors: "Tesson et al.",
    venue: "Nature Communications · 2022",
    title: "Systematic and quantitative view of the antiviral arsenal of prokaryotes",
    doi: "https://doi.org/10.1038/s41467-022-30269-9",
    source: "https://github.com/mdmparis/defense-finder",
    note: "The web workflow uses HMM/system rules and historical model release 1.2.3. It does not assess the complete known defensome.",
  },
  {
    tool: "MacSyFinder",
    authors: "Néron et al.",
    venue: "Peer Community Journal · 2023",
    title:
      "MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes",
    doi: "https://doi.org/10.24072/pcjournal.250",
    source: "https://github.com/gem-pasteur/macsyfinder",
    note: "System-rule engine used by DefenseFinder.",
  },
  {
    tool: "RNAalifold",
    authors: "Bernhart et al.",
    venue: "BMC Bioinformatics · 2008",
    title: "RNAalifold: improved consensus structure prediction for RNA alignments",
    doi: "https://doi.org/10.1186/1471-2105-9-474",
    source: "https://github.com/ViennaRNA/ViennaRNA",
    note:
      "Credit: ViennaRNA authors and the Institute for Theoretical Chemistry, University of Vienna. " +
      "Comparative folding of a researcher-defined homologous family; descriptive covariation is not a statistical significance test.",
  },
  {
    tool: "MAFFT",
    authors: "Katoh and Standley",
    venue: "Molecular Biology and Evolution · 2013",
    title:
      "MAFFT multiple sequence alignment software version 7: improvements in performance and usability",
    doi: "https://doi.org/10.1093/molbev/mst010",
    source: "https://mafft.cbrc.jp/alignment/software/",
    note: "Alignment method used before RNA-family consensus folding.",
  },
  {
    tool: "CCTK",
    authors: "Collins et al.",
    venue: "The CRISPR Journal · 2023",
    title: "CRISPR comparison toolkit",
    doi: "https://doi.org/10.1089/crispr.2022.0080",
    note: "Observed spacer order and set-based sharing; comparison is distinct from evolutionary reconstruction.",
    source: "https://github.com/Alan-Collins/CRISPR_comparison_toolkit",
  },
  {
    tool: "SpacePHARER",
    authors: "Zhang et al.",
    venue: "Bioinformatics · 2021",
    title:
      "SpacePHARER: sensitive identification of phages from CRISPR spacers in prokaryotic hosts",
    doi: "https://doi.org/10.1093/bioinformatics/btab222",
    source: "https://github.com/soedinglab/spacepharer",
    note: "Group–target association evidence from translated and nucleotide searches, with empirical target/control calibration.",
  },
  {
    tool: "CRISPRidentify",
    authors: "Mitrofanov et al.",
    venue: "Nucleic Acids Research · 2021",
    title: "CRISPRidentify: identification of CRISPR arrays using machine learning approach",
    doi: "https://doi.org/10.1093/nar/gkaa1158",
    source: "https://github.com/BackofenLab/CRISPRidentify",
  },
  {
    tool: "CasAndra",
    authors: "Makarova et al.",
    venue: "Nature Microbiology · 2025",
    title: "An updated evolutionary classification of CRISPR–Cas systems including rare variants",
    doi: "https://doi.org/10.1038/s41564-025-02180-8",
    note: "CRISPR–Cas classification reference used by CasAndra.",
  },
  {
    tool: "CRISPRtracrRNA",
    authors: "Mitrofanov et al.",
    venue: "Bioinformatics · 2022",
    title: "CRISPRtracrRNA: robust approach for CRISPR tracrRNA detection",
    doi: "https://doi.org/10.1093/bioinformatics/btac466",
    source: "https://github.com/BackofenLab/CRISPRtracrRNA",
  },
  {
    tool: "CRISPRleader",
    authors: "Alkhnbashi et al.",
    venue: "Bioinformatics · 2016",
    title: "Characterizing leader sequences of CRISPR loci",
    doi: "https://doi.org/10.1093/bioinformatics/btw454",
    source: "https://github.com/BackofenLab/CRISPRleader",
  },
  {
    tool: "CRISPRrepeat",
    authors: "Lorenz et al.",
    venue: "Algorithms for Molecular Biology · 2011",
    title: "ViennaRNA Package 2.0",
    doi: "https://doi.org/10.1186/1748-7188-6-26",
    source: "https://github.com/ViennaRNA/ViennaRNA",
    sourceLabel: "ViennaRNA source",
    note: "Folding engine used by CRISPRrepeat and the RNAalifold family workflow.",
  },
  {
    tool: "RepeatTyper / CRISPRCasTyper",
    authors: "Russel et al.",
    venue: "The CRISPR Journal · 2020",
    title:
      "CRISPRCasTyper: Automated Identification, Annotation, and Classification of CRISPR-Cas Loci",
    doi: "https://doi.org/10.1089/crispr.2020.0059",
    source: "https://github.com/Alexander-Mitrofanov/RepeatTyper",
    note: "Standalone repeat classifier with the preserved 37-label model. Closed-set scores are not calibrated confidence.",
  },
  {
    tool: "CRISPRmap",
    authors: "Lange et al.",
    venue: "Nucleic Acids Research · 2013",
    title:
      "CRISPRmap: an automated classification of repeat conservation in prokaryotic adaptive immune systems",
    doi: "https://doi.org/10.1093/nar/gkt606",
    source: "https://github.com/Alexander-Mitrofanov/CRISPRmap-v2",
    sourceLabel: "Independent mapper source",
    note: "The paper describes the original CRISPRmap. This independent historical-reference mapper does not recreate its classification pipeline.",
  },
  {
    tool: "CRISPRspacer",
    authors: "Alkhnbashi et al.",
    venue: "Nucleic Acids Research · 2021",
    title: "CRISPRloci: comprehensive and accurate annotation of CRISPR–Cas systems",
    doi: "https://doi.org/10.1093/nar/gkab456",
    source: "https://github.com/BackofenLab/CRISPRloci",
    note: "Original web-server and spacer-analysis reference.",
  },
  {
    tool: "BLAST+",
    authors: "Camacho et al.",
    venue: "BMC Bioinformatics · 2009",
    title: "BLAST+: architecture and applications",
    doi: "https://doi.org/10.1186/1471-2105-10-421",
    note: "Sequence search method used by CRISPRspacer and CRISPRtracrRNA.",
  },
  {
    tool: "SpacerPlacer",
    authors: "Fehrenbach et al.",
    venue: "Nucleic Acids Research · 2024",
    title:
      "SpacerPlacer: ancestral reconstruction of CRISPR arrays reveals the evolutionary dynamics of spacer deletions",
    doi: "https://doi.org/10.1093/nar/gkae772",
    source: "https://github.com/fbaumdicker/SpacerPlacer",
  },
  {
    tool: "CRISPR-evOr",
    authors: "Fehrenbach et al.",
    venue: "PLOS Computational Biology · 2025",
    title: "An evolutionary approach to predict the orientation of CRISPR arrays",
    doi: "https://doi.org/10.1371/journal.pcbi.1013706",
    source: "https://github.com/fbaumdicker/SpacerPlacer",
  },
];
</script>

<template>
  <section
    id="references"
    class="references"
    aria-labelledby="references-heading"
  >
    <details ref="disclosure">
      <summary id="references-heading">References</summary>
      <div class="reference-grid">
        <article
          v-for="item in citations"
          :key="item.tool"
        >
          <span>{{ item.authors }} · {{ item.venue }}</span>
          <h3>{{ item.tool }}</h3>
          <p>{{ item.title }}</p>
          <p v-if="item.note">{{ item.note }}</p>
          <div>
            <a
              :href="item.doi"
              target="_blank"
              rel="noopener noreferrer"
              >Publication
              <AppIcon
                name="external"
                :size="15" /></a
            ><a
              v-if="item.source"
              :href="item.source"
              target="_blank"
              rel="noopener noreferrer"
              >{{ item.sourceLabel || "Original source" }}
              <AppIcon
                name="external"
                :size="15"
            /></a>
          </div>
        </article>
      </div>
    </details>
  </section>
</template>
