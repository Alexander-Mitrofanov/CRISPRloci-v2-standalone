import publicApiContract from "./contracts/public-api-v1.json";

export const CATEGORY_POLICIES = Object.freeze([...publicApiContract.enums.category_policies]);
export const DEFAULT_CATEGORY_POLICY = CATEGORY_POLICIES[0];

export const ANALYSIS_MODES = [
  {
    id: "locus_compare",
    title: "Compare gene neighborhoods",
    short: "clinker · source-verified proteins and gene order",
    description:
      "Compare 2–8 related genomic neighborhoods, up to 50,000 nt each. " +
      "Inspect all called genes, verified translations and protein similarity alongside Cas evidence. " +
      "Similarity does not establish orthology or evolutionary history.",
    minimumRecords: 2,
    tools: ["CasAndra", "clinker"],
    requiresAdvertisement: true,
  },
  {
    id: "prophage_context",
    title: "Inspect prophage and spacer context",
    short: "Phigaro · candidate regions and source-spacer matches",
    description:
      "Locate candidate prophage regions and inspect exact matches from detected spacers within the submitted sequences. " +
      "Submit one coherent genome/assembly or closely related strains: Prodigal trains jointly across eligible records in single mode, genetic code 11. " +
      "Unrelated genomes can change predictions for unchanged records. " +
      "The 20,000 nt cutoff is eligibility, not a validated accuracy guarantee. Sequence matches and region overlap do not establish active targeting.",
    minimumRecords: 1,
    tools: ["Phigaro", "CRISPRidentify v2"],
    requiresAdvertisement: true,
  },
  {
    id: "insertion_context",
    title: "Inspect insertion-sequence context",
    short: "ISEScan · mobile-element boundary evidence",
    description:
      "Locate insertion-sequence candidates and inspect their relationship to locus organization. " +
      "Partial elements and missing boundary evidence remain explicit; a predicted overlap does not establish gene inactivation.",
    minimumRecords: 1,
    tools: ["ISEScan", "CasAndra"],
    requiresAdvertisement: true,
  },
  {
    id: "defense_context",
    title: "Survey the defense repertoire",
    short: "DefenseFinder · system rules and component evidence",
    description:
      "Inspect modeled defense systems and individual profile hits on the submitted sequences. " +
      "Model release 1.2.3 provides a defined reference scope; predictions do not establish activity or a complete genome-wide repertoire.",
    minimumRecords: 1,
    tools: ["DefenseFinder", "CasAndra"],
    requiresAdvertisement: true,
  },
  {
    id: "rna_family",
    title: "Compare an RNA family",
    short: "MAFFT · RNAalifold consensus structure",
    description:
      "Supply 2–50 homologous RNA or DNA sequences, 30–500 nt each, in matching 5′→3′ orientation. " +
      "Align the family and inspect its predicted consensus fold and descriptive covariation. Homology and transcript boundaries are user supplied.",
    minimumRecords: 2,
    tools: ["MAFFT", "RNAalifold", "ViennaRNA"],
    requiresAdvertisement: true,
  },
  {
    id: "array_compare",
    title: "Compare observed arrays",
    short: "CCTK · ordered spacers and sharing networks",
    description:
      "Reuse detection arrays.json to compare observed spacer order and shared DNA identities. " +
      "Preserve duplicates, deletion slots and unknown orientation. The service record limit counts observed spacer occurrences.",
    minimumRecords: 1,
    tools: ["CCTK"],
    requiresAdvertisement: true,
  },
  {
    id: "spacer_association",
    title: "Associate spacer groups with phages",
    short: "SpacePHARER · translated and nucleotide evidence",
    description:
      "Compare explicitly grouped spacers with an installed phage reference database. " +
      "Combine translated and nucleotide matches; inspect native calibration and occurrence support.",
    minimumRecords: 1,
    tools: ["SpacePHARER"],
    requiresAdvertisement: true,
  },
  {
    id: "repeat_type",
    title: "Predict repeat subtypes",
    short: "RepeatTyper · repeat-derived subtype evidence",
    description:
      "Classify up to 1000 DNA repeats, at most 200 bases each, with a frozen 37-class model. Inspect uncertain and unsupported sequences alongside native model scores.",
    minimumRecords: 1,
    tools: ["RepeatTyper"],
    requiresAdvertisement: true,
  },
  {
    id: "repeat_map",
    title: "Map repeats to reference families",
    short: "CRISPRmap independent historical-reference mapper · annotated reference evidence",
    description:
      "Compare up to 100 DNA repeats, at most 200 bases each, with a frozen annotated reference. " +
      "Recover exact labels and inspect near matches within three global edits, with conflicting or missing labels retained.",
    minimumRecords: 1,
    tools: ["CRISPRmap independent historical-reference mapper"],
    requiresAdvertisement: true,
  },
  {
    id: "viral_search",
    title: "Find viruses for my spacers",
    short: "Spacer sequences against viral RefSeq",
    description:
      "Search supplied CRISPR spacers against viral reference sequences, including phages. Results retain candidate matches and their sequence evidence.",
    minimumRecords: 1,
    tools: ["CRISPRspacer", "NCBI BLAST+"],
    requiresAdvertisement: true,
  },
  {
    id: "loci",
    number: "",
    title: "Annotate a CRISPR locus",
    short: "Arrays, Cas systems, tracrRNA & leader context",
    description:
      "Combine array detection, Cas cassette classification, tracrRNA candidates and leader context windows. Leader prediction is not yet available.",
    minimumRecords: 1,
    tools: ["CRISPRidentify v2", "CasAndra", "CRISPRtracrRNA v3", "CRISPRleader v2"],
  },
  {
    id: "leader",
    number: "",
    title: "Extract leader context",
    short: "CRISPRleader v2 · context only",
    description:
      "Extract source-verified windows on both sides of accepted CRISPR arrays. Leader prediction is unavailable in this development release.",
    minimumRecords: 1,
    tools: ["CRISPRidentify v2", "CRISPRleader v2"],
  },
  {
    id: "cas",
    number: "",
    title: "Find Cas systems",
    short: "CasAndra",
    description:
      "Detect Cas proteins and classify nearby cassettes. Report exact genomic coordinates and model evidence.",
    minimumRecords: 1,
    tools: ["CasAndra"],
  },
  {
    id: "tracrrna",
    number: "",
    title: "Screen tracrRNA models",
    short: "CRISPRtracrRNA v3",
    description:
      "Screen Type II or Type V-K covariance models. These model-supported candidates do not establish complete transcript boundaries.",
    minimumRecords: 1,
    tools: ["CRISPRtracrRNA v3"],
  },
  {
    id: "detection",
    number: "",
    title: "Detect arrays",
    short: "CRISPRidentify v2",
    description:
      "Locate and categorize candidate CRISPR arrays in one or more DNA records. No evolutionary history is inferred.",
    minimumRecords: 1,
    tools: ["CRISPRidentify v2"],
  },
  {
    id: "repeat_context",
    number: "",
    title: "Compare repeats in array context",
    short: "Observed repeats · isolated and local RNA folding",
    description:
      "Detect arrays and compare each observed repeat's predicted pairing in isolation and in array context. Both DNA orientation hypotheses are retained.",
    minimumRecords: 1,
    tools: ["CRISPRidentify v2", "CRISPRrepeat", "ViennaRNA"],
  },
  {
    id: "repeats",
    number: "",
    title: "Analyze repeat sequences",
    short: "Repeat-only DNA or RNA · structural evidence",
    description:
      "Fold supplied repeats and inspect their predicted pair support. DNA retains both orientation hypotheses; " +
      "RNA is supplied in transcribed 5′→3′ order. No array context is inferred.",
    minimumRecords: 1,
    tools: ["CRISPRrepeat", "ViennaRNA"],
  },
  {
    id: "protospacer",
    number: "",
    title: "Find spacers in my sequence",
    short: "Exact matches to reference spacers",
    description:
      "Compare uploaded DNA with the installed 2017 research spacer collection hosted by NCBI. " +
      "Find full-length exact matches in both orientations; matches alone do not identify a host or establish targeting.",
    minimumRecords: 1,
    tools: ["CRISPRspacer"],
    requiresAdvertisement: true,
  },
  {
    id: "reconstruction",
    number: "",
    title: "Evolution & reconstruction",
    short: "SpacerPlacer",
    description:
      "Detect arrays, form comparable groups, and reconstruct ancestral spacer content on an estimated distance tree.",
    minimumRecords: 2,
    tools: ["CRISPRidentify v2", "SpacerPlacer"],
  },
  {
    id: "orientation",
    number: "",
    title: "Orientation-aware evolution",
    short: "CRISPR-evOr",
    badge: "Recommended for related isolates",
    description:
      "Compare input and reversed array order under the evolutionary model, including the supported—or unresolved default—SpacerPlacer reconstruction.",
    minimumRecords: 2,
    tools: ["CRISPRidentify v2", "CRISPR-evOr", "SpacerPlacer"],
  },
];

export const STAGES = [
  { id: "queued", label: "Queued", detail: "Waiting for an analysis worker" },
  { id: "validate_input", label: "Validate input", detail: "Checking records and analysis policy" },
  { id: "cctk", label: "Compare observed arrays", detail: "CCTK · ordered spacer sharing" },
  {
    id: "clinker",
    label: "Compare gene neighborhoods",
    detail: "clinker · verified protein alignments",
  },
  {
    id: "phigaro",
    label: "Inspect prophage context",
    detail: "Phigaro · candidate regions and source matches",
  },
  {
    id: "isescan",
    label: "Find insertion sequences",
    detail: "ISEScan · boundary and transposase evidence",
  },
  {
    id: "defensefinder",
    label: "Survey defense systems",
    detail: "DefenseFinder · pinned model rules",
  },
  { id: "rnaalifold", label: "Compare RNA family", detail: "MAFFT · ViennaRNA consensus folding" },
  { id: "detect_arrays", label: "Detect arrays", detail: "CRISPRidentify v2" },
  { id: "adapt_arrays", label: "Prepare arrays", detail: "Normalizing accepted calls" },
  { id: "preflight_groups", label: "Preflight groups", detail: "Testing evolutionary eligibility" },
  { id: "reconstruct_spacer_histories", label: "Reconstruct histories", detail: "SpacerPlacer" },
  { id: "compare_orientations", label: "Compare orientations", detail: "CRISPR-evOr" },
  { id: "annotate_cas", label: "Annotate Cas systems", detail: "CasAndra" },
  { id: "predict_tracrrna", label: "Predict tracrRNA candidates", detail: "CRISPRtracrRNA v3" },
  {
    id: "extract_leader_context",
    label: "Extract leader context",
    detail: "CRISPRleader v2 · both sides, no prediction",
  },
  {
    id: "repeattyper",
    label: "Classify repeat subtypes",
    detail: "RepeatTyper · preserved model evidence",
  },
  {
    id: "crisprmap",
    label: "Map repeat references",
    detail: "CRISPRmap independent historical-reference mapper",
  },
  { id: "crisprrepeat", label: "Analyze repeat evidence", detail: "CRISPRrepeat · ViennaRNA" },
  {
    id: "crisprspacer",
    label: "Search spacer references",
    detail: "CRISPRspacer · validated sequence matches",
  },
  {
    id: "spacepharer",
    label: "Associate spacer groups",
    detail: "SpacePHARER · target/control calibration",
  },
  { id: "package_results", label: "Package results", detail: "Reports, provenance, and archive" },
];

export const TERMINAL_STATUSES = new Set(publicApiContract.enums.terminal_statuses);

export function stagesForMode(mode) {
  const annotationStages = {
    locus_compare: ["queued", "validate_input", "annotate_cas", "clinker", "package_results"],
    prophage_context: ["queued", "validate_input", "detect_arrays", "phigaro", "package_results"],
    insertion_context: ["queued", "validate_input", "annotate_cas", "isescan", "package_results"],
    defense_context: [
      "queued",
      "validate_input",
      "annotate_cas",
      "defensefinder",
      "package_results",
    ],
    rna_family: ["queued", "validate_input", "rnaalifold", "package_results"],
    array_compare: ["queued", "validate_input", "cctk", "package_results"],
    spacer_association: ["queued", "validate_input", "spacepharer", "package_results"],
    repeat_type: ["queued", "validate_input", "repeattyper", "package_results"],
    repeat_map: ["queued", "validate_input", "crisprmap", "package_results"],
    viral_search: ["queued", "validate_input", "crisprspacer", "package_results"],
    protospacer: ["queued", "validate_input", "crisprspacer", "package_results"],
    repeats: ["queued", "validate_input", "crisprrepeat", "package_results"],
    repeat_context: [
      "queued",
      "validate_input",
      "detect_arrays",
      "crisprrepeat",
      "package_results",
    ],
    leader: [
      "queued",
      "validate_input",
      "detect_arrays",
      "extract_leader_context",
      "package_results",
    ],
    cas: ["queued", "validate_input", "annotate_cas", "package_results"],
    tracrrna: ["queued", "validate_input", "predict_tracrrna", "package_results"],
    loci: [
      "queued",
      "validate_input",
      "detect_arrays",
      "annotate_cas",
      "predict_tracrrna",
      "extract_leader_context",
      "package_results",
    ],
  }[mode];
  if (annotationStages) return STAGES.filter((stage) => annotationStages.includes(stage.id));
  const excluded =
    {
      detection: new Set([
        "adapt_arrays",
        "preflight_groups",
        "reconstruct_spacer_histories",
        "compare_orientations",
      ]),
      reconstruction: new Set(["compare_orientations"]),
      orientation: new Set(),
    }[mode] || new Set();
  return STAGES.filter(
    (stage) =>
      !excluded.has(stage.id) &&
      ![
        "clinker",
        "phigaro",
        "isescan",
        "defensefinder",
        "rnaalifold",
        "cctk",
        "annotate_cas",
        "predict_tracrrna",
        "extract_leader_context",
        "spacepharer",
        "repeattyper",
        "crisprmap",
        "crisprrepeat",
        "crisprspacer",
      ].includes(stage.id),
  );
}

export function analysisModeAvailable(mode, service) {
  const selected = ANALYSIS_MODES.find((item) => item.id === mode);
  return Boolean(
    selected &&
    (!selected.requiresAdvertisement ||
      (service?.state === "online" && service.modes?.includes(mode))),
  );
}

export function orientationLabel(value) {
  const normalized = String(value || "").toLowerCase();
  if (["forward", "input", "input_order", "input-order"].includes(normalized)) {
    return "Input order supported";
  }
  if (["reverse", "reversed", "reverse_order", "reverse-order"].includes(normalized)) {
    return "Reverse input order supported";
  }
  return "Unresolved";
}

export function categoryClass(value) {
  return String(value || "unknown")
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-");
}
