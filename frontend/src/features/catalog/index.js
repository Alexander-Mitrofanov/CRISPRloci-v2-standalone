const location = ["accession", "Genome", "genome"];
export { sequenceSearchIssue, sequenceComparison } from "./sequenceSearch.js";
export { useSequenceSearch } from "./useSequenceSearch.js";
const interval = ["interval", "Coordinates"];
const sequence = ["sequence", "Sequence", "sequence"];
const tool = ["tool", "Caller"];

export const CATALOG_ENTITIES = [
  {
    id: "assemblies",
    label: "Genomes",
    description:
      "Assembly records retain their analysis status. Failed analyses are not negative biological results.",
    columns: [
      location,
      ["status", "Analysis status"],
      ["sequence_count", "Source contigs"],
      ["array_count", "Arrays"],
      ["cas_gene_count", "Cas genes"],
      ["cas_system_count", "Cas systems"],
    ],
  },
  {
    id: "arrays",
    label: "Arrays",
    description:
      "Predicted arrays with repeat consensus and caller evidence. Open an array to inspect its repeat and spacer occurrences.",
    columns: [
      ["call_id", "Array", "array"],
      location,
      ["sequence_id", "Contig"],
      interval,
      ["category", "Category"],
      ["repeat_consensus", "Repeat consensus", "sequence"],
    ],
  },
  {
    id: "repeats",
    label: "Repeats",
    description:
      "Individual repeat occurrences, including sequence variation. These are not counts of distinct repeat families.",
    columns: [
      ["array_id", "Array", "array"],
      location,
      ["ordinal", "Position in array"],
      interval,
      sequence,
    ],
  },
  {
    id: "spacers",
    label: "Spacers",
    description:
      "Individual spacer occurrences. Deletions retain their interbase boundary and are separate from observed sequences.",
    columns: [
      ["array_id", "Array", "array"],
      location,
      ["ordinal", "Position in array"],
      ["kind", "Kind"],
      interval,
      sequence,
    ],
  },
  {
    id: "cas_genes",
    label: "Cas genes",
    description:
      "Predicted genes and caller annotations. A Cas system in the same genome does not establish an association with an array.",
    columns: [
      ["call_id", "Gene"],
      location,
      ["annotation", "Annotation"],
      ["subtype", "Subtype"],
      interval,
      ["strand", "Strand"],
      tool,
      ["gene_evidence", "Evidence", "evidence"],
    ],
  },
  {
    id: "cas_systems",
    label: "Cas systems",
    description:
      "Predicted Cas systems, with classification as reported by the caller. Missing confidence remains unknown.",
    columns: [
      ["call_id", "System"],
      location,
      ["cas_type", "Type"],
      ["subtype", "Subtype"],
      interval,
      ["confidence", "Caller confidence"],
      tool,
      ["system_evidence", "Evidence", "evidence"],
    ],
  },
];

export function catalogEvidence(row, kind) {
  const fields =
    kind === "gene_evidence"
      ? [
          ["sequence_id", "Contig"],
          ["tool_version", "Caller version"],
          ["profile", "Profile"],
          ["cas_class", "Class"],
          ["cas_type", "Type"],
          ["score", "Caller score"],
          ["score_kind", "Score kind"],
          ["protein_length", "Protein length (aa)"],
          ["translation_table", "Translation table"],
          ["is_cas", "Caller Cas flag"],
          ["cas_system_id", "Predicted system record"],
        ]
      : [
          ["sequence_id", "Contig"],
          ["tool_version", "Caller version"],
          ["cas_class", "Class"],
          ["confidence_kind", "Confidence kind"],
          ["experimental", "Caller experimental flag"],
          ["missing_gene_count", "Missing gene count"],
        ];
  return fields.map(([key, label]) => ({ label, value: catalogValue(row[key]) }));
}

export function catalogNumber(value) {
  return typeof value === "number" && Number.isFinite(value)
    ? value.toLocaleString("en-US")
    : "Not reported";
}

export function catalogValue(value) {
  return value === null || value === undefined || value === "" ? "Not reported" : String(value);
}

export function catalogInterval(row) {
  if (row.kind === "deletion") {
    return row.boundary_position === null || row.boundary_position === undefined
      ? "Deletion; boundary not reported"
      : `Boundary ${catalogNumber(row.boundary_position)} (0-based)`;
  }
  if (row.start === null || row.start === undefined || row.end === null || row.end === undefined)
    return "Not reported";
  return `${catalogNumber(row.start)}–${catalogNumber(row.end)}`;
}

export function catalogCell(row, key) {
  if (key === "interval") return catalogInterval(row);
  if (key === "strand" && !["+", "-"].includes(row.strand)) return "Unknown";
  if (key === "sequence" && row.kind === "deletion") return "No sequence (deletion)";
  return Number.isInteger(row[key]) ? catalogNumber(row[key]) : catalogValue(row[key]);
}

export function catalogCsv(rows, columns) {
  // Avoid spreadsheet formula execution while preserving the displayed biological values.
  const quote = (value) => {
    const text = String(value ?? "");
    return `"${(/^[=+\-@\t\r]/u.test(text) ? `'${text}` : text).replaceAll('"', '""')}"`;
  };
  return [
    columns.map((column) => quote(column[1])).join(","),
    ...rows.map((row) =>
      columns.map(([key]) => quote(key === "interval" ? catalogInterval(row) : row[key])).join(","),
    ),
  ].join("\r\n");
}

export function catalogFasta(rows, entity) {
  return rows
    .flatMap((row) => {
      const sequence = entity === "arrays" ? row.repeat_consensus : row.sequence;
      if (row.kind === "deletion" || typeof sequence !== "string" || !/^[A-Za-z]+$/u.test(sequence))
        return [];
      const header = [entity, row.id, row.accession, row.sequence_id, row.ordinal]
        .filter((value) => value !== null && value !== undefined)
        .join("|")
        .replace(/[\r\n>]/gu, "_");
      return [`>${header}${entity === "arrays" ? " repeat_consensus" : ""}\n${sequence}\n`];
    })
    .join("");
}
