import { createHash } from "node:crypto";

export function genomeFixture() {
  const records = [
    { id: "alpha", sequence: "ACGT".repeat(500) },
    { id: "beta", sequence: "TGCA".repeat(500) },
  ];
  const sources = records.map((row) => ({
    id: row.id,
    length: row.sequence.length,
    sha256: createHash("sha256").update(row.sequence).digest("hex"),
    topology: "unknown",
  }));
  const forward = "1-based-inclusive-source-forward";
  const data = {
    sources,
    detection: {
      coordinate_system: forward,
      arrays: sources.map((source, index) => ({
        id: "array-1",
        source_id: source.id,
        interval: { start: 422, end: 580 },
        category: "Possible",
        model_score: index === 0 ? 0 : null,
        strand: null,
        boundary: { exact: false },
        repeats: [{ id: "repeat-1", ordinal: 1, start: 422, end: 451 }],
        spacers: [
          { id: "spacer-1", ordinal: 1, start: 452, end: 485 },
          { id: "spacer-2", ordinal: 2, start: 516, end: 549 },
        ],
      })),
    },
    cas: {
      coordinate_system: forward,
      genes: [
        {
          id: "cas-1",
          source_id: "alpha",
          start: 700,
          end: 999,
          strand: "-",
          profile: "Cas9",
          score_margin: null,
        },
      ],
      cassettes: [],
    },
    tracrrna: {
      coordinate_system: forward,
      candidates: [
        {
          id: "rna-1",
          source_id: "alpha",
          start: 600,
          end: 640,
          strand: "+",
          score: 0,
          evalue: null,
        },
      ],
    },
    flanks: {
      coordinate_system: "zero_based_half_open",
      contexts: [
        {
          id: "left",
          array_id: "array-1",
          source_id: "alpha",
          side: "left",
          strand: "+",
          segments: [{ start: 384, end: 421 }],
        },
        {
          id: "right",
          array_id: "array-1",
          source_id: "alpha",
          side: "right",
          strand: "-",
          segments: [{ start: 580, end: 617 }],
        },
      ],
    },
  };
  return {
    result: {
      schema: { name: "crisprloci.module-result", version: "1.0.0" },
      module: "loci",
      complete: true,
      status: "completed",
      data,
    },
    input: {
      schema: { name: "crisprloci.sequence-set", version: "1.0.0" },
      kind: "target",
      molecule: "DNA",
      complete: true,
      records,
    },
  };
}

export const genomeArtifacts = [
  { name: "result.json", artifact_id: "result-artifact", size_bytes: 5000 },
  { name: "input-sequences.json", artifact_id: "input-artifact", size_bytes: 5000 },
];
