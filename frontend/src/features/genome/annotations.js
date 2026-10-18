import { coordinates, rows } from "./contract.js";

export function arrayTracks(data, c) {
  coordinates(data);
  const arrays = c.track("arrays", "CRISPR arrays", "#0f766e");
  const repeats = c.track("repeats", "Repeat occurrences", "#65758b");
  const spacers = c.track("spacers", "Spacer occurrences", "#a76716");
  for (const row of rows(data.arrays)) {
    c.add(arrays, row, {
      start: row.interval?.start,
      end: row.interval?.end,
      details: [
        ["Category", row.category],
        ["Raw model score", row.model_score],
        ["Exact boundary", row.boundary?.exact === true ? "Reported" : "Not established"],
      ],
    });
    for (const [key, target] of [
      ["repeats", repeats],
      ["spacers", spacers],
    ]) {
      for (const unit of rows(row[key])) {
        const deletion = unit.status === "deletion";
        c.add(
          target,
          { ...unit, source_id: row.source_id },
          {
            start: deletion ? unit.boundary_position : unit.start,
            end: deletion ? unit.boundary_position : unit.end,
            halfOpen: deletion,
            point: deletion,
            parentId: row.id,
            strand: null,
            label: `${row.id} · ${deletion ? "deletion slot" : key === "repeats" ? "repeat" : "spacer"} ${unit.ordinal}`,
            details: [
              ["Array", row.id],
              ["Ordinal", unit.ordinal],
              ["Occurrence", unit.id],
              [
                "Evidence",
                deletion
                  ? "Deletion slot; no observed DNA"
                  : unit.status === "not_source_verified"
                    ? "Occurrence sequence is not source-verified"
                    : "Reported occurrence in source-forward order",
              ],
              ["Sequence status", unit.status?.replaceAll("_", " ") ?? null],
            ],
          },
        );
      }
    }
  }
}

export function casTracks(data, c) {
  coordinates(data);
  const systems = c.track("cassettes", "Cas cassettes", "#738ab4");
  const genes = c.track("cas", "Cas genes", "#3158a5");
  for (const row of rows(data.cassettes))
    c.add(systems, row, {
      details: [
        ["Type", row.type],
        ["Subtype", row.subtype],
        ["Class", row.class],
        ["Model confidence", row.confidence],
      ],
    });
  for (const row of rows(data.genes))
    c.add(genes, row, {
      label: row.profile ? `${row.profile} · ${row.id}` : row.id,
      details: [
        ["Gene", row.id],
        ["Profile", row.profile],
        ["Raw score margin", row.score_margin],
      ],
    });
}

export function tracrTracks(data, c) {
  coordinates(data);
  const track = c.track("tracrrna", "tracrRNA candidates", "#9b5a21");
  for (const row of rows(data.candidates))
    c.add(track, row, {
      details: [
        ["Evidence", row.evidence_class],
        ["Boundary", row.boundary_kind],
        ["Score kind", row.score_kind],
        ["Raw score", row.score],
        ["E-value", row.evalue],
      ],
    });
}

export function flankTracks(data, c, includeArrays = false) {
  coordinates(data, "zero_based_half_open");
  if (includeArrays) {
    const arrays = c.track("arrays", "CRISPR arrays", "#0f766e");
    for (const row of rows(data.arrays))
      for (const segment of rows(row.segments)) {
        c.add(arrays, row, {
          ...segment,
          halfOpen: true,
          strand: null,
          details: [
            ["Category", row.category],
            ["Exact boundary", "Not established"],
          ],
        });
      }
  }
  const track = c.track("flanks", "Adjacent context windows", "#a989b4");
  for (const row of rows(data.contexts))
    for (const segment of rows(row.segments)) {
      c.add(track, row, {
        ...segment,
        halfOpen: true,
        strand: null,
        parentId: row.array_id,
        label: `${row.array_id} · ${row.side} flank`,
        details: [
          ["Array", row.array_id],
          ["Side", row.side],
          ["Sequence extraction strand", row.strand],
          ["Clipping", row.truncated_reason],
          ["Interpretation", "Analysis window; leader prediction not performed"],
        ],
      });
    }
}
