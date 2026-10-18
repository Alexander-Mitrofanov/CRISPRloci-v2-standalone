import { asList, asRecord } from "./values.js";

const count = (value) => (Number.isInteger(value) && value >= 0 ? value : null);
const text = (value) => (typeof value === "string" ? value : null);
const segments = (value) =>
  asList(value).filter(
    (row) =>
      Number.isInteger(row?.start) &&
      row.start >= 0 &&
      Number.isInteger(row.end) &&
      row.end > row.start,
  );

export function normalizeLeader(value) {
  const report = asRecord(value) || {};
  return {
    status: report.status || "not_requested",
    array_count: count(report.array_count),
    context_count: count(report.context_count),
    flank_length: count(report.flank_length),
    truncated: report.truncated === true,
    contexts: asList(report.contexts)
      .filter(asRecord)
      .map((row) => ({
        id: text(row.id),
        source_id: text(row.source_id),
        array_id: text(row.array_id),
        side: ["left", "right"].includes(row.side) ? row.side : null,
        strand: ["+", "-"].includes(row.strand) ? row.strand : null,
        segments: segments(row.segments),
        available_length: count(row.available_length),
        requested_length: count(row.requested_length),
        truncated_reason: text(row.truncated_reason),
        sequence: text(row.sequence),
        observed_repeat: text(row.observed_repeat),
      })),
  };
}
