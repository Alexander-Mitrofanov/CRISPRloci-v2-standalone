import { asList, firstDefined, firstRecord } from "./values.js";

function normalizeUnit(value, index) {
  const row = firstRecord(value) || {};
  const interval = firstRecord(row.source_interval) || row;
  return {
    ordinal: firstDefined(row, "ordinal", "index") ?? index + 1,
    sequence:
      row.kind === "deletion" ? null : typeof value === "string" ? value : (row.sequence ?? null),
    sequence_status: row.kind === "deletion" ? "deletion" : (row.sequence_status ?? null),
    start: interval.start ?? null,
    end: interval.end ?? null,
  };
}

export function normalizeDetectionArray(value) {
  const row = firstRecord(value) || {};
  return {
    ...row,
    array_id: firstDefined(row, "array_id", "id", "Name"),
    category: firstDefined(row, "category", "Category"),
    end: firstDefined(row, "end", "End"),
    model_score: firstDefined(row, "model_score", "score", "Confidence score"),
    source_id: firstDefined(row, "source_id", "record_id", "sequence_id", "Name"),
    spacer_count: firstDefined(row, "spacer_count", "Number of spacers"),
    start: firstDefined(row, "start", "Start"),
    strand: firstDefined(row, "strand", "Strand"),
    repeats: asList(firstDefined(row, "repeats", "repeat_sequences")).map(normalizeUnit),
    spacers: asList(firstDefined(row, "spacers", "spacer_sequences")).map(normalizeUnit),
    sequences_truncated: row.sequences_truncated === true,
  };
}

export function normalizeDetection(summary) {
  const source = firstRecord(summary?.detection, summary) || {};
  const arrays = asList(firstDefined(source, "arrays", "detected_arrays"));
  const categoryCounts = firstRecord(firstDefined(source, "category_counts", "categories"));
  return {
    ...source,
    arrays: arrays.map(normalizeDetectionArray),
    category_counts: categoryCounts || {},
  };
}
