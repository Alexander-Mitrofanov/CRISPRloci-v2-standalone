import { asList, asRecord } from "./values.js";

export const text = (value) => (typeof value === "string" ? value.slice(0, 2000) : null);
export const count = (value) => (Number.isSafeInteger(value) && value >= 0 ? value : null);
export const metric = (value) =>
  typeof value === "number" && Number.isFinite(value) ? value : null;
export const probability = (value) =>
  metric(value) !== null && value >= 0 && value <= 1 ? value : null;
export const strand = (value) => (["+", "-"].includes(value) ? value : null);
export const flag = (value) => (typeof value === "boolean" ? value : null);
export const rows = (value, limit = 100) => asList(value).slice(0, limit).filter(asRecord);
export const labels = (value, limit = 120) =>
  asList(value)
    .slice(0, limit)
    .map(text)
    .filter((value) => value !== null);
export const counts = (value, keys) =>
  Object.fromEntries(keys.map((key) => [key, count(value?.[key])]));
export const interval = (value) => ({ start: count(value.start), end: count(value.end) });
export const available = (value) => value.available === true && value.execution_completed !== false;
export const casRelations = (value) =>
  rows(value).map((row) => ({
    intervalId: text(row.interval_id),
    sourceId: text(row.source_id),
    nearestCassetteCount: count(row.nearest_cassette_count),
    nearestCassettesTruncated: row.nearest_cassettes_truncated === true,
    nearestCassettes: rows(row.nearest_cassettes, 32).map((cassette) => ({
      id: text(cassette.cassette_id),
      ...interval(cassette),
      overlapBases: count(cassette.overlap_bases),
      interveningBases: count(cassette.intervening_bases),
    })),
  }));
