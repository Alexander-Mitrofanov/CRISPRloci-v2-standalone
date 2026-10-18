import { asList, asRecord } from "./values.js";

const text = (value) => (typeof value === "string" && value ? value : null);
const count = (value) => (Number.isSafeInteger(value) && value >= 0 ? value : null);
const position = (value) => (count(value) > 0 ? value : null);
const hash = (value) => (typeof value === "string" && /^[a-f0-9]{64}$/u.test(value) ? value : null);

function sourceUrl(value) {
  try {
    const url = new URL(value);
    return url.protocol === "https:" && !url.username && !url.password ? url.href : null;
  } catch {
    return null;
  }
}

export function normalizeProtospacer(value) {
  const report = asRecord(value) || {};
  const reference = asRecord(report.reference) || {};
  const counts = asRecord(report.counts) || {};
  const rows = asList(report.matches).filter(asRecord);
  return {
    available: report.available === true,
    search_complete: typeof report.search_complete === "boolean" ? report.search_complete : null,
    outcome: text(report.outcome),
    reference: {
      name: text(reference.name),
      source_url: sourceUrl(reference.source_url),
      source_release: text(reference.source_release),
      source_file_sha256: hash(reference.source_file_sha256),
      database_sha256: hash(reference.database_sha256),
      record_count: count(reference.record_count),
      eligible_record_count: count(reference.eligible_record_count),
      skipped_ambiguous: count(reference.skipped_ambiguous),
      skipped_length: count(reference.skipped_length),
      min_length: count(reference.min_length),
      max_length: count(reference.max_length),
    },
    counts: Object.fromEntries(
      [
        "matches",
        "matched_spacers",
        "targets",
        "target_bases",
        "seed_windows",
        "ambiguous_seed_windows",
        "candidates_checked",
      ].map((key) => [key, count(counts[key])]),
    ),
    matches: rows.slice(0, 100).map((row) => ({
      spacer_id: text(row.spacer_id),
      spacer_sequence: text(row.spacer_sequence),
      target_id: text(row.target_id),
      target_start: position(row.target_start),
      target_end: position(row.target_end),
      relative_strand: ["+", "-", "both"].includes(row.relative_strand)
        ? row.relative_strand
        : null,
      query_length: count(row.query_length),
      identities: count(row.identities),
    })),
    matches_truncated:
      report.matches_truncated === true ||
      rows.length > 100 ||
      (count(counts.matches) != null && counts.matches > rows.length),
  };
}
