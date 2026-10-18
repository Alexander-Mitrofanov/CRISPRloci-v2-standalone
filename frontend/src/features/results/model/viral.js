import { asList, asRecord } from "./values.js";
const text = (v) => (typeof v === "string" ? v : null);
const count = (v) => (Number.isSafeInteger(v) && v >= 0 ? v : null);
export function normalizeViral(value) {
  const r = asRecord(value) || {},
    reference = asRecord(r.reference) || {},
    counts = asRecord(r.counts) || {};
  return {
    available: r.available === true,
    execution_completed: typeof r.execution_completed === "boolean" ? r.execution_completed : null,
    retrieval_completeness: text(r.retrieval_completeness),
    outcome: text(r.outcome),
    reference: Object.fromEntries(
      ["name", "source_release", "reference_sha256"]
        .map((k) => [k, text(reference[k])])
        .concat([["record_count", count(reference.record_count)]]),
    ),
    counts: Object.fromEntries(
      [
        "query_occurrences",
        "eligible_query_occurrences",
        "matched_query_occurrences",
        "matched_subject_accessions",
        "accepted_hsps",
        "raw_hsps",
      ].map((k) => [k, count(counts[k])]),
    ),
    max_substitutions: count(asRecord(r.policy)?.max_substitutions),
    matches: asList(r.matches)
      .filter(asRecord)
      .slice(0, 100)
      .map((h) => ({
        ...Object.fromEntries(
          ["query_id", "query_sequence", "subject_accession", "subject_title"].map((k) => [
            k,
            text(h[k]),
          ]),
        ),
        ...Object.fromEntries(
          ["subject_start", "subject_end", "substitutions", "query_length"].map((k) => [
            k,
            count(h[k]),
          ]),
        ),
        relative_strand: ["+", "-"].includes(h.relative_strand) ? h.relative_strand : null,
      })),
    candidate_viruses: asList(r.candidate_viruses)
      .filter(asRecord)
      .slice(0, 100)
      .map((h) => ({
        accession: text(h.accession),
        title: text(h.title),
        unique_query_sequence_count: count(h.unique_query_sequence_count),
        query_occurrence_count: count(h.query_occurrence_count),
        location_count: count(h.location_count),
      })),
    queries: asList(r.queries)
      .filter(asRecord)
      .slice(0, 100)
      .map((q) => ({
        id: text(q.original_id),
        eligible: typeof q.eligible === "boolean" ? q.eligible : null,
        reasons: asList(q.skip_reasons).filter((v) => typeof v === "string"),
      })),
    truncated:
      r.matches_truncated === true ||
      r.candidates_truncated === true ||
      r.queries_truncated === true ||
      asList(r.matches).length > 100,
  };
}
