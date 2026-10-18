import { asList, asRecord } from "./values.js";
const text = (v) => (typeof v === "string" ? v : null);
const count = (v) => (Number.isSafeInteger(v) && v >= 0 ? v : null);
const label = (v) => {
  const r = asRecord(v) || {};
  return {
    status: text(r.status),
    value: text(r.value),
    values: asList(r.values).filter((v) => typeof v === "string"),
    missing: count(r.missing_reference_count),
  };
};
export function normalizeRepeatMap(value) {
  const r = asRecord(value) || {},
    reference = asRecord(r.reference) || {};
  return {
    available: r.available === true && r.execution_completed === true,
    namespace: text(reference.namespace),
    maxDistance: count(r.policy?.max_distance),
    topK: count(r.policy?.top_k),
    referenceSha256: text(reference.file_sha256),
    truncated: r.results_truncated === true,
    counts: Object.fromEntries(
      ["query_count", "exact_count", "near_count", "no_hit_count", "unsupported_count"].map((k) => [
        k,
        count(asRecord(r.counts)?.[k]),
      ]),
    ),
    results: asList(r.results)
      .filter(asRecord)
      .slice(0, 100)
      .map((row) => ({
        id: text(row.id),
        sequence: text(row.sequence),
        status: text(row.status),
        distance: count(row.best_distance),
        hitCount: count(row.best_hit_count),
        reportedHitCount: count(row.reported_hit_count),
        truncated: row.hits_truncated === true,
        previewTruncated:
          asList(row.hits).length > 100 ||
          (count(row.reported_hit_count) !== null &&
            row.reported_hit_count > Math.min(asList(row.hits).length, 100)),
        family: label(asRecord(row.annotations)?.family),
        motif: label(asRecord(row.annotations)?.motif),
        hits: asList(row.hits)
          .filter(asRecord)
          .slice(0, 100)
          .map((h) =>
            Object.fromEntries(
              ["reference_id", "orientation", "query_aligned", "reference_aligned"].map((k) => [
                k,
                text(h[k]),
              ]),
            ),
          ),
      })),
  };
}
