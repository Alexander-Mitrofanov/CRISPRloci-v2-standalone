import { asList, asRecord } from "./values.js";
const text = (value) => (typeof value === "string" ? value : null);
const count = (value) => (Number.isSafeInteger(value) && value >= 0 ? value : null);
export function normalizeRepeatType(value) {
  const doc = asRecord(value) || {};
  const model = asRecord(doc.model) || {};
  return {
    available: doc.available === true && doc.execution_completed === true,
    modelId: text(model.id),
    modelSha256: text(asRecord(model.files)?.["xgb_repeats.json"]),
    labels: asList(model.labels)
      .filter((v) => typeof v === "string")
      .slice(0, 100),
    threshold: Number.isFinite(doc.policy?.threshold) ? doc.policy.threshold : null,
    truncated: doc.results_truncated === true,
    counts: Object.fromEntries(
      ["records", "predicted", "uncertain", "unsupported"].map((k) => [k, count(doc.counts?.[k])]),
    ),
    results: asList(doc.results)
      .filter(asRecord)
      .slice(0, 1000)
      .map((row) => ({
        id: text(row.id),
        sequence: text(row.sequence),
        status: text(row.status),
        prediction: text(row.prediction),
        winner: text(row.winner),
        reason: text(row.reason),
        nativeScore:
          typeof row.native_score === "number" &&
          Number.isFinite(row.native_score) &&
          row.native_score >= 0 &&
          row.native_score <= 1
            ? row.native_score
            : null,
      })),
  };
}
