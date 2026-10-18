import { asList, asRecord } from "./values.js";
const text = (v) => (typeof v === "string" ? v : null);
const count = (v) => (Number.isSafeInteger(v) && v >= 0 ? v : null);
export function normalizeArrayCompare(value) {
  const doc = asRecord(value) || {};
  const arrays = asList(doc.arrays)
    .filter(asRecord)
    .slice(0, 50)
    .map((a) => ({
      id: text(a.id),
      sourceId: text(a.source_id),
      arrayId: text(a.array_id),
      strand: ["+", "-"].includes(a.strand) ? a.strand : null,
      observedCount: count(a.observed_spacer_count),
      evidenceStatus: text(a.evidence_status),
      truncated: a.spacers_truncated === true,
      spacers: asList(a.spacers)
        .filter(asRecord)
        .slice(0, 50)
        .map((s) => ({
          id: text(s.id),
          ordinal: count(s.ordinal),
          sequence: text(s.sequence),
          symbol:
            typeof s.symbol === "string" && /^s_[a-f0-9]{64}$/.test(s.symbol) ? s.symbol : null,
          deletion: s.status === "deletion",
          start: count(s.start),
          end: count(s.end),
        })),
    }));
  const ids = new Set(arrays.map((a) => a.id));
  const edges = asList(doc.edges)
    .filter(asRecord)
    .slice(0, 100)
    .filter(
      (e) =>
        ids.has(e.a) &&
        ids.has(e.b) &&
        count(e.nshared) > 0 &&
        typeof e.jaccard === "number" &&
        Number.isFinite(e.jaccard) &&
        e.jaccard > 0 &&
        e.jaccard <= 1,
    )
    .map((e) => ({ a: e.a, b: e.b, shared: e.nshared, jaccard: e.jaccard }));
  return {
    available: doc.available === true && doc.execution_completed === true,
    coordinateSystem: text(doc.coordinate_system),
    scoreSemantics: text(doc.score_semantics),
    arrays,
    edges,
    isolates: asList(doc.isolates).filter((id) => typeof id === "string"),
    counts: Object.fromEntries(
      ["arrays", "spacer_occurrences", "deletions", "edges", "components", "isolates"].map((k) => [
        k,
        count(doc.counts?.[k]),
      ]),
    ),
    identity: ["exact", "reverse-complement"].includes(doc.options?.identity)
      ? doc.options.identity
      : null,
    minShared: count(doc.options?.min_shared),
    truncated:
      doc.arrays_truncated === true ||
      doc.edges_truncated === true ||
      arrays.some((a) => a.truncated),
  };
}
