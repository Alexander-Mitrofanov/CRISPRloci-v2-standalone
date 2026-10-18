import { asRecord } from "./values.js";
import {
  available,
  casRelations,
  counts,
  interval,
  metric,
  rows,
  strand,
  text,
} from "./contextValues.js";

export function normalizeInsertionContext(value) {
  const doc = asRecord(value) || {};
  return {
    available: available(doc),
    counts: counts(doc.counts, ["sources", "elements"]),
    casContext: casRelations(doc.cas_context),
    truncated: doc.elements_truncated === true || doc.truncated === true,
    elements: rows(doc.elements).map((element) => ({
      id: text(element.id),
      source_id: text(element.source_id),
      ...interval(element),
      family: text(element.family),
      native_type: ["c", "p"].includes(element.native_type) ? element.native_type : null,
      contig_edge: element.contig_edge === true,
      evalue_this_copy:
        metric(element.evalue_this_copy) !== null && element.evalue_this_copy >= 0
          ? element.evalue_this_copy
          : null,
      transposase_orf: asRecord(element.transposase_orf)
        ? { ...interval(element.transposase_orf), strand: strand(element.transposase_orf.strand) }
        : null,
      terminal_inverted_repeats: asRecord(element.terminal_inverted_repeats),
    })),
  };
}
