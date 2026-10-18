import { asRecord } from "./values.js";
import {
  available,
  casRelations,
  count,
  counts,
  interval,
  labels,
  metric,
  probability,
  rows,
  strand,
  text,
} from "./contextValues.js";

export function normalizeDefenseContext(value) {
  const doc = asRecord(value) || {};
  return {
    available: available(doc),
    casContext: casRelations(doc.cas_context),
    counts: counts(doc.counts, ["sources", "genes", "systems", "components", "hits"]),
    systemsTruncated: doc.systems_truncated === true,
    componentsTruncated: doc.components_truncated === true,
    hitsTruncated: doc.hits_truncated === true,
    systems: rows(doc.systems).map((row) => ({
      id: text(row.id),
      source_id: text(row.source_id),
      type: text(row.type),
      subtype: text(row.subtype),
      model: text(row.model),
      gene_ids: labels(row.gene_ids),
      wholeness: probability(row.wholeness),
      score: metric(row.score),
      ...interval(row),
    })),
    genes: rows(doc.genes, 500).map((row) => ({
      id: text(row.id),
      source_id: text(row.source_id),
      ...interval(row),
      strand: row.strand === 1 ? "+" : row.strand === -1 ? "-" : strand(row.strand),
      position: count(row.position),
      translation_table: count(row.translation_table),
    })),
    components: rows(doc.components, 500).map((row) => ({
      system_id: text(row.system_id),
      gene_id: text(row.gene_id),
      profile: text(row.profile),
      status: text(row.status),
      score: metric(row.score),
      i_evalue: metric(row.i_evalue),
      profile_coverage: probability(row.profile_coverage),
    })),
    hits: rows(doc.hits, 500).map((row) => ({
      gene_id: text(row.gene_id),
      profile: text(row.profile),
      i_evalue: metric(row.i_evalue),
      score: metric(row.score),
      profile_coverage: probability(row.profile_coverage),
      sequence_coverage: probability(row.sequence_coverage),
      protein_start: count(row.protein_start),
      protein_end: count(row.protein_end),
    })),
  };
}
